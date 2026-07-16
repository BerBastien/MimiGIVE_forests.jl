# Integration tests.
#   Part 1 (always runs): a standalone 3-component chain + a Mimi MarginalModel,
#           exercising the marginal-damage plumbing WITHOUT the MimiGIVE stack.
#   Part 2 (opt-in): full MimiGIVE build, enabled with ENV["FOREST_RUN_FULL"]=="true".

using Test, Mimi

for f in ("forest_spatial_socioeconomics.jl", "forest_area_response.jl", "forest_ecosystem_services.jl")
    sym = f == "forest_spatial_socioeconomics.jl" ? :SpatialSocioeconomics :
          f == "forest_area_response.jl"          ? :ForestAreaResponse :
                                                    :ForestEcosystemServices
    isdefined(@__MODULE__, sym) || include(joinpath(@__DIR__, "..", "src", f))
end

# Build a small standalone chain: Spatial -> Area -> ES, wired like the real model.
function _build_standalone(gmst::Vector{Float64})
    years = 2020:2035
    nt = length(years); nc = 2
    m = Model()
    set_dimension!(m, :time, years)
    set_dimension!(m, :country, ["X", "Y"])
    set_dimension!(m, :impact_region, ["A", "B", "C"])   # A,B in X ; C in Y

    add_comp!(m, SpatialSocioeconomics, :Spatial; first = 2020)
    add_comp!(m, ForestAreaResponse, :Area; first = 2020, after = :Spatial)
    add_comp!(m, ForestEcosystemServices, :ES; first = 2020, after = :Area)

    # --- Spatial ---
    popc = [10.0 + 0.1 * (t - 1) + (c - 1) * 10.0 for t in 1:nt, c in 1:nc]
    gdpc = [100.0 + (t - 1) + (c - 1) * 100.0 for t in 1:nt, c in 1:nc]
    update_param!(m, :Spatial, :population_country, popc)
    update_param!(m, :Spatial, :gdp_country, gdpc)
    update_param!(m, :Spatial, :population_share, [0.6, 0.4, 1.0])
    update_param!(m, :Spatial, :gdp_share, [0.7, 0.3, 1.0])
    update_param!(m, :Spatial, :impact_region_country_index, [1, 1, 2])

    # --- Area ---
    update_param!(m, :Area, :gmst, gmst)
    update_param!(m, :Area, :beta1, [-2.0, -1.0, -3.0])
    update_param!(m, :Area, :beta2, [0.0, 0.0, 0.0])
    update_param!(m, :Area, :intercept, [0.0, 0.0, 0.0])
    update_param!(m, :Area, :baseline_forest_area, [10.0, 5.0, 8.0])
    update_param!(m, :Area, :total_impact_region_area, [0.0, 0.0, 0.0])
    update_param!(m, :Area, :delta_gmst_fit_min, [-Inf, -Inf, -Inf])
    update_param!(m, :Area, :delta_gmst_fit_max, [Inf, Inf, Inf])
    update_param!(m, :Area, :baseline_year, 2022)
    update_param!(m, :Area, :change_units_code, 1)
    update_param!(m, :Area, :apply_lower_clip, true)
    update_param!(m, :Area, :apply_upper_clip, false)

    # --- ES ---
    connect_param!(m, :ES => :projected_forest_area, :Area => :projected_forest_area)
    connect_param!(m, :ES => :population_ir, :Spatial => :population_ir)
    connect_param!(m, :ES => :gdp_ir, :Spatial => :gdp_ir)
    update_param!(m, :ES, :baseline_forest_area, [10.0, 5.0, 8.0])
    update_param!(m, :ES, :impact_region_country_index, [1, 1, 2])
    update_param!(m, :ES, :es_value_scale, 1.0)
    update_param!(m, :ES, :es_value_spec_code, 1)
    return m
end

@testset "Marginal-model plumbing (standalone)" begin
    base_gmst  = [1.0 + 0.05 * (t - 1) for t in 1:16]
    pulse_gmst = copy(base_gmst)
    for t in 7:16                      # divergence starts at year index 7 (== 2026)
        pulse_gmst[t] += 0.2
    end

    base = _build_standalone(base_gmst)
    mm = Mimi.create_marginal_model(base, 1.0)
    update_param!(mm.modified, :Area, :gmst, pulse_gmst)
    run(mm)

    # base model: at 2022 (index 3), dT = 0 so projected == baseline
    @test Float64.(mm.base[:Area, :projected_forest_area][3, :]) ≈ [10.0, 5.0, 8.0]

    md_area = mm[:Area, :projected_forest_area]     # (modified - base)/delta
    md_dmg  = mm[:ES, :es_damage_country]
    md_glob = mm[:ES, :es_damage_global]

    # identical outcomes BEFORE the pulse diverges (indices 1..6)
    @test all(abs.(md_area[1:6, :]) .< 1e-10)
    @test all(abs.(md_dmg[1:6, :]) .< 1e-10)

    # small but nonzero marginal forest change AFTER divergence
    @test any(abs.(md_area[7:end, :]) .> 0.0)
    @test any(abs.(md_dmg[7:end, :]) .> 0.0)

    # country marginal damages sum to global marginal damage, every year
    for ti in 1:16
        @test sum(md_dmg[ti, :]) ≈ md_glob[ti]
    end
end

# ---------------------------------------------------------------------------
# Opt-in full MimiGIVE integration (heavy: needs the full stack + instantiate).
# ---------------------------------------------------------------------------
if get(ENV, "FOREST_RUN_FULL", "false") == "true"
    @testset "Full MimiGIVE integration" begin
        forest_root = normpath(joinpath(@__DIR__, ".."))
        processed = joinpath(forest_root, "data", "processed")
        if !isfile(joinpath(processed, "impact_region_country_crosswalk.csv"))
            include(joinpath(forest_root, "scripts", "make_synthetic_fixtures.jl"))
        end
        include(joinpath(forest_root, "ForestGIVE.jl"))
        M = getfield(@__MODULE__, :ForestGIVE)

        built = M.build_forest_give_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245")
        run(built.model)

        i2022 = findfirst(==(2022), collect(1750:2300))
        @test built.model[:ForestAreaResponse, :delta_gmst][i2022] == 0.0
        proj2022 = Float64.(built.model[:ForestAreaResponse, :projected_forest_area][i2022, :])
        @test proj2022 ≈ built.inputs.baseline_forest_area                 # dT=0 => baseline
        @test all(built.model[:ForestEcosystemServices, :es_damage_ir][i2022, :] .== 0.0)

        scc = M.compute_forest_scc(built; year = 2030, gas = :CO2, verbose = false)
        @test isfinite(scc.scc_global)
        @test sum(scc.scc_country) ≈ scc.scc_global
    end
else
    @info "Skipping full MimiGIVE integration test. Set ENV[\"FOREST_RUN_FULL\"]=\"true\" to enable."
end
