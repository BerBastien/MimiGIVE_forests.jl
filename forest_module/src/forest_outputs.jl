using Mimi, DataFrames, CSVFiles, FileIO

# =============================================================================
# forest_outputs.jl  --  tidy tabular outputs + diagnostics
# -----------------------------------------------------------------------------
# All outputs are plain CSVs (long/tidy format), easy to read back into R with
# readr::read_csv(). Geometry-dependent outputs (GeoPackage, maps) are produced
# separately by scripts/make_diagnostic_maps.py, which joins these CSVs back onto
# the impact-region polygons (keeping geometry out of the Mimi runtime).
# =============================================================================

_f(x) = x === missing ? NaN : Float64(x)

function _year_indices(m, years)
    tk = collect(Mimi.dim_keys(m, :time))
    idx = [findfirst(==(y), tk) for y in years]
    any(isnothing, idx) && error("Requested output years not all in model time index: $years")
    return Int.(idx)
end

_model(built) = built isa NamedTuple ? built.model : built

"""
    write_impact_region_output(built; years=2020:10:2100, output_dir=...) -> path

Tidy impact-region table. NOTE: full annual output can be very large
(n_years * n_impact_regions rows); the default subsamples years. Pass `years` to
control it (e.g. a single year for a map layer).
"""
function write_impact_region_output(built; years = 2020:10:2100,
                                    output_dir::String = ForestConfig.OUTPUT_DIRECTORY)
    mkpath(output_dir)
    m = _model(built); inp = built.inputs
    tk = collect(Mimi.dim_keys(m, :time)); idx = _year_indices(m, years)

    gmst   = m[:temperature, :T]
    dgmst  = m[:ForestAreaResponse, :delta_gmst]
    fchg   = m[:ForestAreaResponse, :forest_change]
    proj   = m[:ForestAreaResponse, :projected_forest_area]
    pop    = m[:ForestSpatialSocioeconomics, :population_ir]
    gdp    = m[:ForestSpatialSocioeconomics, :gdp_ir]
    gdppc  = m[:ForestSpatialSocioeconomics, :gdppc_ir]
    esref  = m[:ForestEcosystemServices, :es_value_reference_ir]
    esclim = m[:ForestEcosystemServices, :es_value_climate_ir]
    esdmg  = m[:ForestEcosystemServices, :es_damage_ir]

    nir = length(inp.impact_region_ids)
    rows = length(idx) * nir
    df = DataFrame(
        year = Vector{Int}(undef, rows),
        impact_region_id = Vector{String}(undef, rows),
        give_country_id  = Vector{String}(undef, rows),
        gmst = Vector{Float64}(undef, rows),
        delta_gmst = Vector{Float64}(undef, rows),
        forest_change = Vector{Float64}(undef, rows),
        baseline_forest_area = Vector{Float64}(undef, rows),
        projected_forest_area = Vector{Float64}(undef, rows),
        population_ir = Vector{Float64}(undef, rows),
        gdp_ir = Vector{Float64}(undef, rows),
        gdppc_ir = Vector{Float64}(undef, rows),
        es_value_reference = Vector{Float64}(undef, rows),
        es_value_climate = Vector{Float64}(undef, rows),
        es_damage = Vector{Float64}(undef, rows),
    )
    k = 0
    for (yi, ti) in enumerate(idx)
        y = tk[ti]
        for r in 1:nir
            k += 1
            df.year[k] = y
            df.impact_region_id[k] = inp.impact_region_ids[r]
            df.give_country_id[k]  = inp.give_country_ids[inp.impact_region_country_index[r]]
            df.gmst[k] = _f(gmst[ti])
            df.delta_gmst[k] = _f(dgmst[ti])
            df.forest_change[k] = _f(fchg[ti, r])
            df.baseline_forest_area[k] = inp.baseline_forest_area[r]
            df.projected_forest_area[k] = _f(proj[ti, r])
            df.population_ir[k] = _f(pop[ti, r])
            df.gdp_ir[k] = _f(gdp[ti, r])
            df.gdppc_ir[k] = _f(gdppc[ti, r])
            df.es_value_reference[k] = _f(esref[ti, r])
            df.es_value_climate[k] = _f(esclim[ti, r])
            df.es_damage[k] = _f(esdmg[ti, r])
        end
    end
    path = joinpath(output_dir, "forest_impact_region_output.csv")
    save(path, df)
    return path
end

"""
    write_country_output(built; years=2020:2300, country_names=nothing, output_dir=...) -> path
"""
function write_country_output(built; years = 2020:2300,
                              country_names = nothing,
                              output_dir::String = ForestConfig.OUTPUT_DIRECTORY)
    mkpath(output_dir)
    m = _model(built); inp = built.inputs
    tk = collect(Mimi.dim_keys(m, :time)); idx = _year_indices(m, years)
    countries = String.(Mimi.dim_keys(m, :country))
    names_lut = country_names === nothing ? Dict(c => c for c in countries) : country_names

    ref  = m[:ForestEcosystemServices, :es_value_reference_country]
    clim = m[:ForestEcosystemServices, :es_value_climate_country]
    dmg  = m[:ForestEcosystemServices, :es_damage_country]

    nc = length(countries)
    rows = length(idx) * nc
    df = DataFrame(
        year = Vector{Int}(undef, rows),
        give_country_id = Vector{String}(undef, rows),
        country_name = Vector{String}(undef, rows),
        es_value_reference = Vector{Float64}(undef, rows),
        es_value_climate = Vector{Float64}(undef, rows),
        es_damage = Vector{Float64}(undef, rows),
    )
    k = 0
    for ti in idx
        y = tk[ti]
        for c in 1:nc
            k += 1
            df.year[k] = y
            df.give_country_id[k] = countries[c]
            df.country_name[k] = get(names_lut, countries[c], countries[c])
            df.es_value_reference[k] = _f(ref[ti, c])
            df.es_value_climate[k] = _f(clim[ti, c])
            df.es_damage[k] = _f(dmg[ti, c])
        end
    end
    path = joinpath(output_dir, "forest_country_output.csv")
    save(path, df)
    return path
end

"""
    write_marginal_country_output(scc_res; output_dir=...) -> path

Uses the base/pulse forest damages and discount factors from `compute_forest_scc`.
"""
function write_marginal_country_output(scc_res; output_dir::String = ForestConfig.OUTPUT_DIRECTORY)
    mkpath(output_dir)
    countries = scc_res.countries
    tk = collect(1750:2300)
    year = scc_res.year; last_year = scc_res.last_year
    yi = findfirst(==(year), tk); lyi = findfirst(==(last_year), tk)
    years = tk[yi:lyi]
    df_disc = scc_res.discount_factors    # aligned to years

    base  = scc_res.forest_damage_base
    pulse = scc_res.forest_damage_pulse
    md    = scc_res.marginal_damage_country

    nc = length(countries)
    rows = length(years) * nc
    df = DataFrame(
        year = Vector{Int}(undef, rows),
        give_country_id = Vector{String}(undef, rows),
        forest_damage_base = Vector{Float64}(undef, rows),
        forest_damage_pulse = Vector{Float64}(undef, rows),
        forest_marginal_damage = Vector{Float64}(undef, rows),
        discount_factor = Vector{Float64}(undef, rows),
        discounted_forest_marginal_damage = Vector{Float64}(undef, rows),
    )
    k = 0
    for (j, ti) in enumerate(yi:lyi)
        for c in 1:nc
            k += 1
            df.year[k] = years[j]
            df.give_country_id[k] = countries[c]
            df.forest_damage_base[k] = _f(base[ti, c])
            df.forest_damage_pulse[k] = _f(pulse[ti, c])
            df.forest_marginal_damage[k] = _f(md[ti, c])
            df.discount_factor[k] = df_disc[j]
            df.discounted_forest_marginal_damage[k] = df_disc[j] * _f(md[ti, c])
        end
    end
    path = joinpath(output_dir, "forest_marginal_country_output.csv")
    save(path, df)
    return path
end

"""
    write_summary_output(built, scc_res; output_dir=...) -> (summary_path, contrib_path)
"""
function write_summary_output(built, scc_res; output_dir::String = ForestConfig.OUTPUT_DIRECTORY)
    mkpath(output_dir)
    m = _model(built)
    tk = collect(Mimi.dim_keys(m, :time))

    refg  = m[:ForestEcosystemServices, :es_value_reference_global]
    climg = m[:ForestEcosystemServices, :es_value_climate_global]
    dmgg  = m[:ForestEcosystemServices, :es_damage_global]

    dyears = 2020:2300
    idx = _year_indices(m, dyears)
    global_df = DataFrame(
        year = collect(dyears),
        es_value_reference_global = [_f(refg[i]) for i in idx],
        es_value_climate_global   = [_f(climg[i]) for i in idx],
        es_damage_global          = [_f(dmgg[i]) for i in idx],
    )
    p1 = joinpath(output_dir, "forest_global_summary.csv")
    save(p1, global_df)

    contrib = DataFrame(
        give_country_id = scc_res.countries,
        forest_scc_contribution = scc_res.scc_country,
    )
    sort!(contrib, :forest_scc_contribution, rev = true)
    p2 = joinpath(output_dir, "forest_scc_country_contribution.csv")
    save(p2, contrib)

    # one-line headline file
    headline = DataFrame(
        metric = ["forest_scc_global_placeholder_per_tonne", "pulse_year", "gas", "prtp", "eta"],
        value  = [string(scc_res.scc_global), string(scc_res.year), string(scc_res.gas),
                  string(scc_res.prtp), string(scc_res.eta)],
    )
    p3 = joinpath(output_dir, "forest_scc_headline.csv")
    save(p3, headline)

    return (p1, p2, p3)
end

"""
    write_temperature_diagnostic(built; output_dir=...) -> path

year, gmst, gmst_2022, delta_gmst, delta_gmst_squared
"""
function write_temperature_diagnostic(built; output_dir::String = ForestConfig.OUTPUT_DIRECTORY)
    mkpath(output_dir)
    m = _model(built)
    tk = collect(Mimi.dim_keys(m, :time))
    gmst = m[:temperature, :T]
    dgmst = m[:ForestAreaResponse, :delta_gmst]
    dgmst2 = m[:ForestAreaResponse, :delta_gmst_squared]
    gmst2022 = m[:ForestAreaResponse, :gmst_baseline]   # scalar
    dyears = 2020:2300
    idx = _year_indices(m, dyears)
    df = DataFrame(
        year = collect(dyears),
        gmst = [_f(gmst[i]) for i in idx],
        gmst_2022 = fill(_f(gmst2022), length(idx)),
        delta_gmst = [_f(dgmst[i]) for i in idx],
        delta_gmst_squared = [_f(dgmst2[i]) for i in idx],
    )
    path = joinpath(output_dir, "forest_temperature_diagnostic.csv")
    save(path, df)
    return path
end

"""
    write_clipping_diagnostic(built; years=2020:2300, output_dir=...) -> path

Summarises physical-constraint clipping across the run.
"""
function write_clipping_diagnostic(built; years = 2020:2300,
                                   output_dir::String = ForestConfig.OUTPUT_DIRECTORY)
    mkpath(output_dir)
    m = _model(built); inp = built.inputs
    idx = _year_indices(m, years)
    lower = m[:ForestAreaResponse, :forest_area_clipped_lower]
    upper = m[:ForestAreaResponse, :forest_area_clipped_upper]
    extrap = m[:ForestAreaResponse, :extrapolation_flag]
    nir = length(inp.impact_region_ids)

    per_ir_lower = zeros(Int, nir); per_ir_upper = zeros(Int, nir); per_ir_extrap = zeros(Int, nir)
    for ti in idx, r in 1:nir
        per_ir_lower[r]  += _f(lower[ti, r]) == 1 ? 1 : 0
        per_ir_upper[r]  += _f(upper[ti, r]) == 1 ? 1 : 0
        per_ir_extrap[r] += _f(extrap[ti, r]) == 1 ? 1 : 0
    end
    df = DataFrame(
        impact_region_id = inp.impact_region_ids,
        give_country_id  = [inp.give_country_ids[inp.impact_region_country_index[r]] for r in 1:nir],
        n_years_clipped_lower = per_ir_lower,
        n_years_clipped_upper = per_ir_upper,
        n_years_extrapolated  = per_ir_extrap,
    )
    sort!(df, :n_years_clipped_lower, rev = true)
    path = joinpath(output_dir, "forest_clipping_diagnostic.csv")
    save(path, df)
    return path
end

"""
    write_join_diagnostics(built; output_dir=...) -> path
"""
function write_join_diagnostics(built; output_dir::String = ForestConfig.OUTPUT_DIRECTORY)
    mkpath(output_dir)
    path = joinpath(output_dir, "forest_join_diagnostics.csv")
    save(path, built.inputs.diagnostics)
    return path
end
