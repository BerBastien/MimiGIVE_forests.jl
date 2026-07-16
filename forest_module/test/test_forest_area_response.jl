# Unit tests for Component B: ForestAreaResponse
# These build a TINY standalone Mimi model (no MimiGIVE stack needed).

using Test, Mimi

# Include just the component (avoids pulling in the whole MimiGIVE stack).
if !isdefined(@__MODULE__, :ForestAreaResponse)
    include(joinpath(@__DIR__, "..", "src", "forest_area_response.jl"))
end

@testset "ForestAreaResponse" begin

    years = 2020:2025
    i2022 = findfirst(==(2022), collect(years))   # == 3
    i2024 = findfirst(==(2024), collect(years))   # == 5

    m = Model()
    set_dimension!(m, :time, years)
    set_dimension!(m, :impact_region, ["A", "B"])
    add_comp!(m, ForestAreaResponse; first = 2020)

    # GMST: value at 2022 is 1.2; at 2024 is 1.7  => dT(2024) = 0.5
    update_param!(m, :ForestAreaResponse, :gmst, [1.0, 1.1, 1.2, 1.4, 1.7, 2.0])
    update_param!(m, :ForestAreaResponse, :beta1, [2.0, 0.0])    # A: b1=2 ; B: no response
    update_param!(m, :ForestAreaResponse, :beta2, [1.0, 0.0])    # A: b2=1
    update_param!(m, :ForestAreaResponse, :intercept, [0.0, 0.0])
    update_param!(m, :ForestAreaResponse, :baseline_forest_area, [100.0, 50.0])
    update_param!(m, :ForestAreaResponse, :total_impact_region_area, [0.0, 0.0])
    update_param!(m, :ForestAreaResponse, :delta_gmst_fit_min, [-Inf, -Inf])
    update_param!(m, :ForestAreaResponse, :delta_gmst_fit_max, [Inf, Inf])
    update_param!(m, :ForestAreaResponse, :baseline_year, 2022)
    update_param!(m, :ForestAreaResponse, :change_units_code, 1)   # :percent
    update_param!(m, :ForestAreaResponse, :apply_lower_clip, true)
    update_param!(m, :ForestAreaResponse, :apply_upper_clip, false)

    run(m)

    dgmst  = m[:ForestAreaResponse, :delta_gmst]
    dgmst2 = m[:ForestAreaResponse, :delta_gmst_squared]
    fchg   = m[:ForestAreaResponse, :forest_change]
    proj   = m[:ForestAreaResponse, :projected_forest_area]
    gbase  = m[:ForestAreaResponse, :gmst_baseline]

    @testset "temperature baseline" begin
        @test gbase == 1.2
        @test dgmst[i2022] == 0.0
        @test dgmst2[i2022] == 0.0
        # no intercept => forest change 0 and projected == baseline in 2022
        @test fchg[i2022, 1] == 0.0
        @test proj[i2022, 1] == 100.0
        @test proj[i2022, 2] == 50.0
    end

    @testset "coefficient formula (percent units)" begin
        @test dgmst[i2024] ≈ 0.5
        @test dgmst2[i2024] ≈ 0.25
        # change = 2*0.5 + 1*0.5^2 = 1.25
        @test fchg[i2024, 1] ≈ 1.25
        # area = 100 * (1 + 1.25/100) = 101.25
        @test proj[i2024, 1] ≈ 101.25
        # region B has zero coefficients -> unchanged
        @test all(proj[:, 2] .== 50.0)
    end
end

@testset "ForestAreaResponse: lower clip & zero-response" begin
    years = 2020:2024
    m = Model()
    set_dimension!(m, :time, years)
    set_dimension!(m, :impact_region, ["A"])
    add_comp!(m, ForestAreaResponse; first = 2020)

    # Big negative response so raw area goes below zero and must be clipped.
    update_param!(m, :ForestAreaResponse, :gmst, [0.0, 0.0, 0.0, 5.0, 10.0])  # dT huge later
    update_param!(m, :ForestAreaResponse, :beta1, [-50.0])
    update_param!(m, :ForestAreaResponse, :beta2, [0.0])
    update_param!(m, :ForestAreaResponse, :intercept, [0.0])
    update_param!(m, :ForestAreaResponse, :baseline_forest_area, [10.0])
    update_param!(m, :ForestAreaResponse, :total_impact_region_area, [0.0])
    update_param!(m, :ForestAreaResponse, :delta_gmst_fit_min, [-Inf])
    update_param!(m, :ForestAreaResponse, :delta_gmst_fit_max, [Inf])
    update_param!(m, :ForestAreaResponse, :baseline_year, 2022)
    update_param!(m, :ForestAreaResponse, :change_units_code, 1)
    update_param!(m, :ForestAreaResponse, :apply_lower_clip, true)
    update_param!(m, :ForestAreaResponse, :apply_upper_clip, false)
    run(m)

    proj = m[:ForestAreaResponse, :projected_forest_area]
    raw  = m[:ForestAreaResponse, :projected_forest_area_raw]
    clip = m[:ForestAreaResponse, :forest_area_clipped_lower]
    @test all(proj .>= 0.0)                 # never negative
    @test any(raw .< 0.0)                    # raw did go negative
    @test any(x -> x == 1, clip)             # clipping was flagged

    # zero-response check: rebuild with all coefficients zero -> area == baseline
    m2 = Model()
    set_dimension!(m2, :time, years)
    set_dimension!(m2, :impact_region, ["A"])
    add_comp!(m2, ForestAreaResponse; first = 2020)
    update_param!(m2, :ForestAreaResponse, :gmst, [0.0, 0.0, 0.0, 5.0, 10.0])
    update_param!(m2, :ForestAreaResponse, :beta1, [0.0])
    update_param!(m2, :ForestAreaResponse, :beta2, [0.0])
    update_param!(m2, :ForestAreaResponse, :intercept, [0.0])
    update_param!(m2, :ForestAreaResponse, :baseline_forest_area, [10.0])
    update_param!(m2, :ForestAreaResponse, :total_impact_region_area, [0.0])
    update_param!(m2, :ForestAreaResponse, :delta_gmst_fit_min, [-Inf])
    update_param!(m2, :ForestAreaResponse, :delta_gmst_fit_max, [Inf])
    update_param!(m2, :ForestAreaResponse, :baseline_year, 2022)
    update_param!(m2, :ForestAreaResponse, :change_units_code, 1)
    update_param!(m2, :ForestAreaResponse, :apply_lower_clip, true)
    update_param!(m2, :ForestAreaResponse, :apply_upper_clip, false)
    run(m2)
    @test all(m2[:ForestAreaResponse, :projected_forest_area] .== 10.0)
end
