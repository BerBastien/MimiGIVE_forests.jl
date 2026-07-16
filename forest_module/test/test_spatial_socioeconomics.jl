# Unit tests for Component A: SpatialSocioeconomics (downscaling conservation)

using Test, Mimi

if !isdefined(@__MODULE__, :SpatialSocioeconomics)
    include(joinpath(@__DIR__, "..", "src", "forest_spatial_socioeconomics.jl"))
end

@testset "SpatialSocioeconomics downscaling" begin
    years = 2020:2022
    m = Model()
    set_dimension!(m, :time, years)
    set_dimension!(m, :country, ["X", "Y"])
    set_dimension!(m, :impact_region, ["A", "B", "C"])   # A,B in X ; C in Y
    add_comp!(m, SpatialSocioeconomics; first = 2020)

    # time x country
    popc = [10.0 20.0; 11.0 21.0; 12.0 22.0]
    gdpc = [100.0 200.0; 110.0 210.0; 120.0 220.0]
    update_param!(m, :SpatialSocioeconomics, :population_country, popc)
    update_param!(m, :SpatialSocioeconomics, :gdp_country, gdpc)
    update_param!(m, :SpatialSocioeconomics, :population_share, [0.6, 0.4, 1.0])
    update_param!(m, :SpatialSocioeconomics, :gdp_share, [0.7, 0.3, 1.0])
    update_param!(m, :SpatialSocioeconomics, :impact_region_country_index, [1, 1, 2])
    run(m)

    pop_ir = m[:SpatialSocioeconomics, :population_ir]
    gdp_ir = m[:SpatialSocioeconomics, :gdp_ir]
    gpc_ir = m[:SpatialSocioeconomics, :gdppc_ir]

    for (ti, _) in enumerate(years)
        # conservation: sum over impact regions within a country == country total
        @test pop_ir[ti, 1] + pop_ir[ti, 2] ≈ popc[ti, 1]   # X
        @test pop_ir[ti, 3] ≈ popc[ti, 2]                   # Y
        @test gdp_ir[ti, 1] + gdp_ir[ti, 2] ≈ gdpc[ti, 1]
        @test gdp_ir[ti, 3] ≈ gdpc[ti, 2]
        # gdppc = gdp/pop*1e3 (billion$/million -> $/person)
        @test gpc_ir[ti, 1] ≈ gdp_ir[ti, 1] / pop_ir[ti, 1] * 1e3
    end

    # exact values, first timestep
    @test pop_ir[1, 1] ≈ 6.0    # 10 * 0.6
    @test gdp_ir[1, 2] ≈ 30.0   # 100 * 0.3
    @test gdp_ir[1, 3] ≈ 200.0  # 200 * 1.0
end
