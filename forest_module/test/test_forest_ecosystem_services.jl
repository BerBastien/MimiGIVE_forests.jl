# Unit tests for Component C: ForestEcosystemServices
# (placeholder valuation, reference vs climate, damage, aggregation)

using Test, Mimi

if !isdefined(@__MODULE__, :ForestEcosystemServices)
    include(joinpath(@__DIR__, "..", "src", "forest_ecosystem_services.jl"))
end

@testset "ecosystem_service_value kernel" begin
    @test ecosystem_service_value(2.0, 3.0, 4.0; scale = 1.0, spec_code = 1) == 24.0
    @test ecosystem_service_value(2.0, 3.0, 4.0; scale = 0.5, spec_code = 1) == 12.0
    @test_throws ErrorException ecosystem_service_value(1.0, 1.0, 1.0; spec_code = 99)
end

@testset "ForestEcosystemServices aggregation" begin
    years = 2020:2022
    m = Model()
    set_dimension!(m, :time, years)
    set_dimension!(m, :country, ["X", "Y"])
    set_dimension!(m, :impact_region, ["A", "B", "C"])   # A,B in X ; C in Y
    add_comp!(m, ForestEcosystemServices; first = 2020)

    baseline = [10.0, 5.0, 8.0]
    # projected: identical to baseline in 2020, forest LOSS by 2022
    proj = [10.0 5.0 8.0;    # 2020 == baseline
            9.0  5.0 8.0;    # 2021 A loses
            8.0  4.0 7.0]    # 2022 A,B,C lose
    pop = [1.0 2.0 3.0; 1.0 2.0 3.0; 1.0 2.0 3.0]
    gdp = [4.0 5.0 6.0; 4.0 5.0 6.0; 4.0 5.0 6.0]

    update_param!(m, :ForestEcosystemServices, :baseline_forest_area, baseline)
    update_param!(m, :ForestEcosystemServices, :projected_forest_area, proj)
    update_param!(m, :ForestEcosystemServices, :population_ir, pop)
    update_param!(m, :ForestEcosystemServices, :gdp_ir, gdp)
    update_param!(m, :ForestEcosystemServices, :impact_region_country_index, [1, 1, 2])
    update_param!(m, :ForestEcosystemServices, :es_value_scale, 1.0)
    update_param!(m, :ForestEcosystemServices, :es_value_spec_code, 1)
    run(m)

    ref_ir  = m[:ForestEcosystemServices, :es_value_reference_ir]
    dmg_ir  = m[:ForestEcosystemServices, :es_damage_ir]
    dmg_c   = m[:ForestEcosystemServices, :es_damage_country]
    dmg_g   = m[:ForestEcosystemServices, :es_damage_global]
    ref_c   = m[:ForestEcosystemServices, :es_value_reference_country]
    ref_g   = m[:ForestEcosystemServices, :es_value_reference_global]

    # reference uses BASELINE area (independent of climate), so equal across time
    @test ref_ir[1, 1] ≈ baseline[1] * gdp[1, 1] * pop[1, 1]   # 10*4*1 = 40
    @test ref_ir[3, 1] ≈ ref_ir[1, 1]

    # in 2020, projected == baseline => zero damage everywhere
    @test all(dmg_ir[1, :] .== 0.0)
    @test all(dmg_c[1, :] .== 0.0)
    @test dmg_g[1] == 0.0

    # in 2022, A lost 2 units of area => damage_A = (10-8)*gdp_A*pop_A = 2*4*1 = 8
    @test dmg_ir[3, 1] ≈ 8.0
    # B lost 1 => 1*5*2 = 10 ; C lost 1 => 1*6*3 = 18
    @test dmg_ir[3, 2] ≈ 10.0
    @test dmg_ir[3, 3] ≈ 18.0

    # country aggregation: X = A+B ; Y = C
    @test dmg_c[3, 1] ≈ dmg_ir[3, 1] + dmg_ir[3, 2]     # 18
    @test dmg_c[3, 2] ≈ dmg_ir[3, 3]                    # 18
    # global == sum of countries == sum of impact regions
    @test dmg_g[3] ≈ dmg_c[3, 1] + dmg_c[3, 2]
    @test dmg_g[3] ≈ sum(dmg_ir[3, :])

    # reference country/global aggregation consistency
    @test ref_c[2, 1] ≈ ref_ir[2, 1] + ref_ir[2, 2]
    @test ref_g[2] ≈ sum(ref_ir[2, :])
end
