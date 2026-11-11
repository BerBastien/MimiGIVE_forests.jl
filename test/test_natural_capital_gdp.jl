@testitem "TestNaturalCapitalGDP" setup=[DataDepsSetup] begin

    import MimiGIVE: get_model, set_natural_capital!

    ##------------------------------------------------------------------------------
    ## Natural capital GDP: reference, switch, GDP chain and 2022-2300 emulator
    ##------------------------------------------------------------------------------

    args = (socioeconomics_source = :SSP, SSP_scenario = "SSP245", forest_emulator = :jules)

    # Switch off (default): GDP is the Socioeconomic GDP even though it goes through
    # NaturalCapitalGDP, the component every GDP consumer reads from.
    m_off = get_model(; args...)
    run(m_off)

    # Switch on: GDP is adjusted by the natural capital ratio with elasticity 0.1.
    m_on = get_model(; args..., natural_capital = true, natural_capital_elasticity = 0.1)
    run(m_on)

    num(x) = Float64.(coalesce.(x, 0.0))

    years = collect(Mimi.time_labels(m_on))
    ti    = findall(y -> 2020 <= y <= 2300, years)   # components start in 2020
    t2020 = findfirst(==(2020), years)
    t2022 = findfirst(==(2022), years)
    t2100 = findfirst(==(2100), years)
    t2300 = findfirst(==(2300), years)

    gdp_ssp = num(m_off[:Socioeconomic, :gdp])
    gdp_off = num(m_off[:NaturalCapitalGDP, :gdp])
    gdp_on  = num(m_on[:NaturalCapitalGDP, :gdp])

    ##------------------------------------------------------------------------------
    ## 1. The switch off reproduces the Socioeconomic GDP exactly
    ##------------------------------------------------------------------------------

    @test gdp_off[ti, :] ≈ gdp_ssp[ti, :] rtol = 1e-12
    @test m_off[:NaturalCapitalGDP, :use_nc] == false

    ##------------------------------------------------------------------------------
    ## 2. Reference natural capital is the no-climate counterfactual
    ##------------------------------------------------------------------------------

    @test m_on[:NaturalCapitalGDP, :use_nc] == true
    S_hab = num(m_on[:NaturalCapitalGDP, :S_hab])
    S_ref = num(m_on[:NaturalCapitalGDP, :S_ref])
    hab_nocl = num(m_on[:ForestEmulator, :global_hab_value_nocl])

    @test S_ref[t2022:end] ≈ hab_nocl[t2022:end] rtol = 1e-12
    @test S_hab[t2022] ≈ S_ref[t2022] rtol = 1e-12     # ratio is 1 in the reference year
    @test all(S_ref[t2022:end] .> 0)
    # Before 2022 (first year of the emulator) the connection backup keeps both at zero,
    # so GDP is unadjusted in 2020-2021.
    @test S_hab[t2020] == 0.0
    @test gdp_on[t2020, :] ≈ gdp_ssp[t2020, :] rtol = 1e-12

    ##------------------------------------------------------------------------------
    ## 3. The adjustment follows gdp_ssp * (S_hab / S_ref) ^ b
    ##------------------------------------------------------------------------------

    b = num(m_on[:NaturalCapitalGDP, :b])
    @test all(b .== 0.1)
    for t in (2050, 2100, 2300)
        i = findfirst(==(t), years)
        ratio = S_hab[i] / S_ref[i]
        @test gdp_on[i, :] ≈ gdp_ssp[i, :] .* ratio .^ b rtol = 1e-12
    end

    # The adjustment is not a no-op: at some point it moves GDP
    @test !isapprox(gdp_on[ti, :], gdp_ssp[ti, :]; rtol = 1e-6)

    ##------------------------------------------------------------------------------
    ## 4. Every GDP consumer reads the adjusted GDP
    ##------------------------------------------------------------------------------

    pop = num(m_on[:Socioeconomic, :population])
    @test num(m_on[:PerCapitaGDP, :gdp])[ti, :] ≈ gdp_on[ti, :] rtol = 1e-12
    @test num(m_on[:PerCapitaGDP, :pc_gdp])[t2100, :] ≈ gdp_on[t2100, :] ./ pop[t2100, :] .* 1e3 rtol = 1e-12

    # VSL is downstream of PerCapitaGDP, so it reacts to the adjusted GDP
    @test !isapprox(num(m_on[:VSL, :vsl])[t2100, :],
                    num(m_off[:VSL, :vsl])[t2100, :]; rtol = 1e-6)

    # Every other GDP consumer takes the adjusted GDP
    for (comp, par) in ((:country_netconsumption, :gdp), (:global_netconsumption, :gdp),
                        (:DamageAggregator, :gdp), (:energy_damages, :gdp),
                        (:dice2016R2_damage, :gdp), (:hs_damage, :gdp),
                        (:Agriculture_aggregator_gdp, :input),
                        (:AgricultureDamagesDisaggregator, :gdp_country))
        @test num(m_on[comp, par])[ti, :] ≈ gdp_on[ti, :] rtol = 1e-12
    end

    ##------------------------------------------------------------------------------
    ## 5. The switch can be flipped on an already built model
    ##------------------------------------------------------------------------------

    set_natural_capital!(m_on; use = false)
    run(m_on)
    @test num(m_on[:NaturalCapitalGDP, :gdp])[ti, :] ≈ gdp_ssp[ti, :] rtol = 1e-12

    set_natural_capital!(m_on; use = true)
    run(m_on)
    @test num(m_on[:NaturalCapitalGDP, :gdp])[ti, :] ≈ gdp_on[ti, :] rtol = 1e-12

    # An elasticity of zero leaves GDP unchanged even with the switch on
    set_natural_capital!(m_on; elasticity = 0.0)
    run(m_on)
    @test num(m_on[:NaturalCapitalGDP, :gdp])[ti, :] ≈ gdp_ssp[ti, :] rtol = 1e-12

    ##------------------------------------------------------------------------------
    ## 6. The forest emulator runs to 2300 (2022-2100 regression + extrapolation)
    ##------------------------------------------------------------------------------

    forest_Mha = num(m_on[:ForestEmulator, :global_forest_Mha])
    @test forest_Mha[t2300] > 0
    @test num(m_on[:ForestEmulator, :global_hab_value])[t2300] > 0
    @test num(m_on[:ForestEmulator, :global_co2_flow])[t2300] != 0
    # Forest area stays within its physical bounds after the calibrated period
    @test 0 < forest_Mha[t2100] < sum(num(m_on[:ForestEmulator, :A_Mha]))

    ##------------------------------------------------------------------------------
    ## 7. Natural capital needs the forest emulator
    ##------------------------------------------------------------------------------

    @test_throws ErrorException get_model(; socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                                          natural_capital = true)

end # module
