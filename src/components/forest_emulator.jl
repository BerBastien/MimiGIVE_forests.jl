# Forest and ecosystem services emulator.
#
# Maps global mean temperature to forest area and ecosystem service values
# at the impact-region (IR) level. It is optional:
# when it is added, NaturalCapitalGDP uses its habitat value as the natural
# capital stock, and DamageAggregator can add its market damages to consumption.
#
# Equations. With x = T(t) - T(2022) and services s = rec, hab, nwfp, wat:
#
#   F(r,t)   = min(A(r), max(0, F0(r) + b1(r)*x + b2(r)*x^2))  [Mha]
#   C(r,t)   = max(0, C0(r) + c1(r)*x + c2(r)*x^2)             [tC]
#   mv(r,s,t)= mv(r,s,2022) * (GDP(c,t)/GDP(c,2022))^e_gdp(s)
#   V(r,s,t) = mv(r,s,t) * F(r,t)                             [1e6 USD2020/yr]
#
# The value per hectare changes according to real GDP growth only.
# Population does not affect ecosystem service valuation.
#
# Time span. The regional regressions are estimated on 2022-2100. The model
# carries the emulator to 2300 by evaluating the same quadratic after 2100
# (extrapolation), with the forest area capped at the region area A and the
# vegetation carbon floored at zero, so both stay in a physically plausible
# range. Set `last` to 2100 to run the calibrated period only.
#
# The emulator starts in 2022, so the first timestep is the calibration
# reference: x = 0 and GDP is normalized to its 2022 value.

@defcomp ForestEmulator begin
    country        = Index()
    forest_regions = Index()
    services       = Index()     # 1=rec, 2=hab, 3=nwfp, 4=wat

    # --- Parameters ---
    temperature    = Parameter(index=[time], unit="degC")
    population     = Parameter(index=[time, country])

    # Unadjusted real GDP: the emulator feeds NaturalCapitalGDP, so reading
    # the adjusted GDP here would create a cycle.
    gdp = Parameter(index=[time, country], unit="billion US\$2005/yr")

    region_country = Parameter{Int}(index=[forest_regions])
    A_Mha = Parameter(index=[forest_regions])
    F0    = Parameter(index=[forest_regions])
    b1    = Parameter(index=[forest_regions])
    b2    = Parameter(index=[forest_regions])
    C0    = Parameter(index=[forest_regions])
    c1    = Parameter(index=[forest_regions])
    c2    = Parameter(index=[forest_regions])

    # Regression uncertainty per region (standard errors from the CSV)
    b1_se = Parameter(index=[forest_regions])
    b2_se = Parameter(index=[forest_regions])
    c1_se = Parameter(index=[forest_regions])
    c2_se = Parameter(index=[forest_regions])

    # Standardized shocks ~ Normal(0, 1), sampled by the Monte Carlo.
    # At 0 (default) the model is identical to the deterministic run.
    z_b1 = Parameter(default = 0.0)
    z_b2 = Parameter(default = 0.0)
    z_c1 = Parameter(default = 0.0)
    z_c2 = Parameter(default = 0.0)

    mv = Parameter(index=[forest_regions, services])   # USD2020/ha/yr

    # GDP elasticity of ecosystem service value
    e_gdp = Parameter(index=[services])

    # Retained for compatibility with existing model parameterization.
    # Population elasticity is no longer used in ecosystem service valuation.
    e_pop = Parameter(index=[services])

    # Elasticity uncertainty: robust standard errors of the meta-regression
    # (CWON 2024 Non-Wood Forest Technical Report, Table 5.1)
    e_gdp_se = Parameter(index=[services])
    e_pop_se = Parameter(index=[services])

    # One standardized shock per service ~ Normal(0, 1), sampled by the Monte Carlo
    # (each service is a separate regression, so the shocks are independent)
    z_egdp_rec = Parameter(default = 0.0)
    z_egdp_hab = Parameter(default = 0.0)
    z_egdp_nwfp = Parameter(default = 0.0)
    z_egdp_wat = Parameter(default = 0.0)

    # Retained for compatibility; these shocks are no longer used.
    z_epop_rec = Parameter(default = 0.0)
    z_epop_hab = Parameter(default = 0.0)
    z_epop_nwfp = Parameter(default = 0.0)
    z_epop_wat = Parameter(default = 0.0)

    # Forest damages to consumption (market part of the services)
    # damage_weight[s] = share of service s counted as a market loss (0 = none, 1 = all)
    damage_weight  = Parameter(index=[services])
    usd2020_to_2005 = Parameter(default = 0.773)   # US GDP deflator 2005/2020

    pop_share = Parameter(index=[forest_regions])
    gdp_share = Parameter(index=[forest_regions])

    # --- Variables ---
    T_ref     = Variable()
    cveg_prev = Variable(index=[forest_regions])

    # Reference real GDP in 2022
    gdp0      = Variable(index=[country])

    global_forest_Mha  = Variable(index=[time])
    global_co2_flow    = Variable(index=[time])

    country_forest_Mha = Variable(index=[time, country])
    country_co2_flow   = Variable(index=[time, country])

    country_value      = Variable(index=[time, country, services])

    # Counterfactual: forest area fixed at 2022, while unit values evolve
    # according to real GDP growth.
    country_value_nocl = Variable(index=[time, country, services])

    ir_forest_Mha = Variable(index=[time, forest_regions])
    ir_value      = Variable(index=[time, forest_regions, services])

    # Climate-driven loss of market value:
    # (value with forest fixed at 2022) - (value with climate)
    country_damage = Variable(index=[time, country], unit="US\$2005/yr")

    ir_pop        = Variable(index=[time, forest_regions])
    ir_gdppc      = Variable(index=[time, forest_regions])

    # Global value of habitat services in the run and in its no-climate
    # counterfactual (forest fixed at 2022, GDP evolves).
    # Their ratio is the natural capital stock index used by NaturalCapitalGDP.
    global_hab_value      = Variable(index=[time])
    global_hab_value_nocl = Variable(index=[time])

    function run_timestep(p, v, d, t)

        pcg(t, c) = p.gdp[t, c] / p.population[t, c] * 1e3

        # The emulator starts in 2022, so the first timestep is the calibration
        # reference: x = T(t) - T(2022), and GDP is normalized to its 2022 value.
        if is_first(t)
            v.T_ref = p.temperature[t]

            for c in d.country
                v.gdp0[c] = p.gdp[t, c]
            end
        end

        v.country_forest_Mha[t, :]    .= 0.0
        v.country_co2_flow[t, :]      .= 0.0
        v.country_value[t, :, :]      .= 0.0
        v.country_value_nocl[t, :, :] .= 0.0
        v.country_damage[t, :]        .= 0.0

        v.ir_forest_Mha[t, :] .= 0.0
        v.ir_value[t, :, :]   .= 0.0
        v.ir_pop[t, :]        .= 0.0
        v.ir_gdppc[t, :]      .= 0.0

        # GDP elasticities for this run = estimate + shock * standard error.
        # Population elasticities are not included in ecosystem service values.
        z_gdp = (p.z_egdp_rec, p.z_egdp_hab,
                 p.z_egdp_nwfp, p.z_egdp_wat)

        e_gdp = [p.e_gdp[s] + z_gdp[s] * p.e_gdp_se[s]
                 for s in 1:4]

        # GDP growth factor relative to 2022.
        # The factor is calculated separately for each country and service.
        gdp_factor = zeros(length(d.country), 4)

        for c in d.country, s in 1:4
            gdp_ratio = p.gdp[t, c] / v.gdp0[c]

            gdp_factor[c, s] =
                (isfinite(gdp_ratio) && gdp_ratio > 0) ?
                gdp_ratio^e_gdp[s] : 0.0
        end

        gF = 0.0
        gQ = 0.0

        x = p.temperature[t] - v.T_ref

        for r in d.forest_regions

            # Coefficient for this run = estimate + shock * standard error
            b1 = p.b1[r] + p.z_b1 * p.b1_se[r]
            b2 = p.b2[r] + p.z_b2 * p.b2_se[r]
            c1 = p.c1[r] + p.z_c1 * p.c1_se[r]
            c2 = p.c2[r] + p.z_c2 * p.c2_se[r]

            # Quadratic in x. After 2100 this is an extrapolation of the
            # regression estimated on 2022-2100; A and the zero floor keep it
            # within a plausible range.
            F = min(p.A_Mha[r], max(0.0, p.F0[r] + b1*x + b2*x^2))
            C = max(0.0, p.C0[r] + c1*x + c2*x^2)

            flow = is_first(t) ? 0.0 : (C - v.cveg_prev[r]) * 44 / 12

            v.cveg_prev[r] = C

            gF += F
            gQ += flow

            c = p.region_country[r]

            if c > 0

                v.country_forest_Mha[t, c] += F
                v.country_co2_flow[t, c]   += flow

                v.ir_forest_Mha[t, r] = F
                v.ir_pop[t, r]   = p.population[t, c] * p.pop_share[r]
                v.ir_gdppc[t, r] = p.pop_share[r] > 0 ?
                    pcg(t, c) * p.gdp_share[r] / p.pop_share[r] : 0.0

                for s in 1:4

                    # Unit value updated according to real GDP growth
                    mv_t = p.mv[r, s] * gdp_factor[c, s]

                    # Ecosystem service value with projected forest area
                    val = mv_t * F

                    v.country_value[t, c, s] += val

                    # No-climate counterfactual:
                    # same GDP-adjusted unit value, but forest area fixed at F0.
                    v.country_value_nocl[t, c, s] += mv_t * p.F0[r]

                    v.ir_value[t, r, s] = val
                end
            end
        end

        v.global_forest_Mha[t] = gF
        v.global_co2_flow[t]   = gQ

        v.global_hab_value[t] =
            sum(v.country_value[t, :, 2])

        v.global_hab_value_nocl[t] =
            sum(v.country_value_nocl[t, :, 2])

        # Market damage by country: 1e6 USD2020/yr -> USD2005/yr
        for c in d.country, s in 1:4
            loss = v.country_value_nocl[t, c, s] -
                   v.country_value[t, c, s]

            v.country_damage[t, c] +=
                p.damage_weight[s] * loss * 1e6 * p.usd2020_to_2005
        end
    end
end
