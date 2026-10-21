@defcomp ForestEmulator begin
    country        = Index()
    forest_regions = Index()
    services       = Index()     # 1=rec, 2=hab, 3=nwfp, 4=wat

    # --- Parámetros ---
    temperature    = Parameter(index=[time], unit="degC")
    population     = Parameter(index=[time, country])
    pc_gdp         = Parameter(index=[time, country])
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
    mv    = Parameter(index=[forest_regions, services])   # USD2020/ha/año
    e_gdp = Parameter(index=[services])
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
    z_epop_rec = Parameter(default = 0.0)
    z_epop_hab = Parameter(default = 0.0)
    z_epop_nwfp = Parameter(default = 0.0)
    z_epop_wat = Parameter(default = 0.0)
    # Forest damages to consumption (market part of the services)
    # damage_weight[s] = share of service s counted as a market loss (0 = none, 1 = all)
    damage_weight  = Parameter(index=[services])
    usd2020_to_2005 = Parameter(default = 0.773)   # US GDP deflator 2005/2020 (to verify)
    pop_share = Parameter(index=[forest_regions])
    gdp_share = Parameter(index=[forest_regions])

    # --- Variables ---
    T_ref     = Variable()
    cveg_prev = Variable(index=[forest_regions])
    pop0      = Variable(index=[country])
    gdp0      = Variable(index=[country])
    global_forest_Mha  = Variable(index=[time])
    global_co2_flow    = Variable(index=[time])                  # t CO2/año
    country_forest_Mha = Variable(index=[time, country])
    country_co2_flow   = Variable(index=[time, country])         # t CO2/año
    country_value      = Variable(index=[time, country, services])  # USD/ha × Mha = 1e6 USD/año
    country_value_nocl = Variable(index=[time, country, services])  # bosque fijo, socio evolutiva
    ir_forest_Mha = Variable(index=[time, forest_regions])
    ir_value      = Variable(index=[time, forest_regions, services])
    # Climate-driven loss of market value: (value with forest fixed at 2022) − (value with climate)
    country_damage = Variable(index=[time, country], unit="US\$2005/yr")
    ir_pop        = Variable(index=[time, forest_regions])
    ir_gdppc      = Variable(index=[time, forest_regions])

    function run_timestep(p, v, d, t)
        if is_first(t)
            v.T_ref = p.temperature[t]
            for c in d.country
                v.pop0[c] = p.population[t, c]
                v.gdp0[c] = p.pc_gdp[t, c]
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

        # Elasticities for this run = estimate + shock × standard error
        z_gdp = (p.z_egdp_rec, p.z_egdp_hab, p.z_egdp_nwfp, p.z_egdp_wat)
        z_pop = (p.z_epop_rec, p.z_epop_hab, p.z_epop_nwfp, p.z_epop_wat)
        e_gdp = [p.e_gdp[s] + z_gdp[s] * p.e_gdp_se[s] for s in 1:4]
        e_pop = [p.e_pop[s] + z_pop[s] * p.e_pop_se[s] for s in 1:4]

        fac = zeros(length(d.country), 4)
        for c in d.country, s in 1:4
            pr = p.population[t, c] / v.pop0[c]
            gr = p.pc_gdp[t, c] / v.gdp0[c]
            fac[c, s] = (isfinite(pr) && isfinite(gr)) ? gr^e_gdp[s] * pr^e_pop[s] : 0.0
        end

        gF = 0.0
        gQ = 0.0
        x = p.temperature[t] - v.T_ref
        for r in d.forest_regions
            # Coefficient for this run = estimate + shock × standard error
            b1 = p.b1[r] + p.z_b1 * p.b1_se[r]
            b2 = p.b2[r] + p.z_b2 * p.b2_se[r]
            c1 = p.c1[r] + p.z_c1 * p.c1_se[r]
            c2 = p.c2[r] + p.z_c2 * p.c2_se[r]
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
                    p.pc_gdp[t, c] * p.gdp_share[r] / p.pop_share[r] : 0.0

                for s in 1:4
                    val = p.mv[r, s] * fac[c, s] * F
                    v.country_value[t, c, s]      += val
                    v.country_value_nocl[t, c, s] += p.mv[r, s] * fac[c, s] * p.F0[r]
                    v.ir_value[t, r, s] = val
                end
            end
        end
        v.global_forest_Mha[t] = gF
        v.global_co2_flow[t]   = gQ

        # Market damage by country: 1e6 USD2020/yr → USD2005/yr
        for c in d.country, s in 1:4
            loss = v.country_value_nocl[t, c, s] - v.country_value[t, c, s]
            v.country_damage[t, c] += p.damage_weight[s] * loss * 1e6 * p.usd2020_to_2005
        end
    end
end