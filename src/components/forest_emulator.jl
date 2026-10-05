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
    mv    = Parameter(index=[forest_regions, services])   # USD2020/ha/año
    e_gdp = Parameter(index=[services])
    e_pop = Parameter(index=[services])
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

        v.ir_forest_Mha[t, :] .= 0.0
        v.ir_value[t, :, :]   .= 0.0
        v.ir_pop[t, :]        .= 0.0
        v.ir_gdppc[t, :]      .= 0.0

        fac = zeros(length(d.country), 4)
        for c in d.country, s in 1:4
            pr = p.population[t, c] / v.pop0[c]
            gr = p.pc_gdp[t, c] / v.gdp0[c]
            fac[c, s] = (isfinite(pr) && isfinite(gr)) ? gr^p.e_gdp[s] * pr^p.e_pop[s] : 0.0
        end

        gF = 0.0
        gQ = 0.0
        x = p.temperature[t] - v.T_ref
        for r in d.forest_regions
            F = min(p.A_Mha[r], max(0.0, p.F0[r] + p.b1[r]*x + p.b2[r]*x^2))
            C = max(0.0, p.C0[r] + p.c1[r]*x + p.c2[r]*x^2)
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
    end
end