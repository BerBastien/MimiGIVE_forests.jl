@defcomp ForestEmulator begin
    country        = Index()
    forest_regions = Index()

    temperature    = Parameter(index=[time], unit="degC")
    region_country = Parameter{Int}(index=[forest_regions])
    A_Mha = Parameter(index=[forest_regions])
    F0    = Parameter(index=[forest_regions])
    b1    = Parameter(index=[forest_regions])
    b2    = Parameter(index=[forest_regions])
    C0    = Parameter(index=[forest_regions])
    c1    = Parameter(index=[forest_regions])
    c2    = Parameter(index=[forest_regions])

    T_ref     = Variable()
    cveg_prev = Variable(index=[forest_regions])
    global_forest_Mha = Variable(index=[time])
    global_co2_flow   = Variable(index=[time])               # t CO2/año, todas las regiones
    country_forest_Mha = Variable(index=[time, country])
    country_co2_flow   = Variable(index=[time, country])     # t CO2/año

    function run_timestep(p, v, d, t)
        if is_first(t)
            v.T_ref = p.temperature[t]
        end
        v.country_forest_Mha[t, :] .= 0.0
        v.country_co2_flow[t, :]   .= 0.0
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
            end
        end
        v.global_forest_Mha[t] = gF
        v.global_co2_flow[t]   = gQ
    end
end