

using Mimi, CSV, DataFrames

@defcomp forest_emulator begin
    regions   = Index()
    countries = Index()

    dT             = Parameter(index=[time])            # °C vs 2022
    region_country = Parameter{Int}(index=[regions])    # índice de país
    A_Mha  = Parameter(index=[regions])
    F0     = Parameter(index=[regions])                 # Mha
    b1     = Parameter(index=[regions])
    b2     = Parameter(index=[regions])
    C0     = Parameter(index=[regions])                 # t C
    c1     = Parameter(index=[regions])
    c2     = Parameter(index=[regions])

    forest_area = Variable(index=[time, regions])       # Mha
    cveg        = Variable(index=[time, regions])       # t C
    country_forest_Mha = Variable(index=[time, countries])
    country_co2_flow   = Variable(index=[time, countries])  # t CO2/año

    function run_timestep(p, v, d, t)
        for c in d.countries
            v.country_forest_Mha[t, c] = 0.0
            v.country_co2_flow[t, c]   = 0.0
        end
        x = p.dT[t]
        for r in d.regions
            F = min(p.A_Mha[r], max(0.0, p.F0[r] + p.b1[r]*x + p.b2[r]*x^2))
            C = max(0.0, p.C0[r] + p.c1[r]*x + p.c2[r]*x^2)
            v.forest_area[t, r] = F
            v.cveg[t, r] = C
            flow = is_first(t) ? 0.0 : (C - v.cveg[t-1, r]) * 44 / 12
            c = p.region_country[r]
            v.country_forest_Mha[t, c] += F
            v.country_co2_flow[t, c]   += flow
        end
    end
end

function build_forest_model(regions_csv, fair_csv; ssp = "SSP2")
    reg  = CSV.read(regions_csv, DataFrame)
    fair = CSV.read(fair_csv, DataFrame)
    fair = sort(filter(r -> r.ssp == ssp && 2022 <= r.year <= 2100, fair), :year)

    isos = sort(unique(reg.iso))
    idx  = Dict(i => k for (k, i) in enumerate(isos))

    m = Model()
    set_dimension!(m, :time, 2022:2100)
    set_dimension!(m, :regions, nrow(reg))
    set_dimension!(m, :countries, Symbol.(isos))
    add_comp!(m, forest_emulator)

    update = (n, v) -> update_param!(m, :forest_emulator, n, v)
    set_param!(m, :forest_emulator, :dT, fair.delta_gmst_vs_2022_C)
    set_param!(m, :forest_emulator, :region_country, [idx[i] for i in reg.iso])
    set_param!(m, :forest_emulator, :A_Mha, reg.region_area_ha ./ 1e6)
    set_param!(m, :forest_emulator, :F0, reg.forest_Mha_2022)
    set_param!(m, :forest_emulator, :b1, reg.forest_b1_Mha_per_C)
    set_param!(m, :forest_emulator, :b2, reg.forest_b2_Mha_per_C2)
    set_param!(m, :forest_emulator, :C0, reg.cveg_tC_2022)
    set_param!(m, :forest_emulator, :c1, reg.cveg_b1_tC_per_C)
    set_param!(m, :forest_emulator, :c2, reg.cveg_b2_tC_per_C2)
    return m
end

# --- Validación contra la app ---
m = build_forest_model("C:/Users/USER/export_forests/regions_jules.csv",
                       "C:/Users/USER/MimiGIVE_forests.jl/data/Forest/fair_gmst_2022_2100.csv")
run(m)

mine = getdataframe(m, :forest_emulator, :country_co2_flow)
rename!(mine, :time => :year, :countries => :iso3, :country_co2_flow => :mine)
mine.iso3 = string.(mine.iso3)

app = CSV.read("C:/Users/USER/export_forests/validation_carbon_jules.csv", DataFrame)
app = filter(r -> r.ssp == "SSP2" && r.year >= 2023, app)

cmp = innerjoin(mine, app, on = [:year, :iso3])
cmp.diff = cmp.mine .- cmp.carbon_flux_tCO2
println("Filas comparadas: ", nrow(cmp), " de ", nrow(app))
println("Máx. diferencia absoluta (tCO2): ", maximum(abs, skipmissing(cmp.diff)))