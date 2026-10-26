# Check before adding the CES utility channel: how large is E (non-market
# forest services: rec, hab, wat) relative to consumption C?
#
#   julia --project=. scripts/check_utility_share.jl
#
# 1) Share E0/C0 in 2022, for Mexico and for the world. GIVE's compute_scc
#    (default, no equity weights) discounts with global consumption per capita,
#    so the global share is the one that matters for the discount factor.
# 2) Discount correction. With U = [(1-a)C^θ + aE^θ]^(1/θ) and W = U^(1-η)/(1-η),
#        dW/dC = (1-a) C^(θ-1) U^(1-θ-η),
#    while GIVE uses C^(-η). The ratio relative to 2022 is
#        R_t = (U_t/C_t)^(1-θ-η) / (U_0/C_0)^(1-θ-η),
#    with a/(1-a) = (E0/C0)^(1-θ). R_t close to 1 means GIVE's discount is fine.
#    E_t is taken at constant 2022 prices (mv × forest area), as the CES channel
#    would use it.
#
# Units: country_value and mv × Mha are 1e6 USD2020/yr; gdp is billion USD2005/yr.
# GIVE consumption = GDP (net of damages), no savings rate.

using MimiGIVE, Mimi

m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                       forest_emulator = :jules)
run(m)

yrs   = collect(Mimi.time_labels(m))
i22   = findfirst(==(2022), yrs)
i2100 = findfirst(==(2100), yrs)
cn    = string.(Mimi.dim_keys(m, :country))
imex  = findfirst(==("MEX"), cn)
defl  = m[:ForestEmulator, :usd2020_to_2005]
nm    = [1, 2, 4]                                  # rec, hab, wat

gdp = m[:Socioeconomic, :gdp]                      # [time, country], bn USD2005
V   = m[:ForestEmulator, :country_value]           # [time, country, service], 1e6 USD2020

# E at constant 2022 prices, Mexico (all forest regions map to their country)
mv  = m[:ForestEmulator, :mv]                      # [region, service], USD2020/ha
F   = m[:ForestEmulator, :ir_forest_Mha]           # [time, region], Mha
rc  = m[:ForestEmulator, :region_country]
inm = rc .== imex
Econst(i) = sum(skipmissing(F[i, inm] .* sum(mv[inm, nm], dims = 2)[:])) / 1e3 * defl  # bn USD2005

E0_mex  = sum(skipmissing(V[i22, imex, nm])) / 1e3 * defl
E0_all  = sum(skipmissing(V[i22, :, nm])) / 1e3 * defl
C0_mex  = gdp[i22, imex]
C0_glob = sum(gdp[i22, :])

println("2022, billion USD2005/yr")
println("  E0 Mexico (rec+hab+wat)  = ", round(E0_mex, digits = 1))
println("  E0 check, constant price = ", round(Econst(i22), digits = 1))
println("  C0 Mexico                = ", round(C0_mex, digits = 1))
println("  C0 world                 = ", round(C0_glob, digits = 1))
println("  E0/C0 Mexico             = ", round(100 * E0_mex / C0_mex, digits = 2), " %")
println("  E0/C0 world              = ", round(100 * E0_all / C0_glob, digits = 4), " %")
println("  E/C Mexico 2100 (const.) = ", round(100 * Econst(i2100) / gdp[i2100, imex], digits = 2), " %")

# Discount correction R_t in 2100, Mexico (country-level, equity-weight case)
# and world (GIVE default)
η = 1.45                                           # compute_scc default
function R(E0, C0, Et, Ct, θ)
    k  = (E0 / C0)^(1 - θ)                         # a/(1-a)
    a  = k / (1 + k)
    UC(E, C) = ((1 - a) + a * (E / C)^θ)^(1 / θ)  # U/C
    (UC(Et, Ct) / UC(E0, C0))^(1 - θ - η)
end
C_mex_2100  = gdp[i2100, imex]
C_glob_2100 = sum(gdp[i2100, :])
E_mex_2100  = Econst(i2100)
println("\nDiscount correction R_2100 (1 = no correction needed), η = $η")
println(rpad("θ", 8), rpad("Mexico", 12), "World")
for θ in (1.0, 0.41, -1.02)
    rm = R(E0_mex, C0_mex, E_mex_2100, C_mex_2100, θ)
    rw = R(E0_all, C0_glob, E_mex_2100 * E0_all / E0_mex, C_glob_2100, θ)
    println(rpad(θ, 8), rpad(round(rm, digits = 4), 12), round(rw, digits = 6))
end
