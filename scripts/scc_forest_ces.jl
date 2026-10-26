# SCC 2020 with non-market forest damages through the CES utility.
#
#   julia --project=. scripts/scc_forest_ces.jl
#
# 1) Control: CES with θ = 1 and the emulator prices (ces_income_price = true) must
#    reproduce the additive case (rec, hab, wat added to consumption damages): same
#    SCC (22.64) and same damages.
# 2) Table: nwfp in consumption + rec, hab, wat through the CES at constant 2022
#    prices, for several θ. θ = 1 − income elasticity of the price:
#      θ = 1      perfect substitutes, no relative price change
#      θ = 0.41   Drupp et al. (2025), elasticity 0.59
#      θ = 0      Cobb-Douglas
#      θ = −1.02  CWON water elasticity 2.02
# SSP245, JULES, deterministic. Damages are 0 after 2100 (emulator span).

using MimiGIVE, Mimi

base = (socioeconomics_source = :SSP, SSP_scenario = "SSP245", forest_emulator = :jules,
        forest_damages = true)
nonmkt = [1.0, 1.0, 0.0, 1.0]       # rec, hab, wat
nwfp   = [0.0, 0.0, 1.0, 0.0]

function total_damage_2100(m)
    yrs = collect(Mimi.time_labels(m))
    i = findfirst(==(2100), yrs)
    sum(m[:ForestEmulator, :country_damage][i, :]) / 1e9   # billion USD2005/yr
end

# --- 1) Control ---------------------------------------------------------------
m_add = MimiGIVE.get_model(; base..., forest_damage_weight = nonmkt)
m_ces = MimiGIVE.get_model(; base..., forest_damage_weight = zeros(4),
                           forest_utility_weight = nonmkt, forest_ces_theta = 1.0,
                           forest_ces_income_price = true)
run(m_add); run(m_ces)
dmax = maximum(abs.(m_add[:ForestEmulator, :country_damage] .- m_ces[:ForestEmulator, :country_damage]))
scc_add = MimiGIVE.compute_scc(m_add, year = 2020)
scc_ces = MimiGIVE.compute_scc(m_ces, year = 2020)
println("Control θ = 1")
println("  SCC additive = ", round(scc_add, digits = 2), "   SCC CES = ", round(scc_ces, digits = 2))
println("  max |damage difference| = ", dmax, " USD2005/yr (should be ≈ 0)")

# --- 2) Table -----------------------------------------------------------------
m0 = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245")
println("\n", rpad("Case", 46), rpad("SCC 2020", 12), "Forest damage 2100 (bn USD2005/yr)")
println(rpad("No forest impacts", 46), rpad(round(MimiGIVE.compute_scc(m0, year = 2020), digits = 2), 12), "—")
m = MimiGIVE.get_model(; base..., forest_damage_weight = nwfp)
run(m)
println(rpad("Consumption only (nwfp)", 46), rpad(round(MimiGIVE.compute_scc(m, year = 2020), digits = 2), 12),
        round(total_damage_2100(m), digits = 2))
for (lab, θ) in [("θ = 1 (constant price)", 1.0), ("θ = 0.41 (Drupp)", 0.41),
                 ("θ = 0 (Cobb-Douglas)", 0.0), ("θ = −1.02 (CWON water)", -1.02)]
    m = MimiGIVE.get_model(; base..., forest_damage_weight = nwfp,
                           forest_utility_weight = nonmkt, forest_ces_theta = θ)
    run(m)
    println(rpad("nwfp + CES, " * lab, 46), rpad(round(MimiGIVE.compute_scc(m, year = 2020), digits = 2), 12),
            round(total_damage_2100(m), digits = 2))
end
