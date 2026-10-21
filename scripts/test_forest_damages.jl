# Test of forest market damages in GIVE (consumption channel).
#
#   julia --project=. scripts/test_forest_damages.jl
#
# Forest damage of country c in year t (USD2005/yr):
#   D_c,t = Σ_s w_s · (V_nocl − V)_c,s,t · 1e6 · deflator
#     V_nocl  value with forest fixed at 2022 (socioeconomics evolve)
#     V       value with the climate-driven forest
#     w_s     share of service s counted as a market loss
#             default [rec, hab, nwfp, wat] = [0, 0, 1, 0]
# D enters DamageAggregator like agriculture or energy and lowers net
# consumption, so it reaches the SCC. The emulator runs 2022–2100; damage is 0
# outside that span.
#
# Checks
#   1) forest_damages = false → SCC identical to the model without the change.
#   2) forest_damages = true  → forest damages > 0 and SCC a bit higher.

using MimiGIVE, Mimi

args = (socioeconomics_source = :SSP, SSP_scenario = "SSP245", forest_emulator = :jules)

m0 = MimiGIVE.get_model(; args..., forest_damages = false)
m1 = MimiGIVE.get_model(; args..., forest_damages = true)
run(m0); run(m1)

yrs = collect(Mimi.time_labels(m1))
fd  = m1[:DamageAggregator, :forest_damage]
gdp = sum(m1[:Socioeconomic, :gdp], dims = 2)[:] .* 1e9          # USD2005/yr
for y in (2030, 2050, 2100)
    i = findfirst(==(y), yrs)
    println("Forest damage $y: ", round(fd[i] / 1e9, digits = 2), " billion USD2005/yr (",
            round(100 * fd[i] / gdp[i], digits = 4), "% of GDP)")
end

i = findfirst(==(2100), yrs)
println("Total damage 2100 without / with forest: ",
        round(m0[:DamageAggregator, :total_damage][i] / 1e9, digits = 1), " / ",
        round(m1[:DamageAggregator, :total_damage][i] / 1e9, digits = 1), " billion USD2005/yr")

scc0 = MimiGIVE.compute_scc(m0, year = 2020)
scc1 = MimiGIVE.compute_scc(m1, year = 2020)
println("SCC 2020 without forest damages: ", round(scc0, digits = 2), " USD2005/tCO2")
println("SCC 2020 with forest damages:    ", round(scc1, digits = 2), " USD2005/tCO2")
