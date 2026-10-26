# SCC 2020 with forest impacts through consumption, "utility" and both.
#
#   julia --project=. scripts/scc_forest_table.jl
#
# Provisional utility case: perfect substitution (θ = 1), U = C + E, which is the
# same as adding the non-market services (rec, hab, wat) to consumption damages.
# To be replaced by the CES utility U = [(1−α)C^θ + αE^θ]^(1/θ) when ready.
# Water uses the default price ceiling (M = 3). SSP245, JULES, deterministic.

using MimiGIVE, Mimi

base = (socioeconomics_source = :SSP, SSP_scenario = "SSP245", forest_emulator = :jules)
cases = [("No forest impacts",                 false, [0.0, 0.0, 0.0, 0.0]),
         ("Consumption only (nwfp)",           true,  [0.0, 0.0, 1.0, 0.0]),
         ("Utility only, θ = 1 (rec, hab, wat)", true, [1.0, 1.0, 0.0, 1.0]),
         ("Both (all 4 services)",             true,  [1.0, 1.0, 1.0, 1.0])]

println(rpad("Case", 40), "SCC 2020 (USD2005/tCO2)")
for (name, dam, w) in cases
    m = MimiGIVE.get_model(; base..., forest_damages = dam, forest_damage_weight = w)
    scc = MimiGIVE.compute_scc(m, year = 2020)
    println(rpad(name, 40), round(scc, digits = 2))
end
