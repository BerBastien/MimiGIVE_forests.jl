# Water: price ceiling M = 3 vs. using the Drupp et al. (2025) elasticity directly.
#
#   julia --project=. scripts/compare_water_elasticity.jl
#
# Cases (SSP245, JULES, 2100, deterministic; damages count water only):
#   CWON e = 2.02, no ceiling      (original)
#   CWON e = 2.02, M = 3           (current default)
#   Drupp e = 0.59, no ceiling     (ecosystem services overall)
#   Drupp e = 0.40, no ceiling     (water purification)
# Only the water GDP elasticity changes; the other services keep CWON values.

using MimiGIVE, Mimi

e_cwon = [0.67, 0.79, -1.11, 2.02]
cases = [("CWON 2.02, no ceiling", 2.02, Inf),
         ("CWON 2.02, M = 3",      2.02, 3.0),
         ("Drupp 0.59, no ceiling", 0.59, Inf),
         ("Drupp 0.40, no ceiling", 0.40, Inf)]

println(rpad("Case", 26), "| water value 2100 (bn USD2020/yr) | water damage 2100 (bn USD2005/yr)")
for (name, e_wat, M) in cases
    m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                           forest_emulator = :jules, forest_damages = true,
                           forest_damage_weight = [0.0, 0.0, 0.0, 1.0],
                           forest_price_cap_mult = [Inf, Inf, Inf, M])
    update_param!(m, :ForestEmulator, :e_gdp, [e_cwon[1:3]; e_wat])
    run(m)
    i = findfirst(==(2100), collect(Mimi.time_labels(m)))
    wv = sum(skipmissing(m[:ForestEmulator, :country_value][i, :, 4])) / 1e3
    wd = m[:DamageAggregator, :forest_damage][i] / 1e9
    println(rpad(name, 26), "| ", rpad(round(wv, digits = 1), 34), "| ", round(wd, digits = 2))
end
