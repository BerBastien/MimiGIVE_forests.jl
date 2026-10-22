# Test of the ceiling calibrated to 2022 (by country).
#
#   julia --project=. scripts/test_forest_value_cap_calibrated.jl
#
# The water value's share of each country's GDP may rise to at most
# mult × its own 2022 share. Below the 2022 share nothing changes, so 2022
# values are exact; above, the share saturates smoothly toward mult × S0.
#
# Checks (SSP245, JULES)
#   1) 2022 water value identical with and without the ceiling.
#   2) 2100 water value and damage shrink as mult goes down.

using MimiGIVE, Mimi

println("mult | water value 2022 | water value 2100 (bn USD2020/yr) | water damage 2100 (bn USD2005/yr)")
for mult in (Inf, 5.0, 3.0, 2.0)
    m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                           forest_emulator = :jules, forest_damages = true,
                           forest_damage_weight = [0.0, 0.0, 0.0, 1.0],
                           forest_value_cap_mult = [Inf, Inf, Inf, mult])
    run(m)
    yrs = collect(Mimi.time_labels(m))
    i22, i100 = findfirst(==(2022), yrs), findfirst(==(2100), yrs)
    v = m[:ForestEmulator, :country_value]
    w22  = sum(skipmissing(v[i22, :, 4])) / 1e3
    w100 = sum(skipmissing(v[i100, :, 4])) / 1e3
    wd   = m[:DamageAggregator, :forest_damage][i100] / 1e9
    println(rpad(mult, 5), "| ", rpad(round(w22, digits = 1), 17), "| ",
            rpad(round(w100, digits = 1), 33), "| ", round(wd, digits = 2))
end
