# Test of the asymptotic ceiling on forest service value.
#
#   julia --project=. scripts/test_forest_value_cap.jl
#
# Each country's value of service s is bent toward a ceiling:
#   V = V_raw / (1 + V_raw / (cap_s · GDP_c))
#   ≈ V_raw when V_raw is small, → cap_s · GDP_c (never above) when large.
# Set with get_model(forest_value_cap = [rec, hab, nwfp, wat]); Inf = no ceiling
# (default, so results are unchanged unless a cap is given).
#
# Checks (SSP245, JULES, year 2100)
#   1) cap = Inf  → same water value and damage as before.
#   2) smaller caps → water value and water damage shrink; other services untouched.

using MimiGIVE, Mimi

yr = 2100
println("cap (share of GDP) | water value (bn USD2020/yr) | water damage (bn USD2005/yr)")
for cap in (Inf, 0.05, 0.02, 0.01)
    m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                           forest_emulator = :jules, forest_damages = true,
                           forest_damage_weight = [0.0, 0.0, 0.0, 1.0],
                           forest_value_cap = [Inf, Inf, Inf, cap])
    run(m)
    i = findfirst(==(yr), collect(Mimi.time_labels(m)))
    wv = sum(skipmissing(m[:ForestEmulator, :country_value][i, :, 4])) / 1e3
    wd = m[:DamageAggregator, :forest_damage][i] / 1e9
    println(rpad(cap, 19), "| ", rpad(round(wv, digits = 1), 28), "| ", round(wd, digits = 2))
end
