# Test of the ceiling on the per-hectare price (income effect), calibrated to 2022.
#
#   julia --project=. scripts/test_forest_price_cap.jl
#
# Per-hectare price factor:  g = (GDPpc_t / GDPpc_2022)^e_gdp
# With e_gdp > 1 the price grows faster than income. The ceiling lets it grow at
# most M times faster than income (smoothly), and leaves it as is when it grows
# slower. Value = mv · g · (pop factor) · F stays proportional to forest area F,
# so losing forest keeps its full cost. Set with
# get_model(forest_price_cap_mult = [rec, hab, nwfp, wat]); Inf = off.
# Default: [Inf, Inf, Inf, 3.0] (water capped at M = 3, central case).
#
# Checks (SSP245, JULES)
#   1) 2022 water value identical in every row (price factor = 1 in 2022).
#   2) 2100 water value shrinks as M goes down.
#   3) water damage stays positive and shrinks roughly in proportion to the value
#      (no collapse or sign flip).

using MimiGIVE, Mimi

println("M    | water value 2022 | water value 2100 (bn USD2020/yr) | water damage 2100 (bn USD2005/yr) | damage/value")
for M in (Inf, 5.0, 3.0, 2.0)
    m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                           forest_emulator = :jules, forest_damages = true,
                           forest_damage_weight = [0.0, 0.0, 0.0, 1.0],
                           forest_price_cap_mult = [Inf, Inf, Inf, M])
    run(m)
    yrs = collect(Mimi.time_labels(m))
    i22, i100 = findfirst(==(2022), yrs), findfirst(==(2100), yrs)
    v = m[:ForestEmulator, :country_value]
    w22  = sum(skipmissing(v[i22, :, 4])) / 1e3
    w100 = sum(skipmissing(v[i100, :, 4])) / 1e3
    wd   = m[:DamageAggregator, :forest_damage][i100] / 1e9
    println(rpad(M, 5), "| ", rpad(round(w22, digits = 1), 17), "| ",
            rpad(round(w100, digits = 1), 33), "| ", rpad(round(wd, digits = 2), 34),
            "| ", round(100 * wd / (w100 * 0.773), digits = 2), "%")
end
