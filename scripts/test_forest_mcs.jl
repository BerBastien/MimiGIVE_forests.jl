# Prueba de la incertidumbre del emulador de bosques en el Monte Carlo de GIVE.
#
#   julia --project=. scripts/test_forest_mcs.jl
#
# 1) Corrida determinística: debe dar lo mismo que antes del cambio
#    (≈ 4,282.5 Mha en 2022 y ≈ 4,374.5 Mha en 2100 con SSP2-4.5 y JULES).
# 2) Monte Carlo de 100 corridas: los choques z deben aparecer como dados
#    y el rango del bosque en 2100 debe ser más ancho que sin ellos (~27 Mha).

using MimiGIVE, Mimi, DataFrames, Statistics

m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                       forest_emulator = :jules)
run(m)
F = m[:ForestEmulator, :global_forest_Mha]
yrs = collect(Mimi.time_labels(m))
println("Determinística 2022: ", round(F[findfirst(==(2022), yrs)], digits = 1), " Mha")
println("Determinística 2100: ", round(F[findfirst(==(2100), yrs)], digits = 1), " Mha")

res = MimiGIVE.run_mcs(trials = 100, m = m, save_trials = true,
                       save_list = [(:ForestEmulator, :global_forest_Mha), (:temperature, :T)],
                       output_dir = "output/mcs_bosques_incertidumbre")

df = getdataframe(res, :ForestEmulator, :global_forest_Mha)
f2100 = df[df.time .== 2100, :global_forest_Mha]
println("Monte Carlo 2100 (p5, mediana, p95): ", round.(quantile(f2100, [0.05, 0.5, 0.95]), digits = 1))
println("Ancho p5–p95: ", round(quantile(f2100, 0.95) - quantile(f2100, 0.05), digits = 1), " Mha")
println("Los valores sorteados de z quedan en output/mcs_bosques_incertidumbre/trials.csv (columnas rv_forest_z_*)")
