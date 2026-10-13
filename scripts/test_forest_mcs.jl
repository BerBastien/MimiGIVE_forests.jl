# Test of the forest emulator uncertainty in the GIVE Monte Carlo.
#
#   julia --project=. scripts/test_forest_mcs.jl
#
# 1) Deterministic run: must match the result before this change
#    (≈ 4,282.5 Mha in 2022 and ≈ 4,374.7 Mha in 2100 with SSP2-4.5 and JULES).
# 2) 100-trial Monte Carlo: the z shocks must appear as random variables
#    and the 2100 forest range must be wider than without them (~27 Mha).

using MimiGIVE, Mimi, DataFrames, Statistics

m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                       forest_emulator = :jules)
run(m)
F = m[:ForestEmulator, :global_forest_Mha]
yrs = collect(Mimi.time_labels(m))
println("Deterministic 2022: ", round(F[findfirst(==(2022), yrs)], digits = 1), " Mha")
println("Deterministic 2100: ", round(F[findfirst(==(2100), yrs)], digits = 1), " Mha")

res = MimiGIVE.run_mcs(trials = 100, m = m, save_trials = true,
                       save_list = [(:ForestEmulator, :global_forest_Mha), (:temperature, :T)],
                       output_dir = "output/mcs_forest_uncertainty")

df = getdataframe(res, :ForestEmulator, :global_forest_Mha)
f2100 = df[df.time .== 2100, :global_forest_Mha]
println("Monte Carlo 2100 (p5, median, p95): ", round.(quantile(f2100, [0.05, 0.5, 0.95]), digits = 1))
println("p5–p95 width: ", round(quantile(f2100, 0.95) - quantile(f2100, 0.05), digits = 1), " Mha")
println("Sampled z values are in output/mcs_forest_uncertainty/trials.csv (columns rv_forest_z_*)")
