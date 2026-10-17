# Test of the elasticity uncertainty in the GIVE Monte Carlo.
#
#   julia --project=. scripts/test_elasticity_mcs.jl
#
# Elasticities of the forest service values (CWON 2024 Non-Wood Forest
# Technical Report, Table 5.1, robust SE in parentheses):
#
#              GDP pc          Population
#   rec     0.67 (0.25)      0.72 (0.10)
#   hab     0.79 (0.33)      0.85 (0.20)
#   nwfp   -1.11 (0.17)      0.57 (0.08)
#   wat     2.02 (0.63)      —
#
# Each run uses  e = ê + z · SE,  z ~ Normal(0, 1), one z per service.
#
# 1) Deterministic run (all z = 0): service values must match the result
#    before this change.
# 2) 100-trial Monte Carlo: the 2100 value range per service must be wider
#    than with forest coefficient uncertainty alone, most of all for water.

using MimiGIVE, Mimi, DataFrames, Statistics

services = ["rec", "hab", "nwfp", "wat"]

m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                       forest_emulator = :jules)
run(m)
yrs = collect(Mimi.time_labels(m))
i2100 = findfirst(==(2100), yrs)
V = m[:ForestEmulator, :country_value]          # [time, country, service], 1e6 USD/yr
for s in 1:4
    println("Deterministic 2100 ", services[s], ": ",
            round(sum(skipmissing(V[i2100, :, s])) / 1e3, digits = 1), " billion USD/yr")
end

res = MimiGIVE.run_mcs(trials = 100, m = m, save_trials = true,
                       save_list = [(:ForestEmulator, :country_value)],
                       output_dir = "output/mcs_elasticity_uncertainty")

df = getdataframe(res, :ForestEmulator, :country_value)
df = df[df.time .== 2100, :]
g = combine(groupby(df, [:trialnum, :services]),
            :country_value => (x -> sum(skipmissing(x)) / 1e3) => :value_bn)
for sub in groupby(g, :services)
    q = round.(quantile(sub.value_bn, [0.05, 0.5, 0.95]), digits = 1)
    println("Monte Carlo 2100 ", sub.services[1], " (p5, median, p95): ", q, " billion USD/yr")
end
println("Sampled shocks: output/mcs_elasticity_uncertainty/trials.csv (columns rv_forest_z_e*)")
