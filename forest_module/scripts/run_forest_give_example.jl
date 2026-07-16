# =============================================================================
# run_forest_give_example.jl  --  end-to-end deterministic example
# -----------------------------------------------------------------------------
# Requires the full MimiGIVE stack to be installed & instantiated (it builds a real
# GIVE model). Run from the MimiGIVE.jl repo root with its project active:
#
#     julia --project=. forest_module/scripts/run_forest_give_example.jl
#
# If no processed inputs exist yet, synthetic fixtures are generated so the example
# runs out-of-the-box. Replace them with real preprocessed data (see the two
# preprocess_*.py scripts) for a real run.
# =============================================================================

const HERE = @__DIR__
const FOREST_ROOT = normpath(joinpath(HERE, ".."))
const PROCESSED = joinpath(FOREST_ROOT, "data", "processed")

# 1. Ensure inputs exist ------------------------------------------------------
if !isfile(joinpath(PROCESSED, "impact_region_country_crosswalk.csv"))
    @info "No processed forest inputs found — generating synthetic fixtures."
    include(joinpath(HERE, "make_synthetic_fixtures.jl"))
end

# 2. Load the module ----------------------------------------------------------
include(joinpath(FOREST_ROOT, "ForestGIVE.jl"))
using .ForestGIVE

# 3. Build a deterministic GIVE model with the forest module attached ---------
@info "Building GIVE + forest model (SSP245)…"
built = build_forest_give_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245")

# 4. Run the base (deterministic) model ---------------------------------------
@info "Running base model…"
run(built.model)

# quick sanity peek
let m = built.model
    println("delta_gmst @2022 (should be 0): ",
            m[:ForestAreaResponse, :delta_gmst][findfirst(==(2022), collect(1750:2300))])
end

# 5. Forest SCC via a marginal (pulse) model ----------------------------------
@info "Computing forest SCC (pulse in 2030)…"
scc = compute_forest_scc(built; year = 2030, gas = :CO2, pulse_size = 1.0,
                         prtp = 0.015, eta = 1.45)

# 6. Write tidy outputs -------------------------------------------------------
@info "Writing outputs…"
write_join_diagnostics(built)
write_temperature_diagnostic(built)
write_impact_region_output(built; years = [2030, 2050, 2100])
write_country_output(built; years = 2020:10:2100)
write_marginal_country_output(scc)
write_summary_output(built, scc)
write_clipping_diagnostic(built)

println("\nForest SCC (PLACEHOLDER units / tonne CO2): ", scc.scc_global)
println("Outputs written to: ", ForestGIVE.ForestConfig.OUTPUT_DIRECTORY)
println("\nNOTE: the ecosystem-service value is a placeholder (area*GDP*pop). The SCC")
println("here is for architecture testing only — NOT a policy-relevant number.")
