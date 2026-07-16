# =============================================================================
# make_synthetic_fixtures.jl
# -----------------------------------------------------------------------------
# Writes small, VALID synthetic versions of the five compact input CSVs into
# forest_module/data/processed/ so the whole pipeline (build -> run -> SCC ->
# outputs -> tests) is runnable BEFORE the real spatial data is plugged in.
#
# The synthetic impact regions are assigned to REAL GIVE ISO3 countries (read from
# MimiGIVE's Dimension_countries.csv) so build_forest_give_model can attach them.
#
# Run:  julia --project=. forest_module/scripts/make_synthetic_fixtures.jl
# =============================================================================

using DataFrames, CSVFiles, FileIO, Random

const HERE = @__DIR__
const FOREST_ROOT = normpath(joinpath(HERE, ".."))
const REPO_ROOT   = normpath(joinpath(FOREST_ROOT, ".."))
const PROCESSED   = joinpath(FOREST_ROOT, "data", "processed")
mkpath(PROCESSED)

Random.seed!(42)

# Real GIVE countries (ISO3) so the crosswalk matches the model dimension.
countries = string.(DataFrame(load(joinpath(REPO_ROOT, "data", "Dimension_countries.csv"))).CountryISO)

# Use the first 8 GIVE countries; give each 2–4 synthetic impact regions.
sel_countries = countries[1:8]
ir_ids   = String[]
ir_ctry  = String[]
for c in sel_countries
    nreg = rand(2:4)
    for k in 1:nreg
        push!(ir_ids, string(c, "_IR", k))
        push!(ir_ctry, c)
    end
end
n = length(ir_ids)

# ---- crosswalk --------------------------------------------------------------
save(joinpath(PROCESSED, "impact_region_country_crosswalk.csv"),
     DataFrame(impact_region_id = ir_ids, give_country_id = ir_ctry,
               country_name = ir_ctry))

# ---- coefficients (percent units): mild forest loss with warming -----------
# change% = b1*dT + b2*dT^2 ; choose b1<0 (loss), small b2.
save(joinpath(PROCESSED, "forest_coefficients.csv"),
     DataFrame(impact_region_id = ir_ids,
               beta_delta_gmst = round.(-2.0 .- 1.5 .* rand(n), digits=3),
               beta_delta_gmst_squared = round.(-0.2 .* rand(n), digits=3),
               intercept = zeros(n),
               delta_gmst_fit_min = fill(-1.0, n),
               delta_gmst_fit_max = fill( 5.0, n)))

# ---- baseline forest area (Mha) + total region area ------------------------
base_area  = round.(1.0 .+ 9.0 .* rand(n), digits=3)      # 1–10 Mha
total_area = round.(base_area .* (1.5 .+ rand(n)), digits=3)
save(joinpath(PROCESSED, "baseline_forest_area.csv"),
     DataFrame(impact_region_id = ir_ids,
               baseline_forest_area = base_area,
               total_impact_region_area = total_area))

# ---- population & gdp shares (sum to 1 within each country) -----------------
function shares_for(ids, ctry)
    pop = zeros(length(ids)); gdp = zeros(length(ids))
    for c in unique(ctry)
        idx = findall(==(c), ctry)
        wp = rand(length(idx)); wp ./= sum(wp)
        wg = rand(length(idx)); wg ./= sum(wg)
        pop[idx] .= wp; gdp[idx] .= wg
    end
    return pop, gdp
end
pop_share, gdp_share = shares_for(ir_ids, ir_ctry)
save(joinpath(PROCESSED, "population_shares.csv"),
     DataFrame(impact_region_id = ir_ids, give_country_id = ir_ctry,
               population_share = round.(pop_share, digits=6)))
save(joinpath(PROCESSED, "gdp_shares.csv"),
     DataFrame(impact_region_id = ir_ids, give_country_id = ir_ctry,
               gdp_share = round.(gdp_share, digits=6)))

println("Wrote synthetic fixtures for $n impact regions across $(length(sel_countries)) countries to:")
println("  ", PROCESSED)
