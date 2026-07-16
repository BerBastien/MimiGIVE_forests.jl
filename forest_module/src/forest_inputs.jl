using DataFrames, CSVFiles, FileIO, Statistics

# =============================================================================
# forest_inputs.jl  --  Load & validate the compact, preprocessed CSVs
# -----------------------------------------------------------------------------
# The heavy geospatial work happens once, in Python (scripts/preprocess_*.py),
# which writes small tidy CSVs into PROCESSED_DIR. Here we:
#   * read them,
#   * align every table to a single canonical, ordered list of impact regions,
#   * build the impact_region -> GIVE-country integer index (position-based, as
#     Mimi requires),
#   * validate share sums (== 1 per country), normalising only tiny numerical
#     discrepancies and erroring on large ones,
#   * emit a join-diagnostics report.
#
# R analogy: this is the `readr::read_csv()` + `dplyr` join/validate block that
# produces clean, index-aligned vectors ready to hand to the model.
# =============================================================================

struct ForestInputs
    impact_region_ids::Vector{String}
    give_country_ids::Vector{String}            # == model :country dimension keys
    impact_region_country_index::Vector{Int}    # 1-based index into give_country_ids
    beta1::Vector{Float64}
    beta2::Vector{Float64}
    intercept::Vector{Float64}
    baseline_forest_area::Vector{Float64}
    total_impact_region_area::Vector{Float64}
    delta_gmst_fit_min::Vector{Float64}
    delta_gmst_fit_max::Vector{Float64}
    population_share::Vector{Float64}
    gdp_share::Vector{Float64}
    diagnostics::DataFrame
    summary::Dict{String,Any}
end

_readcsv(path) = DataFrame(load(path))

_hascol(df, name) = name in names(df)

"""Reindex `df` (keyed by `key`) onto `id_order`, returning a Float64 vector of
`col` for each id; error listing any missing ids."""
function _align(df::DataFrame, key::String, id_order::Vector{String}, col::String; what::String)
    _hascol(df, key) || error("Column `$key` not found in $what (columns: $(names(df))).")
    _hascol(df, col) || error("Column `$col` not found in $what (columns: $(names(df))).")
    lut = Dict{String,Float64}()
    for r in eachrow(df)
        lut[string(r[key])] = Float64(r[col])
    end
    missing_ids = [id for id in id_order if !haskey(lut, id)]
    isempty(missing_ids) || error("$what is missing $(length(missing_ids)) impact regions, e.g. $(first(missing_ids, min(5,length(missing_ids)))).")
    return [lut[id] for id in id_order]
end

_align_optional(df, key, id_order, col; default) =
    _hascol(df, col) ? _align(df, key, id_order, col; what="optional column $col") :
                       fill(default, length(id_order))

"""
    load_forest_inputs(processed_dir, model_countries; share_tol=1e-6,
                       share_hard_tol=1e-3, verbose=true)

Load and validate all forest-module inputs from `processed_dir`, aligned to
`model_countries` (the GIVE `:country` dimension, in order). Returns `ForestInputs`.
"""
function load_forest_inputs(processed_dir::String, model_countries::Vector{String};
                            share_tol::Float64 = 1e-6,
                            share_hard_tol::Float64 = 1e-3,
                            verbose::Bool = true)

    xwalk_path = joinpath(processed_dir, "impact_region_country_crosswalk.csv")
    coeff_path = joinpath(processed_dir, "forest_coefficients.csv")
    area_path  = joinpath(processed_dir, "baseline_forest_area.csv")
    pop_path   = joinpath(processed_dir, "population_shares.csv")
    gdp_path   = joinpath(processed_dir, "gdp_shares.csv")

    for (p, nm) in [(xwalk_path,"crosswalk"), (coeff_path,"coefficients"),
                    (area_path,"baseline area"), (pop_path,"population shares"),
                    (gdp_path,"gdp shares")]
        isfile(p) || error("Missing $nm file: $p\nRun the preprocessing scripts (or make_synthetic_fixtures.jl) first.")
    end

    # ---- crosswalk defines the canonical impact-region ordering -------------
    xwalk = _readcsv(xwalk_path)
    _hascol(xwalk, "impact_region_id") || error("crosswalk must have column impact_region_id")
    _hascol(xwalk, "give_country_id")  || error("crosswalk must have column give_country_id")

    ir_ids  = string.(xwalk.impact_region_id)
    ir_ctry = string.(xwalk.give_country_id)

    # ---- join diagnostics ---------------------------------------------------
    dup_ir = [id for id in unique(ir_ids) if count(==(id), ir_ids) > 1]
    isempty(dup_ir) || error("Duplicated impact_region_id in crosswalk: $(first(dup_ir, min(10,length(dup_ir)))).")

    country_pos = Dict(c => i for (i, c) in enumerate(model_countries))
    unmatched_ir = [ir_ids[i] for i in eachindex(ir_ids) if !haskey(country_pos, ir_ctry[i])]
    if !isempty(unmatched_ir)
        bad = unique([ir_ctry[i] for i in eachindex(ir_ids) if !haskey(country_pos, ir_ctry[i])])
        error("$(length(unmatched_ir)) impact regions map to country IDs not present in GIVE: $(first(bad, min(10,length(bad)))). " *
              "Check that crosswalk give_country_id uses GIVE ISO3 codes.")
    end

    ir_country_index = [country_pos[c] for c in ir_ctry]

    # per-country impact-region counts (incl. GIVE countries with none)
    counts = Dict{String,Int}()
    for c in ir_ctry; counts[c] = get(counts, c, 0) + 1; end
    countries_no_ir = [c for c in model_countries if get(counts, c, 0) == 0]

    # ---- coefficients -------------------------------------------------------
    coeff = _readcsv(coeff_path)
    beta1 = _align(coeff, "impact_region_id", ir_ids, "beta_delta_gmst"; what="coefficients")
    beta2 = _align(coeff, "impact_region_id", ir_ids, "beta_delta_gmst_squared"; what="coefficients")
    intercept = _align_optional(coeff, "impact_region_id", ir_ids, "intercept"; default=0.0)
    dtmin = _align_optional(coeff, "impact_region_id", ir_ids, "delta_gmst_fit_min"; default=-Inf)
    dtmax = _align_optional(coeff, "impact_region_id", ir_ids, "delta_gmst_fit_max"; default= Inf)

    # ---- baseline forest area ----------------------------------------------
    area = _readcsv(area_path)
    baseline_area = _align(area, "impact_region_id", ir_ids, "baseline_forest_area"; what="baseline area")
    total_area = _align_optional(area, "impact_region_id", ir_ids, "total_impact_region_area"; default=0.0)

    # ---- shares -------------------------------------------------------------
    popdf = _readcsv(pop_path)
    gdpdf = _readcsv(gdp_path)
    pop_share = _align(popdf, "impact_region_id", ir_ids, "population_share"; what="population shares")
    gdp_share = _align(gdpdf, "impact_region_id", ir_ids, "gdp_share"; what="gdp shares")

    # ---- validate & normalise share sums per country ------------------------
    pop_share, pop_share_report = _validate_shares!(pop_share, ir_country_index, model_countries,
                                                    "population"; tol=share_tol, hard_tol=share_hard_tol)
    gdp_share, gdp_share_report = _validate_shares!(gdp_share, ir_country_index, model_countries,
                                                    "gdp"; tol=share_tol, hard_tol=share_hard_tol)

    # ---- diagnostics table (per impact region) ------------------------------
    diagnostics = DataFrame(
        impact_region_id = ir_ids,
        give_country_id  = ir_ctry,
        country_index    = ir_country_index,
        beta_delta_gmst  = beta1,
        beta_delta_gmst_squared = beta2,
        intercept        = intercept,
        baseline_forest_area = baseline_area,
        population_share = pop_share,
        gdp_share        = gdp_share,
    )

    summary = Dict{String,Any}(
        "n_impact_regions" => length(ir_ids),
        "n_give_countries" => length(model_countries),
        "n_countries_with_impact_regions" => length(model_countries) - length(countries_no_ir),
        "countries_with_no_impact_regions" => countries_no_ir,
        "duplicate_impact_region_ids" => dup_ir,
        "population_share_max_abs_discrepancy" => pop_share_report.max_abs,
        "gdp_share_max_abs_discrepancy" => gdp_share_report.max_abs,
        "population_shares_sum_to_one" => pop_share_report.ok,
        "gdp_shares_sum_to_one" => gdp_share_report.ok,
        "impact_region_count_by_country" => counts,
    )

    if verbose
        println("── Forest inputs: index-alignment diagnostics ───────────────")
        println("  impact regions            : ", summary["n_impact_regions"])
        println("  GIVE countries            : ", summary["n_give_countries"])
        println("  countries with >=1 IR     : ", summary["n_countries_with_impact_regions"])
        println("  countries with 0 IRs      : ", length(countries_no_ir),
                length(countries_no_ir) == 0 ? "" : "  e.g. $(first(countries_no_ir, min(8,length(countries_no_ir))))")
        println("  duplicate IR ids          : ", length(dup_ir))
        println("  pop shares sum to 1       : ", pop_share_report.ok,
                "  (max |sum-1| = ", round(pop_share_report.max_abs, sigdigits=3), ")")
        println("  gdp shares sum to 1       : ", gdp_share_report.ok,
                "  (max |sum-1| = ", round(gdp_share_report.max_abs, sigdigits=3), ")")
        println("─────────────────────────────────────────────────────────────")
    end

    return ForestInputs(ir_ids, model_countries, ir_country_index, beta1, beta2, intercept,
                        baseline_area, total_area, dtmin, dtmax, pop_share, gdp_share,
                        diagnostics, summary)
end

"""Normalise per-country share sums to 1 when within `hard_tol`; error otherwise.
Returns the (possibly normalised) shares and a report NamedTuple."""
function _validate_shares!(shares::Vector{Float64}, country_index::Vector{Int},
                           model_countries::Vector{String}, label::String;
                           tol::Float64, hard_tol::Float64)
    ncountry = length(model_countries)
    sums = zeros(ncountry)
    for i in eachindex(shares)
        sums[country_index[i]] += shares[i]
    end

    # only consider countries that actually have impact regions (sum > 0)
    active = findall(>(0.0), sums)
    discrepancies = [abs(sums[c] - 1.0) for c in active]
    max_abs = isempty(discrepancies) ? 0.0 : maximum(discrepancies)

    # hard failure: any active country off by more than hard_tol
    bad = [model_countries[c] for c in active if abs(sums[c] - 1.0) > hard_tol]
    if !isempty(bad)
        worst = active[argmax(discrepancies)]
        error("$label shares do not sum to 1 for $(length(bad)) countries (worst: " *
              "$(model_countries[worst]) sums to $(round(sums[worst], sigdigits=6))). " *
              "This exceeds the hard tolerance $hard_tol; refusing to silently normalise invalid weights.")
    end

    # soft normalisation for tiny rounding discrepancies
    normalised = false
    for c in active
        if abs(sums[c] - 1.0) > tol
            normalised = true
        end
    end
    if normalised
        for i in eachindex(shares)
            s = sums[country_index[i]]
            s > 0 && (shares[i] = shares[i] / s)
        end
        @warn("$label shares were re-normalised within-country (max original |sum-1| = $(round(max_abs, sigdigits=4))).")
    end

    ok = max_abs <= hard_tol
    return shares, (ok=ok, max_abs=max_abs)
end
