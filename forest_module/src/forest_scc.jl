using Mimi, MimiGIVE

# =============================================================================
# forest_scc.jl  --  forest ecosystem-service SCC via a GIVE marginal model
# -----------------------------------------------------------------------------
# Reuses GIVE's own marginal-model and discounting machinery so the forest sector
# is consistent with the rest of GIVE:
#   * MimiGIVE.get_marginal_model(m; year, gas, pulse_size) builds base + pulse and
#     sets the MarginalModel `delta` = pulse_size * pulse-unit-conversion, so
#     mm[:comp, :var] returns marginal damage PER UNIT of the pulse gas.
#   * We multiply by the molecular-mass conversion (12/44 for CO2, etc.) exactly as
#     _compute_scc does, giving marginal damage per TONNE of the gas.
#   * We discount with GIVE's non-equity-weighted Ramsey convention using the base
#     model's global net consumption per capita:
#       df_i = (cpc[year]/cpc[i])^eta * 1/(1+prtp)^(t-year)
#
# ############################################################################
# ##  UNITS WARNING: es_damage is in PLACEHOLDER units (Mha * billion$ *      ##
# ##  million-persons * scale). The resulting "forest SCC" is therefore in    ##
# ##  placeholder-units per tonne of gas, NOT dollars per tonne. It exists to  ##
# ##  prove the pipeline, not to be reported.                                  ##
# ############################################################################
# =============================================================================

const _FOREST_MODEL_YEARS = collect(1750:2300)

"""
    compute_forest_scc(built; year, gas=:CO2, pulse_size=1.0, prtp=0.015,
                       eta=1.45, last_year=2300, verbose=true)

`built` is either the NamedTuple returned by `build_forest_give_model` or a bare
GIVE `Model` that already has the forest components attached.

Returns a NamedTuple with, among others:
  scc_global                :: Float64                (placeholder-units / tonne gas)
  scc_country               :: Vector{Float64}        (per GIVE country)
  countries                 :: Vector{String}
  marginal_damage_country   :: Matrix (time x country)  (already * molecular conv.)
  forest_damage_base        :: Matrix (time x country)
  forest_damage_pulse       :: Matrix (time x country)
  discount_factors          :: Vector{Float64}
"""
function compute_forest_scc(built; year::Int, gas::Symbol = :CO2,
                            pulse_size::Float64 = 1.0,
                            prtp::Float64 = 0.015, eta::Float64 = 1.45,
                            last_year::Int = 2300, verbose::Bool = true)

    m = built isa NamedTuple ? built.model : built

    year in _FOREST_MODEL_YEARS      || error("year $year not in model time index.")
    last_year in _FOREST_MODEL_YEARS || error("last_year $last_year not in model time index.")

    mm = MimiGIVE.get_marginal_model(m; year = year, gas = gas, pulse_size = pulse_size)
    run(mm)

    year_index      = findfirst(==(year), _FOREST_MODEL_YEARS)
    last_year_index = findfirst(==(last_year), _FOREST_MODEL_YEARS)
    conv            = MimiGIVE.scc_gas_molecular_conversions[gas]

    countries   = String.(Mimi.dim_keys(mm.base, :country))
    dmg_base    = mm.base[:ForestEcosystemServices, :es_damage_country]      # [time, country]
    dmg_pulse   = mm.modified[:ForestEcosystemServices, :es_damage_country]

    # mm[...] = (modified - base) / delta ; * molecular conversion -> per tonne.
    # Pre-first-year cells are `missing`; coerce to a plain Float64 matrix over the
    # damages horizon (1750..2019 are irrelevant and not used in the SCC sum).
    md_raw = mm[:ForestEcosystemServices, :es_damage_country] .* conv         # [time, country]
    md_country = coalesce.(md_raw, 0.0)                                        # dense Float64, shape preserved

    # marginal damages before the emissions year should be zero
    if year_index > 1
        md_country[1:year_index-1, :] .= 0.0
    end

    # --- GIVE non-equity-weighted discounting --------------------------------
    cpc = mm.base[:global_netconsumption, :net_cpc]
    df = Float64[ (cpc[year_index] / cpc[i])^eta * 1 / (1 + prtp)^(t - year)
                  for (i, t) in enumerate(_FOREST_MODEL_YEARS) if year <= t <= last_year ]

    md_slice    = md_country[year_index:last_year_index, :]     # [T, country]
    scc_country = vec(sum(df .* md_slice, dims = 1))            # sum_t df_t * md_{t,c}
    scc_global  = sum(scc_country)

    if verbose
        println("── Forest SCC (PLACEHOLDER units / tonne $gas) ──────────────")
        println("  pulse year            : ", year)
        println("  discounting           : prtp=$prtp, eta=$eta (GIVE Ramsey, global cpc)")
        println("  forest SCC (global)   : ", scc_global)
        top = sortperm(scc_country, rev=true)[1:min(5, length(scc_country))]
        println("  top-5 contributing countries:")
        for c in top
            println("      ", countries[c], "  ", scc_country[c])
        end
        println("─────────────────────────────────────────────────────────────")
    end

    return (scc_global = scc_global, scc_country = scc_country, countries = countries,
            marginal_damage_country = md_country,
            forest_damage_base = dmg_base, forest_damage_pulse = dmg_pulse,
            discount_factors = df, year = year, gas = gas, pulse_size = pulse_size,
            prtp = prtp, eta = eta, last_year = last_year, mm = mm)
end

"""
    compute_total_scc_with_forest(built; year, gas=:CO2, kwargs...)

Prototype total-SCC-including-forest. Computes GIVE's official SCC via
`MimiGIVE.compute_scc` and ADDS the forest sector's SCC after converting the
placeholder units to USD2005 with `forest_damage_to_usd2005`. Guarded by the
`INCLUDE_FOREST_IN_TOTAL_SCC` switch and prints a prominent warning: the placeholder
valuation is NOT calibrated, so the combined number is for architecture testing only.

This is intentionally non-destructive — it does not modify GIVE's DamageAggregator.
"""
function compute_total_scc_with_forest(built; year::Int, gas::Symbol = :CO2,
                                       prtp::Float64 = 0.015, eta::Float64 = 1.45,
                                       last_year::Int = 2300,
                                       forest_damage_to_usd2005::Float64 = ForestConfig.FOREST_DAMAGE_TO_USD2005,
                                       verbose::Bool = true)
    m = built isa NamedTuple ? built.model : built
    include_flag = built isa NamedTuple ? built.include_forest_in_total_scc :
                                          ForestConfig.INCLUDE_FOREST_IN_TOTAL_SCC
    if !include_flag
        error("compute_total_scc_with_forest called but INCLUDE_FOREST_IN_TOTAL_SCC is false. " *
              "Set include_forest_in_total_scc=true in build_forest_give_model to opt in.")
    end

    @warn("""
    Building a TOTAL SCC that includes the PLACEHOLDER forest valuation.
    This number is NOT policy-relevant. Replace `ecosystem_service_value` with a
    defensible valuation before interpreting any combined SCC.
    """)

    official = MimiGIVE.compute_scc(m; year = year, gas = gas, prtp = prtp, eta = eta,
                                    last_year = last_year)
    forest = compute_forest_scc(built; year = year, gas = gas, prtp = prtp, eta = eta,
                                last_year = last_year, verbose = verbose)

    forest_scc_usd = forest.scc_global * forest_damage_to_usd2005
    return (total_scc = official + forest_scc_usd,
            official_scc = official,
            forest_scc_placeholder = forest.scc_global,
            forest_scc_usd = forest_scc_usd,
            forest_damage_to_usd2005 = forest_damage_to_usd2005)
end
