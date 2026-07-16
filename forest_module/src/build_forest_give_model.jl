using Mimi, MimiGIVE

# =============================================================================
# build_forest_give_model.jl  --  wire the forest components into a GIVE model
# -----------------------------------------------------------------------------
# Verified GIVE interfaces (from MimiGIVE v2.2.1-DEV source inspection):
#   * time dimension              : 1750:2300 (annual); damages_first = 2020
#   * country dimension           : 183 ISO3 codes (data/Dimension_countries.csv)
#   * GMST                        : component :temperature, variable :T  (degC)
#   * country population          : :Socioeconomic => :population  ([time,country], million)
#   * country GDP                 : :Socioeconomic => :gdp         ([time,country], billion US$2005/yr)
#   * last GIVE component         : :country_netconsumption
#
# The three forest components are added AFTER :country_netconsumption so that
# :temperature and :Socioeconomic have already run when they execute.
# =============================================================================

const FOREST_FIRST_YEAR = 2020   # == GIVE damages_first

"""
    build_forest_give_model(; kwargs...) -> (model, inputs)

Build a GIVE model with the forest ecosystem-services module attached.

Key keyword arguments (defaults come from ForestConfig):
  socioeconomics_source :: Symbol   (:SSP or :RFF)
  SSP_scenario          :: String   ("SSP245", ...) — used when source == :SSP
  RFFSPsample           :: Union{Int,Nothing}
  coefficient_scenario  :: String   — label only; validated against SSP_scenario
  processed_dir         :: String   — where the compact CSVs live
  forest_change_units   :: Symbol   (:percent / :proportion / :log_change)
  gmst_baseline_year    :: Int      (default 2022)
  es_value_scale        :: Float64
  include_forest_in_total_scc :: Bool
"""
function build_forest_give_model(;
        socioeconomics_source::Symbol = :SSP,
        SSP_scenario::Union{String,Nothing} = "SSP245",
        RFFSPsample::Union{Int,Nothing} = nothing,
        coefficient_scenario::Union{String,Nothing} = nothing,
        processed_dir::String = ForestConfig.PROCESSED_DIR,
        forest_change_units::Symbol = ForestConfig.FOREST_CHANGE_UNITS,
        gmst_baseline_year::Int = ForestConfig.GMST_BASELINE_YEAR,
        gmst_baseline_fallback::Symbol = ForestConfig.GMST_BASELINE_FALLBACK,
        apply_lower_clip::Bool = ForestConfig.APPLY_LOWER_CLIP,
        apply_upper_clip::Bool = ForestConfig.APPLY_UPPER_CLIP,
        es_value_scale::Float64 = ForestConfig.ES_VALUE_SCALE,
        es_value_spec::Symbol = ForestConfig.ES_VALUE_SPEC,
        include_forest_in_total_scc::Bool = ForestConfig.INCLUDE_FOREST_IN_TOTAL_SCC,
        verbose::Bool = true)

    # --- SSP <-> coefficient-scenario consistency (section 16) ---------------
    if coefficient_scenario !== nothing && SSP_scenario !== nothing
        cs = uppercase(coefficient_scenario); ss = uppercase(SSP_scenario)
        (startswith(ss, cs[1:min(4,length(cs))]) || startswith(cs, ss[1:4])) ||
            @warn("coefficient_scenario ($coefficient_scenario) does not obviously match SSP_scenario ($SSP_scenario). Ensure the forest coefficients were fitted for the selected socioeconomic trajectory.")
    end

    # --- base GIVE model -----------------------------------------------------
    m = MimiGIVE.get_model(; socioeconomics_source = socioeconomics_source,
                             SSP_scenario = SSP_scenario, RFFSPsample = RFFSPsample)

    model_countries = String.(Mimi.dim_keys(m, :country))
    model_years     = collect(Mimi.dim_keys(m, :time))

    # --- validate GMST baseline year (section 17) ----------------------------
    baseline_year = _resolve_baseline_year(gmst_baseline_year, FOREST_FIRST_YEAR,
                                           model_years, gmst_baseline_fallback)

    # --- load & validate preprocessed inputs ---------------------------------
    inputs = load_forest_inputs(processed_dir, model_countries; verbose = verbose)

    # --- new dimension -------------------------------------------------------
    set_dimension!(m, :impact_region, inputs.impact_region_ids)

    # --- add components (after the last GIVE component) -----------------------
    add_comp!(m, SpatialSocioeconomics, :ForestSpatialSocioeconomics;
              first = FOREST_FIRST_YEAR, after = :country_netconsumption)
    add_comp!(m, ForestAreaResponse, :ForestAreaResponse;
              first = FOREST_FIRST_YEAR, after = :ForestSpatialSocioeconomics)
    add_comp!(m, ForestEcosystemServices, :ForestEcosystemServices;
              first = FOREST_FIRST_YEAR, after = :ForestAreaResponse)

    # --- Component A: SpatialSocioeconomics ----------------------------------
    connect_param!(m, :ForestSpatialSocioeconomics => :population_country, :Socioeconomic => :population)
    connect_param!(m, :ForestSpatialSocioeconomics => :gdp_country,        :Socioeconomic => :gdp)
    update_param!(m, :ForestSpatialSocioeconomics, :population_share, inputs.population_share)
    update_param!(m, :ForestSpatialSocioeconomics, :gdp_share,        inputs.gdp_share)
    update_param!(m, :ForestSpatialSocioeconomics, :impact_region_country_index, inputs.impact_region_country_index)

    # --- Component B: ForestAreaResponse -------------------------------------
    connect_param!(m, :ForestAreaResponse => :gmst, :temperature => :T)
    update_param!(m, :ForestAreaResponse, :beta1, inputs.beta1)
    update_param!(m, :ForestAreaResponse, :beta2, inputs.beta2)
    update_param!(m, :ForestAreaResponse, :intercept, inputs.intercept)
    update_param!(m, :ForestAreaResponse, :baseline_forest_area, inputs.baseline_forest_area)
    update_param!(m, :ForestAreaResponse, :total_impact_region_area, inputs.total_impact_region_area)
    update_param!(m, :ForestAreaResponse, :delta_gmst_fit_min, inputs.delta_gmst_fit_min)
    update_param!(m, :ForestAreaResponse, :delta_gmst_fit_max, inputs.delta_gmst_fit_max)
    update_param!(m, :ForestAreaResponse, :baseline_year, baseline_year)
    update_param!(m, :ForestAreaResponse, :change_units_code, ForestConfig.forest_change_units_code(forest_change_units))
    update_param!(m, :ForestAreaResponse, :apply_lower_clip, apply_lower_clip)
    update_param!(m, :ForestAreaResponse, :apply_upper_clip, apply_upper_clip)

    # --- Component C: ForestEcosystemServices --------------------------------
    connect_param!(m, :ForestEcosystemServices => :projected_forest_area, :ForestAreaResponse => :projected_forest_area)
    connect_param!(m, :ForestEcosystemServices => :population_ir, :ForestSpatialSocioeconomics => :population_ir)
    connect_param!(m, :ForestEcosystemServices => :gdp_ir,        :ForestSpatialSocioeconomics => :gdp_ir)
    update_param!(m, :ForestEcosystemServices, :baseline_forest_area, inputs.baseline_forest_area)
    update_param!(m, :ForestEcosystemServices, :impact_region_country_index, inputs.impact_region_country_index)
    update_param!(m, :ForestEcosystemServices, :es_value_scale, es_value_scale)
    update_param!(m, :ForestEcosystemServices, :es_value_spec_code, ForestConfig.es_value_spec_code(es_value_spec))

    if include_forest_in_total_scc
        @warn("""
        ============================================================================
        INCLUDE_FOREST_IN_TOTAL_SCC = true
        The forest ecosystem-service value uses the PLACEHOLDER equation
            ESV = scale * forest_area * GDP * population
        which is NOT scientifically calibrated and has NO welfare interpretation.
        Any 'total SCC' that includes it is for architecture testing only and must
        NOT be reported as a policy-relevant value.
        ============================================================================
        """)
    end

    return (model = m, inputs = inputs, baseline_year = baseline_year,
            include_forest_in_total_scc = include_forest_in_total_scc)
end

"""Resolve the GMST baseline year against the model time index and the forest
component's first year, honouring the configured fallback policy."""
function _resolve_baseline_year(baseline_year::Int, forest_first::Int,
                                model_years::Vector{<:Integer}, fallback::Symbol)
    in_range(y) = (y in model_years) && (y >= forest_first)
    if in_range(baseline_year)
        return baseline_year
    end
    if fallback === :error
        error("GMST baseline year $baseline_year is not an available model timestep " *
              ">= $forest_first. Set GMST_BASELINE_FALLBACK to :nearest to override.")
    elseif fallback === :nearest
        candidates = [y for y in model_years if y >= forest_first]
        nearest = candidates[argmin(abs.(candidates .- baseline_year))]
        @warn("GMST baseline year $baseline_year unavailable; using nearest = $nearest.")
        return nearest
    else
        error("Interpolation fallback ($fallback) is not supported by the component " *
              "(it matches an exact timestep). Use :error or :nearest.")
    end
end
