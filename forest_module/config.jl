# =============================================================================
# config.jl  --  Configuration for the MimiGIVE Forest Ecosystem-Services module
# =============================================================================
#
# This file centralises every path and every tunable knob for the forest module.
# EDIT THE PATHS BELOW to point at your real data before running the preprocessing
# scripts.  Nothing in the runtime Julia components reads these files directly --
# the Python preprocessing step (scripts/preprocess_forest_inputs.py) converts the
# spatial data into small, validated CSVs that Julia loads once at model-build time
# (see section 14 of the spec: "do not do expensive geospatial reads inside
# run_timestep").
#
# R analogy: think of this like an `config.R` / `here::here()` block at the top of
# an analysis --- one place to set paths and options, sourced by everything else.
# =============================================================================

module ForestConfig

# -----------------------------------------------------------------------------
# 1. RAW INPUT PATHS  (consumed by the Python preprocessing scripts, not Julia)
# -----------------------------------------------------------------------------
# These can be CSV, GeoPackage (.gpkg), Shapefile, GeoParquet, or NetCDF; the
# preprocessing script inspects each file before assuming its structure.

const FOREST_COEFFICIENTS_PATH = "/Users/dr.bastien/Library/CloudStorage/GoogleDrive-lab.ecoclim.unam@gmail.com/My Drive/PFT patterns/Claude_pft_patterns/emulator_pipeline_JULES/impact_regions/impact_region_regression.gpkg"
const BASELINE_FOREST_PATH     = "/Users/dr.bastien/Documents/GitHub/cil-impact-region-ssp-downscaling/data/processed/impact_region_spatial_weights.gpkg"
const POPULATION_PATTERN_PATH  = "/Users/dr.bastien/Documents/GitHub/cil-impact-region-ssp-downscaling/data/processed/impact_region_spatial_weights.gpkg"
const GDP_PATTERN_PATH         = "/Users/dr.bastien/Documents/GitHub/cil-impact-region-ssp-downscaling/data/processed/impact_region_spatial_weights.gpkg"

# The impact_region -> give_country crosswalk.  If you do not have one, leave this
# as "" (empty) and run scripts/build_impact_region_country_crosswalk.py, which
# derives it by spatially joining the impact-region polygons to a world map and
# writes it to `PROCESSED_DIR/impact_region_country_crosswalk.csv`.
const IMPACT_REGION_COUNTRY_CROSSWALK_PATH = ""

# -----------------------------------------------------------------------------
# 2. PROCESSED (compact) DATA DIRECTORY  --  what Julia actually reads
# -----------------------------------------------------------------------------
# The preprocessing scripts write standardized CSVs here; the Julia loader in
# build_forest_give_model.jl reads them.  make_synthetic_fixtures.jl also writes
# here so the module is runnable/testable before real data is plugged in.
const PROCESSED_DIR = normpath(joinpath(@__DIR__, "data", "processed"))

# -----------------------------------------------------------------------------
# 3. OUTPUT DIRECTORY
# -----------------------------------------------------------------------------
const OUTPUT_DIRECTORY = normpath(joinpath(@__DIR__, "outputs"))

# -----------------------------------------------------------------------------
# 4. SCIENTIFIC / UNIT CONFIGURATION
# -----------------------------------------------------------------------------

# Units produced by the fitted forest-response coefficients BEFORE conversion to
# area.  Confirmed by user = :percent  ("5" means 5%).
#   :percent      ->  A_t = A_2022 * (1 + change/100)
#   :proportion   ->  A_t = A_2022 * (1 + change)
#   :log_change   ->  A_t = A_2022 * exp(change)
const FOREST_CHANGE_UNITS = :percent

# Integer code passed to the Mimi component (Mimi parameters are numeric).
forest_change_units_code(u::Symbol) =
    u === :percent    ? 1 :
    u === :proportion ? 2 :
    u === :log_change ? 3 :
    error("Unknown FOREST_CHANGE_UNITS = $u. Use :percent, :proportion, or :log_change.")

# Units of the baseline forest-area column.  Confirmed by user = :Mha
# (millions of hectares).  This is a *label* used for documentation and output
# metadata; the placeholder valuation is unit-agnostic, but keeping it explicit
# prevents silently mixing ha / km^2 / Mha.
const FOREST_AREA_UNITS = :Mha

# The year the forest coefficients were centred on (delta T = T_t - T_baseline).
const GMST_BASELINE_YEAR = 2022

# Behaviour if GMST_BASELINE_YEAR is not an explicit model timestep.
#   :error       -> stop with an informative error (default / safest)
#   :nearest     -> use nearest available year
#   :interpolate -> linear interpolation between neighbours
const GMST_BASELINE_FALLBACK = :error

# -----------------------------------------------------------------------------
# 5. PHYSICAL CONSTRAINTS
# -----------------------------------------------------------------------------
const APPLY_LOWER_CLIP = true    # forest area cannot be negative (A >= 0)
const APPLY_UPPER_CLIP = false   # cap at total impact-region area (only if supplied)

# -----------------------------------------------------------------------------
# 6. PLACEHOLDER ECOSYSTEM-SERVICE VALUATION
# -----------------------------------------------------------------------------
# ESV = ES_VALUE_SCALE * forest_area * gdp * population   (PLACEHOLDER ONLY)
# The resulting units depend entirely on the input units and have NO welfare
# interpretation until a defensible valuation function replaces the placeholder.
const ES_VALUE_SCALE = 1.0
const ES_VALUE_SPEC  = :placeholder    # only :placeholder is implemented for now
es_value_spec_code(s::Symbol) = s === :placeholder ? 1 :
    error("Unknown ES_VALUE_SPEC = $s. Only :placeholder is implemented.")

# -----------------------------------------------------------------------------
# 7. SCC INTEGRATION SWITCH
# -----------------------------------------------------------------------------
# When false (default) we compute and export forest marginal damages and a forest
# SCC, but DO NOT alter GIVE's official total SCC.  When true, forest damages are
# added to the total (post-hoc, non-destructively) AND a prominent warning is
# printed because the placeholder valuation is not scientifically calibrated.
const INCLUDE_FOREST_IN_TOTAL_SCC = false

# Conversion factor applied to placeholder es_damage BEFORE adding to GIVE's
# dollar-denominated damages, used only when INCLUDE_FOREST_IN_TOTAL_SCC = true.
# There is no defensible value for the placeholder; it exists so the plumbing is
# explicit.  Must be set deliberately by the user.
const FOREST_DAMAGE_TO_USD2005 = 1.0

end # module ForestConfig
