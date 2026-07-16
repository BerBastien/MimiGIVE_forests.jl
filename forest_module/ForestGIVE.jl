# =============================================================================
# ForestGIVE.jl  --  entry point for the MimiGIVE forest ecosystem-services module
# -----------------------------------------------------------------------------
# This is a NON-DESTRUCTIVE extension of MimiGIVE. It does not modify any file in
# MimiGIVE's own src/. To use it:
#
#   include("forest_module/ForestGIVE.jl")
#   using .ForestGIVE
#
# from a Julia session whose active project has MimiGIVE (and its deps) installed.
#
# R analogy: think of this file as an R package's top-level file that `source()`s
# every function definition and re-exports the public API.
# =============================================================================

module ForestGIVE

using Mimi
using MimiGIVE
using DataFrames, CSVFiles, FileIO
using Statistics

# Configuration (paths + knobs) as a submodule.
include("config.jl")

# The three Mimi components.
include("src/forest_spatial_socioeconomics.jl")
include("src/forest_area_response.jl")
include("src/forest_ecosystem_services.jl")

# Data loading, model construction, SCC, and outputs.
include("src/forest_inputs.jl")
include("src/build_forest_give_model.jl")
include("src/forest_scc.jl")
include("src/forest_outputs.jl")

export
    # components
    SpatialSocioeconomics, ForestAreaResponse, ForestEcosystemServices,
    ecosystem_service_value,
    # data + model
    ForestInputs, load_forest_inputs, build_forest_give_model,
    # scc
    compute_forest_scc, compute_total_scc_with_forest,
    # outputs
    write_impact_region_output, write_country_output, write_marginal_country_output,
    write_summary_output, write_temperature_diagnostic, write_clipping_diagnostic,
    write_join_diagnostics,
    # config submodule
    ForestConfig

end # module ForestGIVE
