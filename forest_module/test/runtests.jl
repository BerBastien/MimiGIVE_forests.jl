# Forest module test suite.
#
# Fast component + marginal-plumbing tests need only Mimi:
#     julia --project=. forest_module/test/runtests.jl
#
# To ALSO run the heavy full-MimiGIVE integration test:
#     FOREST_RUN_FULL=true julia --project=. forest_module/test/runtests.jl

using Test

@testset "ForestGIVE module" begin
    include("test_spatial_socioeconomics.jl")
    include("test_forest_area_response.jl")
    include("test_forest_ecosystem_services.jl")
    include("test_forest_give_integration.jl")
end
