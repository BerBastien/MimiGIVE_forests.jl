using Mimi

# =============================================================================
# Component C:  ForestEcosystemServices
# -----------------------------------------------------------------------------
# Computes a PLACEHOLDER ecosystem-service value (ESV) at the impact-region level,
# a no-forest-damage reference value, the forest ecosystem-service damage, and the
# aggregation of all three from impact regions up to GIVE countries and the globe.
#
#   ESV^reference_{r,t} = f(A_{r,2022}, GDP_{r,t}, POP_{r,t})   (forest held at baseline)
#   ESV^climate_{r,t}   = f(A_{r,t},    GDP_{r,t}, POP_{r,t})   (forest responds to climate)
#   Damage^forest_{r,t} = ESV^reference_{r,t} - ESV^climate_{r,t}
#       (+) => loss vs. maintaining baseline forest area; (-) => gain
#
# Country aggregation:  sum over impact regions r within country c.
# Global aggregation:   sum over all countries.
#
# ############################################################################
# ##  SCIENTIFIC LIMITATION                                                  ##
# ##  f(A, GDP, POP) = scale * A * GDP * POP  is a PLACEHOLDER only.          ##
# ##  Its units depend on the input units and it has NO welfare              ##
# ##  interpretation. It exists to exercise the spatial + IAM plumbing and   ##
# ##  is designed to be swapped out (see `ecosystem_service_value` below)    ##
# ##  without touching any other part of the pipeline.                       ##
# ############################################################################
# =============================================================================

"""
    ecosystem_service_value(forest_area, gdp, population; scale=1.0, spec_code=1)

Isolated valuation kernel. Replace the body of the `spec_code == 1` branch (or add
a new spec) with a defensible ecosystem-service valuation function WITHOUT changing
the temperature response, the socioeconomic downscaling, the impact-region
aggregation, the country outputs, or the SCC machinery.

`spec_code`: 1 = :placeholder (the only implemented specification).
"""
function ecosystem_service_value(forest_area, gdp, population; scale::Float64 = 1.0, spec_code::Int = 1)
    if spec_code == 1
        # PLACEHOLDER. Not a scientifically valid valuation.
        return scale * forest_area * gdp * population
    else
        error("Unknown ecosystem-service valuation specification code $spec_code.")
    end
end

@defcomp ForestEcosystemServices begin

    country       = Index()
    impact_region = Index()

    # ---- Inputs -------------------------------------------------------------
    baseline_forest_area  = Parameter(index=[impact_region])
    projected_forest_area = Parameter(index=[time, impact_region])   # from ForestAreaResponse
    population_ir         = Parameter(index=[time, impact_region], unit="million")
    gdp_ir                = Parameter(index=[time, impact_region], unit="billion US\$2005/yr")
    impact_region_country_index = Parameter{Int}(index=[impact_region])

    es_value_scale     = Parameter(default=1.0)
    es_value_spec_code = Parameter{Int}(default=1)

    # ---- Outputs: impact region --------------------------------------------
    es_value_reference_ir = Variable(index=[time, impact_region])
    es_value_climate_ir   = Variable(index=[time, impact_region])
    es_damage_ir          = Variable(index=[time, impact_region])

    # ---- Outputs: GIVE country ---------------------------------------------
    es_value_reference_country = Variable(index=[time, country])
    es_value_climate_country   = Variable(index=[time, country])
    es_damage_country          = Variable(index=[time, country])

    # ---- Outputs: global ----------------------------------------------------
    es_value_reference_global = Variable(index=[time])
    es_value_climate_global   = Variable(index=[time])
    es_damage_global          = Variable(index=[time])

    function run_timestep(p, v, d, t)
        scale = p.es_value_scale
        spec  = p.es_value_spec_code

        # reset country accumulators for this timestep
        for c in d.country
            v.es_value_reference_country[t, c] = 0.0
            v.es_value_climate_country[t, c]   = 0.0
            v.es_damage_country[t, c]          = 0.0
        end

        ref_global = 0.0
        clim_global = 0.0

        for ir in d.impact_region
            c = p.impact_region_country_index[ir]

            ref = ecosystem_service_value(p.baseline_forest_area[ir],
                                          p.gdp_ir[t, ir], p.population_ir[t, ir];
                                          scale = scale, spec_code = spec)
            clim = ecosystem_service_value(p.projected_forest_area[t, ir],
                                           p.gdp_ir[t, ir], p.population_ir[t, ir];
                                           scale = scale, spec_code = spec)
            dmg = ref - clim

            v.es_value_reference_ir[t, ir] = ref
            v.es_value_climate_ir[t, ir]   = clim
            v.es_damage_ir[t, ir]          = dmg

            v.es_value_reference_country[t, c] += ref
            v.es_value_climate_country[t, c]   += clim
            v.es_damage_country[t, c]          += dmg

            ref_global  += ref
            clim_global += clim
        end

        v.es_value_reference_global[t] = ref_global
        v.es_value_climate_global[t]   = clim_global
        v.es_damage_global[t]          = ref_global - clim_global
    end
end
