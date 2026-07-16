using Mimi

# =============================================================================
# Component A:  SpatialSocioeconomics
# -----------------------------------------------------------------------------
# Downscales GIVE country-level population and GDP to impact regions using
# externally supplied, per-country spatial pattern-scaling shares.
#
#   POP_{r,t} = POP_{c,t} * w^pop_{r,c}
#   GDP_{r,t} = GDP_{c,t} * w^gdp_{r,c}
#   GDPpc_{r,t} = GDP_{r,t} / POP_{r,t}     (never downscale GDP-per-capita directly)
#
# where r is an impact region and c = impact_region_country_index[r] is the GIVE
# country that contains it.
#
# Units follow GIVE exactly:
#   population : millions of persons          (GIVE :Socioeconomic => :population)
#   gdp        : billion US$2005 / yr         (GIVE :Socioeconomic => :gdp)
#   gdppc      : US$2005 / yr / person
#
# R analogy: this is a per-country `group_by(country) %>% mutate(pop_ir = pop_c *
# share)` --- but expressed as a position-indexed array operation so it is cheap
# to run every timestep.
# =============================================================================

@defcomp SpatialSocioeconomics begin

    country       = Index()
    impact_region = Index()

    # ---- Inputs (connected to GIVE :Socioeconomic) --------------------------
    population_country = Parameter(index=[time, country], unit="million")
    gdp_country        = Parameter(index=[time, country], unit="billion US\$2005/yr")

    # ---- Inputs (set once at build time from preprocessed shares) -----------
    population_share = Parameter(index=[impact_region])          # w^pop_{r,c}
    gdp_share        = Parameter(index=[impact_region])          # w^gdp_{r,c}
    # Which GIVE country (1-based index into the :country dimension) each impact
    # region belongs to.  Every impact region maps to exactly one country.
    impact_region_country_index = Parameter{Int}(index=[impact_region])

    # ---- Outputs ------------------------------------------------------------
    population_ir = Variable(index=[time, impact_region], unit="million")
    gdp_ir        = Variable(index=[time, impact_region], unit="billion US\$2005/yr")
    gdppc_ir      = Variable(index=[time, impact_region], unit="US\$2005/yr/person")

    function run_timestep(p, v, d, t)
        for ir in d.impact_region
            c = p.impact_region_country_index[ir]

            pop = p.population_country[t, c] * p.population_share[ir]
            gdp = p.gdp_country[t, c]        * p.gdp_share[ir]

            v.population_ir[t, ir] = pop
            v.gdp_ir[t, ir]        = gdp

            # billion$ / million persons * 1e3 = $/person  (matches GIVE PerCapitaGDP)
            v.gdppc_ir[t, ir] = pop > 0 ? gdp / pop * 1e3 : NaN
        end
    end
end
