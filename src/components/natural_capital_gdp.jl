# GDP adjusted by natural capital.
#
# This component is the single GDP source of the model whenever the forest
# emulator is on. Every component that consumes GDP (PerCapitaGDP, VSL, the
# agriculture aggregators, energy damages, DICE2016R2, Howard-Sterner, the
# damage aggregator and the net consumption components) reads `gdp` from here,
# so switching `use_nc` moves all of them between the standard GIVE GDP and the
# natural-capital-adjusted GDP at once.
#
#   gdp(t, c) = gdp_ssp(t, c) * (S_hab(t) / S_ref(t)) ^ b(c)
#
#   S_hab  global value of habitat services in the current run (ForestEmulator)
#   S_ref  the same variable in the run without a pulse, i.e. the counterfactual
#          without climate-driven ecosystem service loss (set_natural_capital!)
#   b      elasticity of GDP to natural capital, per country
#
# When the switch is off, or when either value is missing (2020-2021, before
# the first year of ForestEmulator), the ratio is set to 1 and GDP is exactly
# the Socioeconomic GDP.
@defcomp NaturalCapitalGDP begin
    country = Index()

    gdp_ssp = Parameter(index=[time, country], unit="billion US\$2005/yr")
    S_hab   = Parameter(index=[time])          # global habitat value, ForestEmulator
    S_ref   = Parameter(index=[time])          # S_hab of the base run (no pulse)
    b       = Parameter(index=[country])       # elasticity of GDP to natural capital
    use_nc  = Parameter{Bool}(default = false) # on/off switch

    gdp = Variable(index=[time, country], unit="billion US\$2005/yr")

    function run_timestep(p, v, d, t)
        active = p.use_nc && p.S_ref[t] > 0 && p.S_hab[t] > 0
        ratio  = active ? p.S_hab[t] / p.S_ref[t] : 1.0
        for c in d.country
            v.gdp[t, c] = p.gdp_ssp[t, c] * ratio^p.b[c]
        end
    end
end
