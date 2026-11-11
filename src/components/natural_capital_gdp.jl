@defcomp NaturalCapitalGDP begin
    country = Index()

    gdp_ssp = Parameter(index=[time, country], unit="billion US\$2005/yr")
    S_hab   = Parameter(index=[time])          # global_hab_value de ForestEmulator
    S_ref   = Parameter(index=[time])          # S_hab de una corrida base sin pulso
    b       = Parameter(index=[country])       # elasticidad del PIB a S (CSV del paper)
    use_nc  = Parameter{Bool}(default = false) # interruptor

    gdp = Variable(index=[time, country], unit="billion US\$2005/yr")

    function run_timestep(p, v, d, t)
        active = p.use_nc && p.S_ref[t] > 0 && p.S_hab[t] > 0
        ratio  = active ? p.S_hab[t] / p.S_ref[t] : 1.0
        for c in d.country
            v.gdp[t, c] = p.gdp_ssp[t, c] * ratio^p.b[c]
        end
    end
end