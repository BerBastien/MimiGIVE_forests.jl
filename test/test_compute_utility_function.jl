using MimiGIVE, Mimi, Test

const YEAR  = 2020
const ETA   = 1.45
const PRTP  = 0.015
const RHO   = log(1 + PRTP)     # equivalente continuo exacto de PRTP
const SSPS  = ["SSP119", "SSP245", "SSP585"]

# PLACEHOLDER: reemplazar cuando tengas s y θ confirmados
const S     = 0.3
const THETA = 0.5

const PAR_MOORE = (s = S, θ = THETA)
const PAR_DERIV = (s = S, θ = THETA, s2 = 0.0, θ2 = 0.5, α = ETA, ρ = RHO, i_pc = 1.0)

function model_for(scenario::String)
    m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = scenario)
    update_param!(m, :DamageAggregator, :include_slr, false)
    return m
end

@testset "Regresión" begin
    m = model_for("SSP245")
    base = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP)

    # Moore con s = 0 y sin bosques debe reproducir el original
    @test MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP, method = :moore,
        params = (s = 0.0, θ = THETA, forest = false)) ≈ base

    # Derivadas con s = s2 = 0 y sin bosques debe reproducir el original
    @test MimiGIVE.compute_scc(m, year = YEAR, method = :derivatives,
        params = (s = 0.0, θ = THETA, s2 = 0.0, θ2 = 0.5, α = ETA, ρ = RHO,
                  i_pc = 1.0, forest = false)) ≈ base
end

println("\nSCC (USD2005/tCO2), año $YEAR")
println(rpad("escenario", 10), rpad("original", 14), rpad("moore", 14), "derivatives")

for sc in SSPS
    m = model_for(sc)
    o = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP)
    mo = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
                              method = :moore, params = PAR_MOORE)
    de = MimiGIVE.compute_scc(m, year = YEAR, method = :derivatives, params = PAR_DERIV)
    println(rpad(sc, 10), rpad(round(o, digits = 2), 14),
            rpad(round(mo, digits = 2), 14), round(de, digits = 2))
end