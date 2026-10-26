@testitem "TestComputeUtilityFunction" setup=[DataDepsSetup] begin
    using MimiGIVE
    using Mimi
    using Test

    const YEAR  = 2020
    const ETA   = 1.45
    const PRTP  = 0.015
    const SSPS  = ["SSP119", "SSP245", "SSP585"]

    # Central CES parameters used in the tests and in the scenario table.
    const S     = 0.10
    const THETA = 0.55

    const PAR_MOORE = (forest = true, s = S, θ = THETA)
    const PAR_DERIV = (forest = true, s = S, θ = THETA,
                       s2 = 0.0, θ2 = 0.5, i_pc = 1.0)

    function model_for(scenario::String; forest::Bool = false)
        m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = scenario,
                               forest_emulator = forest ? :jules : nothing)
        update_param!(m, :DamageAggregator, :include_slr, false)
        return m
    end

    @testset "Regresión" begin
        m = model_for("SSP245")
        base = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP)

        # Moore with s = 0 and no forests must reproduce the original SCC.
        @test MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
            method = :moore,
            params = (forest = false, s = 0.0, θ = THETA)) ≈ base

        # Derivatives with s = s2 = 0 and no forests must reproduce the original SCC.
        @test MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
            method = :derivatives,
            params = (forest = false, s = 0.0, θ = THETA,
                      s2 = 0.0, θ2 = 0.5, i_pc = 1.0)) ≈ base
    end

    @testset "Bosques" begin
        m = model_for("SSP245"; forest = true)
        base = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP)

        # s = 0 switches the forest channel off and must reproduce the original.
        mo0 = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
            method = :moore, params = (forest = true, s = 0.0, θ = THETA))
        @test mo0 ≈ base

        # With s > 0 the SCC must be finite and include the ecosystem service cost.
        moore = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
            method = :moore, params = PAR_MOORE)
        @test isfinite(moore)
        @test moore > base

        # With s2 = 0 the nested CES collapses to U, so both methods must agree.
        deriv = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
            method = :derivatives, params = PAR_DERIV)
        @test isfinite(deriv)
        @test moore ≈ deriv

        # A negative θ (elasticity of substitution > 1) must remain finite.
        neg = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
            method = :moore, params = (forest = true, s = S, θ = -0.16))
        @test isfinite(neg)
    end

    println("\nSCC (USD2005/tCO2), año $YEAR")
    println(rpad("escenario", 10), rpad("original", 14), rpad("moore", 14), "derivatives")

    for sc in SSPS
        m = model_for(sc; forest = true)
        o  = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP)
        mo = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
                                  method = :moore, params = PAR_MOORE)
        de = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
                                  method = :derivatives, params = PAR_DERIV)
        println(rpad(sc, 10), rpad(round(o, digits = 2), 14),
                rpad(round(mo, digits = 2), 14), round(de, digits = 2))
    end
end
