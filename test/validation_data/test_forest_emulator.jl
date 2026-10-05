@testitem "TestForestEmulator" setup=[DataDepsSetup] begin

    using DataFrames, Statistics

    import MimiGIVE: get_model

    ##------------------------------------------------------------------------------
    ## Configuración
    ##------------------------------------------------------------------------------
    # Eje de tiempo de GIVE: 1750-2300. El emulador solo corre en 2022-2100,
    # el resto de filas queda en `missing` (esperado).

    years = collect(1750:2300)
    ti    = findall(y -> 2022 <= y <= 2100, years)   # filas 273:351
    t22   = findfirst(==(2022), years)
    t100  = findfirst(==(2100), years)
    svc   = [:rec, :hab, :nwfp, :wat]                # orden de mv_cols en get_model

    has_forest(m) = try
        m[:ForestEmulator, :global_forest_Mha]; true
    catch
        false
    end

    ##------------------------------------------------------------------------------
    ## 1. El componente existe solo cuando se pide
    ##------------------------------------------------------------------------------

    m_0 = get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245")
    m   = get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245",
                    forest_emulator = :jules)
    run(m_0)
    run(m)

    @test !has_forest(m_0)
    @test has_forest(m)

    ##------------------------------------------------------------------------------
    ## 2. Salidas: dimensiones, sin missing/NaN en 2022-2100
    ##------------------------------------------------------------------------------

    rc    = m[:ForestEmulator, :region_country]          # país de cada IR (0 = sin país)
    iv    = m[:ForestEmulator, :ir_value]                # [time, IR, servicio]
    iF    = m[:ForestEmulator, :ir_forest_Mha]           # [time, IR]
    iPop  = m[:ForestEmulator, :ir_pop]
    cv    = m[:ForestEmulator, :country_value]           # [time, país, servicio]
    cvn   = m[:ForestEmulator, :country_value_nocl]
    cF    = m[:ForestEmulator, :country_forest_Mha]
    cflow = m[:ForestEmulator, :country_co2_flow]

    nreg, ncou = size(iF, 2), size(cv, 2)

    @test size(iv) == (length(years), nreg, 4)
    @test length(rc) == nreg

    for x in (iv, iF, iPop, cv, cvn, cF)
        @test !any(ismissing, x[ti, ntuple(_ -> :, ndims(x) - 1)...])
        @test !any(v -> isnan(v), skipmissing(x[ti, ntuple(_ -> :, ndims(x) - 1)...]))
    end
    @test all(skipmissing(iv[ti, :, :]) .>= 0)

    ##------------------------------------------------------------------------------
    ## 3. Bosque: 0 <= F <= A_Mha, y en 2022 F = F0 (x = 0)
    ##------------------------------------------------------------------------------

    A   = Float64.(m[:ForestEmulator, :A_Mha])
    F0  = Float64.(m[:ForestEmulator, :F0])
    b1  = Float64.(m[:ForestEmulator, :b1])
    b2  = Float64.(m[:ForestEmulator, :b2])
    mv  = Float64.(m[:ForestEmulator, :mv])
    eg  = Float64.(m[:ForestEmulator, :e_gdp])
    ep  = Float64.(m[:ForestEmulator, :e_pop])

    F = Float64.(iF[ti, :])
    @test all(F .>= 0)
    @test all(F .<= reshape(A, 1, :) .+ 1e-12)
    @test isapprox(Float64.(iF[t22, :]), min.(A, max.(0.0, F0)); rtol = 1e-9, atol = 1e-12)

    # En 2022 no hay cambio de carbono (flujo = 0) y el valor sin clima = valor con clima
    @test all(isapprox.(Float64.(cflow[t22, :]), 0.0; atol = 1e-12))
    @test isapprox(Float64.(cv[t22, :, :]), Float64.(cvn[t22, :, :]); rtol = 1e-9, atol = 1e-12)

    ##------------------------------------------------------------------------------
    ## 4. Coherencia IR vs país
    ##------------------------------------------------------------------------------

    pop = m[:ForestEmulator, :population]
    ps  = Float64.(m[:ForestEmulator, :pop_share])

    for t in ti
        sumF = zeros(ncou)
        sumV = zeros(ncou, 4)
        sumP = zeros(ncou)
        for r in 1:nreg
            c = rc[r]
            c > 0 || continue
            sumF[c] += iF[t, r]
            sumP[c] += iPop[t, r]
            for s in 1:4
                sumV[c, s] += iv[t, r, s]
            end
        end
        @test isapprox(sumF, Float64.(cF[t, :]);   rtol = 1e-9, atol = 1e-9)
        @test isapprox(sumV, Float64.(cv[t, :, :]); rtol = 1e-9, atol = 1e-9)

        # población: Σ pop IR = pop país × Σ pop_share de sus IR
        sumS = zeros(ncou)
        for r in 1:nreg
            rc[r] > 0 && (sumS[rc[r]] += ps[r])
        end
        popt = Float64.(coalesce.(pop[t, :], 0.0))
        @test isapprox(sumP, popt .* sumS; rtol = 1e-9, atol = 1e-9)
    end

    ##------------------------------------------------------------------------------
    ## 5. Trayectoria de temperatura -> hectáreas -> valor de los 4 servicios
    ##    Se recalcula a mano para unas pocas IR y se compara con el componente.
    ##------------------------------------------------------------------------------

    T    = Float64.(coalesce.(m[:ForestEmulator, :temperature], NaN))
    pcg  = m[:ForestEmulator, :pc_gdp]
    Tref = T[t22]

    # IR con país y bosque en 2022: 5 repartidas en la lista
    cand   = findall(r -> rc[r] > 0 && iF[t22, r] > 0 &&
                          pop[t22, rc[r]] > 0 && pcg[t22, rc[r]] > 0, 1:nreg)
    @test !isempty(cand)
    sample = cand[unique(round.(Int, range(1, length(cand), length = min(5, length(cand)))))]

    for r in sample, t in ti
        x     = T[t] - Tref
        F_exp = min(A[r], max(0.0, F0[r] + b1[r] * x + b2[r] * x^2))
        @test isapprox(Float64(iF[t, r]), F_exp; rtol = 1e-9, atol = 1e-12)

        c = rc[r]
        for s in 1:4
            pr  = pop[t, c] / pop[t22, c]
            gr  = pcg[t, c] / pcg[t22, c]
            fac = (isfinite(pr) && isfinite(gr)) ? gr^eg[s] * pr^ep[s] : 0.0
            @test isapprox(Float64(iv[t, r, s]), mv[r, s] * fac * F_exp; rtol = 1e-8, atol = 1e-12)
        end
    end

    # En 2022 el valor por hectárea es el valor unitario mv (fac = 1)
    for r in sample, s in 1:4
        @test isapprox(Float64(iv[t22, r, s]) / Float64(iF[t22, r]), mv[r, s]; rtol = 1e-8)
    end

    # Vistazo humano: temperatura, hectáreas y 4 servicios para una IR (se imprime al correr)
    r0   = sample[1]
    show_years = [2022, 2040, 2060, 2080, 2100]
    tv   = [findfirst(==(y), years) for y in show_years]
    tab  = DataFrame(year = show_years,
                     T    = T[tv],
                     x    = T[tv] .- Tref,
                     F_Mha = Float64.(iF[tv, r0]),
                     rec  = Float64.(iv[tv, r0, 1]),
                     hab  = Float64.(iv[tv, r0, 2]),
                     nwfp = Float64.(iv[tv, r0, 3]),
                     wat  = Float64.(iv[tv, r0, 4]))
    @info "IR #$r0 (país #$(rc[r0])), valores en 1e6 USD/año" tab

    ##------------------------------------------------------------------------------
    ## 6. Comparación entre trials (opcional: requiere un `res` de Monte Carlo)
    ##------------------------------------------------------------------------------
    # Con SSP la socioeconomía es la misma en todos los trials, así que:
    #   - el valor por hectárea (country_value / country_forest_Mha) debe ser igual,
    #   - la temperatura y las hectáreas deben diferir.
    # Uso: pasa el resultado de tu corrida (n = 2) con save_list que incluya
    #   (:temperature, :T), (:ForestEmulator, :country_value),
    #   (:ForestEmulator, :country_forest_Mha)

    function compare_trials(res; year = 2100)
        dv = dropmissing(getdataframe(res, :ForestEmulator, :country_value))
        dF = dropmissing(getdataframe(res, :ForestEmulator, :country_forest_Mha))
        dT = getdataframe(res, :temperature, :T)

        dv = DataFrame(country = String.(dv.country), services = dv.services,
                       time = Int.(dv.time), trialnum = dv.trialnum,
                       value = Float64.(dv.country_value))
        dF = DataFrame(country = String.(dF.country), time = Int.(dF.time),
                       trialnum = dF.trialnum, F = Float64.(dF.country_forest_Mha))

        dj = innerjoin(dv, dF, on = [:country, :time, :trialnum])
        filter!(r -> r.time == year && r.F > 0, dj)
        dj.usd_ha = dj.value ./ dj.F

        # socioeconomía igual entre trials -> diferencia en USD/ha ≈ 0
        spread_ha = combine(groupby(dj, [:country, :services]),
                            :usd_ha => (x -> maximum(x) - minimum(x)) => :spread)
        @test maximum(spread_ha.spread) < 1e-6 * maximum(dj.usd_ha)

        # la temperatura sí cambia entre trials
        Ty = filter(r -> r.time == year, dT)
        @test std(Float64.(Ty.T)) > 0
        return dj
    end

    # Descomenta y ajusta a tu wrapper de Monte Carlo:
    # res = <tu llamada con n = 2 y save_list>
    # compare_trials(res)

end