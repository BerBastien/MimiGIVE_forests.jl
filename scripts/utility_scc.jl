# ------------------------------------------------------------------------------
# CES utility SCC with forest ecosystem services
#
# Reproduces the numbers and figures of docs/utility_function.tex:
#   - SCC 2020 by scenario: original, :moore and :derivatives
#   - sensitivity of the :moore SCC to s and θ (SSP245)
#   - time paths of C, E, U, marginal utilities and service values
#
# Usage:
#   julia --project=. scripts/utility_scc.jl
#
# Outputs:
#   docs/utility_scc_results.csv
#   docs/utility_scc_sensitivity.csv
#   docs/figures/utility_paths.png / .pdf
#   docs/figures/scc_scenarios.png / .pdf
#   docs/figures/scc_sensitivity.png / .pdf
#
# Set UTILITY_SCC_REUSE=true to reuse the CSVs instead of recomputing the SCCs
# (useful when only the figures need to be regenerated).
# ------------------------------------------------------------------------------

using MimiGIVE, Mimi, CSV, DataFrames, CairoMakie, Statistics

const ROOT    = normpath(joinpath(@__DIR__, ".."))
const DOCS    = joinpath(ROOT, "docs")
const FIGDIR  = joinpath(DOCS, "figures")
const RESULT  = joinpath(DOCS, "utility_scc_results.csv")
const SENS    = joinpath(DOCS, "utility_scc_sensitivity.csv")
const REUSE   = get(ENV, "UTILITY_SCC_REUSE", "false") == "true"

const YEAR    = 2020
const ETA     = 1.45
const PRTP    = 0.015
const SCENARIOS = ["SSP119", "SSP245", "SSP585"]
const S_CENTRAL = 0.10
const THETA_CENTRAL = 0.55

# Services in the order used by ForestEmulator: rec, hab, nwfp, wat.
const SERVICE_LABELS = ["Recreación", "Hábitat", "Prod. no maderables", "Agua"]
const SERVICE_COLORS = ["#0072B2", "#009E73", "#CC79A7", "#D55E00"]
const NONMARKET = [1, 2, 4]

mkpath(FIGDIR)

function model_for(scenario; forest = true, include_slr = true)
    m = MimiGIVE.get_model(socioeconomics_source = :SSP,
                           SSP_scenario = scenario,
                           forest_emulator = forest ? :jules : nothing)
    update_param!(m, :DamageAggregator, :include_slr, include_slr)
    return m
end

moore_params(s = S_CENTRAL, θ = THETA_CENTRAL) = (forest = true, s = s, θ = θ)

deriv_params(s = S_CENTRAL, θ = THETA_CENTRAL) =
    (forest = true, s = s, θ = θ, s2 = 0.10, θ2 = 0.58, i_pc = 1.0)

# ------------------------------------------------------------------------------
# 1. SCC by scenario and method

function compute_scc_table()
    rows = NamedTuple[]
    for sc in SCENARIOS
        @info "SCC table: $sc"
        m = model_for(sc)
        base = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP)
        mo = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
                                  method = :moore, params = moore_params())
        de = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
                                  method = :derivatives, params = deriv_params())
        for (label, value) in (("original", base), ("moore", mo), ("derivatives", de))
            push!(rows, (scenario = sc, method = label, scc = value,
                         pct_change = 100 * (value - base) / base))
        end
    end
    return DataFrame(rows)
end

# ------------------------------------------------------------------------------
# 2. Sensitivity of the :moore SCC to s and θ (SSP245)

function compute_sensitivity()
    @info "SCC sensitivity: SSP245"
    m = model_for("SSP245")
    base = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP)
    rows = NamedTuple[]
    for s in (0.05, 0.10, 0.15), θ in (-0.16, 0.55, 0.86)
        v = MimiGIVE.compute_scc(m, year = YEAR, eta = ETA, prtp = PRTP,
                                 method = :moore, params = moore_params(s, θ))
        push!(rows, (scenario = "SSP245", s = s, θ = θ, scc = v,
                     pct_change = 100 * (v - base) / base))
    end
    return DataFrame(rows)
end

# ------------------------------------------------------------------------------
# 3. Time paths of the CES formula (SSP245, central parameters)

function compute_paths(; scenario = "SSP245", s = S_CENTRAL, θ = THETA_CENTRAL)
    @info "Utility paths: $scenario"
    m = model_for(scenario)
    mm = MimiGIVE.get_marginal_model(deepcopy(m), year = YEAR, gas = :CO2, pulse_size = 1.0)
    run(mm)

    params = moore_params(s, θ)
    E, _ = MimiGIVE._forest_components(mm, :CO2, params)

    i0 = findfirst(==(2020), MimiGIVE._model_years)
    idx = i0:lastindex(MimiGIVE._model_years)
    years = MimiGIVE._model_years[idx]

    C = Float64.(coalesce.(mm.base[:global_netconsumption, :net_cpc][idx], NaN))
    E = E[idx]

    U  = ((1 - s) .* C .^ θ .+ s .* E .^ θ) .^ (1 / θ)
    μC = (1 - s) .* C .^ (θ - 1) .* U .^ (1 - θ)
    μE = s .* E .^ (θ - 1) .* U .^ (1 - θ)

    share_C = μC .* C ./ U
    share_E = μE .* E ./ U

    # Global value of each service, in billion USD2005/yr, 2022-2100.
    defl = 1 / MimiGIVE.pricelevel_2005_to_2020
    cv = mm.base[:ForestEmulator, :country_value]
    syears = collect(2022:2100)
    sidx = [findfirst(==(y), MimiGIVE._model_years) for y in syears]
    service_values = zeros(length(syears), 4)
    for (k, i) in enumerate(sidx), sv in 1:4
        service_values[k, sv] = sum(skipmissing(cv[i, :, sv])) * 1e6 * defl / 1e9
    end

    return (
        years = years,
        C = C,
        E = E,
        U = U,
        μC = μC,
        μE = μE,
        share_C = share_C,
        share_E = share_E,
        service_years = syears,
        service_values = service_values,
        params = params,
    )
end

# ------------------------------------------------------------------------------
# 4. Figures

function fig_paths(p)
    fig = Figure(size = (1200, 850), fontsize = 13)

    # (a) stacked value of non-market services
    ax1 = Axis(fig[1, 1], title = "(a) Valor de los servicios no comerciales",
               xlabel = "Año", ylabel = "miles de millones USD2005/año")
    vals = p.service_values[:, NONMARKET]
    cum = cumsum(vals, dims = 2)
    lower = zeros(size(vals, 1))
    for (j, sv) in enumerate(NONMARKET)
        upper = cum[:, j]
        band!(ax1, p.service_years, lower, upper,
              color = (SERVICE_COLORS[sv], 0.55), label = SERVICE_LABELS[sv])
        lower = upper
    end
    axislegend(ax1, position = :lt, framevisible = false)

    # (b) E / C
    ax2 = Axis(fig[1, 2], title = "(b) Servicios ecosistémicos / consumo",
               xlabel = "Año", ylabel = "E / C (%)")
    vspan!(ax2, 2022, 2100, color = (:gray, 0.08))
    lines!(ax2, p.years, 100 .* p.E ./ p.C, color = "#0072B2", linewidth = 2.5)
    text!(ax2, 2061, maximum(100 .* p.E ./ p.C) * 0.95,
          text = "periodo del emulador", align = (:center, :top), fontsize = 11)

    # (c) per-capita levels
    ax3 = Axis(fig[2, 1], title = "(c) Niveles per cápita",
               xlabel = "Año", ylabel = "USD2005 / persona / año",
               yscale = log10)
    vspan!(ax3, 2022, 2100, color = (:gray, 0.08))
    lines!(ax3, p.years, p.C, color = "black", linewidth = 2.5, label = "Consumo C")
    lines!(ax3, p.years, p.E, color = "#D55E00", linewidth = 2.5,
           label = "Servicios E (rec+hab+wat)")
    axislegend(ax3, position = :rb, framevisible = false)

    # (d) ecosystem share in the CES composite (consumption is the complement)
    ax4 = Axis(fig[2, 2], title = "(d) Participación de E en la utilidad",
               xlabel = "Año", ylabel = "porcentaje de U")
    vspan!(ax4, 2022, 2100, color = (:gray, 0.08))
    sh = 100 .* p.share_E
    lines!(ax4, p.years, sh, color = "#D55E00", linewidth = 2.5,
           label = "Servicios E")
    ylims!(ax4, 0, max(0.5, 1.15 * maximum(sh)))
    text!(ax4, 2023, 0.93 * max(0.5, 1.15 * maximum(sh)),
          text = "el resto es consumo", align = (:left, :top), fontsize = 10)
    axislegend(ax4, position = :rb, framevisible = false)

    save(joinpath(FIGDIR, "utility_paths.png"), fig, px_per_unit = 2)
    save(joinpath(FIGDIR, "utility_paths.pdf"), fig)
    return fig
end

function fig_scenarios(tab)
    fig = Figure(size = (1200, 520), fontsize = 13)
    methods = ["original", "moore", "derivatives"]
    colors  = ["#4D4D4D", "#0072B2", "#D55E00"]
    xt = (1:length(SCENARIOS), SCENARIOS)

    ax1 = Axis(fig[1, 1], title = "(a) SCC 2020", xticks = xt,
               ylabel = "USD2005 / tCO2")
    offsets = [-0.25, 0.0, 0.25]
    for (k, method) in enumerate(methods)
        sub = tab[tab.method .== method, :]
        xs = collect(1:length(SCENARIOS)) .+ offsets[k]
        barplot!(ax1, xs, sub.scc, width = 0.23, color = colors[k], label = method)
    end
    axislegend(ax1, position = :rt, framevisible = false)

    ax2 = Axis(fig[1, 2], title = "(b) Cambio vs original", xticks = xt,
               ylabel = "%")
    offsets = [-0.13, 0.13]
    for (k, method) in enumerate(["moore", "derivatives"])
        sub = tab[tab.method .== method, :]
        xs = collect(1:length(SCENARIOS)) .+ offsets[k]
        barplot!(ax2, xs, sub.pct_change, width = 0.22,
                 color = colors[k + 1], label = method)
    end
    hlines!(ax2, [0.0], color = (:black, 0.4), linestyle = :dash)
    axislegend(ax2, position = :rt, framevisible = false)

    save(joinpath(FIGDIR, "scc_scenarios.png"), fig, px_per_unit = 2)
    save(joinpath(FIGDIR, "scc_scenarios.pdf"), fig)
    return fig
end

function fig_sensitivity(sens)
    fig = Figure(size = (900, 520), fontsize = 13)
    ax = Axis(fig[1, 1], title = "Sensibilidad del SCC a s y θ (SSP245)",
              xlabel = "θ", ylabel = "cambio del SCC vs original (%)")
    hlines!(ax, [0.0], color = (:black, 0.4), linestyle = :dash)
    colors = ["#1B9E77", "#0072B2", "#D95F02"]
    for (k, s) in enumerate((0.05, 0.10, 0.15))
        sub = sort(sens[sens.s .== s, :], :θ)
        lines!(ax, sub.θ, sub.pct_change, color = colors[k], linewidth = 2.5,
               label = "s = $s")
        scatter!(ax, sub.θ, sub.pct_change, color = colors[k], markersize = 9)
    end
    axislegend(ax, position = :lb, framevisible = false)
    save(joinpath(FIGDIR, "scc_sensitivity.png"), fig, px_per_unit = 2)
    save(joinpath(FIGDIR, "scc_sensitivity.pdf"), fig)
    return fig
end

# ------------------------------------------------------------------------------
# Main

function main()
    if REUSE && isfile(RESULT) && isfile(SENS)
        @info "Reusing $RESULT and $SENS"
        tab = CSV.read(RESULT, DataFrame)
        sens = CSV.read(SENS, DataFrame)
    else
        tab = compute_scc_table()
        sens = compute_sensitivity()
        CSV.write(RESULT, tab)
        CSV.write(SENS, sens)
    end

    println("\nSCC (USD2005/tCO2), año $YEAR, JULES, parámetros centrales")
    println(rpad("escenario", 10), rpad("original", 14), rpad("moore", 14),
            rpad("derivatives", 16), "Δ moore")
    for sc in SCENARIOS
        sub = tab[tab.scenario .== sc, :]
        getv(m) = first(sub.scc[sub.method .== m])
        getp(m) = first(sub.pct_change[sub.method .== m])
        println(rpad(sc, 10), rpad(round(getv("original"), digits = 4), 14),
                rpad(round(getv("moore"), digits = 4), 14),
                rpad(round(getv("derivatives"), digits = 4), 16),
                string(round(getp("moore"), digits = 2), " %"))
    end

    println("\nSensibilidad SSP245 (moore):")
    for row in eachrow(sort(sens, [:s, :θ]))
        println("  s=", row.s, "  θ=", row.θ, "  SCC=",
                round(row.scc, digits = 4), "  Δ=", round(row.pct_change, digits = 2), "%")
    end

    p = compute_paths()
    fig_paths(p)
    fig_scenarios(tab)
    fig_sensitivity(sens)

    println("\nFiguras en $FIGDIR")
    return (table = tab, sensitivity = sens, paths = p)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
