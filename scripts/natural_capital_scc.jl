# ------------------------------------------------------------------------------
# Natural capital GDP and SCC combinations
#
# Usage:
#   julia --project=. scripts/natural_capital_scc.jl
#
# For SSP245 with the JULES forest emulator, this script writes:
#
#   docs/natural_capital_gdp.csv           world GDP per year, normal and with
#                                          natural capital
#   docs/natural_capital_gdp_countries.csv the same by country in 2030, 2050,
#                                          2100, 2200 and 2300
#   docs/natural_capital_scc.csv           SCC 2020 for the 12 combinations of
#                                          SCC formula, forest market damages and
#                                          GDP definition
#   docs/figures/natural_capital/*.png     12 figures, one per combination
#
# The 12 combinations are
#
#   method           = :original, :forest, :forest_nested
#   forest_damages   = false, true        (market forest damages in consumption)
#   natural_capital  = false, true        (standard GDP vs natural capital GDP)
#
# Elasticity of GDP to natural capital b: by macroregion (OECD, LAM, MAF, ASIA,
# REF) from Bastien-Olvera et al. (Nature 2024), see scripts/b_por_region.jl.
# ------------------------------------------------------------------------------

using MimiGIVE, Mimi, CSV, DataFrames, CairoMakie, Printf

include(joinpath(@__DIR__, "b_por_region.jl"))

const ROOT   = normpath(joinpath(@__DIR__, ".."))
const DOCS   = joinpath(ROOT, "docs")
const FIGDIR = joinpath(DOCS, "figures", "natural_capital")

const SCENARIO   = "SSP245"
const EMULATOR   = :jules
const B_LABEL    = "paper_regional"          # etiqueta en las salidas
const B_VEC      = b_vector(ROOT)            # un valor por país, orden de get_model
const PULSE_YEAR = 2020
const ETA  = 1.45
const PRTP = 0.015

const BASE_ARGS = (socioeconomics_source = :SSP, SSP_scenario = SCENARIO,
                   forest_emulator = EMULATOR)

const METHODS = (:original, :forest, :forest_nested)

const FOREST_PARAMS        = (forest = true, s = 0.10, θ = 0.55)
const FOREST_NESTED_PARAMS = (forest = true, s = 0.10, θ = 0.55, s2 = 0.10, θ2 = 0.58, i_pc = 1.0)

params_for(method) = method == :forest ? FOREST_PARAMS :
                     method == :forest_nested ? FOREST_NESTED_PARAMS : nothing

label(x::Bool) = x ? "on" : "off"
gdp_label(x::Bool) = x ? "natural capital" : "standard"

mkpath(FIGDIR)

# ------------------------------------------------------------------------------
# Models: one per (forest damages, natural capital) combination

function build_models()
    models = Dict{Tuple{Bool, Bool}, Mimi.Model}()
    for fd in (false, true), nc in (false, true)
        @info "Building model: forest_damages = $fd, natural_capital = $nc"
        m = MimiGIVE.get_model(; BASE_ARGS...,
                               forest_damages = fd,
                               natural_capital = nc,
                               natural_capital_elasticity = B_VEC)
        run(m)
        models[(fd, nc)] = m
    end
    return models
end

# ------------------------------------------------------------------------------
# 1. GDP: normal vs natural capital

function write_gdp_tables(models)
    m = models[(false, true)]      # forest damages do not change GDP
    years = collect(Mimi.time_labels(m))
    ti = findall(y -> 2020 <= y <= 2300, years)
    model_countries = String.(Mimi.dim_keys(m, :country))

    gdp_ssp = m[:Socioeconomic, :gdp]
    gdp_nc  = m[:NaturalCapitalGDP, :gdp]

    world = DataFrame(year = Int[], gdp_normal = Float64[],
                      gdp_natural_capital = Float64[], ratio = Float64[])
    for i in ti
        g0 = sum(skipmissing(gdp_ssp[i, :]))
        g1 = sum(skipmissing(gdp_nc[i, :]))
        push!(world, (year = years[i], gdp_normal = g0, gdp_natural_capital = g1, ratio = g1 / g0))
    end
    CSV.write(joinpath(DOCS, "natural_capital_gdp.csv"), world)

    show_years = [2030, 2050, 2100, 2200, 2300]
    rows = DataFrame(year = Int[], iso3 = String[], gdp_normal = Float64[],
                     gdp_natural_capital = Float64[])
    for y in show_years
        i = findfirst(==(y), years)
        for (c, iso) in enumerate(model_countries)
            push!(rows, (year = y, iso3 = iso,
                         gdp_normal = coalesce(gdp_ssp[i, c], 0.0),
                         gdp_natural_capital = coalesce(gdp_nc[i, c], 0.0)))
        end
    end
    CSV.write(joinpath(DOCS, "natural_capital_gdp_countries.csv"), rows)

    return world
end

# ------------------------------------------------------------------------------
# 2. SCC for the 12 combinations

function compute_scc_table(models)
    rows = DataFrame(method = String[], forest_damages = Bool[],
                     natural_capital = Bool[], elasticity = String[],
                     scc_usd2005_per_tco2 = Float64[])
    for method in METHODS
        for fd in (false, true), nc in (false, true)
            @info "SCC: $method, forest damages $(label(fd)), GDP $(gdp_label(nc))"
            m = models[(fd, nc)]
            scc = MimiGIVE.compute_scc(m; year = PULSE_YEAR, eta = ETA, prtp = PRTP,
                                       method = method, params = params_for(method))
            push!(rows, (method = String(method), forest_damages = fd,
                         natural_capital = nc,
                         elasticity = nc ? B_LABEL : "not applicable",
                         scc_usd2005_per_tco2 = scc))
        end
    end
    CSV.write(joinpath(DOCS, "natural_capital_scc.csv"), rows)
    return rows
end

# ------------------------------------------------------------------------------
# 3. Marginal damage path of a 2020 CO2 pulse (non-SLR channels). It is the same
#    for the three SCC formulas, which only differ in how it is aggregated.

function marginal_paths(models)
    paths = Dict{Tuple{Bool, Bool}, NamedTuple}()
    for fd in (false, true), nc in (false, true)
        @info "Marginal damages: forest damages $(label(fd)), GDP $(gdp_label(nc))"
        mm = MimiGIVE.get_marginal_model(models[(fd, nc)]; year = PULSE_YEAR,
                                         gas = :CO2, pulse_size = 1.0)
        run(mm)
        conv = MimiGIVE.scc_gas_molecular_conversions[:CO2]
        total  = Float64.(coalesce.(mm[:DamageAggregator, :total_damage], 0.0)) .* conv
        forest = Float64.(coalesce.(mm[:DamageAggregator, :forest_damage], 0.0)) .* conv
        paths[(fd, nc)] = (years = collect(Mimi.time_labels(mm.base)),
                           total = total, forest = forest)
    end
    return paths
end

# ------------------------------------------------------------------------------
# 4. Figures: one per combination

function plot_combo(method, fd, nc, path, scc)
    keep = findall(y -> PULSE_YEAR <= y <= 2300, path.years)
    yrs = path.years[keep]
    total = path.total[keep]
    forest = path.forest[keep]

    fig = Figure(size = (1200, 520), fontsize = 13)

    ax1 = Axis(fig[1, 1],
               title = "(a) Total marginal damages of a 2020 CO2 pulse\n" *
                       "SCC 2020 = " * @sprintf("%.2f", scc) * " USD2005/tCO2",
               xlabel = "Year", ylabel = "USD2005 / tCO2")
    lines!(ax1, yrs, total, color = "#0072B2", linewidth = 2.5)
    hlines!(ax1, [0.0], color = (:black, 0.4), linestyle = :dash)

    ax2 = Axis(fig[1, 2],
               title = fd ? "(b) Forest market damages (in consumption)" :
                           "(b) Forest market damages (not in consumption)",
               xlabel = "Year", ylabel = "USD2005 / tCO2")
    lines!(ax2, yrs, forest, color = "#D55E00", linewidth = 2)
    hlines!(ax2, [0.0], color = (:black, 0.4), linestyle = :dash)

    Label(fig[0, :], "formula: $(method) · forest damages: $(label(fd)) · GDP: $(gdp_label(nc))",
          fontsize = 15, font = :bold)

    fname = "scc_$(method)_forest_$(label(fd))_gdp_$(nc ? "natural_capital" : "standard").png"
    save(joinpath(FIGDIR, fname), fig, px_per_unit = 2)
    return fname
end

# ------------------------------------------------------------------------------
# Main

function main()
    models = build_models()
    world = write_gdp_tables(models)
    sccs = compute_scc_table(models)
    paths = marginal_paths(models)

    files = String[]
    for row in eachrow(sccs)
        method = Symbol(row.method)
        push!(files, plot_combo(method, row.forest_damages, row.natural_capital,
                                paths[(row.forest_damages, row.natural_capital)],
                                row.scc_usd2005_per_tco2))
    end

    println("\nWorld GDP (billion USD2005/yr), normal vs natural capital")
    println(rpad("year", 8), rpad("normal", 16), rpad("natural capital", 18), "ratio")
    for y in (2020, 2050, 2100, 2200, 2300)
        r = world[world.year .== y, :][1, :]
        println(rpad(y, 8), rpad(round(r.gdp_normal, digits = 1), 16),
                rpad(round(r.gdp_natural_capital, digits = 1), 18),
                round(r.ratio, digits = 5))
    end

    println("\nSCC 2020 (USD2005/tCO2), SSP245 + $(EMULATOR), b = $(B_LABEL)")
    println(rpad("formula", 15), rpad("forest dmg", 12), rpad("GDP", 18), "SCC")
    for row in eachrow(sccs)
        println(rpad(row.method, 15), rpad(label(row.forest_damages), 12),
                rpad(gdp_label(row.natural_capital), 18),
                round(row.scc_usd2005_per_tco2, digits = 4))
    end

    println("\nWrote docs/natural_capital_gdp.csv, docs/natural_capital_gdp_countries.csv,")
    println("docs/natural_capital_scc.csv and $(length(files)) figures in $FIGDIR")
    return (models = models, world_gdp = world, sccs = sccs, figures = files)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end