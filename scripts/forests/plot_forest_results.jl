# ------------------------------------------------------------------------------
# Article figures from the CSVs written by runforest_mcs.jl
#
#   Fig. 1  Global temperature by scenario (mean and p05-p95)
#   Fig. 2  Total forest area (Mha): one panel per scenario,
#           one color per vegetation model (mean and p05-p95)
#   Fig. 3  Total value of each ecosystem service: rows = services,
#           columns = scenarios, one color per vegetation model
#
# Usage:
#   include("scripts/plot_forest_results.jl")
#   make_figures(outdir = "forest_one", figdir = "figures_test")
#   make_figures(mode = :bars)      # error bars instead of shaded band
#
# NOTE on Fig. 2 and 3: totals sum over IRs. The mean of the sum is the sum
# of the means (exact). Quantiles of the sum are approximated by summing the
# quantiles of each IR, which is exact only if IRs move together (comonotonic)
# across trials. Temperature is the only random input, so this is usually a
# good approximation, but it should be stated in the paper (or the quantiles
# of the total should be computed per trial).
# ------------------------------------------------------------------------------

using CSV, DataFrames, CairoMakie, Statistics

# ------------------------------------------------------------------------------
# Constants

const SCEN_COLORS = Dict(
    "SSP119" => "#1b9e77", "SSP126" => "#66a61e", "SSP245" => "#e6ab02",
    "SSP370" => "#d95f02", "SSP585" => "#b2182b")
const VEG_COLORS = Dict("jules" => "#0072B2", "lpjml" => "#D55E00", "mc2" => "#009E73")
const VEG_ORDER  = ["jules", "lpjml", "mc2"]
const SERVICES   = ["rec", "hab", "nwfp", "wat"]
const SERVICE_LABELS = Dict(
    "rec"  => "Recreation",
    "hab"  => "Habitat",
    "nwfp" => "Non-wood forest products",
    "wat"  => "Water regulation")

# units of the services (not stored in the CSVs)
const VALUE_UNIT = "US\$/ha/year"

# ------------------------------------------------------------------------------
# Helpers

function load_results(outdir)

    # read and stack every forest_ir_*.csv in the folder
    files = filter(f -> occursin(r"^forest_ir_.*\.csv$", f), readdir(outdir))
    isempty(files) && error("No forest_ir_*.csv files in $outdir")

    return reduce(vcat, [CSV.read(joinpath(outdir, f), DataFrame) for f in files])
end

scen_order(df)  = sort(filter(s -> s in unique(df.scenario), collect(keys(SCEN_COLORS))))
veg_present(df) = filter(v -> v in unique(df.veg_model), VEG_ORDER)

# draw mean + p05-p95 range (shaded band or error bars)
function draw_range!(ax, x, lo, m, hi, color; mode = :band, dodge = 0.0, label = nothing)

    if mode == :band
        band!(ax, x, lo, hi, color = (color, 0.2))
        lines!(ax, x, m, color = color, linewidth = 2, label = label)
    else
        xs = x .+ dodge                                  # shift x so models do not overlap
        rangebars!(ax, xs, lo, hi, color = color, whiskerwidth = 5, linewidth = 1.3)
        lines!(ax, xs, m, color = color, linewidth = 2, label = label)
        scatter!(ax, xs, m, color = color, markersize = 6)
    end
end

# sum over IRs by scenario, model and year (mean, p05, p95)
function totals(df, base)

    combine(groupby(df, [:scenario, :veg_model, :year]),
            Symbol(base * "_p05")  => sum => :lo,
            Symbol(base * "_mean") => sum => :m,
            Symbol(base * "_p95")  => sum => :hi)
end

function model_legend_elements(vegs, mode)

    return [[LineElement(color = VEG_COLORS[v], linewidth = 2),
             MarkerElement(color = VEG_COLORS[v], marker = :circle, markersize = 6)] for v in vegs]
end

# x offset of the k-th of nv series (error-bar mode only)
model_dodge(k, nv) = (k - (nv + 1) / 2) * 0.8

# ------------------------------------------------------------------------------
# Fig. 1: temperature

function fig_temperature(df, figdir; mode = :band)

    # T is the same for every IR and model: one value per (scenario, year)
    t = combine(groupby(df, [:scenario, :year]),
                :T_mean => first => :m, :T_p05 => first => :lo, :T_p95 => first => :hi)
    scens = scen_order(df)

    fig = Figure(size = (620, 420))
    ax  = Axis(fig[1, 1], xlabel = "Year", ylabel = "Global temperature (°C above pre-industrial)")

    # one range per scenario
    for (i, s) in enumerate(scens)
        d = sort(t[t.scenario .== s, :], :year)
        draw_range!(ax, d.year, d.lo, d.m, d.hi, SCEN_COLORS[s];
                    mode = mode, dodge = model_dodge(i, length(scens)), label = s)
    end
    axislegend(ax, position = :lt, framevisible = false)

    save(joinpath(figdir, "fig1_temperature.png"), fig, px_per_unit = 3)
    save(joinpath(figdir, "fig1_temperature.pdf"), fig)

    return fig
end

# ------------------------------------------------------------------------------
# Fig. 2: forest area, one panel per scenario

function fig_forest_area(df, figdir; mode = :band)

    t     = totals(df, "F_Mha")
    scens = scen_order(df)
    vegs  = veg_present(df)
    ns    = length(scens)

    fig = Figure(size = (max(520, 300 * ns + 140), 380))
    axs = Axis[]

    # one panel per scenario, one color per vegetation model
    for (j, s) in enumerate(scens)
        ax = Axis(fig[1, j], title = s, xlabel = "Year",
                  ylabel = j == 1 ? "Forest area (Mha)" : "")
        push!(axs, ax)
        for (k, v) in enumerate(vegs)
            d = sort(t[(t.scenario .== s) .& (t.veg_model .== v), :], :year)
            isempty(d) && continue
            draw_range!(ax, d.year, d.lo, d.m, d.hi, VEG_COLORS[v];
                        mode = mode, dodge = model_dodge(k, length(vegs)))
        end
        j > 1 && hideydecorations!(ax, grid = false)
    end
    linkyaxes!(axs...)                                   # same y scale in all panels
    Legend(fig[1, ns + 1], model_legend_elements(vegs, mode), vegs, "Model", framevisible = false)

    save(joinpath(figdir, "fig2_forest_area.png"), fig, px_per_unit = 3)
    save(joinpath(figdir, "fig2_forest_area.pdf"), fig)

    return fig
end

# ------------------------------------------------------------------------------
# Fig. 3: ecosystem services (rows) x scenarios (columns)

function fig_services(df, figdir; mode = :band)

    scens = scen_order(df)
    vegs  = veg_present(df)
    ns    = length(scens)

    fig = Figure(size = (max(560, 300 * ns + 160), 220 * length(SERVICES) + 80))

    # one row per service
    for (i, svc) in enumerate(SERVICES)
        t   = totals(df, svc)
        row = Axis[]

        # one column per scenario
        for (j, s) in enumerate(scens)
            ax = Axis(fig[i, j],
                      title  = i == 1 ? s : "",
                      xlabel = i == length(SERVICES) ? "Year" : "",
                      ylabel = j == 1 ? "$(SERVICE_LABELS[svc])\n$VALUE_UNIT" : "")
            push!(row, ax)
            for (k, v) in enumerate(vegs)
                d = sort(t[(t.scenario .== s) .& (t.veg_model .== v), :], :year)
                isempty(d) && continue
                draw_range!(ax, d.year, d.lo, d.m, d.hi, VEG_COLORS[v];
                            mode = mode, dodge = model_dodge(k, length(vegs)))
            end
            i < length(SERVICES) && hidexdecorations!(ax, grid = false)
            j > 1 && hideydecorations!(ax, grid = false)
        end
        linkyaxes!(row...)                               # same y scale within each service
    end
    Legend(fig[:, ns + 1], model_legend_elements(vegs, mode), vegs, "Model", framevisible = false)

    save(joinpath(figdir, "fig3_services.png"), fig, px_per_unit = 3)
    save(joinpath(figdir, "fig3_services.pdf"), fig)

    return fig
end

# ------------------------------------------------------------------------------
# Run everything

function make_figures(; outdir = "forest_mcs_output", figdir = "figures", mode = :band)

    mkpath(figdir)
    set_theme!(Theme(fontsize = 12, Axis = (xgridvisible = false, ygridvisible = false)))

    df = load_results(outdir)
    @info "Loaded: $(nrow(df)) rows, scenarios = $(scen_order(df)), models = $(veg_present(df))"

    fig_temperature(df, figdir; mode = mode)
    fig_forest_area(df, figdir; mode = mode)
    fig_services(df, figdir; mode = mode)

    @info "Figures saved in $figdir/ (PNG at 300 dpi and vector PDF)"
end