# ------------------------------------------------------------------------------
# Monte Carlo of temperature by SSP -> forest area and ecosystem services by IR
#
# Uses the built-in MimiGIVE.ForestEmulator (formula is not copied).
# Output: one CSV per scenario with mean and quantiles across trials of T,
# forest hectares and the 4 services, by IR, year and vegetation model.
#
# Usage:
#   include("scripts/runforest_mcs.jl")
#   run_all(n = 10_000)
# ------------------------------------------------------------------------------

using Mimi, MimiGIVE, DataFrames, CSV, Statistics, Random

# ------------------------------------------------------------------------------
# Constants

const SCENARIOS   = ["SSP119", "SSP126", "SSP245", "SSP370", "SSP585"]
const VEG_MODELS  = [:jules, :lpjml, :mc2]
const SERVICES    = [:rec, :hab, :nwfp, :wat]   # same order as mv_cols in get_model
const QS          = [0.05, 0.25, 0.5, 0.75, 0.95]
const QNAMES      = ["p05", "p25", "p50", "p75", "p95"]
const YEARS_ALL   = 2022:2100                   # years the emulator runs
const TI          = 273:351                     # those years as rows of the 1750-2300 axis
const OUT_YEARS   = [2022, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100]
const N_FAIR_SETS = 2237                        # number of FAIR parameter sets

stats(v) = (mean(v), quantile(v, QS)...)        # mean + quantiles
f64(x)   = Float64.(coalesce.(x, NaN))          # missing -> NaN

# ------------------------------------------------------------------------------
# 1. Temperature: GIVE Monte Carlo per scenario (unique FAIR sets only)

function temperature_trials(scen; ids, outdir)

    n = length(ids)
    m = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = scen)

    # run the MCS on the unique FAIR sets, save temperature only
    res = MimiGIVE.run_mcs(trials = n, m = m,
            fair_parameter_set = :deterministic, fair_parameter_set_ids = ids,
            save_list = [(:temperature, :T)],
            output_dir = joinpath(outdir, "mcs_raw_$scen"),
            results_in_memory = true)

    # reshape to a years x trials matrix
    df = getdataframe(res, :temperature, :T)
    T  = fill(NaN, length(YEARS_ALL), n)
    for r in eachrow(df)
        if 2022 <= r.time <= 2100 && !ismissing(r.T)
            T[r.time - 2021, r.trialnum] = r.T
        end
    end
    @assert !any(isnan, T) "Missing temperature values (did the variable name change?)"

    return T
end

# cache temperature per scenario so the MCS is not repeated if a later step fails
function cached_temperature(scen, uids, outdir, tag)

    f = joinpath(outdir, "partial", "T_$(scen)_$(tag).csv")

    # read from cache if it exists
    if isfile(f)
        @info "$scen: temperature already computed, reading $f"
        return Matrix{Float64}(CSV.read(f, DataFrame))
    end

    # otherwise run the MCS and save it
    @info "$scen: temperature Monte Carlo ($(length(uids)) unique FAIR sets)"
    T = temperature_trials(scen; ids = uids, outdir = outdir)
    CSV.write(f, DataFrame(T, :auto))

    return T
end

# ------------------------------------------------------------------------------
# 2. Minimal model with the real ForestEmulator (for a subset of IRs)

# m_full: full GIVE model with forest_emulator, already run
# rows: indices of the IRs to include (default: all)
function build_forest_model(m_full, rows = nothing)

    g(v) = m_full[:ForestEmulator, v]
    ncou = size(g(:population), 2)
    nreg = length(g(:region_country))
    rows === nothing && (rows = 1:nreg)

    # dimensions and component
    m = Mimi.Model()
    set_dimension!(m, :time, YEARS_ALL)
    set_dimension!(m, :country, ncou)
    set_dimension!(m, :forest_regions, length(rows))
    set_dimension!(m, :services, 4)
    add_comp!(m, MimiGIVE.ForestEmulator, :ForestEmulator)

    # region parameters (subset to the selected IRs)
    update_param!(m, :ForestEmulator, :region_country, Int.(g(:region_country))[rows])
    for p in (:A_Mha, :F0, :b1, :b2, :C0, :c1, :c2, :pop_share, :gdp_share)
        update_param!(m, :ForestEmulator, p, f64(g(p))[rows])
    end
    update_param!(m, :ForestEmulator, :mv,    f64(g(:mv))[rows, :])

    # socioeconomic parameters
    update_param!(m, :ForestEmulator, :e_gdp, f64(g(:e_gdp)))
    update_param!(m, :ForestEmulator, :e_pop, f64(g(:e_pop)))
    update_param!(m, :ForestEmulator, :population, f64(g(:population))[TI, :])
    update_param!(m, :ForestEmulator, :pc_gdp,     f64(g(:pc_gdp))[TI, :])

    # placeholder, replaced by each temperature trajectory
    update_param!(m, :ForestEmulator, :temperature, zeros(length(YEARS_ALL)))

    return m
end

function load_inputs(scen, veg)

    # full GIVE model for this scenario and vegetation model
    m_full = MimiGIVE.get_model(socioeconomics_source = :SSP, SSP_scenario = scen,
                                forest_emulator = veg)
    run(m_full)

    # region -> country map and region table
    rc  = Int.(m_full[:ForestEmulator, :region_country])
    reg = CSV.read(joinpath(pkgdir(MimiGIVE), "data", "Forest", "regions_$(veg).csv"), DataFrame)
    @assert nrow(reg) == length(rc) "CSV and model have a different number of IRs"

    return (; m_full, rc, reg)
end

# the minimal model must reproduce the component inside full GIVE
function validate(inp)

    (; m_full) = inp
    m_forest = build_forest_model(m_full)

    # feed the same temperature and run
    T = f64(m_full[:ForestEmulator, :temperature])[TI]
    update_param!(m_forest, :ForestEmulator, :temperature, T)
    run(m_forest)

    # compare outputs
    for v in (:ir_forest_Mha, :ir_value)
        full = m_full[:ForestEmulator, v]
        a = f64(full[TI, ntuple(_ -> :, ndims(full) - 1)...])
        b = f64(m_forest[:ForestEmulator, v])
        @assert isapprox(a, b; rtol = 1e-8, atol = 1e-10) "Minimal model does not reproduce $v"
    end

    return true
end

# ------------------------------------------------------------------------------
# 3. Trajectories -> mean and quantiles by IR and year, in blocks of IRs
#    Tu: unique temperatures (years x unique); inv[k] = column of Tu for trial k

function stats_table(inp, scen, veg, Tu, inv, out_years; block = 500)

    (; m_full, rc, reg) = inp
    keep = findall(>(0), rc)                    # IRs without a country are dropped
    yi   = [y - 2021 for y in out_years]
    ny, nu, nk = length(yi), size(Tu, 2), length(keep)

    # output columns
    statnames = ["mean"; QNAMES]
    vars  = vcat(["F_Mha"], string.(SERVICES))
    order = ["$(v)_$(sn)" for v in vcat(["T"], vars) for sn in statnames]
    N     = nk * ny
    cols  = Dict(nm => zeros(N) for nm in order)

    # temperature stats (same for every IR)
    Tfull = Tu[:, inv]                          # years x n trials
    Tst   = [stats(Tfull[y, :]) for y in yi]

    # loop over blocks of IRs to bound memory
    for b0 in 1:block:nk

        idx = b0:min(b0 + block - 1, nk)
        mb  = build_forest_model(m_full, keep[idx])
        arr = Dict(v => Array{Float32}(undef, nu, ny, length(idx)) for v in vars)

        # run the emulator once per unique trajectory
        for k in 1:nu
            update_param!(mb, :ForestEmulator, :temperature, Tu[:, k])
            run(mb)
            iF = f64(mb[:ForestEmulator, :ir_forest_Mha])   # years x IRs in block
            iv = f64(mb[:ForestEmulator, :ir_value])        # years x IRs x service
            arr["F_Mha"][k, :, :] = iF[yi, :]
            for (s, sv) in enumerate(SERVICES)
                arr[string(sv)][k, :, :] = iv[yi, :, s]
            end
        end

        # expand to the n trials and compute stats
        for (jl, jr) in enumerate(idx), j in 1:ny
            row = (jr - 1) * ny + j
            for (i, sn) in enumerate(statnames)
                cols["T_$sn"][row] = Tst[j][i]
            end
            for v in vars
                st = stats(arr[v][inv, j, jl])
                for (i, sn) in enumerate(statnames)
                    cols["$(v)_$sn"][row] = st[i]
                end
            end
        end
        @info "$scen / $veg: IR $(idx.stop) of $nk"
    end

    # assemble the output table
    df = DataFrame(scenario  = fill(scen, N),
                   veg_model = fill(string(veg), N),
                   region_id = repeat(string.(reg.region_id[keep]), inner = ny),
                   iso       = repeat(string.(reg.iso[keep]), inner = ny),
                   year      = repeat(collect(out_years), outer = nk))
    for nm in order
        df[!, nm] = round.(cols[nm]; sigdigits = 6)
    end
    @info "$scen / $veg: $nk IRs with country, $(length(rc) - nk) without (omitted)"

    return df
end

# ------------------------------------------------------------------------------
# 4. Full run: one CSV per scenario (resumable)

function run_all(; n = nothing, seed = 1234, out_years = OUT_YEARS,
                   outdir = "forest_mcs_output",
                   scenarios = SCENARIOS, veg_models = VEG_MODELS,
                   block = 500)

    mkpath(joinpath(outdir, "partial"))

    # draw n FAIR sets with replacement (or use all of them if n is nothing)
    ids  = n === nothing ? collect(1:N_FAIR_SETS) :
                           rand(MersenneTwister(seed), 1:N_FAIR_SETS, n)
    ntr  = length(ids)
    uids = sort(unique(ids))                    # unique sets: the only ones actually run
    inv  = Int.(indexin(ids, uids))             # trial k -> column of the unique temperatures
    tag  = "n$(ntr)_s$(seed)"
    @info "Trials per scenario: $ntr ($(length(uids)) distinct FAIR sets)"

    for scen in scenarios

        parts = [joinpath(outdir, "partial", "$(scen)_$(veg)_$(tag).csv") for veg in veg_models]

        # temperature is only computed if some model is missing
        Tu = nothing
        for (veg, part) in zip(veg_models, parts)
            isfile(part) && (@info "$scen / $veg: already exists, skipping"; continue)
            Tu === nothing && (Tu = cached_temperature(scen, uids, outdir, tag))
            inp = load_inputs(scen, veg)
            validate(inp)
            CSV.write(part, stats_table(inp, scen, veg, Tu, inv, out_years; block = block))
        end

        # merge the vegetation models into one file per scenario
        file = joinpath(outdir, "forest_ir_$(scen).csv")
        CSV.write(file, vcat([CSV.read(p, DataFrame) for p in parts]...))
        @info "Saved $file"
    end
end