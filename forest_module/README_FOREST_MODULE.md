# Forest Ecosystem-Services Damage Module for MimiGIVE

A spatially explicit, modular extension that adds a forest ecosystem-services
damage sector to MimiGIVE. It takes global mean surface temperature (GMST) from
GIVE's climate module, projects impact-region forest area, downscales GIVE
country population and GDP to impact regions, computes a **placeholder**
ecosystem-service value and damage, aggregates back to GIVE countries, and
produces a forest Social Cost of Carbon (SCC) using GIVE's own marginal-model and
discounting machinery.

> **This module is non-destructive.** It lives entirely under `forest_module/`
> and does not modify any file in MimiGIVE's own `src/`. You attach it to a normal
> GIVE model at build time.

> **Scientific limitation, up front.** The valuation is a placeholder,
> `ESV = scale · forest_area · GDP · population`. It exists to prove the spatial +
> IAM pipeline end-to-end. Its units depend on the inputs and it has **no welfare
> interpretation**. Do **not** report any number from this module as a
> policy-relevant ecosystem-service value or SCC until the placeholder is replaced
> with a defensible valuation function (one function, `ecosystem_service_value`,
> is the only thing that must change).

---

## 1. Scientific objective

Establish a valid, modular pipeline:

```
GIVE GMST
  → impact-region forest response      (temperature → forest area)
  → impact-region socioeconomic exposure (downscaled population & GDP)
  → impact-region ecosystem-service value (placeholder)
  → country damages                     (aggregate impact regions → GIVE countries)
  → marginal damages                    (pulse model − base model)
  → forest SCC                          (discounted, GIVE convention)
```

The whole chain is designed so the placeholder valuation can be swapped out without
touching the temperature response, the downscaling, the aggregation, the country
outputs, or the SCC code.

## 2. Two spatial units — do not confuse them

| Term | Meaning | Identifiers used in code |
|------|---------|--------------------------|
| **GIVE country** | Country-level unit GIVE already uses for population/GDP (183 ISO3 codes) | `give_country`, `country`, `country_index` |
| **Impact region** | Finer unit used **only** inside this module | `impact_region`, `impact_region_id`, `ir` |

Every impact region belongs to **exactly one** GIVE country. The authoritative
`impact_region_id → give_country_id` crosswalk is built by
`scripts/build_impact_region_country_crosswalk.py` (a spatial join of the
impact-region polygons to a world map, mapped onto GIVE ISO3 codes).

*R analogy:* impact regions are like fine `sf` polygons; GIVE countries are the
coarse `group_by()` key you sum up to.

## 3. GMST 2022 baseline

The forest coefficients were fitted on **ΔT relative to 2022**, so the module
re-centres GMST on 2022:

```
ΔT_t = GMST_t − GMST_2022
```

Critically, `GMST_2022` is captured **inside each model run** from that run's own
FaIR trajectory (using the same "store-the-baseline-year-value" trick as GIVE's
`GlobalTempNorm`). The base model and the pulse model therefore each use their
**own** 2022 value — never an external observational number. GIVE's annual time
index (1750–2300) always contains 2022; if you change `GMST_BASELINE_YEAR` to a
year outside the model, the build errors unless `GMST_BASELINE_FALLBACK` is set to
`:nearest`.

## 4. Forest-response equation

For each impact region *r* and year *t*:

```
change_{r,t} = α_r + β1_r · ΔT_t + β2_r · ΔT_t²
```

`α_r` (intercept) is **only** used if your coefficient file supplies it; it defaults
to 0. Conversion to area depends on `FOREST_CHANGE_UNITS`:

| `FOREST_CHANGE_UNITS` | Area formula |
|-----------------------|--------------|
| `:percent` (**default, confirmed**) | `A_{r,t} = A_{r,2022} · (1 + change/100)` |
| `:proportion` | `A_{r,t} = A_{r,2022} · (1 + change)` |
| `:log_change` | `A_{r,t} = A_{r,2022} · exp(change)` |

## 5. Coefficient units

**Percentage points** (`5` means 5%). Set in `config.jl` as
`FOREST_CHANGE_UNITS = :percent`. Change it there if your fitted coefficients are
proportions or log-changes.

## 6. Forest-area units

Baseline forest area is in **Mha (millions of hectares)**
(`FOREST_AREA_UNITS = :Mha`). This is a documentation/label field — the placeholder
valuation is unit-agnostic — but keeping it explicit prevents silently mixing ha,
km², and m². Do not mix units in the baseline file.

## 7. Physical constraints

* `A_{r,t} ≥ 0` (lower clip, on by default).
* Optional `A_{r,t} ≤ total_impact_region_area_r` (upper clip, off by default; only
  applied if you supply a total-area column).

The raw, unconstrained projection is kept in `projected_forest_area_raw` for
diagnostics. Clipping and out-of-fit-range (extrapolation) events are counted in
`forest_clipping_diagnostic.csv`.

## 8. Socioeconomic downscaling

```
POP_{r,t} = POP_{c,t} · w^pop_{r,c}      (millions of persons)
GDP_{r,t} = GDP_{c,t} · w^gdp_{r,c}      (billion US$2005/yr)
GDPpc_{r,t} = GDP_{r,t} / POP_{r,t}      (never downscaled directly)
```

where `c` is the GIVE country containing region `r`. The shares `w^pop`, `w^gdp`
must sum to 1 within each country; the loader validates this, normalises tiny
rounding discrepancies within-country, and **errors** (refuses to normalise) if any
country is off by more than 1e-3.

## 9. Ecosystem-service reference and climate values

```
ESV^reference_{r,t} = f(A_{r,2022}, GDP_{r,t}, POP_{r,t})   # forest held at 2022 baseline
ESV^climate_{r,t}   = f(A_{r,t},    GDP_{r,t}, POP_{r,t})   # forest responds to climate
Damage^forest_{r,t} = ESV^reference − ESV^climate
```

Positive damage = loss relative to keeping baseline forest; negative = gain. Value
and damage are kept as separate variables — total value is never labelled "damage".
`f` is the placeholder `scale · A · GDP · POP` (see §Limitation).

## 10. Country aggregation

For each GIVE country `c`, sum over its impact regions; then sum countries for the
global total. `es_damage_global` equals the country sum within numerical tolerance
(unit-tested).

## 11. Marginal model logic

Forest SCC reuses GIVE's marginal model:

```
mm = MimiGIVE.get_marginal_model(m; year, gas, pulse_size)   # base + pulse, delta = pulse units
MD^forest_{c,t} = ( Damage^forest,pulse_{c,t} − Damage^forest,base_{c,t} ) / delta · molecular_conv
```

Because the forest components run inside both base and pulse models, and both derive
their own 2022 baseline, marginal forest damage is **exactly zero before the pulse
perturbs the climate** and small-but-nonzero afterwards (unit-tested on a standalone
chain).

## 12. Integration with the SCC

Discounting follows GIVE's non-equity-weighted Ramsey convention, using the base
model's global net consumption per capita `cpc`:

```
df_i = (cpc[year]/cpc[i])^eta · 1/(1+prtp)^(t−year)
SCC^forest_c = Σ_t df_t · MD^forest_{c,t}
SCC^forest   = Σ_c SCC^forest_c
```

`compute_forest_scc` returns global and per-country forest SCC plus the full
marginal-damage matrix.

**By default the forest sector does NOT alter GIVE's official SCC**
(`INCLUDE_FOREST_IN_TOTAL_SCC = false`). When you set it `true`,
`compute_total_scc_with_forest` adds the forest SCC (converted from placeholder
units via `FOREST_DAMAGE_TO_USD2005`) to GIVE's official SCC — non-destructively,
without editing `DamageAggregator` — and prints a prominent warning that the number
is not calibrated.

## 13. Placeholder valuation limitation

Repeated because it matters: `ESV = scale · forest_area · GDP · population` is **not**
a scientific valuation. Replace the body of `ecosystem_service_value` (in
`src/forest_ecosystem_services.jl`) — nothing else in the pipeline needs to change.

## 14. Required input schemas

Preprocessing (Python) writes these compact CSVs into `data/processed/`, which Julia
reads once at build time:

| File | Required columns | Optional columns |
|------|------------------|------------------|
| `impact_region_country_crosswalk.csv` | `impact_region_id`, `give_country_id` | `country_name` |
| `forest_coefficients.csv` | `impact_region_id`, `beta_delta_gmst`, `beta_delta_gmst_squared` | `intercept`, `delta_gmst_fit_min`, `delta_gmst_fit_max` |
| `baseline_forest_area.csv` | `impact_region_id`, `baseline_forest_area` | `total_impact_region_area` |
| `population_shares.csv` | `impact_region_id`, `give_country_id`, `population_share` | |
| `gdp_shares.csv` | `impact_region_id`, `give_country_id`, `gdp_share` | |

`give_country_id` must be GIVE **ISO3** codes (matching `data/Dimension_countries.csv`).
Coefficients with extra dimensions (SSP, year, biome, model) must be collapsed to
one row per impact region for the selected scenario **before** this step.

## 15. Commands to run

From the MimiGIVE.jl repo root, with the MimiGIVE project active:

```bash
# 0. (once) install geopandas for the preprocessing scripts
pip install geopandas matplotlib

# 1. Build the impact_region → GIVE-country crosswalk from polygons + a world map
python forest_module/scripts/build_impact_region_country_crosswalk.py \
    --regions ".../impact_region_spatial_weights.gpkg" \
    --give-countries data/Dimension_countries.csv
    # add --world ne_10m_admin_0_countries.shp if the geopandas built-in world is unavailable

# 2. Standardize the coefficient + weights files into the compact CSVs
python forest_module/scripts/preprocess_forest_inputs.py \
    --coefficients ".../impact_region_regression.gpkg" \
    --weights      ".../impact_region_spatial_weights.gpkg"
    # pass --coef-b1/--coef-b2/--pop-share-col/... to override auto-detected columns

# 3. Run the deterministic end-to-end example (uses synthetic fixtures if no data yet)
julia --project=. forest_module/scripts/run_forest_give_example.jl

# 4. Run the tests
julia --project=. forest_module/test/runtests.jl
# add the heavy full-GIVE integration test:
FOREST_RUN_FULL=true julia --project=. forest_module/test/runtests.jl

# 5. (optional) diagnostic maps for a chosen year
python forest_module/scripts/make_diagnostic_maps.py \
    --regions ".../impact_region_spatial_weights.gpkg" \
    --ir-output forest_module/outputs/forest_impact_region_output.csv \
    --year 2100
```

Programmatic use:

```julia
include("forest_module/ForestGIVE.jl"); using .ForestGIVE
built = build_forest_give_model(socioeconomics_source=:SSP, SSP_scenario="SSP245")
run(built.model)
scc = compute_forest_scc(built; year=2030, gas=:CO2, pulse_size=1.0)
write_country_output(built); write_summary_output(built, scc)
```

*Julia vs R note:* `include("…jl")` is like `source("…R")`; `using .ForestGIVE`
brings the exported functions into scope like `library()`. `build_forest_give_model`
returns a NamedTuple (`built.model`, `built.inputs`, …) — like a named `list()`.

## 16. Output files (written to `outputs/`)

* `forest_impact_region_output.csv` — per year × impact region (GMST, ΔT, forest change, areas, pop, GDP, GDPpc, ESV reference/climate, damage).
* `forest_country_output.csv` — per year × GIVE country (ESV reference/climate/damage).
* `forest_marginal_country_output.csv` — base/pulse/marginal forest damage, discount factor, discounted marginal damage.
* `forest_global_summary.csv`, `forest_scc_country_contribution.csv`, `forest_scc_headline.csv` — summaries.
* `forest_temperature_diagnostic.csv` — year, gmst, gmst_2022, ΔT, ΔT².
* `forest_clipping_diagnostic.csv` — clipping/extrapolation counts per region.
* `forest_join_diagnostics.csv` — the aligned per-region input table.
* Maps + a GeoPackage for a selected year (from `make_diagnostic_maps.py`).

## 17. Tests

`test/runtests.jl` runs, with only `Mimi` needed:

* **Temperature** — ΔT(2022)=0, ΔT²(2022)=0, and (no intercept) forest change=0 / area=baseline in 2022.
* **Coefficient** — β1=2, β2=1, ΔT=0.5 ⇒ change=1.25% ⇒ area = baseline·1.0125.
* **Clipping / zero-response** — area never negative; all-zero coefficients ⇒ area=baseline.
* **Downscaling** — Σ regional pop = country pop, Σ regional GDP = country GDP.
* **Aggregation** — Σ impact-region damage = country damage; Σ country = global.
* **Marginal model** — base ≡ pulse before divergence; nonzero after; country marginals sum to global marginal.

The opt-in `FOREST_RUN_FULL=true` test builds a real GIVE model and checks the 2022
invariants plus a finite, additive forest SCC.

## 18. Component architecture

| Component | Inputs | Outputs |
|-----------|--------|---------|
| `SpatialSocioeconomics` | `population_country`, `gdp_country` (from GIVE), shares, IR→country index | `population_ir`, `gdp_ir`, `gdppc_ir` |
| `ForestAreaResponse` | `gmst` (from `:temperature=>:T`), β1, β2, intercept, baseline area, (total area) | `delta_gmst`, `forest_change`, `projected_forest_area_raw`, `projected_forest_area`, clip/extrapolation flags |
| `ForestEcosystemServices` | baseline area, `projected_forest_area`, `population_ir`, `gdp_ir`, IR→country index | ESV reference/climate/damage at IR, country, global |

Components are added **after** `:country_netconsumption` so `:temperature` and
`:Socioeconomic` have already run. Heavy geospatial reads happen once in Python; the
runtime components operate only on numeric arrays.

## 19. Known limitations / unresolved issues

* **Placeholder valuation** — not calibrated; units undefined; SCC not policy-relevant.
* **Coefficient collapsing** — if the coefficient file has SSP/year/biome dimensions, the Python step must reduce it to one row per impact region for the chosen scenario; the Julia loader takes the last row per duplicated id.
* **Crosswalk quality** — derived by point-in-polygon + nearest fallback; territories/dependencies whose ISO3 is not in GIVE are reported in `crosswalk_unmatched_iso.csv` and must be remapped or dropped.
* **Column auto-detection** — the preprocessing scripts guess column names; always check the printed schema and override with the CLI flags.
* **`INCLUDE_FOREST_IN_TOTAL_SCC=true`** combines placeholder units with dollars via a user-set factor with no defensible value — for plumbing tests only.
* **Units of ΔT fit-range** — the extrapolation flag compares ΔT (not absolute GMST) to the supplied `delta_gmst_fit_min/max`.

## 20. File map

```
forest_module/
  ForestGIVE.jl                         # module entry (include + using .ForestGIVE)
  config.jl                             # paths + all tunable knobs (ForestConfig)
  src/
    forest_spatial_socioeconomics.jl    # Component A
    forest_area_response.jl             # Component B
    forest_ecosystem_services.jl        # Component C + ecosystem_service_value()
    forest_inputs.jl                    # load + validate compact CSVs, join diagnostics
    build_forest_give_model.jl          # wire components into a real GIVE model
    forest_scc.jl                       # forest marginal damages + SCC
    forest_outputs.jl                   # tidy CSV writers + diagnostics
  scripts/
    build_impact_region_country_crosswalk.py
    preprocess_forest_inputs.py
    make_diagnostic_maps.py
    make_synthetic_fixtures.jl          # valid fake inputs so it runs before real data
    run_forest_give_example.jl          # end-to-end deterministic example
  test/
    runtests.jl + 4 test files
  data/processed/                       # compact CSVs land here (generated)
  outputs/                              # model outputs land here
  README_FOREST_MODULE.md
```
