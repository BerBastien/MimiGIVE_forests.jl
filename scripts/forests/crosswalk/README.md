# Impact regions ↔ INEGI crosswalk (Mexico)

Links the Climate Impact Lab impact regions used by the forest emulator to
INEGI states and municipios (Marco Geoestadístico), so regional results can be
reported by state or municipio. Run the scripts from the repository root.

| Script | What it does | Output |
|---|---|---|
| `ir_inegi_crosswalk.R` | Intersects impact regions with INEGI municipios (EPSG:6372) | `data/Forest/crosswalk/ir_inegi_*.csv`, check map |
| `forest_by_state.R` | Forest-area change by state: boxplot of regions, state totals with p05–p95 | `figures/crosswalk/` |
| `forest_by_municipio.R` | Forest-area change by municipio for one state (default Tabasco) | `figures/crosswalk/` |

Inputs not in the repo (large): CIL `agglomerated-world-new.shp`
(github.com/ClimateImpactLab/data-impact-region-shapefile) and INEGI `00ent.shp`,
`00mun.shp`. Set `dir_inputs` / `ruta_mun` at the top of each script.
`forest_by_*` read the Monte Carlo exports in `forest_mcs_output/forest_ir_<SSP>.csv`.

Aggregation: `value_municipio = Σ_regions value_region × pct_of_region / 100`
(area-weighted; assumes the value is spread evenly inside each region).

Known issues
- CIL region codes use GADM state numbering, not INEGI codes; match by geometry, not by code.
- About 20 of 500 regions fall mostly in a state other than the one in their code
  (Baja California / BCS naming and GADM vs INEGI borders).
- State p05/p95 add up regional percentiles: an envelope, wider than the true percentile.
