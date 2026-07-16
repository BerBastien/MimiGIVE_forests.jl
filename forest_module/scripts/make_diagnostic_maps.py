#!/usr/bin/env python3
"""
make_diagnostic_maps.py
=======================
Join the Julia model's impact-region CSV output back onto the impact-region
polygons and render the diagnostic maps (section 20) for a selected year, plus a
single GeoPackage layer with all mapped variables.

The Mimi model runs in Julia and stays geometry-free; this script does the
geometry join/rendering afterwards.

Maps produced (PNG) for --year:
    1 baseline_forest_area
    2 forest_change (%)
    3 projected_forest_area
    4 population_ir
    5 gdp_ir
    6 es_damage
    (marginal ES damage map, #7, uses forest_marginal_country_output.csv joined at
     country level if you pass --marginal.)

Example:
    python make_diagnostic_maps.py \
        --regions ".../impact_region_spatial_weights.gpkg" \
        --ir-output "../outputs/forest_impact_region_output.csv" \
        --year 2100 --out ../outputs/maps

Requires: geopandas, matplotlib.
"""
import argparse
import os
import sys

try:
    import geopandas as gpd
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:  # pragma: no cover
    sys.exit(f"Needs geopandas + matplotlib. `pip install geopandas matplotlib`. ({e})")

IR_ID_CANDIDATES = ["impact_region_id", "hierid", "region", "region_id", "GID", "id"]
MAP_VARS = [
    ("baseline_forest_area", "Baseline forest area"),
    ("forest_change", "Projected forest change (coef units)"),
    ("projected_forest_area", "Projected forest area"),
    ("population_ir", "Impact-region population"),
    ("gdp_ir", "Impact-region GDP"),
    ("es_damage", "Ecosystem-service damage (placeholder)"),
]


def pick(cols, cands):
    lower = {c.lower(): c for c in cols}
    for c in cands:
        if c in cols: return c
        if c.lower() in lower: return lower[c.lower()]
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--layer", default=None)
    ap.add_argument("--ir-output", required=True, help="forest_impact_region_output.csv")
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--ir-id-col", default=None)
    ap.add_argument("--out", default="maps")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    regions = gpd.read_file(args.regions, layer=args.layer) if args.layer \
        else gpd.read_file(args.regions)
    ir_col = args.ir_id_col or pick(regions.columns, IR_ID_CANDIDATES)
    if ir_col is None:
        sys.exit(f"Could not detect IR id column. Use --ir-id-col. Columns: {list(regions.columns)}")
    regions = regions.rename(columns={ir_col: "impact_region_id"})
    regions["impact_region_id"] = regions["impact_region_id"].astype(str)

    df = pd.read_csv(args.ir_output)
    df = df[df["year"] == args.year].copy()
    if df.empty:
        sys.exit(f"No rows for year {args.year} in {args.ir_output}. "
                 f"Available: {sorted(pd.read_csv(args.ir_output)['year'].unique())}")
    df["impact_region_id"] = df["impact_region_id"].astype(str)

    gdf = regions.merge(df, on="impact_region_id", how="left")

    # GeoPackage with all mapped variables for the selected year
    gpkg_path = os.path.join(args.out, f"forest_impact_region_{args.year}.gpkg")
    gdf.to_file(gpkg_path, driver="GPKG")
    print(f"[out] {gpkg_path}")

    for col, title in MAP_VARS:
        if col not in gdf.columns:
            print(f"  (skip {col}: not in output)")
            continue
        fig, ax = plt.subplots(1, 1, figsize=(11, 6))
        gdf.plot(column=col, ax=ax, legend=True, cmap="viridis",
                 missing_kwds={"color": "lightgrey"})
        ax.set_title(f"{title} — {args.year}")
        ax.set_axis_off()
        png = os.path.join(args.out, f"map_{col}_{args.year}.png")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[out] {png}")

    print("Done.")


if __name__ == "__main__":
    main()
