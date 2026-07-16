#!/usr/bin/env python3
"""
build_impact_region_country_crosswalk.py
========================================
Derive the impact_region_id -> give_country_id crosswalk by spatially joining the
impact-region polygons to a world country map, then mapping the world's ISO-A3
codes onto the GIVE country list (data/Dimension_countries.csv).

Why: the user does not have an authoritative crosswalk. Each impact region sits
inside exactly one country; a point-in-polygon join of the region's representative
point against a world map gives us that country. Regions that fall outside all
polygons (small islands, coastline slivers) are assigned to the NEAREST country.

Output (written to --out, default forest_module/data/processed/):
    impact_region_country_crosswalk.csv     impact_region_id, give_country_id, country_name, match_type
    crosswalk_diagnostics.csv               per-region join diagnostics
    crosswalk_unmatched_iso.csv             world ISO codes not present in GIVE (if any)

Example:
    python build_impact_region_country_crosswalk.py \
        --regions "/path/impact_region_spatial_weights.gpkg" \
        --give-countries "../../data/Dimension_countries.csv" \
        --world "/path/ne_10m_admin_0_countries.shp"     # optional; else auto

Requires: geopandas (pulls in fiona/pyogrio, shapely, pandas). Install with
    pip install geopandas
"""
import argparse
import os
import sys

try:
    import geopandas as gpd
    import pandas as pd
except ImportError as e:  # pragma: no cover
    sys.exit("This script needs geopandas + pandas. Install with `pip install geopandas`.\n"
             f"Import error: {e}")

# Candidate column names, in priority order.
IR_ID_CANDIDATES = ["impact_region_id", "hierid", "region", "region_id", "GID",
                    "gadmid", "ISO_hierid", "id", "OBJECTID"]
ISO3_CANDIDATES = ["ISO_A3_EH", "ISO_A3", "ADM0_A3", "ADM0_A3_US", "iso_a3",
                   "SOV_A3", "GU_A3", "BRK_A3"]
NAME_CANDIDATES = ["ADMIN", "NAME", "NAME_EN", "SOVEREIGNT", "name"]


def pick_column(cols, candidates, what):
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand in cols:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def load_world(world_path):
    if world_path:
        print(f"[world] reading {world_path}")
        return gpd.read_file(world_path)
    # Try geopandas' built-in dataset (removed in geopandas >= 1.0).
    try:
        path = gpd.datasets.get_path("naturalearth_lowres")
        print(f"[world] using geopandas built-in naturalearth_lowres: {path}")
        return gpd.read_file(path)
    except Exception as e:  # pragma: no cover
        sys.exit(
            "No --world path given and the geopandas built-in world dataset is not "
            "available in this geopandas version.\n"
            "Download Natural Earth 'Admin 0 – Countries' (1:10m or 1:110m) from\n"
            "  https://www.naturalearthdata.com/downloads/\n"
            "and pass it with --world /path/to/ne_10m_admin_0_countries.shp\n"
            f"(underlying error: {e})"
        )


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_out = os.path.normpath(os.path.join(here, "..", "data", "processed"))
    default_give = os.path.normpath(os.path.join(here, "..", "..", "data", "Dimension_countries.csv"))

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True, help="Impact-region polygon file (.gpkg/.shp/.parquet)")
    ap.add_argument("--layer", default=None, help="Layer name if the file has multiple layers")
    ap.add_argument("--world", default=None, help="World country polygons (Natural Earth). Optional.")
    ap.add_argument("--give-countries", default=default_give, help="GIVE Dimension_countries.csv")
    ap.add_argument("--ir-id-col", default=None, help="Impact-region id column (auto-detected if omitted)")
    ap.add_argument("--iso-col", default=None, help="World ISO-A3 column (auto-detected if omitted)")
    ap.add_argument("--out", default=default_out, help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ---- impact regions ----------------------------------------------------
    print(f"[regions] reading {args.regions}")
    regions = gpd.read_file(args.regions, layer=args.layer) if args.layer \
        else gpd.read_file(args.regions)
    print(f"[regions] {len(regions)} features; columns: {list(regions.columns)}")

    ir_col = args.ir_id_col or pick_column(regions.columns, IR_ID_CANDIDATES, "impact region id")
    if ir_col is None:
        sys.exit(f"Could not auto-detect the impact-region id column. "
                 f"Pass --ir-id-col. Available: {list(regions.columns)}")
    print(f"[regions] using impact_region_id column: {ir_col}")

    regions = regions[[ir_col, regions.geometry.name]].rename(columns={ir_col: "impact_region_id"})
    regions["impact_region_id"] = regions["impact_region_id"].astype(str)
    if regions.crs is None:
        print("[regions] WARNING: no CRS set; assuming EPSG:4326")
        regions = regions.set_crs(4326)
    regions = regions.to_crs(4326)

    # ---- world -------------------------------------------------------------
    world = load_world(args.world).to_crs(4326)
    iso_col = args.iso_col or pick_column(world.columns, ISO3_CANDIDATES, "world ISO3")
    name_col = pick_column(world.columns, NAME_CANDIDATES, "world name")
    if iso_col is None:
        sys.exit(f"Could not detect an ISO-A3 column in the world file. "
                 f"Pass --iso-col. Available: {list(world.columns)}")
    print(f"[world] {len(world)} countries; ISO col = {iso_col}; name col = {name_col}")

    keep = [iso_col, world.geometry.name] + ([name_col] if name_col else [])
    world = world[keep].rename(columns={iso_col: "world_iso3"})
    if name_col:
        world = world.rename(columns={name_col: "world_name"})
    else:
        world["world_name"] = world["world_iso3"]
    world["world_iso3"] = world["world_iso3"].astype(str).str.upper()

    # ---- point-in-polygon on representative points -------------------------
    pts = regions.copy()
    pts["geometry"] = regions.representative_point()

    joined = gpd.sjoin(pts, world, how="left", predicate="within")
    joined = joined[~joined.index.duplicated(keep="first")]  # guard against overlaps
    joined["match_type"] = joined["world_iso3"].notna().map({True: "within", False: "unmatched"})

    # ---- nearest fallback for points outside all polygons ------------------
    missing = joined["world_iso3"].isna()
    n_missing = int(missing.sum())
    if n_missing:
        print(f"[join] {n_missing} regions had no containing country; assigning nearest.")
        near = gpd.sjoin_nearest(pts.loc[missing, ["impact_region_id", "geometry"]],
                                 world, how="left")
        near = near[~near.index.duplicated(keep="first")]
        joined.loc[missing, "world_iso3"] = near["world_iso3"].values
        joined.loc[missing, "world_name"] = near["world_name"].values
        joined.loc[missing, "match_type"] = "nearest"

    # ---- map world ISO3 -> GIVE ISO3 ---------------------------------------
    give = pd.read_csv(args.give_countries)
    give_iso = set(give.iloc[:, 0].astype(str).str.upper())
    print(f"[give] {len(give_iso)} GIVE countries")

    result = pd.DataFrame({
        "impact_region_id": joined["impact_region_id"].values,
        "give_country_id": joined["world_iso3"].astype(str).str.upper().values,
        "country_name": joined.get("world_name", pd.Series(joined["world_iso3"])).values,
        "match_type": joined["match_type"].values,
    })

    in_give = result["give_country_id"].isin(give_iso)
    unmatched_iso = sorted(set(result.loc[~in_give, "give_country_id"]))
    if unmatched_iso:
        print(f"[give] WARNING: {int((~in_give).sum())} regions map to ISO codes not in GIVE: "
              f"{unmatched_iso[:15]}{' ...' if len(unmatched_iso) > 15 else ''}")
        print("       These are usually territories/dependencies. Add a manual remap or drop them.")

    # ---- write -------------------------------------------------------------
    xwalk_path = os.path.join(args.out, "impact_region_country_crosswalk.csv")
    result[["impact_region_id", "give_country_id", "country_name"]].to_csv(xwalk_path, index=False)
    print(f"[out] crosswalk -> {xwalk_path}  ({len(result)} rows)")

    diag = result.copy()
    diag["in_give"] = in_give.values
    diag_path = os.path.join(args.out, "crosswalk_diagnostics.csv")
    diag.to_csv(diag_path, index=False)
    print(f"[out] diagnostics -> {diag_path}")

    if unmatched_iso:
        pd.DataFrame({"world_iso3_not_in_give": unmatched_iso}).to_csv(
            os.path.join(args.out, "crosswalk_unmatched_iso.csv"), index=False)

    # summary
    counts = result["match_type"].value_counts().to_dict()
    print(f"[summary] match types: {counts}")
    print(f"[summary] regions mapped to a valid GIVE country: {int(in_give.sum())}/{len(result)}")


if __name__ == "__main__":
    main()
