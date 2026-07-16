#!/usr/bin/env python3
"""
preprocess_forest_inputs.py
===========================
One-time preprocessing: read the raw spatial inputs, INSPECT their schema, and
write the small, validated CSVs that the Julia module loads
(forest_module/data/processed/). Do NOT let Julia read the geospatial files at
runtime — this is the "expensive read once" step.

It produces (into --out):
    forest_coefficients.csv        impact_region_id, beta_delta_gmst, beta_delta_gmst_squared[, intercept, delta_gmst_fit_min, delta_gmst_fit_max]
    baseline_forest_area.csv       impact_region_id, baseline_forest_area[, total_impact_region_area]
    population_shares.csv          impact_region_id, give_country_id, population_share
    gdp_shares.csv                 impact_region_id, give_country_id, gdp_share
    (crosswalk is produced by build_impact_region_country_crosswalk.py, or reused)
    preprocess_report.txt          human-readable schema + validation report

Because the real column names are unknown, the script auto-detects using common
name patterns and lets you override every column via CLI flags. It always PRINTS
the detected schema and refuses to guess silently on the critical share columns if
their per-country sums are far from 1.

Supported formats: .csv/.tsv (pandas), .gpkg/.shp/.geojson/.parquet/.geoparquet
(geopandas), .nc (xarray, optional).

Example (single weights gpkg holding shares + baseline; separate coeff gpkg):
    python preprocess_forest_inputs.py \
        --coefficients ".../impact_region_regression.gpkg" \
        --weights      ".../impact_region_spatial_weights.gpkg" \
        --coef-b1 beta_dT --coef-b2 beta_dT2 \
        --pop-share-col pop_share --gdp-share-col gdp_share \
        --baseline-col forest_area_2022 --country-col ISO3
"""
import argparse
import os
import sys

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    sys.exit("Needs pandas (and geopandas for spatial files). `pip install geopandas pandas`.")

IR_ID_CANDIDATES = ["impact_region_id", "hierid", "region", "region_id", "GID", "id"]
COUNTRY_CANDIDATES = ["give_country_id", "ISO3", "iso3", "ISO_A3", "country", "country_id", "GID_0"]
B1_CANDIDATES = ["beta_delta_gmst", "beta_dgmst", "beta1", "b1", "beta_dT", "slope"]
B2_CANDIDATES = ["beta_delta_gmst_squared", "beta_dgmst2", "beta2", "b2", "beta_dT2", "quad"]
INTERCEPT_CANDIDATES = ["intercept", "alpha", "b0", "beta0", "const"]
BASELINE_CANDIDATES = ["baseline_forest_area", "forest_area", "forest_area_2022",
                       "forest_2022", "area", "baseline_area"]
POP_SHARE_CANDIDATES = ["population_share", "pop_share", "pop_weight", "population_weight", "w_pop"]
GDP_SHARE_CANDIDATES = ["gdp_share", "gdp_weight", "w_gdp", "income_share"]
TOTAL_AREA_CANDIDATES = ["total_impact_region_area", "region_area", "total_area", "area_total"]
DTMIN_CANDIDATES = ["delta_gmst_fit_min", "dT_min", "gmst_min", "fit_min"]
DTMAX_CANDIDATES = ["delta_gmst_fit_max", "dT_max", "gmst_max", "fit_max"]


def read_any(path, layer=None):
    """Read a tabular OR spatial file into a (geo)pandas DataFrame."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".csv", ".tsv"):
        return pd.read_csv(path, sep="\t" if ext == ".tsv" else ",")
    if ext == ".nc":
        import xarray as xr
        return xr.open_dataset(path).to_dataframe().reset_index()
    import geopandas as gpd
    df = gpd.read_file(path, layer=layer) if layer else gpd.read_file(path)
    # drop geometry for the compact tabular outputs
    if hasattr(df, "geometry") and df.geometry.name in df.columns:
        df = pd.DataFrame(df.drop(columns=df.geometry.name))
    return df


def pick(cols, candidates, override=None):
    if override:
        if override not in cols:
            sys.exit(f"Requested column '{override}' not found. Available: {list(cols)}")
        return override
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand in cols:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def validate_shares(df, ir_col, country_col, share_col, label, report, hard_tol=1e-3, soft_tol=1e-6):
    sums = df.groupby(country_col)[share_col].sum()
    off = (sums - 1.0).abs()
    max_abs = float(off.max()) if len(off) else 0.0
    bad = sums[off > hard_tol]
    report.append(f"[{label}] countries={len(sums)}  max|sum-1|={max_abs:.3e}")
    if len(bad):
        report.append(f"[{label}] HARD FAIL: {len(bad)} countries off by > {hard_tol}. "
                      f"Worst: {bad.abs().sub(1).abs().idxmax()} -> {float(bad.iloc[0]):.4f}")
        print(f"  !! {label}: {len(bad)} countries have shares not summing to 1 (max {max_abs:.3e}).")
        print("     Refusing to normalise invalid weights. Fix upstream or inspect these countries:")
        print("     ", list(bad.index[:15]))
        raise SystemExit(1)
    if max_abs > soft_tol:
        report.append(f"[{label}] normalising within-country (tiny rounding, max {max_abs:.3e}).")
        df[share_col] = df[share_col] / df[country_col].map(sums)
    return df, max_abs


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_out = os.path.normpath(os.path.join(here, "..", "data", "processed"))

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coefficients", required=True)
    ap.add_argument("--weights", required=True,
                    help="File holding baseline forest area and pop/gdp shares (may be same as --coefficients)")
    ap.add_argument("--coef-layer", default=None)
    ap.add_argument("--weights-layer", default=None)
    ap.add_argument("--out", default=default_out)
    # column overrides
    for flag in ["ir-id-col", "country-col", "coef-b1", "coef-b2", "intercept-col",
                 "baseline-col", "total-area-col", "pop-share-col", "gdp-share-col",
                 "dtmin-col", "dtmax-col"]:
        ap.add_argument(f"--{flag}", default=None)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    report = []

    # ---- coefficients ------------------------------------------------------
    coef = read_any(args.coefficients, args.coef_layer)
    print(f"[coefficients] {len(coef)} rows; columns: {list(coef.columns)}")
    report.append(f"[coefficients] file={args.coefficients} rows={len(coef)} columns={list(coef.columns)}")
    ir_c = pick(coef.columns, IR_ID_CANDIDATES, args.ir_id_col)
    b1 = pick(coef.columns, B1_CANDIDATES, args.coef_b1)
    b2 = pick(coef.columns, B2_CANDIDATES, args.coef_b2)
    if None in (ir_c, b1, b2):
        sys.exit(f"Could not detect coefficient columns (ir_id={ir_c}, b1={b1}, b2={b2}). "
                 f"Use --ir-id-col/--coef-b1/--coef-b2. Columns: {list(coef.columns)}")
    intc = pick(coef.columns, INTERCEPT_CANDIDATES, args.intercept_col)
    dtmin = pick(coef.columns, DTMIN_CANDIDATES, args.dtmin_col)
    dtmax = pick(coef.columns, DTMAX_CANDIDATES, args.dtmax_col)
    print(f"[coefficients] ir_id={ir_c}  b1={b1}  b2={b2}  intercept={intc}  dt_range=({dtmin},{dtmax})")

    coef_out = pd.DataFrame({
        "impact_region_id": coef[ir_c].astype(str),
        "beta_delta_gmst": coef[b1].astype(float),
        "beta_delta_gmst_squared": coef[b2].astype(float),
    })
    if intc:  coef_out["intercept"] = coef[intc].astype(float)
    if dtmin: coef_out["delta_gmst_fit_min"] = coef[dtmin].astype(float)
    if dtmax: coef_out["delta_gmst_fit_max"] = coef[dtmax].astype(float)
    coef_out = coef_out.drop_duplicates("impact_region_id")
    coef_out.to_csv(os.path.join(args.out, "forest_coefficients.csv"), index=False)
    print(f"[out] forest_coefficients.csv ({len(coef_out)} rows)")

    # ---- weights (baseline area + shares) ----------------------------------
    w = read_any(args.weights, args.weights_layer)
    print(f"[weights] {len(w)} rows; columns: {list(w.columns)}")
    report.append(f"[weights] file={args.weights} rows={len(w)} columns={list(w.columns)}")
    ir_w = pick(w.columns, IR_ID_CANDIDATES, args.ir_id_col)
    ctry = pick(w.columns, COUNTRY_CANDIDATES, args.country_col)
    base = pick(w.columns, BASELINE_CANDIDATES, args.baseline_col)
    tot = pick(w.columns, TOTAL_AREA_CANDIDATES, args.total_area_col)
    ps = pick(w.columns, POP_SHARE_CANDIDATES, args.pop_share_col)
    gs = pick(w.columns, GDP_SHARE_CANDIDATES, args.gdp_share_col)
    print(f"[weights] ir_id={ir_w} country={ctry} baseline={base} total_area={tot} "
          f"pop_share={ps} gdp_share={gs}")
    if ir_w is None or base is None:
        sys.exit("Could not detect impact_region_id and/or baseline area column in weights file. "
                 "Use --ir-id-col / --baseline-col.")

    # baseline area
    barea = pd.DataFrame({"impact_region_id": w[ir_w].astype(str),
                          "baseline_forest_area": w[base].astype(float)})
    if tot:
        barea["total_impact_region_area"] = w[tot].astype(float)
    barea = barea.drop_duplicates("impact_region_id")
    barea.to_csv(os.path.join(args.out, "baseline_forest_area.csv"), index=False)
    print(f"[out] baseline_forest_area.csv ({len(barea)} rows)")

    # shares require a country column
    if ctry is None:
        print("  !! No country column found in weights; cannot validate/write shares.")
        print("     Provide --country-col, or ensure the crosswalk is built separately.")
    else:
        if ps:
            popdf = pd.DataFrame({"impact_region_id": w[ir_w].astype(str),
                                  "give_country_id": w[ctry].astype(str).str.upper(),
                                  "population_share": w[ps].astype(float)})
            popdf, _ = validate_shares(popdf, "impact_region_id", "give_country_id",
                                       "population_share", "population", report)
            popdf.to_csv(os.path.join(args.out, "population_shares.csv"), index=False)
            print(f"[out] population_shares.csv ({len(popdf)} rows)")
        else:
            print("  !! population share column not found; pass --pop-share-col.")
        if gs:
            gdpdf = pd.DataFrame({"impact_region_id": w[ir_w].astype(str),
                                  "give_country_id": w[ctry].astype(str).str.upper(),
                                  "gdp_share": w[gs].astype(float)})
            gdpdf, _ = validate_shares(gdpdf, "impact_region_id", "give_country_id",
                                       "gdp_share", "gdp", report)
            gdpdf.to_csv(os.path.join(args.out, "gdp_shares.csv"), index=False)
            print(f"[out] gdp_shares.csv ({len(gdpdf)} rows)")
        else:
            print("  !! gdp share column not found; pass --gdp-share-col.")

    with open(os.path.join(args.out, "preprocess_report.txt"), "w") as fh:
        fh.write("\n".join(report) + "\n")
    print(f"[out] preprocess_report.txt")
    print("Done. If a country column was present you may still need the polygon-based "
          "crosswalk (build_impact_region_country_crosswalk.py) unless give_country_id "
          "here already uses GIVE ISO3 codes.")


if __name__ == "__main__":
    main()
