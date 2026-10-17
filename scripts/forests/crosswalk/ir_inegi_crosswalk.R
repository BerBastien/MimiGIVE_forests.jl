# =============================================================================
#  Crosswalk: Climate Impact Lab impact regions  <->  INEGI states and municipios
# =============================================================================
#
#  What it does
#    Intersects the impact regions of Mexico with INEGI's Marco Geoestadístico
#    (00ent = states, 00mun = municipios) and records, for every piece of
#    overlap, how much of the region and of the municipio it covers.
#
#  Inputs (change the paths in section 0)
#    - Impact regions: the original Climate Impact Lab shapefile (default) or
#      the simplified dashboard layer; only the id, country and geometry are used.
#    - INEGI 00ent.shp and 00mun.shp (Marco Geoestadístico, CRS: MEXICO ITRF2008 LCC).
#
#  Outputs (folder set in section 0)
#    ir_inegi_crosswalk.csv   one row per region × municipio overlap
#    ir_inegi_region.csv      one row per impact region, with its dominant
#                             state and municipio
#    ir_inegi_municipio.csv   one row per municipio, with its dominant region
#    ir_inegi_check_map.png   map to verify the result
#
#  Column dictionary (ir_inegi_crosswalk.csv)
#    region_id        impact region code (CIL)
#    CVE_ENT          INEGI state code (2 digits, text)
#    estado           state name
#    CVE_MUN          INEGI municipio code within the state (3 digits, text)
#    CVEGEO           INEGI geostatistical key = CVE_ENT + CVE_MUN (5 digits)
#    municipio        municipio name
#    area_km2         area of the overlap (km²)
#    pct_of_region    % of the impact region's area inside this municipio
#    pct_of_municipio % of the municipio's area inside this impact region
#
#  How to use it to aggregate results:
#    value by municipio = Σ over regions ( value of region × pct_of_region / 100 )
#    (area-weighted; assumes the value is spread evenly inside each region)
#
#  Areas are computed in INEGI's projection (EPSG:6372, Lambert conformal
#  conic, metres), so they are true areas, not degrees.
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(ggplot2)
})

# --- 0. Paths ----------------------------------------------------------------
# Original CIL shapefile (recommended, full resolution, 24,378 regions):
#   https://github.com/ClimateImpactLab/data-impact-region-shapefile
#   CIL_impact_regions/agglomerated-world-new.shp (+ .shx, .dbf); columns hierid, ISO
# Run from the repository root. Inputs are large, so they stay outside the repo:
dir_inputs <- "path/to/inputs"   # <- change: folder with the CIL and INEGI files
ruta_ir  <- file.path(dir_inputs, "CIL_impact_regions/agglomerated-world-new.shp")
capa_ir  <- NULL          # NULL for a shapefile; "impact_regions" for the dashboard .gpkg
col_id   <- "hierid"      # "region_id" for the dashboard .gpkg
col_iso  <- "ISO"         # "iso" for the dashboard .gpkg
ruta_ent <- file.path(dir_inputs, "INEGI/00ent.shp")
ruta_mun <- file.path(dir_inputs, "INEGI/00mun.shp")
salida   <- "data/Forest/crosswalk"

dir.create(salida, recursive = TRUE, showWarnings = FALSE)
sf_use_s2(FALSE)                    # planar geometry: we work in a projected CRS
crs_mx <- 6372                      # Mexico ITRF2008 / LCC (INEGI), metres

# --- 1. Read and prepare -----------------------------------------------------
ir <- if (is.null(capa_ir)) st_read(ruta_ir, quiet = TRUE) else st_read(ruta_ir, layer = capa_ir, quiet = TRUE)
ir <- ir[ir[[col_iso]] == "MEX", c(col_id, col_iso)]
names(ir)[1:2] <- c("region_id", "iso")
if (is.na(st_crs(ir))) st_crs(ir) <- 4326   # the CIL shapefile ships without .prj (WGS84 lon/lat)
ir <- st_make_valid(st_transform(ir, crs_mx))

ent <- st_read(ruta_ent, quiet = TRUE)
mun <- st_read(ruta_mun, quiet = TRUE)
ent <- st_make_valid(st_transform(ent, crs_mx))
mun <- st_make_valid(st_transform(mun, crs_mx))

# state names onto municipios
mun <- merge(mun[, c("CVEGEO", "CVE_ENT", "CVE_MUN", "NOMGEO")],
             st_drop_geometry(ent)[, c("CVE_ENT", "NOMGEO")],
             by = "CVE_ENT", suffixes = c("_mun", "_ent"))
names(mun)[names(mun) == "NOMGEO_mun"] <- "municipio"
names(mun)[names(mun) == "NOMGEO_ent"] <- "estado"

ir$area_region_km2  <- as.numeric(st_area(ir)) / 1e6
mun$area_mun_km2    <- as.numeric(st_area(mun)) / 1e6

# --- 2. Intersect ------------------------------------------------------------
x <- st_intersection(ir[, c("region_id", "area_region_km2")],
                     mun[, c("CVE_ENT", "estado", "CVE_MUN", "CVEGEO", "municipio", "area_mun_km2")])
x$area_km2 <- as.numeric(st_area(x)) / 1e6

cw <- as.data.table(st_drop_geometry(x))
cw <- cw[area_km2 > 0.01]           # drop slivers < 1 ha from border mismatch
cw[, pct_of_region    := 100 * area_km2 / area_region_km2]
cw[, pct_of_municipio := 100 * area_km2 / area_mun_km2]
setorder(cw, region_id, -area_km2)
crosswalk <- cw[, .(region_id, CVE_ENT, estado, CVE_MUN, CVEGEO, municipio,
                    area_km2 = round(area_km2, 3),
                    pct_of_region = round(pct_of_region, 3),
                    pct_of_municipio = round(pct_of_municipio, 3))]
fwrite(crosswalk, file.path(salida, "ir_inegi_crosswalk.csv"), bom = TRUE)

# --- 3. One row per impact region (dominant state and municipio) ------------
by_state <- cw[, .(area_km2 = sum(area_km2)), by = .(region_id, CVE_ENT, estado, area_region_km2)]
by_state[, pct := 100 * area_km2 / area_region_km2]
dom_state <- by_state[order(-pct), .SD[1], by = region_id]
n_states  <- by_state[, .(n_estados = .N), by = region_id]
dom_mun   <- cw[order(-pct_of_region), .SD[1], by = region_id]
n_muns    <- cw[, .(n_municipios = .N, pct_cubierto = sum(pct_of_region)), by = region_id]

region <- Reduce(function(a, b) merge(a, b, by = "region_id"), list(
  dom_state[, .(region_id, CVE_ENT, estado, pct_estado = round(pct, 2))],
  n_states,
  dom_mun[, .(region_id, CVEGEO_principal = CVEGEO, municipio_principal = municipio,
              pct_municipio_principal = round(pct_of_region, 2))],
  n_muns[, .(region_id, n_municipios, pct_cubierto = round(pct_cubierto, 2))]
))
fwrite(region, file.path(salida, "ir_inegi_region.csv"), bom = TRUE)

# --- 4. One row per municipio (dominant impact region) -----------------------
municipio <- cw[order(-pct_of_municipio), .SD[1], by = CVEGEO][
  , .(CVEGEO, CVE_ENT, estado, CVE_MUN, municipio, region_id_principal = region_id,
      pct_region_principal = round(pct_of_municipio, 2))]
n_reg <- cw[, .(n_regiones = .N), by = CVEGEO]
municipio <- merge(municipio, n_reg, by = "CVEGEO")
fwrite(municipio, file.path(salida, "ir_inegi_municipio.csv"), bom = TRUE)

# --- 5. Checks -----------------------------------------------------------------
cat("Impact regions of Mexico:          ", nrow(ir), "\n")
cat("Regions with a match:              ", uniqueN(cw$region_id), "\n")
cat("Municipios with a match:           ", uniqueN(cw$CVEGEO), "of", nrow(mun), "\n")
cat("Region coverage by municipios (%): median",
    round(median(region$pct_cubierto), 2), "| min", round(min(region$pct_cubierto), 2), "\n")
cat("Regions inside a single state:     ", sum(region$n_estados == 1), "\n")
cat("Regions crossing state lines:      ", sum(region$n_estados > 1), "\n")
cat("Municipios split among >1 regions: ", sum(municipio$n_regiones > 1), "\n")

# --- 6. Check map: each region coloured by its dominant state -----------------
ir_map <- merge(ir, region[, .(region_id, estado)], by = "region_id")
p <- ggplot() +
  geom_sf(data = ir_map, aes(fill = estado), colour = "white", linewidth = 0.05, show.legend = FALSE) +
  geom_sf(data = ent, fill = NA, colour = "black", linewidth = 0.25) +
  labs(title = "Impact regions coloured by dominant INEGI state",
       subtitle = "Black lines: state borders (Marco Geoestadístico)") +
  theme_void()
ggsave(file.path(salida, "ir_inegi_check_map.png"), p, width = 9, height = 6, dpi = 200, bg = "white")
cat("Outputs written to", salida, "\n")
