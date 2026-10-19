# =============================================================================
#  Zoom to municipios: forest-area change for one state (INEGI municipios)
# =============================================================================
#
#  Spreads impact-region results to municipios with the crosswalk
#  (area share), then maps them for the chosen state.
#
#    F_m,t = Σ_r  w_r,m · F_r,t      w_r,m = pct_of_region / 100
#    Δ%_m  = 100 · (F_m,year / F_m,2022 − 1)
#
#  Assumes forest is spread evenly inside each impact region, so municipios
#  that share one region get the same % change: the map shows the resolution
#  of the impact regions, drawn on municipio borders.
#
#  Outputs: fig_municipios_<estado>_<scen>_<model>_<year>.png
#           municipio_forest_<estado>_<scen>_<model>_<year>.csv
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(ggplot2)
})

# --- 0. Settings -------------------------------------------------------------
dir_mc   <- "forest_mcs_output"
dir_cw   <- "data/Forest/crosswalk"
ruta_mun <- "path/to/inputs/INEGI/00mun.shp"   # <- change
salida   <- "figures/crosswalk"

cve_ent   <- "27"        # INEGI state code: 27 = Tabasco, 07 = Chiapas, 31 = Yucatán …
escenario <- "SSP245"
modelo    <- "jules"
anio      <- 2100
min_ha    <- 100         # municipios with less forest than this in 2022 shown in grey

dir.create(salida, recursive = TRUE, showWarnings = FALSE)

# --- 1. Read -----------------------------------------------------------------
mc <- fread(file.path(dir_mc, paste0("forest_ir_", escenario, ".csv")),
            select = c("veg_model", "region_id", "iso", "year", "F_Mha_p05", "F_Mha_p50", "F_Mha_p95"))
mc <- mc[iso == "MEX" & veg_model == modelo & year %in% c(2022L, anio)]

cw  <- fread(file.path(dir_cw, "ir_inegi_crosswalk.csv"),
             colClasses = list(character = c("CVE_ENT", "CVE_MUN", "CVEGEO")))
cw  <- cw[CVE_ENT == cve_ent]
estado <- cw$estado[1]

mun <- st_read(ruta_mun, quiet = TRUE)
mun <- mun[mun$CVE_ENT == cve_ent, c("CVEGEO", "NOMGEO")]

# --- 2. Regions -> municipios ------------------------------------------------
reg <- dcast(mc, region_id ~ year, value.var = c("F_Mha_p05", "F_Mha_p50", "F_Mha_p95"))
setnames(reg, c("region_id", "p05_0", "p05_1", "p50_0", "p50_1", "p95_0", "p95_1"))
x <- merge(cw[, .(region_id, CVEGEO, municipio, w = pct_of_region / 100)], reg, by = "region_id")
m <- x[, .(F0 = sum(w * p50_0), F_p05 = sum(w * p05_1), F_p50 = sum(w * p50_1), F_p95 = sum(w * p95_1)),
       by = .(CVEGEO, municipio)]
m[, `:=`(chg_p50 = 100 * (F_p50 / F0 - 1), chg_p05 = 100 * (F_p05 / F0 - 1), chg_p95 = 100 * (F_p95 / F0 - 1))]
m[F0 * 1e6 < min_ha, c("chg_p05", "chg_p50", "chg_p95") := NA_real_]

tag <- sprintf("%s_%s_%s_%d", gsub(" ", "_", tolower(estado)), tolower(escenario), modelo, anio)
fwrite(m[order(chg_p50), .(CVEGEO, municipio, forest_2022_ha = round(F0 * 1e6), forest_ha = round(F_p50 * 1e6),
                          chg_p05 = round(chg_p05, 2), chg_p50 = round(chg_p50, 2), chg_p95 = round(chg_p95, 2))],
       file.path(salida, paste0("municipio_forest_", tag, ".csv")), bom = TRUE)

# --- 3. Map --------------------------------------------------------------------
mapa <- merge(mun, m[, .(CVEGEO, chg_p50)], by = "CVEGEO", all.x = TRUE)
lim  <- max(5, ceiling(max(abs(m$chg_p50), na.rm = TRUE)))
p <- ggplot(mapa) +
  geom_sf(aes(fill = chg_p50), colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#A4463A", mid = "#EFE7DA", high = "#2D7652", midpoint = 0,
                       limits = c(-lim, lim), na.value = "#C9CDC9", name = "Forest vs 2022 (%)") +
  labs(title = sprintf("Forest-area change by municipio · %s", estado),
       subtitle = sprintf("%s, %s, %d. Median of Monte Carlo trials. Grey: < %d ha of forest in 2022.",
                          escenario, toupper(modelo), anio, min_ha),
       caption = "Impact-region results spread to municipios by area share (INEGI Marco Geoestadístico; Climate Impact Lab regions).") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), plot.caption = element_text(colour = "grey40", hjust = 0),
        legend.position = "right")
ggsave(file.path(salida, paste0("fig_municipios_", tag, ".png")), p, width = 9, height = 6.5, dpi = 300, bg = "white")

cat(estado, ": ", nrow(m), " municipios | median change ", round(median(m$chg_p50, na.rm = TRUE), 1), " %\n", sep = "")
