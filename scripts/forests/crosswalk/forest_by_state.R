# =============================================================================
#  Forest-area change by Mexican state (INEGI), from impact-region results
# =============================================================================
#
#  Uses the Monte Carlo exports of the forest emulator (forest_ir_<SCENARIO>.csv)
#  and the crosswalk ir_inegi_crosswalk.csv / ir_inegi_region.csv.
#
#  Figure 1  fig_state_boxplot_<scen>_<model>_<year>.png
#            One box per state: distribution of the % change in forest area
#            (median across Monte Carlo trials) of the impact regions whose
#            dominant state it is. Box = 25th–75th percentile, line = median,
#            whiskers = 1.5 × IQR, dots = regions beyond that.
#
#  Figure 2  fig_state_totals_<scen>_<year>.png
#            State totals: forest area change (%) with p05–p95 whiskers, one
#            colour per vegetation model. Regional values are split among
#            states with pct_of_region (area share), assuming forest is spread
#            evenly inside each region.
#
#  Table     state_forest_<scen>_<year>.csv  (state × model totals)
#
#  Symbols
#    F_r,t   forest area of region r in year t (Mha); F_r,2022 is the start
#    Δ%_r    100 · (F_r,year / F_r,2022 − 1)   % change of region r
#    w_r,s   pct_of_region / 100: share of region r's area inside state s
#    F_s,t   Σ_r w_r,s · F_r,t                  forest of state s
#
#  Caveat: state p05/p95 add up regional percentiles, so they are an envelope
#  (wider than the true percentile of the state total).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# --- 0. Settings -------------------------------------------------------------
dir_mc   <- "forest_mcs_output"  # forest_ir_*.csv
dir_cw   <- "data/Forest/crosswalk"                         # crosswalk CSVs
salida   <- "figures/crosswalk"

escenario  <- "SSP245"   # SSP119 SSP126 SSP245 SSP370 SSP585
modelo     <- "jules"    # jules lpjml mc2 (figure 1)
anio       <- 2100       # 2030 … 2100
min_ha     <- 1000       # regions with less forest than this in 2022 are left out of the boxes

dir.create(salida, recursive = TRUE, showWarnings = FALSE)

# --- 1. Read -----------------------------------------------------------------
cols <- c("scenario", "veg_model", "region_id", "iso", "year",
          "F_Mha_mean", "F_Mha_p05", "F_Mha_p50", "F_Mha_p95")
mc <- fread(file.path(dir_mc, paste0("forest_ir_", escenario, ".csv")), select = cols)
mc <- mc[iso == "MEX" & year %in% c(2022L, anio)]

cw     <- fread(file.path(dir_cw, "ir_inegi_crosswalk.csv"), colClasses = list(character = c("CVE_ENT", "CVE_MUN", "CVEGEO")))
region <- fread(file.path(dir_cw, "ir_inegi_region.csv"),    colClasses = list(character = c("CVE_ENT", "CVEGEO_principal")))

# short state names for the axis
corto <- c("Coahuila de Zaragoza" = "Coahuila", "Michoacán de Ocampo" = "Michoacán",
           "Veracruz de Ignacio de la Llave" = "Veracruz")
abrev <- function(x) ifelse(x %in% names(corto), corto[x], x)

# --- 2. Regional change (wide: 2022 vs target year) ---------------------------
base <- mc[year == 2022L, .(veg_model, region_id, F0 = F_Mha_p50)]
fin  <- mc[year == anio,  .(veg_model, region_id, F_p05 = F_Mha_p05, F_p50 = F_Mha_p50, F_p95 = F_Mha_p95)]
reg  <- merge(base, fin, by = c("veg_model", "region_id"))
reg[, cambio_pct := 100 * (F_p50 / F0 - 1)]

missing <- setdiff(unique(reg$region_id), region$region_id)
if (length(missing)) message(length(missing), " regions without crosswalk (left out): ", paste(head(missing, 5), collapse = ", "))

# --- 3. Figure 1: boxplot by state (dominant state of each region) -----------
box <- merge(reg[veg_model == modelo & F0 >= min_ha / 1e6],
             region[, .(region_id, estado)], by = "region_id")
box[, estado := abrev(estado)]
orden <- box[, .(m = median(cambio_pct)), by = estado][order(m), estado]
box[, estado := factor(estado, levels = orden)]

p1 <- ggplot(box, aes(x = cambio_pct, y = estado)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_boxplot(fill = "#CFE3D6", colour = "#2D5A43", outlier.size = 0.8, width = 0.65) +
  labs(
    title = sprintf("Forest-area change by state · %s, %s, %d", escenario, toupper(modelo), anio),
    subtitle = sprintf("Each box: impact regions of the state (median of the Monte Carlo trials).\nRegions with < %s ha of forest in 2022 left out.",
                       format(min_ha, big.mark = ",")),
    x = "Forest area vs 2022 (%)", y = NULL,
    caption = "Source: MimiGIVE_forests.jl forest emulator; impact regions: Climate Impact Lab; states: INEGI Marco Geoestadístico."
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey40", hjust = 0))
ggsave(file.path(salida, sprintf("fig_state_boxplot_%s_%s_%d.png", tolower(escenario), modelo, anio)),
       p1, width = 8, height = 8, dpi = 300, bg = "white")

# --- 4. State totals with p05–p95 whiskers (all models) -----------------------
w <- cw[, .(region_id, CVE_ENT, estado, w = pct_of_region / 100)]
tot <- merge(reg, w, by = "region_id", allow.cartesian = TRUE)
tot <- tot[, .(F0 = sum(w * F0), F_p05 = sum(w * F_p05), F_p50 = sum(w * F_p50), F_p95 = sum(w * F_p95)),
           by = .(veg_model, CVE_ENT, estado)]
tot[, `:=`(chg_p05 = 100 * (F_p05 / F0 - 1), chg_p50 = 100 * (F_p50 / F0 - 1), chg_p95 = 100 * (F_p95 / F0 - 1))]
tot <- tot[F0 > 0]
setorder(tot, CVE_ENT, veg_model)
fwrite(tot[, .(CVE_ENT, estado, veg_model, F0_Mha = round(F0, 4), F_p05_Mha = round(F_p05, 4),
               F_p50_Mha = round(F_p50, 4), F_p95_Mha = round(F_p95, 4),
               chg_p05 = round(chg_p05, 2), chg_p50 = round(chg_p50, 2), chg_p95 = round(chg_p95, 2))],
       file.path(salida, sprintf("state_forest_%s_%d.csv", tolower(escenario), anio)), bom = TRUE)

tot[, estado := abrev(estado)]
orden2 <- tot[veg_model == modelo][order(chg_p50), estado]
tot[, estado := factor(estado, levels = orden2)]
colores <- c(jules = "#0072B2", lpjml = "#D55E00", mc2 = "#009E73")

p2 <- ggplot(tot, aes(y = estado, colour = veg_model)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = chg_p05, xmax = chg_p95), height = 0, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(aes(x = chg_p50), size = 1.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = colores, name = "Vegetation model") +
  labs(
    title = sprintf("State forest area in %d vs 2022 · %s", anio, escenario),
    subtitle = "Point: median; whiskers: 5th–95th percentile (sum of regional percentiles, an envelope).",
    x = "Forest area vs 2022 (%)", y = NULL,
    caption = "Regions split among states by area share (pct_of_region)."
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey40", hjust = 0))
ggsave(file.path(salida, sprintf("fig_state_totals_%s_%d.png", tolower(escenario), anio)),
       p2, width = 8, height = 9, dpi = 300, bg = "white")

cat("Regions in boxplot:", nrow(box), "| states:", uniqueN(box$estado), "\n")
cat("Outputs written to", salida, "\n")
