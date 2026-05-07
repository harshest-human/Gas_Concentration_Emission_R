######################################################################
# 01_relabel_NE_SW.R   (v8 wrappers for YBENG-D-26-00430 revision)
# ---------------------------------------------------------------------
# Pure helper functions, no I/O. Sourced by 02..05.
#
# Reviewer R2.18 asks for "north" / "south" to read "NE" / "SW".
# Reviewer R1.12 asks for concentration y-axes clipped at 0.
# Both apply to figures, captions and table headers, NEVER to source
# data files in clean_data/Version_6/.
######################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ---- string-level relabel maps -------------------------------------
v8_loc_map <- c(
  "North background" = "NE background",
  "South background" = "SW background",
  "north background" = "NE background",
  "south background" = "SW background"
)

# ---- relabel column names by suffix --------------------------------
# Renames only suffixes ending in _N or _S so we don't clobber
# columns like "n", "n_dairycows_in", or anything inside a value.
relabel_ne_sw_cols <- function(df) {
  nms  <- names(df)
  nms2 <- ifelse(grepl("_N$", nms), sub("_N$", "_NE", nms), nms)
  nms2 <- ifelse(grepl("_S$", nms2), sub("_S$", "_SW", nms2), nms2)
  names(df) <- nms2
  df
}

# ---- relabel string values in a tibble -----------------------------
# Walks every character / factor column and applies v8_loc_map.
relabel_ne_sw_strings <- function(df) {
  df %>% mutate(across(where(~ is.character(.x) || is.factor(.x)), function(x) {
    if (is.factor(x)) {
      lv <- levels(x)
      lv2 <- ifelse(lv %in% names(v8_loc_map), v8_loc_map[lv], lv)
      factor(as.character(x), levels = lv) -> tmp                 # preserve order
      lv_new <- ifelse(levels(tmp) %in% names(v8_loc_map),
                       v8_loc_map[levels(tmp)], levels(tmp))
      levels(tmp) <- lv_new
      tmp
    } else {
      ifelse(x %in% names(v8_loc_map), v8_loc_map[x], x)
    }
  }))
}

# ---- combined: relabel columns + values ----------------------------
relabel_ne_sw <- function(df) {
  relabel_ne_sw_strings(relabel_ne_sw_cols(df))
}

# ---- v8 cosmetics applied to a ggplot ------------------------------
# clip_y0 = TRUE:  add ylim(0, NA) so negative readings are suppressed
# panel_spacing_pt: row spacing in pt for free-y faceted plots
# Returns the original object if it is not a ggplot (e.g. patchwork).
apply_v8_theme <- function(p, clip_y0 = FALSE, panel_spacing_pt = 20) {
  add_layers <- function(g) {
    if (clip_y0) {
      g <- g + coord_cartesian(ylim = c(0, NA))
    }
    g <- g + theme(panel.spacing.y = unit(panel_spacing_pt, "pt"),
                   panel.spacing.x = unit(8, "pt"))
    g
  }
  if (inherits(p, "patchwork")) {
    # patchwork: apply to every panel
    if (!is.null(p$patches) && !is.null(p$patches$plots)) {
      p$patches$plots <- lapply(p$patches$plots, function(q) {
        if (inherits(q, "ggplot")) add_layers(q) else q
      })
    }
    if (inherits(p, "ggplot")) p <- add_layers(p)
    return(p)
  }
  if (inherits(p, "ggplot")) return(add_layers(p))
  p
}

# ---- helper: rewrite axis / strip / legend labels via relabeller ----
# For plots faceted on a "location" column, applies the NE/SW rename
# to facet strips. Returns a function suitable for facet_wrap(labeller=).
v8_loc_labeller <- function(...) {
  function(value) {
    chr <- as.character(value)
    ifelse(chr %in% names(v8_loc_map), v8_loc_map[chr], chr)
  }
}

# ---- helper: write_excel_csv shim with mkdir -----------------------
v8_write_csv <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_excel_csv(df, path)
  message("v8: wrote ", basename(path), "  (", nrow(df), " rows)")
}

# ---- per-variable parseable math expression --------------------------
# Used by label_parsed for strip text and by scale_*_discrete labels.
v8_var_expr <- c(
  "CO2_mgm3"   = "c[CO[2]]~~'(mg m'^-3*')'",
  "CH4_mgm3"   = "c[CH[4]]~~'(mg m'^-3*')'",
  "NH3_mgm3"   = "c[NH[3]]~~'(mg m'^-3*')'",
  "delta_CO2"  = "Delta*c[CO[2]]~~'(mg m'^-3*')'",
  "delta_CH4"  = "Delta*c[CH[4]]~~'(mg m'^-3*')'",
  "delta_NH3"  = "Delta*c[NH[3]]~~'(mg m'^-3*')'",
  "Q_vent"     = "Q~~'(m'^3*~h^-1*~'LU'^-1*')'",
  "e_CH4_ghLU" = "e[CH[4]]~~'(g h'^-1*~'LU'^-1*')'",
  "e_NH3_ghLU" = "e[NH[3]]~~'(g h'^-1*~'LU'^-1*')'"
)
# Short version (no units) used on the Lin's CCC y-axis tick labels
v8_var_expr_short <- c(
  "CO2_mgm3"   = "c[CO[2]]",
  "CH4_mgm3"   = "c[CH[4]]",
  "NH3_mgm3"   = "c[NH[3]]",
  "delta_CO2"  = "Delta*c[CO[2]]",
  "delta_CH4"  = "Delta*c[CH[4]]",
  "delta_NH3"  = "Delta*c[NH[3]]",
  "Q_vent"     = "Q",
  "e_CH4_ghLU" = "e[CH[4]]",
  "e_NH3_ghLU" = "e[NH[3]]"
)

# ---- v8 boxplot per variable (rows) x location (cols) ---------------
# Box+whisker with jittered raw points, colored by analyzer.
# Row strip text is parsed math expression for the variable + units.
# Returns a ggplot. df_long must have columns: var, value, location,
# analyzer; locations must be a character vector (e.g. NE/SW backgrounds
# +/- "Barn inside"); vars must be the variable names to stack as rows.
v8_box_plot <- function(df_long, vars, locations,
                        analyzer_levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4",
                                            "CRDS.1","CRDS.2","CRDS.3","baseline"),
                        analyzer_colors = c(
                          "FTIR.1"="#1b9e77","FTIR.2"="#d95f02","FTIR.3"="#7570b3",
                          "FTIR.4"="#e7298a","CRDS.1"="#66a61e","CRDS.2"="#e6ab02",
                          "CRDS.3"="#a6761d","baseline"="black")) {
  # Use the parseable expression as the FACTOR LEVEL itself, so
  # labeller = label_parsed reads it directly without an outer
  # as_labeller / default lookup that can blank out the strip text.
  expr_levels <- unname(v8_var_expr[vars])
  d <- df_long %>%
    dplyr::filter(.data$var %in% vars, .data$location %in% locations) %>%
    dplyr::mutate(
      var_expr = factor(unname(v8_var_expr[as.character(var)]),
                        levels = expr_levels),
      location = factor(location, levels = locations),
      analyzer = factor(analyzer, levels = analyzer_levels)
    )
  ggplot(d, aes(x = analyzer, y = value, color = analyzer, fill = analyzer)) +
    geom_jitter(width = 0.25, alpha = 0.45, size = 0.7, shape = 16) +
    geom_boxplot(alpha = 0.10, outlier.shape = NA, linewidth = 0.4) +
    facet_grid(rows = vars(var_expr), cols = vars(location),
               scales = "free_y", switch = "y",
               labeller = labeller(var_expr = label_parsed,
                                   location = label_value)) +
    scale_color_manual(values = analyzer_colors, drop = FALSE) +
    scale_fill_manual(values = analyzer_colors, drop = FALSE) +
    labs(x = NULL, y = NULL, color = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text.y.left = element_text(angle = 90, size = 11),
      strip.text.x      = element_text(size = 11),
      axis.title.y      = element_blank(),
      axis.text.x       = element_text(angle = 30, hjust = 1),
      panel.spacing.y   = unit(20, "pt"),
      panel.spacing.x   = unit(8, "pt"),
      legend.position   = "bottom"
    )
}

message("v8 setup loaded: relabel_ne_sw, apply_v8_theme available.")
