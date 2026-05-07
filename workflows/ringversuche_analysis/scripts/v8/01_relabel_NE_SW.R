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

message("v8 setup loaded: relabel_ne_sw, apply_v8_theme available.")
