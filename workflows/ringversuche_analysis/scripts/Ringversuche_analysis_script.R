################################################################################
##  Ringversuche analysis ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â Version 10
##  Builds on V9 (scripts/Ringversuche_analysis_script.R) and addresses the
##  first-round reviewer/editor comments. Each section header lists the
##  comment IDs it answers (see First_review/.../reviewer_comments_shortened.docx).
##
##  Key V10 changes vs V9:
##    * MIN_DELTA_CO2 guard on Q_vent ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â drops divide-by-near-zero spikes.
##    * FTIR.2 NH3 offset retained for diagnosis at concentration level because
##      it largely cancels in deltas, Q and e.
##    * Paired NE vs SW outdoor-line test, then conditioned on wind sector
##      relative to neighbour sources (W dairy, SW feed store, N/NW manure
##      tank + biogas plant) ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â used to drop the contamination-prone line.
##    * Downstream delta/Q/e analysis runs on the single retained outdoor.
##    * Pairwise Deming + Lin's CCC extended to concentration c and ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc
##      (was Q/e only in V9).
##    * Tukey HSD pairwise tables (ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§16): Tables 3, 4, 6 of the manuscript ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â
##      concentration c, ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc and derived quantities, using Tukey-adjusted pairwise
##      p-values within each panel
##      and ns/*/**/*** significance code.
##    * Delta plots no longer floor at 0 (sign carries meaning); concentration
##      plots still floor.
##    * Sampling-cycle figure for ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§2.2 (M&M).
##    * Text reports written to result_data/text_reports/Version_10 so the
##      manuscript can pull numbers directly.
################################################################################

#### Libraries                                                                ####
library(tidyverse)
library(lubridate)
library(scales)
library(patchwork)
# DescTools::CCC used inline.

#### Shared labels and aesthetics  (unchanged from V9)                        ####
VAR_LABELS_UNITS <- c(
        "CO2_mgm3"    = "c[CO2]~'(mg '*m^-3*')'",
        "CH4_mgm3"    = "c[CH4]~'(mg '*m^-3*')'",
        "NH3_mgm3"    = "c[NH3]~'(mg '*m^-3*')'",
        "r_CH4/CO2"   = "c[CH4]/c[CO2]~'('*'%'*')'",
        "r_NH3/CO2"   = "c[NH3]/c[CO2]~'('*'%'*')'",
        "delta_CO2"   = "Delta*c[CO2]~'(mg '*m^-3*')'",
        "delta_CH4"   = "Delta*c[CH4]~'(mg '*m^-3*')'",
        "delta_NH3"   = "Delta*c[NH3]~'(mg '*m^-3*')'",
        "Q_vent"      = "Q~'('*m^3~h^-1~LU^-1*')'",
        "e_CH4_gh"    = "e[CH4]~'(g '*h^-1*')'",
        "e_NH3_gh"    = "e[NH3]~'(g '*h^-1*')'",
        "e_CH4_ghLU"  = "e[CH4]~'(g '*h^-1~LU^-1*')'",
        "e_NH3_ghLU"  = "e[NH3]~'(g '*h^-1~LU^-1*')'",
        "temp"        = "Temperature~(degree*C)",
        "RH"          = "Relative~Humidity~('%')",
        "ws_mst"      = "Wind~Speed~(m~s^-1)",
        "wd_mst"      = "Wind~Direction~(degree)",
        "n_dairycows" = "'Number of Cows'"
)
VAR_LABELS_PLAIN <- c(
        "CO2_mgm3"    = "c[CO2]",
        "CH4_mgm3"    = "c[CH4]",
        "NH3_mgm3"    = "c[NH3]",
        "delta_CO2"   = "Delta*c[CO2]",
        "delta_CH4"   = "Delta*c[CH4]",
        "delta_NH3"   = "Delta*c[NH3]",
        "Q_vent"      = "Q",
        "e_CH4_ghLU"  = "e[CH4]",
        "e_NH3_ghLU"  = "e[NH3]"
)
LOC_LABELS <- c("Indoor" = "Indoor",
                "Outdoor_NE" = "Outdoor^NE",
                "Outdoor_SW" = "Outdoor^SW")
LOC_DATASET_LABELS <- c("Outdoor_NE" = "Outdoor^NE~dataset",
                        "Outdoor_SW" = "Outdoor^SW~dataset")
loc_from_suffix <- c("in" = "Indoor", "N" = "Outdoor_NE", "S" = "Outdoor_SW")

ANALYSER_COLORS <- c(
        "FTIR.1"="#1b9e77","FTIR.2"="#d95f02","FTIR.3"="#7570b3",
        "FTIR.4"="#e7298a","FTIR.4_old"="#e41a1c",
        "CRDS.1"="#66a61e","CRDS.2"="#e6ab02","CRDS.3"="#a6761d"
)
ANALYSER_SHAPES <- c(
        "FTIR.1"=0,"FTIR.2"=1,"FTIR.3"=2,"FTIR.4"=5,"FTIR.4_old"=13,
        "CRDS.1"=15,"CRDS.2"=19,"CRDS.3"=17,"Input"=16
)
ANALYSER_LEGEND_ORDER <- c("CRDS.1", "CRDS.2", "CRDS.3",
                           "FTIR.1", "FTIR.2", "FTIR.3", "FTIR.4", "FTIR.4_old")
analyser_aes <- function(x, lookup, default) {
        present <- unique(as.character(x))
        out     <- setNames(rep(default, length(present)), present)
        known   <- intersect(present, names(lookup))
        out[known] <- lookup[known]
        out
}

#### V10 ANALYSIS CONFIG                                                       ####
# CO2-balance method needs a non-trivial indoor-outdoor CO2 gradient to be
# well-conditioned. Below this floor, Q = PCO2/dCO2 blows up (V9 produced
# Q ~ 6e7 m^3/h/LU spikes ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â see V9 q_e_trend_plot.png). 50 mg/m^3 ~= 25 ppm,
# i.e. several sigma above analyser noise. R2.24 (Reviewer #2) requires this
# kind of plausibility guard.
MIN_DELTA_CO2 <- 50      # mg m-3
MAX_Q_VENT    <- 5000    # m^3 h-1 LU-1, physical sanity ceiling for NVDB
MIN_Q_VENT    <- 50      # m^3 h-1 LU-1, physical sanity floor

# Reviewer 2 (R2.7) ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â FTIR.2 had a multiplicative-style offset on NH3
# NH3 only. The offset cancels in differences (delta NH3 -> Q -> e_NH3), so
# FTIR.2 showed an NH3 concentration offset. FTIR.4_old is retained only as a
# diagnostic series after Section 14 (R2.7 closure).

# Neighbour-source geometry of the LVAT site (manuscript Fig. 1):
#   W       -> additional dairy barns (cow emissions: CO2, CH4, NH3)
#   SW      -> COVERED feed-storage facility (negligible gas-phase emission)
#   N, NW   -> open manure tank + biogas plant (NH3, CH4)
#   NE,E,SE,S -> open land (clean reference)
#
# V10b revision (2026-05-28). Empirical observation across the campaign:
# Outdoor_NE concentrations were systematically HIGHER than Outdoor_SW for
# all three gases (CO2 +1-3 %, CH4 +17-45 %, NH3 +1-100 %), including the
# easterly sectors that the v3 geometry classed as "NE clean". This is
# inconsistent with the v3 geometry, which had NE clean under E/NE wind
# and SW downwind of LVAT under the same winds ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â under that scheme SW
# should have been HIGHER, not lower.
#
# Reconciliation: the LVAT plume reaches the NE sampling line under all
# southerly and most easterly sectors (ridge venting + side-curtain
# dispersion, not just the simple lee side of the barn), and the SW line
# enjoys near-clean status under the southerly and westerly sectors that
# v3 (incorrectly) coded as contaminated. The covered feed-storage facility
# SW of the barn is not a gas-phase source. Revised assignment:
#   NE line downwind of LVAT under: S, SW, W, SE, E (5 sectors)
#   SW line downwind of LVAT under: N, NE, E         (3 sectors, unchanged)
#   NE neighbour-source under:      N, NW             (unchanged)
#   SW neighbour-source under:      (none ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â feed store covered, W dairies
#                                    are upwind of the open WSW-NE land)
#   NE truly clean under:           NE only
#   SW truly clean under:           S, SE, SW, W
# Under this geometry the ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§11 auto-decision selects Outdoor_SW (lower
# contam_time) for the present E-prevailing campaign, which also matches
# the empirical SW < NE concentration ordering.
NE_CONTAM_LVAT     <- c("S","SW","W","SE","E")
SW_CONTAM_LVAT     <- c("N","NE","E")
NE_CONTAM_NEIGHBOUR <- c("N","NW")   # manure tank / biogas reach NE line
SW_CONTAM_NEIGHBOUR <- c()           # covered feed store; no chemical source
# Truly clean (open-land upwind) sectors for each line:
NE_CLEAN_SECTORS <- c("NE")
SW_CLEAN_SECTORS <- c("S","SE","SW","W")

#### Data preparation                                                         ####
# Indirect CO2 balance ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â V10 adds MIN_DELTA_CO2 guard and drops negative
# delta_CO2 (outdoor > indoor is physically inconsistent with the balance).
indirect.CO2.balance <- function(df, min_delta_co2 = MIN_DELTA_CO2,
                                  q_lo = MIN_Q_VENT, q_hi = MAX_Q_VENT) {
        ppm_to_mgm3 <- function(ppm, molar_mass) {
                T_K <- 273.15; P <- 101325; R <- 8.314472
                (ppm * 1e-6) * molar_mass * 1e3 * P / (R * T_K)
        }
        df <- df %>%
                mutate(
                        hour      = as.numeric(format(DATE.TIME, "%H")),
                        a         = 0.22,
                        h_min     = 2.9,
                        phi       = 5.6 * m_weight^0.75 + 22 * Y1_milk_prod + 1.6e-5 * p_pregnancy_day^3,
                        t_factor  = 1 + 4e-5 * (20 - temp_in)^3,
                        phi_T_cor = phi * t_factor,
                        A_cor     = 1 - a * sin((2*pi/24) * (hour + 6 - h_min)),
                        hpu_T_A_cor_per_cow = phi_T_cor * A_cor,
                        PCO2      = (0.185 * hpu_T_A_cor_per_cow) * 1000,

                        CO2_mgm3_in = ppm_to_mgm3(CO2_ppm_in, 44.01),
                        CO2_mgm3_N  = ppm_to_mgm3(CO2_ppm_N,  44.01),
                        CO2_mgm3_S  = ppm_to_mgm3(CO2_ppm_S,  44.01),
                        NH3_mgm3_in = ppm_to_mgm3(NH3_ppm_in, 17.031),
                        NH3_mgm3_N  = ppm_to_mgm3(NH3_ppm_N,  17.031),
                        NH3_mgm3_S  = ppm_to_mgm3(NH3_ppm_S,  17.031),
                        CH4_mgm3_in = ppm_to_mgm3(CH4_ppm_in, 16.04),
                        CH4_mgm3_N  = ppm_to_mgm3(CH4_ppm_N,  16.04),
                        CH4_mgm3_S  = ppm_to_mgm3(CH4_ppm_S,  16.04),

                        delta_CO2_N = CO2_mgm3_in - CO2_mgm3_N,
                        delta_CO2_S = CO2_mgm3_in - CO2_mgm3_S,
                        delta_NH3_N = NH3_mgm3_in - NH3_mgm3_N,
                        delta_NH3_S = NH3_mgm3_in - NH3_mgm3_S,
                        delta_CH4_N = CH4_mgm3_in - CH4_mgm3_N,
                        delta_CH4_S = CH4_mgm3_in - CH4_mgm3_S,

                        # MIN_DELTA_CO2 guard: only compute Q where the CO2
                        # gradient is well above analyser noise.
                        Q_vent_N = ifelse(delta_CO2_N >= min_delta_co2,
                                          PCO2 / delta_CO2_N, NA_real_),
                        Q_vent_S = ifelse(delta_CO2_S >= min_delta_co2,
                                          PCO2 / delta_CO2_S, NA_real_),

                        # V10b: the [q_lo, q_hi] physical-sanity bound was
                        # removed because it silently dropped rows without
                        # documenting the loss; it is replaced by an explicit
                        # per-analyser 1.5*IQR filter on Q_vent applied in
                        # ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§11.5 with a dropout-count CSV.

                        e_NH3_gh_N = (delta_NH3_N * Q_vent_N / 1000) * n_dairycows_in,
                        e_CH4_gh_N = (delta_CH4_N * Q_vent_N / 1000) * n_dairycows_in,
                        e_NH3_gh_S = (delta_NH3_S * Q_vent_S / 1000) * n_dairycows_in,
                        e_CH4_gh_S = (delta_CH4_S * Q_vent_S / 1000) * n_dairycows_in,

                        e_NH3_ghLU_N = (e_NH3_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_N = (e_CH4_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_NH3_ghLU_S = (e_NH3_gh_S * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_S = (e_CH4_gh_S * 500) / (n_dairycows_in * m_weight)
                )
        for (gas in c("CH4", "NH3")) {
                for (loc in c("in", "N", "S")) {
                        col_gas   <- paste0(gas, "_mgm3_", loc)
                        col_CO2   <- paste0("CO2_mgm3_", loc)
                        ratio_col <- paste0("r_", gas, "/CO2_", loc)
                        df <- df %>%
                                mutate(!!ratio_col := ifelse(.data[[col_CO2]] != 0,
                                                             .data[[col_gas]] / .data[[col_CO2]] * 100,
                                                             NA_real_))
                }
        }
        df
}

# reshaper
reshaper <- function(df) {
        meta_cols <- c("DATE.TIME", "analyser")
        measure_cols <- df %>% select(-all_of(meta_cols)) %>%
                select(where(is.numeric)) %>% names()

        df_long <- df %>%
                pivot_longer(cols = all_of(measure_cols),
                             names_to = "var", values_to = "value") %>%
                mutate(
                        location = case_when(
                                str_detect(var, "_in$") ~ "Indoor",
                                str_detect(var, "_N$")  ~ "Outdoor_NE",
                                str_detect(var, "_S$")  ~ "Outdoor_SW",
                                TRUE ~ NA_character_
                        ),
                        var = str_remove(var, "_(in|N|S)$"),
                        DATE.TIME = as.POSIXct(DATE.TIME),
                        day  = factor(as.Date(DATE.TIME)),
                        hour = factor(format(DATE.TIME, "%H:%M"))
                ) %>%
                select(DATE.TIME, day, hour, location, analyser, var, value) %>%
                arrange(DATE.TIME, var, analyser, location) %>%
                mutate(analyser = case_when(
                        var %in% c("temp", "RH")                           ~ "HOBO",
                        var %in% c("wd_mst", "ws_mst", "wd_trv", "ws_trv") ~ "USA",
                        var %in% c("n_dairycows")                          ~ "RGB",
                        TRUE                                               ~ analyser
                ))

        df_long <- df_long %>%
                mutate(analyser = factor(analyser,
                                         levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4","FTIR.4_old",
                                                    "CRDS.1","CRDS.2","CRDS.3",
                                                    "HOBO","USA","RGB"))) %>%
                arrange(DATE.TIME, location, var)
        df_long
}

#### Statistics helpers                                                       ####
deg_to_compass8 <- function(deg) {
        labs <- c("N","NE","E","SE","S","SW","W","NW")
        factor(labs[floor(((deg %% 360) + 22.5) / 45) %% 8 + 1], levels = labs)
}
deming_fit <- function(x, y) {
        ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
        if (length(x) < 3) return(c(intercept = NA_real_, slope = NA_real_))
        sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
        slope <- (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy)
        c(intercept = mean(y) - slope * mean(x), slope = slope)
}
ba_relative_stats <- function(x, y) {
        ok <- is.finite(x) & is.finite(y) & ((x + y) != 0)
        x <- x[ok]; y <- y[ok]; n <- length(x)
        if (n < 5L) return(NULL)
        d <- 100 * (y - x) / ((x + y) / 2)
        bias <- mean(d); sd_d <- sd(d); se_b <- sd_d / sqrt(n)
        tibble(n = n, bias_pct = bias,
               bias_lo95 = bias - 1.96 * se_b, bias_hi95 = bias + 1.96 * se_b,
               loa_low = bias - 1.96 * sd_d, loa_high = bias + 1.96 * sd_d,
               sd_pct = sd_d)
}

#### Plotting functions  (V10: signed axes for deltas)                        ####
# delta vars carry sign; concentration values and Q/e are non-negative.
.is_signed_var <- function(v) any(grepl("^delta_", v))

emitrendplot <- function(data, y = NULL, location_filter = NULL,
                         plot_err = FALSE, x = "DATE.TIME") {
        if (!is.null(location_filter)) data <- data %>% filter(location %in% location_filter)
        if (!is.null(y))               data <- data %>% filter(var %in% y)
        value_col <- if (plot_err && "pct_err" %in% names(data)) "pct_err" else "value"
        summary_data <- data %>%
                group_by(.data[[x]], analyser, location, var) %>%
                summarise(mean_val = mean(.data[[value_col]], na.rm = TRUE),
                          sd_val   = sd(.data[[value_col]], na.rm = TRUE),
                          .groups  = "drop") %>%
                filter(!is.na(mean_val))
        if (nrow(summary_data) == 0) stop("No data for emitrendplot")
        labels <- if (plot_err) VAR_LABELS_PLAIN else VAR_LABELS_UNITS
        summary_data <- summary_data %>%
                mutate(facet_label = factor(labels[as.character(var)], levels = labels[y]))
        col_vals   <- analyser_aes(summary_data$analyser, ANALYSER_COLORS, "black")
        shape_vals <- analyser_aes(summary_data$analyser, ANALYSER_SHAPES, 16)
        signed <- .is_signed_var(y)
        p <- ggplot(summary_data,
                    aes(x = .data[[x]], y = mean_val,
                        color = analyser, shape = analyser, group = analyser)) +
                geom_line(linewidth = 0.5, alpha = 0.6, na.rm = TRUE) +
                geom_point(size = 1.1, na.rm = TRUE) +
                geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                              width = 0.25, linewidth = 0.35, na.rm = TRUE) +
                scale_color_manual(values = col_vals) +
                scale_shape_manual(values = shape_vals) +
                # When only one location is shown (e.g. retained outdoor only),
                # drop the column facet entirely so the redundant strip header
                # disappears. Otherwise face by location as usual.
                (if (length(unique(summary_data$location)) <= 1)
                        facet_grid(facet_label ~ ., scales = "free_y", switch = "y",
                                   labeller = labeller(facet_label = label_parsed))
                 else
                        facet_grid(facet_label ~ location, scales = "free_y", switch = "y",
                                   labeller = labeller(
                                           facet_label = label_parsed,
                                           location    = as_labeller(LOC_LABELS, label_parsed)))) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                                   labels = scales::label_number(accuracy = 0.1, big.mark = "")) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
                      axis.text.y = element_text(hjust = 1, size = 12),
                      strip.text.y.left = element_text(angle = 0, size = 13, vjust = 0.5),
                      strip.text.x = element_text(size = 12),
                      panel.border = element_rect(color = "black", fill = NA),
                      legend.position = "bottom", legend.title = element_blank(),
                      plot.title = element_text(hjust = 0.5)) +
                guides(color = guide_legend(nrow = 1))
        if (signed) p <- p + geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3)
        else        p <- p + coord_cartesian(ylim = c(0, NA))
        if (x == "DATE.TIME") {
                x_breaks <- seq(floor_date(min(summary_data$DATE.TIME, na.rm = TRUE), unit = "day"),
                                ceiling_date(max(summary_data$DATE.TIME, na.rm = TRUE), unit = "day"),
                                by = "1 day")
                p <- p + scale_x_datetime(breaks = x_breaks, date_labels = "%d-%b")
        }
        p
}

emiboxplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE) {
        if (!is.null(location_filter)) data <- data %>% filter(location %in% location_filter)
        if (!is.null(y))               data <- data %>% filter(var %in% y)
        value_col <- if (plot_err) "pct_err" else "value"
        data <- data %>%
                group_by(location, analyser, var) %>%
                mutate(Q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
                       Q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
                       IQR = Q3 - Q1,
                       is_outlier = (.data[[value_col]] < (Q1 - 1.5 * IQR)) |
                                    (.data[[value_col]] > (Q3 + 1.5 * IQR))) %>%
                filter(!is_outlier) %>% ungroup() %>%
                select(-Q1, -Q3, -IQR, -is_outlier)
        labels <- if (plot_err) VAR_LABELS_PLAIN else VAR_LABELS_UNITS
        data <- data %>%
                mutate(variable_label = factor(labels[as.character(var)], levels = labels[y]))
        col_vals <- analyser_aes(data$analyser, ANALYSER_COLORS, "black")
        signed <- .is_signed_var(y)
        p <- ggplot(data, aes(x = analyser, y = .data[[value_col]], color = analyser)) +
                geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.6) +
                geom_jitter(width = 0.12, alpha = 0.22, size = 0.8) +
                scale_color_manual(values = col_vals) +
                # Hide the redundant column strip when only one location is shown.
                (if (length(unique(data$location)) <= 1)
                        facet_grid(variable_label ~ ., scales = "free_y", switch = "y",
                                   labeller = labeller(variable_label = label_parsed))
                 else
                        facet_grid(variable_label ~ location, scales = "free_y", switch = "y",
                                   labeller = labeller(
                                           variable_label = label_parsed,
                                           location       = as_labeller(LOC_LABELS, label_parsed)))) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                                   labels = scales::label_number(accuracy = 0.1, big.mark = "")) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
                      axis.text.y = element_text(size = 12),
                      strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 13),
                      strip.text.x = element_text(size = 12),
                      panel.border = element_rect(color = "black", fill = NA),
                      legend.position = "bottom", legend.title = element_blank(),
                      plot.title = element_text(hjust = 0.5)) +
                guides(color = guide_legend(nrow = 1))
        if (signed) p <- p + geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3)
        else        p <- p + coord_cartesian(ylim = c(0, NA))
        p
}

hourly_dataset_summary_plot <- function(emission_csv, output_file, var_order,
                                        include_ftir4_old = FALSE) {
        hour_order <- 0:23

        raw_emission <- readr::read_csv(emission_csv, show_col_types = FALSE)
        if (!include_ftir4_old) {
                raw_emission <- raw_emission %>% filter(analyser != "FTIR.4_old")
        }

        derived_long <- raw_emission %>%
                mutate(hour_num = as.integer(hour)) %>%
                pivot_longer(
                        cols = matches(paste0("^(", paste(var_order, collapse = "|"), ")_(N|S)$")),
                        names_to = c("var", "suffix"),
                        names_pattern = "^(.*)_([NS])$",
                        values_to = "value"
                ) %>%
                mutate(
                        dataset = recode(suffix, "N" = "Outdoor_NE", "S" = "Outdoor_SW"),
                        value = as.character(value)
                ) %>%
                select(hour_num, analyser, dataset, var, value)

        shared_inputs <- readr::read_csv(emission_csv, show_col_types = FALSE) %>%
                mutate(hour_num = as.integer(hour)) %>%
                distinct(DATE.TIME, hour_num, n_dairycows_in, temp_N, wd_mst, ws_mst) %>%
                rename(n_dairycows = n_dairycows_in, temp = temp_N) %>%
                tidyr::crossing(dataset = c("Outdoor_NE", "Outdoor_SW")) %>%
                mutate(analyser = "Input") %>%
                pivot_longer(
                        cols = c(n_dairycows, temp, wd_mst, ws_mst),
                        names_to = "var",
                        values_to = "value"
                ) %>%
                mutate(value = as.character(value)) %>%
                select(hour_num, analyser, dataset, var, value)

        all_long <- bind_rows(derived_long, shared_inputs) %>%
                filter(var %in% var_order, hour_num %in% hour_order)

        numeric_summary <- all_long %>%
                mutate(value_num = suppressWarnings(as.numeric(value))) %>%
                filter(is.finite(value_num)) %>%
                group_by(dataset, var, analyser, hour_num) %>%
                summarise(
                        mean_val = mean(value_num, na.rm = TRUE),
                        se_val   = sd(value_num, na.rm = TRUE) / sqrt(sum(is.finite(value_num))),
                        .groups  = "drop"
                )

        wd_summary <- tibble()

        summary_data <- bind_rows(numeric_summary, wd_summary) %>%
                mutate(
                        hour_lab = sprintf("%02d:00", hour_num),
                        var = factor(var, levels = var_order),
                        facet_label = factor(VAR_LABELS_UNITS[as.character(var)],
                                             levels = VAR_LABELS_UNITS[var_order]),
                        dataset_label = factor(LOC_DATASET_LABELS[dataset],
                                               levels = LOC_DATASET_LABELS[c("Outdoor_NE", "Outdoor_SW")]),
                        analyser = factor(analyser,
                                          levels = c(ANALYSER_LEGEND_ORDER, "Input"))
                )

        colour_vals <- analyser_aes(summary_data$analyser,
                                    c(ANALYSER_COLORS, "Input" = "black"),
                                    "black")
        shape_vals <- analyser_aes(summary_data$analyser, ANALYSER_SHAPES, 16)

        numeric_data <- summary_data

        p <- ggplot() +
                geom_line(
                        data = numeric_data,
                        aes(x = hour_num, y = mean_val,
                            colour = analyser, group = analyser),
                        linewidth = 0.45, alpha = 0.75, na.rm = TRUE
                ) +
                geom_point(
                        data = numeric_data,
                        aes(x = hour_num, y = mean_val,
                            colour = analyser, shape = analyser, group = analyser),
                        size = 1.0, na.rm = TRUE
                ) +
                geom_errorbar(
                        data = numeric_data,
                        aes(x = hour_num, ymin = mean_val - se_val, ymax = mean_val + se_val,
                            colour = analyser, group = analyser),
                        width = 0.18, linewidth = 0.30, na.rm = TRUE
                ) +
                facet_grid(facet_label ~ dataset_label,
                           scales = "free_y",
                           switch = "y",
                           labeller = labeller(facet_label = label_parsed,
                                               dataset_label = label_parsed)) +
                scale_colour_manual(values = colour_vals) +
                scale_shape_manual(values = shape_vals) +
                scale_x_continuous(
                        breaks = hour_order,
                        labels = sprintf("%02d:00", hour_order),
                        expand = expansion(mult = c(0.01, 0.01))
                ) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 5),
                        labels = scales::label_number(accuracy = 0.1, big.mark = "")
                ) +
                labs(x = NULL, y = NULL, colour = NULL, shape = NULL) +
                theme_classic(base_size = 12) +
                theme(
                        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                        panel.grid = element_blank(),
                        strip.background = element_rect(fill = "grey96", colour = "black"),
                        strip.text.y.left = element_text(angle = 0, size = 11, vjust = 0.5),
                        strip.text.x = element_text(size = 11),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                        axis.text.y = element_text(size = 8),
                        legend.position = "bottom",
                        legend.text = element_text(size = 10),
                        legend.key.width = unit(1.1, "lines"),
                        plot.margin = margin(8, 8, 8, 8)
                ) +
                guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

        plot_height <- max(6.2, 1.25 * length(var_order) + 1.8)
        ggsave(output_file, p, width = 10.2, height = plot_height, units = "in", dpi = 300, bg = "white")
        p
}

hourly_concentration_summary_plot <- function(emission_csv, output_file,
                                              location_keep = c("Indoor", "Outdoor_NE", "Outdoor_SW"),
                                              include_ftir4_old = TRUE) {
        hour_order <- 0:23
        var_order <- c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3")

        conc_long <- readr::read_csv(emission_csv, show_col_types = FALSE)
        if (!include_ftir4_old) {
                conc_long <- conc_long %>% filter(analyser != "FTIR.4_old")
        }

        conc_long <- conc_long %>%
                mutate(hour_num = as.integer(hour)) %>%
                pivot_longer(
                        cols = c(CO2_mgm3_in, CO2_mgm3_N, CO2_mgm3_S,
                                 CH4_mgm3_in, CH4_mgm3_N, CH4_mgm3_S,
                                 NH3_mgm3_in, NH3_mgm3_N, NH3_mgm3_S),
                        names_to = c("var", "suffix"),
                        names_pattern = "^(.*)_([A-Za-z]+)$",
                        values_to = "value"
                ) %>%
                mutate(
                        location = recode(suffix, "in" = "Indoor", "N" = "Outdoor_NE", "S" = "Outdoor_SW"),
                        value = as.numeric(value)
                ) %>%
                filter(location %in% location_keep,
                       var %in% var_order,
                       is.finite(value),
                       hour_num %in% hour_order) %>%
                group_by(location, var, analyser, hour_num) %>%
                summarise(
                        mean_val = mean(value, na.rm = TRUE),
                        se_val = sd(value, na.rm = TRUE) / sqrt(sum(is.finite(value))),
                        .groups = "drop"
                ) %>%
                mutate(
                        facet_label = factor(VAR_LABELS_UNITS[as.character(var)],
                                             levels = VAR_LABELS_UNITS[var_order]),
                        location_label = factor(LOC_LABELS[location],
                                                levels = LOC_LABELS[location_keep]),
                        analyser = factor(analyser,
                                          levels = ANALYSER_LEGEND_ORDER)
                )

        colour_vals <- analyser_aes(conc_long$analyser, ANALYSER_COLORS, "black")
        shape_vals <- analyser_aes(conc_long$analyser, ANALYSER_SHAPES, 16)

        p <- ggplot(conc_long,
                    aes(x = hour_num, y = mean_val,
                        colour = analyser, shape = analyser, group = analyser)) +
                geom_line(linewidth = 0.45, alpha = 0.75, na.rm = TRUE) +
                geom_point(size = 1.0, na.rm = TRUE) +
                geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                              width = 0.18, linewidth = 0.30, na.rm = TRUE) +
                (if (length(unique(conc_long$location_label)) <= 1)
                        facet_grid(facet_label ~ location_label,
                                   scales = "free_y", switch = "y",
                                   labeller = labeller(facet_label = label_parsed,
                                                       location_label = label_parsed))
                 else
                        facet_grid(facet_label ~ location_label,
                                   scales = "free_y",
                                   switch = "y",
                                   labeller = labeller(facet_label = label_parsed,
                                                       location_label = label_parsed))) +
                scale_colour_manual(values = colour_vals) +
                scale_shape_manual(values = shape_vals) +
                scale_x_continuous(
                        breaks = hour_order,
                        labels = sprintf("%02d:00", hour_order),
                        expand = expansion(mult = c(0.01, 0.01))
                ) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 5),
                        labels = scales::label_number(accuracy = 0.1, big.mark = "")
                ) +
                labs(x = NULL, y = NULL, colour = NULL, shape = NULL) +
                theme_classic(base_size = 12) +
                theme(
                        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                        panel.grid = element_blank(),
                        strip.background = element_rect(fill = "grey96", colour = "black"),
                        strip.text.y.left = element_text(angle = 0, size = 11, vjust = 0.5),
                        strip.text.x = element_text(size = 11),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                        axis.text.y = element_text(size = 8),
                        legend.position = "bottom",
                        legend.text = element_text(size = 10),
                        legend.key.width = unit(1.1, "lines"),
                        plot.margin = margin(8, 8, 8, 8)
                ) +
                guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

        plot_width <- if (length(unique(conc_long$location_label)) <= 1) 8.8 else 10.2
        ggsave(output_file, p, width = plot_width, height = 6.2, units = "in", dpi = 300, bg = "white")
        p
}

hourly_weather_input_plot <- function(emission_csv, output_file) {
        hour_order <- 0:23
        var_order <- c("n_dairycows", "temp", "wd_mst", "ws_mst")

        weather_long <- readr::read_csv(emission_csv, show_col_types = FALSE) %>%
                mutate(hour_num = as.integer(hour)) %>%
                distinct(DATE.TIME, hour_num, n_dairycows_in, temp_N, wd_mst, ws_mst) %>%
                rename(n_dairycows = n_dairycows_in, temp = temp_N) %>%
                pivot_longer(
                        cols = c(n_dairycows, temp, wd_mst, ws_mst),
                        names_to = "var",
                        values_to = "value"
                ) %>%
                mutate(value = as.numeric(value)) %>%
                filter(var %in% var_order, is.finite(value), hour_num %in% hour_order) %>%
                group_by(var, hour_num) %>%
                summarise(
                        mean_val = mean(value, na.rm = TRUE),
                        se_val = sd(value, na.rm = TRUE) / sqrt(sum(is.finite(value))),
                        .groups = "drop"
                ) %>%
                mutate(
                        facet_label = factor(VAR_LABELS_UNITS[as.character(var)],
                                             levels = VAR_LABELS_UNITS[var_order]),
                        analyser = factor("Input", levels = "Input")
                )

        p <- ggplot(weather_long,
                    aes(x = hour_num, y = mean_val, group = 1)) +
                geom_line(linewidth = 0.45, alpha = 0.75, colour = "black", na.rm = TRUE) +
                geom_point(size = 1.0, colour = "black", na.rm = TRUE) +
                geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                              width = 0.18, linewidth = 0.30, colour = "black", na.rm = TRUE) +
                facet_grid(facet_label ~ ., scales = "free_y", switch = "y",
                           labeller = labeller(facet_label = label_parsed)) +
                scale_x_continuous(
                        breaks = hour_order,
                        labels = sprintf("%02d:00", hour_order),
                        expand = expansion(mult = c(0.01, 0.01))
                ) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 5),
                        labels = scales::label_number(accuracy = 0.1, big.mark = "")
                ) +
                labs(x = NULL, y = NULL) +
                theme_classic(base_size = 12) +
                theme(
                        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                        panel.grid = element_blank(),
                        strip.background = element_rect(fill = "grey96", colour = "black"),
                        strip.text.y.left = element_text(angle = 0, size = 11, vjust = 0.5),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                        axis.text.y = element_text(size = 8),
                        plot.margin = margin(8, 8, 8, 8),
                        legend.position = "none"
                )

        ggsave(output_file, p, width = 7.2, height = 7.4, units = "in", dpi = 300, bg = "white")
        p
}

bland_altman_plot <- function(data, var_filter, analyser_pair,
                              location_filter = NULL, x = "DATE.TIME") {
        var_label_expr <- parse(text = VAR_LABELS_UNITS[[var_filter]])[[1]]
        df <- data %>% filter(var == var_filter, analyser %in% analyser_pair)
        if (!is.null(location_filter)) df <- df %>% filter(location %in% location_filter)
        df_wide <- df %>%
                select(all_of(c(x, "location", "analyser", "value"))) %>%
                pivot_wider(names_from = analyser, values_from = value)
        a1 <- analyser_pair[1]; a2 <- analyser_pair[2]
        df_ba <- df_wide %>%
                mutate(mean_val = (.data[[a1]] + .data[[a2]]) / 2,
                       diff_pct = 100 * (.data[[a2]] - .data[[a1]]) / mean_val) %>%
                filter(is.finite(diff_pct))
        bias   <- mean(df_ba$diff_pct, na.rm = TRUE)
        sd_d   <- sd(df_ba$diff_pct, na.rm = TRUE)
        loa_hi <- bias + 1.96 * sd_d; loa_lo <- bias - 1.96 * sd_d
        # Proportional-bias trend (R2.8): regress diff vs mean
        slope_fit <- lm(diff_pct ~ mean_val, data = df_ba)
        slope <- coef(slope_fit)[2]; intercept_fit <- coef(slope_fit)[1]
        subtitle_expr <- if (!is.null(location_filter)) {
                loc_expr <- parse(text = LOC_LABELS[location_filter])[[1]]
                bquote(.(var_label_expr) ~ "|" ~ .(loc_expr))
        } else var_label_expr
        ggplot(df_ba, aes(x = mean_val, y = diff_pct)) +
                geom_point(alpha = 0.5, size = 1) +
                geom_hline(yintercept = 0, color = "grey60") +
                geom_hline(yintercept = bias,   color = "blue", linetype = "dashed", linewidth = 0.9) +
                geom_hline(yintercept = loa_hi, color = "red",  linetype = "dotted", linewidth = 0.9) +
                geom_hline(yintercept = loa_lo, color = "red",  linetype = "dotted", linewidth = 0.9) +
                geom_abline(slope = slope, intercept = intercept_fit,
                            color = "darkgreen", linetype = "longdash", linewidth = 0.7) +
                annotate("text", x = Inf, y = bias, hjust = 1.05, vjust = -0.4, size = 3,
                         label = sprintf("bias = %.1f%% | slope = %.2g %%/unit",
                                         bias, slope)) +
                scale_x_continuous(breaks = pretty_breaks(n = 6),
                                   labels = scales::label_number(accuracy = 0.1)) +
                labs(subtitle = subtitle_expr,
                     x = bquote("Mean of" ~ .(a1) ~ "and" ~ .(a2)),
                     y = bquote("Relative difference (" ~ .(a2) - .(a1) ~ ") %")) +
                theme_classic() +
                theme(plot.subtitle = element_text(hjust = 0.5),
                      plot.title    = element_text(hjust = 0.5))
}

#### 0.  Paths, output directories                                            ####
base_dir   <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
data_dir   <- file.path(base_dir, "clean_data/Version_9/long_format")
clean_dir  <- file.path(base_dir, "clean_data/Version_9")
meta_dir   <- file.path(base_dir, "meta_data")
tables_dir <- file.path(base_dir, "result_data/tables/Version_11")
plots_dir  <- file.path(base_dir, "result_data/plots/Version_11")
report_dir <- file.path(base_dir, "result_data/text_reports/Version_11")
for (d in c(tables_dir, plots_dir, report_dir))
        dir.create(d, showWarnings = FALSE, recursive = TRUE)

start_time <- as.POSIXct("2025-04-08 12:00:00", tz = "UTC")
end_time   <- as.POSIXct("2025-04-14 12:00:00", tz = "UTC")
gases <- c("CO2", "CH4", "NH3")

parse_mixed_datetime_utc <- function(x) {
        parse_date_time(x,
                        orders = c("d/m/Y H:M",
                                   "Y/m/d H:M:S",
                                   "Y-m-d H:M:S",
                                   "d/m/Y H:M:S"),
                        tz = "UTC",
                        exact = FALSE)
}

read_gas_long_csv <- function(path) {
        df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
        bad_dt_names <- grep("DATE\\.TIME$", names(df), value = TRUE)
        if (!("DATE.TIME" %in% names(df)) && length(bad_dt_names) > 0) {
                names(df)[names(df) == bad_dt_names[1]] <- "DATE.TIME"
        }
        names(df) <- names(df) %>%
                str_replace_all("Analyzer", "Analyser") %>%
                str_replace_all("analyzer", "analyser")
        df
}

standardise_analyser_columns <- function(df) {
        names(df) <- names(df) %>%
                str_replace_all("Analyzer", "Analyser") %>%
                str_replace_all("analyzer", "analyser")
        df
}

# Helper for writing one analysis block of the text report.
write_section <- function(file, title, body) {
        cat(strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n",
            body, "\n\n", file = file, append = TRUE, sep = "")
}
report_file <- function(name) file.path(report_dir, name)
init_report <- function(name, header) {
        f <- report_file(name); cat(header, "\n\n", file = f, sep = "")
        invisible(f)
}

parse_lab_analyser_from_file <- function(path, prefix_regex) {
        bn <- basename(path)
        m  <- regmatches(bn, regexec(prefix_regex, bn))[[1]]
        if (length(m) < 3) return(NULL)
        tibble(lab = m[2], analyser = m[3], file = path)
}

resolution_minutes <- function(x) {
        x <- sort(unique(as.POSIXct(x, tz = "UTC")))
        if (length(x) < 2) return(NA_real_)
        round(as.numeric(median(diff(x), na.rm = TRUE), units = "mins"), 1)
}

#### 1.  Sampling-cycle figure ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â clock-style, 60 min = two 30-min cycles      ####
# Two 30-min cycles laid out around a 60-min clock face. Each step is 7.5 min:
# 3.0 min flush (grey, labelled "flushing") + 4.5 min average (white, labelled
# with the sampling location, NE/SW rendered as superscripts via plotmath).
# Ticks placed at every step boundary (0, 3, 7.5, 10.5, 15, 18, 22.5, 25.5,
# 30, 33, 37.5, 40.5, 45, 48, 52.5, 55.5) ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â i.e. each flush end + each
# average end. Box ring is intentionally thin (radial width 0.3).

cycle_segments <- tibble(
        cycle    = rep(c(1, 2), each = 8),
        step     = rep(1:4, each = 2, times = 2),
        segment  = rep(c("flush", "average"), times = 8),
        location = rep(c("Indoor", "Outdoor_NE", "Indoor", "Outdoor_SW"),
                       each = 2, times = 2)
) %>%
        mutate(step_start = (cycle - 1) * 30 + (step - 1) * 7.5,
               start_min  = step_start + ifelse(segment == "flush", 0, 3),
               end_min    = step_start + ifelse(segment == "flush", 3, 7.5),
               midpoint   = (start_min + end_min) / 2,
               # plotmath expressions ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â parse=TRUE in geom_text renders them.
               # atop() stacks "measuring" above the location label inside
               # each averaged box; NE/SW become superscripts via ^"...".
               text_label = case_when(
                       segment  == "flush"      ~ "flushing",
                       location == "Indoor"     ~ "atop('measuring','Indoor')",
                       location == "Outdoor_NE" ~ "atop('measuring',Outdoor^'NE')",
                       location == "Outdoor_SW" ~ "atop('measuring',Outdoor^'SW')"
               ))

# Major tick positions = step boundaries (these get numeric "X minutes" labels).
tick_breaks <- c(0, 3, 7.5, 10.5, 15, 18, 22.5, 25.5,
                 30, 33, 37.5, 40.5, 45, 48, 52.5, 55.5)

# Ruler/protractor minor ticks every 0.5 min ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â short radial lines just outside
# the box ring, no numeric label.
major_ticks <- tibble(x = tick_breaks)
minor_ticks <- tibble(x = setdiff(seq(0, 59.5, by = 0.5), tick_breaks))

# Inner step ring ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â one box per 7.5-min step (4 steps per cycle x 2 cycles).
step_blocks <- tibble(
        cycle = rep(c(1, 2), each = 4),
        step  = rep(1:4, times = 2)
) %>%
        mutate(start_min = (cycle - 1) * 30 + (step - 1) * 7.5,
               end_min   = start_min + 7.5,
               midpoint  = (start_min + end_min) / 2,
               label     = paste0("Step ", step))

centre_label <- tibble(x = 0, y = 0, label = "1 h\n= 2 cycles")

cycle_plot <- ggplot(cycle_segments) +
        # Inner ring ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â one box per 7.5-min step labelled "Step 1..4".
        geom_rect(data = step_blocks,
                  aes(xmin = start_min, xmax = end_min,
                      ymin = 0.80, ymax = 0.97),
                  fill = "grey97", color = "grey55", linewidth = 0.25) +
        geom_text(data = step_blocks,
                  aes(x = midpoint, y = 0.885, label = label),
                  size = 2.8) +
        # Main ring ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â flush (grey) and average (white) sub-segments.
        geom_rect(data = filter(cycle_segments, segment == "flush"),
                  aes(xmin = start_min, xmax = end_min,
                      ymin = 1.0, ymax = 1.3),
                  fill = "grey85", color = "black", linewidth = 0.3) +
        geom_rect(data = filter(cycle_segments, segment == "average"),
                  aes(xmin = start_min, xmax = end_min,
                      ymin = 1.0, ymax = 1.3),
                  fill = "white", color = "black", linewidth = 0.3) +
        geom_text(aes(x = midpoint, y = 1.15, label = text_label),
                  parse = TRUE, size = 2.8, lineheight = 0.85) +
        # Ruler-style minor ticks (lines only ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â no labels).
        geom_segment(data = minor_ticks,
                     aes(x = x, xend = x, y = 1.30, yend = 1.34),
                     color = "grey55", linewidth = 0.25) +
        # Ruler-style major ticks (longer black lines, labels via scale).
        geom_segment(data = major_ticks,
                     aes(x = x, xend = x, y = 1.30, yend = 1.42),
                     color = "black", linewidth = 0.5) +
        geom_text(data = centre_label,
                  aes(x = x, y = y, label = label),
                  size = 4.0, fontface = "bold", lineheight = 0.9) +
        coord_polar(theta = "x", start = 0, direction = 1, clip = "off") +
        scale_x_continuous(limits = c(0, 60),
                           breaks = tick_breaks,
                           labels = function(x) paste0(x, " minutes")) +
        scale_y_continuous(limits = c(0, 1.6), expand = c(0, 0)) +
        labs(title    = "60-min sampling protocol: two 30-min cycles",
             subtitle = "Each step: 3.0 min flushing (discarded) + 4.5 min measuring (retained)") +
        theme_minimal(base_size = 16) +
        theme(axis.title         = element_blank(),
              axis.text.y        = element_blank(),
              axis.ticks.y       = element_blank(),
              axis.text.x        = element_text(size = 13, face = "bold"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.title         = element_text(hjust = 0.5, size = 17),
              plot.subtitle      = element_text(hjust = 0.5, size = 12),
              plot.margin        = margin(0.4, 2.5, 0.4, 2.5, unit = "cm"))
ggsave(file.path(plots_dir, "sampling_cycle.png"), cycle_plot,
       width = 10, height = 8.5, dpi = 300, bg = "white")

# ----- Hourly-averaging methodology paragraph (for manuscript ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§2.2) ----------
hourly_averaging_paragraph <- paste(
        "Each analyser produces eight averaged 4.5-min samples per hour, one per step",
        "of the four-step cycle repeated twice. The hourly mean at each sampling",
        "location is computed as the arithmetic mean of the corresponding samples:",
        "four samples per hour for the indoor ring line (Step 1 and Step 3 of both",
        "cycles) and two samples per hour for each outdoor line (Step 2 of both cycles",
        "for Outdoor_NE; Step 4 of both cycles for Outdoor_SW). Step-to-location",
        "mapping is implemented in three ways depending on each analyser's data",
        "delivery format. Line-based FTIR analysers (ATB FTIR.1: Line 1 = NE,",
        "2 = Indoor, 3 = SW; LUFA FTIR.2: Messstelle 1 = Indoor, 2 = SW, 3 = NE) use",
        "the explicit line identifier embedded in the raw results file. MPV-based CRDS",
        "analysers (ATB CRDS.1: positions 1/2/3 = NE/Indoor/SW; UB CRDS.3: positions",
        "8/1/9 = NE/Indoor/SW; LUFA CRDS.2: positions 3/1/2 = NE/Indoor/SW) treat each",
        "contiguous run of identical MPV position as one step. Time-cycle analysers",
        "(MBBM FTIR.3 and ANECO FTIR.4) have no embedded line identifier; the step is",
        "inferred from the clock time within the cycle, anchored at the campaign",
        "start (2025-04-08 12:00 UTC) with a fixed (Indoor, Outdoor_NE, Indoor,",
        "Outdoor_SW) rotation every 7.5 min. After step assignment, the first 3 min",
        "of each step is discarded (flush) and the remaining 4.5 min is averaged.",
        "The retained averages are then aggregated to hourly means by location for",
        "downstream ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc, Q and emission calculations.",
        sep = " ")
writeLines(hourly_averaging_paragraph,
           file.path(report_dir, "07_hourly_averaging_method.txt"))

#### 2.  Read & merge gas datasets                                            ####
gas_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
gas_data <- map_dfr(gas_files, read_gas_long_csv) %>%
        select(-any_of("MPVPosition")) %>%
        mutate(DATE.TIME = parse_mixed_datetime_utc(DATE.TIME)) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        pivot_wider(id_cols     = c(DATE.TIME, lab, analyser),
                    names_from  = location,
                    values_from = c(CO2, CH4, NH3, H2O),
                    values_fn   = mean) %>%
        rename(CO2_ppm_in = CO2_in, CO2_ppm_N = CO2_N, CO2_ppm_S = CO2_S,
               CH4_ppm_in = CH4_in, CH4_ppm_N = CH4_N, CH4_ppm_S = CH4_S,
               NH3_ppm_in = NH3_in, NH3_ppm_N = NH3_N, NH3_ppm_S = NH3_S) %>%
        arrange(DATE.TIME, lab, analyser)

#### 3.  Animal / climate / wind data                                         ####
animal_data <- read.csv(file.path(meta_dir, "RGB_Animal_count/20250408-15_LVAT_Animal_data.csv"),
                        stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = dmy_hm(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        select(-hour) %>%
        rename(n_dairycows_in = n_animals) %>% distinct()

T_RH_HOBO <- read.csv(file.path(meta_dir, "HOBO_Temp_RH/2025/20250408-20250630_HOBO_Temp_RH_1hour.csv"),
                      stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = mdy_hms(paste(date, time), tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        rename(temp_N  = T_outside, RH_N = RH_out,
               temp_in = T_inside,  RH_in = RH_inside) %>%
        select(DATE.TIME, temp_N, RH_N, temp_in, RH_in)

wind_data <- read.csv(file.path(meta_dir, "USA_mast_wind/20240101_20250825_USA_mast_16_hourly_uvw_wd_ws.csv"),
                      stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = as.POSIXct(datetime_hour, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        rename(wd_mst = wd, ws_mst = ws) %>%
        select(DATE.TIME, wd_mst, ws_mst)

#### 4.  Combined inputs + emissions                                          ####
input_combined <- gas_data %>%
        left_join(animal_data, by = "DATE.TIME") %>%
        left_join(T_RH_HOBO,   by = "DATE.TIME") %>%
        left_join(wind_data,   by = "DATE.TIME") %>%
        arrange(DATE.TIME, lab, analyser)

emission_result <- indirect.CO2.balance(input_combined)
emission_reshaped <- reshaper(emission_result) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

write_excel_csv(input_combined,    file.path(tables_dir, "20250408-15_input_combined.csv"))
write_excel_csv(emission_result,   file.path(tables_dir, "20250408-15_emission_result.csv"))
write_excel_csv(emission_reshaped, file.path(tables_dir, "20250408-15_ringversuche_emission_reshaped.csv"))

#### 4.2 Campaign overview + pre-hourly data-loss summaries                   ####
gas_cycle_files <- list.files(clean_dir, pattern = "^20250408-15_7\\.5_avg_.+\\.csv$", full.names = TRUE)
gas_cycle_campaign <- map_dfr(gas_cycle_files, function(f) {
        meta <- parse_lab_analyser_from_file(f, "^20250408-15_7\\.5_avg_([^_]+)_(.+)\\.csv$")
        if (is.null(meta)) return(NULL)
        d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) %>%
                mutate(DATE.TIME = parse_mixed_datetime_utc(DATE.TIME)) %>%
                filter(DATE.TIME >= start_time & DATE.TIME <= end_time)
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O","N2O"), names(d))
        if (length(gas_cols) == 0 || !"location" %in% names(d)) return(NULL)
        d %>%
                pivot_longer(all_of(gas_cols), names_to = "variable", values_to = "value") %>%
                mutate(lab = meta$lab, analyser = meta$analyser,
                       location = recode(location,
                                         "in" = "Indoor", "N" = "Outdoor_NE", "S" = "Outdoor_SW",
                                         .default = as.character(location)))
})

stage1_files <- list.files(clean_dir, pattern = "^stage1_dropouts_.+\\.csv$", full.names = TRUE)
stage1_loss <- if (length(stage1_files) > 0) {
        map_dfr(stage1_files, ~ read.csv(.x, stringsAsFactors = FALSE) %>% standardise_analyser_columns()) %>%
                rename(variable = Variable,
                       n_before = Before,
                       n_after  = After,
                       n_removed = Removed)
} else {
        tibble()
}

cycle_avg_files <- list.files(clean_dir, pattern = "^20250408-15_7\\.5_avg_.+\\.csv$", full.names = TRUE)
cycle_loss_fallback <- map_dfr(cycle_avg_files, function(f) {
        meta <- parse_lab_analyser_from_file(f, "^20250408-15_7\\.5_avg_([^_]+)_(.+)\\.csv$")
        if (is.null(meta)) return(NULL)
        raw_file <- file.path(data_dir, paste0("20250408-15_long_", meta$lab, "_", meta$analyser, ".csv"))
        if (!file.exists(raw_file)) return(NULL)
        raw_d <- read.csv(raw_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
                mutate(DATE.TIME = parse_mixed_datetime_utc(DATE.TIME)) %>%
                filter(DATE.TIME >= start_time & DATE.TIME <= end_time)
        avg_d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) %>%
                mutate(DATE.TIME = parse_mixed_datetime_utc(DATE.TIME)) %>%
                filter(DATE.TIME >= start_time & DATE.TIME <= end_time)
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O","N2O"), intersect(names(raw_d), names(avg_d)))
        if (length(gas_cols) == 0) return(NULL)
        tibble(
                lab = meta$lab,
                analyser = meta$analyser,
                variable = gas_cols,
                n_before = map_int(gas_cols, ~ sum(is.finite(raw_d[[.x]]))),
                n_after  = map_int(gas_cols, ~ sum(is.finite(avg_d[[.x]])))
        ) %>%
                mutate(n_removed = n_before - n_after) %>%
                filter(n_removed >= 0)
})

pre_hourly_loss <- bind_rows(
        stage1_loss %>% mutate(source = "stage1_dropouts"),
        anti_join(cycle_loss_fallback, stage1_loss %>% select(lab, analyser, variable),
                  by = c("lab", "analyser", "variable")) %>% mutate(source = "derived_from_long_vs_7.5_avg")
) %>%
        arrange(lab, analyser, variable)
write_excel_csv(pre_hourly_loss, file.path(tables_dir, "pre_hourly_outlier_loss.csv"))

gas_campaign_overview <- gas_cycle_campaign %>%
        group_by(source = "gas_concentration", lab, analyser, location, variable) %>%
        summarise(total_observations = n(),
                  resolution_min = resolution_minutes(DATE.TIME),
                  n_missing = sum(!is.finite(value)),
                  mean = mean(value, na.rm = TRUE),
                  min = min(value, na.rm = TRUE),
                  max = max(value, na.rm = TRUE),
                  .groups = "drop") %>%
        left_join(pre_hourly_loss %>% select(lab, analyser, variable, n_removed),
                  by = c("lab", "analyser", "variable")) %>%
        mutate(n_removed = replace_na(n_removed, 0L))

meta_campaign_overview <- bind_rows(
        animal_data %>%
                transmute(source = "animal", lab = "meta", analyser = "RGB",
                          location = "Indoor", variable = "n_dairycows_in",
                          DATE.TIME, value = n_dairycows_in),
        T_RH_HOBO %>%
                pivot_longer(-DATE.TIME, names_to = "variable", values_to = "value") %>%
                mutate(source = "climate", lab = "meta", analyser = "HOBO",
                       location = case_when(
                               variable %in% c("temp_in", "RH_in") ~ "Indoor",
                               variable %in% c("temp_N", "RH_N")   ~ "Outdoor_NE",
                               TRUE ~ "meta")) %>%
                select(source, lab, analyser, location, variable, DATE.TIME, value),
        wind_data %>%
                pivot_longer(-DATE.TIME, names_to = "variable", values_to = "value") %>%
                mutate(source = "wind", lab = "meta", analyser = "USA", location = "Mast")
)

meta_campaign_overview <- meta_campaign_overview %>%
        group_by(source, lab, analyser, location, variable) %>%
        summarise(total_observations = n(),
                  resolution_min = resolution_minutes(DATE.TIME),
                  n_missing = sum(!is.finite(value)),
                  mean = mean(value, na.rm = TRUE),
                  min = min(value, na.rm = TRUE),
                  max = max(value, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(n_removed = 0L)

campaign_overview <- bind_rows(gas_campaign_overview, meta_campaign_overview) %>%
        mutate(across(c(mean, min, max), ~ round(.x, 3))) %>%
        arrange(source, lab, analyser, location, variable)
write_excel_csv(campaign_overview, file.path(tables_dir, "campaign_data_overview.csv"))

#### 4.5 Raw 7.5-min descriptive stats aggregation  (handover ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§2.3 / Annex A1) ####
# Single consolidated CSV for the manuscript Annex (per-Annex-table-1
# question ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â option 2). Per-analyser raw_stats_<lab>_<analyser>.csv files
# are produced by the cleaners (FTIR/CRDS_ringversuche_cleaning.R) on the
# pre-Stage-1 7.5-min data; if any analyser is missing its file (e.g. the
# cleaners have not been re-run since the new logic was added), the
# fallback computes the stats from the existing post-Stage-1 7.5_avg CSVs
# in clean_dir and flags those rows in a `source` column.

rs_files <- list.files(clean_dir,
                       pattern = "^raw_stats_.+\\.csv$",
                       full.names = TRUE)
raw_stats_pre <- if (length(rs_files) > 0) {
        map_dfr(rs_files, function(f) {
                d <- read.csv(f, stringsAsFactors = FALSE)
                d$source <- "pre-Stage-1 (from cleaner)"
                d
        })
} else {
        tibble()
}
known_pairs <- if (nrow(raw_stats_pre) > 0) {
        unique(paste(raw_stats_pre$lab, raw_stats_pre$analyser))
} else {
        character(0)
}

cycle_files <- list.files(clean_dir,
                          pattern = "^20250408-15_7\\.5_avg_.+\\.csv$",
                          full.names = TRUE)
raw_stats_fallback <- map_dfr(cycle_files, function(f) {
        bn <- basename(f)
        m  <- regmatches(bn, regexec("^20250408-15_7\\.5_avg_([^_]+)_(.+)\\.csv$", bn))[[1]]
        if (length(m) < 3) return(NULL)
        lab_x <- m[2]; ana_x <- m[3]
        if (paste(lab_x, ana_x) %in% known_pairs) return(NULL)
        d <- read.csv(f, stringsAsFactors = FALSE)
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O","N2O"), names(d))
        if (length(gas_cols) == 0 || !"location" %in% names(d)) return(NULL)
        d %>%
                pivot_longer(any_of(gas_cols), names_to = "gas", values_to = "value") %>%
                filter(!is.na(value), is.finite(value)) %>%
                group_by(location, gas) %>%
                summarise(n = n(), mean = mean(value),
                          median = median(value), min = min(value),
                          max = max(value), sd = sd(value), .groups = "drop") %>%
                mutate(lab = lab_x, analyser = ana_x,
                       source = "post-Stage-1 (from 7.5_avg CSV)",
                       .before = 1)
})

raw_stats_all <- bind_rows(raw_stats_pre, raw_stats_fallback) %>%
        arrange(lab, analyser, location, gas) %>%
        mutate(across(c(mean, median, min, max, sd), ~ round(.x, 3))) %>%
        relocate(lab, analyser, location, gas, n,
                 mean, median, min, max, sd, source)
write_excel_csv(raw_stats_all,
                file.path(tables_dir, "raw_stats_all_analysers.csv"))

#### 5.  Animal-count diagnostic  (R2.22 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â were enough cows present?)         ####
animal_summary <- animal_data %>%
        summarise(n_hours       = n(),
                  n_mean        = mean(n_dairycows_in, na.rm = TRUE),
                  n_median      = median(n_dairycows_in, na.rm = TRUE),
                  n_min         = min(n_dairycows_in, na.rm = TRUE),
                  n_max         = max(n_dairycows_in, na.rm = TRUE),
                  n_hours_lt_30 = sum(n_dairycows_in < 30, na.rm = TRUE),
                  n_hours_lt_40 = sum(n_dairycows_in < 40, na.rm = TRUE))
write_excel_csv(animal_summary, file.path(tables_dir, "animal_count_summary.csv"))


#### 6.  Concentration plots WITH FTIR.4_old (justifies its drop)             ####
c_trend_plot <- emitrendplot(emission_reshaped, y = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"))
c_boxplot    <- emiboxplot  (emission_reshaped, y = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"))
ggsave(file.path(plots_dir, "c_trend_plot.png"), c_trend_plot, width = 14, height = 6.8, dpi = 300)
ggsave(file.path(plots_dir, "c_boxplot.png"),    c_boxplot,    width = 14, height = 6.8, dpi = 300)

#### 7.  FTIR.4 vs FTIR.4_old ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â spectral-library re-evaluation check (R2.7)   ####
ftir_old <- input_combined %>% filter(analyser == "FTIR.4_old")
ftir_new <- input_combined %>% filter(analyser == "FTIR.4")
lib_pair <- function(g, l) {
        col <- paste0(g, "_ppm_", l)
        inner_join(tibble(DATE.TIME = ftir_old$DATE.TIME, v_old = ftir_old[[col]]),
                   tibble(DATE.TIME = ftir_new$DATE.TIME, v_new = ftir_new[[col]]),
                   by = "DATE.TIME") %>%
                filter(complete.cases(v_old, v_new))
}
lib_compare <- map_dfr(gases, function(g) {
        map_dfr(c("in", "N", "S"), function(l) {
                j <- lib_pair(g, l)
                if (nrow(j) < 3) return(NULL)
                dem <- deming_fit(j$v_old, j$v_new)
                tibble(gas = g, location = unname(loc_from_suffix[l]),
                       n   = nrow(j),
                       mean_old = mean(j$v_old), mean_new = mean(j$v_new),
                       mean_diff = mean(j$v_new - j$v_old),
                       RPD_pct = 100 * mean(j$v_new - j$v_old) /
                                 mean((j$v_old + j$v_new) / 2),
                       RMSD = sqrt(mean((j$v_new - j$v_old)^2)),
                       deming_intercept = unname(dem["intercept"]),
                       deming_slope     = unname(dem["slope"]),
                       CCC = DescTools::CCC(j$v_old, j$v_new)$rho.c[, "est"],
                       p_wilcox = suppressWarnings(
                               wilcox.test(j$v_old, j$v_new, paired = TRUE)$p.value))
        })
})
write_excel_csv(lib_compare, file.path(tables_dir, "FTIR4_old_vs_new_comparison.csv"))



#### 9.  Drop FTIR.4_old; rebuild working datasets  (R2.7 step 1)              ####
emission_result_v  <- emission_result %>% filter(analyser != "FTIR.4_old")

#### 9.5  Stage 2 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â mask Q/e/ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc cells when any ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc is negative (handover ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§2.6c) ####
# Steady-state mass balance: c_indoor >= c_outdoor for emitting gases.
# A negative ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc indicates either a measurement artefact (analyser noise,
# adsorption lag) or a transient mass-balance violation (e.g. brief
# outdoor plume from the western dairies blowing across the indoor ring
# line). Per handover ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§2.6c the mask is per-(analyser, hour, outdoor):
# when any of CO2/CH4/NH3 has ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc < 0 for an outdoor line in a given row,
# only that outdoor's ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc / Q / e cells are masked. Other analysers' rows
# at the same hour stay intact; the same row's _N data are unaffected
# when only _S has a negative ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc and vice versa.

stage2_dropouts <- emission_result_v %>%
        select(DATE.TIME, lab, analyser,
               matches("^delta_(CO2|CH4|NH3)_(N|S)$")) %>%
        pivot_longer(matches("^delta_"),
                     names_to      = c("gas", "outdoor"),
                     names_pattern = "^delta_(CO2|CH4|NH3)_(N|S)$",
                     values_to     = "delta") %>%
        mutate(outdoor = recode(outdoor,
                                "N" = "Outdoor_NE", "S" = "Outdoor_SW")) %>%
        group_by(analyser, gas, outdoor) %>%
        summarise(n_total      = n(),
                  n_negative   = sum(delta < 0, na.rm = TRUE),
                  pct_negative = round(100 * n_negative / n_total, 2),
                  .groups      = "drop")
write_excel_csv(stage2_dropouts, file.path(tables_dir, "stage2_dropouts.csv"))

mask_cols_N <- c("delta_CO2_N",   "delta_CH4_N", "delta_NH3_N",
                 "Q_vent_N",
                 "e_CH4_gh_N",    "e_NH3_gh_N",
                 "e_CH4_ghLU_N",  "e_NH3_ghLU_N")
mask_cols_S <- gsub("_N$", "_S", mask_cols_N)

emission_result_v <- emission_result_v %>%
        mutate(.neg_N = (delta_CO2_N < 0) | (delta_CH4_N < 0) | (delta_NH3_N < 0),
               .neg_S = (delta_CO2_S < 0) | (delta_CH4_S < 0) | (delta_NH3_S < 0)) %>%
        mutate(across(any_of(mask_cols_N),
                      ~ ifelse(.neg_N %in% TRUE, NA_real_, .x)),
               across(any_of(mask_cols_S),
                      ~ ifelse(.neg_S %in% TRUE, NA_real_, .x))) %>%
        select(-.neg_N, -.neg_S)

emission_reshaped_v <- reshaper(emission_result_v) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
        mutate(analyser = fct_drop(analyser))
input_combined_v   <- input_combined %>% filter(analyser != "FTIR.4_old")
analysers          <- sort(unique(input_combined_v$analyser))
conc_vars          <- c("CO2_mgm3","CH4_mgm3","NH3_mgm3")
delta_vars         <- c("delta_CO2","delta_CH4","delta_NH3")
qe_vars            <- c("Q_vent","e_CH4_ghLU","e_NH3_ghLU")

#### 9.7 Delta-quality diagnostics                                           ####
delta_quality <- emission_result_v %>%
        select(DATE.TIME, lab, analyser,
               delta_CO2_N, delta_CO2_S,
               delta_CH4_N, delta_CH4_S,
               delta_NH3_N, delta_NH3_S) %>%
        pivot_longer(-c(DATE.TIME, lab, analyser),
                     names_to = c("variable", "suffix"),
                     names_pattern = "^(delta_[A-Z0-9]+)_([NS])$",
                     values_to = "value") %>%
        mutate(location = recode(suffix, "N" = "Outdoor_NE", "S" = "Outdoor_SW")) %>%
        group_by(analyser, variable, location) %>%
        group_modify(~{
                x <- .x$value[is.finite(.x$value)]
                if (length(x) == 0) {
                        return(tibble(n_total = 0L, n_negative = 0L, n_near_zero = 0L,
                                      n_extreme_high = 0L, median_positive = NA_real_,
                                      near_zero_cutoff = NA_real_, high_cutoff = NA_real_))
                }
                x_pos <- x[x > 0]
                med_pos <- if (length(x_pos)) median(x_pos, na.rm = TRUE) else NA_real_
                near_cut <- if (is.finite(med_pos)) 0.1 * med_pos else NA_real_
                q3 <- suppressWarnings(quantile(x, 0.75, na.rm = TRUE, names = FALSE))
                iqr_x <- IQR(x, na.rm = TRUE)
                high_cut <- if (is.finite(q3) && is.finite(iqr_x)) q3 + 1.5 * iqr_x else NA_real_
                tibble(
                        n_total = length(x),
                        n_negative = sum(x < 0, na.rm = TRUE),
                        n_near_zero = sum(x > 0 & is.finite(near_cut) & x < near_cut, na.rm = TRUE),
                        n_extreme_high = sum(is.finite(high_cut) & x > high_cut, na.rm = TRUE),
                        median_positive = med_pos,
                        near_zero_cutoff = near_cut,
                        high_cutoff = high_cut
                )
        }) %>%
        ungroup() %>%
        mutate(across(c(median_positive, near_zero_cutoff, high_cutoff), ~ round(.x, 3)))
write_excel_csv(delta_quality, file.path(tables_dir, "delta_quality_summary.csv"))

#### 10. NE vs SW outdoor-line comparison  (R2.24 / R2.8 / E.1)                ####
# Outdoor-line choice should rely on the most credible background measurements.
# FTIR.4 is excluded from this decision due to unresolved plausibility issues,
# and FTIR.2 is excluded for NH3 only due to its known concentration offset.
OUTDOOR_SELECTION_EXCLUDE_DEFAULT <- c("FTIR.4")
OUTDOOR_SELECTION_EXCLUDE_BY_GAS <- list(
        "NH3" = c("FTIR.2", "FTIR.4")
)

get_outdoor_selection_analysers <- function(gas_name) {
        unique(c(OUTDOOR_SELECTION_EXCLUDE_DEFAULT,
                 OUTDOOR_SELECTION_EXCLUDE_BY_GAS[[gas_name]])) %>%
                purrr::discard(is.null) -> excl
        setdiff(as.character(analysers), excl)
}

# Step A: paired test per gas & analyser (campaign-wide).
ne_sw_tests <- map_dfr(gases, function(g) {
        sel_analysers <- get_outdoor_selection_analysers(g)
        map_dfr(sel_analysers, function(a) {
                d  <- input_combined_v %>% filter(analyser == a)
                ne <- d[[paste0(g, "_ppm_N")]]; sw <- d[[paste0(g, "_ppm_S")]]
                ok <- complete.cases(ne, sw); ne <- ne[ok]; sw <- sw[ok]
                if (length(ne) < 3) return(NULL)
                tibble(gas = g, analyser = a, n = length(ne),
                       mean_NE = mean(ne), mean_SW = mean(sw),
                       mean_diff = mean(ne - sw),
                       RPD_pct = 100 * mean(ne - sw) / mean((ne + sw) / 2),
                       p_ttest  = t.test(ne, sw, paired = TRUE)$p.value,
                       p_wilcox = suppressWarnings(
                               wilcox.test(ne, sw, paired = TRUE)$p.value))
        })
}) %>%
        mutate(p_wilcox_holm = p.adjust(p_wilcox, method = "holm"),
               significant   = ifelse(p_wilcox_holm < 0.05, "yes", "no"))
write_excel_csv(ne_sw_tests, file.path(tables_dir, "outdoor_NE_vs_SW_tests.csv"))

# Step B: wind-sector-conditioned NE vs SW. For each gas and 8-sector wind
# bin, average NE and SW concentrations across the retained analyser subset.
ws_long_ne_sw <- input_combined_v %>%
        filter(!is.na(wd_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst)) %>%
        select(DATE.TIME, analyser, wind_sector,
               matches("^(CO2|CH4|NH3)_ppm_(N|S)$")) %>%
        pivot_longer(cols = matches("_ppm_"),
                     names_to = c("gas", "loc"),
                     names_pattern = "(.+)_ppm_(.+)",
                     values_to = "ppm") %>%
        mutate(location = unname(loc_from_suffix[loc]),
               gas      = as.character(gas)) %>%
        rowwise() %>%
        filter(!(analyser %in% unique(c(OUTDOOR_SELECTION_EXCLUDE_DEFAULT,
                                        OUTDOOR_SELECTION_EXCLUDE_BY_GAS[[gas]])))) %>%
        ungroup() %>%
        mutate(gas = factor(gas, levels = gases))

ne_sw_by_sector <- ws_long_ne_sw %>%
        group_by(gas, location, wind_sector) %>%
        summarise(mean_ppm = mean(ppm, na.rm = TRUE),
                  sd_ppm   = sd(ppm,   na.rm = TRUE),
                  n        = n(), .groups = "drop")
write_excel_csv(ne_sw_by_sector, file.path(tables_dir, "outdoor_NE_vs_SW_by_wind_sector.csv"))

# Annotate each sector with the contamination class of each line.
sector_class <- tibble(
        wind_sector = factor(c("N","NE","E","SE","S","SW","W","NW"),
                             levels = c("N","NE","E","SE","S","SW","W","NW")),
        NE_class = case_when(
                wind_sector %in% NE_CONTAM_NEIGHBOUR ~ "neighbour_source",
                wind_sector %in% NE_CONTAM_LVAT      ~ "downwind_LVAT",
                wind_sector %in% NE_CLEAN_SECTORS    ~ "clean",
                TRUE                                 ~ "other"),
        SW_class = case_when(
                wind_sector %in% SW_CONTAM_NEIGHBOUR ~ "neighbour_source",
                wind_sector %in% SW_CONTAM_LVAT      ~ "downwind_LVAT",
                wind_sector %in% SW_CLEAN_SECTORS    ~ "clean",
                TRUE                                 ~ "other"))
write_excel_csv(sector_class, file.path(tables_dir, "wind_sector_contamination_class.csv"))

# Step C: per-line contamination summary ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â mean concentration in each class.
line_class_summary <- ne_sw_by_sector %>%
        left_join(sector_class, by = "wind_sector") %>%
        mutate(class_for_line = ifelse(location == "Outdoor_NE", NE_class, SW_class)) %>%
        group_by(gas, location, class_for_line) %>%
        summarise(mean_ppm = mean(mean_ppm, na.rm = TRUE),
                  n_sectors = n_distinct(wind_sector),
                  .groups = "drop")
write_excel_csv(line_class_summary, file.path(tables_dir, "outdoor_line_class_summary.csv"))

# Step E: bar plot of RPD vs wind sector to make the contamination pattern
# easy to read for the manuscript.
ne_sw_rpd_sector <- ne_sw_by_sector %>%
        select(gas, location, wind_sector, mean_ppm) %>%
        pivot_wider(names_from = location, values_from = mean_ppm) %>%
        mutate(RPD_pct = 100 * (Outdoor_NE - Outdoor_SW) /
                        ((Outdoor_NE + Outdoor_SW) / 2)) %>%
        left_join(sector_class, by = "wind_sector")
rpd_sector_plot <- ggplot(ne_sw_rpd_sector, aes(x = wind_sector, y = RPD_pct, fill = gas)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7,
                 color = "black", linewidth = 0.2) +
        geom_hline(yintercept = 0, color = "grey50") +
        facet_wrap(~ gas, ncol = 1, scales = "free_y") +
        scale_fill_brewer(palette = "Dark2", guide = "none") +
        labs(title = "Relative percentage difference NE - SW by wind sector",
             subtitle = "Positive = NE higher than SW; negative = SW higher than NE",
             x = "Wind sector (mast)", y = "RPD (%)") +
        theme_bw(base_size = 12)
ggsave(file.path(plots_dir, "outdoor_NE_vs_SW_RPD_by_sector.png"), rpd_sector_plot,
       width = 10, height = 9, dpi = 300, bg = "white")

# Outdoor concentration plot by wind sector and analyser (mg m-3), excluding
# the indoor line and the superseded FTIR.4_old diagnostic series.
sector_conc_by_analyser <- emission_result_v %>%
        filter(!is.na(wd_mst),
               analyser %in% c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4")) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst)) %>%
        select(DATE.TIME, analyser, wind_sector,
               matches("^(CO2|CH4|NH3)_mgm3_(N|S)$")) %>%
        pivot_longer(cols = matches("_mgm3_"),
                     names_to = c("gas", "loc"),
                     names_pattern = "(.+)_mgm3_(.+)",
                     values_to = "mgm3") %>%
        mutate(
                gas = factor(gas, levels = c("CO2", "CH4", "NH3")),
                location = factor(unname(loc_from_suffix[loc]),
                                  levels = c("Outdoor_NE", "Outdoor_SW")),
                analyser = factor(as.character(analyser),
                                  levels = c("CRDS.1","CRDS.2","CRDS.3",
                                             "FTIR.1","FTIR.2","FTIR.3","FTIR.4")),
                wind_sector = factor(wind_sector,
                                     levels = c("N","NE","E","SE","S","SW","W","NW"))
        ) %>%
        group_by(gas, analyser, location, wind_sector) %>%
        summarise(mean_mgm3 = mean(mgm3, na.rm = TRUE),
                  sd_mgm3   = sd(mgm3, na.rm = TRUE),
                  n         = n(),
                  .groups = "drop")
write_excel_csv(sector_conc_by_analyser,
                file.path(tables_dir, "outdoor_sector_concentration_by_analyser.csv"))

sector_conc_plot <- ggplot(
        sector_conc_by_analyser,
        aes(x = wind_sector, y = mean_mgm3, fill = location)
) +
        geom_col(
                position = position_dodge(width = 0.8),
                width = 0.72,
                colour = "black",
                linewidth = 0.22
        ) +
        facet_grid(
                rows = vars(gas),
                cols = vars(analyser),
                scales = "free_y",
                switch = "y",
                labeller = labeller(
                        gas = as_labeller(c(
                                "CO2" = "c[CO2]~'(mg '*m^-3*')'",
                                "CH4" = "c[CH4]~'(mg '*m^-3*')'",
                                "NH3" = "c[NH3]~'(mg '*m^-3*')'"
                        ), label_parsed)
                )
        ) +
        scale_fill_manual(
                values = c("Outdoor_NE" = "grey75", "Outdoor_SW" = "grey40"),
                labels = c("Outdoor_NE" = "Outdoor^NE", "Outdoor_SW" = "Outdoor^SW")
        ) +
        labs(
                x = "Wind sector",
                y = NULL,
                fill = "Outdoor line"
        ) +
        theme_bw(base_size = 15) +
        theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey95", colour = "black"),
                strip.text = element_text(size = 14),
                strip.text.y.left = element_text(angle = 0, size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
                axis.text.y = element_text(size = 11),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_blank(),
                legend.position = "bottom",
                legend.title = element_text(size = 13),
                legend.text = element_text(size = 12)
        )
ggsave(file.path(plots_dir, "outdoor_sector_concentration_by_analyser.png"),
       sector_conc_plot, width = 18, height = 8.5, dpi = 300, bg = "white")

#### 11. Decide which outdoor line to retain                                  ####
# Decision rule:
#   * The retained line is the one whose mean over its "clean" sectors is
#     closer to its mean over the "downwind_LVAT" sectors AND that has the
#     lower neighbour-source contamination, weighted by the time fraction
#     each line spends in each class during this campaign.
#   * In a tie, prefer the line on the upwind side of the prevailing wind.
# Practically this is read off the line_class_summary + wind-rose
# frequencies below; the actual drop is written here so all downstream
# tables stay consistent.

wind_freq <- input_combined_v %>%
        distinct(DATE.TIME, wd_mst) %>%
        filter(!is.na(wd_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst)) %>%
        count(wind_sector, .drop = FALSE) %>%
        mutate(freq = n / sum(n))
write_excel_csv(wind_freq, file.path(tables_dir, "wind_sector_frequency.csv"))

# Time spent with each line in each class.
line_time_class <- wind_freq %>%
        left_join(sector_class, by = "wind_sector") %>%
        pivot_longer(cols = c(NE_class, SW_class),
                     names_to = "line", values_to = "class") %>%
        mutate(line = recode(line, "NE_class" = "Outdoor_NE", "SW_class" = "Outdoor_SW")) %>%
        group_by(line, class) %>% summarise(time_frac = sum(freq), .groups = "drop")
write_excel_csv(line_time_class, file.path(tables_dir, "outdoor_line_time_in_class.csv"))

contam_score <- line_time_class %>%
        filter(class %in% c("downwind_LVAT", "neighbour_source")) %>%
        group_by(line) %>% summarise(contam_time = sum(time_frac), .groups = "drop") %>%
        arrange(contam_time)
RETAINED_OUTDOOR <- contam_score$line[1]
DROPPED_OUTDOOR  <- contam_score$line[contam_score$line != RETAINED_OUTDOOR][1]

# Persist the decision for the manuscript.
cat(sprintf("retained_outdoor: %s\ndropped_outdoor:  %s\n",
            RETAINED_OUTDOOR, DROPPED_OUTDOOR),
    file = file.path(tables_dir, "outdoor_line_choice.txt"))

#### 11.5 Per-analyser IQR filter on Q_vent + dropout diagnostic              ####
# Replaces the V10 [MIN_Q_VENT, MAX_Q_VENT] sanity bound that was previously
# applied inside indirect.CO2.balance(). Reviewer-requested change: caps were
# silent and conflated genuine low-gradient hours (already handled by
# MIN_DELTA_CO2 in ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§2.3) with per-analyser anomalies (especially FTIR.4).
# The new filter applies 1.5*IQR per analyser to the surviving Q_vent values
# and records the number of rows lost at every stage so the loss can be
# reported in Results ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§3.2 and Discussion ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§4.4.

# Count inputs and MIN_DELTA_CO2 dropouts per (analyser, outdoor line).
# An input row is a row of emission_result_v with a finite delta_CO2_*
# (i.e. the row had both indoor and outdoor CO2 measurements). It is
# dropped by MIN_DELTA_CO2 when delta_CO2_* < MIN_DELTA_CO2.
qvent_dropouts_min_delta <- emission_result_v %>%
        group_by(analyser) %>%
        summarise(
                n_input_N        = sum(is.finite(delta_CO2_N)),
                n_dropped_minD_N = sum(is.finite(delta_CO2_N) & delta_CO2_N <  MIN_DELTA_CO2),
                n_after_minD_N   = sum(is.finite(delta_CO2_N) & delta_CO2_N >= MIN_DELTA_CO2),
                n_input_S        = sum(is.finite(delta_CO2_S)),
                n_dropped_minD_S = sum(is.finite(delta_CO2_S) & delta_CO2_S <  MIN_DELTA_CO2),
                n_after_minD_S   = sum(is.finite(delta_CO2_S) & delta_CO2_S >= MIN_DELTA_CO2),
                .groups = "drop")

# Apply per-analyser 1.5*IQR removal on Q_vent_N and Q_vent_S.
iqr_drop <- function(x) {
        if (sum(is.finite(x)) < 4) return(x)
        q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
        H <- 1.5 * IQR(x, na.rm = TRUE)
        ifelse(x < (q[1] - H) | x > (q[2] + H), NA_real_, x)
}
emission_result_v <- emission_result_v %>%
        group_by(analyser) %>%
        mutate(
                Q_vent_N_pre_iqr = Q_vent_N,
                Q_vent_S_pre_iqr = Q_vent_S,
                Q_vent_N         = iqr_drop(Q_vent_N),
                Q_vent_S         = iqr_drop(Q_vent_S)
        ) %>%
        ungroup() %>%
        # Re-propagate the IQR-filtered Q into the derived emission columns.
        mutate(
                e_NH3_gh_N   = (delta_NH3_N * Q_vent_N / 1000) * n_dairycows_in,
                e_CH4_gh_N   = (delta_CH4_N * Q_vent_N / 1000) * n_dairycows_in,
                e_NH3_gh_S   = (delta_NH3_S * Q_vent_S / 1000) * n_dairycows_in,
                e_CH4_gh_S   = (delta_CH4_S * Q_vent_S / 1000) * n_dairycows_in,
                e_NH3_ghLU_N = (e_NH3_gh_N * 500) / (n_dairycows_in * m_weight),
                e_CH4_ghLU_N = (e_CH4_gh_N * 500) / (n_dairycows_in * m_weight),
                e_NH3_ghLU_S = (e_NH3_gh_S * 500) / (n_dairycows_in * m_weight),
                e_CH4_ghLU_S = (e_CH4_gh_S * 500) / (n_dairycows_in * m_weight)
        )

qvent_dropouts_iqr <- emission_result_v %>%
        group_by(analyser) %>%
        summarise(
                n_after_minD_N   = sum(is.finite(Q_vent_N_pre_iqr)),
                n_dropped_iqr_N  = sum(is.finite(Q_vent_N_pre_iqr) & !is.finite(Q_vent_N)),
                n_kept_N         = sum(is.finite(Q_vent_N)),
                n_after_minD_S   = sum(is.finite(Q_vent_S_pre_iqr)),
                n_dropped_iqr_S  = sum(is.finite(Q_vent_S_pre_iqr) & !is.finite(Q_vent_S)),
                n_kept_S         = sum(is.finite(Q_vent_S)),
                .groups = "drop")

qvent_dropouts <- qvent_dropouts_min_delta %>%
        select(analyser, n_input_N, n_dropped_minD_N,
               n_input_S, n_dropped_minD_S) %>%
        left_join(qvent_dropouts_iqr, by = "analyser") %>%
        mutate(
                pct_kept_N = 100 * n_kept_N / pmax(n_input_N, 1),
                pct_kept_S = 100 * n_kept_S / pmax(n_input_S, 1)
        ) %>%
        select(analyser,
               n_input_N, n_dropped_minD_N, n_dropped_iqr_N, n_kept_N, pct_kept_N,
               n_input_S, n_dropped_minD_S, n_dropped_iqr_S, n_kept_S, pct_kept_S)

write_excel_csv(qvent_dropouts,
                file.path(tables_dir, "qvent_filter_dropouts.csv"))

#### 12. Build single-outdoor working dataset                                 ####
# Map RETAINED -> suffix column for *_N / *_S in emission_result_v
retain_suffix  <- c("Outdoor_NE" = "N", "Outdoor_SW" = "S")[[RETAINED_OUTDOOR]]
drop_suffix    <- c("Outdoor_NE" = "N", "Outdoor_SW" = "S")[[DROPPED_OUTDOOR]]
# Keep only the retained-outdoor delta/Q/e columns in emission_result_single
emission_result_single <- emission_result_v %>%
        select(-matches(paste0("_", drop_suffix, "$"))) %>%
        rename_with(~ str_replace(.x, paste0("_", retain_suffix, "$"), ""),
                    .cols = matches(paste0("_", retain_suffix, "$")))

# Reshape ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â single outdoor, so location is Indoor + retained outdoor only.
# We feed reshaper a frame with three "location-coded" CO2/CH4/NH3 columns
# but with delta/Q/e already collapsed; pivot manually.
single_long <- emission_result_v %>%
        select(DATE.TIME, lab, analyser,
               CO2_mgm3_in, CH4_mgm3_in, NH3_mgm3_in,
               any_of(paste0(c("CO2_mgm3","CH4_mgm3","NH3_mgm3"), "_", retain_suffix)),
               any_of(paste0(c("delta_CO2","delta_CH4","delta_NH3",
                                "Q_vent","e_CH4_gh","e_NH3_gh",
                                "e_CH4_ghLU","e_NH3_ghLU"), "_", retain_suffix))) %>%
        rename_with(~ str_replace(.x, paste0("_", retain_suffix, "$"), "_OUT"),
                    .cols = matches(paste0("_", retain_suffix, "$"))) %>%
        pivot_longer(cols = -c(DATE.TIME, lab, analyser),
                     names_to = "raw_var", values_to = "value") %>%
        mutate(location = case_when(
                str_detect(raw_var, "_in$")  ~ "Indoor",
                str_detect(raw_var, "_OUT$") ~ RETAINED_OUTDOOR,
                TRUE                         ~ "Combined"),
               var = str_remove(raw_var, "_(in|OUT)$")) %>%
        select(DATE.TIME, lab, analyser, location, var, value)

single_long <- single_long %>%
        mutate(analyser = factor(analyser,
                                 levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4",
                                            "CRDS.1","CRDS.2","CRDS.3")))
write_excel_csv(single_long, file.path(tables_dir,
                                       sprintf("emission_long_single_outdoor_%s.csv",
                                               RETAINED_OUTDOOR)))

fit_repeated_measures_model <- function(df_long, vars, location_filter) {
        if (!requireNamespace("nlme", quietly = TRUE)) {
                stop("Package 'nlme' is required for repeated-measures mixed models.")
        }

        model_input <- df_long %>%
                filter(location == location_filter, var %in% vars) %>%
                transmute(
                        DATE.TIME = as.POSIXct(DATE.TIME, tz = "UTC"),
                        time_id = factor(DATE.TIME),
                        analyser = factor(analyser,
                                          levels = c("CRDS.1","CRDS.2","CRDS.3",
                                                     "FTIR.1","FTIR.2","FTIR.3","FTIR.4")),
                        var,
                        value
                ) %>%
                filter(is.finite(value)) %>%
                mutate(analyser = forcats::fct_drop(analyser))

        anova_rows  <- list()
        fixed_rows  <- list()
        random_rows <- list()

        for (v in vars) {
                dat_v <- model_input %>% filter(var == v)
                if (nrow(dat_v) == 0) next

                fit_v <- nlme::lme(
                        fixed = value ~ analyser,
                        random = ~1 | time_id,
                        data = dat_v,
                        na.action = na.omit,
                        control = nlme::lmeControl(returnObject = TRUE)
                )

                aov_v <- as.data.frame(nlme::anova.lme(fit_v))
                aov_v$term <- rownames(aov_v)
                rownames(aov_v) <- NULL
                names(aov_v) <- gsub("^p-value$", "p_value", names(aov_v))
                anova_rows[[v]] <- aov_v %>%
                        as_tibble() %>%
                        mutate(variable = v, location = location_filter, .before = 1) %>%
                        select(variable, location, term, numDF, denDF, `F-value`, p_value)

                fixed_rows[[v]] <- summary(fit_v)$tTable %>%
                        as.data.frame() %>%
                        tibble::rownames_to_column("term") %>%
                        as_tibble() %>%
                        mutate(variable = v, location = location_filter, .before = 1) %>%
                        rename(
                                estimate = Value,
                                std_error = `Std.Error`,
                                df = DF,
                                t_value = `t-value`,
                                p_value = `p-value`
                        ) %>%
                        select(variable, location, term, estimate, std_error, df, t_value, p_value)

                vc_v <- nlme::VarCorr(fit_v)
                random_rows[[v]] <- tibble(
                        variable = v,
                        location = location_filter,
                        component = c("time_intercept_sd", "residual_sd"),
                        estimate = c(
                                as.numeric(vc_v[1, "StdDev"]),
                                as.numeric(vc_v[nrow(vc_v), "StdDev"])
                        ),
                        n_obs = nrow(dat_v),
                        n_time = n_distinct(dat_v$time_id),
                        n_analysers = n_distinct(dat_v$analyser)
                )
        }

        list(
                anova = bind_rows(anova_rows),
                fixed = bind_rows(fixed_rows),
                random = bind_rows(random_rows)
        )
}

mixed_model_sw <- fit_repeated_measures_model(
        df_long = single_long,
        vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
        location_filter = "Outdoor_SW"
)

write_excel_csv(
        mixed_model_sw$anova,
        file.path(tables_dir, "mixed_model_analyser_effect_Outdoor_SW_anova.csv")
)
write_excel_csv(
        mixed_model_sw$fixed,
        file.path(tables_dir, "mixed_model_analyser_effect_Outdoor_SW_fixed_effects.csv")
)
write_excel_csv(
        mixed_model_sw$random,
        file.path(tables_dir, "mixed_model_analyser_effect_Outdoor_SW_random_effects.csv")
)

mixed_model_sw_report <- mixed_model_sw$anova %>%
        filter(term == "analyser") %>%
        mutate(
                F_value = round(`F-value`, 3),
                p_value = scales::pvalue(p_value, accuracy = 0.001)
        ) %>%
        transmute(
                line = sprintf(
                        "%s at %s: analyser fixed effect -> F(%s, %s) = %s, p = %s",
                        recode(variable,
                               "Q_vent" = "Q",
                               "e_CH4_ghLU" = "eCH4",
                               "e_NH3_ghLU" = "eNH3"),
                        recode(location, "Outdoor_SW" = "Outdoor^SW"),
                        numDF, denDF, F_value, p_value
                )
        ) %>%
        pull(line)
write_lines(
        c(
                "Repeated-measures mixed model: value ~ analyser + (1 | time_id)",
                "Dataset: retained single-outdoor workflow (Outdoor^SW)",
                mixed_model_sw_report
        ),
        file.path(report_dir, "mixed_model_analyser_effect_Outdoor_SW.txt")
)

mixed_model_sw_effects <- mixed_model_sw$fixed %>%
        filter(term != "(Intercept)") %>%
        mutate(
                analyser = str_remove(term, "^analyser"),
                conf_low = estimate - 1.96 * std_error,
                conf_high = estimate + 1.96 * std_error
        ) %>%
        select(variable, location, analyser, estimate, std_error, conf_low, conf_high, p_value)

mixed_model_sw_effects <- bind_rows(
        tibble(
                variable = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                location = "Outdoor_SW",
                analyser = "CRDS.1",
                estimate = 0,
                std_error = 0,
                conf_low = 0,
                conf_high = 0,
                p_value = NA_real_
        ),
        mixed_model_sw_effects
) %>%
        mutate(
                analyser = factor(analyser,
                                  levels = c("CRDS.1","CRDS.2","CRDS.3",
                                             "FTIR.1","FTIR.2","FTIR.3","FTIR.4")),
                variable = factor(variable,
                                  levels = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"))
        )

mixed_model_sw_axis_guides <- mixed_model_sw_effects %>%
        group_by(variable) %>%
        summarise(
                max_abs = max(abs(c(conf_low, conf_high, estimate)), na.rm = TRUE),
                .groups = "drop"
        ) %>%
        mutate(
                max_abs = ifelse(max_abs == 0, 1, max_abs),
                max_abs = max_abs * 1.12,
                analyser = factor("CRDS.1",
                                  levels = c("CRDS.1","CRDS.2","CRDS.3",
                                             "FTIR.1","FTIR.2","FTIR.3","FTIR.4"))
        ) %>%
        select(variable, analyser, max_abs)

write_excel_csv(
        mixed_model_sw_effects,
        file.path(tables_dir, "mixed_model_analyser_effect_Outdoor_SW_effect_plot_data.csv")
)

mixed_model_effect_plot <- ggplot(
        mixed_model_sw_effects,
        aes(x = analyser, y = estimate, colour = analyser)
) +
        geom_blank(data = mixed_model_sw_axis_guides, aes(y = max_abs), inherit.aes = FALSE) +
        geom_blank(data = mixed_model_sw_axis_guides, aes(y = -max_abs), inherit.aes = FALSE) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
        geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.16, linewidth = 0.8) +
        geom_point(size = 2.6) +
        facet_wrap(
                ~ variable,
                ncol = 3,
                scales = "free_y",
                labeller = labeller(variable = as_labeller(VAR_LABELS_UNITS, label_parsed))
        ) +
        scale_colour_manual(values = ANALYSER_COLORS[names(ANALYSER_COLORS) != "FTIR.4_old"]) +
        scale_y_continuous(
                breaks = function(lims) {
                        lim <- max(abs(lims), na.rm = TRUE)
                        seq(-lim, lim, length.out = 9)
                },
                labels = function(x) {
                        lim <- max(abs(x), na.rm = TRUE)
                        if (lim >= 10) {
                                format(round(x, 0), trim = TRUE, scientific = FALSE, nsmall = 0)
                        } else if (lim >= 1) {
                                format(round(x, 2), trim = TRUE, scientific = FALSE, nsmall = 2)
                        } else {
                                format(round(x, 2), trim = TRUE, scientific = FALSE, nsmall = 2)
                        }
                }
        ) +
        labs(x = NULL, y = "Estimated analyser effect relative to CRDS.1") +
        theme_bw(base_size = 14) +
        theme(
                panel.grid = element_blank(),
                strip.background = element_rect(fill = "grey96", colour = "black"),
                strip.text = element_text(size = 14, face = "plain"),
                axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13),
                legend.position = "none"
        )

ggsave(
        file.path(plots_dir, "mixed_model_analyser_effect_Outdoor_SW.png"),
        mixed_model_effect_plot,
        width = 13.8,
        height = 4.9,
        dpi = 300,
        bg = "white"
)

#### 13. Delta / Q / e plots ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â single outdoor                                 ####
single_for_plot <- single_long %>%
        mutate(day  = factor(as.Date(DATE.TIME)),
               hour = factor(format(DATE.TIME, "%H:%M")))

make_plot_frame_for_outdoor <- function(df_long, outdoor_label, suffix_tag) {
        out <- df_long %>%
                filter(location %in% c("Indoor", outdoor_label)) %>%
                mutate(day  = factor(as.Date(DATE.TIME)),
                       hour = factor(format(DATE.TIME, "%H:%M")))
        d_trend <- emitrendplot(out,
                                y = c("delta_CO2","delta_CH4","delta_NH3"),
                                location_filter = outdoor_label)
        d_box <- emiboxplot(out,
                            y = c("delta_CO2","delta_CH4","delta_NH3"),
                            location_filter = outdoor_label)
        qe_trend <- emitrendplot(out,
                                 y = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                                 location_filter = outdoor_label)
        qe_box <- emiboxplot(out,
                             y = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                             location_filter = outdoor_label)
        ggsave(file.path(plots_dir, paste0("d_trend_plot_", suffix_tag, ".png")), d_trend,
               width = 10.5, height = 6.2, dpi = 300)
        ggsave(file.path(plots_dir, paste0("d_boxplot_", suffix_tag, ".png")), d_box,
               width = 10.5, height = 6.2, dpi = 300)
        ggsave(file.path(plots_dir, paste0("q_e_trend_plot_", suffix_tag, ".png")), qe_trend,
               width = 10.5, height = 6.2, dpi = 300)
        ggsave(file.path(plots_dir, paste0("q_e_boxplot_", suffix_tag, ".png")), qe_box,
               width = 10.5, height = 6.2, dpi = 300)
}

make_plot_frame_for_outdoor(emission_reshaped_v, "Outdoor_NE", "Outdoor_NE")
make_plot_frame_for_outdoor(single_for_plot, RETAINED_OUTDOOR, "Outdoor_SW")

hourly_dataset_summary_plot(
        emission_csv = file.path(tables_dir, "20250408-15_emission_result.csv"),
        output_file  = file.path(plots_dir, "hourly_dataset_summary_delta.png"),
        var_order    = c("delta_CO2", "delta_CH4", "delta_NH3"),
        include_ftir4_old = FALSE
)

hourly_dataset_summary_plot(
        emission_csv = file.path(tables_dir, "20250408-15_emission_result.csv"),
        output_file  = file.path(plots_dir, "hourly_dataset_summary_q_e.png"),
        var_order    = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
        include_ftir4_old = FALSE
)

hourly_concentration_summary_plot(
        emission_csv = file.path(tables_dir, "20250408-15_emission_result.csv"),
        output_file  = file.path(plots_dir, "hourly_concentration_summary_indoor.png"),
        location_keep = "Indoor",
        include_ftir4_old = TRUE
)

hourly_concentration_summary_plot(
        emission_csv = file.path(tables_dir, "20250408-15_emission_result.csv"),
        output_file  = file.path(plots_dir, "hourly_concentration_summary_outdoors.png"),
        location_keep = c("Outdoor_NE", "Outdoor_SW"),
        include_ftir4_old = TRUE
)

hourly_weather_input_plot(
        emission_csv = file.path(tables_dir, "20250408-15_emission_result.csv"),
        output_file  = file.path(plots_dir, "hourly_weather_input_summary.png")
)

#### 14. Bland-Altman, shared-line analyser pairs                             ####
# Relative Bland-Altman analysis is applied to analyser pairs sharing the same
# sampling line within a laboratory (Lab A and Lab B).
ba_pairs <- list(AnalyserA = c("FTIR.1", "CRDS.1"),
                 AnalyserB = c("FTIR.2", "CRDS.2"))

save_ba_family <- function(data, vars, locations, analyser_pair, tag, file, ncol = 2,
                           width = 10, height = 8) {
        panels <- map(locations, function(loc) {
                map(vars, function(v) {
                        bland_altman_plot(data, var_filter = v,
                                          analyser_pair = analyser_pair,
                                          location_filter = loc) +
                                theme(plot.margin = margin(8, 8, 8, 8))
                })
        }) %>% unlist(recursive = FALSE)
        ggsave(file.path(plots_dir, file),
               wrap_plots(panels, ncol = ncol),
               width = width, height = height, units = "in", dpi = 300)
}

save_ba_by_location <- function(data, vars, locations, analyser_pair, tag, prefix,
                                width = 12, height = 7) {
        walk(locations, function(loc) {
                file <- paste0(prefix, "_", loc, "_", tag, ".png")
                save_ba_family(data,
                               vars = vars,
                               locations = loc,
                               analyser_pair = analyser_pair,
                               tag = tag,
                               file = file,
                               ncol = length(vars),
                               width = width,
                               height = height)
        })
}

save_ba_by_location(emission_reshaped_v,
                    vars = conc_vars,
                    locations = c("Indoor", "Outdoor_NE", "Outdoor_SW"),
                    analyser_pair = ba_pairs$AnalyserA,
                    tag = "AnalyserA",
                    prefix = "c_BlandAltman",
                    width = 12, height = 4.5)
save_ba_by_location(emission_reshaped_v,
                    vars = conc_vars,
                    locations = c("Indoor", "Outdoor_NE", "Outdoor_SW"),
                    analyser_pair = ba_pairs$AnalyserB,
                    tag = "AnalyserB",
                    prefix = "c_BlandAltman",
                    width = 12, height = 4.5)

save_ba_by_location(emission_reshaped_v,
                    vars = delta_vars,
                    locations = c("Outdoor_NE", "Outdoor_SW"),
                    analyser_pair = ba_pairs$AnalyserA,
                    tag = "AnalyserA",
                    prefix = "d_BlandAltman",
                    width = 12, height = 4.5)
save_ba_by_location(emission_reshaped_v,
                    vars = delta_vars,
                    locations = c("Outdoor_NE", "Outdoor_SW"),
                    analyser_pair = ba_pairs$AnalyserB,
                    tag = "AnalyserB",
                    prefix = "d_BlandAltman",
                    width = 12, height = 4.5)

save_ba_by_location(emission_reshaped_v,
                    vars = qe_vars,
                    locations = c("Outdoor_NE", "Outdoor_SW"),
                    analyser_pair = ba_pairs$AnalyserA,
                    tag = "AnalyserA",
                    prefix = "qe_BlandAltman",
                    width = 12, height = 4.5)
save_ba_by_location(emission_reshaped_v,
                    vars = qe_vars,
                    locations = c("Outdoor_NE", "Outdoor_SW"),
                    analyser_pair = ba_pairs$AnalyserB,
                    tag = "AnalyserB",
                    prefix = "qe_BlandAltman",
                    width = 12, height = 4.5)

ba_specs <- tribble(
        ~group_name,      ~data_name,            ~vars,        ~locations,
        "concentration",  "emission_reshaped_v", conc_vars,    c("Indoor", "Outdoor_NE", "Outdoor_SW"),
        "delta",          "emission_reshaped_v", delta_vars,   c("Outdoor_NE", "Outdoor_SW"),
        "derived",        "emission_reshaped_v", qe_vars,      c("Outdoor_NE", "Outdoor_SW")
)

ba_data_lookup <- list(emission_reshaped_v = emission_reshaped_v)
ba_table <- map_dfr(names(ba_pairs), function(tag) {
        pr <- ba_pairs[[tag]]
        map_dfr(seq_len(nrow(ba_specs)), function(i) {
                vars_i <- ba_specs$vars[[i]]
                locs_i <- ba_specs$locations[[i]]
                dat_i  <- ba_data_lookup[[ba_specs$data_name[[i]]]]
                map_dfr(locs_i, function(loc) {
                        map_dfr(vars_i, function(v) {
                                w <- dat_i %>%
                                        filter(var == v, location == loc, analyser %in% pr) %>%
                                        select(DATE.TIME, analyser, value) %>%
                                        pivot_wider(names_from = analyser, values_from = value)
                                if (!all(pr %in% names(w))) return(NULL)
                                st <- ba_relative_stats(w[[pr[1]]], w[[pr[2]]])
                                if (is.null(st)) return(NULL)
                                bind_cols(tibble(pair = tag,
                                                 group_name = ba_specs$group_name[[i]],
                                                 variable = v,
                                                 location = loc),
                                          st)
                        })
                })
        })
})
write_excel_csv(ba_table, file.path(tables_dir, "bland_altman_relative_all.csv"))

#### 15. Pairwise method comparison ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â extended to concentration c + ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc       ####
# R2.9 / R1.13: report Deming slope+intercept + Lin's CCC + Pearson r for
# every analyser pair, not just for Q/e. Also exposes ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âc agreement.
pairwise_compare <- function(df_long, var_sel, loc_sel) {
        w <- df_long %>%
                filter(var == var_sel, location == loc_sel,
                       !analyser %in% c("HOBO","USA","RGB")) %>%
                select(DATE.TIME, analyser, value) %>%
                mutate(analyser = as.character(analyser)) %>%
                pivot_wider(names_from = analyser, values_from = value,
                            values_fn = ~ mean(.x, na.rm = TRUE))
        ana <- setdiff(names(w), "DATE.TIME")
        if (length(ana) < 2) return(NULL)
        map_dfr(combn(ana, 2, simplify = FALSE), function(pr) {
                x <- w[[pr[1]]]; y <- w[[pr[2]]]
                ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
                if (length(x) < 5) return(NULL)
                dem <- deming_fit(x, y)
                lm_fit <- lm(y ~ x)
                tibble(var = var_sel, location = loc_sel,
                       analyser_x = pr[1], analyser_y = pr[2], n = length(x),
                       pearson_r        = cor(x, y),
                       ccc              = DescTools::CCC(x, y)$rho.c[, "est"],
                       deming_slope     = unname(dem["slope"]),
                       deming_intercept = unname(dem["intercept"]),
                       ols_slope        = unname(coef(lm_fit)[2]),
                       ols_intercept    = unname(coef(lm_fit)[1]),
                       r_squared        = summary(lm_fit)$r.squared)
        })
}
# Inputs: emission_reshaped_v carries both location-coded concentrations and deltas
conc_vars <- c("CO2_mgm3","CH4_mgm3","NH3_mgm3")
delta_vars    <- c("delta_CO2","delta_CH4","delta_NH3")
qe_vars       <- c("Q_vent","e_CH4_ghLU","e_NH3_ghLU")
locs_abs      <- c("Indoor","Outdoor_NE","Outdoor_SW")
locs_delta_qe <- c("Outdoor_NE","Outdoor_SW")

pairwise_abs_source <- emission_reshaped %>%
        filter(analyser != "FTIR.4") %>%
        mutate(analyser = ifelse(analyser == "FTIR.4_old", "FTIR.4_old", as.character(analyser)))

pairwise_abs <- map_dfr(conc_vars, function(v)
        map_dfr(locs_abs, function(l) pairwise_compare(pairwise_abs_source, v, l)))
pairwise_delta <- map_dfr(delta_vars, function(v)
        map_dfr(locs_delta_qe, function(l) pairwise_compare(emission_reshaped_v, v, l)))
pairwise_qe <- map_dfr(qe_vars, function(v)
        map_dfr(locs_delta_qe, function(l) pairwise_compare(emission_reshaped_v, v, l)))
pairwise_tbl <- bind_rows(pairwise_abs, pairwise_delta, pairwise_qe)
write_excel_csv(pairwise_tbl, file.path(tables_dir, "pairwise_regression_ccc.csv"))

# CCC heatmaps for each variable group
heatmap_ccc_conc_tri <- function(tbl, var_set, plain_labels, title, file) {
        analyser_order <- c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4_old")
        d <- tbl %>%
                filter(var %in% var_set,
                       analyser_x %in% analyser_order,
                       analyser_y %in% analyser_order) %>%
                mutate(facet_label = factor(plain_labels[var], levels = plain_labels[var_set]),
                       analyser_x = factor(analyser_x, levels = analyser_order),
                       analyser_y = factor(analyser_y, levels = rev(analyser_order)))
        p <- ggplot(d, aes(x = analyser_x, y = analyser_y, fill = ccc)) +
                geom_tile(color = "white", linewidth = 0.4) +
                geom_text(aes(label = sprintf("%.2f", ccc)), size = 3.4, color = "white") +
                scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#1a9850",
                                     midpoint = 0.5, limits = c(0, 1), name = "Lin's CCC") +
                scale_y_discrete(position = "right") +
                facet_grid(location ~ facet_label,
                           switch = "y",
                           labeller = labeller(facet_label = label_parsed,
                                               location    = as_labeller(LOC_LABELS, label_parsed))) +
                coord_fixed() +
                labs(x = NULL, y = NULL) +
                theme_bw(base_size = 11) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      axis.text.y.right = element_text(angle = 45, hjust = 0.5),
                      axis.text.y.left = element_blank(),
                      axis.ticks.y.left = element_blank(),
                      strip.placement = "outside",
                      legend.position = "bottom",
                      panel.grid = element_blank())
        ggsave(file.path(plots_dir, file), p, width = 12, height = 9, dpi = 300, bg = "white")
}

heatmap_ccc <- function(tbl, var_set, plain_labels, title, file, location_filter = NULL) {
        analyser_order <- c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4")
        d <- tbl %>%
                filter(var %in% var_set) %>%
                { if (!is.null(location_filter)) filter(., location %in% location_filter) else . } %>%
                filter(analyser_x %in% analyser_order,
                       analyser_y %in% analyser_order) %>%
                mutate(facet_label = factor(plain_labels[var], levels = plain_labels[var_set]),
                       analyser_x = factor(analyser_x, levels = analyser_order),
                       analyser_y = factor(analyser_y, levels = rev(analyser_order)))
        p <- ggplot(d, aes(x = analyser_x, y = analyser_y, fill = ccc)) +
                geom_tile(color = "white") +
                geom_text(aes(label = sprintf("%.2f", ccc)), size = 2.6) +
                scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#1a9850",
                                     midpoint = 0.5, limits = c(0, 1), name = "Lin's CCC") +
                scale_y_discrete(position = "right") +
                facet_grid(location ~ facet_label,
                           switch = "y",
                           labeller = labeller(facet_label = label_parsed,
                                               location    = as_labeller(LOC_LABELS, label_parsed))) +
                labs(x = NULL, y = NULL) +
                theme_bw(base_size = 11) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      axis.text.y.right = element_text(angle = 45, hjust = 0.5),
                      axis.text.y.left = element_blank(),
                      axis.ticks.y.left = element_blank(),
                      strip.placement = "outside",
                      legend.position = "bottom")
        out_width  <- if (!is.null(location_filter)) 10 else 12
        out_height <- if (!is.null(location_filter)) 4.8 else 9
        ggsave(file.path(plots_dir, file), p, width = out_width, height = out_height, dpi = 300, bg = "white")
}
heatmap_ccc_conc_tri(pairwise_abs, conc_vars, VAR_LABELS_PLAIN,
                     "Pairwise agreement, concentrations (includes FTIR.4_old; Lin's CCC)",
                     "pairwise_ccc_concentration_c.png")
heatmap_ccc(pairwise_delta, delta_vars, VAR_LABELS_PLAIN,
            "Pairwise agreement, ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â concentrations - Outdoor^NE (Lin's CCC)",
            "pairwise_ccc_delta_c_Outdoor_NE.png",
            location_filter = "Outdoor_NE")
heatmap_ccc(pairwise_delta, delta_vars, VAR_LABELS_PLAIN,
            "Pairwise agreement, ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â concentrations - Outdoor^SW (Lin's CCC)",
            "pairwise_ccc_delta_c_Outdoor_SW.png",
            location_filter = "Outdoor_SW")
heatmap_ccc(pairwise_qe, qe_vars, VAR_LABELS_PLAIN,
            "Pairwise agreement, ventilation + emissions - Outdoor^NE (Lin's CCC)",
            "pairwise_ccc_q_e_Outdoor_NE.png",
            location_filter = "Outdoor_NE")
heatmap_ccc(pairwise_qe, qe_vars, VAR_LABELS_PLAIN,
            "Pairwise agreement, ventilation + emissions - Outdoor^SW (Lin's CCC)",
            "pairwise_ccc_q_e_Outdoor_SW.png",
            location_filter = "Outdoor_SW")

regression_scatter_panels <- function(df_long, pair_tbl, vars, location_filter,
                                      file, width = 10, height = 8, ncol = 6) {
        spec_tbl <- pair_tbl %>%
                filter(var %in% vars, location == location_filter) %>%
                distinct(var, location, analyser_x, analyser_y,
                         pearson_r, ccc, deming_slope, deming_intercept, r_squared) %>%
                mutate(panel_label = paste0("atop(",
                                            VAR_LABELS_PLAIN[var],
                                            ", '", analyser_x, " vs ", analyser_y, "')"))
        if (!nrow(spec_tbl)) return(invisible(NULL))

        plot_df <- purrr::pmap_dfr(spec_tbl[, c("var", "location", "analyser_x", "analyser_y", "panel_label")],
                                   function(var, location, analyser_x, analyser_y, panel_label) {
                w <- df_long %>%
                        filter(var == .env$var,
                               location == .env$location,
                               analyser %in% c(.env$analyser_x, .env$analyser_y)) %>%
                        select(DATE.TIME, analyser, value) %>%
                        pivot_wider(names_from = analyser, values_from = value,
                                    values_fn = ~ mean(.x, na.rm = TRUE))
                if (!all(c(analyser_x, analyser_y) %in% names(w))) return(NULL)
                out <- w %>%
                        transmute(x = .data[[analyser_x]], y = .data[[analyser_y]]) %>%
                        filter(complete.cases(x, y))
                if (nrow(out) < 5) return(NULL)
                out %>%
                        mutate(var = var, location = location,
                               analyser_x = analyser_x, analyser_y = analyser_y,
                               panel_label = panel_label)
        })

        ann <- spec_tbl %>%
                mutate(label = sprintf("y = %.2f + %.2f x\nr = %.2f\nCCC = %.2f\nR^2 = %.2f",
                                       deming_intercept, deming_slope, pearson_r, ccc, r_squared))
        if (!nrow(plot_df) || !nrow(ann)) return(invisible(NULL))

        p <- ggplot(plot_df, aes(x = x, y = y)) +
                geom_point(shape = 16, size = 0.9, alpha = 0.45, color = "#4d4d4d") +
                geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                            linewidth = 0.45, color = "grey55") +
                geom_abline(data = ann,
                            aes(slope = deming_slope, intercept = deming_intercept),
                            color = "#1a9850", linewidth = 0.7, inherit.aes = FALSE) +
                geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
                          hjust = -0.05, vjust = 1.1, size = 2.8, inherit.aes = FALSE) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        text = element_text(size = 14),
                        panel.border = element_rect(color = "black", fill = NA),
                        panel.grid = element_blank(),
                        strip.text = element_text(face = "bold", size = 11),
                        plot.title = element_blank()
                ) +
                facet_wrap(~ panel_label, scales = "free", ncol = ncol,
                           labeller = label_parsed)
        ggsave(file.path(plots_dir, file), p, width = width, height = height, dpi = 300, bg = "white")
}
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("delta_CO2"), "Outdoor_NE",
                          "pairwise_regression_delta_CO2_Outdoor_NE.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("delta_CO2"), "Outdoor_SW",
                          "pairwise_regression_delta_CO2_Outdoor_SW.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("delta_CH4"), "Outdoor_NE",
                          "pairwise_regression_delta_CH4_Outdoor_NE.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("delta_CH4"), "Outdoor_SW",
                          "pairwise_regression_delta_CH4_Outdoor_SW.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("delta_NH3"), "Outdoor_NE",
                          "pairwise_regression_delta_NH3_Outdoor_NE.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("delta_NH3"), "Outdoor_SW",
                          "pairwise_regression_delta_NH3_Outdoor_SW.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("Q_vent"), "Outdoor_NE",
                          "pairwise_regression_Q_vent_Outdoor_NE.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("Q_vent"), "Outdoor_SW",
                          "pairwise_regression_Q_vent_Outdoor_SW.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("e_CH4_ghLU"), "Outdoor_NE",
                          "pairwise_regression_e_CH4_Outdoor_NE.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("e_CH4_ghLU"), "Outdoor_SW",
                          "pairwise_regression_e_CH4_Outdoor_SW.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("e_NH3_ghLU"), "Outdoor_NE",
                          "pairwise_regression_e_NH3_Outdoor_NE.png",
                          width = 18, height = 6, ncol = 7)
regression_scatter_panels(emission_reshaped_v, pairwise_tbl, c("e_NH3_ghLU"), "Outdoor_SW",
                          "pairwise_regression_e_NH3_Outdoor_SW.png",
                          width = 18, height = 6, ncol = 7)

#### 16. Tukey HSD pairwise tables  (Tables 3, 4, 6 of the manuscript)        ####
# One-way ANOVA per (variable, location) on the analyser factor, followed by
# Tukey's honestly significant difference post-hoc. Each pair is reported with
# the per-analyser means, the relative difference (RD %), the Tukey-adjusted
# p-value, and a significance code (ns / * / ** / ***).
# RD follows the manuscript equation for analyser-vs-analyser comparison:
#   pairwise RD = 100 * (mean_1 - mean_2) / mean_2
#
# Three tables are produced:
#   Table 3 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â concentrations  (c_CO2, c_CH4, c_NH3) x 3 locations
#   Table 4 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â concentration differences (delta_CO2, delta_CH4, delta_NH3) x 2 outdoors
#   Table 6 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â derived quantities       (Q_vent, e_CH4_ghLU, e_NH3_ghLU) x 2 outdoors
#
# Inputs: emission_reshaped_v (post-FTIR.4_old drop). FTIR.2 is kept because
# the NH3 concentration offset cancels in deltas ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â see Section 8 diagnostics.

# Significance code helper (ns / * / ** / ***)
sig_code <- function(p) {
        ifelse(is.na(p), "",
        ifelse(p < 0.001, "***",
        ifelse(p < 0.01,  "**",
        ifelse(p < 0.05,  "*",  "ns"))))
}

# Run Tukey HSD for one (variable, location) panel and return long-form rows.
tukey_panel <- function(df_long, var_sel, loc_sel) {
        x <- df_long %>%
                filter(var == var_sel, location == loc_sel,
                       !analyser %in% c("HOBO","USA","RGB"),
                       !is.na(value), is.finite(value)) %>%
                mutate(analyser = as.character(analyser))
        if (nrow(x) < 10 || length(unique(x$analyser)) < 2) return(NULL)
        an  <- factor(x$analyser)
        fit <- aov(value ~ an, data = data.frame(value = x$value, an = an))
        tk  <- TukeyHSD(fit, conf.level = 0.95)$an
        means <- x %>% group_by(analyser) %>%
                summarise(m = mean(value, na.rm = TRUE), n = n(), .groups = "drop")
        pairs <- rownames(tk)
        map_dfr(seq_along(pairs), function(i) {
                ab <- strsplit(pairs[i], "-", fixed = TRUE)[[1]]
                ax <- ab[1]; ay <- ab[2]
                mx <- means$m[means$analyser == ax]
                my <- means$m[means$analyser == ay]
                if (length(mx) == 0 || length(my) == 0) return(NULL)
                rd <- ifelse(my != 0, 100 * (mx - my) / my, NA_real_)
                tibble(variable = var_sel, location = loc_sel,
                       analyser_1 = ax, analyser_2 = ay,
                       n_1     = means$n[means$analyser == ax],
                       n_2     = means$n[means$analyser == ay],
                       mean_1  = round(mx, 3), mean_2 = round(my, 3),
                       diff    = round(tk[i, "diff"], 3),
                       lwr     = round(tk[i, "lwr"], 3),
                       upr     = round(tk[i, "upr"], 3),
                       RD_pct  = round(rd, 1),
                       p_tukey = signif(tk[i, "p adj"], 3))
        })
}

# Walk the (variable x location) grid, apply Tukey, and add significance codes.
tukey_table <- function(df_long, vars, locs) {
        out <- map_dfr(vars, function(v) {
                map_dfr(locs, function(l) tukey_panel(df_long, v, l))
        })
        if (nrow(out) == 0) return(out)
        out %>%
                mutate(sig = sig_code(p_tukey)) %>%
                arrange(variable, location, analyser_1, analyser_2)
}

build_pairwise_plot_data <- function(tukey_tbl, df_long, vars, locs,
                                     analyser_order) {
        panel_keys <- tidyr::expand_grid(variable = vars, location = locs)
        out <- map2_dfr(panel_keys$variable, panel_keys$location, function(v, l) {
                tuk <- tukey_tbl %>%
                        filter(variable == v, location == l) %>%
                        transmute(variable, location,
                                  row_analyser = analyser_1,
                                  col_analyser = analyser_2,
                                  cell_type = "upper",
                                  label = sig,
                                  p_value = p_tukey,
                                  RD_pct = NA_real_) %>%
                        bind_rows(
                                tukey_tbl %>%
                                        filter(variable == v, location == l) %>%
                                        transmute(variable, location,
                                                  row_analyser = analyser_2,
                                                  col_analyser = analyser_1,
                                                  cell_type = "lower",
                                                  label = formatC(RD_pct, digits = 1, format = "f"),
                                                  p_value = NA_real_,
                                                  RD_pct = RD_pct)
                        )

                diag_cells <- tibble(
                        variable = v,
                        location = l,
                        row_analyser = analyser_order,
                        col_analyser = analyser_order,
                        cell_type = "diag",
                        label = "-",
                        p_value = NA_real_,
                        RD_pct = NA_real_
                )

                bind_rows(tuk, diag_cells) %>%
                        mutate(
                                row_analyser = factor(row_analyser, levels = analyser_order),
                                col_analyser = factor(col_analyser, levels = analyser_order)
                        )
        })
        out
}

plot_pairwise_matrix <- function(plot_df, vars, locs, analyser_order, file_name,
                                 width, height, fill_title = "Tukey HSD",
                                 fixed_aspect = TRUE) {
        if (!nrow(plot_df)) return(invisible(NULL))

        full_grid <- tidyr::expand_grid(
                variable = vars,
                location = locs,
                row_analyser = factor(analyser_order, levels = analyser_order),
                col_analyser = factor(analyser_order, levels = analyser_order)
        )

        annot_df <- plot_df %>%
                filter(cell_type == "upper", !is.na(label), nzchar(label))

        rd_df <- plot_df %>% filter(cell_type == "lower", !is.na(RD_pct))
        diag_df <- plot_df %>% filter(cell_type == "diag")

        p <- ggplot(full_grid, aes(x = col_analyser, y = row_analyser)) +
                geom_tile(fill = "white", color = "white", linewidth = 0.45) +
                geom_tile(data = annot_df, fill = "white", color = "white", linewidth = 0.45) +
                geom_tile(data = rd_df, fill = "white", color = "white", linewidth = 0.45) +
                geom_tile(data = diag_df, fill = "grey95", color = "white", linewidth = 0.45) +
                geom_text(data = annot_df, aes(label = label), size = 4.9, family = "", color = "black") +
                geom_text(data = rd_df, aes(label = label), size = 4.5, family = "", color = "black") +
                geom_text(data = diag_df, aes(label = label), size = 4.9, family = "", color = "black") +
                facet_grid(
                        rows = vars(factor(location, levels = locs, labels = LOC_LABELS[locs])),
                        cols = vars(factor(variable, levels = vars, labels = VAR_LABELS_PLAIN[vars])),
                        switch = "y",
                        labeller = labeller(.rows = label_parsed, .cols = label_parsed)
                ) +
                labs(x = NULL, y = NULL) +
                theme_bw(base_size = 15) +
                theme(
                        panel.grid = element_blank(),
                        panel.spacing = unit(0.16, "lines"),
                        strip.background = element_rect(fill = "grey96", colour = "black"),
                        strip.text = element_text(size = 14, face = "plain"),
                        strip.text.y.left = element_text(angle = 0),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                        axis.text.y.right = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 13),
                        axis.text.y.left = element_blank(),
                        axis.ticks.y.left = element_blank(),
                        axis.ticks.y.right = element_line(),
                        axis.text.y = element_blank(),
                        legend.position = "none"
                ) +
                scale_y_discrete(position = "right")

        if (fixed_aspect) {
                p <- p + coord_fixed()
        }

        ggsave(file.path(plots_dir, file_name), p,
               width = width, height = height, dpi = 300, bg = "white")
        invisible(NULL)
}

tukey_abs   <- tukey_table(emission_reshaped_v,
                           vars = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
                           locs = c("Indoor","Outdoor_NE","Outdoor_SW"))
tukey_delta <- tukey_table(emission_reshaped_v,
                           vars = c("delta_CO2","delta_CH4","delta_NH3"),
                           locs = c("Outdoor_NE","Outdoor_SW"))
tukey_qe    <- tukey_table(emission_reshaped_v,
                           vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                           locs = c("Outdoor_NE","Outdoor_SW"))

write_excel_csv(tukey_abs,   file.path(tables_dir, "tukey_concentrations_long.csv"))
write_excel_csv(tukey_delta, file.path(tables_dir, "tukey_delta_concentrations_long.csv"))
write_excel_csv(tukey_qe,    file.path(tables_dir, "tukey_ventilation_emission_long.csv"))

# Remove legacy matrix-style table exports so Version_11 keeps the PNG panels as
# the manuscript-facing output format for these comparisons.
legacy_tukey_files <- c(
        "tukey_concentrations.csv",
        "tukey_delta_concentrations.csv",
        "tukey_delta_concentrations_Outdoor_NE.csv",
        "tukey_delta_concentrations_Outdoor_SW.csv",
        "tukey_ventilation_emission.csv",
        "tukey_ventilation_emission_Outdoor_NE.csv",
        "tukey_ventilation_emission_Outdoor_SW.csv"
)
walk(file.path(tables_dir, legacy_tukey_files), ~ if (file.exists(.x)) file.remove(.x))
legacy_matrix_dir <- file.path(tables_dir, "triangular_pairwise_matrices")
if (dir.exists(legacy_matrix_dir)) unlink(legacy_matrix_dir, recursive = TRUE, force = TRUE)

legacy_matrix_plots <- c(
        "pairwise_tukey_rd_concentrations_all.png",
        "pairwise_tukey_rd_concentrations_Indoor.png",
        "pairwise_tukey_rd_concentrations_Outdoor_NE.png",
        "pairwise_tukey_rd_concentrations_Outdoor_SW.png",
        "pairwise_tukey_rd_delta_all.png",
        "pairwise_tukey_rd_delta.png",
        "pairwise_tukey_rd_qe_all.png",
        "pairwise_tukey_rd_qe.png"
)
walk(file.path(plots_dir, legacy_matrix_plots), ~ if (file.exists(.x)) file.remove(.x))

conc_plot_df <- build_pairwise_plot_data(
        tukey_tbl = tukey_abs,
        df_long = emission_reshaped_v,
        vars = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
        locs = c("Indoor","Outdoor_NE","Outdoor_SW"),
        analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4")
)

delta_ne_plot_df <- build_pairwise_plot_data(
        tukey_tbl = tukey_delta,
        df_long = emission_reshaped_v,
        vars = c("delta_CO2","delta_CH4","delta_NH3"),
        locs = c("Outdoor_NE"),
        analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4")
)

delta_sw_plot_df <- build_pairwise_plot_data(
        tukey_tbl = tukey_delta,
        df_long = emission_reshaped_v,
        vars = c("delta_CO2","delta_CH4","delta_NH3"),
        locs = c("Outdoor_SW"),
        analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4")
)

qe_ne_plot_df <- build_pairwise_plot_data(
        tukey_tbl = tukey_qe,
        df_long = emission_reshaped_v,
        vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
        locs = c("Outdoor_NE"),
        analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4")
)

qe_sw_plot_df <- build_pairwise_plot_data(
        tukey_tbl = tukey_qe,
        df_long = emission_reshaped_v,
        vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
        locs = c("Outdoor_SW"),
        analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4")
)

plot_pairwise_matrix(
        conc_plot_df,
        vars = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
        locs = c("Indoor","Outdoor_NE","Outdoor_SW"),
        analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4"),
        file_name = "pairwise_tukey_rd_concentrations.png",
        width = 22, height = 12,
        fixed_aspect = FALSE
)
plot_pairwise_matrix(delta_ne_plot_df %>% filter(location == "Outdoor_NE"),
                     vars = c("delta_CO2","delta_CH4","delta_NH3"),
                     locs = c("Outdoor_NE"),
                     analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4"),
                     file_name = "pairwise_tukey_rd_delta_Outdoor_NE.png",
                     width = 22, height = 6.2,
                     fixed_aspect = FALSE)
plot_pairwise_matrix(delta_sw_plot_df %>% filter(location == "Outdoor_SW"),
                     vars = c("delta_CO2","delta_CH4","delta_NH3"),
                     locs = c("Outdoor_SW"),
                     analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4"),
                     file_name = "pairwise_tukey_rd_delta_Outdoor_SW.png",
                     width = 22, height = 6.2,
                     fixed_aspect = FALSE)
plot_pairwise_matrix(qe_ne_plot_df %>% filter(location == "Outdoor_NE"),
                     vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                     locs = c("Outdoor_NE"),
                     analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4"),
                     file_name = "pairwise_tukey_rd_qe_Outdoor_NE.png",
                     width = 22, height = 6.2,
                     fixed_aspect = FALSE)
plot_pairwise_matrix(qe_sw_plot_df %>% filter(location == "Outdoor_SW"),
                     vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                     locs = c("Outdoor_SW"),
                     analyser_order = c("CRDS.1","CRDS.2","CRDS.3","FTIR.1","FTIR.2","FTIR.3","FTIR.4"),
                     file_name = "pairwise_tukey_rd_qe_Outdoor_SW.png",
                     width = 22, height = 6.2,
                     fixed_aspect = FALSE)
# Long-form CSV layout (one row per analyser pair) is still kept for reporting:
#   variable, location, analyser_1, analyser_2, n_1, n_2, mean_1, mean_2,
#   diff (mean_1 - mean_2), lwr/upr (95% CI on diff from Tukey),
#   RD_pct = 100 * (mean_1 - mean_2) / mean_2,
#   p_tukey (Tukey-adjusted within the (variable, location) panel),
#   sig (ns / * / ** / ***).

#### 17. Wind characterisation ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â daily campaign + seasonal historical         ####
# Two polar-bar plots. In both, bars are stacked by wind-speed bin and each
# sector is annotated (in blue) with its hourly count and per-panel % share.
#
#   wind_rose_daily.png    ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â 7 panels, one per campaign day
#                            (08.04.2025 to 14.04.2025; 15.04 deliberately
#                            excluded because the campaign ends 14.04 12:00)
#   wind_rose_seasonal.png ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â 4 panels, historical seasons in calendar order
#                            (Summer 2024 -> Autumn 2024 -> Winter 2024/25
#                            -> Spring 2025 up to 07.04.2025)

# Shared wind-speed bins + colour palette (Figure-2 reference style:
# burgundy at low speeds -> blue at high speeds).
WS_BREAKS  <- c(0, 0.5, 1.0, 2.0, 3.0, 4.0, Inf)
WS_LABELS  <- c("0.0-0.5","0.5-1.0","1.0-2.0","2.0-3.0","3.0-4.0",">=4.0")
WS_COLOURS <- c("0.0-0.5" = "#8C1838",   # dark wine
                "0.5-1.0" = "#E04344",   # red
                "1.0-2.0" = "#F49649",   # orange
                "2.0-3.0" = "#E3E394",   # pale yellow-green
                "3.0-4.0" = "#4DBC8E",   # teal-green
                ">=4.0"   = "#3E81BA")   # blue

# Read the FULL wind series (campaign window + historical year).
wind_full <- read.csv(file.path(meta_dir, "USA_mast_wind/20240101_20250825_USA_mast_16_hourly_uvw_wd_ws.csv"),
                      stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = as.POSIXct(datetime_hour, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
        rename(wd_mst = wd, ws_mst = ws) %>%
        filter(!is.na(wd_mst), !is.na(ws_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst),
               ws_bin      = cut(ws_mst, breaks = WS_BREAKS, labels = WS_LABELS,
                                 include.lowest = TRUE, right = FALSE))

# --- (a) Daily wind roses for the 7 campaign days ----------------------------
daily_data <- wind_full %>%
        filter(DATE.TIME >= start_time, DATE.TIME <= end_time) %>%
        mutate(day_label = format(as.Date(DATE.TIME), "%d.%m.%Y")) %>%
        group_by(day_label) %>% mutate(day_total = n()) %>% ungroup()

daily_stack <- daily_data %>%
        group_by(day_label, day_total, wind_sector, ws_bin) %>%
        summarise(n_bin = n(), .groups = "drop") %>%
        mutate(pct = 100 * n_bin / day_total)

daily_totals <- daily_stack %>%
        group_by(day_label, wind_sector) %>%
        summarise(total_n   = sum(n_bin),
                  total_pct = sum(pct),
                  .groups   = "drop")
# Annotate ONLY the most prevalent sector per day, with the share as "X.X%".
daily_top <- daily_totals %>%
        group_by(day_label) %>%
        slice_max(total_pct, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(label = sprintf("%.1f%%", total_pct))

# Order panels chronologically (dd.mm.yyyy string sort wouldn't be chronological).
day_levels   <- daily_data %>% distinct(day = as.Date(DATE.TIME), day_label) %>%
        arrange(day) %>% pull(day_label)
daily_stack <- daily_stack %>% mutate(day_label = factor(day_label, levels = day_levels))
daily_top   <- daily_top   %>% mutate(day_label = factor(day_label, levels = day_levels))

# Ring labels ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â placed at the S compass position so they sit on each ring
# radius going outward from the centre. One copy per panel.
RING_STEPS_DAILY <- c(20, 40, 60, 80)
ring_labels_daily <- expand.grid(day_label = factor(day_levels, levels = day_levels),
                                 y_value   = RING_STEPS_DAILY,
                                 stringsAsFactors = FALSE) %>%
        mutate(label       = sprintf("%d%%", y_value),
               wind_sector = factor("S", levels = c("N","NE","E","SE","S","SW","W","NW")))

daily_rose <- ggplot(daily_stack, aes(x = wind_sector, y = pct, fill = ws_bin)) +
        geom_col(color = "black", linewidth = 0.2, width = 0.95) +
        geom_text(data = ring_labels_daily,
                  aes(x = wind_sector, y = y_value, label = label),
                  color = "grey30", size = 3.3, inherit.aes = FALSE) +
        geom_text(data = daily_top,
                  aes(x = wind_sector,
                      y = pmin(total_pct + 4, 88),
                      label = label),
                  color = "#1F4E79", size = 5.0,
                  fontface = "bold", inherit.aes = FALSE) +
        coord_polar(start = -pi / 8) +
        facet_wrap(~ day_label, nrow = 1) +
        scale_fill_manual(values = WS_COLOURS, drop = FALSE,
                          name = expression("Wind speed (m s"^-1*")"),
                          guide = guide_legend(nrow = 1, byrow = TRUE)) +
        scale_y_continuous(limits = c(0, 80),
                           breaks = RING_STEPS_DAILY) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 14) +
        theme(legend.position  = "bottom",
              strip.text       = element_text(face = "bold", size = 14),
              strip.background = element_rect(fill = "grey95"),
              axis.text.x      = element_text(size = 13),
              axis.text.y      = element_blank(),
              axis.ticks.y     = element_blank(),
              panel.grid.minor = element_blank(),
              legend.text      = element_text(size = 13),
              legend.title     = element_text(size = 14))
ggsave(file.path(plots_dir, "wind_rose_daily.png"), daily_rose,
       width = 24, height = 4.8, dpi = 300, bg = "white")

# --- (b) Seasonal wind roses for the year preceding the campaign -------------
SEASON_RANGES <- tibble::tribble(
        ~season_label,                              ~start,       ~end,
        "Summer 2024\n01.06.2024 - 31.08.2024",     "2024-06-01", "2024-08-31",
        "Autumn 2024\n01.09.2024 - 30.11.2024",     "2024-09-01", "2024-11-30",
        "Winter 2024/25\n01.12.2024 - 28.02.2025",  "2024-12-01", "2025-02-28",
        "Spring 2025\n01.03.2025 - 07.04.2025",     "2025-03-01", "2025-04-07"
)

seasonal_data <- map_dfr(seq_len(nrow(SEASON_RANGES)), function(i) {
        rng <- SEASON_RANGES[i, ]
        wind_full %>%
                filter(DATE.TIME >= as.POSIXct(paste(rng$start, "00:00:00"), tz = "UTC"),
                       DATE.TIME <= as.POSIXct(paste(rng$end,   "23:59:59"), tz = "UTC")) %>%
                mutate(season = rng$season_label)
}) %>%
        mutate(season = factor(season, levels = SEASON_RANGES$season_label))

seasonal_stack <- seasonal_data %>%
        group_by(season) %>% mutate(season_total = n()) %>% ungroup() %>%
        group_by(season, season_total, wind_sector, ws_bin) %>%
        summarise(n_bin = n(), .groups = "drop") %>%
        mutate(pct = 100 * n_bin / season_total)

seasonal_totals <- seasonal_stack %>%
        group_by(season, wind_sector) %>%
        summarise(total_n   = sum(n_bin),
                  total_pct = sum(pct),
                  .groups   = "drop")
# Annotate ONLY the most prevalent sector per season, with the share as "X.X%".
seasonal_top <- seasonal_totals %>%
        group_by(season) %>%
        slice_max(total_pct, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(label = sprintf("%.1f%%", total_pct))

RING_STEPS_SEASONAL <- c(10, 20, 30, 40)
ring_labels_seasonal <- expand.grid(
        season  = factor(SEASON_RANGES$season_label,
                         levels = SEASON_RANGES$season_label),
        y_value = RING_STEPS_SEASONAL,
        stringsAsFactors = FALSE) %>%
        mutate(label       = sprintf("%d%%", y_value),
               wind_sector = factor("S", levels = c("N","NE","E","SE","S","SW","W","NW")))

seasonal_rose <- ggplot(seasonal_stack, aes(x = wind_sector, y = pct, fill = ws_bin)) +
        geom_col(color = "black", linewidth = 0.2, width = 0.95) +
        geom_text(data = ring_labels_seasonal,
                  aes(x = wind_sector, y = y_value, label = label),
                  color = "grey30", size = 3.0, inherit.aes = FALSE) +
        geom_text(data = seasonal_top,
                  aes(x = wind_sector,
                      y = pmin(total_pct + 3, 88),
                      label = label),
                  color = "#1F4E79", size = 4.6,
                  fontface = "bold", inherit.aes = FALSE) +
        coord_polar(start = -pi / 8) +
        facet_wrap(~ season, nrow = 1) +
        scale_fill_manual(values = WS_COLOURS, drop = FALSE,
                          name = expression("Wind speed (m s"^-1*")"),
                          guide = guide_legend(nrow = 1, byrow = TRUE)) +
        scale_y_continuous(limits = c(0, 40),
                           breaks = RING_STEPS_SEASONAL) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 14) +
        theme(legend.position  = "bottom",
              strip.text       = element_text(face = "bold", size = 14),
              strip.background = element_rect(fill = "grey95"),
              axis.text.x      = element_text(size = 13),
              axis.text.y      = element_blank(),
              axis.ticks.y     = element_blank(),
              panel.grid.minor = element_blank(),
              legend.text      = element_text(size = 13),
              legend.title     = element_text(size = 14))
ggsave(file.path(plots_dir, "wind_rose_seasonal.png"), seasonal_rose,
       width = 20, height = 6.5, dpi = 300, bg = "white")

#### 18. Headline numbers ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â campaign mean Q and e on retained outdoor         ####
# These are the numbers that go into the manuscript's ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§3.3/3.4 headlines.
headline_stats <- single_for_plot %>%
        filter(var %in% c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
               !analyser %in% c("HOBO","USA","RGB"),
               is.finite(value)) %>%
        group_by(var, analyser) %>%
        summarise(n = n(), mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE), .groups = "drop")
headline_consortium <- headline_stats %>%
        group_by(var) %>%
        summarise(n_analysers       = n(),
                  n_obs             = sum(n),
                  consortium_mean   = mean(mean),
                  consortium_median = median(mean),
                  sd_across         = sd(mean),
                  se_mean           = sd_across / sqrt(n_analysers),
                  ci_lo95           = consortium_mean - 1.96 * se_mean,
                  ci_hi95           = consortium_mean + 1.96 * se_mean,
                  ci_half_pct       = 1.96 * se_mean / consortium_mean * 100,
                  .groups = "drop")
# Also a version excluding FTIR.4 from the consortium ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â FTIR.4 has 2-3x higher
# SD in derived quantities than the others even after library correction, so a
# robust headline must be reported in both forms.
headline_consortium_noF4 <- headline_stats %>%
        filter(analyser != "FTIR.4") %>%
        group_by(var) %>%
        summarise(n_analysers       = n(),
                  n_obs             = sum(n),
                  consortium_mean   = mean(mean),
                  consortium_median = median(mean),
                  sd_across         = sd(mean),
                  se_mean           = sd_across / sqrt(n_analysers),
                  ci_lo95           = consortium_mean - 1.96 * se_mean,
                  ci_hi95           = consortium_mean + 1.96 * se_mean,
                  ci_half_pct       = 1.96 * se_mean / consortium_mean * 100,
                  .groups = "drop")
write_excel_csv(headline_stats,           file.path(tables_dir, "headline_per_analyser.csv"))
write_excel_csv(headline_consortium,      file.path(tables_dir, "headline_consortium.csv"))
write_excel_csv(headline_consortium_noF4, file.path(tables_dir, "headline_consortium_noFTIR4.csv"))

#### 19. Text reports ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â copy-paste material for the manuscript               ####
fmt_num <- function(x, d = 1) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "NA")

# --- Report 1: V10 analysis configuration ---
r1 <- init_report("01_config.txt",
                  "Ringversuche V10 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â analysis configuration")
write_section(r1, "Thresholds",
              sprintf("MIN_DELTA_CO2 = %g mg m-3\nMIN_Q_VENT = %g, MAX_Q_VENT = %g m3 h-1 LU-1",
                      MIN_DELTA_CO2, MIN_Q_VENT, MAX_Q_VENT))
write_section(r1, "Outdoor-line choice",
              sprintf("Retained: %s\nDropped:  %s", RETAINED_OUTDOOR, DROPPED_OUTDOOR))

# --- Report 2: animal-count summary (R2.22) ---
r2 <- init_report("02_animal_counts.txt",
                  "Campaign in-barn animal counts (R2.22)")
write_section(r2, "Counts",
              paste(sprintf("n_hours = %d", animal_summary$n_hours),
                    sprintf("mean    = %.1f cows", animal_summary$n_mean),
                    sprintf("median  = %.1f cows", animal_summary$n_median),
                    sprintf("min     = %.0f cows", animal_summary$n_min),
                    sprintf("max     = %.0f cows", animal_summary$n_max),
                    sprintf("hours with <30 cows: %d", animal_summary$n_hours_lt_30),
                    sprintf("hours with <40 cows: %d", animal_summary$n_hours_lt_40),
                    sep = "\n"))

# --- Report 3: NE vs SW outdoor-line decision (R2.24) ---
r3 <- init_report("03_outdoor_line_choice.txt",
                  "NE vs SW outdoor-line analysis (R2.24)")
write_section(r3, "Paired-test summary (campaign-wide)",
              paste(capture.output(print(
                      ne_sw_tests %>% select(gas, analyser, n, mean_NE, mean_SW,
                                             RPD_pct, p_wilcox_holm, significant) %>%
                              as.data.frame(),
                      row.names = FALSE)),
                    collapse = "\n"))
write_section(r3, "Wind-sector contamination class for each line",
              paste(capture.output(print(as.data.frame(sector_class), row.names = FALSE)),
                    collapse = "\n"))
write_section(r3, "Time fraction each line spends in each contamination class",
              paste(capture.output(print(
                      line_time_class %>% as.data.frame(), row.names = FALSE)),
                    collapse = "\n"))
write_section(r3, "Decision",
              sprintf("Retained outdoor: %s (less time downwind of LVAT + neighbour sources)\nDropped outdoor:  %s",
                      RETAINED_OUTDOOR, DROPPED_OUTDOOR))

# --- Report 4: Q and emission headline numbers ---
r4 <- init_report("04_q_e_headlines.txt",
                  "Headline Q + emission numbers, retained outdoor only")
write_section(r4, "Consortium mean and 95% CI (all 7 analysers)",
              paste(capture.output(print(as.data.frame(headline_consortium),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r4, "Consortium robust headline excluding FTIR.4 (high SD even post-library)",
              paste(capture.output(print(as.data.frame(headline_consortium_noF4),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r4, "Per-analyser means",
              paste(capture.output(print(as.data.frame(headline_stats),
                                          row.names = FALSE)),
                    collapse = "\n"))

# --- Report 5: BA + pairwise summary ---
r5 <- init_report("05_pairwise_and_BA.txt",
                  "Bland-Altman + pairwise regression / Lin's CCC summary")
write_section(r5, "Bland-Altman (intra-lab pairs, retained outdoor)",
              paste(capture.output(print(as.data.frame(ba_table),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r5, "Pairwise regression + Lin's CCC (mean per variable group)",
              paste(capture.output(print(
                      pairwise_tbl %>%
                              group_by(var) %>%
                              summarise(mean_pearson = mean(pearson_r, na.rm=TRUE),
                                        mean_ccc     = mean(ccc, na.rm=TRUE),
                                        mean_slope   = mean(deming_slope, na.rm=TRUE),
                                        .groups = "drop") %>%
                              as.data.frame(), row.names = FALSE)),
                    collapse = "\n"))

# --- Report 6: FTIR.4_old library check ---
r6 <- init_report("06_FTIR4_library_check.txt",
                  "FTIR.4_old vs FTIR.4 library re-evaluation (R2.7)")
write_section(r6, "Per gas x location",
              paste(capture.output(print(as.data.frame(lib_compare),
                                          row.names = FALSE)),
                    collapse = "\n"))

# --- Report 7: campaign overview + preprocessing loss ---
r7 <- init_report("07_campaign_overview_and_loss.txt",
                  "Campaign overview and preprocessing loss")
write_section(r7, "Campaign data overview",
              paste(capture.output(print(as.data.frame(campaign_overview),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r7, "Pre-hourly outlier-removal loss",
              paste(capture.output(print(as.data.frame(pre_hourly_loss),
                                          row.names = FALSE)),
                    collapse = "\n"))

# --- Report 8: plausibility diagnostics ---
r8 <- init_report("08_plausibility_diagnostics.txt",
                  "Plausibility diagnostics")
write_section(r8, "FTIR.4_old versus corrected FTIR.4",
              paste(capture.output(print(as.data.frame(lib_compare),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r8, "FTIR.2 NH3 offset note",
              paste("FTIR.2 showed a positive NH3 concentration offset.",
                    "The offset was not corrected directly because it largely cancelled in delta NH3 and",
                    "therefore did not propagate in the same way to Q and eNH3. See plots",
                    "."))

# --- Report 9: delta-quality screening ---
r9 <- init_report("09_delta_quality_screening.txt",
                  "Delta-quality screening")
write_section(r9, "Negative-delta masking summary",
              paste(capture.output(print(as.data.frame(stage2_dropouts),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r9, "Delta-quality summary",
              paste(capture.output(print(as.data.frame(delta_quality),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r9, "Q-filter dropout summary",
              paste(capture.output(print(as.data.frame(qvent_dropouts),
                                          row.names = FALSE)),
                    collapse = "\n"))

# --- Report 10: H0,1 analyser-only agreement ---
r10 <- init_report("10_H01_analyser_only_agreement.txt",
                   "H0,1 analyser-only agreement")
write_section(r10, "Bland-Altman summary",
              paste(capture.output(print(as.data.frame(ba_table),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r10, "Selected pairwise regression and concordance metrics",
              paste(capture.output(print(
                      pairwise_tbl %>%
                              filter((analyser_x == "FTIR.1" & analyser_y == "CRDS.1") |
                                     (analyser_x == "FTIR.2" & analyser_y == "CRDS.2")) %>%
                              as.data.frame(),
                      row.names = FALSE)),
                    collapse = "\n"))

# --- Report 11: H0,2 full-configuration agreement ---
r11 <- init_report("11_H02_full_configuration_agreement.txt",
                   "H0,2 full-configuration agreement")
write_section(r11, "Tukey HSD concentrations",
              paste(capture.output(print(as.data.frame(tukey_abs),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r11, "Tukey HSD delta concentrations",
              paste(capture.output(print(as.data.frame(tukey_delta),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r11, "Tukey HSD ventilation and emissions",
              paste(capture.output(print(as.data.frame(tukey_qe),
                                          row.names = FALSE)),
                    collapse = "\n"))

#### 20. Console summary                                                      ####
cat("\n===== V10 run summary =====\n")
cat("Retained outdoor line: ", RETAINED_OUTDOOR, "\n", sep = "")
cat("Dropped outdoor line:  ", DROPPED_OUTDOOR,  "\n", sep = "")
cat("Tables:  ", tables_dir, "\n", sep = "")
cat("Plots:   ", plots_dir,  "\n", sep = "")
cat("Reports: ", report_dir, "\n", sep = "")

