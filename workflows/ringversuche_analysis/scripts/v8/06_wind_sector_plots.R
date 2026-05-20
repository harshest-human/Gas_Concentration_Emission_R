######################################################################
# 06_wind_sector_plots.R   (v8)
# ---------------------------------------------------------------------
# Generates new gas concentration graphs for NE and SW sampling points,
# sorted and categorized by wind direction sectors (8 compass sectors).
######################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(ggplot2)
})

v8_base       <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
v8_v6_tables  <- file.path(v8_base, "result_data/tables/Version_6")
v8_plots_out  <- file.path(v8_base, "result_data/plots/Version_8")
dir.create(v8_plots_out, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Load Reshaped Data -----------------------------------------
message("Loading concentration data...")
cr <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_concentration_reshaped.csv"),
               show_col_types = FALSE, guess_max = 50000)

message("Loading emission data to extract wind direction...")
er <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_emission_reshaped.csv"),
               show_col_types = FALSE, guess_max = 50000)

# ---- 2) Extract Wind Data & Merge ----------------------------------
# Extract mast wind direction (wd_mst)
wd_data <- er %>%
  filter(var == "wd_mst") %>%
  select(DATE.TIME, wd_mst = value) %>%
  distinct() %>%
  filter(!is.na(wd_mst))

# Merge wind data into concentration data
cr_wind <- cr %>%
  inner_join(wd_data, by = "DATE.TIME")

# ---- 3) Define Wind Sectors ----------------------------------------
cr_wind <- cr_wind %>%
  mutate(wind_sector = case_when(
    is.na(wd_mst) ~ NA_character_,
    (wd_mst >= 337.5 | wd_mst < 22.5) ~ "N",
    (wd_mst >= 22.5 & wd_mst < 67.5) ~ "NE",
    (wd_mst >= 67.5 & wd_mst < 112.5) ~ "E",
    (wd_mst >= 112.5 & wd_mst < 157.5) ~ "SE",
    (wd_mst >= 157.5 & wd_mst < 202.5) ~ "S",
    (wd_mst >= 202.5 & wd_mst < 247.5) ~ "SW",
    (wd_mst >= 247.5 & wd_mst < 292.5) ~ "W",
    (wd_mst >= 292.5 & wd_mst < 337.5) ~ "NW"
  )) %>%
  mutate(wind_sector = factor(wind_sector, levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")))

# Relabel locations to match V8 naming convention
cr_wind <- cr_wind %>%
  mutate(location = case_when(
    location == "North background" ~ "NE background",
    location == "South background" ~ "SW background",
    TRUE ~ location
  ))

# ---- 4) Filter for NE & SW points and specific variables -----------
target_vars <- c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3")

plot_data <- cr_wind %>%
  filter(location %in% c("NE background", "SW background"),
         var %in% target_vars,
         !is.na(wind_sector))

# Clean up variables for plot faceting
facet_labels_expr <- c(
  "CO2_mgm3" = "c[CO2]~'(mg '*m^-3*')'",
  "CH4_mgm3" = "c[CH4]~'(mg '*m^-3*')'",
  "NH3_mgm3" = "c[NH3]~'(mg '*m^-3*')'"
)

plot_data <- plot_data %>%
  mutate(facet_label = factor(facet_labels_expr[as.character(var)], levels = facet_labels_expr[target_vars]))

# Define aesthetics matching the rest of the manuscript
analyzer_colors <- c(
  "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
  "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
  "CRDS.3" = "#a6761d", "baseline" = "black"
)

# Ensure all analyzers have a color mapping
all_analyzers <- unique(plot_data$analyzer)
analyzer_colors_full <- setNames(rep("black", length(all_analyzers)), all_analyzers)
analyzer_colors_full[names(analyzer_colors)] <- analyzer_colors

# ---- 5) Generate Plot ----------------------------------------------
message("Generating boxplots by wind sector...")

p <- ggplot(plot_data, aes(x = wind_sector, y = value, fill = analyzer)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8, lwd = 0.4) +
  scale_fill_manual(values = analyzer_colors_full) +
  facet_grid(facet_label ~ location, scales = "free_y", switch = "y", 
             labeller = labeller(facet_label = label_parsed, location = label_value)) +
  theme_classic() +
  labs(x = "Wind Direction Sector", y = NULL, fill = "Analyzer") +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text.y.left = element_text(size = 14, vjust = 0.5),
    strip.placement = "outside",
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom"
  )

out_file <- file.path(v8_plots_out, "c_boxplot_by_wind.png")
ggsave(out_file, p, width = 14, height = 9, dpi = 300)
message("Plot saved to: ", out_file)

# Produce a secondary plot showing counts of observations per sector
p_count <- ggplot(plot_data %>% group_by(location, wind_sector, var) %>% count() %>% filter(var == "CO2_mgm3"), 
                  aes(x = wind_sector, y = n, fill = location)) +
  geom_col(position = "dodge") +
  theme_classic() +
  labs(x = "Wind Direction Sector", y = "Number of Observations", fill = "Location") +
  theme(text = element_text(size = 14))

out_file_count <- file.path(v8_plots_out, "wind_sector_counts.png")
ggsave(out_file_count, p_count, width = 10, height = 6, dpi = 300)
message("Count plot saved to: ", out_file_count)

message("Script complete.")
