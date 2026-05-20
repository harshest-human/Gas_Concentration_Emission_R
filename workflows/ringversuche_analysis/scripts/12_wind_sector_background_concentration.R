library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

# Set up paths
data_file <- "../result_data/tables/Version_8/20250408-15_ringversuche_input_combined_data.csv"
output_dir <- "../result_data/plots/wind_sector_background"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load data
cat("Loading combined dataset...\n")
data <- read.csv(data_file)
data$DATE.TIME <- ymd_hms(data$DATE.TIME)

# Filter for the 7-day campaign period
filtered_data <- data %>%
  filter(DATE.TIME >= ymd_hms("2025-04-08 12:00:00"),
         DATE.TIME <= ymd_hms("2025-04-14 12:00:00"))

# 2. Extract wind sector
deg_to_compass8 <- function(direction_deg) {
  labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  width <- 360 / length(labels)
  idx <- floor(((direction_deg %% 360) + width / 2) / width) %% length(labels) + 1
  factor(labels[idx], levels = labels)
}

processed_data <- filtered_data %>%
  filter(!is.na(wd_mst)) %>%
  mutate(wind_sector = deg_to_compass8(wd_mst))

# 3. Aggregate data over the entire 7-day period by wind sector
sector_summary <- processed_data %>%
  group_by(wind_sector) %>%
  summarise(
    CO2_N = mean(CO2_ppm_NE, na.rm = TRUE),
    CO2_S = mean(CO2_ppm_SW, na.rm = TRUE),
    CH4_N = mean(CH4_ppm_NE, na.rm = TRUE),
    CH4_S = mean(CH4_ppm_SW, na.rm = TRUE),
    NH3_N = mean(NH3_ppm_NE, na.rm = TRUE),
    NH3_S = mean(NH3_ppm_SW, na.rm = TRUE),
    .groups = "drop"
  )

# 4. Pivot data for plotting
plot_data <- sector_summary %>%
  pivot_longer(
    cols = -wind_sector,
    names_to = c("Gas", "Location"),
    names_sep = "_",
    values_to = "Concentration"
  ) %>%
  mutate(
    Gas = factor(Gas, levels = c("CO2", "CH4", "NH3")),
    Location = factor(Location, levels = c("N", "S"))
  )

# 5. Create Plot
cat("Generating faceted plot for background concentrations by wind sector...\n")

p <- ggplot(plot_data, aes(x = wind_sector, y = Concentration, fill = Location)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  facet_wrap(~Gas, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("N" = "#377EB8", "S" = "#E41A1C"), 
                    labels = c("N" = "Background North (NE)", "S" = "Background South (SW)")) +
  labs(
    title = "Average Background Gas Concentrations by Wind Sector",
    subtitle = "Campaign period: 08.04.2025 12:00 to 14.04.2025 12:00",
    x = "Wind Sector",
    y = "Average Concentration (ppm)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80")
  )

# Save the plot
file_name <- file.path(output_dir, "background_concentration_by_wind_sector.png")
ggsave(file_name, plot = p, width = 10, height = 12, dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", file_name))
