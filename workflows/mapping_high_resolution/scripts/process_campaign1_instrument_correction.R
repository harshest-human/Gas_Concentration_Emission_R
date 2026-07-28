###############################################################################
##### CAMPAIGN 1: HARMONISE FTIR CONCENTRATIONS TO THE CRDS SCALE
###############################################################################
#
# During the 24-hour co-located field comparison, FTIR concentrations were:
#   - 6% higher than CRDS for CO2;
#   - 6% higher than CRDS for CH4; and
#   - 9% higher than CRDS for NH3.
#
# Therefore, FTIR measurements are divided by 1.06 (CO2 and CH4) or 1.09
# (NH3). This is an instrument-scale harmonisation, not an absolute calibration.
#
# The original concentrations are preserved in *_raw columns. Corrected values
# are stored in *_corr columns. CRDS values are unchanged, and H2O is unchanged
# for every analyser because no H2O comparison factor was established.
#
# Location "s" is retained in the output and represents location 52, the
# south-background sampling point. It must be excluded later from analyses of
# spatial heterogeneity among the 51 internal barn locations.
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

library(data.table)

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)

campaign1_dir <- file.path(workflow_dir, "clean_data", "1_campaign")

input_path <- file.path(
  campaign1_dir,
  "20241001.000232_20241024.143858_FTIR1_FTIR2_CRDS8_CRDS9_combined.csv"
)

output_path <- file.path(
  campaign1_dir,
  paste0(
    "20241001.000232_20241024.143858_",
    "campaign1_raw_and_CRDS_scale_corrected.csv"
  )
)

campaign1 <- fread(
  input_path,
  colClasses = list(
    character = c("DATE.TIME", "analyser", "location")
  ),
  showProgress = FALSE
)

required_columns <- c(
  "DATE.TIME", "analyser", "location", "CO2", "CH4", "NH3", "H2O"
)
missing_columns <- setdiff(required_columns, names(campaign1))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

campaign1[, DATE.TIME := as.POSIXct(
  DATE.TIME,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)]

if (anyNA(campaign1[["DATE.TIME"]])) {
  stop("One or more DATE.TIME values could not be parsed.")
}

if (!all(campaign1[["analyser"]] %in% c("FTIR1", "FTIR2", "CRDS8", "CRDS9"))) {
  stop(
    "Unexpected analyser(s): ",
    paste(
      setdiff(unique(campaign1[["analyser"]]), c("FTIR1", "FTIR2", "CRDS8", "CRDS9")),
      collapse = ", "
    )
  )
}

# Preserve the source concentrations.
setnames(
  campaign1,
  old = c("CO2", "CH4", "NH3", "H2O"),
  new = c("CO2_raw", "CH4_raw", "NH3_raw", "H2O_raw")
)

# Begin with corrected values equal to the source values.
campaign1[, `:=`(
  CO2_corr = CO2_raw,
  CH4_corr = CH4_raw,
  NH3_corr = NH3_raw,
  H2O_corr = H2O_raw
)]

# Harmonise only FTIR1 and FTIR2 to the CRDS measurement scale.
campaign1[
  analyser %in% c("FTIR1", "FTIR2"),
  `:=`(
    CO2_corr = CO2_raw / 1.06,
    CH4_corr = CH4_raw / 1.06,
    NH3_corr = NH3_raw / 1.09
  )
]

campaign1[, correction_applied := as.integer(
  analyser %in% c("FTIR1", "FTIR2")
)]

campaign1[, correction_basis := fifelse(
  correction_applied == 1L,
  "FTIR harmonised to CRDS: CO2/1.06; CH4/1.06; NH3/1.09",
  "CRDS retained without correction"
)]

setcolorder(
  campaign1,
  c(
    "DATE.TIME", "analyser", "location",
    "CO2_raw", "CH4_raw", "NH3_raw", "H2O_raw",
    "CO2_corr", "CH4_corr", "NH3_corr", "H2O_corr",
    "correction_applied", "correction_basis"
  )
)
setorder(campaign1, DATE.TIME, analyser, location)

if (!"s" %in% campaign1[["location"]]) {
  stop("South-background location 's' is missing from Campaign 1.")
}

fwrite(
  campaign1,
  output_path,
  quote = TRUE,
  dateTimeAs = "write.csv"
)

message("Campaign 1 rows written: ", nrow(campaign1))
message(
  "FTIR rows corrected: ",
  campaign1[["correction_applied"]] |> sum(na.rm = TRUE)
)
message(
  "South-background rows retained: ",
  sum(campaign1[["location"]] == "s", na.rm = TRUE)
)
message("Output: ", output_path)
