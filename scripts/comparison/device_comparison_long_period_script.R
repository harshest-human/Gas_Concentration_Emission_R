####### libraries ######
library(tidyverse)
library(reshape2)
library(hablar)
library(lubridate)
library(psych)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)
library(readr)
library(scales)
library(purrr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
project_dir <- if (length(script_path) > 0) {
        normalizePath(file.path(dirname(script_path[1]), "..", ".."), winslash = "/", mustWork = FALSE)
} else {
        normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source(file.path(project_dir, "scripts", "utils", "remove_outliers_function.R"))
source(file.path(project_dir, "scripts", "utils", "indirect.CO2.balance function.R"))

result_dir <- file.path(project_dir, "workflows", "device_comparison", "result_data")
plot_dir <- file.path(result_dir, "plots")
table_dir <- file.path(result_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

####### Import CRDS dataset #######
# List all crds files
crds_files <- list.files(file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg"),
                                       pattern = "\\.csv$",
                                       full.names = TRUE,
                                       recursive = TRUE)

                                       
# Read, Filter Date Time and Calculate Delta
crds_data <- purrr::map_dfr(crds_files,read.csv) %>% 
        mutate(DATE.HOUR = ymd_hms(DATE.HOUR),        
               delta_CO2 = CO2_in - CO2_S,       
               delta_CH4 = CH4_in - CH4_S,
               delta_NH3 = NH3_in - NH3_S,
               analyzer = "crds") %>%
        filter(DATE.HOUR >= "2025-10-01 00:00:00", DATE.HOUR <= "2026-03-30 23:00:00") %>%
        select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) %>% 
        remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"),
                        group_cols = c("DATE.HOUR"))


####### Import logas_ndir dataset #######
# List all logas_ndir files
logas_ndir_files <- list.files(
        path = file.path(project_dir, "workflows", "device_comparison", "raw_data", "logas_ndir_raw", "Messdaten"),
        pattern = "^differenzmessung_.*\\.txt$",
        full.names = TRUE)

# Function to safely read one file
logas_ndir_read_log_file <- function(file) {
        message("Reading: ", basename(file))
        df <- tryCatch(
                read.table(
                        file,
                        header = TRUE,
                        sep = "\t",
                        dec = ",",
                        check.names = FALSE,
                        stringsAsFactors = FALSE,
                        fill = TRUE,
                        comment.char = "",
                        colClasses = "character"   
                ),
                error = function(e) {
                        message("Error reading ", basename(file), ": ", e$message)
                        return(NULL)
                }
        )
        return(df)
}

# Read and combine all files safely
logas_ndir_data  <- logas_ndir_files %>%
        lapply(logas_ndir_read_log_file) %>%
        bind_rows() %>%
        mutate(DATE.TIME = dmy_hms(`Datum Uhrzeit`)) %>%
        mutate(DATE.HOUR = floor_date(DATE.TIME, "hour")) %>%
        group_by(DATE.HOUR) %>%
        summarise(delta_CH4 = mean(as.numeric(gsub(",", ".", `CH4 in ppm`)), na.rm = TRUE),
                  delta_CO2 = mean(as.numeric(gsub(",", ".", `CO2 in ppm`)), na.rm = TRUE),
                  delta_N2O = mean(as.numeric(gsub(",", ".", `N2O in ppm`)), na.rm = TRUE),
                  delta_NH3 = mean(as.numeric(gsub(",", ".", `NH3 in ppm`)), na.rm = TRUE)) %>%
        mutate(analyzer = "logas_ndir") %>%
        filter(DATE.HOUR >= "2025-10-01 00:00:00", DATE.HOUR <= "2026-03-30 23:00:00") %>%
        select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) %>% 
        remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"),
                        group_cols = c("DATE.HOUR"))



####### Import logas_tdlas dataset #########
# List all logas_tdlas files
logas_tdlas_files <- list.files(path = file.path(project_dir, "workflows", "device_comparison", "raw_data", "logas_tdlas_raw"), 
                    pattern = "^[^~].*\\.xlsx$", 
                    full.names = TRUE)

# Read all files and bind them row-wise
logas_tdlas_data <- logas_tdlas_files %>%
        map_df(~ read_excel(.x)) %>%
        mutate(Time = ymd_hms(Time),
               DATE.HOUR = floor_date(Time, "hour")) %>%
        pivot_longer(cols = any_of(c("CH4", "NH3", "CO2")),
                     names_to = "gas", values_to = "value") %>%
        mutate(gas_name = paste0(gas, ifelse(Type == 1, "_in", "_S"))) %>%
        group_by(DATE.HOUR, gas_name) %>%
        summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>% 
        pivot_wider(names_from = gas_name, values_from = value) %>%
        arrange(DATE.HOUR) %>%
        mutate(DATE.HOUR = as.POSIXct(DATE.HOUR),        
               delta_CO2 = CO2_in - CO2_S,       
               delta_CH4 = CH4_in - CH4_S,
               delta_NH3 = NH3_in - NH3_S,
               analyzer = "logas_tdlas")  %>%
        filter(DATE.HOUR >= "2025-10-01 00:00:00", DATE.HOUR <= "2026-03-30 23:00:00") %>%
        select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) %>% 
        remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"),
                        group_cols = c("DATE.HOUR"))

####### Combine crds, logas_ndir, and logas_tdlas #######
Gas_data <- bind_rows(logas_ndir_data, crds_data, logas_tdlas_data) %>% arrange(DATE.HOUR)
Gas_data <- remove_outliers(Gas_data, exclude_cols = c("DATE.HOUR", "analyzer"), group_cols = NULL)

####### Calculate ventilation rate and Emissions ######
# Import animal data
animal_data <- read_excel(file.path(project_dir, "shared_data", "clean_data", "animal_clean", "animal_data_2025-01-10_2026-01-19.xlsx")) 

# Import temperature data
temp_data <- list.files(file.path(project_dir, "shared_data", "clean_data", "temp_rh_clean"),
                        pattern = "\\.csv$",
                        full.names = TRUE,
                        recursive = TRUE) %>%
        purrr::map_dfr(~ readr::read_csv(.x, show_col_types = FALSE)) %>%
        dplyr::select(Date, T_inside) %>%
        dplyr::rename(temp_in = T_inside) %>%
        dplyr::mutate(DATE.TIME = lubridate::mdy_hms(sub(" \\+0000$", "", Date)),
                      DATE.HOUR = lubridate::floor_date(DATE.TIME, "hour")) %>%
        dplyr::group_by(DATE.HOUR) %>%
        dplyr::summarise(temp_in = mean(temp_in, na.rm = TRUE),
                         .groups = "drop")

# Combine the Input data
input_data <- Gas_data %>%
        left_join(animal_data, by = "DATE.HOUR", relationship = "many-to-many") %>%
        left_join(temp_data,   by = "DATE.HOUR", relationship = "many-to-many") %>%
        filter(DATE.HOUR >= ymd_hms("2026-03-10 00:00:00"),
               DATE.HOUR <= ymd_hms("2026-03-30 23:00:00")) %>%
        rename("DATE.TIME" = "DATE.HOUR")

# Calculate Emissions
emission_data <- indirect.CO2.balance(input_data)

emission_reshaped <-reshaper(emission_data)

#write_excel_csv(emission_reshaped, "emission_reshaped_20251209_20251222.csv")

#write_excel_csv(emission_data, "emission_data_20251209_20251222.csv")

####### Data Visualization ########
d_errorbarplot <- emierrorbarplot(emission_reshaped, y = c("delta_CO2", "delta_CH4", "delta_NH3"))

d_trendplot <- emitrendplot(emission_reshaped,
                            y = c("delta_CO2", "delta_CH4", "delta_NH3"))

q_e_errorbarplot <- emierrorbarplot(data = emission_reshaped,
                                    y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"))

q_e_trend_plot <- emitrendplot(data = emission_reshaped,
                               y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"))

####### Statistical Data Analysis #######
emi_daily <- read.csv(file.path(table_dir, "emission_data_20251209_20251222.csv")) %>%
        mutate(
                DATE.TIME = ymd_hms(DATE.TIME),
                DATE = as.Date(DATE.TIME)
        ) %>%
        group_by(DATE, analyzer) %>%
        summarise(
                Q_vent = mean(Q_vent, na.rm = TRUE),
                e_CH4_ghLU = mean(e_CH4_ghLU, na.rm = TRUE),
                e_NH3_ghLU = mean(e_NH3_ghLU, na.rm = TRUE),
                .groups = "drop"
        )


emi_compare <- emi_daily %>%
        pivot_wider(
                names_from = analyzer,
                values_from = c(Q_vent, e_CH4_ghLU, e_NH3_ghLU)
        ) %>%
        mutate(
                RE_Qvent_logas_tdlas   = (Q_vent_logas_tdlas   - Q_vent_crds) / Q_vent_crds * 100,
                RE_Qvent_logas_ndir = (Q_vent_logas_ndir - Q_vent_crds) / Q_vent_crds * 100,
                
                RE_CH4_logas_tdlas     = (e_CH4_ghLU_logas_tdlas    - e_CH4_ghLU_crds) / e_CH4_ghLU_crds * 100,
                RE_CH4_logas_ndir   = (e_CH4_ghLU_logas_ndir  - e_CH4_ghLU_crds) / e_CH4_ghLU_crds * 100,
                
                RE_NH3_logas_tdlas     = (e_NH3_ghLU_logas_tdlas    - e_NH3_ghLU_crds) / e_NH3_ghLU_crds * 100,
                RE_NH3_logas_ndir   = (e_NH3_ghLU_logas_ndir  - e_NH3_ghLU_crds) / e_NH3_ghLU_crds * 100
        )


emi_analyzer <- emi_compare %>%
        summarise(Q_vent_crds = mean(Q_vent_crds, na.rm = TRUE),
                  Q_vent_logas_tdlas = mean(Q_vent_logas_tdlas, na.rm = TRUE),
                  Q_vent_logas_ndir = mean(Q_vent_logas_ndir, na.rm = TRUE),
                  e_CH4_ghLU_crds = mean(e_CH4_ghLU_crds, na.rm = TRUE),
                  e_CH4_ghLU_logas_tdlas = mean(e_CH4_ghLU_logas_tdlas, na.rm = TRUE),
                  e_CH4_ghLU_logas_ndir = mean(e_CH4_ghLU_logas_ndir, na.rm = TRUE),
                  e_NH3_ghLU_crds = mean(e_NH3_ghLU_crds, na.rm = TRUE),
                  e_NH3_ghLU_logas_tdlas = mean(e_NH3_ghLU_logas_tdlas, na.rm = TRUE),
                  e_NH3_ghLU_logas_ndir = mean(e_NH3_ghLU_logas_ndir, na.rm = TRUE),
                  RE_Qvent_logas_tdlas = mean(RE_Qvent_logas_tdlas),
                  RE_Qvent_logas_ndir = mean(RE_Qvent_logas_ndir),
                  RE_CH4_logas_tdlas = mean(RE_CH4_logas_tdlas),
                  RE_CH4_logas_ndir = mean(RE_CH4_logas_ndir),
                  RE_NH3_logas_tdlas = mean(RE_NH3_logas_tdlas),
                  RE_NH3_logas_ndir = mean(RE_NH3_logas_ndir),
                  .groups = "drop"
        )
                


library(broom)

models <- list(
        Qvent_logas_tdlas   = lm(Q_vent_logas_tdlas ~ Q_vent_crds, data = emi_compare),
        Qvent_logas_ndir = lm(Q_vent_logas_ndir ~ Q_vent_crds, data = emi_compare),
        CH4_logas_tdlas     = lm(e_CH4_ghLU_logas_tdlas ~ e_CH4_ghLU_crds, data = emi_compare),
        CH4_logas_ndir   = lm(e_CH4_ghLU_logas_ndir ~ e_CH4_ghLU_crds, data = emi_compare),
        NH3_logas_tdlas     = lm(e_NH3_ghLU_logas_tdlas ~ e_NH3_ghLU_crds, data = emi_compare),
        NH3_logas_ndir   = lm(e_NH3_ghLU_logas_ndir ~ e_NH3_ghLU_crds, data = emi_compare)
)

sapply(models, function(x) summary(x)$r.squared)


