round_to_interval <- function(datetime, interval_sec = 450) {
        # Step 1: round to nearest interval
        dt_num <- as.numeric(datetime)
        rounded <- round(dt_num / interval_sec) * interval_sec
        rounded <- as.POSIXct(rounded, origin = "1970-01-01", tz = tz(datetime))
        
        # Step 2: if exactly at full hour, shift back by 1 sec
        full_hour <- format(rounded, "%M:%S") == "00:00"
        rounded[full_hour] <- rounded[full_hour] - 1
        
        return(rounded)
}
