# Manuscript 3, Campaign 1 recovered-data analysis (version 1).
# Describes concentration/ratio heterogeneity, sampling-density error, temporal
# aggregation sensitivity, and a four-bin Shannon-entropy adaptation of Chen et al.

library(data.table)
library(ggplot2)

workflow <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/mapping_high_resolution"
input <- file.path(workflow, "clean_data", "1_campaign", "recovered_crds_v01",
  "campaign1_recovered_four_analyser_CRDS_scale_corrected_v02.csv")
out <- file.path(workflow, "clean_data", "manuscript3_campaign1_recovered_v01")
fig <- file.path(workflow, "plots", "manuscript3_campaign1_recovered_v01")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

x <- fread(input, colClasses = c(DATE.TIME = "character", location = "character"))
x[, DATE.TIME := as.POSIXct(DATE.TIME, "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")]
x[, location_n := suppressWarnings(as.integer(location))]
x <- x[!is.na(location_n) & location_n %between% c(1L, 51L)]
x[CO2_corr < 300, CO2_corr := NA_real_]
x[CH4_corr <= 0, CH4_corr := NA_real_]
x[NH3_corr <= 0, NH3_corr := NA_real_]
x[, height := factor(c("top", "mid", "bottom")[(location_n - 1L) %% 3L + 1L],
                    levels = c("bottom", "mid", "top"))]
x[, horizontal_position := ceiling(location_n / 3)]
x[, `:=`(CH4_CO2 = CH4_corr / CO2_corr,
         NH3_CO2 = NH3_corr / CO2_corr,
         NH3_CH4 = NH3_corr / CH4_corr)]

metrics <- c("CO2_corr", "CH4_corr", "NH3_corr", "CH4_CO2", "NH3_CO2", "NH3_CH4")
labels <- c(CO2_corr="CO2", CH4_corr="CH4", NH3_corr="NH3",
            CH4_CO2="CH4/CO2", NH3_CO2="NH3/CO2", NH3_CH4="NH3/CH4")

# Location summaries retain robust and conventional dispersion statistics.
long <- melt(x, id.vars = c("DATE.TIME", "analyser", "location_n", "height",
  "horizontal_position"), measure.vars = metrics, variable.name = "metric",
  value.name = "value")
summary_location <- long[is.finite(value), .(
  n = .N, mean = mean(value), median = median(value), sd = sd(value),
  mad = mad(value), q25 = quantile(value, .25), q75 = quantile(value, .75),
  cv_pct = 100 * sd(value) / abs(mean(value))
), by = .(metric, location_n, height, horizontal_position)]
summary_location[, metric_label := labels[metric]]
fwrite(summary_location, file.path(out, "Table_location_mean_median_SD_MAD_CV.csv"))

# Equal-width, four-bin Shannon entropy across two-hour observed blocks.
# H is normalised by log2(4), making its range 0-1. This adapts Chen et al.'s
# cross-case CFD entropy to repeated empirical time blocks.
x[, block_2h := as.POSIXct(floor(as.numeric(DATE.TIME) / 7200) * 7200,
  origin = "1970-01-01", tz = "Europe/Berlin")]
block <- x[, lapply(.SD, mean, na.rm = TRUE),
           by = .(block_2h, location_n, height, horizontal_position), .SDcols = metrics]
complete_blocks <- block[, .(n_locations = uniqueN(location_n)), by = block_2h][
  n_locations == 51L, block_2h]
block <- block[block_2h %in% complete_blocks]
block_long <- melt(block, id.vars = c("block_2h", "location_n", "height",
  "horizontal_position"), measure.vars = metrics, variable.name = "metric",
  value.name = "value")
entropy4 <- function(v) {
  v <- v[is.finite(v)]; if (length(v) < 8 || diff(range(v)) == 0) return(NA_real_)
  b <- cut(v, breaks = seq(min(v), max(v), length.out = 5), include.lowest = TRUE)
  p <- as.numeric(table(b)) / length(v); p <- p[p > 0]
  -sum(p * log2(p)) / 2
}
entropy <- block_long[, .(n_blocks = sum(is.finite(value)),
  entropy_normalised = entropy4(value)), by = .(metric, location_n, height,
                                                horizontal_position)]
entropy[, rank_high_to_low := frank(-entropy_normalised, ties.method = "min"), by = metric]
entropy[, metric_label := labels[metric]]
fwrite(entropy, file.path(out, "Table_Chen_adapted_four_bin_Shannon_entropy.csv"))

# Full-network baseline and progressively reduced spatial configurations.
sets <- list("Full 51"=1:51,
  "Top + bottom 34"=setdiff(1:51, seq(2,50,3)),
  "Bottom 17"=seq(3,51,3), "Every second column, all heights 27"=unlist(
    lapply(seq(1,17,2), function(h) (3*h-2):(3*h))))
base <- block[, lapply(.SD, mean, na.rm = TRUE), by = block_2h, .SDcols = metrics]
cfg <- rbindlist(lapply(names(sets), function(nm) {
  z <- block[location_n %in% sets[[nm]], lapply(.SD, mean, na.rm = TRUE),
             by = block_2h, .SDcols = metrics]
  z[, configuration := nm]; z
}))
cfg_long <- melt(cfg, id.vars = c("block_2h", "configuration"),
                 measure.vars = metrics, variable.name = "metric", value.name = "estimate")
base_long <- melt(base, id.vars = "block_2h", measure.vars = metrics,
                  variable.name = "metric", value.name = "baseline")
cfg_long <- merge(cfg_long, base_long, by = c("block_2h", "metric"))
cfg_long[, relative_error_pct := 100 * (estimate - baseline) / baseline]
fwrite(cfg_long, file.path(out, "campaign1_spatial_configuration_block_results.csv"))
cfg_sum <- cfg_long[, .(n_blocks=.N, mean_RE_pct=mean(relative_error_pct,na.rm=TRUE),
  median_RE_pct=median(relative_error_pct,na.rm=TRUE),
  mean_absolute_RE_pct=mean(abs(relative_error_pct),na.rm=TRUE),
  p95_absolute_RE_pct=quantile(abs(relative_error_pct),.95,na.rm=TRUE)),
  by=.(metric,configuration)]
fwrite(cfg_sum, file.path(out, "Table_spatial_configuration_relative_error.csv"))

# Temporal aggregation sensitivity: spatial CV at 2, 4, 8 and 24 h.
temporal <- rbindlist(lapply(c(2,4,8,24), function(hours) {
  z <- copy(x); z[, time_block := as.POSIXct(floor(as.numeric(DATE.TIME)/(hours*3600)) *
    hours*3600, origin="1970-01-01", tz="Europe/Berlin")]
  q <- z[, lapply(.SD, mean, na.rm=TRUE), by=.(time_block,location_n), .SDcols=metrics]
  ql <- melt(q,id.vars=c("time_block","location_n"),measure.vars=metrics,
             variable.name="metric",value.name="value")
  ql[,.(spatial_cv_pct=100*sd(value,na.rm=TRUE)/abs(mean(value,na.rm=TRUE))),
     by=.(time_block,metric)][,hours:=hours]
}))
fwrite(temporal, file.path(out, "campaign1_temporal_aggregation_spatial_CV.csv"))
fwrite(temporal[,.(median_spatial_CV_pct=median(spatial_cv_pct,na.rm=TRUE),
  q25=quantile(spatial_cv_pct,.25,na.rm=TRUE),q75=quantile(spatial_cv_pct,.75,na.rm=TRUE)),
  by=.(metric,hours)], file.path(out,"Table_temporal_aggregation_sensitivity.csv"))

# Publication figures.
p1 <- ggplot(summary_location, aes(horizontal_position, mean, colour=height,
  group=height)) + geom_line() + geom_point(aes(size=pmin(cv_pct,100))) +
  facet_wrap(~metric_label, scales="free_y", ncol=3) +
  scale_colour_manual(values=c(bottom="steelblue1",mid="green3",top="orange")) +
  scale_size_continuous("CV (%)", range=c(1.5,5)) + theme_bw() +
  labs(x="Horizontal position",y="Campaign mean",colour="Height")
ggsave(file.path(fig,"Fig1_spatial_means_with_CV.png"),p1,width=12,height=7,dpi=300)

p2 <- ggplot(entropy, aes(horizontal_position, height, fill=entropy_normalised)) +
  geom_tile(colour="white") + geom_text(aes(label=location_n),size=2.6) +
  facet_wrap(~metric_label,ncol=3) + scale_fill_viridis_c(limits=c(0,1)) +
  theme_bw() + labs(x="Horizontal position",y="Height",
    fill="Normalised\nShannon entropy")
ggsave(file.path(fig,"Fig2_Chen_adapted_entropy_map.png"),p2,width=12,height=7,dpi=300)

p3 <- ggplot(cfg_long[configuration!="Full 51"],
  aes(configuration,relative_error_pct,fill=configuration)) +
  geom_hline(yintercept=0,lty=2)+geom_boxplot(outlier.alpha=.15)+
  facet_wrap(~labels[metric],scales="free_y",ncol=3)+theme_bw()+
  theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="none")+
  labs(x=NULL,y="Relative error against full 51-location mean (%)")
ggsave(file.path(fig,"Fig3_sampling_configuration_relative_error.png"),p3,
       width=12,height=7,dpi=300)

agree <- cfg_long[configuration == "Top + bottom 34" & metric %in%
                    c("CO2_corr","CH4_corr","NH3_corr")]
p4 <- ggplot(agree, aes(baseline, estimate)) + geom_abline(slope=1,intercept=0) +
  geom_point(alpha=.45,size=.8) + facet_wrap(~labels[metric],scales="free",nrow=1) +
  theme_bw() + labs(x="Full 51-location mean (ppm)",
    y="Top-and-bottom 34-location mean (ppm)")
ggsave(file.path(fig,"Fig4_full_vs_top_bottom_agreement.png"),p4,
       width=12,height=4.5,dpi=300)

writeLines(c(
  "Campaign 1 recovered-data analysis v01",
  "CRDS cycles: ordered 1:16 sequences only; first 60 s flushed.",
  "Four-minute dwell: remaining 180 s averaged; two 180-s exceptions retained for sensitivity.",
  "FTIR correction: CO2/1.06, CH4/1.06, NH3/1.09.",
  "Entropy: four equal-width bins, Shannon bits divided by log2(4).",
  "This is an empirical time-block adaptation of Chen et al., not CFD case entropy.",
  "Relative error = 100 * (configuration mean - full-network mean) / full-network mean."
), file.path(out,"analysis_readme.txt"))
