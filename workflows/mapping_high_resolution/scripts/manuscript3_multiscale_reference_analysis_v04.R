###############################################################################
# Manuscript 3 v04: multiscale spatial persistence and reference-location study
###############################################################################

library(data.table)
library(ggplot2)
library(lme4)
library(mgcv)
library(cluster)
library(rpart)

set.seed(20260803)
tz <- "Europe/Berlin"
workflow <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/mapping_high_resolution"
out <- file.path(workflow, "clean_data", "manuscript3_multiscale_reference_v04")
fig <- file.path(workflow, "plots", "manuscript3_multiscale_reference_v04")
models <- file.path(out, "models")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(fig, recursive=TRUE, showWarnings=FALSE)
dir.create(models, recursive=TRUE, showWarnings=FALSE)

gases <- c("CO2", "CH4", "NH3")
ratios <- c("CH4_CO2", "NH3_CO2", "NH3_CH4")
responses <- c(gases, ratios)
top <- seq(1L,49L,3L); mid <- seq(2L,50L,3L); bottom <- seq(3L,51L,3L)
height_cols <- c(top="orange",mid="green3",bottom="steelblue1")

parse_local <- function(z) as.POSIXct(as.character(z), "%Y-%m-%d %H:%M:%S", tz=tz)
height_of <- function(z) fifelse(z%in%top,"top",fifelse(z%in%mid,"mid","bottom"))
floor_time <- function(z, hours) as.POSIXct(floor(as.numeric(z)/(hours*3600))*
  hours*3600, origin="1970-01-01", tz=tz)

###############################################################################
##### DATA
###############################################################################

months <- c("2024-06","2024-07","2024-08")
c1_files <- file.path(workflow,"clean_data","1_campaign",
  "recovered_crds_june_august_v03",months,
  paste0("campaign1_recovered_four_analyser_CRDS_scale_corrected_",months,"_v03.csv"))
c1 <- rbindlist(lapply(c1_files, fread,
  colClasses=c(DATE.TIME="character",location="character")),fill=TRUE)
c1[,DATE.TIME:=parse_local(DATE.TIME)]
c1[,campaign:="Campaign 1"]

c2_file <- file.path(workflow,"clean_data","2_campaign",
  "20241116.000647_20241231.235710_CRDS8_CRDS9_campaign2_combined.csv")
c2 <- fread(c2_file,colClasses=c(DATE.TIME="character",location="character"))
c2[,DATE.TIME:=parse_local(DATE.TIME)]
c2[,campaign:="Campaign 2"]
for(g in gases) c2[,(paste0(g,"_corr")):=get(g)]

keep <- c("DATE.TIME","campaign","analyser","location",
          paste0(gases,"_corr"))
x <- rbindlist(list(c1[,..keep],c2[,..keep]),fill=TRUE)
x[,location_n:=suppressWarnings(as.integer(location))]
x <- x[!is.na(location_n)]
x[CO2_corr<300,CO2_corr:=NA_real_]
x[CH4_corr<=0,CH4_corr:=NA_real_]
x[NH3_corr<=0,NH3_corr:=NA_real_]
x[,`:=`(CO2=CO2_corr,CH4=CH4_corr,NH3=NH3_corr)]
x[,`:=`(CH4_CO2=CH4/CO2,NH3_CO2=NH3/CO2,NH3_CH4=NH3/CH4,
        height=factor(height_of(location_n),levels=c("bottom","mid","top")),
        horizontal_position=ceiling(location_n/3),date=as.Date(DATE.TIME,tz=tz),
        hour=as.numeric(format(DATE.TIME,"%H"))+
             as.numeric(format(DATE.TIME,"%M"))/60)]

validation <- x[,.(rows=.N,first=min(DATE.TIME),last=max(DATE.TIME),
  locations=uniqueN(location_n),dates=uniqueN(date)),by=campaign]
fwrite(validation,file.path(out,"Table_01_data_validation.csv"))

###############################################################################
##### MULTISCALE DESCRIPTIVES
###############################################################################

summarise_scale <- function(hours,label){
  z <- copy(x); z[,time_block:=floor_time(DATE.TIME,hours)]
  b <- z[,lapply(.SD,mean,na.rm=TRUE),by=.(campaign,time_block,location_n,height,
       horizontal_position),.SDcols=responses]
  long <- melt(b,id.vars=c("campaign","time_block","location_n","height",
    "horizontal_position"),measure.vars=responses,variable.name="response",
    value.name="value")
  long[is.finite(value),.(n=.N,mean=mean(value),median=median(value),
    minimum=min(value),p05=quantile(value,.05),p95=quantile(value,.95),
    maximum=max(value),sd=sd(value),mad=mad(value),cv_pct=100*sd(value)/abs(mean(value))),
    by=.(campaign,location_n,height,horizontal_position,response)][,scale:=label][]
}
descriptive <- rbindlist(list(summarise_scale(1,"hour"),summarise_scale(24,"day"),
  summarise_scale(168,"week")))

# Whole-campaign statistics use every valid observation at each location. They
# are deliberately calculated as a fourth scale rather than inferred from the
# weekly summaries.
campaign_descriptive <- melt(x,id.vars=c("campaign","location_n","height",
  "horizontal_position"),measure.vars=responses,variable.name="response",
  value.name="value")[is.finite(value),.(n=.N,mean=mean(value),median=median(value),
  minimum=min(value),p05=quantile(value,.05),p95=quantile(value,.95),
  maximum=max(value),sd=sd(value),mad=mad(value),
  cv_pct=100*sd(value)/abs(mean(value))),
  by=.(campaign,location_n,height,horizontal_position,response)][,scale:="campaign"]
descriptive <- rbind(descriptive,campaign_descriptive,fill=TRUE)
fwrite(descriptive,file.path(out,"Table_02_multiscale_location_descriptives.csv"))

scale_levels <- c("hour","day","week","campaign")
descriptive[,scale:=factor(scale,levels=scale_levels)]
scale_comparison <- descriptive[,.(locations=.N,
  median_of_location_means=median(mean),spatial_sd_of_location_means=sd(mean),
  median_of_location_medians=median(median),
  spatial_sd_of_location_medians=sd(median),
  median_location_SD=median(sd),median_location_MAD=median(mad),
  median_location_CV_pct=median(cv_pct),
  median_robust_range=median(p95-p05),
  median_full_range=median(maximum-minimum)),by=.(campaign,response,scale)]
fwrite(scale_comparison,file.path(out,"Table_02b_hour_day_week_campaign_comparison.csv"))

# Repeated-location Friedman tests ask whether the descriptive statistic changes
# systematically with aggregation scale. Campaign values are one observation
# per location and are directly comparable with the location summaries at the
# other three scales.
friedman_results <- rbindlist(lapply(c("sd","mad","cv_pct"),function(stat){
  descriptive[is.finite(get(stat)),{
    z<-dcast(.SD,location_n~scale,value.var=stat)
    z<-z[complete.cases(z)]
    if(nrow(z)<3) list(statistic=NA_real_,df=NA_real_,p_value=NA_real_,locations=nrow(z))
    else {ft<-friedman.test(as.matrix(z[,-1]));list(statistic=unname(ft$statistic),
      df=unname(ft$parameter),p_value=ft$p.value,locations=nrow(z))}
  },by=.(campaign,response)][,statistic_name:=stat]
}))
fwrite(friedman_results,file.path(out,"Table_02c_scale_Friedman_tests.csv"))

###############################################################################
##### BLOCK RELATIVE ERRORS AND CONVERGENCE
###############################################################################

make_block <- function(hours){
  z <- copy(x); z[,time_block:=floor_time(DATE.TIME,hours)]
  b <- z[,lapply(.SD,mean,na.rm=TRUE),by=.(campaign,time_block,location_n,height,
       horizontal_position,date=as.Date(time_block)),.SDcols=responses]
  bl <- melt(b,id.vars=c("campaign","time_block","date","location_n","height",
    "horizontal_position"),measure.vars=responses,variable.name="response",
    value.name="value")
  bl <- bl[is.finite(value)]
  bl[,barn_mean:=mean(value),by=.(campaign,time_block,response)]
  bl[,`:=`(RE_pct=100*(value-barn_mean)/barn_mean,
           abs_RE_pct=abs(100*(value-barn_mean)/barn_mean),
           location_rank=frank(value,ties.method="average")),
     by=.(campaign,time_block,response)]
  bl[,hours:=hours]
  bl
}
hourly <- make_block(1)
daily <- make_block(24)
convergence <- rbindlist(lapply(c(1,2,4,8,12,24,72,168),make_block))
fwrite(hourly,file.path(out,"campaign_hourly_location_relative_errors.csv"))
fwrite(daily,file.path(out,"campaign_daily_location_relative_errors.csv"))
conv_summary <- convergence[,.(blocks=uniqueN(time_block),
  median_abs_RE_pct=median(abs_RE_pct),p95_abs_RE_pct=quantile(abs_RE_pct,.95),
  spatial_RE_sd=sd(RE_pct)),by=.(campaign,response,hours)]
fwrite(conv_summary,file.path(out,"Table_03_temporal_convergence.csv"))

# Power-law convergence: median absolute RE = a * duration^b.
decay <- conv_summary[median_abs_RE_pct>0,{
  fit <- lm(log(median_abs_RE_pct)~log(hours))
  .(intercept=coef(fit)[1],exponent=coef(fit)[2],R2=summary(fit)$r.squared,
    predicted_persistent_direction=ifelse(coef(fit)[2]<0,"converging","not converging"))
},by=.(campaign,response)]
fwrite(decay,file.path(out,"Table_04_power_law_convergence.csv"))

###############################################################################
##### RANK STABILITY
###############################################################################

rank_results <- daily[,{
  wide <- dcast(.SD,time_block~location_n,value.var="location_rank")
  m <- as.matrix(wide[,-1]); m <- m[complete.cases(m),,drop=FALSE]
  if(nrow(m)<2) list(days=nrow(m),median_pairwise_spearman=NA_real_,kendall_W=NA_real_)
  else {
    cors <- cor(t(m),method="spearman",use="pairwise.complete.obs")
    R <- colSums(m); n <- ncol(m); k <- nrow(m)
    W <- 12*sum((R-mean(R))^2)/(k^2*(n^3-n))
    list(days=k,median_pairwise_spearman=median(cors[upper.tri(cors)],na.rm=TRUE),
         kendall_W=W)
  }
},by=.(campaign,response)]
fwrite(rank_results,file.path(out,"Table_05_daily_rank_stability.csv"))

rank_location <- daily[,.(median_rank=median(location_rank),rank_sd=sd(location_rank),
  probability_top_decile=mean(location_rank>=.9*max(location_rank),na.rm=TRUE),
  probability_bottom_decile=mean(location_rank<=.1*max(location_rank),na.rm=TRUE)),
  by=.(campaign,response,location_n)]
fwrite(rank_location,file.path(out,"Table_06_location_rank_persistence.csv"))

# Whole-campaign location ranks and deviations from the campaign spatial mean.
campaign_location <- x[,lapply(.SD,mean,na.rm=TRUE),
  by=.(campaign,location_n,height,horizontal_position),.SDcols=responses]
campaign_rank <- melt(campaign_location,id.vars=c("campaign","location_n","height",
  "horizontal_position"),measure.vars=responses,variable.name="response",
  value.name="value")[is.finite(value)]
campaign_rank[,campaign_barn_mean:=mean(value),by=.(campaign,response)]
campaign_rank[,`:=`(
  campaign_RE_pct=100*(value-campaign_barn_mean)/campaign_barn_mean,
  campaign_rank=frank(value,ties.method="average")),by=.(campaign,response)]
fwrite(campaign_rank,file.path(out,"Table_06b_whole_campaign_location_ranks.csv"))

###############################################################################
##### MIXED-EFFECTS VARIANCE DECOMPOSITION
###############################################################################

mixed_input <- hourly[value>0]
mixed_input[,hour_num:=as.integer(format(time_block,"%H"))]
mixed_tables <- list()
for(resp in responses){
  d <- mixed_input[response==resp]
  d[,`:=`(log_value=log(value),location_f=factor(location_n),date_f=factor(date))]
  # An additive campaign term is used because Campaign 2 has no middle-height
  # observations; a full campaign-by-three-height interaction is not estimable.
  fit <- lmer(log_value~campaign+height+factor(hour_num)+
                (1|location_f)+(1|date_f),data=d,REML=TRUE,
              control=lmerControl(optimizer="bobyqa"))
  saveRDS(fit,file.path(models,paste0("mixed_",resp,".rds")))
  vc <- as.data.table(as.data.frame(VarCorr(fit)))
  total <- sum(vc$vcov)
  vc[,`:=`(response=resp,variance_fraction=vcov/total)]
  mixed_tables[[resp]] <- vc
}
variance_components <- rbindlist(mixed_tables,fill=TRUE)
fwrite(variance_components,file.path(out,"Table_07_mixed_model_variance_components.csv"))

###############################################################################
##### GAM TEMPORAL STRUCTURE
###############################################################################

gam_tests <- list()
for(camp in unique(mixed_input$campaign)) for(resp in responses){
  d <- mixed_input[campaign==camp & response==resp]
  d[,`:=`(log_value=log(value),location_f=factor(location_n),
           date_num=as.numeric(date-min(date))+1,
           hour=as.numeric(format(time_block,"%H")))]
  fit <- bam(log_value~height+s(hour,bs="cc",k=12)+s(date_num,k=15)+
    s(location_f,bs="re"),data=d,method="fREML",discrete=TRUE,
    knots=list(hour=c(0,24)))
  saveRDS(fit,file.path(models,paste0("GAM_",gsub(" ","_",camp),"_",resp,".rds")))
  st <- as.data.table(summary(fit)$s.table,keep.rownames="smooth")
  st[,`:=`(response=resp,campaign=camp,converged=fit$converged)]
  gam_tests[[paste(camp,resp)]] <- st
}
fwrite(rbindlist(gam_tests,fill=TRUE),file.path(out,"Table_08_GAM_smooth_tests.csv"))

###############################################################################
##### PCA AND STABLE CLUSTERING
###############################################################################

pca_loadings <- list(); cluster_rows <- list(); silhouette_rows <- list()
for(camp in unique(hourly$campaign)) for(resp in responses){
  z <- hourly[campaign==camp & response==resp]
  w <- dcast(z,time_block~location_n,value.var="RE_pct")
  m <- as.matrix(w[,-1]); m <- m[,colSums(is.finite(m))>.8*nrow(m),drop=FALSE]
  for(j in seq_len(ncol(m))) m[!is.finite(m[,j]),j] <- median(m[,j],na.rm=TRUE)
  pc <- prcomp(m,center=TRUE,scale.=TRUE)
  load <- as.data.table(pc$rotation[,1:min(3,ncol(pc$rotation)),drop=FALSE],
                        keep.rownames="location_n")
  load[,`:=`(campaign=camp,response=resp,PC1_variance=summary(pc)$importance[2,1])]
  pca_loadings[[paste(camp,resp)]] <- load
  dist_loc <- as.dist(1-cor(m,use="pairwise.complete.obs"))
  hc <- hclust(dist_loc,method="average")
  sil <- rbindlist(lapply(2:min(8,ncol(m)-1),function(k){cl<-cutree(hc,k);
    data.table(k=k,mean_silhouette=mean(silhouette(cl,dist_loc)[,"sil_width"]))}))
  best_k <- sil[which.max(mean_silhouette),k]
  cluster_rows[[paste(camp,resp)]] <- data.table(campaign=camp,response=resp,
    location_n=as.integer(names(cutree(hc,best_k))),cluster=cutree(hc,best_k),best_k=best_k)
  sil[,`:=`(campaign=camp,response=resp)];silhouette_rows[[paste(camp,resp)]]<-sil
}
fwrite(rbindlist(pca_loadings,fill=TRUE),file.path(out,"Table_09_PCA_location_loadings.csv"))
fwrite(rbindlist(cluster_rows),file.path(out,"Table_10_location_clusters.csv"))
fwrite(rbindlist(silhouette_rows),file.path(out,"Table_11_cluster_silhouette.csv"))

###############################################################################
##### ENTROPY, MUTUAL INFORMATION AND REFERENCE-LOCATION SELECTION
###############################################################################

entropy4 <- function(v){v<-v[is.finite(v)];if(length(v)<8||diff(range(v))==0)return(NA_real_)
  p<-table(cut(v,seq(min(v),max(v),length.out=5),include.lowest=TRUE))/length(v)
  -sum(p[p>0]*log2(p[p>0]))/2}
mutual_info4 <- function(a,b){ok<-is.finite(a)&is.finite(b);a<-a[ok];b<-b[ok]
  if(length(a)<20||diff(range(a))==0||diff(range(b))==0)return(NA_real_)
  ca<-cut(a,unique(seq(min(a),max(a),length.out=5)),include.lowest=TRUE)
  cb<-cut(b,unique(seq(min(b),max(b),length.out=5)),include.lowest=TRUE)
  tab<-table(ca,cb);pab<-tab/sum(tab);pa<-rowSums(pab);pb<-colSums(pab)
  sum(pab[pab>0]*log2(pab[pab>0]/as.vector(outer(pa,pb))[pab>0]))}

reference_metrics <- hourly[,.(coverage=.N,bias_pct=mean(RE_pct),
  abs_bias_pct=abs(mean(RE_pct)),MAE_pct=mean(abs_RE_pct),RMSE_pct=sqrt(mean(RE_pct^2)),
  RE_sd=sd(RE_pct),spearman=cor(value,barn_mean,method="spearman"),
  entropy=entropy4(value),mutual_information=mutual_info4(value,barn_mean)),
  by=.(campaign,response,location_n,height,horizontal_position)]

common_locations <- Reduce(intersect,lapply(split(x$location_n,x$campaign),unique))
reference_metrics[,`:=`(
  common_location=location_n%in%common_locations,
  confirmed_problem=fcase(location_n==19L,"fan proximity",location_n==40L,"line leakage",default="none"),
  conservative_fan_corridor=horizontal_position%in%c(7L,15L))]

# Percentile ranks are oriented so that larger is always better.
reference_metrics[,`:=`(
  s_bias=1-frank(abs_bias_pct,ties.method="average")/.N,
  s_mae=1-frank(MAE_pct,ties.method="average")/.N,
  s_rmse=1-frank(RMSE_pct,ties.method="average")/.N,
  s_stability=1-frank(RE_sd,ties.method="average")/.N,
  s_correlation=frank(spearman,ties.method="average")/.N,
  s_information=frank(mutual_information,ties.method="average")/.N),
  by=.(campaign,response)]
reference_metrics[,score_component:=rowMeans(.SD,na.rm=TRUE),
  .SDcols=c("s_bias","s_mae","s_rmse","s_stability","s_correlation","s_information")]
reference_location <- reference_metrics[common_location & confirmed_problem=="none",
  .(statistical_score=mean(score_component[response%in%gases]),
    all_response_sensitivity_score=mean(score_component),
    worst_gas_score=min(score_component[response%in%gases]),
    mean_MAE_pct=mean(MAE_pct[response%in%gases]),
    max_MAE_pct=max(MAE_pct[response%in%gases]),
    mean_spearman=mean(spearman[response%in%gases]),
    mean_mutual_information=mean(mutual_information[response%in%gases]),
    conservative_fan_corridor=first(conservative_fan_corridor),height=first(height),
    horizontal_position=first(horizontal_position)),by=location_n]
reference_location[,rank_all_eligible:=frank(-statistical_score,ties.method="min")]
reference_location[,rank_conservative:=fifelse(!conservative_fan_corridor,
  frank(-statistical_score,ties.method="min",na.last="keep"),NA_integer_)]
reference_location[,sort_key:=fifelse(is.na(rank_conservative),9999L,rank_conservative)]
setorder(reference_location,sort_key,rank_all_eligible)
fwrite(reference_metrics,file.path(out,"Table_12_reference_location_metrics.csv"))
fwrite(reference_location,file.path(out,"Table_13_reference_location_ranking.csv"))

###############################################################################
##### DATE-GROUPED PREDICTIVE BENCHMARKS
###############################################################################

ml <- hourly[response%in%gases]
ml[,`:=`(date_group=as.character(as.Date(time_block)),hour_of_day=as.numeric(format(time_block,"%H")),
         location_factor=factor(location_n))]
cv_results <- list()
for(resp in gases){
  d <- ml[response==resp]
  dates <- sort(unique(d$date_group)); fold_map <- setNames(rep(1:5,length.out=length(dates)),dates)
  d[,fold:=fold_map[date_group]]
  for(k in 1:5){
    train<-d[fold!=k];test<-d[fold==k]
    linear<-lm(RE_pct~campaign+height+horizontal_position+hour_of_day+location_factor,data=train)
    gamfit<-bam(RE_pct~campaign+height+s(horizontal_position,k=10)+
      s(hour_of_day,bs="cc",k=10)+s(location_factor,bs="re"),data=train,
      discrete=TRUE,knots=list(hour_of_day=c(0,24)))
    tree<-rpart(RE_pct~campaign+height+horizontal_position+hour_of_day+location_n,
                data=train,control=rpart.control(cp=.002,minbucket=30))
    preds<-list(linear=predict(linear,test),GAM=predict(gamfit,test),tree=predict(tree,test))
    cv_results[[paste(resp,k)]]<-rbindlist(lapply(names(preds),function(model){p<-preds[[model]];
      data.table(response=resp,fold=k,model=model,n=sum(is.finite(p)),
        RMSE=sqrt(mean((test$RE_pct-p)^2,na.rm=TRUE)),
        MAE=mean(abs(test$RE_pct-p),na.rm=TRUE),
        R2=1-sum((test$RE_pct-p)^2,na.rm=TRUE)/sum((test$RE_pct-mean(test$RE_pct))^2,na.rm=TRUE))}))
  }
}
cv <- rbindlist(cv_results)
fwrite(cv,file.path(out,"Table_14_date_grouped_predictive_CV.csv"))
fwrite(cv[,.(RMSE=mean(RMSE),MAE=mean(MAE),R2=mean(R2)),by=.(response,model)],
       file.path(out,"Table_15_predictive_model_summary.csv"))

###############################################################################
##### EXPLORATORY CHANGE-POINT SCREENING
###############################################################################

daily_state <- daily[,.(barn_mean=first(barn_mean),spatial_cv=100*sd(value)/abs(mean(value))),
  by=.(campaign,response,time_block)]
change_points <- rbindlist(lapply(split(daily_state,by=c("campaign","response")),function(d){
  if(nrow(d)<10)return(NULL)
  d[,time_num:=as.numeric(time_block)]
  rbindlist(lapply(c("barn_mean","spatial_cv"),function(y){
    fit<-rpart(as.formula(paste(y,"~ time_num")),data=d,
               control=rpart.control(cp=.01,minbucket=5))
    splits<-if(is.null(fit$splits))numeric() else sort(unique(fit$splits[,"index"]))
    data.table(campaign=first(d$campaign),response=first(d$response),outcome=y,
      split_time=as.POSIXct(splits,origin="1970-01-01",tz=tz))
  }))
}),fill=TRUE)
fwrite(daily_state,file.path(out,"Table_16_daily_barn_state_for_change_points.csv"))
fwrite(change_points,file.path(out,"Table_17_exploratory_change_points.csv"),
       dateTimeAs="write.csv")

###############################################################################
##### FIGURES
###############################################################################

pconv <- ggplot(conv_summary,aes(hours,median_abs_RE_pct,colour=response))+
  geom_line()+geom_point()+scale_x_log10(breaks=c(1,2,4,8,12,24,72,168))+
  facet_wrap(~campaign)+theme_bw()+labs(x="Averaging duration (h)",
  y="Median absolute location error (%)",colour="Response")
ggsave(file.path(fig,"Fig_01_temporal_convergence.png"),pconv,width=11,height=6,dpi=300)

prank <- ggplot(daily[response%in%gases],aes(time_block,factor(location_n),fill=location_rank))+
  geom_tile()+facet_grid(response~campaign,scales="free_x",space="free_x")+
  scale_fill_viridis_c()+theme_bw()+labs(x="Date",y="Location",fill="Daily rank")+
  theme(axis.text.y=element_text(size=5))
ggsave(file.path(fig,"Fig_02_daily_location_rank_heatmap.png"),prank,width=14,height=10,dpi=300)

pref <- reference_location[rank_all_eligible<=15]
pref_long <- melt(pref,id.vars=c("location_n","conservative_fan_corridor"),
  measure.vars=c("statistical_score","worst_gas_score"),variable.name="score",value.name="value")
ppref <- ggplot(pref_long,aes(reorder(factor(location_n),value),value,fill=score))+
  geom_col(position="dodge")+coord_flip()+theme_bw()+labs(x="Candidate location",y="Score",
  fill=NULL,title="Cross-campaign reference-location candidates")
ggsave(file.path(fig,"Fig_03_reference_location_ranking.png"),ppref,width=8,height=7,dpi=300)

pcv <- ggplot(cv,aes(model,RMSE,fill=model))+geom_boxplot()+facet_wrap(~response,scales="free_y")+
  theme_bw()+theme(legend.position="none")+labs(x=NULL,y="Date-grouped CV RMSE (percentage points)")
ggsave(file.path(fig,"Fig_04_predictive_model_validation.png"),pcv,width=9,height=5,dpi=300)

pscale <- ggplot(scale_comparison,
  aes(scale,median_location_CV_pct,group=response,colour=response))+
  geom_line()+geom_point(size=2)+facet_wrap(~campaign)+theme_bw()+
  scale_x_discrete(drop=FALSE)+labs(x="Aggregation level",
    y="Median location CV (%)",colour="Response")
ggsave(file.path(fig,"Fig_04b_hour_day_week_campaign_CV.png"),pscale,
       width=11,height=6,dpi=300)

# Engineering check: statistical selection and conservative fan corridors on
# the supplied floor plan. Pixel coordinates come from the existing mapping
# workflow and are not re-estimated here.
if(requireNamespace("png",quietly=TRUE)){
  floor_file <- file.path(workflow,
    "Manuscript_3_Mapping_high_resolution_concentration_Latex_draft",
    "figures","Fig1a_campaign1_sampling_layout_no_ringline.png")
  coord_file <- file.path(workflow,"clean_data","campaign_cv_spatial",
                          "location_coordinate_lookup.csv")
  if(file.exists(floor_file)&&file.exists(coord_file)){
    img <- png::readPNG(floor_file); h <- dim(img)[1]; w <- dim(img)[2]
    coords <- fread(coord_file)
    coords[,plot_y:=h-map_y_pixel]
    coords[,status:=fcase(location_number==25L,"Selected candidate",
      horizontal_position%in%c(7L,15L),"Conservative fan corridor",
      location_number==40L,"Confirmed leakage",default="Other location")]
    pmap <- ggplot()+annotation_raster(img,0,w,0,h)+
      geom_point(data=coords[status!="Other location"],
        aes(map_x_pixel,plot_y,shape=status,colour=status),size=4,stroke=1.3)+
      geom_label(data=coords[location_number==25L],aes(map_x_pixel,plot_y,label="25"),
        nudge_y=55,size=4,fontface="bold")+
      scale_colour_manual(values=c("Selected candidate"="blue",
        "Conservative fan corridor"="red","Confirmed leakage"="black"))+
      scale_shape_manual(values=c("Selected candidate"=8,
        "Conservative fan corridor"=4,"Confirmed leakage"=4))+
      coord_fixed(xlim=c(0,w),ylim=c(0,h),expand=FALSE)+theme_void()+
      theme(legend.position="bottom")+labs(colour=NULL,shape=NULL)
    ggsave(file.path(fig,"Fig_05_reference_location_floorplan.png"),pmap,
           width=13,height=7,dpi=300,bg="white")
  }
}

writeLines(c(
  "MANUSCRIPT 3 MULTISCALE AND REFERENCE-LOCATION ANALYSIS V04",
  "==========================================================",
  "Campaign 1: recovered June-August four-analyser data.",
  "Campaign 2: reduced 32-location CRDS design.",
  "All validation folds are grouped by date.",
  "Reference screening excludes confirmed location 19 fan influence and location 40 leakage.",
  "Primary pragmatic ranking additionally excludes floor-plan fan corridors 19-21 and 43-45.",
  "The final engineering choice requires author confirmation of fan axes.",
  "Regression-tree change points are exploratory screening, not formal causal break tests.",
  paste("R version:",R.version.string)
),file.path(out,"analysis_readme.txt"))

capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
