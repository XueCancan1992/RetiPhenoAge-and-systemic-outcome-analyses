

### 1. PheWAS Association
# ============================================================
ds<-readRDS("your_dataset.rds")

phecodes_file<-read.csv("your_phecode_map.csv")#
head(phecodes_file)
phecodes<-phecodes_file$description_clean



library(survival)

####*******************************
# script about phecode grouping

{
  phecode.grouping <- read.csv("your_phecode_map.csv")
  phecode = phecode.grouping$description_clean
  phecodes[!phecodes %in% phecode]
  phecode.grouping$phecodename = phecode
}


cox.result <- NULL
covariates <- c("RetiPhenoAge_sd","Age", "Sex.0","Ethnicity3")
#covariates <- c("RetiPhenoAge_quartile_sd","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0")
#covariates <- c("RetiPhenoAge_quartile_sd","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0","PhenoAge")


#covariates <- c("RetiPhenoAge_quartile","Age", "Sex.0","Ethnicity3")
#covariates <- c("RetiPhenoAge_quartile","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0")
#covariates <- c("RetiPhenoAge_quartile","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0","PhenoAge")


#covariates <- c("RetiPhenoAge_quartile_num","Age", "Sex.0","Ethnicity3")
#covariates <- c("RetiPhenoAge_quartile_num","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0")
#covariates <- c("RetiPhenoAge_quartile_num","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0","PhenoAge")




for (phegroup in unique(phecode.grouping$group)) {
  phegroup.phecodes = phecode.grouping$phecodename[phecode.grouping$group %in% phegroup]
  dt <- as.data.frame(ds)
  for (phegroup.phecode in phegroup.phecodes) {
    if (paste('all_inpatient_',phegroup.phecode,'_before',sep='') %in% colnames(ds)) {
      dt <- dt[!dt[,paste('all_inpatient_',phegroup.phecode,'_before',sep='')] %in% 1,]
    }
  }
  if (nrow(dt)>0) {
    for (phecode in phegroup.phecodes) {
      try({
        type = 'IMevent'
        fm <- as.formula(paste("Surv(all_IMevent_", phecode, "_year,", "all_IMevent_", phecode, ")~", paste(covariates, collapse = "+"),sep=''))
        m <- coxph(fm, data = dt)
        cox.result <- rbind(cox.result,
                            cbind(phecode,type,covariates,summary(m)$coef))
        rm(fm,m,type)
      },silent=TRUE)
    }
  }
  rm(dt)
}
dim(cox.result)
class(cox.result)

cox.result <- as.data.frame(cox.result)


for (c in 4:8) cox.result[,c] <- as.numeric(as.character(cox.result[,c]))
cox.result$HR.LCL <- exp(cox.result[,'coef'] - 1.96*cox.result[,'se(coef)'])
cox.result$HR.UCL <- exp(cox.result[,'coef'] + 1.96*cox.result[,'se(coef)'])
cox.result <- cox.result[,c(1:3,5,9:10,4,6,8)]
colnames(cox.result)[4:9] <- c('HR','LCL','UCL','Coef','SE','P')

head(cox.result)




library(dplyr)

adjust_fdr <- function(data) {
  data %>%
    group_by(group) %>%
    mutate(fdr_P = p.adjust(P_value, method = "fdr")) %>%
    ungroup() %>%
    as.data.frame()  # Ensure the output is a data.frame
}
results <- adjust_fdr(cox.results)

head(results)


write.csv(result,file='your results.csv')


# ============================================================


#### 2 C-index and cummulative incidence plot
# ============================================================
library(survival)
#library(survminer)
library(pROC)
library(epitools)

coxph.sum <- function(m,alpha=0.05) {
  require(survival)
  # summary: c('Coef','C.index','C.count','GOF')
  {
    # alpha <- 0.05 # 1- confidence level
    ms <- summary(m)
    # Coef
    {
      ms$Coef <- cbind(ms$coefficients,ms$conf.int)
      ms$Coef <- ms$Coef[,c(1,3,6:9,4:5)]
      if (!'matrix' %in% class(ms$Coef)) {
        ms$Coef <- as.data.frame(t(ms$Coef))
        rownames(ms$Coef) <- rownames(ms$coefficients)
      } else {
        ms$Coef <- as.data.frame(ms$Coef)
      }
      colnames(ms$Coef) <- c('Coef','SE','HR','HR-','LCL','UCL','Z','P')
    }
    # C.index
    {
      ms$C.index <- c(ms$concordance,ms$concordance[1]+qnorm(c(alpha/2,1-alpha/2))*ms$concordance[2])
      names(ms$C.index) <- c('C','SE','LCL','UCL')
      ms$C.index <- as.data.frame(t(ms$C.index))
      rownames(ms$C.index) <- 'Concordance'
    }
    # C.count
    {
      ms$C.count <- concordance(m)$count
    }
    # GOF: goodness of fit
    {
      ms$GOF <- as.data.frame(rbind(
        c(ms$loglik,ms$logtest), 
        c(NA,NA,ms$sctest), 
        c(NA,NA,ms$waldtest), 
        cbind(NA,NA,cox.zph(m)$table))) 
      colnames(ms$GOF) <- c('logLik0','logLik1','stat','df','P')
      rownames(ms$GOF) <- c('LRT','Score.test','Wald.test',paste('PH.Chisq.test.',rownames(cox.zph(m)$table),sep=''))
    }
  }
  return(ms)
}



coxph.predict <- function(m,X,time) {
  require(survival)
  cstat <- concordance(m)$concordance
  cstat.se <- sqrt(concordance(m)$var)
  cstat.LCL <- cstat-1.96*cstat.se
  cstat.UCL <- cstat+1.96*cstat.se
  cstat.table <- as.data.frame(cbind(cstat,cstat.se,cstat.LCL,cstat.UCL))
  colnames(cstat.table) <- c('C.stat','C.SE','C.LCL','C.UCL')
  rownames(cstat.table) <- 'C.statistic'
  
  cumhaz <- basehaz(m,centered=FALSE)
  cumhaz.t <- tail(cumhaz$haz[cumhaz$time<=time],1)
  beta <- as.matrix(m$coefficients,ncol=1)
  beta.var <- m$var
  coef.table <- as.data.frame(summary(m)$coef)
  if ('robust se' %in% colnames(coef.table)) {
    colnames(coef.table) <- c('Coef','OR','Coef.SE','Coef.robust.SE','zstat','p.value')
    coef.table$Coef.LCL <- coef.table[,1]-1.96*coef.table[,4]
    coef.table$Coef.UCL <- coef.table[,1]+1.96*coef.table[,4]
    coef.table$OR.LCL <- exp(coef.table$Coef.LCL)
    coef.table$OR.UCL <- exp(coef.table$Coef.UCL)
    coef.table <- coef.table[,c(2,9,10,1,4,7,8,5,6)]
  } else {
    colnames(coef.table) <- c('Coef','OR','Coef.SE','zstat','p.value')
    coef.table$Coef.LCL <- coef.table[,1]-1.96*coef.table[,3]
    coef.table$Coef.UCL <- coef.table[,1]+1.96*coef.table[,3]
    coef.table$OR.LCL <- exp(coef.table$Coef.LCL)
    coef.table$OR.UCL <- exp(coef.table$Coef.UCL)
    coef.table <- coef.table[,c(2,8,9,1,3,6,7,4,5)]
  }
  
  lp <- X %*% beta
  risk <- exp(lp)
  expected <- cumhaz.t*risk
  survival <- exp(-expected)
  
  out <- NULL
  out$model <- m
  out$cstat.table <- cstat.table
  out$cstat <- cstat
  out$cstat.se <- cstat.se
  out$cstat.LCL <- cstat.LCL
  out$cstat.UCL <- cstat.UCL
  out$cumhaz <- cumhaz
  out$cumhaz.t <- cumhaz.t
  out$beta <- beta
  out$beta.var <- beta.var
  out$coef.table <- coef.table
  out$lp <- lp
  out$risk <- risk
  out$expected <- expected
  out$survival <- survival
  return(out)
}

concordance.comp <- function(contrast=c(1,-1),...) {
  require(survival)
  cstat <- concordance(...)
  nc <- length(coef(cstat))
  contrast <- contrast[1:min(nc,length(contrast))]
  contrast <- c(contrast,rep(0,nc-length(contrast)))
  zstat <- (contrast %*% coef(cstat)) / (contrast %*% vcov(cstat) %*% contrast)
  pstat <- pnorm(-abs(zstat),0,1)*2
  cstat <- cbind(cstat$concordance,sqrt(diag(cstat$var)),
                 cstat$concordance-1.96*sqrt(diag(cstat$var)),
                 cstat$concordance+1.96*sqrt(diag(cstat$var)),
                 c(pstat,rep(NA,nc-1)))
  colnames(cstat) <- c('CStat','SE','LCL','UCL','P')
  rownames(cstat) <- paste('model',1:length(contrast))
  return(cstat)
}

roc.sum <- function(formula,data,threshold='youden') {
  require(pROC)
  require(epitools)
  roc <- roc(formula,data,ci=TRUE)
  auc <- roc$ci[c(2,1,3)]
  if (threshold=='youden') {
    roc$youden <- roc$sensitivities + roc$specificities - 1
    roc$prediction <- as.integer(cut(roc$predictor,c(-Inf,roc$thresholds[which.max(roc$youden)],Inf)))-1
  } else {
    roc$prediction <- as.integer(cut(roc$predictor,c(-Inf,threshold,Inf)))-1
  }
  spec <- epitools::binom.approx(sum((roc$response %in% 0)&(roc$prediction %in% 0)),sum(roc$response %in% 0))[3:5]
  sens <- epitools::binom.approx(sum((roc$response %in% 1)&(roc$prediction %in% 1)),sum(roc$response %in% 1))[3:5]
  npv <- epitools::binom.approx(sum((roc$response %in% 0)&(roc$prediction %in% 0)),sum(roc$prediction %in% 0))[3:5]
  ppv <- epitools::binom.approx(sum((roc$response %in% 1)&(roc$prediction %in% 1)),sum(roc$prediction %in% 1))[3:5]
  roc.sum <- as.data.frame(rbind(auc,spec,sens,npv,ppv))
  rownames(roc.sum) <- c('AUC','Spec','Sens','NPV','PPV')
  colnames(roc.sum) <- c('Est','LCL','UCL')
  return(roc.sum)
}

plot.surv.strata <- function(data,time,event,strata,fun=NULL,file,
                             legend.title = strata,
                             legend.position = c("top", "bottom", "left", "right", "none")[1],
                             xscale='none',yscale='none',
                             censor=FALSE,conf.int = FALSE,test.for.trend = FALSE,pval=FALSE,
                             risk.table = TRUE,risk.table.y.text = TRUE,risk.table.title = 'No. at risk',
                             tables.y.text.col = FALSE,cumevents = TRUE,cumcensor = FALSE,cumevents.y.text.col = FALSE,
                             xlab = 'Time (year)',ylab = 'Cumulative events (%)',ylim=NULL,xlim=NULL,break.y.by=NULL,break.x.by=NULL,fs=14,
                             ...) {
  if (tolower(tools::file_ext(file))=='svg') {
    svg(file,...)
  } else if (tolower(tools::file_ext(file))=='png') {
    png(file,...)
  }
  # ggsurvplot
  {
    eval(parse(
      text=paste("fit <- survfit(Surv(", time, ",", event, ")~",strata,",data=",data,")",sep='') ))
    g <- ggsurvplot(fit,data=eval(parse(text=data)),fun=fun,
                    # surv.scale=c('default','percent')[2],
                    censor=censor,
                    conf.int = conf.int,
                    # ylim=ylim,break.y.by=break.y.by,
                    # xlim=xlim,break.x.by=break.x.by,
                    font.x = c(fs,"plain", "black"),
                    font.y = c(fs,"plain", "black"),
                    font.tickslab = c(fs,"plain", "black"),
                    test.for.trend = test.for.trend,pval=pval,
                    # ggtheme = theme_bw(),
                    risk.table = risk.table,
                    cumevents = cumevents,
                    # cumevents.y.text.col = cumevents.y.text.col,
                    cumcensor = cumcensor,
                    risk.table.y.text = risk.table.y.text,
                    tables.y.text.col = tables.y.text.col,
                    risk.table.title = 'No. at risk',
                    tables.theme =
                      theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            legend.position="none",
                            panel.background=element_blank(),
                            panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            plot.background=element_blank()),
                    # palette=c('#0000AA','#AA0000','#00AA00','#DDAA00'),
                    legend.title = legend.title,
                    legend.labs = levels(eval(parse(text=paste(data,'$',strata,sep='')))),
                    legend = legend.position,
                    font.legend = c(fs,"plain", "black"),
                    xlab=xlab,ylab=ylab)
  }
  # g <- ggpar(g,xscale=xscale,yscale=yscale)
  print(g)
  dev.off()
}

# ============================================================

# ============================================================


ds<-readRDS("your_dataset.rds")

phecodes_file<-read.csv("your_phecode_map.csv")
head(phecodes_file)
phecodes<-phecodes_file$description_clean



# YOUR INPUT

timepoint <- 10
covariates0 <- c("Age", "Sex.0")
#covariates0 <- 'Age'
#covariates <- 'RetiPhenoAge'
covariates <- c('RetiPhenoAge',"Age","Sex.0")
#covariates <- c('RetiPhenoAge',"Age") 
strata <- NULL # for numerical retiphenoage (no plotting)
#strata <- 'RetiPhenoAge_quartile' # for categorical retiphenoage (plotting generated)
{
  legend.position = c("top", "bottom", "left", "right", "none")[2]
  ylim=c(.6,1)
  break.y.by=.1
  xlim=c(0,10)
  break.x.by=2
} 
output1 <- 'model.csv'
output2 <- 'AUC.csv'
output3 <- 'C_index.csv'



cox.result <- NULL
cox.c.index <- NULL
cox.auc <- NULL



{

  phecode.grouping <- read.csv("your_phecode_map.csv")
  phecode = phecode.grouping$description_clean
  
  # check is there any phecodes not in phecode
  phecodes[!phecodes %in% phecode]
  
  phecode.grouping$phecodename = phecode
}





for (phegroup in unique(phecode.grouping$group)) {
  phegroup.phecodes = phecode.grouping$phecodename[phecode.grouping$group %in% phegroup]
  dt <- as.data.frame(ds)
  dt <- dt[rowSums(is.na(dt[,covariates]))==0,]
  for (phegroup.phecode in phegroup.phecodes) {
    if (paste('all_inpatient_',phegroup.phecode,'_before',sep='') %in% colnames(ds)) {
      dt <- dt[!dt[,paste('all_inpatient_',phegroup.phecode,'_before',sep='')] %in% 1,]
    }
  }
  if (nrow(dt)>0) {
    DT <- dt
    for (phecode in phegroup.phecodes) {
      try({
        type = 'IMevent'
        timevar <- paste('all_IMevent_',phecode,'_year',sep='')
        eventvar <- paste('all_IMevent_',phecode,sep='')
        dt <- DT
        dt <- dt[!is.na(dt[,paste('all_IMevent_',phecode,sep='')]),]
        dt <- dt[!is.na(dt[,paste('all_IMevent_',phecode,'_year',sep='')]),]
        dt[,paste('all_IMevent_',phecode,sep='')] <- as.integer(factor(dt[,paste('all_IMevent_',phecode,sep='')],0:1))
        fm <- as.formula(paste("Surv(all_IMevent_", phecode, "_year,", "all_IMevent_", phecode, ")~", paste(covariates, collapse = "+"),sep=''))
        fm <- (paste("Surv(all_IMevent_", phecode, "_year,", "all_IMevent_", phecode, ")~", paste(covariates, collapse = "+"),sep=''))
        fm0 <- as.formula(paste("Surv(all_IMevent_", phecode, "_year,", "all_IMevent_", phecode, ")~", paste(covariates0, collapse = "+"),sep=''))
        # coxph
        {
          m <- coxph(as.formula(fm), data = dt)
          m0 <- coxph(as.formula(fm0), data = dt)
          ms <- coxph.sum(m)
          ms0 <- coxph.sum(m0)
          LRT <- pchisq((ms$GOF[1,2] - ms0$GOF[1,2])*2,ms$GOF[1,4]-ms0$GOF[1,4],lower.tail=FALSE)
        }
        # prediction
        {
          dt$outcome_at_timepoint <- ifelse(
            (dt[,paste("all_IMevent_", phecode,sep='')] %in% 0)&
              (dt[,paste("all_IMevent_", phecode,'_year',sep='')] >= timepoint),0,
            ifelse(
              (dt[,paste("all_IMevent_", phecode,sep='')] %in% 1)&
                (dt[,paste("all_IMevent_", phecode,'_year',sep='')] <= timepoint),1,
              ifelse(
                (dt[,paste("all_IMevent_", phecode,sep='')] %in% 0)&
                  (dt[,paste("all_IMevent_", phecode,'_year',sep='')] < timepoint),NA,
                0)))
          mp <- coxph.predict(m,model.matrix(as.formula(fm),dt)[,-1],timepoint)
          dt$risk_at_timepoint <- 1-mp$survival
          roc.result <- roc.sum(outcome_at_timepoint~risk_at_timepoint,dt)
        }
        # summary
        {
          cox.result <- rbind(cox.result,
                              cbind(phecode,type,rownames(ms$Coef),ms$Coef))
          cox.c.index <- rbind(cox.c.index,
                               cbind(phecode,type,
                                     t(concordance.comp(contrast=c(-1,1),m,m0)[1,]),
                                     t(concordance.comp(contrast=c(-1,1),m,m0)[2,]),
                                     LRT))
          cox.auc <- rbind(cox.auc,
                           cbind(phecode,type,rownames(roc.result),roc.result))
        }
        if (!is.null(strata)) {
          plot.surv.strata('dt',timevar,eventvar,strata,fun=NULL,
                           file=paste('KMcurve_',type,'_',phecode,'_by_',strata,'.png',sep=''),
                           legend.title = strata,xscale='none',yscale='none',
                           legend.position = legend.position,
                           censor=FALSE,conf.int = TRUE,test.for.trend = FALSE,pval=FALSE,
                           risk.table = TRUE,risk.table.y.text = TRUE,risk.table.title = 'No. at risk',
                           tables.y.text.col = TRUE,cumevents = FALSE,cumcensor = FALSE,cumevents.y.text.col = FALSE,
                           xlab = 'Time (year)',ylab = 'Cumulative events (%)',ylim=ylim,xlim=xlim,break.y.by=break.y.by,break.x.by=break.x.by,fs=14)
        }
        rm(dt,type,m,ms,fm0,m0,ms0,LRT,mp,roc.result,timevar,eventvar)
      },silent=TRUE)
    }
  }
  rm(dt)
}

# formatting results
{
  colnames(cox.result)[3] <- 'covariate'
  cox.result <- as.data.frame(cox.result)
  for (c in 4:ncol(cox.result)) cox.result[,c] <- as.numeric(as.character(cox.result[,c]))
  
  colnames(cox.auc)[3] <- 'metric'
  cox.auc <- as.data.frame(cox.auc)
  for (c in 4:ncol(cox.auc)) cox.auc[,c] <- as.numeric(as.character(cox.auc[,c]))
  
  colnames(cox.c.index)[8:12] <- paste(colnames(cox.c.index)[8:12],'_H0',sep='')
  cox.c.index <- cox.c.index[,-12]
  cox.c.index <- as.data.frame(cox.c.index)
  for (c in 3:ncol(cox.c.index)) cox.c.index[,c] <- as.numeric(as.character(cox.c.index[,c]))
}


head(cox.result)
head(cox.auc)
head(cox.c.index)

write.csv(cox.result,file=output1)
write.csv(cox.auc,file=output2)
write.csv(cox.c.index,file=output3)


# ============================================================



### 3 calibration analyses

# ============================================================
# Calibration statistics for all outcomes
# ============================================================
df<-ds
library(survival)
library(dplyr)
library(purrr)

# Load outcome CSV
outcome_df <- read.csv("your_outcome_list.csv")

# Change Outcome to your actual column name
outcomes <- outcome_df$Outcome
outcomes <- outcomes[!is.na(outcomes) & outcomes != ""]

# Settings
time_horizon <- 10


covariates <- c(
  "RetiPhenoAge_sd","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0"
)

dir.create("calibration_results", showWarnings = FALSE)

run_calibration_stats <- function(outcome_name, df, time_horizon, covariates) {
  
  event_var <- paste0("all_IMevent_", outcome_name)
  time_var  <- paste0("all_IMevent_", outcome_name, "_year")
  
  needed_vars <- c(event_var, time_var, covariates)
  
  if (!all(needed_vars %in% names(df))) {
    return(data.frame(
      outcome = outcome_name,
      n = NA,
      events = NA,
      calibration_slope = NA,
      calibration_intercept = NA,
      brier_score = NA,
      status = "missing variables"
    ))
  }
  
  dat <- df[, needed_vars]
  dat <- na.omit(dat)
  
  names(dat)[names(dat) == event_var] <- "event"
  names(dat)[names(dat) == time_var]  <- "time"
  
  dat <- dat %>% filter(time > 0)
  
  if (nrow(dat) < 100 || sum(dat$event == 1) < 20) {
    return(data.frame(
      outcome = outcome_name,
      n = nrow(dat),
      events = sum(dat$event == 1),
      calibration_slope = NA,
      calibration_intercept = NA,
      brier_score = NA,
      status = "too few events"
    ))
  }
  
  formula_txt <- paste(
    "Surv(time, event) ~",
    paste(covariates, collapse = " + ")
  )
  
  fit <- coxph(as.formula(formula_txt), data = dat, x = TRUE)
  
  lp <- predict(fit, type = "lp")
  
  # Calibration slope
  slope_fit <- coxph(Surv(time, event) ~ lp, data = dat)
  calibration_slope <- coef(slope_fit)[1]
  
  # Predicted risk at time_horizon
  base_surv <- survfit(fit)
  s0_t <- summary(base_surv, times = time_horizon, extend = TRUE)$surv
  
  pred_risk <- 1 - s0_t ^ exp(lp)
  pred_risk <- pmin(pmax(pred_risk, 1e-6), 1 - 1e-6)
  
  # Observed event by time_horizon
  observed_event <- ifelse(dat$time <= time_horizon & dat$event == 1, 1, 0)
  
  # Brier score
  brier_score <- mean((observed_event - pred_risk)^2, na.rm = TRUE)
  
  # Approximate calibration intercept
  intercept_fit <- glm(
    observed_event ~ offset(qlogis(pred_risk)),
    family = binomial
  )
  
  calibration_intercept <- coef(intercept_fit)[1]
  
  data.frame(
    outcome = outcome_name,
    n = nrow(dat),
    events = sum(dat$event == 1),
    calibration_slope = calibration_slope,
    calibration_intercept = calibration_intercept,
    brier_score = brier_score,
    status = "ok"
  )
}

calibration_stats <- map_dfr(
  outcomes,
  run_calibration_stats,
  df = df,
  time_horizon = time_horizon,
  covariates = covariates
)

write.csv(
  calibration_stats,
  "calibration_results/calibration_statistics_all_outcomes.csv",
  row.names = FALSE
)

getwd()




### 4 calibration plots

# ============================================================
# Calibration plots for selected key outcomes
# ============================================================

library(survival)
library(rms)
library(dplyr)
library(ggplot2)

# Key outcomes only



plot_event_calibration <- function(
    outcome_name,
    df,
    time_horizon = 10,
    covariates = c("RetiPhenoAge_sd","Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0")
) {
  
  event_var <- paste0("all_IMevent_", outcome_name)
  time_var  <- paste0("all_IMevent_", outcome_name, "_year")
  
  needed_vars <- c(event_var, time_var, covariates)
  
  dat <- df[, needed_vars]
  dat <- na.omit(dat)
  
  names(dat)[names(dat) == event_var] <- "event"
  names(dat)[names(dat) == time_var]  <- "time"
  
  dat <- dat %>% filter(time > 0)
  
  formula_txt <- paste(
    "Surv(time, event) ~",
    paste(covariates, collapse = " + ")
  )
  
  fit <- coxph(as.formula(formula_txt), data = dat)
  
  # baseline survival
  basefit <- survfit(fit)
  
  s0_t <- summary(basefit, times = time_horizon, extend = TRUE)$surv
  
  # linear predictor
  lp <- predict(fit, type = "lp")
  
  # predicted event probability
  dat$pred_risk <- 1 - (s0_t ^ exp(lp))
  
  # group into deciles
  dat$decile <- ntile(dat$pred_risk, 10)
  
  # observed event probability via KM
  calibration_df <- dat %>%
    group_by(decile) %>%
    summarise(
      mean_predicted = mean(pred_risk),
      observed = 1 - summary(
        survfit(Surv(time, event) ~ 1, data = cur_data()),
        times = time_horizon,
        extend = TRUE
      )$surv
    )
  
  p <- ggplot(calibration_df,
              aes(x = mean_predicted,
                  y = observed)) +
    geom_point(size = 3) +
    geom_line() +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    labs(
      x = paste0("Predicted ", time_horizon, "-year event probability"),
      y = paste0("Observed ", time_horizon, "-year event probability"),
      title = paste0("Calibration: ", outcome_name)
    ) +
    theme_bw()
  
  print(p)
  
  return(calibration_df)
}





# ============================================================
# Smoothed calibration plots: predicted vs observed event risk
# Cox model, grouped KM observed risk + LOESS smooth
# ============================================================

library(survival)
library(dplyr)
library(ggplot2)



key_outcomes <- c(
  "Stroke",
  "Coronary_atherosclerosis_Atherosclerotic_heart_disease",
  "Hypertension",
  "Ischemic_heart_disease",
  "Hyperlipidemia",
  "Type_2_diabetes",
  "Pneumonia",
  "Chronic_obstructive_pulmonary_disease_COPD",
  "Pleural_effusion",
  "Anemia",
  "Acute_kidney_failure",
  "Chronic_kidney_disease",
  "Dementias_and_cerebral_degeneration",
  "Osteoporosis",
  "Malignant_neoplasm_of_the_digestive_organs"
)
time_horizon <- 10

covariates <- c(
  "RetiPhenoAge_sd"
,"Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0")
  
  


dir.create("calibration_plots", showWarnings = FALSE)

plot_event_calibration <- function(
    outcome_name,
    df,
    time_horizon = 10,
    covariates,
    n_groups = 10
) {
  
  event_var <- paste0("all_IMevent_", outcome_name)
  time_var  <- paste0("all_IMevent_", outcome_name, "_year")
  
  needed_vars <- c(event_var, time_var, covariates)
  
  missing_vars <- setdiff(needed_vars, names(df))
  
  if (length(missing_vars) > 0) {
    message(
      "Skipping ", outcome_name,
      ": missing variables: ",
      paste(missing_vars, collapse = ", ")
    )
    return(NULL)
  }
  
  dat <- df[, needed_vars]
  dat <- na.omit(dat)
  
  names(dat)[names(dat) == event_var] <- "event"
  names(dat)[names(dat) == time_var]  <- "time"
  
  dat <- dat %>%
    filter(time > 0)
  
  if (nrow(dat) < 100 || sum(dat$event == 1) < 20) {
    message("Skipping ", outcome_name, ": too few observations/events")
    return(NULL)
  }
  
  formula_txt <- paste(
    "Surv(time, event) ~",
    paste(covariates, collapse = " + ")
  )
  
  fit <- coxph(
    as.formula(formula_txt),
    data = dat,
    x = TRUE
  )
  
  base_surv <- survfit(fit)
  s0_t <- summary(base_surv, times = time_horizon, extend = TRUE)$surv
  
  lp <- predict(fit, type = "lp")
  
  dat$pred_risk <- 1 - (s0_t ^ exp(lp))
  
  dat <- dat %>%
    mutate(risk_group = ntile(pred_risk, n_groups))
  
  calibration_df <- dat %>%
    group_by(risk_group) %>%
    group_modify(~ {
      km <- survfit(Surv(time, event) ~ 1, data = .x)
      km_sum <- summary(km, times = time_horizon, extend = TRUE)
      
      tibble(
        n = nrow(.x),
        events = sum(.x$event == 1),
        mean_predicted = mean(.x$pred_risk, na.rm = TRUE),
        observed = 1 - km_sum$surv,
        lower = 1 - km_sum$upper,
        upper = 1 - km_sum$lower
      )
    }) %>%
    ungroup()
  
  max_axis <- max(
    calibration_df$mean_predicted,
    calibration_df$observed,
    calibration_df$upper,
    na.rm = TRUE
  )
  
  min_axis <- min(
    calibration_df$mean_predicted,
    calibration_df$observed,
    calibration_df$lower,
    na.rm = TRUE
  )
  
  p <- ggplot(
    calibration_df,
    aes(x = mean_predicted, y = observed)
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.005,
      alpha = 0.6
    ) +
    geom_point(size = 2.5) +
    geom_smooth(
      method = "loess",
      se = FALSE,
      span = 1,
      linewidth = 0.8
    ) +
    coord_equal(
      xlim = c(min_axis, max_axis),
      ylim = c(min_axis, max_axis)
    ) +
    labs(
      title = paste0("Calibration: ", outcome_name),
      x = paste0("Predicted ", time_horizon, "-year event probability"),
      y = paste0("Observed ", time_horizon, "-year event probability")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = paste0(
      "calibration_plots/calibration_event_probability_",
      outcome_name,
      ".svg"
    ),
    plot = p,
    width = 6,
    height = 6
  )
  
  return(list(
    outcome = outcome_name,
    model = fit,
    calibration_data = calibration_df,
    plot = p
  ))
}

calibration_results <- lapply(
  key_outcomes,
  plot_event_calibration,
  df = df,
  time_horizon = time_horizon,
  covariates = covariates,
  n_groups = 10
)



# 5  Proportional hazards assumption check for all outcomes
# Schoenfeld residuals using cox.zph()
# ============================================================

library(survival)
library(dplyr)
library(purrr)

df<-readRDS("your_dataset.rds")
# Outcome list
outcome_df <- read.csv("your_outcome_list.csv")

# Change Outcome to your actual column name
outcomes <- outcome_df$Outcome
outcomes <- outcomes[!is.na(outcomes) & outcomes != ""]
outcomes <- unique(outcomes)



# Model covariates
covariates <- c(
  "RetiPhenoAge_sd"
  ,"Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0")


dir.create("PH_results", showWarnings = FALSE)
dir.create("PH_plots", showWarnings = FALSE)

check_ph_assumption <- function(outcome_name, df, covariates) {
  
  event_var <- paste0("all_IMevent_", outcome_name)
  time_var  <- paste0("all_IMevent_", outcome_name, "_year")
  
  needed_vars <- c(event_var, time_var, covariates)
  missing_vars <- setdiff(needed_vars, names(df))
  
  if (length(missing_vars) > 0) {
    return(data.frame(
      outcome = outcome_name,
      term = NA,
      chisq = NA,
      df = NA,
      p = NA,
      status = paste0("missing variables: ", paste(missing_vars, collapse = ", "))
    ))
  }
  
  dat <- df[, needed_vars]
  dat <- na.omit(dat)
  
  names(dat)[names(dat) == event_var] <- "event"
  names(dat)[names(dat) == time_var]  <- "time"
  
  dat <- dat %>% filter(time > 0)
  
  
  if (sum(dat$event == 1) == 0) {
    return(data.frame(
      outcome = outcome_name,
      term = NA,
      chisq = NA,
      df = NA,
      p = NA,
      status = "no events"
    ))
  }
  
  
  
  
  
  
  formula_txt <- paste(
    "Surv(time, event) ~",
    paste(covariates, collapse = " + ")
  )
  
  fit <- tryCatch(
    coxph(as.formula(formula_txt), data = dat, x = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      outcome = outcome_name,
      term = NA,
      chisq = NA,
      df = NA,
      p = NA,
      status = "cox model failed"
    ))
  }
  
  ph <- tryCatch(
    cox.zph(fit),
    error = function(e) NULL
  )
  
  if (is.null(ph)) {
    return(data.frame(
      outcome = outcome_name,
      term = NA,
      chisq = NA,
      df = NA,
      p = NA,
      status = "cox.zph failed"
    ))
  }
  
  ph_table <- as.data.frame(ph$table)
  ph_table$term <- rownames(ph_table)
  rownames(ph_table) <- NULL
  
  ph_table <- ph_table %>%
    rename(
      chisq = chisq,
      p = p
    ) %>%
    mutate(
      outcome = outcome_name,
      status = "ok"
    ) %>%
    select(outcome, term, chisq, df, p, status)
  
  return(ph_table)
}

# Run PH checks for all outcomes
ph_results <- map_dfr(
  outcomes,
  check_ph_assumption,
  df = df,
  covariates = covariates
)



# Summary by outcome: global PH test
ph_global_summary <- ph_results %>%
  filter(term == "GLOBAL") %>%
  mutate(
    PH_violation_p005 = ifelse(p < 0.05, TRUE, FALSE),
    PH_violation_BH_FDR005 = p.adjust(p, method = "BH") < 0.05
  )

write.csv(
  ph_global_summary,
  "PH_results/proportional_hazards_global_summary.csv",
  row.names = FALSE
)

getwd()



# 6 Restricted cubic spline (RCS) analyses for ALL outcomes
# Cox models + overall P + nonlinearity P
# ============================================================

library(survival)
library(rms)
library(dplyr)
library(purrr)


outcome_df <- read.csv("your_outcome_file.csv")

# Change column name if needed
outcomes <- outcome_df$Outcome

outcomes <- outcomes[!is.na(outcomes)]
outcomes <- outcomes[outcomes != ""]
outcomes <- unique(outcomes)

# Exposure variable
exposure_var <- "RetiPhenoAge_sd"

# Covariates
covariates <- c(
"Age", "Sex.0","Ethnicity3","Smoking.formatted","Alcohol_status.formatted","BMI","Townsend_deprivation_index_at_recruitment.0")


# Number of spline knots
n_knots <- 4

dir.create("RCS_results", showWarnings = FALSE)

# ============================================================
#  Function for spline analysis
# ============================================================

run_rcs_analysis <- function(
    outcome_name,
    df,
    exposure_var,
    covariates,
    n_knots = 4
) {
  
  event_var <- paste0("all_IMevent_", outcome_name)
  time_var  <- paste0("all_IMevent_", outcome_name, "_year")
  
  needed_vars <- c(
    event_var,
    time_var,
    exposure_var,
    covariates
  )
  
  missing_vars <- setdiff(needed_vars, names(df))
  
  if (length(missing_vars) > 0) {
    
    return(data.frame(
      outcome = outcome_name,
      n = NA,
      events = NA,
      overall_p = NA,
      nonlinearity_p = NA,
      status = paste0(
        "missing variables: ",
        paste(missing_vars, collapse = ", ")
      )
    ))
  }
  
  dat <- df[, needed_vars]
  
  dat <- na.omit(dat)
  
  names(dat)[names(dat) == event_var] <- "event"
  names(dat)[names(dat) == time_var]  <- "time"
  
  dat <- dat %>%
    filter(time > 0)
  
  # Skip only if no events
  if (sum(dat$event == 1) == 0) {
    
    return(data.frame(
      outcome = outcome_name,
      n = nrow(dat),
      events = 0,
      overall_p = NA,
      nonlinearity_p = NA,
      status = "no events"
    ))
  }
  
  # datadist required by rms
  dd_tmp <- datadist(dat)
  assign("dd_tmp", dd_tmp, envir = .GlobalEnv)
  options(datadist = "dd_tmp")
  
  # Build formula
  formula_txt <- paste0(
    "Surv(time, event) ~ ",
    "rcs(", exposure_var, ", ", n_knots, ") + ",
    paste(covariates, collapse = " + ")
  )
  
  fit <- tryCatch(
    
    cph(
      as.formula(formula_txt),
      data = dat,
      x = TRUE,
      y = TRUE,
      surv = TRUE
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    
    return(data.frame(
      outcome = outcome_name,
      n = nrow(dat),
      events = sum(dat$event == 1),
      overall_p = NA,
      nonlinearity_p = NA,
      status = "model failed"
    ))
  }
  
  # ANOVA for spline terms
  a <- tryCatch(
    anova(fit),
    error = function(e) NULL
  )
  
  if (is.null(a)) {
    
    return(data.frame(
      outcome = outcome_name,
      n = nrow(dat),
      events = sum(dat$event == 1),
      overall_p = NA,
      nonlinearity_p = NA,
      status = "anova failed"
    ))
  }
  
  a_df <- as.data.frame(a)
  a_df$term <- rownames(a_df)
  
  # Overall exposure P
  overall_row <- a_df %>%
    filter(term == exposure_var)
  
  # Nonlinear component P
  nonlinear_row <- a_df %>%
    filter(grepl("Nonlinear", term))
  
  overall_p <- ifelse(
    nrow(overall_row) > 0,
    overall_row$P[1],
    NA
  )
  
  nonlinear_p <- ifelse(
    nrow(nonlinear_row) > 0,
    nonlinear_row$P[1],
    NA
  )
  
  return(data.frame(
    outcome = outcome_name,
    n = nrow(dat),
    events = sum(dat$event == 1),
    overall_p = overall_p,
    nonlinearity_p = nonlinear_p,
    status = "ok"
  ))
}

# ============================================================
# Run analyses for ALL outcomes
# ============================================================

rcs_results <- map_dfr(
  outcomes,
  run_rcs_analysis,
  df = df,
  exposure_var = exposure_var,
  covariates = covariates,
  n_knots = n_knots
)

# ============================================================
#  Multiple-testing correction
# ============================================================

rcs_results <- rcs_results %>%
  
  mutate(
    
    overall_FDR_p = p.adjust(overall_p, method = "BH"),
    
    nonlinearity_FDR_p = p.adjust(
      nonlinearity_p,
      method = "BH"
    ),
    
    significant_overall = overall_FDR_p < 0.05,
    
    significant_nonlinearity =
      nonlinearity_FDR_p < 0.05
  )

# ============================================================
#  Save results
# ============================================================

write.csv(
  rcs_results,
  "RCS_results/restricted_cubic_spline_all_outcomes.csv",
  row.names = FALSE
)




