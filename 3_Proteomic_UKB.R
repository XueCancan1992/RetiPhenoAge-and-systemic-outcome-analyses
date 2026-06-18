
#1. Data preparation: QC, imputation

# library ----
library(Amelia)
library(bestNormalize)
library(readr)
# load olink data ----
olink = readr::read_csv(
  'olink.csv',
  col_types=paste(
    'iiiccc',
    paste(rep('d',length(7:2929)),collapse=''),
    sep=''))
olink = as.data.frame(olink)
colnames(olink)[7:2929] = gsub('olink([0-9]|[0-9][0-9])_','',colnames(olink)[7:2929])

olink_assay = readr::read_csv('olink_assay.csv',col_types='ccc')
olink_assay = as.data.frame(olink_assay)

olink_qc = readr::read_csv('olink_qc.csv',col_types='ccccidciccciD')
olink_qc = as.data.frame(olink_qc)

# assay grouping in txt file
cat = read.csv('olink_assay.dat.txt',header=TRUE,sep='\t')
cat$Panel = gsub(' II','',cat$Panel)

# missing data filtering (1st layer) ----
ds = olink
# Exclude proteins that are missing in >20% of participants.
ds = ds[,unique( c(1:6, (1:ncol(ds))[colSums(is.na(ds))/nrow(ds) <= .2]) )]
# Exclude participants with >20% missing values.
ds = ds[rowSums(is.na(ds[,7:ncol(ds)]))/(ncol(ds)-6) <= .2,]
dim(ds)

# split olink by panel ----
olink_assay = cat
d = NULL
i = 1
for (p in unique(olink_assay$Panel)) {
  assay = olink_assay$Assay[olink_assay$Panel %in% p]
  assay = sort(colnames(ds)[toupper(colnames(ds)) %in% toupper(gsub('[-]','_',assay))])
  dt = NULL
  dt$panel = p
  dt$data = ds[,assay]
  d[[i]] = dt
  rm(dt)
  i = i+1
}

# count protein by panel ----
{
  cbind(sapply(1:length(d),function(i) d[[i]]$panel),
        sapply(1:length(d),function(i) dim(d[[i]]$data)[2]))
}
# missing data filtering (2nd layer) ----
# Exclude participants with >20% missing values in each panel group
id = rownames(ds)
for (i in 1:length(d)) {
  dt = d[[i]]$data
  # Exclude participants with >20% missing values.
  dt = dt[rowSums(is.na(dt))/ncol(dt) <= .2,]
  d[[i]]$data = dt
  id = id[id %in% rownames(dt)]
  rm(dt)
}
# structure data d ----
for (i in 1:length(d)) {
  dt = d[[i]]$data
  dt = dt[id,]
  d[[i]]$data = dt
  rm(dt)
}
d[[1]]$id = ds[id,1:6]

rm(list=ls()[ls()!='d'])
gc()

sapply(1:length(d),function(i) dim(d[[i]]$data)[1])

# transform data: Yeo-Johnson Normalization
# imputation: joint modeling imputation by multinormal distribution
# inverse transform: inverse Yeo-Johnson transformation
{
  yjparm = NULL
  yjp = data.frame(lambda=NA,mean=NA,sd=NA,norm_stat=NA)[-1,]
  for (i in 1:length(d)) yjparm[[i]] = yjp
  rm(yjp)
  
  seed = 250131
  for (i in 1:length(d)) {
    dt = d[[i]]$data
    x.t = dt
    d[[i]]$imp1 = dt
    d[[i]]$imp2 = dt
    d[[i]]$imp3 = dt
    d[[i]]$imp4 = dt
    d[[i]]$imp5 = dt
    
    gc()
    # Yeo-Johnson Normalization
    yjp = data.frame(lambda=NA,mean=NA,sd=NA,norm_stat=NA)[-1,]
    for (c in 1:ncol(dt)) {
      yj = yeojohnson(dt[,c])
      yjp[c,] = NA
      yjp[c,'lambda'] = yj$lambda
      yjp[c,'mean'] = yj$mean
      yjp[c,'sd'] = yj$sd
      yjp[c,'norm_stat'] = yj$norm_stat
      x.t[,c] = yj$x.t
      yj$x = NULL
      yj$x.t = NULL
      yj$n = NULL
    }
    yjparm[[i]] = yjp
    
    # joint modeling imputation by multinormal distribution
    dt = x.t
    gc()
    set.seed(seed)
    imp = amelia(dt,m=5,p2s=0)$imputations
    
    gc()
    # inverse Yeo-Johnson transformation
    for (c in 1:ncol(dt)) {
      yj$lambda = yjparm[[i]]$lambda[c]
      yj$mean = yjparm[[i]]$mean[c]
      yj$sd = yjparm[[i]]$sd[c]
      yj$norm_stat = yjparm[[i]]$norm_stat[c]
      
      d[[i]]$imp1[,c] = predict(yj,imp$imp1[,c],inverse=TRUE)
      d[[i]]$imp2[,c] = predict(yj,imp$imp2[,c],inverse=TRUE)
      d[[i]]$imp3[,c] = predict(yj,imp$imp3[,c],inverse=TRUE)
      d[[i]]$imp4[,c] = predict(yj,imp$imp4[,c],inverse=TRUE)
      d[[i]]$imp5[,c] = predict(yj,imp$imp5[,c],inverse=TRUE)
    }
    
    seed = seed+1
    rm(dt,x.t,imp,yj,yjp)
  }
}
olink_imputed = d



### 2. Regression analysis

# merge the olink data with RetiPhenoAge and covariates before analyses

# adjusted covariates
zs = c('Age','Sex','Ethnicity',
       'PPP_cohort','Batch','centre','genotyping_array',
       'Fasting_time','Season_blood_sample_collected','Year_Blood_Sample_To_Olink',
       'BMI','Hypertension','Diabetes.comb','Smoking.formatted','Alcohol_status.formatted',
       'Townsend_deprivation_index_at_recruitment')

regression_output = NULL
for (i in 60:2970) { # 60:2970 are the columns of proteomic variables
  x = colnames(d0)[i]
  beta_i = NULL
  for (k in 1:length(imp)) {
    eval(parse(text=
                 paste("beta_i = rbind(beta_i,summary(lm(retiphenoage~",
                       paste(c(x,zs),collapse='+'),
                       ",imp[[",k,"]]))$coef[2,1:2])",sep='')))
  }
  # Rubin's rule: beta_MI = mean beta_i; var(beta_MI)=mean(var(beta_i))+(1+1/M)sum((beta_i-beta_MI)^2/(M-1))
  beta_MI5 = mean(beta_i[1:5,1])
  beta_MI5.se = sqrt(mean(beta_i[1:5,2]^2)+(6/5)*sum((beta_i[1:5,1]-beta_MI5)^2/4))
  beta_MI5.lcl = beta_MI5 + qnorm(.025)*beta_MI5.se
  beta_MI5.ucl = beta_MI5 + qnorm(.975)*beta_MI5.se
  beta_MI5.p = pnorm(-abs(beta_MI5/beta_MI5.se),0,1,lower.tail=TRUE,log.p=FALSE)*2
  
  regression_output = rbind(
    regression_output,
    cbind(x,
          beta_MI5,beta_MI5.se,beta_MI5.lcl,beta_MI5.ucl,beta_MI5.p))
  rm(x,beta_i,
     beta_MI5,beta_MI5.se,beta_MI5.lcl,beta_MI5.ucl,beta_MI5.p)
}

