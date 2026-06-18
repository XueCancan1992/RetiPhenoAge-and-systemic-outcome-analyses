


# 1. Logisic regression

## Results organize 
export_logistic_OR <- function(model, filename = "logistic_results.csv") {
  coefs <- summary(model)$coefficients
  OR <- exp(coefs[, 1])
  CI_lower <- exp(coefs[, 1] - 1.96 * coefs[, 2])
  CI_upper <- exp(coefs[, 1] + 1.96 * coefs[, 2])
  pval <- coefs[, 4]
  results <- data.frame(
    Variable = rownames(coefs),
    OR = OR,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    P_value = pval,
    row.names = NULL
  )
  write.csv(results, file = filename, row.names = FALSE)
  
  message("Saved results to: ", filename)
  return(results)
}



SEED_retiPhenoAge<-read.csv("your_dataset.csv")



### Association




# Diabetes

SEED_dm<-subset(SEED_retiPhenoAge, SEED_retiPhenoAge$dm5==0) #exclude baseline cases


dm_m1 <- glm(dm5_2 ~ 
             Retiphenoage_SD,
           data = SEED_dm,
           family = binomial(link = "logit"))

dm_m2 <- glm(dm5_2 ~ Retiphenoage_SD+age+gender,
           data = SEED_dm,
           family = binomial(link = "logit"))


dm_m3 <- glm(dm5_2  ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi,
             data = SEED_dm,
             family = binomial(link = "logit"))

dm_m4 <- glm(dm5_2 ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi+ethnicity,
             data = SEED_dm,
             family = binomial(link = "logit"))


export_logistic_OR(dm_m1, "dm_m1.csv")
export_logistic_OR(dm_m2, "dm_m2.csv")
export_logistic_OR(dm_m3, "dm_m3.csv")
export_logistic_OR(dm_m4, "dm_m4.csv")




#hypertension

SEED_HT<-subset(SEED_retiPhenoAge, SEED_retiPhenoAge$hypertension==0)#exclude baseline cases




HT_m1 <- glm(hypertension_2~
               Retiphenoage_SD,
             data = SEED_HT,
             family = binomial(link = "logit"))

HT_m2 <- glm(hypertension_2 ~ Retiphenoage_SD+age+gender,
             data = SEED_HT,
             family = binomial(link = "logit"))


HT_m3 <- glm(hypertension_2  ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi,
             data = SEED_HT,
             family = binomial(link = "logit"))

HT_m4 <- glm(hypertension_2 ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi+ethnicity,
             data = SEED_HT,
             family = binomial(link = "logit"))


export_logistic_OR(HT_m1, "HT_m1.csv")
export_logistic_OR(HT_m2, "HT_m2.csv")
export_logistic_OR(HT_m3, "HT_m3.csv")
export_logistic_OR(HT_m4, "HT_m4.csv")






#CKD

SEED_CKD<-subset(SEED_retiPhenoAge, SEED_retiPhenoAge$CKD_EPI ==0)



CKD_m1 <- glm(CKD_EPI_2~
               Retiphenoage_SD,
             data = SEED_CKD,
             family = binomial(link = "logit"))

CKD_m2 <- glm(CKD_EPI_2 ~ Retiphenoage_SD+age+gender,
             data = SEED_CKD,
             family = binomial(link = "logit"))


CKD_m3 <- glm(CKD_EPI_2  ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi,
             data = SEED_CKD,
             family = binomial(link = "logit"))

CKD_m4 <- glm(CKD_EPI_2 ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi+ethnicity,
             data = SEED_CKD,
             family = binomial(link = "logit"))


export_logistic_OR(CKD_m1, "CKD_m1.csv")
export_logistic_OR(CKD_m2, "CKD_m2.csv")
export_logistic_OR(CKD_m3, "CKD_m3.csv")
export_logistic_OR(CKD_m4, "CKD_m4.csv")



# Hyperlipidemia


SEED_hyperlipid<-subset(SEED_retiPhenoAge, SEED_retiPhenoAge$hyperlipidaemia_1 ==0)




hyperlipid_m1 <- glm(hyperlipidaemia_2~
                Retiphenoage_SD,
              data = SEED_hyperlipid,
              family = binomial(link = "logit"))

hyperlipid_m2 <- glm(hyperlipidaemia_2 ~ Retiphenoage_SD+age+gender,
              data = SEED_hyperlipid,
              family = binomial(link = "logit"))


hyperlipid_m3 <- glm(hyperlipidaemia_2  ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi,
              data = SEED_hyperlipid,
              family = binomial(link = "logit"))

hyperlipid_m4 <- glm(hyperlipidaemia_2 ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi+ethnicity,
              data = SEED_hyperlipid,
              family = binomial(link = "logit"))


export_logistic_OR(hyperlipid_m1, "hyperlipid_m1.csv")
export_logistic_OR(hyperlipid_m2, "hyperlipid_m2.csv")
export_logistic_OR(hyperlipid_m3, "hyperlipid_m3.csv")
export_logistic_OR(hyperlipid_m4, "hyperlipid_m4.csv")



##2 test of interaction


dm <- glm(dm5_2 ~ Retiphenoage_SD*ethnicity+age + gender + smkcurr+alc_cat+bmi,
          data = SEED_dm,
          family = binomial(link = "logit"))

HT <- glm(hypertension_2 ~ Retiphenoage_SD*ethnicity+age + gender + smkcurr+alc_cat+bmi,
             data = SEED_HT,
             family = binomial(link = "logit"))

Hyperlipidaemia <- glm(hyperlipidaemia_2 ~ Retiphenoage_SD*ethnicity+age + gender + smkcurr+alc_cat+bmi,
             data = SEED_hyperlipid,
             family = binomial(link = "logit"))


CKD <- glm(CKD_EPI__2  ~ Retiphenoage_SD*ethnicity+age + gender + smkcurr+alc_cat+bmi,
             data = SEED_CKD,
             family = binomial(link = "logit"))


##3 Subgroup analyses for Diabetes



DM_M<-subset(SEED_dm, SEED_dm$ethnicity=="Malay") 
DM_I<-subset(SEED_dm, SEED_dm$ethnicity=="Indian") 
DM_C<-subset(SEED_dm, SEED_dm$ethnicity=="Chinese") 



DM_Malay <- glm(dm5_2  ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi,
                     data = DM_M,
                     family = binomial(link = "logit"))


DM_Indian <- glm(dm5_2  ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi,
                data = DM_I,
                family = binomial(link = "logit"))

DM_Chinese <- glm(dm5_2  ~ Retiphenoage_SD+age + gender + smkcurr+alc_cat+bmi,
                data = DM_C,
                family = binomial(link = "logit"))





