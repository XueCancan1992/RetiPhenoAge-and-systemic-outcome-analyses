
data<-readRDS("your_dataset2.rds")

independent_var <- "RetiPhenoAge_sd"
covariates <- c("Age", "Sex.0", "Ethnicity3", "BMI", "Smoking.formatted", "Alcohol_status.formatted","Townsend_deprivation_index_at_recruitment.0","BP_lowering","statin",
                "Hypertension","Hyperlipidaemia","Diabetes")
results <- data.frame(Variable = character(), Correlation = numeric(), 
                      R2 = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

results <- data.frame(
  Variable = character(), 
  Estimate = numeric(),
  SE = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  R2 = numeric(), 
  P_value = numeric(), 
  stringsAsFactors = FALSE
)



for (i in 3:75) {
  dependent_var_name <- colnames(data)[i]
  
  # Check if dependent_var_name exists and is not the same as independent_var or covariates
  if (dependent_var_name %in% c(independent_var, covariates)) next
  
  # Standardize the dependent variable
  data[[dependent_var_name]] <- scale(data[[dependent_var_name]])
  
  # Create a formula for the regression
  try({
    formula <- as.formula(
      paste(dependent_var_name, "~", independent_var, "+", paste(covariates, collapse = " + "))
    )
    
    # Fit the linear regression model
    model <- lm(formula, data = data)
    
    # Extract the coefficient for the independent variable
    if (independent_var %in% rownames(summary(model)$coefficients)) {
      coefficient <- summary(model)$coefficients[independent_var, "Estimate"]
      se <- summary(model)$coefficients[independent_var, "Std. Error"]
      p_value <- summary(model)$coefficients[independent_var, "Pr(>|t|)"]
      r2 <- summary(model)$r.squared
      
      ci <- confint(model, parm = independent_var)  
      ci_lower <- ci[1]
      ci_upper <- ci[2]
      
      # Append results to the data frame
      results <- rbind(results, data.frame(
        Variable = dependent_var_name,
        Estimate = coefficient,  
        SE = se,
        CI_lower = ci_lower,
        CI_upper = ci_upper,
        R2 = r2,
        P_value = p_value,
        stringsAsFactors = FALSE
      ))
    }
  }, silent = TRUE)  # Skip problematic variables
}   

# Check results
if (nrow(results) == 0) {
  warning("No results were generated. Check if the independent variable and covariates exist in the dataset, and if dependent variables are valid.")
}


#map with group
quanti_map<-read.csv("your_quantiative_traits_map.csv")

head(quanti_map)

results<-merge(quanti_map,results, by="Variable", all.x = TRUE)

head(results)

library(dplyr)

adjust_fdr <- function(data) {
  data %>%
    group_by(group) %>%
    mutate(fdr_P = p.adjust(P_value, method = "fdr")) %>%
    ungroup() %>%
    as.data.frame()  # Ensure the output is a data.frame
}
results <- adjust_fdr(results)

head(results)


write.csv(results, "RetiPhenoAge_multivariable_per_sd.csv", row.names = FALSE)
head(results)


