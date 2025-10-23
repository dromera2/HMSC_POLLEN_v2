targets::tar_load(combined_eval)

format.data <- function(combined_eval, group1, group2, metric) {
  data <- combined_eval %>% 
    filter(group %in% c(group1, group2)) %>% 
    filter(variable == metric) %>% 
    select(short_species, group, value) %>% 
    group_by(short_species, group) %>% 
    summarize(value = mean(value), .groups = "drop") %>% 
    pivot_wider(names_from = group, values_from = value)
  
  colnames(data) <- c("id", "A", "B")
  
  data
}


normality.check <- function(data){
  diff <- data$A - data$B
  shap <- shapiro.test(diff)
  shap
}

exp_AUC <- format.data(combined_eval, "EXP_JSDM", "EXP_SDM", "AUC")
exp_TSS <- format.data(combined_eval, "EXP_JSDM", "EXP_SDM", "TSS")
pred_AUC <- format.data(combined_eval, "PRED_JSDM", "PRED_SDM", "AUC")
pred_TSS <- format.data(combined_eval, "PRED_JSDM", "PRED_SDM", "TSS")


normality.check(exp_AUC)
normality.check(exp_TSS)
normality.check(pred_AUC)
normality.check(pred_TSS)


t.test(exp_AUC$A, exp_AUC$B, paired = TRUE)
t.test(exp_TSS$A, exp_TSS$B, paired = TRUE)
t.test(pred_AUC$A, pred_AUC$B, paired = TRUE)
wilcox.test(pred_TSS$A, pred_TSS$B, paired = TRUE)
