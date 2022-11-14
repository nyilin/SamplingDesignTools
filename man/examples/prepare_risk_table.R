library(SamplingDesignTools)
data("ncc_2")
risk_table <- prepare_risk_table(ncc = ncc_2, t_match_name = "Time", y_name = "Fail",
                                 match_var_names = c("age_cat", "gender"))
head(risk_table)
