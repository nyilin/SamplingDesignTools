library(SamplingDesignTools)
# Load mini cohort
data("mini_cohort")
mini_cohort
# Manually prepare a 1:1 NCC data
ncc <- rbind(
  data.frame(Set = 1, Map = c(1, 5), Time = mini_cohort$t[1], 
             Fail = c(1, 0), t = mini_cohort$t[c(1, 5)]), 
  data.frame(Set = 2, Map = c(3, 4), Time = mini_cohort$t[3], 
             Fail = c(1, 0), t = mini_cohort$t[c(3, 4)]), 
  data.frame(Set = 3, Map = c(4, 10), Time = mini_cohort$t[4], 
             Fail = c(1, 0), t = mini_cohort$t[c(4, 10)]), 
  data.frame(Set = 4, Map = c(6, 7), Time = mini_cohort$t[6], 
             Fail = c(1, 0), t = mini_cohort$t[c(6, 7)]), 
  data.frame(Set = 5, Map = c(8, 9), Time = mini_cohort$t[8], 
             Fail = c(1, 0), t = mini_cohort$t[c(8, 9)]), 
  data.frame(Set = 6, Map = c(9, 10), Time = mini_cohort$t[9], 
             Fail = c(1, 0), t = mini_cohort$t[c(9, 10)]) 
)
rownames(ncc) <- NULL
ncc
# Map the NCC sample to the original cohort, break the matching, identify the
# subjects selected into the NCC, and return this subset with KM type weights
# computed for them.
# First create the sampling and status indicator:
sample_stat <- numeric(nrow(mini_cohort))
sample_stat[unique(ncc$Map[ncc$Fail == 0])] <- 1
sample_stat[ncc$Map[ncc$Fail == 1]] <- 2
# Then find the sampled subset and compute weights:
ncc_nodup <- compute_km_weights(
  cohort = mini_cohort, t_name = "t", y_name = "status",
  sample_stat = sample_stat, n_per_case = 1
)
ncc_nodup
# Alternatively, if the cohort is not available, the weights can be computed
# as long as number of subjects at risk at event times in each strata is
# available elsewhere, and the actual time of event/censoring is available
# for each subject in the NCC.
# Compute the number of subjects at risk from mini_cohort:
risk_table <- compute_risk_table(cohort = mini_cohort, t_name = "t", 
                                 y_name = "status")
risk_table
# The following command computes the same weights as in ncc_nodup:
ncc_nodup_v2 <- compute_km_weights(
  ncc = ncc[, -1], risk_table_manual = risk_table,
  id_name = "Map", t_match_name = "Time", t_name = "t", y_name = "Fail",
  n_per_case = 1
)
ncc_nodup_v2
