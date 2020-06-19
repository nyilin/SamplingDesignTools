#' Simulated cohort 1
#'
#' @format A data frame with 10000 rows and 5 variables:
#' \describe{
#'   \item{id}{Subject ID of each subject.}
#'   \item{y}{Event/censoring indicator. Events are indicated by y=1.}
#'   \item{t}{Event/censoring time (in years). Maximum follow-up time is 25 years.}
#'   \item{age}{Age of subjects (rounded to integers).}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#' }
"cohort_1"

#' Nested case-control (NCC) data 1
#' 
#' @description Time-matched data drawn from \code{\link{cohort_1}}, with one  
#' controls matched to each case.
#' 
#' @format A data frame with 1164 rows and 6 variables:
#' \describe{
#'   \item{Set}{Set ID.}
#'   \item{Map}{Row numbers in \code{cohort_1}.}
#'   \item{Time}{Event time (in years) of the case in each \code{Set}.}
#'   \item{Fail}{Case-control indicator.}
#'   \item{age}{Age of subjects (rounded to integers).}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#' }
"ncc_1"

#' Simulated cohort 2
#'
#' @format A data frame with 100000 rows and 8 variables:
#' \describe{
#'   \item{id}{Subject ID of each subject.}
#'   \item{y}{Event/censoring indicator. Events are indicated by y=1.}
#'   \item{t}{Event/censoring time (in years). Maximum follow-up time is 25 years.}
#'   \item{x}{Binary exposure.}
#'   \item{age}{Age of subjects (rounded to integers and mean-centred).}
#'   \item{age_cat}{Age category of subjects.}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#'   \item{z}{Binary effect modifier.}
#' }
"cohort_2"

#' Nested case-control (NCC) data 2
#' 
#' @description NCC data drawn from \code{\link{cohort_2}} matched on age group
#'   and gender, with five controls matched to each case.
#' 
#' @format A data frame with 16638 rows and 10 variables:
#' \describe{
#'   \item{Set}{Set ID.}
#'   \item{Map}{Row numbers in \code{cohort_2}.}
#'   \item{Time}{Event time (in years) of the case in each \code{Set}.}
#'   \item{Fail}{Case-control indicator.}
#'   \item{age_cat}{Age category of subjects.}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#'   \item{x}{Binary exposure.}
#'   \item{age}{Age of subjects (rounded to integers).}
#'   \item{z}{Binary effect modifier.}
#'   \item{t}{Actual time/censoring (in years) of each subject.}
#' }
"ncc_2"

#' Number of subjects at risk at each time point in \code{ncc_2}
#' 
#' @format A data frame with 2773 rows and 4 variables:
#' \describe{
#'   \item{t}{Event time (in years) of cases in \code{ncc_2}.}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#'   \item{age_cat}{Age category of subjects.}
#'   \item{n.risk}{Number of subjects at risk at each event time in each gender-age strata.}
#' }
#' @details Note how column names in \code{n_at_risk} match those in
#'   \code{ncc_2}. This is critical for computing the KM-type weights for
#'   subjects in \code{ncc_2}.
"n_at_risk"