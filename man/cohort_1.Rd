% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cohort_1}
\alias{cohort_1}
\title{Simulated cohort 1}
\format{
A data frame with 10000 rows and 5 variables:
\describe{
  \item{id}{Subject ID of each subject.}
  \item{y}{Event/censoring indicator. Events are indicated by y=1.}
  \item{t}{Event/censoring time (in years). Maximum follow-up time is 25 years.}
  \item{age}{Age of subjects (rounded to integers).}
  \item{gender}{Gender of subjects (1 for male and 0 for female).}
}
}
\usage{
cohort_1
}
\description{
Survival outcome was simulated with true hazard:
\deqn{\log \{h(t)\} = \log \{h_0\} + \log(1.1)Age + \log(2)Gender.}
}
\examples{
library(survival)
data("cohort_1", package = "SamplingDesignTools")
summary(coxph(Surv(t, y) ~ age + gender, data = cohort_1))

}
\keyword{datasets}
