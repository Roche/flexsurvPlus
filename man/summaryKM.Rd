% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryKM.R
\name{summaryKM}
\alias{summaryKM}
\title{Extract information about non-parametric survival models}
\usage{
summaryKM(
  data,
  time_var,
  event_var,
  weight_var = "",
  strata_var,
  int_name,
  ref_name,
  types = c("survival", "cumhaz", "median", "rmst"),
  t = NULL,
  ci = FALSE,
  se = FALSE,
  ...
)
}
\arguments{
\item{data}{A data frame containing individual patient data for the relevant
time to event outcomes.}

\item{time_var}{Name of time variable in 'data'. Variable must be numerical and >0.}

\item{event_var}{Name of event variable in 'data'. Variable must be
numerical and contain 1's to indicate an event and 0 to indicate a censor.}

\item{weight_var}{Optional name of a variable in "data" containing case weights.}

\item{strata_var}{Name of stratification variable in "data". This is usually
the treatment variable and must be categorical. Not required if only one
arm is being analyzed.}

\item{int_name}{Character to indicate the name of the treatment of interest,
must be a level of the "strata_var" column in "data", used for labelling
the parameters.}

\item{ref_name}{Character to indicate the name of the reference treatment,
must be a level of the "strata_var" column in "data", used for labelling
the parameters. Not required if only one arm is being analyzed.}

\item{types}{A list of statistics to extract - options include "survival",
"cumhaz", "median", and "rmst". For details see the vignette on descriptive
analysis.}

\item{t}{The time points to be used - this only controls the rmst statistic.}

\item{ci}{Should a confidence interval be returned (TRUE or FALSE)}

\item{se}{Should a standard error be returned (TRUE or FALSE)}

\item{...}{Additional arguments passed to \code{\link{survfit}}}
}
\value{
A data frame containing the following values and similar to that returned by \code{\link{summaryPSM}}
\itemize{
\item Model - returned as "Kaplan Meier"
\item ModelF - an ordered factor of Model
\item Dist - returned as "Kaplan Meier"
\item DistF - an ordered factor of Dist
\item distr - returned as "km"
\item Strata - Either Intervention or Reference
\item StrataName - As specified by int_name and ref_name respectively.
\item type - as defined by the types parameter.
\item variable - "est", "lcl", "ucl", "se" respectively
\item time - either NA or the time the statistic is evaluated at
\item value - estimated value
}
}
\description{
Extract information about non-parametric survival models
}
\examples{

require(dplyr)
require(ggplot2)

PFS_data <- sim_adtte(seed = 2020, rho = 0.6) \%>\%
filter(PARAMCD=="PFS") \%>\%
transmute(USUBJID,
            ARMCD,
            PFS_days = AVAL,
            PFS_event = 1- CNSR,
            wt = runif(500,0,1)
)

pfs_info <- summaryKM(
  data = PFS_data,
  time_var = "PFS_days",
  event_var = "PFS_event",
  strata_var = "ARMCD",
  int_name = "A",
  ref_name = "B",
  ci = TRUE,
  t = c(500, 700))

ggplot(data = filter(pfs_info, type == "survival", variable == "est"),
       aes(x = time, y = value, color = StrataName)) +
  geom_step() +
  geom_step(data = filter(pfs_info, type == "survival", variable == "lcl"), linetype = 2) +
  geom_step(data = filter(pfs_info, type == "survival", variable == "ucl"), linetype = 2) +
  geom_point(data = filter(pfs_info, type == "survival", variable == "censored")) +
  xlab("Time") +
  ylab("Survival") +
  ggtitle("KM estimates and 95\% CI")

filter(pfs_info, type == "cumhaz", variable == "est") \%>\%
  ggplot(aes(x = time, y = value, color = StrataName)) +
  geom_step() +
  xlab("Time") +
  ylab("Cumulative hazard") 
 
filter(pfs_info, type == "median") \%>\%
  transmute(StrataName, variable, value)
  
filter(pfs_info, type == "rmst") \%>\%
  transmute(StrataName, variable, time, value)
 
# example with weights
 pfs_info_wt <- summaryKM(
   data = PFS_data,
   time_var = "PFS_days",
   event_var = "PFS_event",
   strata_var = "ARMCD",
   weight_var = "wt",
   int_name = "A",
   ref_name = "B",
   types = "survival"
   )
   
   ggplot(data = filter(pfs_info, type == "survival", variable == "est"),
          aes(x = time, y = value, color = StrataName)) +
     geom_step(aes(linetype = "Original")) +
     geom_step(data = filter(pfs_info_wt, type == "survival", variable == "est"), 
               aes(linetype = "Weighted")) +
     xlab("Time") +
     ylab("Survival") +
     ggtitle("KM estimates and 95\% CI")
   
}
