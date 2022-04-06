#' Extract information about non-parametric survival models
#'
#'
#' @param data A data frame containing individual patient data for the relevant
#'   time to event outcomes.  
#' @param time_var Name of time variable in 'data'. Variable must be numerical and >0.
#' @param event_var Name of event variable in 'data'. Variable must be
#'   numerical and contain 1's to indicate an event and 0 to indicate a censor.
#' @param weight_var Optional name of a variable in "data" containing case weights.
#' @param strata_var Name of stratification variable in "data". This is usually
#'   the treatment variable and must be categorical. Not required if only one 
#'   arm is being analyzed.
#' @param int_name Character to indicate the name of the treatment of interest,
#'   must be a level of the "strata_var" column in "data", used for labelling
#'   the parameters.
#' @param ref_name Character to indicate the name of the reference treatment,
#'    must be a level of the "strata_var" column in "data", used for labelling
#'    the parameters. Not required if only one arm is being analyzed.
#' @param types A list of statistics to extract - options include "survival", 
#'    "cumhaz", "median", and "rmst". For details see the vignette on descriptive
#'     analysis.
#' @param t The time points to be used - this only controls the rmst statistic.
#' @param ci Should a confidence interval be returned (TRUE or FALSE)
#' @param se Should a standard error be returned (TRUE or FALSE)
#' @param ... Additional arguments passed to \code{\link{survfit}}
#'
#' 
#' @return A data frame containing the following values and similar to that returned by \code{\link{summaryPSM}}
#' \itemize{
#'   \item Model - returned as "Kaplan Meier"
#'   \item ModelF - an ordered factor of Model
#'   \item Dist - returned as "Kaplan Meier"
#'   \item DistF - an ordered factor of Dist
#'   \item distr - returned as "km"
#'   \item Strata - Either Intervention or Reference
#'   \item StrataName - As specified by int_name and ref_name respectively.
#'   \item type - as defined by the types parameter.
#'   \item variable - "est", "lcl", "ucl", "se" respectively
#'   \item time - either NA or the time the statistic is evaluated at
#'   \item value - estimated value
#'   }
#'   
#' @examples
#' 
#' require(dplyr)
#' require(ggplot2)
#' 
#' PFS_data <- sim_adtte(seed = 2020, rho = 0.6) %>%
#' filter(PARAMCD=="PFS") %>%
#' transmute(USUBJID,
#'             ARMCD,
#'             PFS_days = AVAL,
#'             PFS_event = 1- CNSR,
#'             wt = runif(500,0,1)
#' )
#' 
#' pfs_info <- summaryKM(
#'   data = PFS_data,
#'   time_var = "PFS_days",
#'   event_var = "PFS_event",
#'   strata_var = "ARMCD",
#'   int_name = "A",
#'   ref_name = "B",
#'   ci = TRUE,
#'   t = c(500, 700))
#'
#' ggplot(data = filter(pfs_info, type == "survival", variable == "est"),
#'        aes(x = time, y = value, color = StrataName)) +
#'   geom_step() +
#'   geom_step(data = filter(pfs_info, type == "survival", variable == "lcl"), linetype = 2) +
#'   geom_step(data = filter(pfs_info, type == "survival", variable == "ucl"), linetype = 2) +
#'   geom_point(data = filter(pfs_info, type == "survival", variable == "censored")) +
#'   xlab("Time") +
#'   ylab("Survival") +
#'   ggtitle("KM estimates and 95% CI")
#'
#' filter(pfs_info, type == "cumhaz", variable == "est") %>%
#'   ggplot(aes(x = time, y = value, color = StrataName)) +
#'   geom_step() +
#'   xlab("Time") +
#'   ylab("Cumulative hazard") 
#'  
#' filter(pfs_info, type == "median") %>%
#'   transmute(StrataName, variable, value)
#'   
#' filter(pfs_info, type == "rmst") %>%
#'   transmute(StrataName, variable, time, value)
#'  
#' # example with weights
#'  pfs_info_wt <- summaryKM(
#'    data = PFS_data,
#'    time_var = "PFS_days",
#'    event_var = "PFS_event",
#'    strata_var = "ARMCD",
#'    weight_var = "wt",
#'    int_name = "A",
#'    ref_name = "B",
#'    types = "survival"
#'    )
#'    
#'    ggplot(data = filter(pfs_info, type == "survival", variable == "est"),
#'           aes(x = time, y = value, color = StrataName)) +
#'      geom_step(aes(linetype = "Original")) +
#'      geom_step(data = filter(pfs_info_wt, type == "survival", variable == "est"), 
#'                aes(linetype = "Weighted")) +
#'      xlab("Time") +
#'      ylab("Survival") +
#'      ggtitle("KM estimates and 95% CI")
#'    
#' @export
summaryKM <- function(data,
                      time_var, event_var, weight_var = "",
                      strata_var,
                      int_name, ref_name,
                      types = c("survival", "cumhaz", 
                                "median", "rmst"),
                      t = NULL,
                      ci = FALSE,
                      se = FALSE,
                      ...
                      ){
 
  # check the data passed
  assertthat::assert_that(
    !missing(data),
    msg = "data parameter must be specified"
  )
  
  assertthat::assert_that(
    !missing(time_var),
    msg = "time_var parameter must be specified"
  )
  
  assertthat::assert_that(
    !missing(event_var),
    msg = "event_var parameter must be specified"
  )
  
  
  # check if we have a single arm or two arm specification
  if (!missing(strata_var) | !missing(ref_name)){
    model.type <- "two.arm"
  } else {
    model.type <- "one.arm"
  }
    
  # prepare the appropriate standard data
  if (model.type == "two.arm"){
    
    # additional call checks
    assertthat::assert_that(
      !missing(strata_var) & !missing(ref_name) & !missing(int_name),
      msg = "For analysis of data with two arms strata_var, ref_name and int_name must be specified."
    )
    
    # prepare and validate data
    data_standard <- Format_data(data = data, time_var = time_var, event_var = event_var, weight_var = weight_var, 
                              strata_var = strata_var, int_name = int_name, ref_name = ref_name)
    model.formula <- survival::Surv(Time, Event==1) ~ ARM
  }
  
  if (model.type == "one.arm"){
    
    # additional call checks
    assertthat::assert_that(
      missing(strata_var) & missing(ref_name) & !missing(int_name),
      msg = "For analysis of data with one arm only int_name must be specified. strata_var and ref_name should not be specified."
    )
    
    # prepare and validate data
    data_standard <- Format_data_onearm(data = data, time_var = time_var, event_var = event_var, weight_var = weight_var, 
                                     int_name = int_name)
    model.formula <- survival::Surv(Time, Event==1) ~ 1
    
    # fake assign for later
    ref_name <- ""
    
  }
  
  # impute t if not provided
  if (is.null(t)){
    t <- seq(0, max(data_standard$Time), length.out = 10)
  }
  
  
  # fix for binding checks
  Weight <- Strata <- variable <- value <- type <- ARM <- time <- NULL
  v2 <- maxtime <- survival <- lag_time <- lag_surv <- area <- cum_area <- NULL
  Model <- Dist <- NULL
  
  # get a survfit object 
  # check if a Weight  was specified
  # as standard format names are fixed
  
  if ("Weight" %in% names(data_standard)) {
    sf <- survival::survfit(model.formula, data = data_standard, weights = Weight, ...)
  } else {
    sf <- survival::survfit(model.formula, data = data_standard, ...)
  }
  
  # add a fake strata so processing is the same for both data types 
  
  if(model.type == "one.arm"){
    sf$strata <- length(sf$surv)
    names(sf$strata) <- "ARM=Int"
  }
    
  # extract initial statistics from the survfit object
  
  # define a structure for return
  rc_struct <- tibble::tibble(
    ARM = "",
    type = "",
    variable = "",
    time = NaN,
    value = NaN) %>%
    dplyr::filter(1 == 2)
  
  # basis return
  rc <- rc_struct
  
  ##############################################
  # survival statistic
  ##############################################
  
  if ("survival" %in% types){
    
    rc_surv <- tibble::tibble(
      ARM = rep(gsub("ARM=","",names(sf$strata), fixed = TRUE), sf$strata),
      type = "survival",
      variable = "est",
      time = sf$time,
      value = sf$surv)
    
    # include censored obs for plotting
    
    smry_cens <- summary(sf, censored = TRUE)
    
    if(model.type == "one.arm"){
      smry_cens$strata <- "ARM=Int"
    }
    
    rc_cens <- tibble::tibble(
      ARM = gsub("ARM=","",smry_cens$strata, fixed = TRUE),
      type = "survival",
      variable = "censored",
      time = smry_cens$time,
      value = smry_cens$surv,
      v2 = smry_cens$n.censor
    ) %>%
      dplyr::filter(v2 > 0) %>%
      dplyr::select(-v2)
  
    
    rc_surv <- rc_surv %>%
      dplyr::bind_rows(rc_cens)

    if (ci == TRUE){
      rc_surv <- 
        rc_surv %>%
        dplyr::filter(variable == "est") %>%
        dplyr::mutate(variable = "lcl",
                      value = sf$lower) %>%
        dplyr::bind_rows(rc_surv)
      
      rc_surv <- 
        rc_surv %>%
        dplyr::filter(variable == "est") %>%
        dplyr::mutate(variable = "ucl",
                      value = sf$upper) %>%
        dplyr::bind_rows(rc_surv)
    }
    
    if (se == TRUE){
      rc_surv <- 
        rc_surv %>%
        dplyr::filter(variable == "est") %>%
        dplyr::mutate(variable = "se",
                      value = sf$std.err) %>%
        dplyr::bind_rows(rc_surv)
      }
  
    # add a 0 time with survival of 1
    
    rc_surv <- rc_surv %>%
      dplyr::filter(variable %in% c("est", "lcl", "ucl")) %>%
      dplyr::group_by(ARM, type, variable) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(time = 0, value = 1) %>%
      dplyr::ungroup() %>%
      dplyr::bind_rows(rc_surv) 
    
    rc <- rc %>%
      dplyr::bind_rows(rc_surv)
    
  } # end survival
  
  ##############################################
  # cum hazard statistic
  ##############################################
  
  if ("cumhaz" %in% types){
    
    rc_ch <- tibble::tibble(
      ARM = rep(gsub("ARM=","",names(sf$strata), fixed = TRUE), sf$strata),
      type = "cumhaz",
      variable = "est",
      time = sf$time,
      value = sf$cumhaz)
    
    # add a 0 time with cum hazard of 0
    
    rc_ch <- rc_ch %>%
      dplyr::filter(variable %in% c("est", "lcl", "ucl")) %>%
      dplyr::group_by(ARM, type, variable) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(time = 0, value = 0) %>%
      dplyr::ungroup() %>%
      dplyr::bind_rows(rc_ch) 
    
    # add to the return dataset
    rc <- rc %>%
      dplyr::bind_rows(rc_ch)
    
  }
  
  
  ##############################################
  # median statistic
  ##############################################
  
  if ("median" %in% types){
  
    qnt <- stats::quantile(sf, probs = c(0.5), conf.int = TRUE)
    
    rc_median <- tibble::tibble(
      ARM = gsub("ARM=","",row.names(qnt$quantile), fixed = TRUE),
      type = "median",
      variable = "est",
      value = as.numeric(qnt$quantile))
    
    if (ci == TRUE){
      rc_median <- 
        rc_median %>%
        dplyr::filter(variable == "est") %>%
        dplyr::mutate(variable = "lcl",
                      value = as.numeric(qnt$lower)) %>%
        dplyr::bind_rows(rc_median)
      
      rc_median <- 
        rc_median %>%
        dplyr::filter(variable == "est") %>%
        dplyr::mutate(variable = "ucl",
                      value = as.numeric(qnt$upper)) %>%
        dplyr::bind_rows(rc_median)
    }
    
    # add to the return dataset
    rc <- rc %>%
      dplyr::bind_rows(rc_median)
    
  }
  
  ##############################################
  # rmst statistic
  # TBD if another package should be used here 
  # or if this should be a standalone function
  ##############################################
  
  if ("rmst" %in% types){
    
    temp_surv <- tibble::tibble(
      ARM = rep(gsub("ARM=","",names(sf$strata), fixed = TRUE), sf$strata),
      time = sf$time,
      survival = sf$surv) 
    
    # add a time 0 with survival of 1
    temp_surv <- temp_surv %>%
      dplyr::group_by(ARM) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(time = 0, survival = 1) %>%
      dplyr::ungroup() %>%
      dplyr::bind_rows(temp_surv) %>%
      dplyr::arrange(ARM, time)
    
 
    for (this.t in t){
      
      rc_rmst <- temp_surv %>%
        dplyr::group_by(ARM) %>%
        dplyr::mutate(maxtime = max(time)) %>%
        dplyr::filter(time < this.t) 
        
      # add a time for rmst
      # if beyond max time then adding NA will cause later area to be also NA
      rc_rmst <- rc_rmst %>%
        dplyr::slice(1) %>%
        dplyr::mutate(time = ifelse(this.t <= maxtime, this.t, NA),
                      survival = NA) %>%
        dplyr::bind_rows(rc_rmst) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ARM, time, dplyr::desc(survival))
        
      # calculate the area
      
      rc_rmst <- rc_rmst %>%
        dplyr::group_by(ARM) %>%
        dplyr::mutate(lag_time = dplyr::lag(time),
                      lag_surv = dplyr::lag(survival),
          diff = time - lag_time,
          area = ifelse(time == 0, 0, diff*lag_surv)) %>%
        dplyr::summarise(
          cum_area = sum(area),
          ) %>%
        dplyr::ungroup() 
      
      # add to return data set
      
      rc <- rc_rmst %>%
        dplyr::transmute(
          ARM,
          type = "rmst",
          variable = "est",
          time = this.t,
          value = cum_area) %>%
        dplyr::bind_rows(rc)
        
    }
    
  }
  
  
  ##################################################
  # decode the strata names and return dataset
  ##################################################
  
  # define factors for models (consistent with other functions)
  Model.levels <- c("Kaplan Meier",
                    "Common shape", "Independent shape",
                    "Separate - Reference", "Separate - Intervention",
                    "One arm - Intervention")
  
  Dist.orig <- c("exp", "weibull", "lnorm", "gamma", "gengamma", 
                 "genf", "llogis", "gompertz", "Kaplan Meier")
  Dist.levels <- c("Exponential", "Weibull", "Log Normal", "Gamma", "Generalized Gamma", 
                   "Generalized F", "Log Logistic", "Gompertz", "Kaplan Meier")
  
  
  rc <- rc %>%
    dplyr::transmute(
      Model = "Kaplan Meier",
      ModelF = factor(Model, levels = Model.levels, ordered = TRUE),
      Dist = "Kaplan Meier",
      DistF = factor(Dist, levels = Dist.levels, ordered = TRUE),
      distr = "km",
      Strata = ifelse(ARM=="Int", "Intervention", "Reference"),
      StrataName = ifelse(Strata == "Intervention", int_name, ref_name),
      type,
      variable,
      time,
      value) %>%
    dplyr::arrange(
      Strata, type, variable, time
    )
  
  return(rc)
  
}





