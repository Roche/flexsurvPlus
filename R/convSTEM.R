#' Convert survival parameters to SAS/STEM parametric forms
#'
#' @param x An object created by calling \code{\link{runPSM}}
#' @param samples An object created by calling \code{\link{boot}} with \code{\link{bootPSM}}
#' @param use an optional character string giving a method for computing covariances in the presence of missing values. See  \code{\link{cov}} 
#'            for details. Option "complete.obs" maybe needed when some bootstrap samples do not converge to estimate covariance only using
#'            those that do. 
#'
#' This function primarily exists for backward compatibility with older excel models where parametric extrapolation was
#' performed with SAS and alternative parametric forms were used for distributions. As such only a subset of models are supported.
#' One or both of \code{x} and \code{samples} must be specified and affect what is returned. For more details please see the 
#' vignette \code{vignette("STEM_compatability", package = "flexsurvPlus")}
#' 
#' Possible distributions include
#' 
#' \itemize{
#'   \item Exponential ('exp')
#'   \item Weibull ('weibull')
#'   \item Gompertz ('gompertz')
#'   \item Log-normal ('lnorm')
#'   \item Log-logistic ('llogis')
#'   \item Generalized gamma ('gengamma')
#'   \item Gamma ('gamma')
#'   }
#'
#' @return a list containing 4 data frames
#' \itemize{
#'   \item stem_param Converted parameter estimates
#'   \item stem_cov Converted covariance matrix (if \code{samples} provided)
#'   \item stem_modsum Converted model summary (if \code{x} provided)
#'   \item stem_boot Converted bootstrap samples (if \code{samples} provided)
#'   }
#' @export
convSTEM <- function(x = NULL, samples = NULL, use = "everything"){
  
  #check that at least one object is provided
  assertthat::assert_that(
    !is.null(x) | !is.null(samples),
    msg = "One of x or samples must be provided"
  )
  
  # check x is an object from runPSM
  assertthat::assert_that(
    is.null(x) | any(names(x) == c("parameters_vector")),
    msg = "x must be created using runPSM"
  )
  
  # check samples is a bootstrap object
  assertthat::assert_that(
    is.null(samples) | class(samples) == "boot",
    msg = "samples must be a boot object"
  )
  
  if (!is.null(samples)){
    assertthat::assert_that(
      all(unlist(lapply(strsplit(names(samples$t0), ".", fixed = TRUE), function(x){x[1]})) %in% 
            c("comshp", "indshp", "sep", "onearm")), 
      msg = "This function only works with bootstrap samples created using boot with statistic = bootPSM"
    )
  }
  
  # check if both x and samples provided and are from same data/models
  if (!is.null(samples) & !is.null(x)) {
    assertthat::assert_that(
      all(x$parameters_vector == samples$t0, na.rm = TRUE),
      msg = "x and samples provided do not match. Please confirm these are from the same models"
    )
  }
  
  # fix no binding checks
  Estimate <- bootid <- BootSample <- Model <- Dist <- Param <- Value <- rowParam <- colParam <- NULL
  ModelF <- DistF <- ParamF <- rowParamF <- colParamF <- NULL
  AIC <- AIC_SAS <- BIC <- BIC_SAS <- Status <- NULL
  
  
  ####################################################################################
  # conversion of the parameters to a data frame
  ####################################################################################
  
  # convert estimates to a data frame
  # either from x or from samples t0
  if (!is.null(x)){
    ests <- cbind(Estimate = "Main", bootid = 0, as.data.frame(t(x$parameters_vector)))
  } else {
    ests <- cbind(Estimate = "Main", bootid = 0, as.data.frame(t(samples$t0)))
  }
  
  # if samples are provided then include these into the estimates dataframe
  if (!is.null(samples)) {
    tsamp <- samples$t
    colnames(tsamp) <- names(samples$t0)
    covests <- tibble::as_tibble(tsamp)
    covests$Estimate <- paste0("Bootstrap")
    covests$bootid <- 1:nrow(covests)
    ests <- dplyr::bind_rows(ests, covests)
  }
  
  # do the conversions
  # this is not pretty code...but hopefully transformations are transparent
  
  # set up a return objects for each parameter to convert
  rc_struct <- ests %>%
    dplyr::transmute(Estimate, bootid, Model = "", Dist = "", Param = "", Value =NaN)
  
  # empty dataframe to append converted parameters too
  rc <- rc_struct %>%
    dplyr::filter(1 == 2)
  
  #############################
  #### exponential models
  #############################
  
  # common shape
  if (any(grepl("comshp.exp.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Exponential",
                    Param = "INTERCEPT", Value = -log(ests$comshp.exp.rate.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Exponential",
                    Param = "TX(Intervention)", Value = -(log(ests$comshp.exp.rate.int) - log(ests$comshp.exp.rate.ref))
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # separate shape
  if (any(grepl("sep.exp.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Exponential",
                    Param = "INTERCEPT", Value = -log(ests$sep.exp.rate.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Exponential",
                    Param = "INTERCEPT", Value = -log(ests$sep.exp.rate.int)
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # independent shape (not really supported by STEM so handled as if separate models)
  if (any(grepl("indshp.exp.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Exponential",
                    Param = "INTERCEPT", Value = -log(ests$indshp.exp.rate.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Exponential",
                    Param = "INTERCEPT", Value = -log(ests$indshp.exp.rate.int)
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # one arm models
  if (any(grepl("onearm.exp.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Exponential",
                    Param = "INTERCEPT", Value = -log(ests$onearm.exp.rate.int)
      ) %>%
      dplyr::bind_rows(rc)
    
  }
  
  #############################
  #### weibull models
  #############################
  
  # common shape
  if (any(grepl("comshp.weibull.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Weibull",
                    Param = "INTERCEPT",
                    Value = log(ests$comshp.weibull.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Weibull",
                    Param = "TX(Intervention)",
                    Value = log(ests$comshp.weibull.scale.int) - log(ests$comshp.weibull.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Weibull",
                    Param = "SCALE",
                    Value = 1/ests$comshp.weibull.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # separate
  if (any(grepl("sep.weibull.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Weibull",
                    Param = "INTERCEPT",
                    Value = log(ests$sep.weibull.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Weibull",
                    Param = "SCALE",
                    Value = 1/ests$sep.weibull.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Weibull",
                    Param = "INTERCEPT",
                    Value = log(ests$sep.weibull.scale.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Weibull",
                    Param = "SCALE",
                    Value = 1/ests$sep.weibull.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # independent shape (not really supported by STEM so handled as if separate models)
  if (any(grepl("indshp.weibull.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Weibull",
                    Param = "INTERCEPT",
                    Value = log(ests$indshp.weibull.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Weibull",
                    Param = "SCALE",
                    Value = 1/ests$indshp.weibull.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Weibull",
                    Param = "INTERCEPT",
                    Value = log(ests$indshp.weibull.scale.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Weibull",
                    Param = "SCALE",
                    Value = 1/ests$indshp.weibull.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # one arm
  if (any(grepl("onearm.weibull.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Weibull",
                    Param = "INTERCEPT",
                    Value = log(ests$onearm.weibull.scale.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Weibull",
                    Param = "SCALE",
                    Value = 1/ests$onearm.weibull.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  #############################
  #### lognormal models
  #############################
  
  # common shape
  if (any(grepl("comshp.lnorm.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Log Normal",
                    Param = "INTERCEPT",
                    Value = ests$comshp.lnorm.meanlog.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Log Normal",
                    Param = "TX(Intervention)",
                    Value = ests$comshp.lnorm.meanlog.int - ests$comshp.lnorm.meanlog.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Log Normal",
                    Param = "SCALE",
                    Value = ests$comshp.lnorm.sdlog.ref
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # separate
  if (any(grepl("sep.lnorm.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Log Normal",
                    Param = "INTERCEPT",
                    Value = ests$sep.lnorm.meanlog.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Log Normal",
                    Param = "SCALE",
                    Value = ests$sep.lnorm.sdlog.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Log Normal",
                    Param = "INTERCEPT",
                    Value = ests$sep.lnorm.meanlog.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Log Normal",
                    Param = "SCALE",
                    Value = ests$sep.lnorm.sdlog.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # independent shape (not really supported by STEM so handled as if separate models)
  if (any(grepl("indshp.lnorm.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(
        Model = "Independent shape - Reference", Dist = "Log Normal",
        Param = "INTERCEPT",
        Value = ests$indshp.lnorm.meanlog.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(
        Model = "Independent shape - Reference", Dist = "Log Normal",
        Param = "SCALE",
        Value = ests$indshp.lnorm.sdlog.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    
    rc <- rc_struct %>%
      dplyr::mutate(
        Model = "Independent shape - Intervention", Dist = "Log Normal",
        Param = "INTERCEPT",
        Value = ests$indshp.lnorm.meanlog.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(
        Model = "Independent shape - Intervention", Dist = "Log Normal",
        Param = "SCALE",
        Value = ests$indshp.lnorm.sdlog.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # One arm
  if (any(grepl("onearm.lnorm.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Log Normal",
             Param = "INTERCEPT",
             Value = ests$onearm.lnorm.meanlog.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Log Normal",
             Param = "SCALE",
             Value = ests$onearm.lnorm.sdlog.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  
  #############################
  #### loglogistic models
  #############################
  
  # common shape
  if (any(grepl("comshp.llogis.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Log Logistic",
             Param = "INTERCEPT",
             Value = log(ests$comshp.llogis.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Log Logistic",
             Param = "TX(Intervention)",
             Value = log(ests$comshp.llogis.scale.int) - log(ests$comshp.llogis.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Log Logistic",
             Param = "SCALE",
             Value = 1/ests$comshp.llogis.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # separate
  if (any(grepl("sep.llogis.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Log Logistic",
             Param = "INTERCEPT",
             Value = log(ests$sep.llogis.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Log Logistic",
             Param = "SCALE",
             Value = 1/ests$sep.llogis.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Log Logistic",
             Param = "INTERCEPT",
             Value = log(ests$sep.llogis.scale.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Log Logistic",
             Param = "SCALE",
             Value = 1/ests$sep.llogis.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # independent shape (not really supported by STEM so handled as if separate models)
  if (any(grepl("indshp.llogis.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Log Logistic",
             Param = "INTERCEPT",
             Value = log(ests$indshp.llogis.scale.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Log Logistic",
             Param = "SCALE",
             Value = 1/ests$indshp.llogis.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Log Logistic",
             Param = "INTERCEPT",
             Value = log(ests$indshp.llogis.scale.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Log Logistic",
             Param = "SCALE",
             Value = 1/ests$indshp.llogis.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # one arm
  if (any(grepl("onearm.llogis.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Log Logistic",
             Param = "INTERCEPT",
             Value = log(ests$onearm.llogis.scale.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Log Logistic",
             Param = "SCALE",
             Value = 1/ests$onearm.llogis.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  #############################
  #### gengamma models
  #############################
  
  # common shape
  if (any(grepl("comshp.gengamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Generalized Gamma",
             Param = "INTERCEPT",
             Value = ests$comshp.gengamma.mu.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Generalized Gamma",
             Param = "TX(Intervention)",
             Value = ests$comshp.gengamma.mu.int - ests$comshp.gengamma.mu.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Generalized Gamma",
             Param = "SCALE",
             Value = ests$comshp.gengamma.sigma.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Generalized Gamma",
             Param = "SHAPE",
             Value = ests$comshp.gengamma.Q.ref
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # Separate
  if (any(grepl("sep.gengamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Generalized Gamma",
             Param = "INTERCEPT",
             Value = ests$sep.gengamma.mu.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Generalized Gamma",
             Param = "SCALE",
             Value = ests$sep.gengamma.sigma.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Generalized Gamma",
             Param = "SHAPE",
             Value = ests$sep.gengamma.Q.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Generalized Gamma",
             Param = "INTERCEPT",
             Value = ests$sep.gengamma.mu.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Generalized Gamma",
             Param = "SCALE",
             Value = ests$sep.gengamma.sigma.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Generalized Gamma",
             Param = "SHAPE",
             Value = ests$sep.gengamma.Q.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # independent shape (not really supported by STEM so handled as if separate models)
  if (any(grepl("indshp.gengamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Generalized Gamma",
             Param = "INTERCEPT",
             Value = ests$indshp.gengamma.mu.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Generalized Gamma",
             Param = "SCALE",
             Value = ests$indshp.gengamma.sigma.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Generalized Gamma",
             Param = "SHAPE",
             Value = ests$indshp.gengamma.Q.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Generalized Gamma",
             Param = "INTERCEPT",
             Value = ests$indshp.gengamma.mu.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Generalized Gamma",
             Param = "SCALE",
             Value = ests$indshp.gengamma.sigma.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Generalized Gamma",
             Param = "SHAPE",
             Value = ests$indshp.gengamma.Q.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # one arm
  if (any(grepl("onearm.gengamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Generalized Gamma",
             Param = "INTERCEPT",
             Value = ests$onearm.gengamma.mu.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Generalized Gamma",
             Param = "SCALE",
             Value = ests$onearm.gengamma.sigma.int
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Generalized Gamma",
             Param = "SHAPE",
             Value = ests$onearm.gengamma.Q.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  #############################
  #### gompertz models
  #############################
  
  # common shape
  if (any(grepl("comshp.gompertz.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Gompertz",
             Param = "INTERCEPT",
             Value = -log(ests$comshp.gompertz.rate.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Gompertz",
             Param = "TX(Intervention)",
             Value = -(log(ests$comshp.gompertz.rate.int) - log(ests$comshp.gompertz.rate.ref))
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Gompertz",
             Param = "SCALE",
             Value = ests$comshp.gompertz.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # separate
  if (any(grepl("sep.gompertz.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Gompertz",
             Param = "INTERCEPT",
             Value = -log(ests$sep.gompertz.rate.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Gompertz",
             Param = "SCALE",
             Value = ests$sep.gompertz.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Gompertz",
             Param = "INTERCEPT",
             Value = -log(ests$sep.gompertz.rate.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Gompertz",
             Param = "SCALE",
             Value = ests$sep.gompertz.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # independent shape (not really supported by STEM so handled as if separate models)
  if (any(grepl("indshp.gompertz.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Gompertz",
             Param = "INTERCEPT",
             Value = -log(ests$indshp.gompertz.rate.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Gompertz",
             Param = "SCALE",
             Value = ests$indshp.gompertz.shape.ref
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Gompertz",
             Param = "INTERCEPT",
             Value = -log(ests$indshp.gompertz.rate.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Gompertz",
             Param = "SCALE",
             Value = ests$indshp.gompertz.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # onearm
  if (any(grepl("onearm.gompertz.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Gompertz",
             Param = "INTERCEPT",
             Value = -log(ests$onearm.gompertz.rate.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Gompertz",
             Param = "SCALE",
             Value = ests$onearm.gompertz.shape.int
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  #############################
  #### gamma models
  #############################
  
  # common shape
  if (any(grepl("comshp.gamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Gamma",
             Param = "INTERCEPT",
             Value = -log(ests$comshp.gamma.rate.ref / ests$comshp.gamma.shape.ref )
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Gamma",
             Param = "TX(Intervention)",
             Value = -(log(ests$comshp.gamma.rate.int) - log(ests$comshp.gamma.rate.ref))
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Common shape", Dist = "Gamma",
             Param = "SCALE",
             Value = ests$comshp.gamma.shape.ref^-0.5
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # separate
  if (any(grepl("sep.gamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Gamma",
             Param = "INTERCEPT",
             Value = -log(ests$sep.gamma.rate.ref  / ests$sep.gamma.shape.ref )
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Reference", Dist = "Gamma",
             Param = "SCALE",
             Value = ests$sep.gamma.shape.ref^-0.5
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Gamma",
             Param = "INTERCEPT",
             Value = -log(ests$sep.gamma.rate.int / ests$sep.gamma.shape.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Separate - Intervention", Dist = "Gamma",
             Param = "SCALE",
             Value = ests$sep.gamma.shape.int^-0.5
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # independent shape (not really supported by STEM so handled as if separate models)
  if (any(grepl("indshp.gamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Gamma",
             Param = "INTERCEPT",
             Value = -log(ests$indshp.gamma.rate.ref / ests$indshp.gamma.shape.ref)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Reference", Dist = "Gamma",
             Param = "SCALE",
             Value = ests$indshp.gamma.shape.ref^-0.5
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Gamma",
             Param = "INTERCEPT",
             Value = -log(ests$indshp.gamma.rate.int / ests$indshp.gamma.shape.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "Independent shape - Intervention", Dist = "Gamma",
             Param = "SCALE",
             Value = ests$indshp.gamma.shape.int^-0.5
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  # onearm
  if (any(grepl("onearm.gamma.", names(ests)))){
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Gamma",
             Param = "INTERCEPT",
             Value = -log(ests$onearm.gamma.rate.int / ests$onearm.gamma.shape.int)
      ) %>%
      dplyr::bind_rows(rc)
    
    rc <- rc_struct %>%
      dplyr::mutate(Model = "One arm - Intervention", Dist = "Gamma",
             Param = "SCALE",
             Value = ests$onearm.gamma.shape.int^-0.5
      ) %>%
      dplyr::bind_rows(rc)
  }
  
  
  
  stemParamDF <- rc
  
  ####################################################################################
  # calculate the covariance matrices from the bootstrap samples if exist
  ####################################################################################
  
  rc.cov.df <- tibble::tibble(Model = "", Dist = "", rowParam = "", colParam = "", Value = 1)%>%
    dplyr::filter(1==2)
  
  if (any(stemParamDF$Estimate == "Bootstrap")){
    
    included_models <- stemParamDF %>%
      dplyr::filter(Estimate == "Bootstrap") %>%
      dplyr::transmute(Model, Dist) %>%
      unique()
    
    for (i in 1:nrow(included_models)) {
      # get covariance
      this.model.cov <- stemParamDF %>%
        dplyr::filter(Model == included_models$Model[i],
               Dist == included_models$Dist[i],
               Estimate == "Bootstrap") %>%
        dplyr::transmute(bootid, Param, Value) %>%
        stats::reshape(direction = "wide", idvar = "bootid", timevar = "Param", v.names = "Value") %>%
        dplyr::select(-bootid) %>%
        stats::cov(use = use)
      
      # test code for matrix dimensioning
      # this.model.cov <- matrix(data = c(1,2,3,4,5,6), ncol = 2)
      # rownames(this.model.cov) <- c("a","b","c")
      # colnames(this.model.cov) <- c("x","y")
      
      this.cov.df <- tibble::tibble(
        Model = included_models$Model[i],
        Dist = included_models$Dist[i],
        rowParam = rep(rownames(this.model.cov), times = dim(this.model.cov)[2]),
        colParam = rep(colnames(this.model.cov), each = dim(this.model.cov)[1]),
        Value = as.numeric(this.model.cov)
      )
      
      rc.cov.df <- rc.cov.df %>%
        dplyr::bind_rows(this.cov.df)
      
    }
    
    #tidy up the names
    
    rc.cov.df <- rc.cov.df %>%
      dplyr::mutate(rowParam = gsub("Value.", "", rowParam, fixed=TRUE),
                    colParam = gsub("Value.", "", colParam, fixed=TRUE)
      )
    
  }
  
  ####################################################################################
  # calculate a summary of estimates if exist
  ####################################################################################
  
  rc.est.df <- stemParamDF %>%
    dplyr::filter(Estimate == "Main") %>%
    dplyr::transmute(Model, Dist, Param, Estimate = Value)
  
  
  ####################################################################################
  # get the model summaries
  # the naming of the models do not match other values
  ####################################################################################
  
  if (!is.null(x)){
    
    rc.modsum <- x$model_summary %>%
      dplyr::mutate(AIC_SAS = NaN,
                    BIC_SAS = NaN)
    
    
    for (i in 1:nrow(rc.modsum)){
      this.mdl <- rc.modsum$flexsurvfit[i]
      mdl <- x$models[[this.mdl]]
      
      # convert the AIC & BIC to match SAS formulations
      # difference is due to choice of log(t) vs t as response in likelihood
      # extract the event times from the flexsurv object
      
      sum_log_event_time <- mdl$data$Y[mdl$data$Y[,"status"]==1, "time"] %>%
        log() %>%
        sum()
      
      rc.modsum$AIC_SAS[i] <- stats::AIC(mdl) - 2 * sum_log_event_time
      rc.modsum$BIC_SAS[i] <- stats::BIC(mdl) - 2 * sum_log_event_time
    }
    
  } else {
    rc.modsum <- tibble::tibble(Model = "", Dist = "", Status = "", AIC = NaN, AIC_SAS = NaN, BIC = NaN, BIC_SAS = NaN) %>%
      dplyr::filter(1==2)
  }
  
  ####################################################################################
  # Add factors for easy sort
  ####################################################################################
  
  STEMModel.levels <- c("Common shape", "Separate - Reference", "Separate - Intervention",
                    "Independent shape - Reference", "Independent shape - Intervention",
                    "Independent shape",
                    "One arm - Intervention")
  
  STEMDist.levels <- c("Exponential", "Weibull", "Log Normal", "Gamma", "Generalized Gamma", "Log Logistic", "Gompertz")
  
  STEMParam.levels <-  c("INTERCEPT", "TX(Intervention)", "SCALE", "SHAPE")
  
  stem_param <- rc.est.df %>%
    dplyr::transmute(Model, ModelF = factor(Model, levels = STEMModel.levels, ordered = TRUE),
              Dist, DistF = factor(Dist, levels = STEMDist.levels, ordered = TRUE),
              Param, ParamF = factor(Param, levels = STEMParam.levels, ordered = TRUE),
              Estimate
    ) %>%
    dplyr::arrange(ModelF, DistF, ParamF)
  
  stem_cov <- rc.cov.df %>%
    dplyr::transmute(Model, ModelF = factor(Model, levels = STEMModel.levels, ordered = TRUE),
              Dist, DistF = factor(Dist, levels = STEMDist.levels, ordered = TRUE),
              rowParam, colParam,
              rowParamF = factor(rowParam, levels = STEMParam.levels, ordered = TRUE),
              colParamF = factor(colParam, levels = STEMParam.levels, ordered = TRUE),
              rowNum = as.numeric(rowParamF),
              colNum = as.numeric(colParamF),
              CovEst = Value
    ) %>%
    dplyr::arrange(ModelF, DistF, rowParamF, colParamF)
  
  stem_modsum <- rc.modsum %>%
    dplyr::transmute(Model, ModelF = factor(Model, levels = STEMModel.levels, ordered = TRUE),
              Dist, DistF = factor(Dist, levels = STEMDist.levels, ordered = TRUE),
              Status,
              AIC, AIC_SAS,
              BIC, BIC_SAS) %>%
    dplyr::filter(!is.na(DistF))
  
  
  stem_boot <- stemParamDF %>%
    dplyr::filter(Estimate == "Bootstrap") %>%
    dplyr::transmute(Model, ModelF = factor(Model, levels = STEMModel.levels, ordered = TRUE),
              Dist, DistF = factor(Dist, levels = STEMDist.levels, ordered = TRUE),
              BootSample = bootid,
              Param, ParamF = factor(Param, levels = STEMParam.levels, ordered = TRUE),
              Estimate = Value) %>%
    dplyr::arrange(ModelF, DistF, BootSample, ParamF)
  
  ####################################################################################
  # return the data
  ####################################################################################
  
  stem_return <- list(stem_param = stem_param,
                      stem_cov = stem_cov,
                      stem_modsum = stem_modsum,
                      stem_boot = stem_boot)
  
  
  return(stem_return)
}