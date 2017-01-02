
#' simulate_pulsets
#'
#'
#'   Simulates hormone concentration per minute for 1 subject.
#'
#'   Pulse simulation function adapted from Ken Horton by reducing to single
#'   subject, restructuring, and commenting. Ken's originally titled "subject
#'   model sim aim1 mua.R". You have two options for passing arguments: 1)
#'   specify arguments individually or 2) pass a named vector of parameters to
#'   parms_list (see function code for naming scheme).
#'
#'
#' @param parms_list Named vector of parameters: T, t, vare, ipimean, ipivar,
#' ipimin, mua, sda, muw, and sdw are required.  Either muh and varh or
#' constant.halflife must be specified.  Same for mub and varb or
#' constant_baseline.
#' @export
#' @keywords pulse simulation
#' simulate_pulsets()
simulate_pulsets <- function(parms_list) {
    

    # Handle vector of parameters
    if (!is.null(parms_list)){
        T       = parms_list[["T"]] # num. samples
        t       = parms_list[["t"]] # sampling interval
        vare    = parms_list[["vare"]]
        ipimean = parms_list[["ipimean"]]
        ipivar  = parms_list[["ipivar"]]
        ipimin  = parms_list[["ipimin"]]
        mua     = parms_list[["mua"]]
        sda    = parms_list[["sda"]]
        muw     = parms_list[["muw"]]
        sdw    = parms_list[["sdw"]]
        if (!is.null(parms_list[["muh"]]) & !is.na(parms_list[["muh"]])) {
            muh     = parms_list[["muh"]]
            varh    = parms_list[["varh"]]
            mub     = parms_list[["mub"]]
            varb    = parms_list[["varb"]]
        } else {
            constant_halflife = parms_list[["constant_halflife"]]
            constant_baseline = parms_list[["constant_baseline"]]
        }
    } else {
      stop("parms_list argument NULL or invalid")

    }


    
    #---------------------------------------
    # Calculate parameters from input params
    #
    #   Using a renewal process, define by interpulsatile interval and variance
    #   then convert gamma parameters
    # 
    # mean = alpha / beta
    # var = alpha / beta^2
    #---------------------------------------
    gammamean <- ipimean - ipimin
    alpha     <- gammamean * gammamean / ipivar
    beta      <- gammamean / ipivar

    
    
    #---------------------------------------
    # Create data structures for generated parameters
    #---------------------------------------
    B      <- 0     # Baseline concentration
    H      <- 0     # Decay rate (in terms of half-life)
    tau    <- rep(0, 25) # 
    ysim   <- vector(mode = "numeric", length = length(seq(10,(T*t),t)))
    
    
    #---------------------------------------
    # Helper functions for drawing hormone 
    #   concentration (w/o error)
    #---------------------------------------
    # Normal CDF 
    erfFn <- function(x){
    	    y <- 2 * stats::pnorm(x * sqrt(2), 0, 1) - 1
    	    y
    	}
    
    # Model concentation over time given pulse parameters  
    # t = seq(10, 1440, by = 10)
    # b = 0 baseline, added later
    # a = pulse masses
    # tau1 = pulse locations
    # lam = decay rate (transformed from half-life)
    # s2p = pulse widths
    # NOTE: based on above and eqn, it looks to me like widths, half-life are on
    # the minutes scale, not sampling units
    meanI <- function(t, b, a, tau1, lam, s2p){
    	    m <- b + (a / 2) * exp((tau1 - t) * lam + 0.5 * lam^2 * s2p) * 
                 (1 + erfFn((t - (tau1 + lam * s2p)) / sqrt(2 * s2p)))
    	    m
    	}
    
    
    #---------------------------------------
    # Simulate subject-specific parameters
    #---------------------------------------
    # Get subject specific parameters
    if (is.null(constant_baseline)) {
        while(B <= 0){
            B <- stats::rnorm(1, mub, sqrt(varb))
        }
    } else {
        B <- constant_baseline
    }
    
    # Sample for decay rate from normal given muh, varh 
    #   H is half-life, H=ln(2)/lambda_x, where lambda_x is the decay constant
    if (is.null(constant_halflife)) {
        while(H <= 8){
            H <- stats::rnorm(1, muh, sqrt(varh))
        }
    } else {
        # set constant half-life if 
        H <- constant_halflife
    }

    
    #---------------------------------------
    # Simulate pulse locations
    #---------------------------------------
    # pulse location vector - max 25 pulses per day
    tau <- rep(0, 25)
    # generate initial pulse location
    tau[1] <- t * (stats::rgamma(1, alpha, beta) - 10)
    
    # generate i+1, i+2,... pulse locations
    i <- 1
    while (tau[i] < (T * t)){
    	  i      <- i + 1
    	  tmp    <- ipimin + stats::rgamma(1, alpha, beta)
    	  tau[i] <- tau[i-1] + (t * tmp)
    }
    
    # reduce pulse location vector to values within time range (<0, 1440)
    tau <- subset(tau, tau < (T * t))
    tau <- subset(tau, tau != 0)
    np  <- length(tau) # No. of pulses - length of vector of pulse locations
                       #    w/in time range
    
    
    #---------------------------------------
    # Generate pulse-specific parameters 
    #---------------------------------------
    # storage vectors for pulse-specific parms
    A     <- rep(0, np)             # pulse mass
    s2p   <- rep(0, np)             # pulse width
    taxis <- seq(10, (T*t), t)      # time axis 144*10
    ytmp  <- rep(0, length(taxis))  # hormone concentration
    
    # Generate parameters
    for (i in 1:np) {
        #while (A[i] <= 0) {
            A[i] <- exp(stats::rnorm(1, mua, sda))
        #}
        #while (s2p[i] <= 0) {
            s2p[i] <- exp(stats::rnorm(1, muw, sdw))
        #}
        # Truncated T /gamma normal mixture
        #kappa1 = rgamma(1, 2, 2)
        #kappa2 = rgamma(1, 2, 2)
        #tvar   = sda^2 / kappa1
        #t2var  = sdw^2 / kappa2

        #while (A[i] < 0.25) {
        #    A[i]   <- rnorm(1, mua, sqrt(tvar))
        #}

        #while (s2p[i] < 0.5) {
    	  #  s2p[i] <- rnorm(1, muw, sqrt(t2var))
        #}
    }

    # Draw mean concentration
    for (i in 1:np) {
       ytmp   <- ytmp + meanI(seq(10, (T*t), t), 0, A[i],
                              tau[i], log(2) / H, s2p[i])
    }
    
    # Create vector of pulse parameters 
    allpulseparms <- cbind("pulse_no" = seq(1, np),
                           "mass"     = A,
                           "width"    = s2p,
                           "location" = tau)
    allpulseparms <- tibble::as_data_frame(allpulseparms)
    
    
    #---------------------------------------
    # Combine into final simulated data
    #---------------------------------------
    # Add baseline to generated hormone concentration
    ysim <- ytmp + B
    
    # Two possible formulations of error variance
    # first is using coefficient of variation, second is using explicit variance
    #error_var <- mean(ysim) * vare
    error_var <- vare
    # Generate and add random noise on log scale
    errors    <- stats::rnorm((length(taxis)), 0, sqrt(error_var))
    ysimerror <- ysim * exp(errors)
    ysim_df   <- cbind("observation"   = 1:length(taxis),
                       "time"          = taxis,
                       "concentration" = ysimerror) 
    ysim_df   <- tibble::as_data_frame(ysim_df)

    return(list("pulse_data" = ysim_df, "pulse_parms" = allpulseparms))

}




#-------------------------------------------------------------------------------
# End of file  
#-------------------------------------------------------------------------------
