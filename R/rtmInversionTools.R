
################################################################################
## Inversion helper functions
## Marco Visser, Gamboa, Panama, February 2014
## updated for Bayesian tools Dec2019
############################################################
## Updates: Marco Visser
## 2014-2019, many locations.
################################################################################


################################################################################
## adminstrative functions
################################################################################

## moving window gelman diagnostic
## to estimate the burnin time automatically 
gelmanMW <- function(data=NULL, window=0.10, step=2,
                     thres=1.1,stepIt=1){

    if(is.null(data)) stop("error data = NULL")

    total <- nrow(data[[1]]$Z)
    window <- floor(total*0.1)
    step <- round(window/step)
    spots <- seq(from=1, to=c(total-window), by=step)
    result <- numeric(length = length(spots))
    
    for(i in 1:length(spots)){

        wind <- lapply(data,function(X)
            coda::as.mcmc(X$Z[spots[i]:(spots[i]+window),]))
        wind <- coda::as.mcmc.list(wind)
        result[i] <- coda::gelman.diag(wind)$mpsrf
    }

    burn <- spots[head(which(result<thres),1)]
    ext <- 0

    ## Convergence was achieved in last stepIt
    ## so returning 50% stepIf
    if(length(burn)==0) {
        ext <- 1
        burn <- 1+total-(stepIt/2)
        warning("gelmanMW: Burnin could not be determined.")
    }
        
    return(list(mpsrf=result,
                burnin=burn,
                nsamples=total,
                exit=ext
                )
              )
}

## clear a cat line
catclr <- function() {

    cat("\33[2K")

}

## Neat concatenate and Print repeatedly
## Takes a given message
## and makes sure it is not too long
## for a repeated call
neatCat <- function(string="Lores Ipsum",nmx=65) {

    if(nmx<30) nmx <- 30
    
    if(exists(".globalNmx")) nmx <- .globalNmx

    if(length(string)>1) paste(string,collapse=" ")
    
    Ncs <- nchar(string)
    if(Ncs>nmx){
        string1  <- substr(string,start=1,stop=(nmx-15))
        string2  <- substr(string,start=Ncs-15,stop=Ncs)
        string <- paste(string1," [...] ", string2)
    }

    cat(string)
    
}
    
## make a fit log
## @param steps number of steps in fit sequence
## @param step current step 
## @param error any errors
## @param fit fit statistic
updateLog <- function(steps,
                      step=NULL,
                      errors=NULL,
                      fit=NULL,
                      clear=FALSE){

    if(clear) {
        rm(.fitLog,envir=.GlobalEnv)
        return(invisible(1))
    } else {
    
        if(!exists(".fitLog")){
            .fitLog <<- vector("list",steps)
        }

        if(!is.null(step)){
            .fitLog[[step]] <<- c(0,NA)
            names(.fitLog[[step]])<<- c("e","fit")
        }
        ## find current step
        if(is.null(step))
            step <- tail(which(!sapply(.fitLog,is.null)),1)
        if(!is.null(errors)){

           .fitLog[[step]]["e"] -> nE 
           .fitLog[[step]]["e"] <<- nE+errors

        }
        if(!is.null(fit)){
            .fitLog[[step]]["fit"] <<- fit
            
        }
        
        
    }
    
}

## function to get the current error rate
getErr <- function(){


    if(!exists(".fitLog")){
        stop("NO LOG EXISTS!")
        }

    step <- tail(which(!sapply(.fitLog,is.null)),1)

    return(.fitLog[[step]]["e"]) 


}



## build parameter data.frame
## builds parameters to feed to package hsdar prospect
## @param prms parameter vector from MCMC
## @param n names of all parameters
## @parm fixed a named vector of fixed values
buildParList <- function(prms,nms=c("N", "Cab","Car", "Cw", "Cm", "Cbrown"),
                         fixed=NULL){

    if(is.null(fixed)){
    param <- prms[1:length(nms)]
    names(param) <- nms
    }
    
    if(!is.null(fixed)){

        if(is.null(names(prms))) stop("parameters must be named!")
        
        prms <- prms[names(prms)%in%nms]
        param <- numeric(length(nms))
        param[!(nms%in%names(fixed))] <- prms
        param[nms%in%names(fixed)] <- fixed
        names(param) <- nms
    }
    
    return(as.data.frame(t(param)))
}


## S3- methods for getting defaults for models
##       - N     = leaf structure parameter
##       - Cab   = chlorophyll a+b content 
##       - Car   = carotenoids content 
##       - Cbrown= brown pigments content in arbitrary units
##       - Cw    = equivalent water thickness
##       - Cm    = dry matter content 
##
getDefaults<-function(modelname=NULL, ...){
    class(modelname) <- modelname
    defaults(modelname)
}

## defualts function
defaults <- function(x, ...){
    UseMethod("defaults",x)
}

## get some default values for the prospect model 
##       - N     = leaf structure parameter
##       - Cab   = chlorophyll a+b content 
##       - Car   = carotenoids content 
##       - Cbrown= brown pigments content in arbitrary units
##       - Cw    = equivalent water thickness
##       - Cm    = dry matter content 
##
defaults.prospect<-function(version="5B"){
    
    
    ## typical values are the rice values
    typical = c(2.698,70.8,20, 0.0117,0.009327,0.001,0.1)
    ## N, Cab,Car, Cw, Cm, Cbrown, resid
    omega <- c(5,150,140,.4,0.45,1,4) 
    alpha <- c(1,10,1,1e-5,1e-5,0,1e-6) 
    
    names(typical)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown","sigma")
    names(omega)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown","sigma")
    names(alpha)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown","sigma")
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega
    
    if(version=="5") def <- def[-6]
    
    return(def)
}

## get some default values for the prosail model 
##       - N     = leaf structure parameter
##       - Cab   = chlorophyll a+b content 
##       - Car   = carotenoids content 
##       - Cbrown= brown pigments content in arbitrary units
##       - Cw    = equivalent water thickness
##       - Cm    = dry matter content psoil: Dry/Wet soil factor
##       - LAI   =  Leaf area index
##       - TypeLidf= Type of leaf angle distribution. 
##       - lidfa = Leaf angle distribution parameter a. 
##       - lidfb = Leaf angle distribution parameter b.
##       - hspot = Hotspot parameter
##       - tts   = Solar zenith angle (degrees)
##       - tto   = Observer zenith angle (degrees)
##       - psi   = Relative azimuth angle (degrees)
defaults.prosail<-function(version="5B"){

    nms <- c("N", "Cab","Car", "Cw", "Cm", "Cbrown",
             "psoil","LAI","TypeLidf","lidfa","lidfb",
             "hspot","tts","tto","psi",
             "sigma")
    
    ## typical values are the rice values
    typical = c(2.698,70.8,20, 0.0117,0.009327,0.001,
                0,4.77,1,-0.35,-0.15,
                0.01,30,10,0,
                0.1)
    
    ## N, Cab,Car, Cw, Cm, Cbrown, resid
    omega <- c(5,150,140,.4,0.45,1,
               1,50,1,1,1,
               1,90,90,180,
               4) 
    alpha <- c(1,0,0,1e-5,1e-5,0,
               0,0.1,1,-1,-1,
               1e-6,0,0,0,        
               1e-8) 
    
    names(typical)<-nms
    names(omega)<-nms
    names(alpha)<-nms
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega
    
    if(version=="5") def <- def[-6,]
    
    return(def)

}

################################################################################
## Likelihood functions
################################################################################

##' Normal difference log-likelihood generator
## Disclaimer: this function works by assigning objects
## in the normDiffLogLik environment, the environment
## is then automatically refferenced
## by the "out function". Note that this is a
## performance poor way of programming functions in R
## (needless copies in memory of observed for instance
## whenever obs is larger this is inefficient)
## 
## A better way would be to open a new environment X
## and save all needed objects in X once
## to which all function can be pointed
normDiffLogLike <- function(nparams,
                            model,
                            observed,
                        penalty=-1e12,
                        verbose = TRUE, ...) {
    
    stopifnot(nparams >= 1, nparams %% 1 == 0,
              is.function(model),
              is.numeric(observed))
    
    n_obs <- length(observed)
    out <- function(x) {

        ## prep model parameters
        params <- x[seq_len(nparams)]
        ## extract residual error
        rsd <- x[nparams]
        
        mod <- model(params, ...)
        
        if (any(c(is.na(mod),is.infinite(mod)))) {
            updateLog(error=1)
            
            if (verbose) {
                catclr()
                neatCat(paste0("\r Flagging event # ",getErr(),
                    "- Errors ecountered in model output (Inf,NA,NaN). ",
                    "Returning penalty loglike = "
                   ,penalty,"\t\t"))
                             }
        
            return(penalty)
        }
        
        err <- mod - observed
        SS <- sum(err * err)
        sigma2 <- as.numeric(rsd * rsd)
        
        ll <- -0.5 * (n_obs * log(sigma2) + SS / sigma2)
    
        if (is.na(ll)|is.infinite(ll)) {
            updateLog(error=1)

        if (verbose) {
            catclr()
            neatCat(paste0("\r Errors encountered in likelihood (NA,Inf). Penalty loglike = ",
                           penalty,"\n",
                "Total error (residual sd): ", SS,"(", rsd,")","\n",
                "LL  = ",ll,"\t\t"))
            
        }
        return(penalty)
    }
    return(ll)
}
return(out)
} 


########################## PROIR functions #########################
#' Default prior parameters for PROSPECT models
#' @param priorparams list of length 2 with named vectors of the first and second
#' parameter for the prior ditributions 
#' @param sdInflate Standard deviation multiplier (default = 1)
#' Use with caution. Note for a lognormal values >1 deflate variance,
#' <1 increase and that changing the sd influences the log normal mean.
#' It is far superior to directly change the sd at the prior definition!
#' @export
priorDefaultsProspect <- function(priorparams=NULL,sdInflate=1) {
    
    if(is.null(priorparams)){
        
        priormean  <- c(N = -0.62, Cab = 3.53, Car = 2.12, Cw = -4.56,
                        Cm = -5.3,Cbrown=-4.6,sigma=-3.68)
        priorsd    <- c(N = .48, Cab = 0.75, Car = 0.59, Cw = 0.42,
                        Cm = 0.48,Cbrown=0.42,sigma=0.25)*sdInflate

        priorname <- rep("dlnorm",length(priormean))
        names(priorname) <- names(priormean)

        
    } else {

        priormean <- priorparams[[1]]    
        priorsd <- priorparams[[2]]    
        priorname <- priorparams[[3]]
        
    }
    
    return(list(mu = priormean, sigma = priorsd,prior=priorname))
} 


#' Default PROSPECT 5 prior density function
#' 
#' @details Assumes lognormal distribution for all parameters. NOTE that prior 
#' on N is shifted by 1.
#' @param pmu Lognormal mu parameter
#' @param psigma Lognormal sigma parameter
#' @export
buildPriorDenProspect <- function() {

if(!exists("dpriors")) {dpriors<<-priorDefaultsProspect()}

  prior <- function(params) {
    if (is.null(names(params))) {
      warning("Parameters are not named.", "\n", 
              "Assuming N Cab (Car) (Cbrown) Cw Cm for priors")
      params[1] <- params[1] - 1
    } else {
      params["N"] <- params["N"] - 1
    }
    priors <- dlnorm(params, meanlog=dpriors$mu,
                     sdlog=dpriors$sigma, log = TRUE)
    return(sum(priors))
  }
  return(prior)
} # priorfunc.prospect


#' Default lognormal prior sampling function for PROSPECT 5/5B
#' 
#' @details Assumes lognormal distribution for all parameters. NOTE that prior 
#' on N is shifted by 1.
#' @param pmu Lognormal mu parameter
#' @param psigma Lognormal sigma parameter
#' @export
buildPriorSampProspect <- function() {
	if (is.null(names(dpriors$mu))) {
      stop("Parameters are not named \n")
      } 
	sampler <- function(n=1) {
    
   samples <- sapply(1:length(dpriors$mu),function(X)
		rlnorm(n, meanlog=dpriors$mu[X],
			 sdlog=dpriors$sigma[X]))

	if(n==1){
      names(samples)<-names(dpriors$mu)
      samples["N"] <- samples["N"]+1
          } else{
      colnames(samples)<-names(dpriors$mu)
      samples[,"N"] <- samples[,"N"]+1
     
      }
    return(samples)
  }
  return(sampler)
}

#' Default prior parameters for PROSAIL models
#'
#' @inheritParams priorDefaultsProspect
#' 
#' @export
priorDefaultsProsail <- function(priorparams=NULL, sdInflate = 1) {

    
    if(is.null(priorparams)){
        
        priormean  <- c(N = -0.62, Cab = 3.53, Car = 2.12, Cw = -4.56,
                        Cm = -5.3,Cbrown=-4.6,
                        psoil=-3,LAI=.9,lidfa=-0.35,lidfb=-0.15,
                        hspot=-4.6,
                        tts=3,
                        tto=2,
                        psi=2,
                        sigma=-3.68)

        priorsd    <- c(N = .48, Cab = 0.75, Car = 0.59, Cw = 0.42,
                        Cm = 0.48,Cbrown=0.42,
                        psoil=0.01,LAI=.65,lidfa=0.65,lidfb=0.65,
                        hspot=2.5,
                        tts=.65,
                        tto=1,
                        psi=1,
                        sigma=.55) * sdInflate

    
        prior <- rep("dlnorm",length(priormean))
        prior[9:10] <- "dnorm"
        names(prior) <- names(priormean)
        
    
    } else {

    priormean <- priorparams[[1]]
    priorsd <- priorparams[[2]]
    prior <- priorparams[[3]]    

    }
    
    return(list(mu = priormean, sigma = priorsd,
                prior=prior))
} 


#' Default PROSAIL 5 prior density function
#' 
#' @details Assumes lognormal distribution for non negative parameters.
#' Otherwise a normal distribution is used. 
#' NOTE that prior on N is shifted by 1.
#' 
#' @param parSel which parameters are optimized?
#' @export
buildPriorDenProsail <- function() {
    
    if(!exists("dpriors")) {dpriors<<-priorDefaultsProsail()}
    
    prior <- function(params) {
    if (is.null(names(params))) {
        stop("Parameters are not named") 
    } else {
        params["N"] <- params["N"] - 1
    }
    
    inc <- which(dpriors$prior=="dlnorm")
    priorsLn <- dlnorm(params[inc],
                       meanlog=dpriors$mu[inc],
                       sdlog=dpriors$sigma[inc], log = TRUE)
    priorsN <- dnorm(params[-inc],
                     mean=dpriors$mu[-inc],
                     sd=dpriors$sigma[-inc], log = TRUE)

    return(sum(c(priorsN,priorsLn)))
    }
    return(prior)
} # priorfunc.prosail


#' Default lognormal prior sampling function for PROSAIL 5/5B
#' 
#' @details Assumes lognormal distribution for non negative parameters.
#' Otherwise a normal distribution is used.
#' @param parSel which parameters are optimized?
#' 
#' @export
buildPriorSampProsail <- function() {
    
    if (is.null(names(dpriors$mu))) {
        stop("Parameters are not named \n")
    } 
    sampler <- function(n=1) {
        
        samples <- sapply(1:length(dpriors$mu),function(X)
            ifelse(dpriors$prior[X]=="dlnorm",
                   rlnorm(n, meanlog=dpriors$mu[X],
                          sdlog=dpriors$sigma[X]),
                   rnorm(n, mean=dpriors$mu[X],
                         sd=dpriors$sigma[X])))
        
        
	if(n==1){
            names(samples)<-names(dpriors$mu)
            samples["N"] <- samples["N"]+1
        } else{
            colnames(samples)<-names(dpriors$mu)
            samples[,"N"] <- samples[,"N"]+1
            
        }
        return(samples)
    }
    return(sampler)
}



###############################################################################
## Functions for fitting 
###############################################################################

#' check for convergence trend
#' function to help decide weather to continue
#' a mcmcm run
#' @param fitStats a time series of fit stats
#' @requires Kendall
convTrend <- function(fitStats,thres=10){

    ## discard earler observations
    discard <- floor(thres/4)
    
    N <- length(fitStats)
    if(N<thres) return(TRUE)
    
    MK <- Kendall::SeasonalMannKendall(ts(fitStats[-(1:discard)]))

    if(N>thres*2) {

        if(MK$tau>=0) return(FALSE)
        if(MK$tau<0&MK$sl<0.1) {
            return(TRUE)
        } else {
            return(FALSE)
                }
        
    } else if(N<=thres*2) {

        if(MK$tau>=0) return(FALSE)
        if(MK$tau<0&MK$sl<0.3) {
            return(TRUE)
        } else {
            return(FALSE)
        }
        
    }
    
}


##' check convergence 
checkConverge <- function(samples,
                          threshold = 1.1,
                          use_CI = TRUE,
                          use_mpsrf = TRUE,
                          verbose=TRUE) {
  i <- ifelse(use_CI, 2, 1)
  gelman <- try(BayesianTools::gelmanDiagnostics(samples))
  if (class(gelman) == 'try-error') {
  if(verbose)  neatCat('\r Error in gelman diagnostic. Assuming no convergence..')
    return(FALSE)
  }
  if (use_mpsrf) {
    gelman_vec <- c(gelman$psrf[,i], mpsrf = gelman$mpsrf)
  } else {
    gelman_vec <- gelman$psrf[,i]
  }

  exceeds <- gelman_vec > threshold

  if (any(exceeds)) {
    exceeds_vec <- gelman_vec[exceeds]
    exceeds_char <- sprintf('%s: %.2f', names(exceeds_vec), exceeds_vec)
    exceeds_str <- paste(exceeds_char, collapse = '; ')
    
    if(nchar(exceeds_str)>100) {
        exceeds_str <- paste0(length(exceeds_vec)-sum(use_mpsrf),
                              " Parameters, mpsrf: ",
                              tail(gelman_vec))
    }
    
    if(verbose) {
        catclr()
        neatCat(paste0('\r Parameters exceding threshold: ', exceeds_str))
        }
    return(FALSE)
  } else {
    return(TRUE)
  }
}


#' Fit and Converge a Bayesian Tools model 
#' 
#' \code{FitConv} - runs a model and checks for convergence
## using gelman stats. It then returns standard output from runMCMC
#' 
#' @param btSetup A setup object from Bayesian tools
#' @param btSampler the desired MCMC algorithm 
#' @param btSettings A settings object from Bayesian tools
#' @param gelmanthr Gelman threshold values to except fit
#' @param steps number of checks to conduct while sampling
#' stops every btSettings$interation/steps to see
#' if MCMC chains have converged
#' @param maxSteps maximum number of steps to take before
#' stopping (maximum number of possible iteration is therefore
#' btSettings$interation/steps * maxSteps). 
#' @param attemps if program determines that no
#' gelman diagnostic improvements are being
#' made, how many simes to restart?
#' Basically, if model fails to converge, then this
#' sets the number of attempts at which the model will be restarted
#' (with new initial values).
#' @param autoburnin attempt find the burnin period
#'  and return samples after this period
#' @param errTol how many errors in the model/likelihood to
#' tolerate before stopping? High error rate could
#' indicate chain is in a unfeasible region
#' @author Marco D. Visser
#' 
#' @export
fitConv <- function(btSetup = bayesianSetup,
                    btSampler = "DEzs",
                    btSettings = settings,
                    gelmanthr=1.1,
                    steps=10,
                    maxSteps=200,
                    attempts=3,
                    autoburnin=FALSE,
                    verbose=TRUE,
                    convThres=10,
                    errTol=1/3,
                    plotTS=FALSE
                    ){

    maxIt <- btSettings$iterations
    stepIt <- ceiling(maxIt/steps)
    maxIt <- maxIt*maxSteps

    
    btSettings$iterations <- stepIt

    btSettings$burnin <- stepIt/2 

    attempt <- 1
    GelmanStat <- vector("list",steps)
    
    CurIt <- 1
    for(j in 1:attempts){

        if(verbose) {
            catclr()
            neatCat(paste0("\r ****  Starting fit attempt # ",j,
                           " **** \t\t\t\t \n"))
           }
        
        ## prepare log
        if(exists(".fitLog")) updateLog(clear=TRUE)
        updateLog(steps=maxSteps,step=1)

        ## run first step
        out <- runMCMC(bayesianSetup = btSetup,
                       sampler=btSampler,
                       settings = btSettings)

        updateLog(fit=gelmanDiagnostics(out)$mpsrf)

        for(i in 2:maxSteps){

            updateLog(step=i)

            if(checkConverge(out,threshold=gelmanthr,
                             verbose=verbose)){

                ## burnin would be CurIt
                if(autoburnin){
                    
                    B <- gelmanMW(out,window=0.1,step=5,
                                  thres=gelmanthr,stepIt=stepIt)
                    
                    samples <- lapply(out,
                                      function(X)
                                          coda::as.mcmc(X$Z[B$burnin:B$nsamples,])
                                      )
                    
                } else{

                    return(out)
                }
                
                if(verbose) cat("\n")
                if(exists(".fitLog")) updateLog(clear=TRUE)
                return(coda::as.mcmc.list(samples))
            }
            
            if(CurIt>=maxIt) {break}
            if(CurIt<maxIt){

                             
                ## run though another iteration
                out <- runMCMC(out)
                updateLog(fit=gelmanDiagnostics(out)$mpsrf)

                ## check error rate
                if((getErr()/stepIt)>errTol) {
                    if(verbose) {
                        catclr()
                        neatCat("\r - !! Error rate > tolerance .. stopping !! - \n")
                        }
                    break
                }

                 
                gelmts <- ts(unlist(sapply(.fitLog,
                                           function(X) X["fit"])))
                if(plotTS&i%%3==0) {
                    par(mfrow=c(2,1))
                    plot(gelmts,main="fit",ylab="mpfrs")
                    abline(h=gelmanthr,col="red")
                    errs <- ts(unlist(sapply(.fitLog,
                                               function(X) X["e"])))
                    
                    plot(errs/stepIt,main="error rate")
                    abline(h=errTol,col="red")
                    }

                if(verbose) {
                    catclr()
                    neatCat(paste0("\r  Fit attemp: ",j,
                                   " - Current mpfrs: ",
                                   tail(gelmts),
                                   "\t\t"))
                }

                ## check if sequence is converging
                if(!convTrend(gelmts,thres=convThres)){
                    catclr()
                 if(verbose) cat("\r - !! failed convergence test .. stopping !! - \n")
                    break
                }
                
                CurIt <- stepIt*i
                
            }
        }


        if(verbose) cat("Convergence not achieved in",j,"attemp(s) - restarting  \n")
    }
    
    if(verbose) cat("Convergence not achieved in",attempts,"attemp(s) - stopping \n")
    if(exists(".fitLog")) updateLog(clear=TRUE)

    return(out)
}


################################################################################
## parallel functions

defaultControls <- function(){

    controls <- data.frame(
        gelmanthr=1.1,
        steps=20,
        maxSteps=60,
        attemps=5,
        errTol=0.25,
        convThres=10,
        autoburnin=TRUE,
        plotTS=FALSE,
        stoponerror=FALSE,
        sinkfile=TRUE,
        verbose=FALSE
            )

    return(controls)
}

#' Function to run a single parallel execution of a fit function for all rows
#' of a dataset (refmat)
#'
#' @param refmat matrix with reflection values and obervations in rows
#' and wavelengths in columns
#' @param refPars parameter information as return by e.g. getDefaults("prosail")
#' @param modelFunc a model formulation that returns reflection values that
#' correspond to refmat
#' @param prior a BayesianTools prior
#' @param settings a BayesianTools settings object
#' @param start if starting somewhere other than 1, this feeds information
#' on where the parallel worker is to a log file
#' @param chainName unique name for each call to fitSection
#' help identify output in log files and fail safes.
#'  If NULL a random string is created. 
#' @param controls list of settings to send to fitConv
#' @export
fitSection <- function(refmat,refPars,modelFunc,prior,settings,
                       strt=1,chainName=NULL,controls=NULL){

    if(is.null(controls)){
        controls <- defaultControls()
    } 
    
    resultsMu  <- array(dim=c(nrow(refmat),nrow(refPars)))
    colnames(resultsMu) <- rownames(refPars)
    resultsL  <- resultsU <- resultsMu

    ## create unique output code to prevent overwrite
    
    if(is.null(chainName)){
        Ucode <- paste(sample(letters,6,replace=TRUE),collapse="")
    } else {
        Ucode <- chainName
    }
    

    if(controls$sinkfile){
        logF <- file(paste0("./objects/messagedumpCore_",Ucode,
                            strt,".txt"), open = "wt")
        
        sink(file=logF,type = "message",append=TRUE)
    }
    
    for(i in 1:nrow(refmat)){
        
        ## build a likelihood function      
        likelihood <- normDiffLogLike(nrow(refPars),model=modelFunc,
                                      observed=refmat[i,],
                                      penalty=-1e12,
                                      verbose = FALSE)
        
        bayesianSetup <- createBayesianSetup(likelihood, prior,
                                             names = rownames(refPars))
        T0 <- Sys.time()
        cat(paste(Ucode," Start: ", (strt-1)+i," ",T0,"\n"),
            append=TRUE,file="./objects/fitlog.txt")

        out <- try(fitConv(btSetup = bayesianSetup,
                           btSampler = "DEzs",
                           btSettings = settings,
                           gelmanthr=controls$gelmanthr,
                           steps=controls$steps,
                           maxSteps=controls$maxSteps,
                           convThres=controls$convThres,
                           attempts=controls$attemps,
                           errTol=controls$errTol,
                           autoburnin=controls$autoburnin,
                           plotTS=controls$plotTS,
                           verbose=controls$verbose),
                   silent=TRUE)

        T1 <- Sys.time()
        dT <- T1-T0

        if(class(out)=="try-error") {
            cat(paste(Ucode," ERROR in fitConv: ", (strt-1)+i," ",T1,"\n"),
                append=TRUE,
                file="./objects/fitlog.txt")

          if(controls$stoponerror) return(out)

            
        } else {
            
            cat(paste(Ucode," Finished: ", (strt-1)+i," ",T1,
                      " run time: ",round(as.numeric(dT[[1]]),2)," ",
                      attributes(dT)$unit,
                      "\n"),append=TRUE,
                file="./objects/fitlog.txt")
            
            if(class(out)=="mcmc.list"){

                full <- summary(out)
                resultsMu[i,] <- full$statistics[,1]
                resultsL[i,] <- full$quantiles[,1]
                resultsU[i,] <- full$quantiles[,5]

                saveRDS(list(resultsMu,resultsL,resultsU),
                        paste0("./objects/resultsdumpCore_",Ucode,strt,".rds"))
                
            } else {

                cat(paste(Ucode," Finished: ", (strt-1)+i,": non mcmc.list output:",
                          class(out),
                      "\n"),append=TRUE,
                file="./objects/fitlog.txt")

                saveRDS(out,
                        paste0("./objects/resultsdump_",Ucode,i,".rds"))


            }
            
        }
        
        gc()
    }

    if(controls$sinkfile) sink()

    ##for(i in seq_len(sink.number("message")))
    #close(logF)
    
    return(list(resultsMu,resultsL,resultsU))
    
}

################################################################################
## functions to help interpret and plot the PROSAIl RT model 
################################################################################
## Gives probabilty density of leaves at different angles
## from flat (horizonal) 0 to fully erect 90.
dladgen <- function(a,b){

    litab <- c(5.,15.,25.,35.,45.,55.,65.,75.,81.,83.,85.,87.,89.)
    #litab <- seq(1,89,1)
    freq <- numeric(length(litab))
                    
    for(i1 in 1:8){
        t <- i1*10
        freq[i1] <- dcum(a,b,t);
    }
    
    for(i2 in 9:12){
        t <- 80.+(i2-8)*2.
    freq[i2] <- dcum(a,b,t)
    }


    freq[13] <- 1
    for(i in 13:2){
        freq[i] <- freq[i]-freq[i-1];
    }
    return(list(freq,litab))
}

dcum <- function(a,b,t){

    rd <- pi/180
    if (a>=1){
        f <- 1-cos(rd*t)
    }
    eps <- 1e-8
    delx <- 1
    x <- 2*rd*t
    p <- x
	while(delx >= eps) {
        y <- a*sin(x)+.5*b*sin(2.*x)
        dx <- .5*(y-x+p)
        x <- x+dx
        delx <- abs(dx)
        }
    f <- (2.*y+p)/pi
}



## SLA to structural parameter N and reverse
## @param SLA (cm^2/mg)
## @param N
## function A from Jacquemoud and Baret (1990)
## function B from Cecato et al  Remote Sensing of Environment 77 (2001) 22-33
## N = (1/(1/(SLA-0.01)))^(1/4)
# REMOTE SENS. ENVIRON. 34:75-91 
SLA2N <- function(SLA=NULL,N=NULL,method="A"){

if(method=="A"){
    if(!is.null(N)){
        return((0.1*N+0.025)/(N-0.9))
    } else if(!is.null(SLA)){

        return((1+36*SLA)/(40*SLA-4))
    } else {
        stop("inputs are NULL")

    }
}


if(method=="B"){
    if(!is.null(N)){
        return((100+N^4)/(100*N^4))
    } else if(!is.null(SLA)){
        return((1/(SLA-0.01))^(1/4))
    } else {
        stop("inputs are NULL")

    }
}

}

