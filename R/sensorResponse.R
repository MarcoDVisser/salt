#' resampleBands
#'
#' Function for spectral response modelling of common
#' sensor platforms
#' @param wl a vector of wavelengths
#' @param spectra a vector of reflectance/radiance values
#' @param sensor sensor bands to simulate / resample
#' @export
resampleBands <- function(wl,spectra,sensor,...){

    dat <- cbind(wl,spectra)
    colnames(dat) <- c("wl","R")
    class(dat) <- sensor
    resampleSensor(dat,...)
}



## S3 resample sensor function
resampleSensor <- function(x){
    UseMethod("resampleSensor",x)
}


## Hyperion has two spectrometers, their range overlaps.
## Hyperion has 220 unique channels, but sensors together 
## have 242[VNIR (70 channels, 356 nm - 1058 nm),
##  SWIR (172 channels, 852 nm - 2577 nm)]
## Here the overlapping channels are removed,
## and replaces with the less sensitive channels.
resampleSensor.Hyperion <- function(x){


    data(response_functions)
    res <- response_functions[["Hyperion"]]

    ## get rid of overlap
    res <- res[-(71:92),]

    final <- array(dim=c(nrow(res),3))
    final[,1] <- 1:nrow(res)
    final[,2] <- res[,"mid"]
    
    for(i in 1:nrow(res)){
        normdens <- dnorm(x[,"wl"],res[i,"mid"],res[i,"sd"])
        final[i,3] <- sum((x[,"R"]*normdens)/sum(normdens))
               
    }

  
    colnames(final) <- c("Band","Wavelength","Response")
  
    return(final)
}

#' Data on sensor response function for various Remote Sensing Platforms
#'
#' @section details:
#' ***********************************************************************
#' response_functions (Jan, 20th 2020)
#'  The dataset contains the following response functions (in a list):
#' ***********************************************************************
#'  \itemize{
#'
#' \item [1] = Hyperion, data gives mean (nw) and sigma of best fitting Gaussuin
#' for each band (note that band differ for each of the 256 spatial pixels!)
#' This data uses averages. Data was originally supplied by the EO-1 Instrument Teams.
#' http://www.eoc.csiro.au/hswww/oz_pi/specresp.htm
#' }
#'
#' @section references:
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name response_functions
#' @usage data(response_functions)
NULL
