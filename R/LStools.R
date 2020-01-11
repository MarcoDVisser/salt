require(tools)
################################################################################
## Landsat data analysis tools
################################################################################
## Nijmegen, March 2018
## 
## 
##                 ,,ggddY"""Ybbgg,,
##            ,agd888b,_ "Y8, ___`""Ybga,
##         ,gdP""88888888baa,.""8b    "888g,
##       ,dP"     ]888888888P'  "Y     `888Yb,
##     ,dP"      ,88888888P"  db,       "8P""Yb,
##    ,8"       ,888888888b, d8888a           "8,
##  ,8'        d88888888888,88P"' a,          `8,
## ,8'         88888888888888PP"  ""           `8,
## d'          I88888888888P"                   `b
## 8           `8"88P""Y8P'                      8
## 8            Y 8[  _ "                        8
## 8              "Y8d8b  "Y a                   8
## 8                 `""8d,   __                 8
## Y,                    `"8bd888b,             ,P
## `8,                     ,d8888888baaa       ,8'
##  `8,                    888888888888'      ,8'
##   `8a                   "8888888888I      a8'
##    `Yba                  `Y8888888P'    adP'
##      "Yba                 `888888P'   adY"
##        `"Yba,             d8888P" ,adP"'  
##           `"Y8baa,      ,d888P,ad8P"'     
##                ``""YYba8888P""''          
################################################################################


## Read key components from the metadata of a set of
## LS images contain in a directory
## @param files list of LS tar.gz files
## @param stats stats to read (see details)
## @param \\dots additional parameters pased to untar
compileInfo <- function(path,
                     stats=c("SPACECRAFT_ID",
                             "SENSOR",
                             "DATE_ACQUIRED",
                             "DATA_TYPE",
                             "MAP_PROJECTION",
                             "DATUM"
                             )
                    ,...){


    path <- file_path_as_absolute(path) # standardize file path
    
    fls <- list.files(path)


    if(any(grepl(".xml",fls))){

        xml <- fls[grepl(".xml",fls)] ## extract potential ARD xml
    }

    ## cloud cover available?
    if(any(grepl("pixel_qa.tif",fls))){
        clouds <- fls[grepl("pixel_qa.tif",fls)]
    } else{


        clouds <- NA
    }

     fls <- fls[grepl("_MTL.txt",fls)] ## subset to metadata
   

    if (!length(fls)){
       stop("path appears empty or contains no metadata elements")
    }


    meta <- lapply(paste0(path,"/",fls) , readLines)

    ## get key information
    REX <- paste(stats,collapse="|")
    metaShort <- meta[[1]][sapply(meta,function(X) grepl(REX,X))]
    tmp <- sapply(metaShort,function(X)
        gsub("\""," ",X,fixed=TRUE))
    tmp <- sapply(tmp,function(X) gsub("[[:space:]]", "", X))
    tmp <- lapply(tmp, function(X) strsplit(X,"="))
    tmp <- t(sapply(tmp,function(X) rbind(c(X[[1]][1],X[[1]][2]))))
    rownames(tmp) <- NULL

    if(length(fls)==1){
    tmp <- rbind(tmp,c("METADATA",fls),c("Pixel_QA",clouds))
    }

    if(exists("xml")){
    if(length(xml)==1){
    tmp <- rbind(tmp,c("XML_METADATA",xml))
    }
    }

    tmp <- data.frame(subject=tmp[,1],value=tmp[,2])
    class(tmp) <- c("quickInfo",class(tmp))
    return(tmp)
}


## Read information from USGS LS file 
## @param filen USGS landsat file name or vector of filenames
getFileInfo <- function(filen){

    filen <- sapply(filen,basename,USE.NAMES=FALSE)

    comp <- sapply(filen,strsplit, split="_",fixed=TRUE)
    
    compn <- c("SATTELITE-SENSOR","Processing correction level",
                     "Acq. date","Proc. date", "Collection number",
                     "Collection category","Product","band.extention")

    ## translate info
    
    return(transinf)

}

## Quickly get satelite and sensor name 
## from a compile info object
## @param info info object returned by compileInfo
getSatName <- function(info){

    a <- tolower(gsub("_","",info[info$subject=="SPACECRAFT_ID",2]))
    b <- tolower(gsub("_","",info[info$subject=="SENSOR_ID",2]))

    sat <- paste0(a,b)
    if(nchar(sat)>11){

        sat <- substr(sat,1,11)
    }
    return(sat)

}


## Return interger vector with band order
## for TC transform
## @param sat satname returned by getSatName
bandOrder <- function(sat){
    
binfo <- structure(list(sat = c("Landsat4TM", "Landsat5TM", "Landsat7ETM", 
"Landsat8OLI", "MODIS"), bands = c("1,2,3,4,5,7", "1,2,3,4,5,7", 
"1,2,3,4,5,7", "2,3,4,5,6,7", "1,2,3,4,5,6,7"), coefficients = c("Crist", 
"Crist", "Huang", "Baig", "Lobser"), data = c(1985L, 1985L, 2002L, 
2014L, 2007L), unit = c("reflectance", "reflectance", "reflectance", 
"reflectance", "reflectance")), .Names = c("sat", "bands", "coefficients", 
"data", "unit"), class = "data.frame", row.names = c(NA, -5L))

bands <- binfo$bands[sat==tolower(binfo$sat)]
as.integer(strsplit(bands,",")[[1]])
}


## plot a landsat images and automatically select bands
## for a set of common band configurations
## @param img a landsat image (as raster object)
## @param meta "ImageMetaData" = landsat metadata object
## or a "quickInfo" object
## @param type type of plot can be one of
## ("RGB","HVI","NDVI","NDWI","NBRI")
plot.lsat <- function(img,meta,type="RGB",saveFile=NULL){


    ## get rgb bands
    if(any("ImageMetaData"%in%class(meta))){

        inc <- as.logical(sapply(gsub("_","",tolower(LSbands$exp)), function(X) grepl(X,tolower(meta$SATELLITE))))
        inc2 <- as.logical(sapply(gsub("_","", tolower(LSbands$sensor)),
                                  function(X) grepl(X,gsub("_","",tolower(meta$SENSOR)))))

        focalbands <- LSbands[inc&inc2,]

    } else if(any("quickInfo"%in%class(meta))){
        focalbands <- LSbands[as.logical(sapply(LSbands$exp,
                 function(X) grepl(X,subset(meta,subject=="SPACECRAFT_ID")$value))),]

                  } else {

                      stop("meta of unrecognized format")
                  }
    

    if(type=='RGB'){

        bandorder <- match(c("red","green","blue"),focalbands$band)

        plotRGB(img,
            r = bandorder[1], g = bandorder[2], b = bandorder[3],
            stretch = "lin",
            axes = TRUE,
            main = "RGB bands\n true color estima")

    } else if(type=='HVI'){
        bandorder <- match(c("red","green","nir"),focalbands$band)

    plotRGB(img,
            r = bandorder[1], g = bandorder[2], b = bandorder[3],
            stretch = "lin",
            axes = TRUE,
            main = "Healthy vegatation image\n Landsat Bands red, green, NIR")

    } else if(type=='NDVI'){

               bandorder <- match(c("red","nir"),focalbands$band)

    inx <- spectralIndices(lsat, red = bandorder[1], nir = bandorder[2], indices = "NDVI")

        plot(inx,
            axes = TRUE,
            main = "NDVI")

        if(!is.null(saveFile)){

            writeRaster(inx,saveFile,options=c('TFW=YES'))
        }


    } else if(type=='NDWI'){

        bandorder <- match(c("green","nir"),focalbands$band)

        inx <- spectralIndices(lsat, green = bandorder[1], nir = bandorder[2], indices = "NDWI")

        plot(inx,
            axes = TRUE,
            main = "NDWI")

        if(!is.null(saveFile)){

            writeRaster(inx,saveFile,options=c('TFW=YES'))
        }


    } else if(type=='NBRI'){

                bandorder <- match(c("swir2","nir"),focalbands$band)

        inx <- spectralIndices(lsat, swir3 = bandorder[1], nir = bandorder[2], indices = "NBRI")

        plot(inx,
            axes = TRUE,
            main = "NBRI")

        if(!is.null(saveFile)){

            writeRaster(inx,saveFile,options=c('TFW=YES'))
        }



    } else{
    stop("unrecognized type")
    }


}    




