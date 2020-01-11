## Obtain temporal distribution of landsat data for a study area
## @param metaDir the directory containing the USGS metadata
## on archived landsat 1,2,3,4,5,7,and 8 scenes
## @param pathRow data.frame containing the paths and rows
## you are interested in (must contain tow columsn called
## "PATH" and "ROW".
queryArchive <- function(metaDir,pathRow){
    metaFiles<-list.files(metaDir)
    
    LSscens<-vector("list",length(metaFiles))
    names(LSscens) <-metaFiles
    
    cat("\r Searching archives: ", 0,"% \r")
    for(i in 1:length(metaFiles)){
        gzmeta <- gzfile(paste0(metaDir,metaFiles[i])) ## open connection to gz file
        metaDat<- read.csv(gzmeta) ## load meta data
        idx<-paste0(metaDat$path,metaDat$row)
        idy<-paste0(pathRow[,"PATH"],pathRow[,"ROW"])
        LSscens[[i]]<-subset(metaDat,idx%in%idy)
        cat("\r Searching archives:", round(i/length(metaFiles),3)*100,"% \r")
    }
    cat("\n")
    
    ## retrieve useful data
    
    colnms<-c("path","row","acquisitionDate",
              "LANDSAT_PRODUCT_ID","sensor",
              "cloudCover","CLOUD_COVER_LAND","dayOrNight")
    
    metaScene<-do.call("rbind",lapply(LSscens,function(X) X[colnms]))
    metaScene$acquisitionDate<-as.Date(metaScene$acquisitionDate)
    metaScene$scode<-as.numeric(factor(metaScene$sensor))
    labelz<-tapply(metaScene$scode,metaScene$sensor,mean)
    metaScene$cloudCover[metaScene$cloudCover<0]<-0
    metaScene$cloudCover[metaScene$cloudCover>100]<-100
    metaScene$dayOrNight<-toupper(metaScene$dayOrNight)

    return(metaScene)
}
