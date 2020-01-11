require(SearchTrees)
require(proxy) ## distance between two sets of points    

## 10 x 10 hash table
makeHash <- function(dat,dx=10,dy=dx){


    ## mesh mid points
    dat$yh <- trunc(dat$gy/dy)
    dat$xh <- trunc(dat$gx/dx)

    ## build hash table reference
    xh <- 0:trunc(ceiling(max(dat$gx))/dx)
    yh <- 0:trunc(ceiling(max(dat$gy))/dy)

    hash <- expand.grid(xh,yh)
    hash <- matrix(c(hash[,1],hash[,2],hash[,1]*dx,hash[,2]*dy)
                  ,ncol=4,byrow=FALSE)
    
    colnames(hash) <- c("xh","yh","x","y")
    

    attr(dat,which="hash") <- hash
    attr(dat,which="hashdxdy") <- c(dx,dy)

    return(dat)
}


hashSubset <- function(dat,x,y,r=10){

    hash <- attr(dat,"hash")
    mesh <- attr(dat,"hashdxdy")
    dx <- mesh[1];dy <- mesh[2]

    ## mid points
    hash[,3:4] <- hash[,3:4]+c(dx/2,dy/2)

    d <- sqrt((x-hash[,3])^2 + (y-hash[,4])^2)
    R <- r+sqrt((dx/2)^2+(dy/2)^2) ## reference diameter
    inc <- hash[d<R,]
    inc2 <- dat$xh%in%inc[,1] & dat$yh%in%inc[,2]

    return(dat[inc2,])
}


normSubset <- function(dat,x,y,r=10){

    
    d <- sqrt((dat$gx-x)^2 + (dat$gy-y)^2)
    inc <- d<r

   return(dat[inc,])   
}


qSubset <- function(dat,qtr,x,y,r=10){

    inc <- rectLookup(qtr,xlims=c(x-r,x+r),
                      ylims=c(y-r,y+r))
    return(dat[inc,])

}




## Make Biomass or height map
## @param dat dataset with x,y tree x and y coordinates
## and statistics h and b arespread out over radius r
## @param dx,dy mesh size
## @param celldx numerical integration size
mapMaker <- function(dat,dx=5,dy=dx,celldx=.1){

    ranger <- c(floor(min(dat$gx)),ceiling(max(dat$gx)),
                floor(min(dat$gy)),ceiling(max(dat$gy)))

    xi <- seq(ranger[1]+(dx/2),ranger[2],dx)
    yi <- seq(ranger[3]+(dy/2),ranger[4],dy)
    cxi <- seq(0+(celldx/2),dx,celldx) - dx/2 ## within cell
    cyi <- seq(0+(celldx/2),dy,celldx) - dy/2 ## within cell
        
    mapgrid <- expand.grid(xi,yi)
    cellgrid <- expand.grid(cxi,cyi)
    mapindex <- expand.grid(1:length(xi),1:length(yi))

    rasth <- matrix(numeric(nrow(mapindex)),ncol=length(xi))
    rastb <- matrix(numeric(nrow(mapindex)),ncol=length(xi))
    
    colnames(mapgrid) <- c("x","y")
    #mapgrid[,1] <- mapgrid[,1] + dx/2
    #mapgrid[,1] <- mapgrid[,1] + dy/2

    Rmax <- max(dat$r) ## max search radius
    qtr <- createTree(dat,treeType = "quad", dataType = "point",columns=1:2) # quandtree
    
    for(i in 1:nrow(mapgrid)){

        xtmp <- mapgrid[i,1]
        ytmp <- mapgrid[i,2]

        ## center focalgrid at xtmp and ytmp
        focalgrid <- cellgrid
        focalgrid[,1] <- focalgrid[,1]+xtmp
        focalgrid[,2] <- focalgrid[,2]+ytmp

        ## get local trees
        tmp <- qSubset(dat,qtr,xtmp,ytmp,sqrt(dx^2+dy^2)+Rmax)

        if(dim(tmp)[1]>0){
        #tmp <- normSubset(dat,xtmp,ytmp,sqrt(dx^2+dy^2)+Rmax)

        # distmat trees and focalgrid
        Rmat <- proxy::dist(focalgrid,tmp[,c("gx","gy")]) 

        ## Grid points within crown radii
        hitMat <- sapply(1:nrow(tmp),function(X) Rmat[,X]<=tmp$r[X])
        Hmat <- matrix(rep(tmp$h,nrow(cellgrid)),ncol=nrow(tmp),byrow=TRUE)
        Bmat <- matrix(rep(tmp$b,nrow(cellgrid)),ncol=nrow(tmp),byrow=TRUE)
        Rmat <- pi*matrix(rep(tmp$r,nrow(cellgrid)),ncol=nrow(tmp),byrow=TRUE)^2
        
        ## find height at each point
        cellvalh <- apply(hitMat*Hmat,1,max)
        cellvalb <- apply(((hitMat*celldx)/Rmat)*Bmat,1,mean)

                ## final height for cell
        
            rasth[mapindex[i,2],mapindex[i,1]] <- mean(cellvalh)
            rastb[mapindex[i,2],mapindex[i,1]] <- mean(cellvalb)
        } else {
            
            rasth[mapindex[i,2],mapindex[i,1]] <- NA
            rastb[mapindex[i,2],mapindex[i,1]] <- NA
        
        }
        
        cat("\r Creating heightmap", 100*round(i/nrow(mapgrid),4),"% \r")
        
                 }
    
    cat("\n")
       
    return(list(rasth,rastb))

}



## Magnetic north to transform
## p1 - p4; corner south west, south east
## north west, north east resp 
ctrans <- function(p1=MnciSub[1,],
                   p2=MnciSub[2,],
                   p3=MnciSub[4,],
                   p4=MnciSub[3,]
                   ){

    ## calculate transformation of y over x
    slYoX <- mean(p2[2]-p1[2],p4[2]-p3[2])
    slXoY <- mean(p3[1]-p1[1],p4[1]-p2[1])
    scalYX <- mean(p2[1]-p1[1],p4[1]-p3[1])
    scalXY <- mean(p3[2]-p1[2],p4[2]-p2[2])
    full <- rbind(p1,p2,p3,p4)
    xmm <- c(min(full[,1]),max(full[,1]))
    ymm <- c(min(full[,2]),max(full[,2]))

    
    return(list("Y"=slYoX/scalYX,"X"=slXoY/scalXY,
                "South-West"=p1,
                "Xmin-Max"=xmm,
                "Ymin-Max"=ymm                
                ))
    
}
