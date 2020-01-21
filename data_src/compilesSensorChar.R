## compile sensor characteristics for some platforms
setwd("~/symlinks/git/salt/src/")

## save list of response functions
N <- 1 ## add more sensors
response_functions <- vector("list",N)

################################################################################
## Land sat
################################################################################
bands <- read.csv("./src/bands.csv")
save(bands,file="./data/bands.rda")

################################################################################
## Hyperion

fwhm <- read.table("./Hyp_band/Hyperion_cen_fwhm.dat") ## full width half maximum
bw <- read.table("./Hyp_band/BandWidthL0.dat") ## bandwidth for each pixel (FWHM)
cwl <- read.table("./Hyp_band/SpectralL0_revA.dat") ## central wavelength in nm
afwhm <- read.table("./Hyperion_cen_fwhm_av.dat",header=TRUE,skip=1) ## averaged fwhm
 


##  create spectral response gaussians for hyperion
## get half of the halfmax width
gm <- (1/sqrt(2*pi)) # guassian max
hm <- abs(qnorm(gm/2)) # right hand halfmax
sig <- (afwhm$fwhm/2)*hm ## response function half max

resFunc <- array(dim=c(nrow(afwhm),2))
colnames(resFunc) <- c("mid","sd")

resFunc[,1] <- afwhm$cwl
resFunc[,2] <- sig

response_functions[[1]]<-resFunc


################################################################################
## finalize and save
################################################################################
names(response_functions) <- c("Hyperion")
save(response_functions,file="../data/response_functions.rda")
