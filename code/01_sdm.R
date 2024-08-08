
### fit distribution models for every species at 810 m resolution


# libraries
options(java.parameters = "-Xmx1g" )
require(dismo)
require(raster)
library(data.table)
library(dplyr)
library(doParallel)


# cluster setup
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)


# load occurrence and background data
occ <- read.csv("data/caBryoOccs_2024_07_05-1.csv")
bg <- readRDS("data/10000_CA_plants_bg_810m_occbias.rdata")


# load climate data
clim_files <- list.files("data/bcm", pattern="filled3x", full.names=T) 
clim_vars <- substr(basename(clim_files), 1, 3)
clim <- lapply(clim_files, readRDS) %>% do.call("stack", .)
names(clim) <- clim_vars
clim <- clim[[sort(names(clim))]]


# load landscape intactness data
intactness <- raster("data/intactness_810m.tif")


# fit models
results <- foreach(pres = split(occ, occ$name), .combine="c") %dopar% {
      
      options(java.parameters = "-Xmx1g" )
      require(dismo)
      require(raster) 
      
      # format presence data
      pres <- as.data.frame(pres)
      coordinates(pres) <- c("Long", "Lat")
      projection(pres) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
      pres <- spTransform(pres, crs(clim))
      
      # drop occurrences that fall in the water
      valid <- !is.na(raster::extract(clim[[1]], pres))
      pres <- pres[valid,]
      
      # fit and project maxent model
      mx <- try(maxent(clim, pres, bg,
                       args=c("-a", "-z", "outputformat=raw", "maximumbackground=10000", 
                              "nothreshold", "nohinge")))
      if(class(mx)=="try-error") return("failure")
      suit <- predict(mx, clim)
      
      # fit and project distance model
      sigma <- 25 # km
      prox <- distanceFromPoints(clim[[1]], pres)
      prox <- mask(prox, clim[[1]])
      gauss <- function(x, sigma) exp(-.5 * (x/sigma) ^ 2)
      prox <- calc(prox, function(x) gauss(x, sigma*1000))
      
      # combine climate, distance, and intactness; save output
      pred <- suit * prox * intactness
      saveRDS(pred, paste0("results/sdm/", pres$name[1], ".rds"))
      return("success")
}

stopCluster(cl)
