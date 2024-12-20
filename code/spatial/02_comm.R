
# generate 15km site-by-taxon matrices for two different distribution data sets:
# - 810m SDM rasters
# - specimen occurrence points


# libraries
library(raster)
library(dplyr)
library(doParallel)
library(ape)
library(phytools)
library(data.table)
select <- dplyr::select


# cluster setup
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)


# load 15km CPPP raster template
template <- stack("data/cpad_cced_raster_15km.tif")[[2]]


# function to calculate 15km probability from 810m probabilities
aggProb <- function(x, ...){
      # probability of a single presence -- too liberal
      # 1 - prod(1 - x, ...)
      # average probability -- reflects proportion of range in target cell
      mean(x, ...)
}


# function to convert SDM raster to vector of occurrence probabilities per coarse grid cell
upscale <- function(x, template, vector = T){
      require(raster)
      r <- readRDS(x)
      d <- as.data.frame(rasterToPoints(r))
      coordinates(d) <- c("x", "y")
      crs(d) <- crs(r)
      u <- rasterize(d, template, field = "layer", fun = aggProb)
      u <- mask(u, template)
      if(vector) return(values(u))
      return(u)
}


# generate site-by-species matrix
f <- list.files("results/sdm", full.names = T)
spp <- sub("\\.rds", "", gsub(" ", "_", basename(f)))
sxsp <- foreach(x = f, .combine = "cbind") %dopar% upscale(x, template)
colnames(sxsp) <- spp
saveRDS(sxsp, "results/comm/site_by_species.rds")


# load species occurrence records
occ <- read.csv("data/caBryoOccs_2024_07_05-1.csv")

# function to convert point locations to vector of cell presences/absences
occ_to_rast <- function(x, template, vector = T){
      require(raster)
      pres <- occ[occ$name == x, ]
      coordinates(pres) <- c("Long", "Lat")
      projection(pres) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
      pres <- spTransform(pres, crs(template))
      pres$occ <- 1
      u <- raster::rasterize(pres, template, "occ")
      u[is.na(u)] <- 0
      u <- mask(u, template)
      if(vector) return(values(u))
      return(u)
}

sxocc <- foreach(x = spp, .combine = "cbind") %dopar% occ_to_rast(x, template)
colnames(sxocc) <- spp
saveRDS(sxocc, "results/comm/site_by_occurrence.rds")


# terminate compute cluster
stopCluster(cl)
