
### convert 810m distribution rasters into 15km site-by-taxon matrix


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
upscale <- function(x, template){
      require(raster)
      r <- readRDS(x)
      d <- as.data.frame(rasterToPoints(r))
      coordinates(d) <- c("x", "y")
      crs(d) <- crs(r)
      u <- rasterize(d, template, field = "layer", fun = aggProb)
      u[template[] == 0] <- NA # remove pixels outside state
      values(u)
}


# generate site-by-species matrix
f <- list.files("results/sdm", full.names = T)
sxsp <- foreach(x = f, .combine = "cbind") %dopar% upscale(x, template)
colnames(sxsp) <- sub("\\.rds", "", gsub(" ", "_", basename(f)))
saveRDS(sxsp, "results/comm/site_by_species.rds")


stopCluster(cl)
