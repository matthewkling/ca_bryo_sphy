
# libraries
library(tidyverse)
library(spatialphy) # devtools::install_github("matthewkling/spatialphy")
library(raster)
library(ape)


# protected areas =====================================

# load and transform existing protected area data
reserves <- raster("data/protection_status.tif")
protected <- reserves %>%
      rasterToPoints() %>%
      as.data.frame()
coordinates(protected) <- c("x", "y")
crs(protected) <- crs(reserves)
template <- stack("data/cpad_cced_raster_15km.tif")[[2]]
protected <- spTransform(protected, crs(template))
protected <- rasterize(protected, template, field = "protection_status", fun = mean)
protectedNA <- protected
protected[is.na(protected)] <- 1


# analyses using SDM inputs =====================================

# function to run spatial phylogenetic analyses for a given phylogeny
sphylo <- function(tree_file = "data/moss_chrono.tree", 
                   comm_file = "results/comm/site_by_species.rds"){
      
      comm <- readRDS(comm_file)
      tree <- read.tree(file = tree_file)#[[2]]
      
      # intersect and clean up
      tree$tip.label <- str_remove_all(tree$tip.label, "-")
      colnames(comm) <- str_remove_all(colnames(comm), "-")
      xcom <- comm[, colnames(comm) %in% tree$tip.label]
      tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(comm)))
      xcom <- xcom[, tree$tip.label]
      xcom[is.na(xcom)] <- 0 # NA values not allowed in sphy functions
      
      # detect whether input data are binary or quantitative
      binary <- all(xcom %in% 0:1)
      
      # construct spatial phylo object
      sp <- sphy(tree, xcom, template)
      
      # alpha diversity measures
      div <- sp %>% sphy_diversity()
      
      # regionalizations
      sp <- sphy_dist(sp, normalize = T, add = T)
      methods <- c("kmeans", "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
      reg <- purrr::map(methods, function(m) sphy_regions(sp, method = m, normalize = TRUE)) %>%
            stack() %>% setNames(paste0("region_", methods))
      rgb <- sphy_rgb(sp, "nmds", rank) %>% setNames(paste0("rgb", 1:3))
      
      # conservation prioritization
      con <- sp %>% sphy_prioritize(protected)
      
      # neo and paleo edemism
      n_rand <- 1000
      n_iter <- 100000
      rnd <- sp %>% sphy_rand(n_rand = n_rand, n_iter = n_iter, n_strata = ifelse(binary, 2, 5), 
                              transform = sqrt, n_cores = parallel::detectCores() - 2)
      
      # construct binary version of data set and run canaper randomizations
      # (note that a reasonable alternative stage to threshold would be after aggregating clade ranges;
      # but doing it before matches the output of canaper functions which has some appeal)
      if(binary){
            sp_binary <- sp
      }else{
            threshold <- .25
            xcom_binary <- apply(xcom, 2, function(x) as.integer(x > (threshold * max(x, na.rm = T))))
            sp_binary <- sphy(tree, xcom_binary, template)
      }
      cpr <- sp_binary %>% sphy_canape(n_reps = n_rand, n_iterations = n_iter)
      names(cpr) <- paste0("canape_", names(cpr))
      
      # combine
      stack(div, rnd, reg, rgb, con, cpr) %>%
            mask(protectedNA)
}

# generate results for mosses and liverworts
moss <- sphylo("data/moss_chrono.tree", "results/comm/site_by_species.rds")
liverwort <- sphylo("data/liverworts_chrono.tree", "results/comm/site_by_species.rds")
combined <- sphylo("data/combined_chrono.tree", "results/comm/site_by_species.rds")

# save to disk
moss %>% writeRaster("results/sphy/moss_sphy.tif", overwrite = T)
liverwort %>% writeRaster("results/sphy/liverwort_sphy.tif", overwrite = T)
combined %>% writeRaster("results/sphy/combined_sphy.tif", overwrite = T)
names(moss) %>% saveRDS("results/sphy/sphy_layer_names.rds")



# analyses using only occurrence records =====================================

# generate results for mosses and liverworts
moss <- sphylo("data/moss_chrono.tree", "results/comm/site_by_occurrence.rds")
liverwort <- sphylo("data/liverworts_chrono.tree", "results/comm/site_by_occurrence.rds")
combined <- sphylo("data/combined_chrono.tree", "results/comm/site_by_occurrence.rds")

# save to disk
moss %>% writeRaster("results/sphy/moss_sphy_occ.tif", overwrite = T)
liverwort %>% writeRaster("results/sphy/liverwort_sphy_occ.tif", overwrite = T)
combined %>% writeRaster("results/sphy/combined_sphy_occ.tif", overwrite = T)



