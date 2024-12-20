
# libraries
library(tidyverse)
library(phylospatial) # devtools::install_github("matthewkling/phylospatial")
library(terra)
library(raster)
library(ape)


# protected areas =====================================

# load and transform existing protected area data
reserves <- rast("data/protection_status.tif")
protected <- reserves %>% as.data.frame(xy = TRUE)
coordinates(protected) <- c("x", "y")
crs(protected) <- crs(reserves)
template <- rast("data/cpad_cced_raster_15km.tif")[[2]]
protected <- spTransform(protected, crs(template))
protected <- rasterize(vect(protected), template, field = "protection_status", fun = mean)
protectedNA <- protected
protected[is.na(protected)] <- 1


# analyses using SDM inputs =====================================

# function to run spatial phylogenetic analyses for a given phylogeny
sphylo <- function(tree_file = "results/chronograms/moss_chrono.tree", 
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
      
      species <- colnames(xcom)
      if(tree_file == "data/combined_chrono.tree") saveRDS(species, "results/sphy/species.rds")
      
      # detect whether input data are binary or quantitative
      binary <- all(xcom %in% 0:1)
      
      # construct spatial phylo object
      ps <- phylospatial(xcom, tree, template, data_type = "prob")
      
      # alpha diversity measures
      div <- ps_diversity(ps)
      
      # beta diversity
      ps <- ps_add_dissim(ps, "sorensen", normalize = TRUE)
      rgb <- ps_rgb(ps, "nmds", rank) %>% setNames(paste0("rgb", 1:3))
      methods <- c("kmeans", "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
      reg <- purrr::map(methods, function(m) ps_regions(ps, method = m, normalize = TRUE)) %>%
            rast() %>% setNames(paste0("region_", methods))
      
      # conservation prioritization
      con <- ps_prioritize(ps, protected)
      
      # neo and paleo edemism
      n_rand <- 1000
      n_iter <- 100000
      rnd <- ps_rand(ps, n_rand = n_rand, n_iter = n_iter, n_strata = ifelse(binary, 2, 5), 
                     transform = sqrt, n_cores = parallel::detectCores() - 2)
      
      # construct binary version of data set and run canaper randomizations
      # (note that a reasonable alternative stage to threshold would be after aggregating clade ranges;
      # but doing it before matches the output of canaper functions which has some appeal)
      if(binary){
            ps_binary <- ps
      }else{
            threshold <- .25
            xcom_binary <- apply(xcom, 2, function(x) as.integer(x > (threshold * max(x, na.rm = T))))
            ps_binary <- phylospatial(xcom_binary, tree, template, data_type = "binary")
      }
      cpr <- ps_binary %>% ps_canape(n_reps = n_rand, n_iterations = n_iter)
      names(cpr) <- paste0("canape_", names(cpr))
      
      # combine
      c(div, rnd, reg, rgb, con, cpr) %>%
            mask(protectedNA)
}

# run conservation prioritization only, using probabilistic method
con <- function(tree_file = "results/chronograms/moss_chrono.tree", 
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
      ps <- phylospatial(tree, xcom, template)
      
      # conservation prioritization
      ps_prioritize(ps, protected, method = "probable", max_iter = 50, n_reps = 1000, n_cores = 8)
}


# primary analyses using continuous SDM probabilities ===============================

# generate results for mosses and liverworts
moss <- sphylo("results/chronograms/moss_chrono.tree", "results/comm/site_by_species.rds")
liverwort <- sphylo("results/chronograms/liverworts_chrono.tree", "results/comm/site_by_species.rds")
combined <- sphylo("results/chronograms/combined_chrono.tree", "results/comm/site_by_species.rds")

# save to disk
moss %>% writeRaster("results/sphy/moss_sphy.tif", overwrite = T)
liverwort %>% writeRaster("results/sphy/liverwort_sphy.tif", overwrite = T)
combined %>% writeRaster("results/sphy/combined_sphy.tif", overwrite = T)
names(moss) %>% saveRDS("results/sphy/ps_layer_names.rds")

combined_con <- con("results/chronograms/combined_chrono.tree", "results/comm/site_by_species.rds")
combined_con %>% writeRaster("results/sphy/combined_prioritization.tif", overwrite = T)
names(combined_con) %>% saveRDS("results/sphy/cp_layer_names.rds")


# binary analyses using only occurrence records =====================================

# generate results for mosses and liverworts
moss <- sphylo("results/chronograms/moss_chrono.tree", "results/comm/site_by_occurrence.rds")
liverwort <- sphylo("results/chronograms/liverworts_chrono.tree", "results/comm/site_by_occurrence.rds")
combined <- sphylo("results/chronograms/combined_chrono.tree", "results/comm/site_by_occurrence.rds")

# save to disk
moss %>% writeRaster("results/sphy/moss_ps_occ.tif", overwrite = T)
liverwort %>% writeRaster("results/sphy/liverwort_ps_occ.tif", overwrite = T)
combined %>% writeRaster("results/sphy/combined_ps_occ.tif", overwrite = T)



# summary stats for manuscript ============================

species <- readRDS("results/sphy/species.rds")
occ <- read.csv("data/caBryoOccs_2024_07_05-1.csv") %>%
      mutate(name = str_remove_all(name, "-")) %>%
      filter(name %in% species)


