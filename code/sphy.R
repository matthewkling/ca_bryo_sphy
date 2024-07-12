
library(tidyverse)
library(spatialphy) # devtools::install_github("matthewkling/spatialphy")
library(raster)
library(ape)



#### data setup ####

# template spatial grid (15 km res)
grid <- raster("data/dem15km.tif")

# occurrence data
occ_data <- read_csv("data/caBryoOccs_2024_07_05-1.csv")
coordinates(occ_data) <- c("Long", "Lat")
crs(occ_data) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# map occurrences to grid cells
occ_data <- spTransform(occ_data, crs(grid))
grid <- extend(grid, occ_data) %>% extend(1)
grid[] <- 1:ncell(grid)
occ_data$cell <- raster::extract(grid, occ_data)
comm <- occ_data %>%
      as.data.frame() %>%
      dplyr::select(name, cell) %>%
      distinct() %>%
      mutate(pres = 1) %>%
      spread(name, pres, fill = 0) %>%
      arrange(cell)
cell <- comm$cell
i <- grid[] %in% cell
comm <- dplyr::select(comm, -cell)

# convert to raster
comm_ml <- apply(comm, 2, function(x){
      y <- grid
      z <- rep(0, ncell(grid))
      z[i] <- x
      y[] <- z
      return(y)}) %>%
      stack()

# existing reserves
reserves <- raster("data/protection_status.tif")
protected <- reserves %>%
      rasterToPoints() %>%
      as.data.frame()
coordinates(protected) <- c("x", "y")
crs(protected) <- crs(reserves)
protected <- spTransform(protected, crs(comm_ml))
protected <- rasterize(protected, comm_ml[[1]], field = "protection_status", fun = mean)
protectedNA <- protected
protected[is.na(protected)] <- 1


#### run analyses ####

sphylo <- function(x){
      
      # load tree and align names
      tree <- read.tree(file = x)[[2]]
      comm <- subset(comm_ml, intersect(tree$tip.label, names(comm_ml)))
      tree <- drop.tip(tree, setdiff(tree$tip.label, names(comm)))
      sp <- sphy(tree, comm)
      
      # run spatial phylogenetic analyses
      div <- sp %>% sphy_diversity()
      rnd <- sp %>% sphy_rand()
      reg <- sp %>% sphy_regions(k = 8)
      con <- sp %>% sphy_prioritize(protected)
      
      stack(div, rnd, reg, con) %>%
            mask(protectedNA)
}

moss <- sphylo("data/mosses_rooted.tree")
liverwort <- sphylo("data/liverworts_rooted.tree")



#### export results and figures ####

# save prioritization rasters
writeRaster(moss$priority, "results/moss_priority.tif")
writeRaster(liverwort$priority, "results/liverwort_priority.tif")

# state bounday data for maps
cali <- map_data("state") %>%
      filter(region == "california")
cali_pts <- select(cali, long, lat)
coordinates(cali_pts) <- c("long", "lat")
crs(cali_pts) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
cali_pts <- spTransform(cali_pts, crs(moss))
cali <- coordinates(cali_pts) %>% bind_cols(cali) %>% rename(x = coords.x1, y = coords.x2)

# function to reformat data for plotting
r2df <- function(var){
      stack(moss[[var]],
            liverwort[[var]]) %>%
            setNames(c("moss", "liverwort")) %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            as_tibble() %>%
            gather(taxon, value, moss, liverwort) %>%
            filter(!is.na(value))
}

# priority
pd <- r2df("priority")
p <- ggplot() +
      facet_wrap(~taxon) +
      geom_raster(data = pd, aes(x, y, fill = value)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_viridis_c(trans = "log10", na.value = "white") +
      theme_void() +
      coord_fixed() +
      theme(legend.position = "bottom") +
      guides(fill = guide_colorbar(barwidth = 12)) +
      labs(fill = "priority")
ggsave("figures/priority.png", p, width = 8, height = 6, units = "in")

# phyloregion
pd <- r2df("phyloregion") %>%
      mutate(value = factor(value))
p <- ggplot() +
      facet_wrap(~taxon) +
      geom_raster(data = pd, aes(x, y, fill = value)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_brewer(type = "qual") +
      theme_void() +
      coord_fixed() +
      theme(legend.position = "bottom") +
      labs(fill = "phyloregion")
ggsave("figures/region.png", p, width = 8, height = 6, units = "in")

# diversity measures
for(x in names(moss)[1:10]){
      pd <- r2df(x) %>%
            mutate(value = ifelse(value == 0, NA, value))
      p <- ggplot() +
            facet_wrap(~taxon) +
            geom_raster(data = pd, aes(x, y, fill = value)) +
            geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
            scale_fill_viridis_c(na.value = "white", trans = "sqrt") +
            theme_void() +
            coord_fixed() +
            theme(legend.position = "bottom") +
            guides(fill = guide_colorbar(barwidth = 12)) +
            labs(fill = x)
      ggsave(paste0("figures/", x, ".png"), p, width = 8, height = 6, units = "in")
}

# rand
for(x in names(moss)[grepl("obs_z", names(moss)) & !grepl("alt", names(moss))]){
      pd <- r2df(x) %>%
            mutate(value = ifelse(value == 0, NA, value))
      p <- ggplot() +
            facet_wrap(~taxon) +
            geom_raster(data = pd, aes(x, y, fill = value)) +
            geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
            scale_fill_gradient2(mid = "gray90", high = "red", low = "blue", na.value = "white") +
            theme_void() +
            coord_fixed() +
            theme(legend.position = "bottom") +
            guides(fill = guide_colorbar(barwidth = 12)) +
            labs(fill = x)
      ggsave(paste0("figures/", x, ".png"), p, width = 8, height = 6, units = "in")
}

# endem type
pd <- r2df("endem_type") %>%
      mutate(value = factor(value, levels = 0:4, labels = c("not-significant", "neo", "paleo", "mixed", "super")))
p <- ggplot() +
      facet_wrap(~taxon) +
      geom_raster(data = pd, aes(x, y, fill = value)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_manual(values = c("gray80", "red", "blue", "violet", "darkorchid4"), na.value = "white") +
      theme_void() +
      coord_fixed() +
      theme(legend.position = "bottom") +
      labs(fill = "CANAPE classification ")
ggsave(paste0("figures/", "canape", ".png"), p, width = 8, height = 6, units = "in")
