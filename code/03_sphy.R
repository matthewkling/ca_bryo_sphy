
# libraries
library(tidyverse)
library(spatialphy) # devtools::install_github("matthewkling/spatialphy")
library(raster)
library(ape)

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


# function to run spatial phylogenetic analyses for a given phylogeny
sphylo <- function(tree_file = "data/mosses_rooted.tree", 
                   comm_file = "results/comm/site_by_species.rds"){
      
      comm <- readRDS(comm_file)
      tree <- read.tree(file = tree_file)[[2]]
      
      # intersect and cleanup
      tree$tip.label <- str_remove_all(tree$tip.label, "-")
      colnames(comm) <- str_remove_all(colnames(comm), "-")
      xcom <- comm[, colnames(comm) %in% tree$tip.label]
      tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(comm)))
      xcom <- xcom[, tree$tip.label]
      xcom[is.na(xcom)] <- 0 # NA values not allowed in sphy functions
      
      # construct spatial phylo object, run analyses
      sp <- sphy(tree, xcom, template)
      div <- sp %>% sphy_diversity()
      reg <- sp %>% sphy_regions(k = 5)
      con <- sp %>% sphy_prioritize(protected)
      
      # construct binary version of data set and run canaper randomizations
      # (note that a reasonable alternative stage to threshold would be after aggregating clade ranges;
      # but doing it before matches the output of canaper functions which has some appeal)
      threshold <- .25
      xcom_binary <- apply(xcom, 2, function(x) as.integer(x > (threshold * max(x, na.rm = T))))
      sp_binary <- sphy(tree, xcom_binary, template)
      rnd <- sp_binary %>% sphy_rand(null_model = "curveball")
      
      # combine
      stack(div, rnd, reg, con) %>%
            mask(protectedNA)
}

# run analyses for mosses and liverworts
moss <- sphylo("data/mosses_rooted.tree")
liverwort <- sphylo("data/liverworts_rooted.tree")

# save results
moss %>% writeRaster("results/sphy/moss_sphy.tif", overwrite = T)
liverwort %>% writeRaster("results/sphy/liverwort_sphy.tif", overwrite = T)
names(moss) %>% saveRDS("results/sphy/sphy_layer_names.rds")


### figures ###

# load results
moss <- stack("results/sphy/moss_sphy.tif") %>% setNames(readRDS("results/sphy/sphy_layer_names.rds"))
liverwort <- stack("results/sphy/liverwort_sphy.tif") %>% setNames(readRDS("results/sphy/sphy_layer_names.rds"))

# state bounday data for maps
select <- dplyr::select
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
for(x in names(moss)[grepl("obs_z|obs_p_upper", names(moss)) & !grepl("alt", names(moss))]){
      pd <- r2df(x) %>%
            mutate(value = ifelse(value == 0, NA, value))
      midpoint <- ifelse(grepl("obs_z", x), 0, .5)
      p <- ggplot() +
            facet_wrap(~taxon) +
            geom_raster(data = pd, aes(x, y, fill = value)) +
            geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
            scale_fill_gradient2(mid = "gray90", high = "red", low = "blue", na.value = "white", midpoint = midpoint) +
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



# SNAPE

library(terra)
library(patchwork)

d <- bind_rows(liverwort %>% rast() %>% as.data.frame(xy = T) %>% mutate(tree = "liverwort"),
                moss %>% rast() %>% as.data.frame(xy = T) %>% mutate(tree = "moss")) %>%
      group_by(tree) %>% 
      mutate(PE = PE / mean(PE), CE = CE / mean(CE)) %>%
      mutate(endem_type = factor(endem_type, levels = 0:4, 
                                 labels = c("not-significant", "neo", "paleo", "mixed", "super"))) 

canape <- ggplot() +
      facet_grid(method ~ tree, switch = "y") +
      geom_raster(data = d %>% mutate(method = "CANAPE"), 
                  aes(x, y, fill = endem_type)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_manual(values = c("gray95", "red", "blue", "violet", "darkorchid4"), na.value = "white", drop = F) +
      theme_void()

snape <- ggplot() +
      facet_grid(method ~ tree, switch = "y") +
      geom_raster(data = d %>% mutate(method = "SNAPE"),
                  aes(x, y, fill = PE / CE, alpha = PE + CE)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_gradientn(colors = c("darkred", "red", "violet", "violet", "blue", "darkblue"),
                           values = c(0, .25, .45, .55, .75, 1),
                           limits = exp(max(abs(log(d$PE / d$CE)), na.rm = T) * c(-1, 1)),
                           na.value = "white",
                           trans = "log10") +
      scale_alpha(range = 0:1) +
      theme_void() +
      theme(strip.text.x = element_blank())

p <- canape / snape
ggsave(paste0("figures/", "snape", ".png"), p, width = 8, height = 8, units = "in")



