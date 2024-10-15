
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
sphylo <- function(tree_file = "data/mosses_chrono.tree", 
                   comm_file = "results/comm/site_by_species.rds"){
      
      comm <- readRDS(comm_file)
      tree <- read.tree(file = tree_file)#[[2]]
      
      # intersect and cleanup
      tree$tip.label <- str_remove_all(tree$tip.label, "-")
      colnames(comm) <- str_remove_all(colnames(comm), "-")
      xcom <- comm[, colnames(comm) %in% tree$tip.label]
      tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(comm)))
      xcom <- xcom[, tree$tip.label]
      xcom[is.na(xcom)] <- 0 # NA values not allowed in sphy functions
      
      # construct spatial phylo object
      sp <- sphy(tree, xcom, template)
      
      # alpha diversity measures
      div <- sp %>% sphy_diversity()
      rnd <- sp %>% sphy_rand(n_rand = 1000, n_strata = 5, transform = sqrt, n_cores = parallel::detectCores() - 2)
      
      # regionalizations
      sp <- sphy_dist(sp, normalize = T, add = T)
      methods <- c("kmeans", "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
      reg <- map(methods, function(m) sphy_regions(sp, method = m, normalize = TRUE)) %>%
            stack() %>% setNames(paste0("region_", methods))
      
      # conservation prioritization
      con <- sp %>% sphy_prioritize(protected)
      
      # construct binary version of data set and run canaper randomizations
      # (note that a reasonable alternative stage to threshold would be after aggregating clade ranges;
      # but doing it before matches the output of canaper functions which has some appeal)
      threshold <- .25
      xcom_binary <- apply(xcom, 2, function(x) as.integer(x > (threshold * max(x, na.rm = T))))
      sp_binary <- sphy(tree, xcom_binary, template)
      cpr <- sp_binary %>% sphy_canape(n_reps = 1000)
      names(cpr) <- paste0("canape_", names(cpr))
      
      # combine
      stack(div, rnd, reg, con, cpr) %>%
            mask(protectedNA)
}

# run analyses for mosses and liverworts
moss <- sphylo("data/moss_chrono.tree")
liverwort <- sphylo("data/liverworts_chrono.tree")
combined <- sphylo("data/combined_chrono.tree")

# save results
moss %>% writeRaster("results/sphy/moss_sphy.tif", overwrite = T)
liverwort %>% writeRaster("results/sphy/liverwort_sphy.tif", overwrite = T)
combined %>% writeRaster("results/sphy/combined_sphy.tif", overwrite = T)
names(moss) %>% saveRDS("results/sphy/sphy_layer_names.rds")


### figures ###

# load results
lnm <- readRDS("results/sphy/sphy_layer_names.rds")
moss <- stack("results/sphy/moss_sphy.tif") %>% setNames(lnm)
liverwort <- stack("results/sphy/liverwort_sphy.tif") %>% setNames(lnm)
combined <- stack("results/sphy/combined_sphy.tif") %>% setNames(lnm)

# state boundary data for maps
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
            liverwort[[var]],
            combined[[var]]) %>%
            setNames(c("moss", "liverwort", "combined")) %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            as_tibble() %>%
            gather(taxon, value, moss, liverwort, combined) %>%
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
ggsave("figures/priority.png", p, width = 9, height = 5, units = "in")

# regions
for(x in names(moss)[grepl("region", names(moss))]){
      pd <- r2df(x) %>%
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
      ggsave(paste0("figures/", x, ".png"), p, width = 9, height = 5, units = "in")
}

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
      ggsave(paste0("figures/", x, ".png"), p, width = 9, height = 5, units = "in")
}

# diversity randomization measures
for(x in names(moss)[11:20]){
      pd <- r2df(x)
      p <- ggplot() +
            facet_wrap(~taxon) +
            geom_raster(data = pd, aes(x, y, fill = value)) +
            geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
            scale_fill_viridis_c(na.value = "white", limits = 0:1) +
            theme_void() +
            coord_fixed() +
            theme(legend.position = "bottom") +
            guides(fill = guide_colorbar(barwidth = 12)) +
            labs(fill = x)
      ggsave(paste0("figures/", x, ".png"), p, width = 9, height = 5, units = "in")
}

# rand
for(x in names(moss)[grepl("obs_p_upper", names(moss)) & !grepl("alt", names(moss))]){
      pd <- r2df(x) %>%
            mutate(value = ifelse(value == 0, NA, value))
      # midpoint <- ifelse(grepl("obs_z", x), 0, .5)
      p <- ggplot() +
            facet_wrap(~taxon) +
            geom_raster(data = pd, aes(x, y, fill = value)) +
            geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
            # scale_fill_gradient2(mid = "gray90", high = "red", low = "blue", na.value = "white", midpoint = midpoint) +
            scale_fill_viridis_c(na.value = "white", limits = 0:1) +
            theme_void() +
            coord_fixed() +
            theme(legend.position = "bottom") +
            guides(fill = guide_colorbar(barwidth = 12)) +
            labs(fill = x)
      ggsave(paste0("figures/", x, ".png"), p, width = 9, height = 5, units = "in")
}

# endem type
pd <- r2df("canape_endem_type") %>%
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
ggsave(paste0("figures/", "canape", ".png"), p, width = 9, height = 5, units = "in")



# SNAPE

library(terra)
library(patchwork)

snape_plot <- function(x, name){
      r <- x %>% rast() %>% as.data.frame(xy = T) %>%
            select(x, y, PE, CE, qPE, qCE, qRPE) %>%
            na.omit()
      
      sn <- snape(r, transform = function(x) x^10, palette = c("gray80", "blue", "red"))
      
      r$color <- sn$snape$color
      
      map <- ggplot() +
            geom_tile(data = r, aes(x, y, fill = color)) +
            geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
            scale_fill_identity() +
            theme_void()
      
      legend <- ggplot(sn$legend, aes(qPCE, qRPE, fill = color)) +
            geom_tile() +
            scale_fill_identity() +
            theme_bw() +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))
      
      scatter <- ggplot(r, aes(PE, CE, color = color)) +
            geom_point() +
            scale_color_identity() +
            theme_bw()
      
      p <- map + legend + scatter +
            plot_layout(design = c("AB
                             AC"), widths = c(2, 1))
      ggsave(paste0("figures/snape_", name, ".png"), p, width = 8, height = 6, units = "in")
}

snape_plot(moss, "moss")
snape_plot(liverwort, "liverwort")
snape_plot(combined, "combined")



