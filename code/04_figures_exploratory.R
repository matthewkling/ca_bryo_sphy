

# libraries
library(tidyverse)
library(spatialphy) # devtools::install_github("matthewkling/spatialphy")
library(raster)
library(ape)



for(dataset in c("occ", "sdm")){
      
      figdir <- paste0("figures/exploratory_", dataset)
      
      files <- c("results/sphy/moss_sphy.tif",
                 "results/sphy/liverwort_sphy.tif",
                 "results/sphy/combined_sphy.tif")
      if(dataset == "occ") files <- str_replace(files, ".tif", "_occ.tif")
      
      # load results
      lnm <- readRDS("results/sphy/sphy_layer_names.rds")
      moss <- stack(files[1]) %>% setNames(lnm)
      liverwort <- stack(files[2]) %>% setNames(lnm)
      combined <- stack(files[3]) %>% setNames(lnm)
      
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
      ggsave(paste0(figdir, "/priority.png"), p, width = 9, height = 5, units = "in")
      
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
            ggsave(paste0(figdir, "/", x, ".png"), p, width = 9, height = 5, units = "in")
      }
      
      # rgb
      pd <- r2df("rgb1") %>% rename(r = value) %>% 
            left_join(r2df("rgb2") %>% rename(g = value)) %>% 
            left_join(r2df("rgb3") %>% rename(b = value))
      pd$color <- rgb(pd$r, pd$g, pd$b)
      
      p <- ggplot() +
            facet_wrap(~taxon) +
            geom_raster(data = pd, aes(x, y, fill = color)) +
            geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
            scale_fill_identity() +
            theme_void() +
            coord_fixed()
      ggsave(paste0(figdir, "/ordination_rgb.png"), p, width = 9, height = 5, units = "in")
      
      
      
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
            ggsave(paste0(figdir, "/", x, ".png"), p, width = 9, height = 5, units = "in")
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
            ggsave(paste0(figdir, "/", x, ".png"), p, width = 9, height = 5, units = "in")
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
            ggsave(paste0(figdir, "/", x, ".png"), p, width = 9, height = 5, units = "in")
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
      ggsave(paste0(figdir, "/canape.png"), p, width = 9, height = 5, units = "in")
      
      
      
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
            ggsave(paste0(figdir, "/snape_", name, ".png"), p, width = 8, height = 6, units = "in")
      }
      
      snape_plot(moss, "moss")
      snape_plot(liverwort, "liverwort")
      snape_plot(combined, "combined")
      
}
