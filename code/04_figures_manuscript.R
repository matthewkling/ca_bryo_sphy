
# libraries
library(tidyverse)
library(spatialphy) # devtools::install_github("matthewkling/spatialphy")
library(raster)
library(ape)
library(ggtree)
library(patchwork)






# methods figure ======================================

# load occurrences and phylogeny
occ <- read.csv("data/caBryoOccs_2024_07_05-1.csv")
tree <- read.tree(file = "data/combined_chrono.tree")

# state boundary data for maps
select <- dplyr::select
cali <- map_data("state") %>%
      filter(region == "california")



p1 <- ggplot(occ, aes(Long, Lat)) +
      geom_polygon(data = cali, aes(long, lat, group = group), fill = "gray80", color = "gray80") +
      geom_point(size = .25) +
      theme_void()

p2 <- ggplot(tree) + 
      geom_tree(linewidth = .5) + 
      theme_tree() +
      coord_polar(theta = "y") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))

p <- p1 + p2 + plot_layout(nrow = 1, widths = c(1, 1.25))
ggsave("figures/manuscript/data.png", p, width = 15, height = 8, units = "in")


# data =======================================================

lnm <- readRDS("results/sphy/sphy_layer_names.rds")
r2df <- function(r) r %>% stack() %>% setNames(lnm) %>% rasterToPoints() %>% as.data.frame() %>% as_tibble()
dco <- r2df("results/sphy/combined_sphy_occ.tif")
dc <- r2df("results/sphy/combined_sphy.tif")
dm <- r2df("results/sphy/moss_sphy.tif")
dl <- r2df("results/sphy/liverwort_sphy.tif")

# projected state boundary
cali_pts <- select(cali, long, lat)
coordinates(cali_pts) <- c("long", "lat")
crs(cali_pts) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
cali_pts <- spTransform(cali_pts, crs(stack("results/sphy/combined_sphy.tif")))
cali <- coordinates(cali_pts) %>% bind_cols(cali) %>% rename(x = coords.x1, y = coords.x2)


# alpha diversity ===========================================

# vanilla diversity and endemism measures
pd <- dc %>%
      select(x, y,
             TR, TE, PD, PE) %>%
      gather(stat, value, -x, -y) %>%
      na.omit() %>%
      mutate(stat = factor(stat, levels = c("TR", "PD", "TE", "PE"),
                           labels = c("terminal\ndiversity", "phylogenetic\ndiversity",
                                      "terminal\nendemism", "phylogenetic\nendemism"))) %>%
      group_by(stat) %>%
      mutate(value = value / max(value))
p1 <- pd %>%
      ggplot(aes(x, y, fill = value)) +
      facet_wrap(~stat, nrow = 2) +
      geom_raster() +
      scale_x_continuous(expand = c(0,0)) +
      scale_fill_viridis_c(option = "B") +
      guides(fill = guide_colorbar(barwidth = 8)) +
      theme_void() +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold", size = 12)) +
      labs(fill = "value as fraction \nof maximum  ")


# randomization results for PD and RPD
pd <- dc %>%
      select(x, y, qPD, qRPD) %>%
      na.omit() %>%
      gather(stat, value, -x, -y) %>%
      na.omit() %>%
      mutate(stat = factor(stat, levels = c("qPD", "qRPD"),
                           labels = c("phylogenetic\ndiversity",
                                      "relative\nphylogenetic diversity")))
p2 <- pd %>%
      mutate(value = pmax(.001, pmin(.999, value))) %>%
      ggplot(aes(x, y, fill = value)) +
      facet_wrap(~stat, nrow = 2) +
      geom_raster() +
      scale_x_continuous(expand = c(0,0)) +
      # scale_fill_viridis_c(option = "B") +
      scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "gray", "forestgreen", "darkgreen"),
                           # values = c(0, .01, .3, .7, .99, 1),
                           breaks = c(.01, .1, .5, .9, .99),
                           labels = c(".01", ".1", ".5", ".9", ".99"),
                           trans = "logit") +
      guides(fill = guide_colorbar(barwidth = 6)) +
      theme_void() +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold", size = 12)) +
      labs(fill = "quantile in null \ndistribution")

p <- p1 + p2 + plot_layout(nrow = 1, widths = c(2, 1))
ggsave("figures/manuscript/alpha.png", 
       p, width = 8, height = 8, units = "in")





# beta diversity ============================================

pd <- dc %>%
      select(rgb1, rgb2, rgb3, x, y) %>%
      na.omit()

clr <- function(order, inversion){
      x <- pd
      x$color <- colors3d::colors3d(select(x, rgb1:rgb3), 
                                    order = order, inversion = inversion)
      x$order <- order
      x$inversion <- inversion
      x
}
pd <- expand_grid(order = 1:6, inversion = 1:8) %>%
      pmap_dfr(clr)

p <- ggplot() +
      facet_grid(order ~ inversion) +
      geom_raster(data = pd, aes(x, y, fill = color)) +
      scale_fill_identity() +
      theme_void() +
      coord_fixed()
ggsave("figures/manuscript/ordination_rgb_menu.png", 
       p, width = 9, height = 8, units = "in")

p1 <- ggplot() +
      geom_point(data = pd %>% filter(order == 3, inversion == 5),
                 aes(rgb1, rgb2, color = color), size = .5) +
      scale_color_identity() +
      theme_bw() +
      coord_fixed() +
      theme(plot.background = element_blank()) +
      labs(x = "NMDS 1\n(z-axis = NMDS 3)", y = "NMDS 2")

p2 <- ggplot() +
      geom_raster(data = pd %>% filter(order == 3, inversion == 5),
                  aes(x, y, fill = color)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_identity() +
      theme_void() +
      coord_fixed()

p <- p2 + inset_element(p1, .49, .49, 1, 1)
# ggsave("figures/manuscript/ordination_rgb.png", 
#        p, width = 6, height = 6, units = "in")


#### dendrogram ==================

# community phylogenetic distance matrix
tree_file <- "data/moss_chrono.tree"
comm_file <- "results/comm/site_by_species.rds"
comm <- readRDS(comm_file)
tree <- read.tree(file = tree_file)
tree$tip.label <- str_remove_all(tree$tip.label, "-")
colnames(comm) <- str_remove_all(colnames(comm), "-")
xcom <- comm[, colnames(comm) %in% tree$tip.label]
tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(comm)))
xcom <- xcom[, tree$tip.label]
xcom[is.na(xcom)] <- 0 # NA values not allowed in sphy functions
template <- stack("data/cpad_cced_raster_15km.tif")[[2]]
sp <- sphy(tree, xcom, template)
sp <- sphy_dist(sp, normalize = T, add = T)

# dendrogram
a <- rowSums(sp$occ) > 0
d <- as.matrix(sp$dist)
rownames(d) <- colnames(d) <- paste("cell", 1:ncol(d))
da <- d[a, a]
da[is.infinite(da)] <- max(da[!is.infinite(da)]) + 1000
da <- as.dist(da)

pd <- pd %>% 
      filter(order == 3, inversion == 5) %>%
      mutate(label = rownames(d)[a])

descendants <- function(tree, node, desc = NULL){
      if(is.null(desc)) desc <- vector()
      daughters <- tree %>% activate(edges) %>% filter(from == node) %>% pull(to)
      desc <- c(desc, daughters)
      w <- which(! tree %>% activate(nodes) %>% slice(daughters) %>% pull(leaf))
      if(length(w)>0) for(i in 1:length(w))
            desc <- descendants(tree, daughters[w[i]], desc)
      return(desc)
}

mean_color <- function(x) x %>% col2rgb() %>% apply(1, mean) %>% "/"(255) %>% as.list() %>% do.call(rgb, .)

interpolate_internal_colors <- function(graph){
      
      # assign colors to internal nodes
      graph <- activate(graph, nodes) %>% mutate(id = 1:length(leaf))
      ids <- graph %>% pull(id)
      colors <- graph %>% pull(color)
      for(i in 1:length(ids)){
            d <- descendants(graph, ids[i])
            if(length(d) > 0) colors[i] <- graph %>% slice(d) %>% pull(color) %>% mean_color()
      }
      graph <- graph %>% mutate(color = colors)
      
      # assign colors to edges
      graph <- graph %>% activate(edges) %>% mutate(color = pull(graph, color)[to],
                                                    height = pull(graph, height)[to])
      graph
}

graph <- da %>%
      hclust("average") %>%
      as.dendrogram() %>%
      reorder(pd$y) %>%
      as_tbl_graph() %>%
      activate(nodes) %>%
      left_join(select(pd, label, color)) %>%
      interpolate_internal_colors()

p3 <- ggraph(graph, 'dendrogram', height = height ^ .5) + 
      geom_edge_elbow(aes(edge_color = color, edge_width = height ^ 2)) +
      scale_edge_color_identity() +
      scale_edge_width_continuous(range = c(.1, 5)) +
      coord_flip() +
      theme_void() +
      theme(legend.position = "none")

pp <- p + p3 + plot_layout(nrow = 1)
ggsave("figures/manuscript/ordination_rgb.png", 
       pp, width = 10, height = 6, units = "in")



# nape ====================================================

# methods comparison for the combined dataset:
# - occ vs sdm
# - curveball vs quantize
# - canape vs snape

snape <- function(rand,
                  palette = c("gray90", "blue", "red"),
                  transform = function(x) x^2){
      
      if(inherits(rand, "RasterBrick")) rand <- values(rand)
      
      snape_colors <- function(peq, rpeq){
            colors <- colorRampPalette(palette[length(palette):2])(101)
            colors <- colors[round(rpeq * 100)+1]
            colors <- cbind(t(col2rgb(colors)), peq)
            colors <- apply(colors, 1, function(x) x[1:3] + (col2rgb(palette[1]) - x[1:3]) * (1 - transform(x[4])))
            colors <- rgb(colors[1,], colors[2,], colors[3,], maxColorValue = 255)
            colors
      }
      
      snp <- data.frame(qPE = rand[,"qPE"],
                        qCE = rand[,"qCE"],
                        qPCE = pmax(rand[,"qPE"], rand[,"qCE"]),
                        qRPE = rand[,"qRPE"])
      snp$color <- snape_colors(snp$qPCE, snp$qRPE)
      
      lgd <- expand.grid(qPCE = seq(.005, .995, .01),
                         qRPE = seq(.005, .995, .01))
      lgd$color <- snape_colors(lgd$qPCE, lgd$qRPE)
      
      return(list(snape = snp,
                  legend = lgd))
}

# non-sig, neo, paleo, mixed/super
pal <- c("gray80", "red", "dodgerblue", "purple")

# endem type
c1 <- bind_rows(dc %>% 
                      select(x, y, value = canape_endem_type) %>%
                      mutate(data = "SDMs", randomization = "curveball", pvalues = "threshold"),
                dco %>% 
                      select(x, y, value = canape_endem_type) %>%
                      mutate(data = "occurrences", randomization = "curveball", pvalues = "threshold"))

# sdm, quantize, canape
c2 <- dc %>%
      mutate(rpe_obs_p_lower = 1 - qRPE) %>%
      select(x, y,
             pe_obs_p_upper = qPE, 
             pe_alt_obs_p_upper = qCE, 
             rpe_obs_p_upper = qRPE, 
             rpe_obs_p_lower) %>% as.data.frame() %>%
      na.omit() %>%
      canaper::cpr_classify_endem() %>%
      mutate(value = as.integer(factor(endem_type, levels = c("not significant", "neo", "paleo", "mixed", "super"))) - 1) %>% 
      mutate(data = "SDMs", randomization = "quantize", pvalues = "threshold") %>%
      select(x, y, value, data, randomization, pvalues)

cp <- bind_rows(c1, c2) %>%
      na.omit() %>%
      mutate(value = ifelse(value == 4, 3, value),
             color = pal[value + 1],
             value = factor(value, levels = 0:3, labels = c("not-significant", "neo", "paleo", "mixed")))

# sdm, quantize, snape
r <- dc %>%
      select(x, y, PE, CE, qPE, qCE, qRPE) %>% as.data.frame() %>%
      na.omit()
sn <- snape(r, transform = function(x) x^10, palette = pal[c(1, 3, 4, 4, 2)])
r1 <- r %>% mutate(color = sn$snape$color,
                   data = "SDMs", randomization = "quantize", pvalues = "smooth")

# occ, quantize, snape
r <- dco %>%
      select(x, y, PE, CE, qPE, qCE, qRPE) %>% as.data.frame() %>%
      na.omit()
sn <- snape(r, transform = function(x) x^10, palette = pal[c(1, 3, 4, 4, 2)])
r2 <- r %>% mutate(color = sn$snape$color,
                   data = "occurrences", randomization = "quantize", pvalues = "smooth")

# sdm, curveball, snape
r <- dc %>%
      select(x, y, 
             PE = canape_pe_obs, 
             CE = canape_pe_alt_obs, 
             qPE = canape_pe_obs_p_upper, 
             qCE = canape_pe_alt_obs_p_upper, 
             qRPE = canape_rpe_obs_p_upper) %>% 
      as.data.frame() %>% na.omit()
sn <- snape(r, transform = function(x) x^10, palette = pal[c(1, 3, 4, 4, 2)])
r3 <- r %>% mutate(color = sn$snape$color,
                   data = "SDMs", randomization = "curveball", pvalues = "smooth")

# occ, curveball, snape
r <- dco %>%
      select(x, y, 
             PE = canape_pe_obs, 
             CE = canape_pe_alt_obs, 
             qPE = canape_pe_obs_p_upper, 
             qCE = canape_pe_alt_obs_p_upper, 
             qRPE = canape_rpe_obs_p_upper) %>% 
      as.data.frame() %>% na.omit()
sn <- snape(r, transform = function(x) x^10, palette = pal[c(1, 3, 4, 4, 2)])
r4 <- r %>% mutate(color = sn$snape$color,
                   data = "occurrences", randomization = "curveball", pvalues = "smooth")


r <- bind_rows(r1, r3, r4)
cpr <- bind_rows(cp, r) %>%
      mutate(randomization = factor(randomization, levels = c("curveball", "quantize")),
             pvalues = factor(pvalues, levels = c("threshold", "smooth")))

labels <- cpr %>% select(randomization, data, pvalues) %>% distinct() %>%
      arrange(pvalues, desc(data), randomization) %>%
      mutate(label = paste0("(", letters[1:6], ")"))

map <- ggplot() +
      facet_grid(pvalues ~ data + randomization, labeller = "label_both") +
      geom_tile(data = cpr, aes(x, y, fill = color)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      geom_text(data = labels, aes(label = label), 
                x = max(cpr$x), y = max(cpr$y), 
                size = 5, fontface = "bold", hjust = 1, vjust = 1) +
      scale_fill_identity() +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_rect(color = "black", fill = "black"),
            strip.text = element_text(color = "white"))

legend <- ggplot(cpr, aes(pmax(qPE, qCE), qRPE, color = color)) +
      geom_point(shape = 15) +
      annotate(geom = "segment",
               x = c(.95, .95), xend = c(.999, .999),
               y = c(.975, .025), yend = c(.975, .025),
               linewidth = .25) +
      geom_vline(xintercept = .95, linewidth = .25) +
      annotate(geom = "label", fill = "white",  size = 2.5, fontface = "bold",
               x = c(rep(.995, 3), .1), 
               y = c(.005, .5, .995, .5),
               label = c("neo", "mixed", "paleo", "non-significant"),,
               color = c(pal[c(2, 4, 3)], "gray40"),
               lineheight = .8, size = 2.75) +
      scale_color_identity() +
      scale_x_continuous(trans = "logit", breaks = c(.001, .05, .5, .95, .999)) +
      scale_y_continuous(trans = "logit", breaks = c(.001, .025, .5, .975, .999)) +
      theme_bw() +
      theme(plot.background = element_blank(),
            text = element_text(size = 8)) +
      labs(x = "max null quantile, PE | CE",
           y = "null quantile, RPE (PE/CE)") +
      coord_fixed()

p <- legend + map + plot_layout(nrow = 1, widths = c(1, 2))
ggsave("figures/manuscript/nape_method_comp.png", 
       p, width = 9, height = 5, units = "in")


# combined vs moss vs liverwort empirical results:
# sdm, quantize, snape
r <- bind_rows(dm %>% select(x, y, PE, CE, qPE, qCE, qRPE) %>% as.data.frame() %>% mutate(clade = "mosses"),
               dl %>% select(x, y, PE, CE, qPE, qCE, qRPE) %>% as.data.frame() %>% mutate(clade = "liverworts"),
               dc %>% select(x, y, PE, CE, qPE, qCE, qRPE) %>% as.data.frame() %>% mutate(clade = "mosses + liverworts")) %>%
      na.omit() %>%
      mutate(clade = factor(clade, levels = c("mosses", "liverworts", "mosses + liverworts")))
sn <- snape(r, transform = function(x) x^10, palette = pal[c(1, 3, 4, 4, 2)])

r$color <- sn$snape$color
r$qPCE <- sn$snape$qPCE

labels <- data.frame(clade = unique(r$clade),
                     label = paste0("(", letters[1:3], ")"))

p <- ggplot() +
      facet_wrap(~ clade, nrow = 2) +
      geom_tile(data = r, aes(x, y, fill = color)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      geom_text(data = labels, aes(label = label), 
                x = max(r$x), y = max(r$y), 
                size = 5, fontface = "bold", hjust = 1, vjust = 1) +
      scale_fill_identity() +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_rect(color = "black", fill = "black"),
            strip.text = element_text(color = "white"))

l <- ggplot(r, aes(qPCE, qRPE, color = color)) +
      geom_point(shape = 15) +
      annotate(geom = "segment",
               x = c(.95, .95), xend = c(.999, .999),
               y = c(.975, .025), yend = c(.975, .025),
               linewidth = .25) +
      geom_vline(xintercept = .95, linewidth = .25) +
      annotate(geom = "label", fill = "white",  size = 2.5, fontface = "bold",
               x = c(rep(.995, 3), .1), 
               y = c(.005, .5, .995, .5),
               label = c("neo", "mixed", "paleo", "non-significant"),,
               color = c(pal[c(2, 4, 3)], "gray40"),
               lineheight = .8, size = 2.75) +
      scale_color_identity() +
      scale_x_continuous(trans = "logit", breaks = c(.001, .05, .5, .95, .999)) +
      scale_y_continuous(trans = "logit", breaks = c(.001, .025, .5, .975, .999)) +
      theme_bw() +
      theme(plot.background = element_blank(),
            text = element_text(size = 8)) +
      labs(x = "max null quantile, PE | CE",
           y = "null quantile, RPE (PE/CE)") +
      coord_fixed()

pp <- p + inset_element(l, .5, 0, 1, .45)
ggsave("figures/manuscript/snape_moss_liverwort.png", 
       pp, width = 5, height = 6, units = "in")



# conservation =============================================

pd <- dc %>%
      select(x, y, value = priority) %>%
      na.omit()

p1 <- ggplot() +
      geom_raster(data = pd, aes(x, y, fill = value)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_viridis_c(option = "D", trans = "log10") +
      coord_fixed() +
      guides(fill = guide_colorbar(barheight = 8, reverse = TRUE)) +
      theme_void() +
      theme(legend.position = c(.6, .79),
            strip.text = element_text(face = "bold", size = 12)) +
      labs(fill = "priority")

p2 <- ggplot() +
      geom_raster(data = pd, aes(x, y, fill = value <= 50)) +
      geom_path(data = cali, aes(x, y, group = group), alpha = .5) +
      scale_fill_manual(values = c("gray80", "darkred")) +
      coord_fixed() +
      theme_void() +
      theme(legend.position = c(.6, .89),
            strip.text = element_text(face = "bold", size = 12)) +
      labs(fill = "top-50 priority")

p <- p1 + p2 + plot_layout(nrow = 1)
ggsave("figures/manuscript/conservation.png", 
       p, width = 10, height = 6, units = "in")
