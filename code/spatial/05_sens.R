
## sensitivity analysis to sigma parameter ##



# 01_sdm.R ===========

sdm(sigma = 10, outdir = "results/sdm10km/")
sdm(sigma = 50, outdir = "results/sdm50km/")
sdm(sigma = 100, outdir = "results/sdm100km/")



# 02_comm.R ==========

# generate site-by-species matrix
make_comm <- function(sigma){
      f <- list.files(paste0("results/sdm", sigma, "km"), full.names = T)
      spp <- sub("\\.rds", "", gsub(" ", "_", basename(f)))
      sxsp <- foreach(x = f, .combine = "cbind", 
                      .export = c("upscale", "template", "aggProb")) %dopar% upscale(x, template)
      colnames(sxsp) <- spp
      saveRDS(sxsp, paste0("results/comm", sigma, "km/site_by_species.rds"))
}

make_comm(10)
make_comm(50)
make_comm(100)



# 03_sphy.R ===============

rand <- function(tree_file = "results/chronograms/combined_chrono.tree", 
                 comm_file = "results/comm/site_by_species.rds", 
                 n_rand = 1000, n_iter = 100000){
      build_ps(tree_file, comm_file) %>%
            ps_rand(fun = "quantize", method = "curveball", 
                    n_rand = n_rand, burnin = n_iter, n_strata = 5, 
                    transform = sqrt, n_cores = parallel::detectCores() - 2)
}

rand(comm_file = "results/comm10km/site_by_species.rds") %>% writeRaster("results/sphy/rand_p10.tif", overwrite = T)
rand(comm_file = "results/comm/site_by_species.rds") %>% writeRaster("results/sphy/rand_p25.tif", overwrite = T)
rand(comm_file = "results/comm50km/site_by_species.rds") %>% writeRaster("results/sphy/rand_p50.tif", overwrite = T)
rand(comm_file = "results/comm100km/site_by_species.rds") %>% writeRaster("results/sphy/rand_p100.tif", overwrite = T)


# 04_figures_manuscript.R ===============

r2df <- function(r) r %>% stack() %>% rasterToPoints() %>% as.data.frame() %>% as_tibble()
r <- bind_rows(r2df("results/sphy/rand_p10.tif") %>% select(x, y, qPE, qCE, qRPE, qPD) %>% as.data.frame() %>% mutate(sigma = "10 km"),
               r2df("results/sphy/rand_p25.tif") %>% select(x, y, qPE, qCE, qRPE, qPD) %>% as.data.frame() %>% mutate(sigma = "25 km"),
               r2df("results/sphy/rand_p50.tif") %>% select(x, y, qPE, qCE, qRPE, qPD) %>% as.data.frame() %>% mutate(sigma = "50 km"),
               r2df("results/sphy/rand_p100.tif") %>% select(x, y, qPE, qCE, qRPE, qPD) %>% as.data.frame() %>% mutate(sigma = "100 km")) %>%
      na.omit() %>%
      mutate(sigma = factor(sigma, levels = c("10 km", "25 km", "50 km", "100 km")))
sn <- snape(r, transform = function(x) x^10, palette = cpal)

r$color <- sn$snape$color
r$qPCE <- sn$snape$qPCE


p <- ggplot() +
      facet_wrap(~ sigma, nrow = 1, labeller = label_both) +
      geom_tile(data = r, aes(x, y, fill = color)) +
      scale_fill_identity() +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_rect(color = "black", fill = "black"),
            strip.text = element_text(color = "white")) +
      coord_fixed()

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
               color = c(pal[c(2, 3, 4)], "gray40"),
               lineheight = .8, size = 2.75) +
      scale_color_identity() +
      scale_x_continuous(trans = "logit", breaks = c(.001, .05, .5, .9, .95, .999)) +
      scale_y_continuous(trans = "logit", breaks = c(.001, .025, .1, .9, .975, .999)) +
      theme_bw() +
      theme(plot.background = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
            text = element_text(size = 8)) +
      labs(x = "max null quantile, PE | CE",
           y = "null quantile, RPE (PE/CE)") +
      coord_fixed()

pp <- p + l + plot_layout(heights = c(1, 1))
ggsave("figures/manuscript/fig_S7.tiff",
       pp, width = 8, height = 6, units = "in")

p <- r %>%
      select(x, y, qPD, qPE, sigma) %>%
      gather(metric, value, -x, -y, -sigma) %>%
      mutate(value = pmax(.001, pmin(.999, value)),
             metric = str_remove(metric, "q")) %>%
      ggplot(aes(x, y, fill = value)) +
      facet_grid(metric~sigma, labeller = label_both) +
      geom_raster() +
      scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "gray", "forestgreen", "darkgreen"),
                           breaks = c(.01, .1, .5, .9, .99),
                           labels = c(".01", ".1", ".5", ".9", ".99"),
                           trans = "logit") +
      guides(fill = guide_colorbar(barwidth = 12)) +
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.text = element_text(color = "white"),
            strip.background = element_rect(color = "black", fill = "black"),
            legend.position = "bottom") +
      coord_fixed() +
      labs(fill = "Significance (quantile in null distribution)")
ggsave("figures/manuscript/fig_S6.tiff",
       p, width = 8, height = 5.5, units = "in")


