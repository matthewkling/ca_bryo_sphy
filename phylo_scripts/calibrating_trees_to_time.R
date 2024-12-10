###############################
# Scaling phylogenies to time #
###############################

library(ape)
library(dplyr)
library(adephylo)

# Let's take each phylogeny and scale it to time using callibrations from other studies

#----------------------------------------------
# First liverworts

# read topology
liverworts <- read.nexus("data/rooted_trees/liverworts.nex")

# check the tree
plot.phylo(liverworts, cex = 0.7)
nodelabels(frame = "none", cex = 0.7, col = "blue")


# nodes to calibrate
# Using dates from Bechteler, J., Peñaloza‐Bojacá, G., Bell, D., Gordon Burleigh, J., McDaniel, S. F., Christine Davis, E., ... & Villarreal A, J. C. (2023). Comprehensive phylogenomic time tree of bryophytes reveals deep relationships and uncovers gene incongruences in the last 500 million years of diversification. American Journal of Botany, 110(11), e16249.
# Sup mat S8

col_names <- c("node", "age.min", "age.max")

# Root
root.liv    <- c(120, 444, 450)
names(root.liv) <- col_names

# Complex thalloids
com_th        <- c(121,  387, 395)
names(com_th) <- col_names

# Sphaerocarpales
sph        <- c(123, 272, 287)
names(sph) <- col_names

# Aytoniaceae
ayto        <- c(134, 169, 186)
names(ayto) <- col_names

# leafy clade

# porellales (including Ptilidium)
porellales        <- c(167, 330, 349)
names(porellales) <- col_names

# Marsupella + Blephanostora
marsup_clade        <- c(219, 213, 225)
names(marsup_clade) <- col_names

# make the calibration table
  # bind rows
calib <- bind_rows(root.liv, com_th, sph,  ayto, porellales, marsup_clade)

  # add last column
calib$soft.bounds <- FALSE

# calibration
liverworts_chrono <- ape::chronos(liverworts,
                                  calibration = calib,
                                  model = "relaxed")
# see the tree
plot.phylo(liverworts_chrono, cex = 0.7)

# save it to output
#write.tree(liverworts_chrono, "output/chronograms/liverworts_chrono.tree")

#----------------------------------------------
# Next moss

# read moss tree
moss <- read.nexus("data/rooted_trees/moss.nex")

# check the tree
plot.phylo(moss, cex = 0.3)
nodelabels(frame = "none", cex = 0.3, col = "blue")

# nodes to calibrate

# colnames
col_names <- c("node", "age.min", "age.max")

# Root
moss.root        <- c(467, 416, 423)
names(moss.root) <- col_names


# Psychomitrium + Discelium clade
psy.dis        <- c(508, 245, 258)
names(psy.dis) <- col_names


# Orthodontium and sister clade
orth        <- c(810, 135, 220)
names(orth) <- col_names


# Plagiomnium + Pohlia
plagio.poh        <- c(712, 58, 116)
names(plagio.poh) <- col_names

# Tortula + Dicranowesia clade
tor.dicr        <- c(617, 155, 213)
names(tor.dicr) <- col_names


# Grimmia + Blindia clade
grim.blin        <- c(540, 163, 193)
names(grim.blin) <- col_names


# make the calibration table
# bind rows
calib.moss <- bind_rows(moss.root, orth, plagio.poh, psy.dis, tor.dicr, grim.blin)

# add last column
calib.moss$soft.bounds <- FALSE

# calibration
moss_chrono <- ape::chronos(moss,
                            calibration = calib.moss,
                            model = "relaxed")
# see the tree
plot.phylo(moss_chrono, cex = 0.5)

# save it to output
#write.tree(moss_chrono, "output/chronograms/moss_chrono.tree")


#----------------------------------------------
# MERGE TWO TREES
# Divergence between mosses and liverworts in that paper is 486

moss.liv.div <- 486

# read the two calibrated trees
chrono_moss <- read.tree("output/chronograms/moss_chrono.tree")
chrono_liv  <- read.tree("output/chronograms/liverworts_chrono.tree")

# the length of the liverworts tree is 450
distRoot(chrono_liv, "Marsupella_bolanderi")
liv_TL <- 450

# the length of the moss tree is 423
distRoot(chrono_moss, "Kindbergia_oregana")
moss_TL <- 423

# add the root lengths to the trees
chrono_moss$root.edge <- moss.liv.div - moss_TL
chrono_liv$root.edge  <- moss.liv.div - liv_TL


# bind the trees
liverwort_moss_chrono <- bind.tree(chrono_liv, chrono_moss, where = "root", position = 36)

plot(test, cex = 0.2)

# checking that both liverwort and mosses have the same age
distRoot(test, "Marsupella_bolanderi")
distRoot(test, "Kindbergia_oregana")


# save it to output
#write.tree(liverwort_moss_chrono, "output/chronograms/combined_chrono.tree")





plot(liverwort_moss_chrono, cex = 0.5)
nodelabels(col = "blue", cex = 0.5)









