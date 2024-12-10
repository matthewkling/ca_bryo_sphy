
library(data.table)
library(dplyr)
library(tidyverse)

# read the data
names <- read.csv("data/CAthesaurus2024.csv", header = T)


# get rid of underscores in the first column
names <- names %>% separate(ACCtrinom, into = c("gen", "sp", "sub"), sep= "_")
names <- names %>% unite("ACCtrinom", c("gen", "sp", "sub"), sep = " ", na.rm = T)

# get rid of underscores in the second column
names <- names %>% separate(SYNtri, into = c("gen", "sp", "sub"), sep= "_")
names <- names %>% unite("SYNtri", c("gen", "sp", "sub"), sep = " ", na.rm = T)


# separate by liverworts and mosses
hep_names  <- names[names$grp == "hep", ]
moss_names <- names[names$grp == "moss", ]


# produce a list in which each item is a list of the synonyms

# for liverworts

  # get a subset of the df with only the names
hep_names_only <- hep_names[ , 1:2]

  # get the names in wider format
hep_wide_names <- hep_names_only %>%
  group_by(ACCtrinom) %>%
  mutate(rn = paste0("syn", row_number())) %>%
  spread(rn, SYNtri, fill = "")

col.order <- c("ACCtrinom", "syn1", "syn2", "syn3", "syn4", "syn5","syn6", "syn7",
               "syn8", "syn9")

hep_wide_names <- hep_wide_names[, col.order]




# for mosses

# get a subset of the df with only the names
moss_names_only <- moss_names[ , 1:2]

# get the names in wider format
moss_wide_names <- moss_names_only %>%
  group_by(ACCtrinom) %>%
  mutate(rn = paste0("syn", row_number())) %>%
  spread(rn, SYNtri, fill = "")

col.order <- c("ACCtrinom", "syn1", "syn2", "syn3", "syn4", "syn5","syn6", "syn7",
               "syn8", "syn9", "syn10", "syn11", "syn12", "syn13", "syn14", "syn15",
               "syn16", "syn17", "syn18")

moss_wide_names <- moss_wide_names[, col.order]
