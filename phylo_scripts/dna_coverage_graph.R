#####################################
# Visualizing gene coverage results #
#####################################

#libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

#Reading in the data from MatrixMaker results
hep_data <- read.csv("output/per_group/liverwort_genbank_accessions.csv", header = T)
colnames(hep_data)[1] <- "taxon"
glimpse(hep_data)

#presence/abscence variables
hep_data <- hep_data %>% mutate(rbcLpa = ifelse(is.na(rbcL)==T, "0", "1"))
hep_data <- hep_data %>% mutate(matKpa = ifelse(is.na(matK)==T, "0", "1"))
hep_data <- hep_data %>% mutate(trnKpa = ifelse(is.na(trnK)==T, "0", "1"))
hep_data <- hep_data %>% mutate(trnLpa = ifelse(is.na(trnL)==T, "0", "1"))
hep_data <- hep_data %>% mutate(ndhFpa = ifelse(is.na(ndhF)==T, "0", "1"))
hep_data <- hep_data %>% mutate(ITSpa  = ifelse(is.na(ITS)==T, "0", "1"))


#extracting presence/abscence columns
columns <- c(1, 8:13)

 #percentage coverage graph
#data rearrangig
d <- hep_data[ , columns]
r<- nrow(d)
d$rbcLpa <- as.numeric(d$rbcLpa)
d$matKpa <- as.numeric(d$matKpa)
d$trnKpa <- as.numeric(d$trnKpa)
d$trnLpa <- as.numeric(d$trnLpa)
d$ndhFpa <- as.numeric(d$ndhFpa)
d$ITSpa  <- as.numeric(d$ITSpa)



#getting the percentages
genes <- colnames(hep_data[, 2:7])

num.seq <- d %>% select(2:7) %>%
  summarise_all(funs(sum)) %>% t() %>%
  as.data.frame()

num.seq <- num.seq %>% mutate(percentage_coverage = V1/r)
num.seq <- num.seq %>% mutate(gene = genes)


#plotting
ggplot(num.seq, aes(genes, percentage_coverage)) +
  geom_col() + geom_text(aes(label=round(percentage_coverage, digits=3)*100),
                         vjust=-0.30) +
  theme_bw()


# coverage graphs
# format for graph
data.pa <- hep_data[ , columns]
data.pa <- pivot_longer(data.pa, cols = c(2:7), names_to = "gene", values_to = "presence")
data.pa$presence <- as.factor(data.pa$presence)

data.pa %>%
  ggplot(aes(gene, taxon), fill=presence) +
  geom_tile(aes(fill=presence)) +
  theme(legend.position = "none", axis.text.x = element_text(size = rel(0.6), angle=90),
        axis.text.y = element_text(size = rel(0.2)))



###How many genes each terminal has?
markers <- data.pa %>% group_by(taxon) %>%
  filter(presence == "1") %>%
  summarise(num_mark = n())

genes_per_taxa<-markers %>% group_by(num_mark) %>%
  summarise(num_taxa = n())

#calculatig how many of the species don't have a marker
n <- data_frame(num_mark = 0, num_taxa = nrow(d) - nrow(markers))
#attaching
genes_per_taxa <- bind_rows(genes_per_taxa, n)

#graph
genes_per_taxa$num_mark = as.factor(genes_per_taxa$num_mark)
ggplot(genes_per_taxa, aes(num_mark, num_taxa)) +
  geom_col(fill="darkslategray4", colour= "gray38") +
  theme_grey() +
  labs(title = "Number of markers per taxon", x = "number of markers",
       y = 'Number of taxa') +
  geom_text(aes(label=num_taxa, vjust=-0.30))

# how many sequences will I have with only rbcl and matk
d_5 <- d[,-6]
nrow(d_5) #total


d_with5m <- d_5 %>% filter(rbcLpa != 0 | matKpa != 0 | trnKpa != 0 | trnLpa != 0 | ITSpa != 0 )





# FOR MOSSES -----------

#Reading in the data from MatrixMaker results
moss_data <- read.csv("output/per_group/moss_genbank_accessions.csv", header = T)
colnames(moss_data)[1] <- "taxon"
glimpse(moss_data)

#presence/abscence variables
moss_data <- moss_data %>% mutate(rbcLpa = ifelse(is.na(rbcL)==T, "0", "1"))
moss_data <- moss_data %>% mutate(matKpa = ifelse(is.na(matK)==T, "0", "1"))
moss_data <- moss_data %>% mutate(trnKpa = ifelse(is.na(trnK)==T, "0", "1"))
moss_data <- moss_data %>% mutate(trnLpa = ifelse(is.na(trnL)==T, "0", "1"))
moss_data <- moss_data %>% mutate(ndhFpa = ifelse(is.na(ndhF)==T, "0", "1"))
moss_data <- moss_data %>% mutate(ITSpa  = ifelse(is.na(ITS)==T, "0", "1"))


#extracting presence/abscence columns
columns <- c(1, 8:13)

#percentage coverage graph
#data rearrangig
d <- moss_data[ , columns]
r<- nrow(d)
d$rbcLpa <- as.numeric(d$rbcLpa)
d$matKpa <- as.numeric(d$matKpa)
d$trnKpa <- as.numeric(d$trnKpa)
d$trnLpa <- as.numeric(d$trnLpa)
d$ndhFpa <- as.numeric(d$ndhFpa)
d$ITSpa  <- as.numeric(d$ITSpa)



#getting the percentages
genes <- colnames(moss_data[, 2:7])

num.seq <- d %>% select(2:7) %>%
  summarise_all(funs(sum)) %>% t() %>%
  as.data.frame()

num.seq <- num.seq %>% mutate(percentage_coverage = V1/r)
num.seq <- num.seq %>% mutate(gene = genes)


#plotting
ggplot(num.seq, aes(genes, percentage_coverage)) +
  geom_col() + geom_text(aes(label=round(percentage_coverage, digits=3)*100),
                         vjust=-0.30) +
  theme_bw()


# coverage graphs
# format for graph
data.pa <- moss_data[ , columns]
data.pa <- pivot_longer(data.pa, cols = c(2:7), names_to = "gene", values_to = "presence")
data.pa$presence <- as.factor(data.pa$presence)

data.pa %>%
  ggplot(aes(gene, taxon), fill=presence) +
  geom_tile(aes(fill=presence)) +
  theme(legend.position = "none", axis.text.x = element_text(size = rel(0.6), angle=90),
        axis.text.y = element_text(size = rel(0.2)))



###How many genes each terminal has?
markers <- data.pa %>% group_by(taxon) %>%
  filter(presence == "1") %>%
  summarise(num_mark = n())

genes_per_taxa<-markers %>% group_by(num_mark) %>%
  summarise(num_taxa = n())

#calculatig how many of the species don't have a marker
n <- data_frame(num_mark = 0, num_taxa = nrow(d) - nrow(markers))
#attaching
genes_per_taxa <- bind_rows(genes_per_taxa, n)

#graph
genes_per_taxa$num_mark = as.factor(genes_per_taxa$num_mark)
ggplot(genes_per_taxa, aes(num_mark, num_taxa)) +
  geom_col(fill="darkslategray4", colour= "gray38") +
  theme_grey() +
  labs(title = "Number of markers per taxon", x = "number of markers",
       y = 'Number of taxa') +
  geom_text(aes(label=num_taxa, vjust=-0.30))



d_with5m <- d_5 %>% filter(rbcLpa != 0 | matKpa != 0 | trnKpa != 0 | trnLpa != 0 | ITSpa != 0 )



