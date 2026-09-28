
#   Description: R script for PGLS comparison of gene number between Brachycera 
#   (true flies) and all other insects
#
#   Biological question: do brachyceran flies have significantly more lysozyme
#   genes than all other insects?
#
#   Author: Lucas Rocha


library(ape)
library(caper)
library(phytools)
library(dplyr)
library(here)


# 1. Load data
tree <- read.tree(here("03_phylo_stats/data/Insecta_v2.NEWICK2.nwk"))
data <- read.table(here("03_phylo_stats/data/insecta_genes-taxa2.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Bug fixed
# 2. Maintain data present in both archives
same <- intersect(tree$tip.label, data$species)
data <- data[data$species %in% same, ]
tree <- keep.tip(tree, same)

# 3. Resolve polytomies and zero branches
tree <- multi2di(tree)
tree$edge.length[tree$edge.length == 0] <- 1e-6

# 4. Data preparation
data$taxa <- as.factor(data$taxa)

# 5. Paired data object
comp_data <- comparative.data(phy = tree, 
                               data = data, 
                               names.col = "species", 
                               vcv = TRUE)

# 6. PGLS
# lambda = "ML" for automatic phylo signal estimation
model_pgls <- pgls(log(lys_genes) ~ taxa, 
                    data = comp_data, 
                    lambda = "ML")

summary(model_pgls)
