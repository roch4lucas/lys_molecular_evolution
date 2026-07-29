
#   Description: R script for PGLS comparison of gene and pseudogene number 
#   between Calliphoridae species  with distinct larval feeding habits
#   (facultative parasites x obligate parasites)
#
#   Biological question: do obligate parasite blowflies have lost lysozyme genes?
#
#   Author: Lucas Rocha


library(ape)
library(caper)
library(here)


# 1. File loading
tree <- read.tree(here("03_phylo_stats/data/speciestree_calliphoridae.nwk"))
table_genes <- read.table(here("03_phylo_stats/data/calliphoridae_genes.tsv"), 
                          header = TRUE, sep = "\t")
table_pseudo <- read.table(here("03_phylo_stats/data/calliphoridae_pseudos.tsv"), 
                           header = TRUE, sep = "\t")

# 2. Tree transformation
tree <- compute.brlen(tree, method = "Grafen") # --> Grafen's method (1989)
tree <- multi2di(tree) # --> binary tree without polytomies

# 3. Data preparation
table_genes$Larval_feeding_habit <- as.factor(table_genes$Larval_feeding_habit)
table_pseudo$Larval_feeding_habit <- as.factor(table_pseudo$Larval_feeding_habit)

# 4. VCV matrix inversion by caper
comp_data_genes <- comparative.data(phy = tree, 
                                     data = table_genes, 
                                     names.col = "Species", 
                                     vcv = TRUE)

comp_data_pseudo <- comparative.data(phy = tree, 
                               data = table_pseudo, 
                               names.col = "Species", 
                               vcv = TRUE)

# 5. PGLS - genes
model_pgls_genes <- pgls(LysBr_genes ~ Larval_feeding_habit, 
                          data = comp_data_genes, 
                          lambda = "ML")

summary(model_pgls_genes)

# 5. PGLS - pseudogenes
model_pgls_pseudo <- pgls(Pseudogenes ~ Larval_feeding_habit, 
                          data = comp_data_pseudo, 
                          lambda = "ML")

summary(model_pgls_pseudo)