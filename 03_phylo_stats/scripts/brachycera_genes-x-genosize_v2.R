##### PGLS: Genes x genome size pgls for Brachycera and other insects V2 #####

# Biological question: if brachyceran flies really have more lysozyme genes
# than all other insects, is this driven only (or mainly) by genome size expansion alone?

library(ape)
library(caper)
library(phytools)
library(dplyr)

# Load data
data <- read.table("insecta-reglinear.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tree <- read.tree("Insecta_v2.NEWICK.nwk")

# Prepare tree
tree$node.label <- NULL 
tree <- multi2di(tree)

# Divide data per taxa
data_brachycera <- subset(data, taxa == "Brachycera")
data_others <- subset(data, taxa != "Brachycera") # Or taxa == "Other"

# Create separate comparative objects

# Brachycera
comp_brachy <- comparative.data(phy = tree, 
                                data = data_brachycera, 
                                names.col = "species", 
                                vcv = TRUE, 
                                na.omit = TRUE)

# Outros
comp_others <- comparative.data(phy = tree, 
                                data = data_others, 
                                names.col = "species", 
                                vcv = TRUE, 
                                na.omit = TRUE)

# PGLS

# Brachycera
print("--- RESULTS: BRACHYCERA ---")
model_brachy <- pgls(lys_genes ~ genome_size, data = comp_brachy, lambda = "ML")
summary(model_brachy)

# Outros
print("--- RESULTS: OTHERS ---")
model_others <- pgls(lys_genes ~ genome_size, data = comp_others, lambda = "ML")
summary(model_outhers)

# Combined plot (for analysis check)

plot(lys_genes ~ genome_size, data = data,
     main = "Genes Lys x Genome (Groups)",
     xlab = "Genome size (Mb)",
     ylab = "Lys gene number",
     type = "n")

# Brachycera data
points(data_brachycera$genome_size, data_brachycera$lys_genes, pch = 19, col = "red")

# Others data
points(data_others$genome_size, data_others$lys_genes, pch = 19, col = "blue")

# PGLS lines
abline(model_brachy, col = "red", lwd = 2)  # Red line
abline(model_others, col = "blue", lwd = 2) # Blue line

legend("topright", legend = c("Brachycera", "Others"), 
       col = c("red", "blue"), pch = 19, lwd = 2)
