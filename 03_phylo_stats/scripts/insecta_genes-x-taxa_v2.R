##### PGLS: Genes x Insecta taxa V2 #####

# Biological question: do brachyceran flies have significantly more lysozyme
# genes than all other insects?

library(ape)
library(caper)
library(phytools)
library(dplyr)

# Load data
tree <- read.tree("Insecta_v2.NEWICK.nwk")
data <- read.table("insecta_genes-taxa2.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Bug fixed
# Maintain data present in both archives
same <- intersect(tree$tip.label, data$species)
data <- data[data$species %in% same, ]
tree <- keep.tip(tree, same)

# Resolve polytomies
tree <- multi2di(tree)

# Data preparation
data$taxa <- as.factor(data$taxa)

# Paired data object
dados_comp <- comparative.data(phy = tree, 
                               data = data, 
                               names.col = "species", 
                               vcv = TRUE)

# PGLS
# lambda = "ML" for automatic phylo signal estimation
model_pgls <- pgls(lys_genes ~ taxa, 
                    data = dados_comp, 
                    lambda = "ML")

summary(model_pgls)

# ==============================================================================
# Visualization 1: Phenogram (Traitgram) - FIXED COLOR MAPPING
# ==============================================================================
cat("\nGenerating Phenogram with mapped ancestral states...\n")

# 1. Extract the discrete trait (Suborder) as a named factor
suborder_trait <- setNames(as.factor(data_matched$diptera_taxa1), rownames(data_matched))

# 2. Map the discrete trait on the phylogeny using stochastic character mapping
# This infers internal node states so the branches can be colored properly.
cat("Running stochastic character mapping (make.simmap) for branch colors...\n")
mapped_tree <- make.simmap(tree_matched, suborder_trait, model = "ER", nsim = 1, message = FALSE)

# 3. Define the color palette mapping EXACTLY to the factor levels
colors_mapped <- setNames(c("#D55E00", "#0072B2"), c("Brachycera", "Nematocera"))

pdf("Phenogram_LysGenes.pdf", width = 10, height = 8)
phenogram(mapped_tree, lys_vector, 
          colors = colors_mapped, 
          spread.labels = TRUE, 
          spread.cost = c(1, 0), 
          fsize = 0.6,
          xlab = "Time (relative)", 
          ylab = "Number of Lysozyme Genes",
          main = "Phenogram of Lysozyme Gene Evolution in Diptera")
legend("topleft", legend = c("Brachycera", "Nematocera"), 
       col = c("#D55E00", "#0072B2"), lty = 1, lwd = 3, bty = "n")
dev.off()
