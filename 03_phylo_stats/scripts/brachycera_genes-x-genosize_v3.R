
#   Description: R script for PGLS analyses and continuous trait model fitting 
#   of gene number and genome size between brachyceran flies and mosquitoes
#   - Includes scatterplot and phenogram ploting
#   - Verbose script
#
#   Biological question: do brachyceran flies have significantly more lysozyme
#   genes than their sister group "nematocera" (mosquitoes)? Is there a
#   correlation between genome size and the number of genes in these groups?
#
#   Author: Lucas Rocha


library(ape)
library(geiger)
library(phytools)
library(caper)
library(ggplot2)
library(dplyr)
library(here)


# 1. Load data and tree --------------------------------------------------------
cat("Loading data and tree...\n")
tree <- read.tree(here("03_phylo_stats/data/speciestree_insecta.nwk"))

tree <- compute.brlen(tree, method = "Grafen") # --> Grafen's method (1989)
tree <- multi2di(tree) # --> binary tree without polytomies

data <- read.table(here("03_phylo_stats/data/diptera_data.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(data) <- data$species

# 2. Match tree and tata & clean numeric nariables -----------------------------
cat("Matching tree tips with dataset rows...\n")
matched_obj <- treedata(tree, data, sort = TRUE)
tree_matched <- matched_obj$phy
data_matched <- as.data.frame(matched_obj$data, stringsAsFactors = FALSE)

# Clean number formats (remove thousand separator dots, replace decimal commas)
data_matched$lys_genes <- as.numeric(data_matched$lys_genes)
clean_genome <- gsub("\\.", "", data_matched$genome_length_Mb) 
clean_genome <- gsub(",", ".", clean_genome)                   
data_matched$genome_length_Mb <- as.numeric(clean_genome)

# Restore species column explicitly
data_matched$species <- rownames(data_matched)
lys_vector <- setNames(data_matched$lys_genes, rownames(data_matched))

# 4. Computes phylogenetic signal & model fitting across different taxa --------
cat("\n--- Running Phylogenetic Signal & Evolutionary Models for Subsets ---\n")

groups_to_test <- c("All", "Brachycera", "Nematocera")
results_list <- list()

for (grp in groups_to_test) {
  cat(sprintf("Processing group: %s...\n", grp))
  
  # Subset data based on the group
  if (grp == "All") {
    sub_data <- data_matched
  } else {
    sub_data <- data_matched[data_matched$diptera_taxa1 == grp, ]
  }
  
  # Match the tree to the specific subset
  matched_sub <- treedata(tree_matched, sub_data, sort = TRUE, warnings = FALSE)
  phy_sub <- matched_sub$phy
  dat_sub <- as.data.frame(matched_sub$data, stringsAsFactors = FALSE)
  
  # Named vector for the trait in this subset
  lys_sub <- setNames(as.numeric(dat_sub$lys_genes), rownames(dat_sub))
  
  # Phylogenetic Signal
  sig_K <- phylosig(phy_sub, lys_sub, method = "K", test = TRUE, nsim = 999)
  sig_L <- phylosig(phy_sub, lys_sub, method = "lambda", test = TRUE)
  
  # Evolutionary Models (Continuous Trait)
  # Suppress warnings for optimization bounds which are common in small subsets
  fit_BM <- suppressWarnings(fitContinuous(phy_sub, lys_sub, model = "BM"))
  fit_OU <- suppressWarnings(fitContinuous(phy_sub, lys_sub, model = "OU"))
  fit_EB <- suppressWarnings(fitContinuous(phy_sub, lys_sub, model = "EB"))
  
  # AICc extraction and best model selection
  aicc_scores <- c(BM = fit_BM$opt$aicc, OU = fit_OU$opt$aicc, EB = fit_EB$opt$aicc)
  best_mod <- names(which.min(aicc_scores))
  
  # Append to results list
  results_list[[grp]] <- data.frame(
    Group = grp,
    N_taxa = Ntip(phy_sub),
    Blombergs_K = round(sig_K$K, 3),
    P_val_K = round(sig_K$P, 4),
    Pagels_Lambda = round(sig_L$lambda, 3),
    P_val_Lambda = round(sig_L$P, 4),
    AICc_BM = round(aicc_scores["BM"], 2),
    AICc_OU = round(aicc_scores["OU"], 2),
    AICc_EB = round(aicc_scores["EB"], 2),
    Best_Model = best_mod,
    stringsAsFactors = FALSE
  )
}

# Combine and export to TSV
final_subset_results <- bind_rows(results_list)
write.table(final_subset_results, here("03_phylo_stats/PhyloSignal_Models_Results.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("-> Exported subset analyses to 'PhyloSignal_Models_Results.tsv'\n")

final_subset_results

# 5. Phylogenetic Generalized Least Squares (pGLS) -----------------------------
cat("\n--- Running PGLS Models ---\n")
comp_data <- comparative.data(phy = tree_matched, 
                              data = data_matched, 
                              names.col = "species", 
                              vcv = TRUE, 
                              na.omit = TRUE)

pgls_mod1 <- pgls(lys_genes ~ genome_length_Mb, data = comp_data, lambda = "ML")
pgls_mod2 <- pgls(lys_genes ~ genome_length_Mb * diptera_taxa1, data = comp_data, lambda = "ML")

cat("\n[PGLS Base Model Summary]\n")
print(summary(pgls_mod1))
cat("\n[PGLS Interaction Model Summary]\n")
print(summary(pgls_mod2))


# ==============================================================================
# Phenogram (Traitgram)
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

#pdf("Phenogram_LysGenes.pdf", width = 10, height = 8)
phenogram(mapped_tree, lys_vector, 
          colors = colors_mapped, 
          spread.labels = TRUE, 
          spread.cost = c(1, 0), 
          fsize = 0.6,
          xlab = "Time (Relative)", 
          ylab = "Number of Lysozyme Genes",
          main = "Phenogram of Lysozyme Gene Evolution in Diptera")
legend("topleft", legend = c("Brachycera", "Nematocera"), 
       col = c("#D55E00", "#0072B2"), lty = 1, lwd = 3, bty = "n")
#dev.off()

# ==============================================================================
# PGLS scatterplot
# ==============================================================================
cat("Generating PGLS Scatter Plot...\n")

coefs <- coef(pgls_mod2)
int_brach <- coefs["(Intercept)"]
slope_brach <- coefs["genome_length_Mb"]

nem_int_term <- "diptera_taxa1Nematocera"
nem_slope_term <- "genome_length_Mb:diptera_taxa1Nematocera"

if (nem_int_term %in% names(coefs) && nem_slope_term %in% names(coefs)) {
  int_nem <- int_brach + coefs[nem_int_term]
  slope_nem <- slope_brach + coefs[nem_slope_term]
} else {
  int_nem <- int_brach
  slope_nem <- slope_brach
}

pgls_plot <- ggplot(data_matched, aes(x = genome_length_Mb, y = lys_genes, color = diptera_taxa1)) +
  geom_point(size = 3, alpha = 0.8) +
  # BUG FIX: Replaced 'size' with 'linewidth' to comply with ggplot2 3.4.0+
  geom_abline(intercept = int_brach, slope = slope_brach, 
              color = "#D55E00", linewidth = 1.2, linetype = "dashed") +
  geom_abline(intercept = int_nem, slope = slope_nem, 
              color = "#0072B2", linewidth = 1.2, linetype = "dashed") +
  scale_color_manual(values = c("Brachycera" = "#D55E00", "Nematocera" = "#0072B2")) +
  theme_minimal(base_size = 14) +
  labs(x = "Genome Length (Mb)",
       y = "Lysozyme Genes",
       color = "Suborder",
       title = "Phylogenetic Generalized Least Squares (PGLS)",
       subtitle = "Relationship between genome size and lysozyme gene count") +
  theme_minimal() +
  theme(
    # Fundo e grades
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.5),
    
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  )

print(pgls_plot)

#ggsave("PGLS_Lys_GenomeSize.pdf", plot = pgls_plot, width = 8, height = 6, dpi = 300)

cat("\nAnalysis complete. Check your working directory for the TSV table and PDF plots.\n")