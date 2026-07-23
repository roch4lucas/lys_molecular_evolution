##### PGLS: Calliphoridae hábitos alimentares x genes LysBr #####

library(ape)
library(caper)

# 1. Carregamento dos arquivos
arvore <- read.tree("speciestree_calliphoridae.nwk")
tabela <- read.table("calliphoridae_pseudos.tsv", header = TRUE, sep = "\t")

# 2. Transformação da árvore
arvore <- compute.brlen(arvore, method = "Grafen") # *Estima ramos baseados na topologia* --> Método de Grafen (1989)!
arvore <- multi2di(arvore) # Garante que a árvore seja binária (sem politomias)

# 3. Preparação dos dados
tabela$Larval_feeding_habit <- as.factor(tabela$Larval_feeding_habit)

# Inverção da matriz VCV pelo pacote caper
dados_comp <- comparative.data(phy = arvore, 
                               data = tabela, 
                               names.col = "Species", 
                               vcv = TRUE)

# 4. PGLS
# Se o erro persistir com lambda="ML", tente lambda=1 (assume sinal filogenético máximo)
modelo_pgls <- pgls(Pseudogenes ~ Larval_feeding_habit, 
                    data = dados_comp, 
                    lambda = "ML")

summary(modelo_pgls)
