#Script para analisar os dados de expressão gênica obtidos no GEO a partir
#do id: GSE162689. São dados de biópsia de pâncreas de pacientes com Diabetes
#tipo 1 e pacientes controle não-diabéticos. Vamos fazer uma análise de expressão
#diferencial de genes usando os dados obtidos.

#Chamada de pacotes
library(dplyr)
library(DESeq2)

#Definir pasta de trabalho
setwd("/Users/Usuario/Desktop/Curso ciências ômicas/")

#Carregar count table
count_data <- read.csv(file = "GSE162689/GSE162689_raw_counts_all.txt.gz",sep = "\t")
count_data <- count_data[-1,]

count_data <- count_data %>%
  filter(!is.na(is_133_001),
         !grepl(pattern = "-Mar|-Sep|-Dec",Title))

#Carregar metadados e manipular as tabelas
accession2id <- read.csv("https://raw.githubusercontent.com/csbl-inovausp/RNAseq_DESeq2_R_tutorial/main/data/accession2id.csv")

sra_run_table <- read.csv("https://raw.githubusercontent.com/csbl-inovausp/RNAseq_DESeq2_R_tutorial/main/data/SraRunTable.txt")

proto_pheno <- sra_run_table %>%
  left_join(y = accession2id,by = c("Sample.Name"="GEO_Accession"))

pheno_data <- proto_pheno %>%
  dplyr::select(sample_id,disease_state,tissue)

pheno_data$disease_state <- factor(x = pheno_data$disease_state,
                                   levels = c("non-diabetic","type 1 diabetic"))

pheno_islet <- pheno_data %>%
  filter(tissue=="Islet tissue")

islet_cols <- c(1,which(colnames(count_data) %in% pheno_islet$sample_id))
count_data_islet <- count_data[,islet_cols]

#Criar objeto dds (DESeq2)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data_islet,
                                      colData = pheno_islet,
                                      tidy = T,
                                      design = ~disease_state)

nrow(dds)

keep_10 <- rowSums(DESeq2::counts(dds) >= 10) >= 5
dds_filtrado_final <- dds[keep_10,]

nrow(dds_filtrado_final)

#Plotar PCA das amostras
vsd <- DESeq2::vst(dds_filtrado_final, blind = FALSE)
DESeq2::plotPCA(vsd, intgroup = "disease_state")  


#Plotar heatmap de distância das amostras
sampleDists <- dist(t(assay(vsd)))

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$disease_state, sep = " - " )
colnames(sampleDistMatrix) <- paste( vsd$disease_state, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#Rodar análise de expressão diferencial




