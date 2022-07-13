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

accession2id <- read.csv("https://raw.githubusercontent.com/csbl-inovausp/RNAseq_DESeq2_R_tutorial/main/accession2id.csv")

sra_run_table <- read.csv("https://raw.githubusercontent.com/csbl-inovausp/RNAseq_DESeq2_R_tutorial/main/SraRunTable.txt")

proto_pheno <- sra_run_table %>%
  left_join(y = accession2id,by = c("Sample.Name"="GEO_Accession"))

pheno_data <- proto_pheno %>%
  dplyr::select(sample_id,disease_state)
  
  