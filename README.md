# Análise de Expressão Diferencial

Este repositório contém o código e a documentação para uma análise de expressão diferencial de dados de adenocarcinoma pulmonar (LUAD). O estudo foi desenvolvido no âmbito da Unidade Curricular de Extração de Conhecimento de Dados Biológicos do Mestrado em Bioinformática, e é resultado do trabalho conjunto de Mariana Oliveira (PG52648), Rui Sousa (PG21019) e Samuel Baptista (PG49130).


## Materiais e Métodos

Utilizamos um dataset público do The Cancer Genome Atlas Program (TCGA), disponibilizado através do NCI Genomic Data Commons (GDC), contendo 600 amostras de pacientes com LUAD. O dataset pode ser encontrado [aqui](https://portal.gdc.cancer.gov/projects/TCGA-LUAD). 

A análise dos dados foi realizada utilizando a ferramenta edgeR.

## Requisitos

Para reproduzir a análise, é recomendada a instalação das seguintes bibliotecas do Bioconductor: 

```R
 if (!require("BiocManager"))
      install.packages("BiocManager")
  
  if (!require("edgeR"))
    BiocManager::install("edgeR")
  
  if (!require("limma"))
    BiocManager::install("limma")
  
  if (!require("Glimma"))
    BiocManager::install(c("Glimma"))
  
  if (!require("gplots"))
    BiocManager::install(c("gplots"))
  
  if (!require("RColorBrewer"))
    BiocManager::install(c("RColorBrewer"))
  
  if (!require("TCGAbiolinks"))
    BiocManager::install("TCGAbiolinks")
  
  if (!require("SummarizedExperiment"))
    BiocManager::install("SummarizedExperiment")
  
  if (!require("biomaRt"))
    BiocManager::install("biomaRt")
  
  if (!require("GSEABase"))
    BiocManager::install(c("GSEABase"))
  
  if (!require("fgsea"))
    BiocManager::install(c("fgsea"))
  
  if (!require("clusterProfiler"))
    BiocManager::install("clusterProfiler")
  
  if (!require("ggplot2"))
    install.packages("ggplot2")
  
  if (!require("dplyr"))
    install.packages("dplyr")
  
  if (!require("annotables")) {
    install.packages("devtools")
    devtools::install_github("stephenturner/annotables")
  }
  
  if (!require("caret"))
    install.packages("caret")
  
  if (!require("EnsDb.Hsapiens.v79"))
    BiocManager::install("EnsDb.Hsapiens.v79")

  if (!require("factoextra"))
    install.packages("factoextra")
  
```

### Transferência e Carregamento do Data Set

O dataset foi transferido do TCGA utilizando o pacote TCGAbiolinks e processado para gerar um objeto SummarizedExperiment. Este objeto contém os dados de expressão genética necessários para a análise.

Para carregar os dados, execute o seguinte comando:

```R
library(TCGAbiolinks)
luad = GDCquery( 
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = 'RNA-Seq',
  workflow.type = "STAR - Counts",
  access = "open"
)
GDCdownload(luad)
luad_data = GDCprepare(luad)
```
