# Análise de Expressão Diferencial

Este repositório contém o código e a documentação para uma análise de expressão diferencial de dados de adenocarcinoma pulmonar (LUAD). O estudo foi desenvolvido no âmbito da Unidade Curricular de Extração de Conhecimento de Dados Biológicos do Mestrado em Bioinformática, e é resultado do trabalho conjunto de Mariana Oliveira (PG52648), Rui Sousa (PG21019) e Samuel Baptista (PG49130).


## Materiais e Métodos

Utilizamos um dataset público do The Cancer Genome Atlas Program (TCGA), disponibilizado através do NCI Genomic Data Commons (GDC), contendo 600 amostras de pacientes com LUAD. O dataset pode ser encontrado [aqui](https://portal.gdc.cancer.gov/projects/TCGA-LUAD). 

A análise dos dados foi realizada utilizando a ferramenta edgeR.

## Requisitos

Para reproduzir a análise, é recomendada a instalação das seguintes bibliotecas do Bioconductor: 

```R
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(TCGAbiolinks)
library(SummarizedExperiment)
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
