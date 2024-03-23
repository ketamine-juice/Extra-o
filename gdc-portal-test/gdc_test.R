# Carregar libs
library(TCGAbiolinks)
library(DESeq2)

# https://rpubs.com/tiagochst/TCGAbiolinks_to_DESEq2

# sugestão: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1708

gdcprojects = getGDCprojects()


BRCA = GDCquery( 
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = 'RNA-Seq',
  workflow.type = "STAR - Counts",
  access = "open",
  barcode = "TCGA-E9-A1RH-01A-21R-A169-07"
          )

out = getResults(BRCA)

GDCdownload(BRCA)


data = GDCprepare(BRCA, summarizedExperiment = TRUE)

metadata(data)

BiocManager::install("curatedTCGAData")
                     