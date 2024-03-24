#
library(TCGAbiolinks)
library(DESeq2)

# Query para chegar ao estudo
luad = GDCquery( 
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = 'RNA-Seq',
  workflow.type = "STAR - Counts",
  access = "open"
)

# Olhar para os resultados / casos
out = getResults(luad)

# O data set tem 600 casos e um tamanho considerávle (2 GB)
# Proponho escolher 210 aleatoriamente

#casos = out$cases[sample(1:600, 210)]

# Rescrevemos a query
luad = GDCquery( 
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = 'RNA-Seq',
  workflow.type = "STAR - Counts",
  access = "open",
  barcode = casos
)

# Buscamos os dados clínicos

clin = GDCquery(
  project = 'TCGA-LUAD',
  data.category = 'Clinical',
  data.type = 'Clinical Supplement',
  data.format = 'BCR Biotab'
)


#GDCdownload(clin)
# GDCdownload(luad)

# Carregamos os dados do estudo e os dados clínicos
luad_data = GDCprepare(luad, summarizedExperiment = TRUE)
clin_data = GDCprepare(clin)

# Carregamos os dados RNASeq 
seqdata = assay(luad_data, 'unstranded')

# Limitar ao tumour stage (alterar se necessário)
data_de = seqdata[,!is.na(luad_data$paper_Tumor.stage)]

# Análise exploratória
dim(data_de)
head(data_de[,1:5])
any(is.na(data_de))

# Calculamos CPM
calccpm = cpm(data_de)

# Aparamos os dados, removendo genes com baixa expressão
# é geralmente aceite a eliminação de genes com CPM inferior a 0.5 em mais do que 2 amostras
thresh = calccpm > 0.5
keep = rowSums(thresh) >= 2
counts.keep = data_de[keep,]
summary(keep)
dim(counts.keep)

# Foram excluídos 28182 genes
dim(data_de) - dim(counts.keep)


# Criar objetos para a análise de expr. diferencial
dgeObj = DGEList(counts.keep)
dgeObj
names(dgeObj)
dgeObj$samples


## distributions - log transform
logcounts <- cpm(dgeObj,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

