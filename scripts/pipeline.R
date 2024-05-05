

# Lista de requisitos para poder correr o código deste notebook corretamente

suppressMessages({
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
  
})




# carregamento de bibliotecas
suppressMessages({
  library(edgeR)
  library(limma)
  library(Glimma)
  library(gplots)
  library(RColorBrewer)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(biomaRt)
  library(GSEABase)
  library(fgsea)
  library(clusterProfiler)
  library(ggplot2)
  library(dplyr)
  library(annotables)
  library(EnsDb.Hsapiens.v79)
  library(caret)
  library(factoextra)
})


# Elaboração da query para conectar ao servidor GDC

#luad = GDCquery( 
#  project = "TCGA-LUAD",
#  data.category = "Transcriptome Profiling",
#  experimental.strategy = 'RNA-Seq',
#  workflow.type = "STAR - Counts",
#  access = "open"
#)


# Transferência dos ficheiros

#GDCdownload(luad)


#' Carregamento dos ficheiros num objeto SummarizedExperiment
#' 
#' O parâmetro save permite guardar o histórico num ficheiro compacto
#' tornando a leitura e carregamento numa próxima utilização mais ágil/rápida

#luad_data = GDCprepare(luad,
#                       save = TRUE,
#                       save.filename = 'luad_data_load.rda',
#                       summarizedExperiment = TRUE)

#' No próximo carregamento utilizar apenas o seguinte comando para não saturar o servidor e tornar o processo mais rápido
#' O ficheiro encontra-se aqui: https://ruisousa.me/luad_dataset/luad_data_load.rda

# Alteramos o working  directory para a raiz, para ser mais fácil utilizar as diferentes pastas (modelos, scripts, etc.)
#setwd('../')

# Carregamento do ficheiro
load('luad_data_load.rda')
luad_data = data


# Tamanho do objeto
dim(luad_data)

# Informação sobre o estudo
metadata(luad_data)

# Informação sobre os tipos de dados de RNASeq
names(assays(luad_data))

# Tipos de metadados associados a cada gene
names(rowData(luad_data))

# Tipos de metadados associados a cada amostra
names(colData(luad_data))




# Carregamento de metadados

gender = luad_data$gender

sample_id = substr(luad_data$sample_id, 1, 10)

expr = luad_data$paper_expression_subtype

meta = data.frame(sample_id = sample_id, gender = gender, expr = expr)





##Género

#Verificação da presença de NAs

any(is.na(meta$gender))

# Converte 'gender' para factor com dois níveis
meta$gender = factor(meta$gender, levels = c("male", "female"))

filter = meta$gender != '[Not Available]'

meta = meta[filter,]

#Dimensões da variável

table(meta$gender)


##Subtipos de expressão

#Verificação da presença de NAs

any(is.na(meta$expr))

#Eliminação das colunas correspondentes aos NAs

meta = meta[!is.na(meta$expr),]

filter = meta$expr != '[Not Available]'

meta = meta[filter,]

length(meta$expr)

#Verificação da eliminação dos NAs

any(is.na(meta$expr))

#Verificação de duplicados

# Verificar se há duplicatas na coluna "sample_id"
duplicados <- any(duplicated(meta$sample_id))

duplicados

# Remover as linhas duplicadas com base na coluna "sample_id"
meta = distinct(meta, sample_id, .keep_all = TRUE)

#Dimensões da variável

table(meta$expr)
length(meta$expr)



seqdata_filter = seqdata[,meta$expr]




# Garantir que as dims estão corretas

dim(seqdata_filter)

dim(meta)

# Corrigir nomes e garantir ordem 

colnames(seqdata_filter) = meta$sample_id

all(names(seqdata_filter) == meta$sample_id)




genero = table(luad_data$gender)

colors <- c("pink", "lightblue")  

max_level <- 400

ylim <- c(0, max_level)

barplot(genero, col = colors, ylim = ylim, main = "Distribuição de Sexo", ylab = "Ocorrências")

legend("topright", legend = names(genero), fill = colors)

for (i in 1:length(genero)) {
  text(i, genero[i], labels = genero[i], pos = 3, cex = 0.8, col = "black")
}

#chi_genero = chisq.test(genero)


ages = na.omit(luad_data$paper_Age.at.diagnosis)
ages = ages[ages != '[Not Available]']
ages = as.numeric(as.character(ages))

summary(ages)

# Discretização
breaks = c(40, 50, 60, 70, 80, Inf)
labels = c("40-49", "50-59", "60-69", "70-79", "80 and over")
ages_int = cut(ages, breaks = breaks, labels = labels, include.lowest = TRUE)

# Create bar plot with title and x-axis label
barplot(table(ages_int), main = "Distribuição de idades", xlab = "Intervalos de idades")
abline(h = mean(ages), col = "red", lty = 2, lwd = 1)
legend("topright", legend = paste("Idade média:", round(mean(ages), 2)), col = "red", lty = 2, lwd = 1)




# Combinar as colunas de estado vital e estágio do tumor
vital_stage <- data.frame(vital_status = luad_data$vital_status, tumor_stage = luad_data$paper_Tumor.stage)

# Substituir "Not Available" por NA na variável tumor_stage
vital_stage$tumor_stage <- gsub("Not Available", NA, vital_stage$tumor_stage)

# Remover valores omissos
vital_stage <- vital_stage[complete.cases(vital_stage$tumor_stage), ]

# Criação da tabela de contingência
contingency_table <- table(vital_stage$vital_status, vital_stage$tumor_stage)

# Realização do teste de qui-quadrado
chi_sq_test <- chisq.test(contingency_table)
chi_sq_test

# Representação gráfica
barplot(contingency_table, beside = TRUE, legend.text = TRUE,
        main = "Vital Status vs. Tumor Stage",
        xlab = "Estágio tumor", ylab = "Frequência",
        col = c("lightgreen", "pink"))




# Cálculo CPM
calccpm = cpm(seqdata_filter)

# Remoção de genes com baixa expressão

thresh = calccpm > 0.5
keep = rowSums(thresh) >= 2
counts_keep = seqdata_filter[keep,]

summary(keep)

dim(counts_keep)




# Criação do objeto para a análise de expr. diferencial

dgeObj = DGEList(counts_keep)

names(dgeObj)

head(dgeObj$samples)




## distributions - log transform
logcounts = cpm(dgeObj,log=TRUE)

# Set up the connection to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_dataset <- useDataset('hsapiens_gene_ensembl', mart = ensembl)

# Remove version numbers from Ensembl gene IDs in logcounts
ensembl_ids <- gsub("\\..*", "", rownames(logcounts))
rownames(logcounts) <- ensembl_ids

# Retrieve gene symbols for the Ensembl IDs
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_ids,
                      mart = ensembl)

# Loop through each Ensembl gene ID in seqdata
for (i in seq_along(ensembl_ids)) {
  # Find the index of the Ensembl gene ID in gene_symbols
  idx <- match(ensembl_ids[i], gene_symbols$ensembl_gene_id)
  # If a corresponding gene symbol is found and it's not an empty string, replace the row name
  if (!is.na(idx) && gene_symbols$external_gene_name[idx] != "") {
    rownames(logcounts)[i] <- gene_symbols$external_gene_name[idx]
  }
}



boxplot(logcounts[,1:50], xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots de logCPMs")

boxplot(logcounts[,51:100], xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots de logCPMs")


boxplot(logcounts[,101:150], xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots de logCPMs")


boxplot(logcounts[,151:200], xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots de logCPMs")

boxplot(logcounts[,201:239], xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots de logCPMs")


# Cálculo da variância de cada gene nos dados de contagem logaritmizada
var_genes = apply(logcounts, 1, var)

# Seleção dos 10 genes com maior variabilidade
select_var = names(sort(var_genes, decreasing=TRUE))[1:30]
select_var

# Seleciona as linhas da matriz 'logcounts' com base nos índices fornecidos em 'select_var'
highly_variable_lcpm = logcounts[select_var,]
#head(highly_variable_lcpm)

# Calcula as dimensões (número de linhas e colunas) da matriz 'highly_variable_lcpm'
dim(highly_variable_lcpm)

# Criação do mapa
mypalette = brewer.pal(9,"RdYlBu")
morecols = colorRampPalette(mypalette)

labels = levels(meta$expr)
col.cell1 = c("purple","orange", "green")[meta$expr]
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 10 genes mais variáveis",
          ColSideColors=col.cell1,scale="row",
          margins = c(5, 12))
legend("left", legend=labels, fill=c("purple","orange","green")) # corrigir etiqueta

labels = levels(meta$expr)
col.cell2 = c("purple","orange", "green")[meta$expr][1:50]
heatmap.2(highly_variable_lcpm[,1:50], 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 10 genes mais variáveis entre 50 amostras aleatórias",
          ColSideColors=col.cell2,scale="row",
          margins = c(5, 12))
legend("left", legend=labels, fill=c("purple","orange","green")) # corrigir etiqueta



dgeObj = calcNormFactors(dgeObj)

# demo

plotMD(logcounts, column = 7)
abline(h=0,col="grey")

plotMD(dgeObj, column = 7)
abline(h=0,col="grey")

dgeObj


# Substituir caracteres especiais nos valores de expr
meta_expr <- gsub("prox.-prolif.", "prox_prolif", meta$expr)
meta_expr <- gsub("prox.-inflam", "prox_inflam", meta$expr)

meta_expr = paste(meta$expr)

meta_expr = as.character(meta_expr)

meta_expr


# Definir a variável de design

design = model.matrix(~meta_expr)


head(design)

dgeObj = estimateCommonDisp(dgeObj)

dgeObj$common.dispersion

dgeObj = estimateGLMTrendedDisp(dgeObj)

dgeObj = estimateTagwiseDisp(dgeObj)

plotBCV(dgeObj)


fit = glmFit(dgeObj, design)

head(fit$coefficients)

lrt.BvsL = glmLRT(fit, coef = 2) 

topTags(lrt.BvsL)

results_PP <- as.data.frame(topTags(lrt.BvsL,n = Inf))
results_PP
dim(results_PP)
summary(de <- decideTestsDGE(lrt.BvsL))

detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.BvsL, de.tags=detags)


signif <- -log10(results_PP$FDR)
plot(results_PP$logFC,signif,pch=16)
points(results_PP[detags,"logFC"],-log10(results_PP[detags,"FDR"]),pch=16,col="red")


lrt.BvsL = glmLRT(fit, coef = 3) 

topTags(lrt.BvsL)

results_TRU <- as.data.frame(topTags(lrt.BvsL,n = Inf))
results_TRU
dim(results_TRU)
summary(de <- decideTestsDGE(lrt.BvsL))

detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.BvsL, de.tags=detags)

signif <- -log10(results_TRU$FDR)
plot(results_TRU$logFC,signif,pch=16)
points(results_TRU[detags,"logFC"],-log10(results_TRU[detags,"FDR"]),pch=16,col="red")


# ML



seqdata = as.data.frame(assay(luad, 'unstranded'))
colnames(seqdata) = luad$sample_id
seqdata = seqdata[,expr]

# Cálculo CPM
calccpm = cpm(seqdata)

# Remoção de genes com baixa expressão
thresh = calccpm > 0.5
keep = rowSums(thresh) >= 2
counts_keep = seqdata[keep,]

#summary(keep)
dim(counts_keep)


# Criação do objeto para a análise de expr. diferencial
dgeObj = DGEList(counts_keep)
#names(dgeObj)
#head(dgeObj$samples)

## distributions - log transform
logcounts = cpm(dgeObj,log=TRUE)
#head(logcounts)

suppressMessages(
  library(biomaRt)
)

# Set up the connection to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_dataset <- useDataset('hsapiens_gene_ensembl', mart = ensembl)

# Remove version numbers from Ensembl gene IDs in logcounts
ensembl_ids <- gsub("\\..*", "", rownames(logcounts))
rownames(logcounts) <- ensembl_ids

# Retrieve gene symbols for the Ensembl IDs
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_ids,
                      mart = ensembl)

# Loop through each Ensembl gene ID in seqdata
for (i in seq_along(ensembl_ids)) {
  # Find the index of the Ensembl gene ID in gene_symbols
  idx <- match(ensembl_ids[i], gene_symbols$ensembl_gene_id)
  # If a corresponding gene symbol is found and it's not an empty string, replace the row name
  if (!is.na(idx) && gene_symbols$external_gene_name[idx] != "") {
    rownames(logcounts)[i] <- gene_symbols$external_gene_name[idx]
  }
}

# Cálculo da variância de cada gene nos dados de contagem logaritmizada
var_genes       = apply(logcounts, 1, var)
var_genes.100   = apply(logcounts, 1, var)
var_genes.500   = apply(logcounts, 1, var)
var_genes.1000  = apply(logcounts, 1, var)
var_genes.5000  = apply(logcounts, 1, var)
var_genes.10000 = apply(logcounts, 1, var)

# Seleção dos 10 genes com maior variabilidade
select_var       = names(sort(var_genes, decreasing=TRUE))
select_var.100   = names(sort(var_genes, decreasing=TRUE))[1:100]
select_var.500   = names(sort(var_genes, decreasing=TRUE))[1:500]
select_var.1000  = names(sort(var_genes, decreasing=TRUE))[1:1000]
select_var.5000  = names(sort(var_genes, decreasing=TRUE))[1:5000]
select_var.10000 = names(sort(var_genes, decreasing=TRUE))[1:10000]

# Seleciona as linhas da matriz 'logcounts' com base nos índices fornecidos em 'select_var'
highly_variable_lcpm       = logcounts[select_var,]
highly_variable_lcpm.100   = logcounts[select_var.100,]
highly_variable_lcpm.500   = logcounts[select_var.500,]
highly_variable_lcpm.1000  = logcounts[select_var.5000,]
highly_variable_lcpm.5000  = logcounts[select_var.1000,]
highly_variable_lcpm.10000 = logcounts[select_var.10000,]



# Visualização de mapas para diferentes perfis para ter uma percepção de possíveis resultados

mypalette = RColorBrewer::brewer.pal(9,"RdYlBu")
morecols = colorRampPalette(mypalette)

# 100

labels = levels(expr)
col.cell2 = c("purple","orange", "green")[expr][1:50]
heatmap.2(highly_variable_lcpm.100[,1:50],
          col=rev(morecols(50)),
          trace="column",
          main="Top 100 genes mais variáveis entre 50 amostras aleatórias",
          ColSideColors=col.cell2,scale="row",
          margins = c(5, 12))
legend("left", legend=labels, fill=c("purple","orange","green")) # corrigir etiqueta

# 500

labels = levels(expr)
col.cell2 = c("purple","orange", "green")[expr][1:50]
heatmap.2(highly_variable_lcpm.500[,1:50],
          col=rev(morecols(50)),
          trace="column",
          main="Top 500 genes mais variáveis entre 50 amostras aleatórias",
          ColSideColors=col.cell2,scale="row",
          margins = c(5, 12))
legend("left", legend=labels, fill=c("purple","orange","green")) # corrigir etiqueta

# 1000

labels = levels(expr)
col.cell2 = c("purple","orange", "green")[expr][1:50]
heatmap.2(highly_variable_lcpm.1000[,1:50],
          col=rev(morecols(50)),
          trace="column",
          main="Top 1000 genes mais variáveis entre 50 amostras aleatórias",
          ColSideColors=col.cell2,scale="row",
          margins = c(5, 12))
legend("left", legend=labels, fill=c("purple","orange","green")) # corrigir etiqueta

# 5000

labels = levels(expr)
col.cell2 = c("purple","orange", "green")[expr][1:50]
heatmap.2(highly_variable_lcpm.5000[,1:50],
          col=rev(morecols(50)),
          trace="column",
          main="Top 5000 genes mais variáveis entre 50 amostras aleatórias",
          ColSideColors=col.cell2,scale="row",
          margins = c(5, 12))
legend("left", legend=labels, fill=c("purple","orange","green")) # corrigir etiqueta

# 10000

labels = levels(expr)
col.cell2 = c("purple","orange", "green")[expr][1:50]
heatmap.2(highly_variable_lcpm.10000[,1:50],
          col=rev(morecols(50)),
          trace="column",
          main="Top 10000 genes mais variáveis entre 50 amostras aleatórias",
          ColSideColors=col.cell2,scale="row",
          margins = c(5, 12))
legend("left", legend=labels, fill=c("purple","orange","green")) # corrigir etiqueta




# Preparar as amostras e o tipo
expr_data = t(highly_variable_lcpm)

# Utilizamos a variável dos meta dados para a usar "fresca"
phenotype = expr

expr_profiles = cbind(expr_data,phenotype)
expr_profiles = as.data.frame(expr_profiles)
labels = c('prox.-inflam','prox.-prolif','TRU')
expr_profiles$phenotype = factor(expr_profiles$phenotype, levels = 1:3 , labels = labels)




# Criação do training set
inTrain = createDataPartition(y = expr_profiles$phenotype , p = 0.80, list = F)
trainData = expr_profiles[inTrain,]
testData = expr_profiles[-inTrain,]

dim(trainData)



# Carregamento dos modelos pré-treinados
#setwd('../')

# Random forests
luad.rf.100   = readRDS('./models/random_forests/luad.rf.100.rds')
luad.rf.500   = readRDS('./models/random_forests/luad.rf.500.rds')
luad.rf.1000  = readRDS('./models/random_forests/luad.rf.1000.rds')
luad.rf.5000  = readRDS('./models/random_forests/luad.rf.5000.rds')
luad.rf.10000 = readRDS('./models/random_forests/luad.rf.10000.rds')
luad.rf.19529 = readRDS('./models/random_forests/luad.rf.19529.rds')

# MLP
luad.mlp.100   = readRDS('./models/multilayer_perceptron/luad.mlp.100.rds')
luad.mlp.500   = readRDS('./models/multilayer_perceptron/luad.mlp.500.rds')
luad.mlp.1000  = readRDS('./models/multilayer_perceptron/luad.mlp.1000.rds')
luad.mlp.5000  = readRDS('./models/multilayer_perceptron/luad.mlp.5000.rds')
luad.mlp.10000 = readRDS('./models/multilayer_perceptron/luad.mlp.10000.rds')
luad.mlp.19529 = readRDS('./models/multilayer_perceptron/luad.mlp.19529.rds')

dim(expr_profiles)






# Resultados do modelo
luad.rf.100$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.rf.100 = predict(luad.rf.100,testData[,1:100])

# Geramos e analisamos a matriz de confusão
cm.rf.100 = caret::confusionMatrix(pred.rf.100, testData$phenotype)
cm.rf.100
cm.rf.100$overall
cm.rf.100$byClass



# Resultados do modelo
luad.rf.500$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.rf.500 = predict(luad.rf.500,testData[,1:500])

# Geramos e analisamos a matriz de confusão
cm.rf.500 = caret::confusionMatrix(pred.rf.500, testData$phenotype)
cm.rf.500$overall
cm.rf.500$byClass



# Resultados do modelo
luad.rf.1000$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.rf.1000 = predict(luad.rf.1000,testData[,1:1000])

# Geramos e analisamos a matriz de confusão
cm.rf.1000 = caret::confusionMatrix(pred.rf.1000, testData$phenotype)
cm.rf.1000
cm.rf.1000$overall
cm.rf.1000$byClass



# Resultados do modelo
luad.rf.5000$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.rf.5000 = predict(luad.rf.5000,testData[,1:5000])

# Geramos e analisamos a matriz de confusão
cm.rf.5000 = caret::confusionMatrix(pred.rf.5000, testData$phenotype)
cm.rf.5000$overall
cm.rf.5000$byClass



# Resultados do modelo
luad.rf.10000$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.rf.10000 = predict(luad.rf.10000,testData[,1:10000])

# Geramos e analisamos a matriz de confusão
cm.rf.10000 = caret::confusionMatrix(pred.rf.10000, testData$phenotype)
cm.rf.10000
cm.rf.10000$overall
cm.rf.10000$byClass



# Resultados do modelo
luad.rf.19529$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.rf.19529 = predict(luad.rf.19529,testData[,1:19529])

# Geramos e analisamos a matriz de confusão
cm.rf.19529 = caret::confusionMatrix(pred.rf.19529, testData$phenotype)
cm.rf.19529$overall
cm.rf.19529$byClass

#### Multi-Layer Perceptron (MLP)


# Resultados do modelo
luad.mlp.100$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.mlp.100 = predict(luad.mlp.100,testData[,interv])

# Geramos e analisamos a matriz de confusão
cm.mlp.100 = caret::confusionMatrix(pred.mlp.100, testData$phenotype)
cm.mlp.100
cm.mlp.100$overall
cm.mlp.100$byClass



# Resultados do modelo
luad.mlp.500$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.mlp.500 = predict(luad.mlp.500,testData[,1:500])

# Geramos e analisamos a matriz de confusão
cm.mlp.500 = caret::confusionMatrix(pred.mlp.500, testData$phenotype)
cm.mlp.500
cm.mlp.500$overall
cm.mlp.500$byClass



# Resultados do modelo
luad.mlp.1000$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.mlp.1000 = predict(luad.mlp.1000,testData[,1:1000])

# Geramos e analisamos a matriz de confusão
cm.mlp.1000 = caret::confusionMatrix(pred.mlp.1000, testData$phenotype)
cm.mlp.1000
cm.mlp.1000$overall
cm.mlp.1000$byClass





# Resultados do modelo
luad.mlp.5000$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.mlp.5000 = predict(luad.mlp.5000,testData[,1:5000])

# Geramos e analisamos a matriz de confusão
cm.mlp.5000 = caret::confusionMatrix(pred.mlp.5000, testData$phenotype)
cm.mlp.5000
cm.mlp.5000$overall
cm.mlp.5000$byClass



# Resultados do modelo
luad.mlp.10000$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.mlp.10000 = predict(luad.mlp.10000,testData[,1:10000])

# Geramos e analisamos a matriz de confusão
cm.mlp.10000 = caret::confusionMatrix(pred.mlp.10000, testData$phenotype)
cm.mlp.10000
cm.mlp.10000$overall
cm.mlp.10000$byClass




# Resultados do modelo
luad.mlp.19529$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.mlp.19529 = predict(luad.mlp.19529,testData[,1:19529])

# Geramos e analisamos a matriz de confusão
cm.mlp.19529 = caret::confusionMatrix(pred.mlp.19529, testData$phenotype)
cm.mlp.19529
cm.mlp.19529$overall
cm.mlp.19529$byClass
