library(caret)

# Preparação de amostras (matrizes não transpostas)
# Alterar variáveis em função do necessário
expr_data = t(highly_variable_lcpm)
phenotype = meta$expr
expr_profiles = cbind(expr_data,phenotype)
expr_profiles = as.data.frame(expr_profiles)
labels = c('prox.-inflam','prox.-prolif','TRU')
expr_profiles$phenotype = factor(expr_profiles$phenotype, levels = 1:3 , labels = labels)
dim(expr_profiles)

# Criação do training set
set.seed(100)

# Preparar as amostras e o tipo
inTrain = createDataPartition(y = t_ex[,101] , p = 0.8, list = F)
trainData = t_ex[inTrain,]
testData = t_ex[-inTrain,]

# Exemplo para random forests num intervalo de 100 colunas

interv = 1:100

luad.rf.100 = caret::train(
  trainData[,1:100],
  trainData$phenotype, 
  method = c('rf'), 
  trControl = trainControl(method = 'boot', number = 5), 
  tuneLength = 5, 
  preProc = c("center", "scale")
)

# Exemplo para MLP para um intervalo de 100 colunas

interv = 1:100

luad.mlp.100 = caret::train(
  trainData[,interv],
  trainData$phenotype, 
  method = 'mlp', 
  trControl = trainControl(method='cv',number='5')
)

# Análise de resultados (alterar em função do nome do modelo)

# Resultados do modelo
luad.rf.100$results

# Utilizamos o modelo de inferência no conjunto de testes
pred.rf.100 = predict(luad.rf.100,testData[,1:100])

# Geramos e analisamos a matriz de confusão
cm.rf.100 = caret::confusionMatrix(pred.rf.100, testData$phenotype)
cm.rf.100
cm.rf.100$overall
cm.rf.100$byClass