# Inspirado nisto: https://www.sciencedirect.com/science/article/pii/S0307904X23002809

library(caret)
library(RSNNS)
library(party)
library(randomForest)

set.seed(100)

inTrain = createDataPartition(y = logcounts , p = 0.8, list = F)
trainData = logcounts[inTrain]
testData = logcounts[-inTrain]

# Usando ridge reg, conventional rf e conditional inference rf, MLP
# Baseado na seguinte lógica: https://chat.openai.com/share/adf05b1e-265c-4b62-bc2f-4cdd51f78b7f

models = c(
'ORFridge',
'cforest',
'rf',
'mlp'
)

luad_models = train(logcounts[,1:4],logcounts[,5], method = models)

preds_luad = predict(luad_models, testData[,1:4])

# agora falta testar o snippet
