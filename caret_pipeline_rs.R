#libs
library(caret)
library(RSNNS)
library(party)
library(randomForest)

# Criação do training set
set.seed(100)

# Preparar as amostras e o tipo
t_hvl = t(highly_variable_lcpm)
t_ex = cbind(t_hvl,meta$expr)
t_ex = as.data.frame(t_ex)
# Garantir que está tudo OK
all(rownames(t_hvl) == meta$sample_id)

inTrain = createDataPartition(y = t_ex[,101] , p = 0.8, list = F)
trainData = t_ex[inTrain,]
testData = t_ex[-inTrain,]

# Usando ridge reg, conventional rf e conditional inference rf, MLP
# Baseado na seguinte lógica: https://chat.openai.com/share/adf05b1e-265c-4b62-bc2f-4cdd51f78b7f

models = c(
  'cforest',
  'rf',
  'mlp'
)

luad_models = caret::train(t_ex[,1:100],t_ex[,101], method = 'rf')
preds_luad = predict(luad_models, testData[,1:100])
preds_luad
confusionMatrix(preds_luad, testData[,101])