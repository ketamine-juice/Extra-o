# datasets

data(iris)
dim(iris)
class(iris$Species)
names(iris)

library(MASS)
data(cpus)
dim(cpus)
names(cpus)
class(cpus$perf)

# error metrics

pecc = function(obs,pred) 	sum(obs==pred)/length(obs)
rmse = function(obs, pred) sqrt(mean((obs-pred)^2))
mad = function(obs, pred) mean(abs(obs-pred))

# train/test sets

set.seed(12345)
ind = sample(2, nrow(iris), replace=TRUE, prob=c(0.7, 0.3))
trainData = iris[ind==1,]
testData = iris[ind==2,]
dim(trainData)
dim(testData)
table(trainData$Species)
table(testData$Species)

set.seed(12345)
ordem = sample(nrow(cpus))
tam_treino = 2/3 * nrow(cpus)
ind_tr = ordem[1:tam_treino]
ind_ts = ordem[(tam_treino+1):nrow(cpus)]
cpuTr = cpus[ind_tr,]
cpuTs = cpus[ind_ts,]
dim(cpuTr)
dim(cpuTs)
mean(cpuTr$perf)
mean(cpuTs$perf)
names(cpus)
head(cpus)

# KNN
library(class)
knn_pred = knn(trainData[,1:4], testData[,1:4], trainData$Species)
knn_pred

t = table(knn_pred, testData$Species)
t
pecc(knn_pred, testData$Species)
vp_versicolor = t[2,2]
vn_versicolor = t[1,1] + t[3,3]
fp_versicolor = t[2,1] + t[2,3]
fn_versicolor = t[1,2] + t[3,2]
sensib_versicolor = vp_versicolor / (vp_versicolor + fn_versicolor)
sensib_versicolor
especif_versicolor = vn_versicolor / (vn_versicolor + fp_versicolor)
especif_versicolor

# naive bayes
library(e1071)
model = naiveBayes(Species ~ ., trainData)
nb_pred = predict(model, testData)
nb_pred
table(nb_pred, testData$Species)
sum(nb_pred==testData$Species)/length(testData$Species)
pecc(testData$Species, nb_pred)

# trees
library(party)
iris_ctree = ctree(Species ~., data=trainData)
print(iris_ctree)
plot(iris_ctree)

testPred = predict(iris_ctree, testData)
testPred[1]
table(testPred, testData$Species)  
pecc(testData$Species, testPred)
#sum(testPred==testData$Species)/length(testData$Species)

formula <- Species ~ Sepal.Length + Sepal.Width + Petal.Length
iris_ctree2 <- ctree(formula , data=trainData)
plot(iris_ctree2)
testPred2 <- predict(iris_ctree2, testData)
table(testPred2, testData$Species)        
#sum(testPred2==testData$Species)/length(testData$Species)
pecc(testData$Species, testPred2)

library(tree)
tree1 <- tree(Species ~ Sepal.Width + Petal.Width, data=iris)
plot(tree1)
text(tree1)

plot(iris$Petal.Width,iris$Sepal.Width,pch=19, col=as.numeric(iris$Species))
partition.tree(tree1,label="Species",add=TRUE)
legend(2,4.5,legend=unique(iris$Species), col=unique(as.numeric(iris$Species)),pch=19)
val_prev = predict(tree1,iris, type="class")

length(iris$Species) - sum(val_prev == iris$Species)
plot(jitter(iris$Petal.Width),jitter(iris$Sepal.Width),pch=19, col=as.numeric(iris$Species))
partition.tree(tree1,label="Species",add=TRUE)
legend(2,4.5,legend=unique(iris$Species), col=unique(as.numeric(iris$Species)),pch=19)

tree2 = prune.tree(tree1, best = 3)
plot(tree2)
text(tree2) 
val_prev = predict(tree2,iris, type="class")
sum(val_prev != iris$Species)

library(rpart)
arvreg = rpart(perf ~ syct + mmin + mmax + chmax + chmin, data = cpuTr)
plot(arvreg)
text(arvreg)
val_prev = predict(arvreg, cpuTs)
val_prev
rmse(val_prev, cpuTs$perf)
mad(val_prev, cpuTs$perf)
