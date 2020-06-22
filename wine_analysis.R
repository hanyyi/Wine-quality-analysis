setwd("E:/first year phd/STA 232A/final project/Project_Data_with_Description/explotary")
library(psych)
library(ggplot2)
#library(devtools)
#install_github("vqv/ggbiplot", force = TRUE)
library(ggbiplot)
library(FactoMineR)
library(FactoMineR)
library(psych)
library(cluster) 
library(fpc)
library(factoextra)
library(car)
library(corrplot)
library(boot)
library(ISLR)
library(leaps)
library(glmnet)
library(pls)
library(e1071)
library(tree)
library(randomForest)
library(car)
library(ROCR)
library(MASS)


data_wine = read.csv("winequality-red.csv", header = T, sep = ";")
head(data_wine)
dim(data_wine)
str(data_wine)
summary(data_wine)
pairs.panels(data_wine)
data_wine = scale(data_wine)
data_wine = as.data.frame(data_wine)

cor(data_wine)
png("pairs.png",width = 4500, height = 3000, res = 72 * 5)
pairs.panels(data_wine)
dev.off()



#PCA
par(mfrow=c(1,1))
data_x=data_wine[,-12]
attach(data_wine)
winepca = princomp(~., data = data_x, cor = TRUE)
summary(winepca, loadings = TRUE)
screeplot(winepca, type = "lines")
biplot(winepca)

png("biplot.png",width = 1500, height = 900, res = 72 * 2)
biplot(winepca, expand = 13, xlim = c(-1.5, 1.5), ylim = c(-1, 1))
dev.off()

PCA(scale(data_x))


fa.parallel(data_x, fa = "pc", n.iter = 100, 
            show.legend = FALSE, main = "Scree plot with parallel analysis")

mydata.pr1 = principal(data_x, nfactors = 4, rotate = "none")
mydata.pr1

#kmeans
data_x_scale = scale(data_x)
#decide the optimal number of clusters 3
fviz_nbclust(data_x_scale, kmeans, method = "wss") + 
  geom_vline(xintercept = 3, linetype = 5,col = "darkred")

k.means = kmeans(data_x_scale, 3, nstart = 25)
k.means

k.means$cluster

clusplot(data_x_scale, k.means$cluster, main = "2D representation of the Cluster", 
         color = TRUE, shade = TRUE,
         labels = 2, lines = 0)

fviz_cluster(object = k.means,
             data = data_x_scale,
             ellipse.type = "norm",
             geom = "point",
             palette = "jco",
             main = "",
             ggtheme = theme_minimal())

data_x_scale=as.data.frame(data_x_scale)
data_x_scale$cluster = k.means$cluster
ggplot(data_x_scale, aes(pH,alcohol,color=as.factor(cluster)))+geom_point()
ggplot(data_x_scale,aes(pH,sulphates,color=as.factor(cluster)))+geom_point()
ggplot(data_x_scale,aes(pH,total.sulfur.dioxide,color=as.factor(cluster)))+geom_point()
ggplot(data_x_scale,aes(alcohol,sulphates,color=as.factor(cluster)))+geom_point()
ggplot(data_x_scale,aes(alcohol,total.sulfur.dioxide,color=as.factor(cluster)))+geom_point()
ggplot(data_x_scale,aes(sulphates,total.sulfur.dioxide,color=as.factor(cluster)))+geom_point()

#########################################################################################
#################################Model Build###################################
#linear model#
set.seed(2)
indexes = sample(1:nrow(data_wine), size=0.75*nrow(data_wine))
Train75 <- data_wine[indexes,]
Test50 <- data_wine[-indexes,]

attach(Train75)
linear_model = lm(quality~., data=Train75)
summary(linear_model)

vif(linear_model)
corrplot(cor(data_x), method = "number")

mean((Test50[,12] - predict(linear_model, Test50[,-12]))^2)
detach(Train75)

#glm_fit = glm(data_wine$quality ~., data=data_wine)
#cv.err = cv.glm(data_wine, glm_fit)



lm.forward.AIC = step(linear_model, k=2,direction="forward",test="F")
summary(lm.forward.AIC)
mean((Test50[,12] - predict(lm.forward.AIC, Test50[,-12]))^2)

lm.backward.AIC = step(linear_model, k=2,direction="backward",test="F")
mean((Test50[,12] - predict(lm.backward.AIC, Test50[,-12]))^2)

lm.both.AIC = step(linear_model,k=2,direction="both",test="F")
mean((Test50[,12] - predict(lm.both.AIC, Test50[,-12]))^2)


###################################################################
data_wine = as.data.frame(scale(read.csv("winequality-red.csv", header = T, sep = ";")))
linear_model_full = lm(quality~. , data = data_wine)
summary(linear_model_full)
bc = boxcox(linear_model_full, lambda=seq(0, 2, by=0.1))
lambda = bc$x[which.max(bc$y)]
y = data_wine[,12]
ylam = (y^lambda - 1)/lambda

data_bc = cbind(data_wine[,-12],ylam)
names(data_bc)[12] = "quality"

linear_model_full_bc = lm(quality~., data=data_bc)
summary(linear_model_full_bc)
plot(linear_model_full_bc)


par(mfrow = c(1, 1))
plot(linear_model_full$residuals)
vif(linear_model_full)


lm.backward.AIC = step(linear_model_full, k=2,direction="backward",test="F")

n = dim(data_wine)[1]
lm.backward.BIC = step(linear_model_full, k=log(n),direction="both",test="F")
par(mfrow = c(2, 2))
plot(lm.backward.BIC)
par(opar)


lm.backward = regsubsets(quality~., data = data_wine, nvmax=11, method = "backward")
summary(lm.backward)
coef(lm.backward, 7)
summary(lm.backward)$bic


lm.summary.backward=summary(lm.backward)
plot(lm.summary.backward$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")
qplot(y=lm.summary.backward$bic,x=c(1:length(lm.summary.backward$bic)),geom="point")
BIC=as.data.frame(cbind(t(t(c(1:length(lm.summary.backward$bic)))),t(t(lm.summary.backward$bic)),t(t(c(0,0,0,0,0,1,0,0,0,0,0)))))
names(BIC) =c("number_of_predict_random_variables","value_of_BIC","color")
ggplot(data = BIC, mapping = aes(x = number_of_predict_random_variables, y = value_of_BIC)) + geom_point(size = 3,color="green")+geom_point(data=BIC[BIC$color == 1,],color="red",size=3)+labs(x="number of predictors",y="value of BIC")

points(which.min(lm.summary.backward$bic), lm.summary.backward$bic[which.min(lm.summary.backward$bic)],
       col = "red", cex = 2, pch = 20)
plot(lm.backward, scale = "bic")


#################################33train and test data###################
set.seed(1)
train = sample(c(TRUE, FALSE), nrow(data_wine), rep = TRUE)
test = !train
x=data_wine[,-12]
y=data_wine[,12]

set.seed(1)
train = sample(1:nrow(x), nrow(x)/2)
test = (-train)
y.test = y[test]

lm.best.train = regsubsets(quality~., data = data_wine[train,], nvmax=11)
test.mat = model.matrix(quality~., data = data_wine[test,])

val.errors = rep(NA, 11)
for (i in 1:11) {
  coefi = coef(lm.best.train, id = i)
  pred = test.mat[, names(coefi)] %*% coefi
  val.errors[i] = mean((data_wine$quality[test] - pred)^2)
}

val.errors
which.min(val.errors)

coef(lm.backward, 8)

##################################k-fold#################
k = 10
set.seed(1)
folds = sample(1:k, nrow(data_wine), replace = TRUE)
cv.errors = matrix(NA, k, 11, dimnames = list(NULL, paste(1:11)))

predict.regsubsets = function(object, newdata, id, ...){
  form = as.formula(object$call[[2]])
  mat = model.matrix(form, newdata)
  coefi = coef(object, id=id)
  xvars = names(coefi)
  mat[,xvars] %*% coefi
}

for (j in 1:k) {
  best.fit = regsubsets(quality~., data = data_wine[folds != j,], nvmax = 11)
  for (i in 1:11) {
    pred = predict(best.fit, data_wine[folds == j,], id = i)
    cv.errors[j,i] = mean((data_wine$quality[folds == j ] - pred)^2)
  }
}
mean.cv.errors = apply(cv.errors, 2, mean)
mean.cv.errors
plot(mean.cv.errors, type="b")
which.min(mean.cv.errors)

coef(lm.backward, 6)
#################################### Ridge Regression and LASSO######
x = model.matrix(quality~., data_wine)[,-1]
y = data_wine$quality
grid = 10^seq(10, -2, length = 100)
ridge.mod = glmnet(x, y, alpha = 0, lambda = grid)

set.seed(1)
train = sample(1:nrow(x), nrow(x)/2)
test = (-train)
y.test = y[test]

set.seed(1)
cv.out = cv.glmnet(x[train,], y[train], alpha = 0)
plot(cv.out)
bestlamda=cv.out$lambda.min
bestlamda

ridge.mod = glmnet(x[train,], y[train], alpha = 0, lambda = grid, thresh = 1e-12)
ridge.pred = predict(ridge.mod, s=bestlamda, newx = x[test,])
mean((ridge.pred-y.test)^2) #0.3954

out = glmnet(x, y, alpha = 0)
predict(out, type="coefficients", s=bestlamda) # last ridge regression model


##########LASSO########
lasso.mod = glmnet(x[train,], y[train], alpha = 1, lambda=grid)
plot(lasso.mod)

set.seed(1)
cv.out.lasso = cv.glmnet(x[train,], y[train], alpha=1)
plot(cv.out.lasso)
bestlamda.lasso = cv.out.lasso$lambda.min
lasso.pred = predict(lasso.mod, s=bestlamda.lasso, newx =x[test,])
mean((lasso.pred - y.test)^2) #0.3947

out.lasso = glmnet(x,y,alpha = 1,lambda = grid)
lasso.coef = predict(out.lasso, type="coefficients", s=bestlamda.lasso)
lasso.coef



############################PCR#################################
set.seed(1)
pcr.fit = pcr(quality~., data=data_wine,subset=train, scale=TRUE, validation="CV")
summary(pcr.fit)

par(mfrow=c(1,1))
validationplot(pcr.fit,val.type = "MSEP")


pcr.pred = predict(pcr.fit,x[test,],ncomp = 3)
mean((pcr.pred-y.test)^2)

pcr.fit = pcr(y~x, scale=TRUE, ncomp = 3)
summary(pcr.fit)

##########################svm#############
#set.seed(100)
#index <- sample(2,nrow(cats),replace = TRUE,prob=c(3/4,1/4))
data_wine = read.csv("winequality-red.csv", header = T, sep = ";")
set.seed(100)
train = sample(1:nrow(x), nrow(x)*3/4)
test = (-train)
y.test = y[test]


dat = data.frame(x=data_wine[train,-12],y=as.factor(data_wine[train,12]))
out = svm(y~., data=dat, kernel="radial",cost = 10)
summary(out)
table(out$fitted, dat$y)

dat.test = data.frame(x=data_wine[test,-12],y=as.factor(data_wine[test,12]))
pred.test = predict(out, newdata = dat.test)
table_test = table(pred.test, dat.test$y)
sum(diag(table_test))/sum(table_test)

factor_wine = data_wine
factor_wine[which(factor_wine[,12] >= 7), 12] = "high"
factor_wine[which(factor_wine[,12]<=5),12]= "low"
factor_wine[which(factor_wine[,12] == 6),12] = "medium"
factor_wine[,12]


dat.factor.train = data.frame(x=factor_wine[train,-12],y=as.factor(factor_wine[train,12]))
out.factor = svm(y~., data = dat.factor.train, kernel="radial",type="C-classification" ,gamma=0.1,cost=10)
summary(out.factor)
table(out.factor$fitted, dat.factor.train$y)

dat.factor.test = data.frame(x=factor_wine[test,-12],y=as.factor(factor_wine[test,12]))
fitted = attributes(predict(out.factor, dat.factor.train,decision.values = TRUE))




pred.factor.test = predict(object=out.factor,newdata = dat.factor.test)
table_factor_test = table(pred.factor.test, dat.factor.test$y)
sum(diag(table_factor_test))/sum(table_factor_test)

modelroc = roc(dat.factor.test$y,pred.factor.test)

roc_plot=function(test_data_y,pred_data_y,flag,notflage){
  test_data_y=as.vector(test_data_y)
  test_data_y[which(test_data_y!=flag)]=notflage
  test_data_y=as.factor(test_data_y)
  
  pred_data_y=as.vector(pred_data_y)
  pred_data_y[which(pred_data_y!=flag)]=notflage
  pred_data_y=as.factor(pred_data_y)
  rpc_test1=roc(test_data_y,as.ordered(pred_data_y))
  return(rpc_test1)
}

notlow = roc_plot(test_data_y = dat.factor.test[,12],pred_data_y = pred.factor.test,flag = "low",notflage = "notlow")

plot(notlow)

notmedium = roc_plot(test_data_y = dat.factor.test[,12],pred_data_y = pred.factor.test,flag = "medium",notflage = "notmedium")
plot(notmedium)

nothigh = roc_plot(test_data_y = dat.factor.test[,12],pred_data_y = pred.factor.test,flag = "high",notflage = "nothigh")
plot(nothigh)

plot(notlow,col="green")
par(mfrow=c(1,1),new = T)
plot(notmedium,col="red")
par(mfrow=c(1,1),new = T)
plot(nothigh,col="blue")
legend('bottomright', c("low","med",'high'),lty=c(1, 2), pch=c(15, 17), col=c("green", 'red',"blue"),cex = 0.8)
#tune.svm(y ~., data =dat.factor.train, gamma = 10^(-100:-1), cost = 10^(0:3))

###############################random forests################################
redwineTree <- tree(dat.factor.train[,12] ~. , data=dat.factor.train[,-12], method="class")
summary(redwineTree)
plot(redwineTree)
text(redwineTree, pretty=0, cex=0.6)
misclass.tree(redwineTree, detail=T)
Treefit1 <- predict(redwineTree, dat.factor.test[,-12], type="class")
table_tree = table(Treefit1, dat.factor.test[,12])
sum(diag(table_tree))/sum(table_tree)


set.seed(3)
cv.redwine.tree=cv.tree(redwineTree,FUN = prune.misclass)
#names(cv.redwine.tree)
cv.redwine.tree

par(mfrow=c(1,2))
plot(cv.redwine.tree$size,cv.redwine.tree$dev,type="b")
plot(cv.redwine.tree$k,cv.redwine.tree$dev,type="b")

prune.redwine.tree = prune.misclass(redwineTree, best=10)
par(mfrow=c(1,1))
plot(prune.redwine.tree)
text(prune.redwine.tree, pretty=0)

prune.pre.tree = predict(prune.redwine.tree,dat.factor.test[,-12], type="class")
table_prune = table(prune.pre.tree, dat.factor.test[,12])
sum(diag(table_prune))/sum(table_prune)







set.seed(1)
rf_wine <- randomForest(dat.factor.train[,12] ~ . , data=dat.factor.train[,-12], mtry=4, importance=T,ntree=500)
plot(c(1:length(rf_wine$err.rate)),rf_wine$err.rate)


rf_predict <- predict(rf_wine, newdata=dat.factor.test[,-12], type="class")


notlow_tree = roc_plot(test_data_y = dat.factor.test[,12],pred_data_y = rf_predict ,flag = "low",notflage = "notlow")
notmedium_tree = roc_plot(test_data_y = dat.factor.test[,12],pred_data_y = rf_predict ,flag = "medium",notflage = "notmedium")
nothigh_tree = roc_plot(test_data_y = dat.factor.test[,12],pred_data_y = rf_predict ,flag = "high",notflage = "nothigh")

plot(notlow_tree,col="green")
par(mfrow=c(1,1),new = T)
plot(notmedium_tree,col="red")
par(mfrow=c(1,1),new = T)
plot(nothigh_tree,col="blue")


table_rf = table(rf_predict, dat.factor.test[,12])
sum(diag(table_rf))/sum(table_rf)
importance(rf_wine)

plot(rf_wine,  main="", cex=0.8)
varImpPlot(rf_wine,  main="", cex=0.8)




boxplot.stats(data_wine[,12])$out
