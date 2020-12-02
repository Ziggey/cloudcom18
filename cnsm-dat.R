library(glmnet)
library(caret)
library(dplyr)
library(ggplot2)
library(leaps)

x = read.delim("periodicx_nc.dat", header = TRUE, sep = ",")
y = read.delim("periodicy_nc.dat", header = TRUE, sep = ",")
tcp = read.delim("periodicoad.dat", header = TRUE, sep = ",")

########## Histogram of TCPSCK
qplot(tcp, geom = "histogram")

ggplot(data=tcp, aes(tcp$X19)) + 
  geom_histogram(breaks=seq(19, 110, by=1), 
                 col="red", 
                 fill="green", 
                 alpha = .2) +
  ylim(c(0,2000))

### Combine x and y; rename the response variable
xy = cbind(y,x)
names(xy)[1]<-"dispframe"

head(xy[, 1:5])
table(is.na(xy))
#################################################
min(y)
max(y)
max(tcp)
min(tcp)

#### lets filter based on tcp values between 24 & 72 where counts are more than 500
tcp_reduced = filter(xy, X19 >= 24 & X19 <= 72)
min(tcp_reduced$X19)
max(tcp_reduced$X19)

################ Subset selection on the 

?regsubsets
######## Linear Regression #### 
lr = lm(dispframe~., data = xy)
summary(lr)
names(tcp)
names(xy)
names(xy[237])
summary(lr)



write.csv(x, "xfile.csv", row.names = FALSE)
graph_range = range(0, tcp, y)

plot(y[1:100,], type = "o", col ="red", ylim = graph_range, axes = FALSE, ann = FALSE)
box()
lines(tcp[1:100,], type = "o", pch=22, lty = 2, col = "blue")
str(y)



#this plot showed the flunctuations between the tcp socket count and the dispframe rate
plot(tcp[1:25000,], type="l", ylim=c(0,600), col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(y[1:25000,], type="l", col = "red")
legend("topright",legend=c("Load - TCPSCK", "Service Metric y"), col=c("blue", "red"),lty=1.2,cex=0.6)

pdf(file = "flunctuation.pdf")


##############################

2^6

########
attach(train)
par(mfrow = c(1, 2))
fit_ridge = glmnet(train, train$dispframe, alpha = 0)
plot(fit_ridge)
plot(fit_ridge, xvar = "lambda", label = TRUE)


###EN#
lambda.grid = 10^seq(2,-2,length=50)
alpha.grid = seq(0,1,length =10)
searchGrd = expand.grid(.alpha=alpha.grid, .lambda=lambda.grid)


# MODEL COMPARISM
##########################################################
## Data partition
set.seed(102)
id = sample(2, nrow(xy), replace = T, prob = c(0.7, 0.3))
xy_train = xy[id == 1,]
xy_test = xy[id == 2,]
#################Custom COntrol Parameters#################
custom = trainControl(method = "repeatedcv",
                      number = 10,
                      repeats = 5,
                      verboseIter = T)
##################### Linear Model ########################
lm = train(dispframe~.,
           xy_train,
           method = 'lm',
           trControl = custom)
lm$results
lm
summary(lm)
##################### Ridge Regression ####################

ridge = train(dispframe ~ .,
              xy_train,
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = 0,
                                     lambda = seq(0.000001,10, length = 5)),
              trControl = custom)
plot(ridge)
ridge
plot(ridge$finalModel, xvar = "lambda", label = T)
plot(ridge$finalModel, xvar = "dev", label = T)
plot(varImp(ridge, scale = F))

################### Lasso Regression #######################

lasso = train(dispframe ~ .,
              xy_train,
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = 1,
                                     lambda = seq(0.0001, 10, length = 5)),
              trControl = custom)
lasso
plot(lasso$finalModel, xvar = "lambda", label = T)
plot(lasso$finalModel, xvar = "dev", label = T)
varImp(lasso)

################# Elastic Net Regression #######################
en = train(dispframe ~ .,
           xy_train,
           method = 'glmnet',
           tuneGrid = expand.grid(alpha = seq(0,1, length = 10),
                                  lambda = seq(0.0001, 10, length = 5)),
           trControl = custom)
plot(en)
plot(en$finalModel, xvar = "lambda", label = T)
plot(en$finalModel, xvar = "dev", label = T)
en$bestTune
best = en$finalModel
coef(best, s = en$bestTune$lambda)

############## Compare Models #################################
model_list = list(LinearModel = lm, Ridge = ridge, Lasso = lasso, ElasticNet = en)
res = resamples(model_list)
summary(res)
xyplot(res, metric = 'RMSE')

############# Predict on the EN Model 
#RMSE for train
p1 = predict(en, xy_train)
rmse_train = sqrt(mean((xy_train$dispframe-p1)^2))
rmse_train
#RMSE for test
p2 = predict(en, xy_test)
rmse_test = sqrt(mean((xy_test$dispframe-p2)^2))
rmse_test

min(xy$X19)
max(xy$X19)
mean(xy$X19)

######## 1000 Time Points Sampling ###
xy_train_la = xy_train[1:1000,]
la_lm = lm(dispframe~X19, data = xy_train_la)
summary(la_lm)
sum(xy$X19>39)
sum(xy$X19<39)
33768 + 16484
51042 - 50252
sum(xy$X19==39)
sum(xy$X19==70)

##### LR prediction using the test set
lr_predic = predict(lm, xy_test)
lr_predic

plot(xy_test$dispframe, type = "l", lty = 1.8, col = "blue")
lines(lr_predic, type = "l", col = "red")
summary(lr_predic)
names(lr_predic)

ridge$lambda
ridge
plot(en)
###### ISLR method
set.seed(234)
x1= xy[1:10000,1:265]

rm(x_islr)
#make x and y
x_islr = model.matrix(dispframe~., x1)[,-1]
y_islr = x1$dispframe
#Ridge
grid = 10^seq(10,-2,length=100)
ridge.mod = glmnet(x_islr,y_islr, alpha = 0, lambda = grid)
lasso.mod = glmnet(x_islr,y_islr, alpha = 1, lambda = grid)
dim(coef(ridge.mod))

ridge.mod$lambda[60]


############## Load Effect
a1 = x1 %>%
  select(dispframe,X1.12,X0.50,X19,X0.122,X29.5,X1664,X0.3,X2003,X0.91,X0.88,X1630) %>%
  filter(X19 <= 30)

a2 = x1 %>%
select(dispframe,X1.12,X0.50,X19,X0.122,X29.5,X1664,X0.3,X2003,X0.91,X0.88,X1630) %>%
  filter(X19 > 30)
varImp(en)
##################
la_en = train(dispframe ~.,
           a1_train,
           method = 'glmnet',
           tuneGrid = expand.grid(alpha = seq(0,1, length = 10),
                                  lambda = seq(0.0001, 10, length = 5)),
           trControl = custom)

la_en
en
set.seed(1000)
id2 = sample(2, nrow(a1), replace = T, prob = c(0.7, 0.3))
a1_train = a1[id2 == 1,]
a2_test = a1[id2 == 2,]

set.seed(1001)
id3 = sample(2, nrow(a2), replace = T, prob = c(0.7, 0.3))
a22_train = a2[id3 == 1,]
a22_test = a2[id3 == 2,]


la_lr = lm(dispframe~X1.12+X0.50+X19+X0.122+X29.5+X1664+X0.3+X2003+X0.91+X0.88+X1630, data = a1_train)
summary(la_lr)
plot(la_en)
plot(ridge.mod)
plot(lasso)
plot(lasso$finalModel, xvar = "lambda", label = T)
plot(lasso$finalModel, xvar = "dev", label = T)
?AIC

AIC(la_lr)
AIC(la_lm)
plot(ridge)
ridge
lasso
en
la_en
la_r = train(dispframe ~.,
              a1_train,
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = 0,
                                     lambda = seq(0.0001, 10, length = 5)),
              trControl = custom)
la_la = train(dispframe ~.,
             a1_train,
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 1,
                                    lambda = seq(0.0001, 10, length = 5)),
             trControl = custom)
#compare LA models
model_list1 = list(LA_Ridge = la_r, LA_Lasso = la_la, LA_EN = la_en)
res1 = resamples(model_list1)
summary(res1)

#RMSE for train
p11 = predict(la_en, a1_train)
r_train = sqrt(mean((a1_train$dispframe-p11)^2))
r_train
#RMSE for test
p22 = predict(la_en, a2_test)
r_test = sqrt(mean((a2_test$dispframe-p22)^2))
r_test
plot(la_r)
plot(la_la)
plot(la_r$finalModel, xvar = "lambda", label = T)
plot(lasso$finalModel, xvar = "dev", label = T)
plot(la_en)

la_en2 = train(dispframe ~.,
             a22_train,
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = seq(0,1, length = 10),
                                    lambda = seq(0.0001, 10, length = 5)),
             trControl = custom)
la_en2
#RMSE for train
p222 = predict(la_en2, a22_train)
r2_train = sqrt(mean((a22_train$dispframe-p222)^2))
r2_train
#RMSE for test
p221 = predict(la_en2, a22_test)
r2_test = sqrt(mean((a22_test$dispframe-p221)^2))
r2_test
plot(la_en)
line(la_en2)

y_s = data.frame(y[1:100, ])
tcp_s = data.frame(tcp[1:100,])


p222 = predict(la_en2, a22_train)
r2_train = sqrt(mean((a22_train$dispframe-p222)^2))
r2_train
#RMSE for test
p221 = predict(la_en2, a22_test)
r2_test = sqrt(mean((a22_test$dispframe-p221)^2))
r2_test

en3 = train(dispframe ~.,
              a1_train,
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = seq(0,1, length = 10),
                                     lambda = seq(0.0001, 200, length = 5)),
              trControl = custom)

en4 = train(dispframe ~.,
              a1_train[,-4],
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = seq(0,1, length = 10),
                                     lambda = seq(0.0001, 200, length = 5)),
              trControl = custom)



#compare LA models
model_list2 = list(la_la_x19 = la_la, la_la_alone = la_la1)
res2 = resamples(model_list2)
summary(res2)


#compare LA models
model_list3 = list(EN3 = en3, EN4 = en4)
res3 = resamples(model_list3)
summary(res3)
en3
en4
################################################################
##################### Linear Model ########################
lm_model = train(dispframe~.,
                 xy_train,
                 method = 'lm',
                 trControl = custom)
summary(lm_model)

###################### Predictions on xy_train ##############
xy_train$pred = predict(lm_model)

rmse_train <- RMSE(xy_train$pred, xy_train$dispframe)

rmse_train #25.42566

###################### Predictions on xy_test ##############
xy_test$pred = predict(lm_model, newdata = xy_test)

rmse_test <- RMSE(xy_test$pred, xy_test$dispframe)
rmse_test  #2727107416

###################### plots of predictions vs actuals
ggplot(xy_train, aes(x = pred, y = dispframe)) + 
  geom_point() +
  geom_abline(color = "blue")


ggplot(xy_test, aes(x = pred, y = dispframe)) + 
  geom_point() +
  geom_abline(color = "blue")


#train = sample(1:nrow(x), nrow(x)*0.70) or nrow(x)/2 if you want to halve it
#test = (-train)
#ytest = y[test]
#a = xy[train,]
#b = xy[test,]
#rm(train, test, a, b)

# get train and test
train = sample(1:nrow(x), nrow(x) * 0.70)
test = (-train)
x_train = x[train,]
x_test = x[test,]
y_train = y[train,]

# x and y for glmnet package
x_mat = model.matrix(dispframe~.-1, data = xy)
y_mat = xy$dispframe

#### Ridge & Lasso - ISLR method
#Ridge
ridge.mod = glmnet(x_mat, y_mat, alpha = 0)
plot(ridge.mod, xvar = "lambda", label = TRUE)
plot(ridge.mod, xvar = "dev", label = TRUE)
cv.ridge = cv.glmnet(x_mat, y_mat, alpha = 0)
plot(cv.ridge)
cv.ridge$lambda.min
#Lasso
lasso.mod = glmnet(x_mat, y_mat)
plot(lasso.mod, xvar = "lambda", label = TRUE)
plot(lasso.mod, xvar = "dev", label = TRUE)
cv.lasso = cv.glmnet(x_mat, y_mat)
plot(cv.lasso)
coef(cv.lasso)
cv.lasso$lambda.min
###use the train index / validation division to select lambda value
lasso.tr = glmnet(x_mat[train,], y_mat[train])
lasso.tr
str(lasso.tr)
lasso.pred = predict(lasso.tr, x_mat[-train])

alist = list(10,20,30,40,50)
multiply = function(x, factor){
  x * factor
}
x3 = lapply(alist, multiply, factor = 5)
x3
unlist(x3)

below_zero <- function(x) {
  return(x[x < 0])
}
bl = c(10, 2, 3, 4, -4, -100)
below_zero(bl)

sort(bl)
Sys.Date()
Sys.timezone()
Sys.time()
