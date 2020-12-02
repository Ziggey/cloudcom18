library(glmnet)
library(caret)
library(dplyr)
library(ggplot2)
library(leaps)
library(MASS)
#library(data.table)

# Experiments are based on the Rerngvit IM dataset on mldata.org

X = read.csv("XY.csv", header = TRUE)
#rm(X1, x0,x0_re)
Y = read.csv("Y.csv", header = TRUE)
names(Y)
X = X[,-1]
Y = Y[,-1]

#Extract the response variables in dataframe format
Y1 = as.data.frame(Y[,c(2,8,12)], drop=FALSE)
names(Y1) #noAudioPlayed, NoRTPPkts, DispFrames

### Combine x and y
XY = cbind(Y1,X)

#remove NAs

names(XY)
table(is.na(XY))
XY = na.omit(XY)

set.seed(102)
id = sample(2, nrow(XY), replace = T, prob = c(0.6, 0.4))
train = XY[id == 1,]
test = XY[id == 2,]
names(train)


get_best_result = function(caret_fit){
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}
#plot RTP versus TCPSCK(Load)
plot(Y$NoRTPPkts[1:25000], type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(X$tcpsck[1:25000], type = "l", col = "red")
legend("topright",legend=c("RTP", "Load"), col=c("blue", "red"),lty=1.2,cex=0.6)


####################### UNADJUSTED ##############################
####### Unadjusted learning using EN for Dispframes ##############

#custom control
en_cv = trainControl(method = "cv", number = 10)

en_D = train(DispFrames ~ .,
            train[,c(-1,-2)],
            method = 'glmnet',
            tuneGrid = expand.grid(alpha = 0:1,
                                   lambda = 10^seq(10, -2, length = 100)),
            trControl = en_cv)

en_D$bestTune
#  alpha    lambda
#113     1  0.01
#####################################################
############# Predict on the EN Model - DispFrames
#RMSE for train
pD1 = predict(en_D, train)
r_D_train_UA = sqrt(mean((pD1-train$DispFrames)^2))
r_D_train_UA  #2.987113

#RMSE for test
pD2 = predict(en_D, test)
r_D_test_UA = sqrt(mean((pD2-test$DispFrames)^2))
r_D_test_UA  #3.038319
###############################################

plot(train$dispframe, p1)
plot(test$dispframe, p2)

####################### UNADJUSTED ##############################
####### Unadjusted learning using EN for Audio Buffer Rate####
en_A = train(noAudioPlayed ~ .,
                 train[,c(-2,-3)],
                 method = 'glmnet',
                 tuneGrid = expand.grid(alpha = 0:1,
                                        lambda = 10^seq(10, -2, length = 100)),
                 trControl = en_cv)

en_A$bestTune
#alpha     lambda
# 1       0.01

#####################################################
############# Predict on the EN Model - Audio
#RMSE for train (Audio Buffer Rate)
pA1 = predict(en_A, train)
r_A_train_UA = sqrt(mean((pA1-train$noAudioPlayed)^2))
r_A_train_UA   #14.02757

#RMSE for test
pA2 = predict(en_A, test)
r_D_test_UA = sqrt(mean((pA2-test$noAudioPlayed)^2))
r_D_test_UA #14.34256
#######################################################
get_best_result(en_audio)
#alpha     lambda     RMSE   Rsquared      MAE    RMSESD RsquaredSD
# 1       0.01747528 19.71877 0.08758169 14.31381 0.3798021 0.01550241
#MAESD
#1 0.2163152

##############################################################
####### Unadjusted learning using EN for NoRTPPkts##
en_R = train(NoRTPPkts ~ .,
                     train[,c(-1,-3)],
                     method = 'glmnet',
                     tuneGrid = expand.grid(alpha = 0:1,
                                            lambda = 10^seq(10, -2, length = 100)),
                     trControl = en_cv)

en_R$bestTune
#alpha    lambda
#1         0.01

#####################################################
############# Predict on the EN Model - NoRTPPkts
#RMSE for train (NoRTPPkts)
pR1 = predict(en_R, train)
r_R_train_UA = sqrt(mean((pR1-train$NoRTPPkts)^2))
r_R_train_UA  #26.06071

#RMSE for test
pR2 = predict(en_R, test)
r_R_test_UA = sqrt(mean((pR2-test$NoRTPPkts)^2))
r_R_test_UA #25.96774
###############################################
#Compare the EN models for Dispframes, Audio & NoRTPPkts
model_list = list(Video = en_disp, Audio = en_audio, RTP = en_rtp)
res = resamples(model_list)
summary(res)

##################################################################################
################### Load Adjusted ##############
max(XY$tcpsck)
min(XY$tcpsck)

t30 = filter(XY, tcpsck == 30)
t40 = filter(XY, tcpsck == 40)
t50 = filter(XY, tcpsck == 50)
t60 = filter(XY, tcpsck == 60)
t70 = filter(XY, tcpsck == 70)

##########################-##############################################
#DispFrames
en30_D = train(DispFrames ~ .,
             t30[,c(-1,-2)],
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 0:1,
                                    lambda = 10^seq(10, -2, length = 100)),
             trControl = en_cv)
p30_D = predict(en30_D, t30)
r_p30_D = sqrt(mean((p30_D-t30$DispFrames)^2))
r_p30_D #3.133531

en40_D = train(DispFrames ~ .,
             t40[,c(-1,-2)],
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 0:1,
                                    lambda = 10^seq(10, -2, length = 100)),
             trControl = en_cv)

p40_D = predict(en40_D, t40)
r_p40_D = sqrt(mean((p40_D-t40$DispFrames)^2))
r_p40_D #3.223687

en50_D = train(DispFrames ~ .,
             t50[,c(-1,-2)],
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 0:1,
                                    lambda = 10^seq(10, -2, length = 100)),
             trControl = en_cv)
p50_D = predict(en50_D, t50)
r_p50_D = sqrt(mean((p50_D-t50$DispFrames)^2))[]
r_p50_D #2.621742

en60_D = train(DispFrames ~ .,
             t30[,c(-1,-2)],
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 0:1,
                                    lambda = 10^seq(10, -2, length = 100)),
             trControl = en_cv)
p60_D = predict(en60_D, t60)
r_p60_D = sqrt(mean((p60_D-t60$DispFrames)^2))
r_p60_D #5.184658

en70_D = train(DispFrames ~ .,
             t70[,c(-1,-2)],
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 0:1,
                                    lambda = 10^seq(10, -2, length = 100)),
             trControl = en_cv)
p70_D = predict(en70_D, t70)
r_p70_D = sqrt(mean((p70_D-t70$DispFrames)^2))
r_p70_D #1.997892

en70_audio$bestTune
###################################################################################

#audio
en30_A = train(noAudioPlayed ~ .,
             t30[,c(-2,-3)],
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 0:1,
                                    lambda = 10^seq(10, -2, length = 100)),
             trControl = en_cv)

p30_A = predict(en30_A, t30)
r_p30_A = sqrt(mean((p30_A-t30$noAudioPlayed)^2))
r_p30_A #10.22653

en40_A = train(noAudioPlayed ~ .,
                   t40[,c(-2,-3)],
                   method = 'glmnet',
                   tuneGrid = expand.grid(alpha = 0:1,
                                          lambda = 10^seq(10, -2, length = 100)),
                   trControl = en_cv)

p40_A = predict(en40_A, t40)
r_p40_A = sqrt(mean((p40_A-t40$noAudioPlayed)^2))
r_p40_A #8.155041

en50_A = train(noAudioPlayed ~ .,
                   t50[,c(-2,-3)],
                   method = 'glmnet',
                   tuneGrid = expand.grid(alpha = 0:1,
                                          lambda = 10^seq(10, -2, length = 100)),
                   trControl = en_cv)

p50_A = predict(en50_A, t50)
r_p50_A = sqrt(mean((p50_A-t50$noAudioPlayed)^2))
r_p50_A #9.903063 or 11.3505

en60_A = train(noAudioPlayed ~ .,
                   t60[,c(-2,-3)],
                   method = 'glmnet',
                   tuneGrid = expand.grid(alpha = 0:1,
                                          lambda = 10^seq(10, -2, length = 100)),
                   trControl = en_cv)

p60_A = predict(en60_A, t60)
r_p60_A = sqrt(mean((p60_A-t60$noAudioPlayed)^2))
r_p60_A #14.85639

en70_A = train(noAudioPlayed ~ .,
                   t70[,c(-2,-3)],
                   method = 'glmnet',
                   tuneGrid = expand.grid(alpha = 0:1,
                                          lambda = 10^seq(10, -2, length = 100)),
                   trControl = en_cv)

p70_A = predict(en70_A, t70)
r_p70_A = sqrt(mean((p70_A-t70$noAudioPlayed)^2))
r_p70_A #12.52547

#################################################################################
#RTP
en30_R = train(NoRTPPkts ~ .,
                   t30[,c(-1,-3)],
                   method = 'glmnet',
                   tuneGrid = expand.grid(alpha = 0:1,
                                          lambda = 10^seq(10, -2, length = 100)),
                   trControl = en_cv)

p30_R = predict(en30_R, t30)
r_p30_R = sqrt(mean((p30_R-t30$NoRTPPkts)^2))
r_p30_R #37.89082

en40_R = train(NoRTPPkts ~ .,
                 t40[,c(-1,-3)],
                 method = 'glmnet',
                 tuneGrid = expand.grid(alpha = 0:1,
                                        lambda = 10^seq(10, -2, length = 100)),
                 trControl = en_cv)

p40_R = predict(en40_R, t40)
r_p40_R = sqrt(mean((p40_R-t40$NoRTPPkts)^2))
r_p40_R #25.24521

en50_R = train(NoRTPPkts ~ .,
                 t50[,c(-1,-3)],
                 method = 'glmnet',
                 tuneGrid = expand.grid(alpha = 0:1,
                                        lambda = 10^seq(10, -2, length = 100)),
                 trControl = en_cv)

p50_R = predict(en50_R, t50)
r_p50_R = sqrt(mean((p50_R-t50$NoRTPPkts)^2))
r_p50_R #20.67602

en60_R = train(NoRTPPkts ~ .,
                 t60[,c(-1,-3)],
                 method = 'glmnet',
                 tuneGrid = expand.grid(alpha = 0:1,
                                        lambda = 10^seq(10, -2, length = 100)),
                 trControl = en_cv)

p60_R = predict(en60_R, t60)
r_p60_R = sqrt(mean((p60_R-t60$NoRTPPkts)^2))
r_p60_R #18.57304

en70_R = train(NoRTPPkts ~ .,
                 t70[,c(-1,-3)],
                 method = 'glmnet',
                 tuneGrid = expand.grid(alpha = 0:1,
                                        lambda = 10^seq(10, -2, length = 100)),
                 trControl = en_cv)

p70_R = predict(en70_R, t70)
r_p70_R = sqrt(mean((p70_R-t70$NoRTPPkts)^2))
r_p70_R #19.98439

get_best_result(en70_rtp)


###############################################################
########### Plot trials and calculations
###############################################################################
#Train and Test predictions

p1-train$DispFrames

plot(errors, type = "l", lty = 1.8, col = "blue")
lines(p2, type = "l", col = "red")
legend("topright",legend=c("Train RMSE", "Test RMSE"), col=c("blue", "red"),lty=1.2,cex=0.6)

train_error = train$DispFrames-p1
test_error = test$DispFrames-p2

par(mfrow = c(2,1))
par(mar=c(1,2,1,1))
plot(test_error[1:500], type = "l", lty = 1.8, col = "blue")
lines(train_error[1:500], type = "l", col = "red")
dev.off()

plot(train_error[1:500], type = "l", lty = 1.8, col = "blue")
lines(p30[1:500], type = "l", col = "red")
legend("topright",legend=c("Train RMSE", "Test RMSE"), col=c("blue", "red"),lty=1.2,cex=0.6)
################################################################################
get_best_result(en60_rtp)
plot(res)
a = c(40.11,32.73,24.34,24.61,23.46)
b = c(38.43,25.48,24.20,19.15,18.91)
plot(a, type = "l", col = "blue")
lines(b, type="l", col="red")

v = c(3.133531,3.223687,2.621742,5.184658,2.04679)
v_mean = mean(v)
v_mean - 3.038319

p2-test$DispFrames

par(mfrow = c(1,1))
vfr_ua = p2-test$DispFrames
vfr_30 = p30-t30$DispFrames
R1 = p40_R-t40$NoRTPPkts
r2 = p30_R-t30$NoRTPPkts
plot(e30[1:500], type = "l", lty = 1.8, col = "blue")
plot(e40[1:500], type = "l", lty = 1.8, col = "blue")
lines(e40[1:10], type = "l", col = "red")
lines(p40_R[1:10], type = "l", col = "black")
legend("topright",legend=c("UA Error", "LA_30 Error"), col=c("blue", "red"),lty=1.2,cex=0.6)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(Y$NoRTPPkts[1:25000], type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(X$tcpsck[1:25000], type = "l", col = "red")
legend("topright",legend=c("RTP", "Load"), col=c("blue", "red"),lty=1.2,cex=0.6)

plot(Y$DispFrames[1:100], type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(X$tcpsck[1:100], type = "l", col = "red")
legend("topright",legend=c("RTP", "Load"), col=c("blue", "red"),lty=1.2,cex=0.6)

############################
e30 = p30_R-t30$NoRTPPkts
e30a = p30_A-t30$noAudioPlayed
e30r = p30_R-t30$NoRTPPkts
la30 = cbind(e30,e30a,e30r)

aa = data.frame(la30)
names(aa) = c("RTP", "ABR", "VFR")
boxplot(aa, ylab="MSE")

#######################
# plot the RTP errors against the load values
e30 = p30_R-t30$NoRTPPkts
e40 = p40_R-t40$NoRTPPkts
e50 = p50_R-t50$NoRTPPkts
e60 = p60_R-t60$NoRTPPkts
e70 = p70_R-t70$NoRTPPkts
rtp_all = data.frame(e30[1:600],e40[1:600],e50[1:600],e60[1:600],e70[1:600])
names(rtp_all)= c("30", "40", "50", "60", "70")
par(mfrow = c(1,1))
boxplot(rtp_all, main="Errors in RTP for each Load Value",ylab="Errors in RTP rate estimates", xlab="Load Values, K")

par(mfrow = c(1,2))
boxplot(e30)
boxplot(e40)
################################

par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
plot(1:10)
mtext("Small Margins", side = 3, line = 1, cex = 1.2)
mtext("par(mar = c(2, 2, 2, 2))", side = 3)

# Second Plot with large margins
par(mar = c(5, 5, 5, 5)) # Set the margin on all sides to 6
plot(1:10)
mtext("Large Margins", side = 3, line = 1, cex = 1.2)
mtext("par(mar = c(5, 5, 5, 5))", side = 3)

par(mfrow = c(2,1))
plot(test$NoRTPPkts[1:1000], type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(p6[1:1000], type = "l", col = "red")
legend("topleft",legend=c("Actual RTP", "LAEN estimates for RTP"), col=c("blue", "red"),lty=1.2,pt.cex=0.6)

plot(test$NoRTPPkts[1:1000], type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(p30_R[1:1000], type = "l", col = "red")
legend("topleft",legend=c("Actual RTP", "LAEN estimates for RTP"), col=c("blue", "red"),lty=1.2,cex=0.2)
###############
dev.off()
e30a = p30_A-t30$noAudioPlayed
e40a = p40_A-t40$noAudioPlayed

e30d = p30-t30$DispFrames
e40d = p40-t40$DispFrames
a_d = data.frame(e30a[1:700],e40a[1:700],e30d[1:700],e40d[1:700])
names(a_d) = c("30a", "40a", "30d", "40d")
min(a_d[,1])
###################################################
#Save entire workspace

save.image(file = "paper_work_ext.RData")

####### 

get_best_result = function(caret_fit){
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}
squares = function(a,b){
  a^2 + b^2
  
}
squares(2,4)

