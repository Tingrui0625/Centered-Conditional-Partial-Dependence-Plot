#loading the Wine quality dataset
data <- read.csv("Case1-Red-wine-quality-master/winequality-red.csv", sep = ';')
data1=data

#removing outliers
for(i in 1:(ncol(data)-1)){
  quartiles <- quantile(data1[,i], probs=c(.25, .75), na.rm = FALSE)
  IQR <- IQR(data1[,i])

  Lower <- quartiles[1] - 1.5*IQR
  Upper <- quartiles[2] + 1.5*IQR
  if(length(which(data1[,i] > Lower & data1[,i] < Upper))!=0){
    data1_no_outlier= subset(data1, data1[,i] > Lower & data1[,i] < Upper)
  }
  data1=data1_no_outlier
}

#tuning hyperparameters with file "1. Tuning and Model Selection.R"
size=11
decay=0.0911569

#saparating the data into train&test parts
set.seed(41)
plot_id=sample(nrow(data1),0.1*nrow(data1))
data_test=data1[plot_id,]
data_train=data1[-plot_id,]

#train the model
task = as_task_regr(data_train,target="quality")
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
model = lrn$train(task = task)

#evaluate the model with test dataset
task_test=as_task_regr(data_test,target="quality")
prediction=model$predict(task_test)
measure=msr("regr.mse")
prediction$score(measure)
# regr.mse 
#0.3278117

#baseline model
prediction_base = lrn("regr.featureless")$train(task)$predict(task_test)
prediction_base$score(measure)
# regr.mse 
# 0.5022465 

#important features
feature=colnames(data)[6]
feature=colnames(data)[7]

#hyperparameters for ccPDPs are tuned with file "3. Tuning the hyperparameters in ccPDPs.R"
xswcpdp <-effect2 (mod = model$model, data = data1, feature = feature, target = "quality",
                   predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 1,kernel.width=2)

xcwcpdp = effect2(mod = model$model, data = data1, feature = feature, target = "quality",
                  gower.power = 0.00133, predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=2)

xcwccpdp = effect2(mod = model$model, data = data1, feature = feature, target = "quality",
                   gower.power = 0.00133, predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=2)

xswccpdp <- effect2 (mod = model$model, data = data1, feature = feature, target = "quality",
                    predict.fun = predict, h = h, method = "xs-wccpdp",gower.power = 0.00133,kernel.width=2)

cpdp = effect2(mod = model$model, data = data1, feature = feature, target = "quality",
               gower.power = 0.00133, predict.fun = predict, h = h, method = "cpdp",kernel.width=2)
pred <- Predictor$new(model=model$model, data = data1, y = "quality")
q<-quantile(data1[,feature], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
#ice <- FeatureEffect$new(pred, feature = feature, method = "ice", grid.points =q)

mpl = mplot(data1, feature, target = "quality", eps = 0.5)

xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
xswccpdp$y=xswccpdp$y-mean(xswccpdp$y,na.rm=TRUE)
xcwcpdp$y=xcwcpdp$y-mean(xcwcpdp$y)
xcwccpdp$y=xcwccpdp$y-mean(xcwccpdp$y,na.rm=TRUE)
cpdp$y=cpdp$y-mean(cpdp$y)

pdp$results$.value=pdp$results$.value-mean(pdp$results$.value)
mpl$mplot=mpl$mplot-mean(mpl$mplot)
p1=ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = xswcpdp, aes(x = x, y = y, col = "xswcpdp"), lty = 1, lwd = 1) +
  geom_line(data = xswccpdp, aes(x = x, y = y, col = "xswccpdp"), lty = 1, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  geom_line(data = xcwcpdp, aes(x = x, y = y, col = "xcwcpdp"), lty = 1, lwd = 1) +
  geom_line(data = xcwccpdp, aes(x = x, y = y, col = "xcwccpdp"), lty = 1, lwd = 1) +
  #geom_line(data = mpl, aes(x, mplot, col = "mplot"),lty = 1, lwd = 1) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  scale_color_manual(name = "Effects", 
                     values = c("ale" = "lightgreen", 
                                "pdp" = "orange1", 
                                "cpdp" = "grey27",
                                "xswcpdp" = "maroon3",
                                "xswccpdp"="palevioletred",
                                "xcwcpdp"="dodgerblue1",
                                "xcwccpdp"="cornflowerblue",
                                "xsxcpdp"="blueviolet",
                                "cxsxcpdp"="mediumpurple4",
                                "mplot"="mistyrose3"
                     ))+
  guides()
ggsave(p1,
       file=paste0("effect2_case1_x7_All.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )
ggplot(xswcpdp,aes(x = x, y = y, col = "xswcpdp"))+
  geom_point() +
  geom_rug(col="steelblue",alpha=0.1, linewidth=1.5)+
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 2) +
  geom_line(data = xswcpdp, aes(x = x, y = y, col = "xswcpdp"), lty = 1, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 7, lwd = 2)
# effect1 -----------------------------------------------------------------


cpdp <- effect1 (mod = model$model, data = data1, feature = feature, target = "quality", 
                 predict.fun = predict, h = h, method = "cpdp",gower.power = 0.00258,kernel.width=2)
#kernel.width=1.22
#fixed_sigma=0.5
#gamma=0.1
xswcpdp <-effect1 (mod = model$model, data = data1, feature = feature, target = "quality",
                   predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 0.00258,kernel.width=2)

xcwcpdp = effect1(mod = model$model, data = data1, feature = feature, target = "quality",
                  gower.power = 0.00258, predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=2)

xcwccpdp = effect1(mod = model$model, data = data1, feature = feature, target = "quality",
                   gower.power = 0.00258, predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=2)

xswccpdp <- effect1 (mod = model$model, data = data1, feature = feature, target = "quality",
                     predict.fun = predict, h = h, method = "xs-wccpdp",gower.power = 0.00258,kernel.width=2)

xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
xswccpdp$y=xswccpdp$y-mean(xswccpdp$y)
xcwcpdp$y=xcwcpdp$y-mean(xcwcpdp$y)
xcwccpdp$y=xcwccpdp$y-mean(xcwccpdp$y)
cpdp$y=cpdp$y-mean(cpdp$y)
pdp$results$.value=pdp$results$.value-mean(pdp$results$.value)
p1=ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = xswcpdp, aes(x = x, y = y, col = "xswcpdp"), lty = 1, lwd = 1) +
  #geom_line(data = xswccpdp, aes(x = x, y = y, col = "xswccpdp"), lty = 1, lwd = 1) +
  #geom_line(data = xcwcpdp, aes(x = x, y = y, col = "xcwcpdp"), lty = 1, lwd = 1) +
  #geom_line(data = xcwccpdp, aes(x = x, y = y, col = "xcwccpdp"), lty = 1, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 8, lwd = 2) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  scale_color_manual(name = "Effects", 
                     values = c("ale" = "lightgreen", 
                                "pdp" = "orange1", 
                                "cpdp" = "grey27",
                                "xswcpdp" = "maroon3",
                                "xswccpdp"="palevioletred",
                                "xcwcpdp"="dodgerblue1",
                                "xcwccpdp"="cornflowerblue",
                                "xsxcpdp"="blueviolet",
                                "cxsxcpdp"="mediumpurple4",
                                "mplot"="mistyrose3"
                     ))+
  guides()
ggsave(p1,
       file=paste0("effect1_case1_x7_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )


# effect4 -----------------------------------------------------------------

xsxcpdp = effect4(mod = model$model, data = data, feature = feature, target = "quality",
                  predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=20)

xsxccpdp = effect4(mod = model$model, data = data, feature = feature, target = "quality",
                  predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma=20)

cpdp = effect4(mod = model$model, data = data, feature = feature, target = "quality",
                predict.fun = predict, h = h, method = "cpdp",gamma=20)

xsxcpdp$y=xsxcpdp$y-mean(xsxcpdp$y)
xsxccpdp$y=xsxccpdp$y-mean(xsxccpdp$y)

p1=ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = xsxcpdp, aes(x = x, y = y, col = "xsxcpdp"), lty = 1, lwd = 1) +
  #geom_line(data = xsxccpdp, aes(x = x, y = y, col = "cxsxcpdp"), lty = 1, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) + 
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 8, lwd = 2) +
  scale_color_manual(name = "Effects", 
                     values = c("ale" = "lightgreen", 
                                "pdp" = "orange1", 
                                "cpdp" = "grey27",
                                "xswcpdp" = "maroon3",
                                "xswccpdp"="palevioletred",
                                "xcwcpdp"="dodgerblue1",
                                "xcwccpdp"="cornflowerblue",
                                "xsxcpdp"="blueviolet",
                                "cxsxcpdp"="mediumpurple4",
                                "mplot"="mistyrose3"
                     ))+
  guides()
ggsave(p1,
       file=paste0("effect4_case1_x7_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )

# effect3 -----------------------------------------------------------------


xsxcpdp = effect3(mod = model$model, data = data, feature = feature, target = "quality",
                   predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=20)

xsxccpdp = effect3(mod = model$model, data = data, feature = feature, target = "quality",
                   predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma=20)

cpdp = effect3(mod = model$model, data = data, feature = feature, target = "quality",
               predict.fun = predict, h = h, method = "cpdp",gamma=20)

xsxcpdp$y=xsxcpdp$y-mean(xsxcpdp$y)
xsxccpdp$y=xsxccpdp$y-mean(xsxccpdp$y)

p1=ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  #geom_line(data = xsxccpdp, aes(x = x, y = y, col = "cxsxcpdp"), lty = 1, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  geom_line(data = xsxcpdp, aes(x = x, y = y, col = "xsxcpdp"), lty = 1, lwd = 1) +
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 1, lwd = 1) +
  scale_color_manual(name = "Effects", 
                     values = c("ale" = "lightgreen", 
                                "pdp" = "orange1", 
                                "cpdp" = "grey27",
                                "xswcpdp" = "maroon3",
                                "xswccpdp"="palevioletred",
                                "xcwcpdp"="dodgerblue1",
                                "xcwccpdp"="cornflowerblue",
                                "xsxcpdp"="blueviolet",
                                "cxsxcpdp"="mediumpurple4",
                                "mplot"="mistyrose3"
                     ))+
  guides()
ggsave(p1,
       file=paste0("effect3_case1_x7_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )


# effect5 -----------------------------------------------------------------
xsxcpdp = effect5(mod = model$model, data = data, feature = feature, target = "quality",kernel.width=2,gower.power=0.00133,
                  predict.fun = predict, h = h, method = "xs-xc-pdp")

xsxccpdp = effect5(mod = model$model, data = data, feature = feature, target = "quality",kernel.width=2,gower.power=0.00133,
                   predict.fun = predict, h = h, method = "xs-xc-cpdp")

cpdp = effect5(mod = model$model, data = data, feature = feature, target = "quality",kernel.width=2,gower.power=0.00133,
               predict.fun = predict, h = h, method = "cpdp")

xsxcpdp$y=xsxcpdp$y-mean(xsxcpdp$y)
xsxccpdp$y=xsxccpdp$y-mean(xsxccpdp$y)

p1=ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  #geom_line(data = xsxccpdp, aes(x = x, y = y, col = "xsxccpdp"), lty = 1, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  geom_line(data = xsxcpdp, aes(x = x, y = y, col = "xsxcpdp"), lty = 1, lwd = 1) +
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 1, lwd = 1) +
  scale_color_manual(name = "Effects", 
                     values = c("ale" = "lightgreen", 
                                "pdp" = "orange1", 
                                "cpdp" = "grey27",
                                "xswcpdp" = "maroon3",
                                "xswccpdp"="palevioletred",
                                "xcwcpdp"="dodgerblue1",
                                "xcwccpdp"="cornflowerblue",
                                "xsxcpdp"="blueviolet",
                                "cxsxcpdp"="mediumpurple4",
                                "mplot"="mistyrose3"
                     ))+
  guides()
ggsave(p1,
       file=paste0("effect5_case1_x7_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )


