# general setting ---------------------------------------------------------

setwd(" ")
core_number = 8
example_name = "example4"
feature = c("x1","x2")
effect_name = "effect1"
set_seed_number = 41
iteration_number = 112

feature="x4"
h = 50
set.seed(41)
data = create_xor_corr(n = 5000)
X = data[,setdiff(names(data),"y")]
task = as_task_regr(data,target="y")
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
model = lrn$train(task = task)

cpdp <- effect1 (mod = model$model, data = data, feature = feature, target = "y", 
                 predict.fun = predict, h = h, method = "cpdp",gower.power = 2,kernel.width=2)

xswcpdp <-effect1 (mod = model$model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 2,kernel.width=2)

xcwcpdp = effect1(mod = model$model, data = data, feature = feature, target = "y",
                  gower.power = 2, predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=2)

xcwccpdp = effect1(mod = model$model, data = data, feature = feature, target = "y",
                   gower.power = 2, predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=2)

xswccpdp <- effect1 (mod = model$model, data = data, feature = feature, target = "y",
                     predict.fun = predict, h = h, method = "xs-wccpdp",gower.power = 2,kernel.width=2)

pred <- Predictor$new(model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)

mpl = mplot(data, feature, target = "y", eps = 0.5)

xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
xswccpdp$y=xswccpdp$y-mean(xswccpdp$y)
xcwcpdp$y=xcwcpdp$y-mean(xcwcpdp$y)
xcwccpdp$y=xcwccpdp$y-mean(xcwccpdp$y)
cpdp$y=cpdp$y-mean(cpdp$y)

pdp$results$.value=pdp$results$.value-mean(pdp$results$.value)
ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 2) +
  #geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 2, lwd = 2) +
  geom_line(data = xswcpdp, aes(x = x, y = y, col = "xswcpdp"), lty = 3, lwd = 2) +
  geom_line(data = xcwcpdp, aes(x = x, y = y, col = "xcwcpdp"), lty = 3, lwd = 2) +
  #geom_line(data = xswccpdp, aes(x = x, y = y, col = "xswccpdp"), lty = 4, lwd = 2) +
  #geom_line(data = xcwcpdp, aes(x = x, y = y, col = "xcwcpdp"), lty = 5, lwd = 2) +
  #geom_line(data = xcwccpdp, aes(x = x, y = y, col = "xcwccpdp"), lty = 6, lwd = 2) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 7, lwd = 2) +
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 8, lwd = 2) +
  guides()

# test for effect 2 -------------------------------------------------------

feature="x2"
h = 20
set.seed(41)
data = create_xor_corr(n = 1000)
X = data[,setdiff(names(data),"y")]
task = as_task_regr(data,target="y")
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
model = lrn$train(task = task)
# fit_xs <- kde1d(data[[feature]])
# fit_xc <- kdevine(data[, setdiff(colnames(X),feature), drop = FALSE])

xswcpdp <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 2,kernel.width=0.6)

# xcwcpdp = effect2(mod = model$model, data = data, feature = feature, target = "y",
#                   gower.power = 0.49, predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=0.1341259)
# 
# xcwccpdp = effect2(mod = model$model, data = data, feature = feature, target = "y",
#                    gower.power = 0.49, predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=0.1341259)
# 
# xswccpdp <- effect2 (mod = model$model, data = data, feature = feature, target = "y",
#                      predict.fun = predict, h = h, method = "xs-wccpdp",gower.power = 0.49,kernel.width=0.1341259)
# 
cpdp = effect2(mod = model$model, data = data, feature = feature, target = "y",
               gower.power = 2, predict.fun = predict, h = h, method = "cpdp",kernel.width=0.6)
pred <- Predictor$new(model=model$model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp+ice", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)

mpl = mplot(data, feature, target = "y", eps = 0.5)

xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
# xswccpdp$y=xswccpdp$y-mean(xswccpdp$y,na.rm=TRUE)
# xcwcpdp$y=xcwcpdp$y-mean(xcwcpdp$y)
# xcwccpdp$y=xcwccpdp$y-mean(xcwccpdp$y,na.rm=TRUE)
cpdp$y=cpdp$y-mean(cpdp$y)
pdp$results$.value=pdp$results$.value-mean(pdp$results$.value)


# library(yaImpute)
# install.packages("yaImpute")
# eff <- FeatureEffect$new(pred, feature = c("x2", "x4"),method="ale")
# eff$plot()
# eff$plot(show.data = TRUE)

ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = xswcpdp, aes(x = x, y = y, col = "xswcpdp"), lty = 3, lwd = 1) +
  #geom_line(data = xswccpdp, aes(x = x, y = y, col = "xswccpdp"), lty = 4, lwd = 2) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  #geom_line(data = xcwcpdp, aes(x = x, y = y, col = "xcwcpdp"), lty = 5, lwd = 2) +
  #geom_line(data = xcwccpdp, aes(x = x, y = y, col = "xcwccpdp"), lty = 6, lwd = 2) +
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 8, lwd = 2) +
  #geom_line(data = curve, aes(x = x, y = y, col = "knockoff"),lty = 8, lwd = 2) +
  # geom_line(data = xsxcpdp, aes(x = x, y = y, col = "xsxcpdp"), lty = 3, lwd = 2) +
  # geom_line(data = xsxccpdp, aes(x = x, y = y, col = "xsxccpdp"), lty = 4, lwd = 2) +
  # geom_line(data = xsxcpdp1, aes(x = x, y = y, col = "xsxcpdp1"), lty = 3, lwd = 2) +
  # geom_line(data = xsxccpdp1, aes(x = x, y = y, col = "xsxccpdp1"), lty = 4, lwd = 2) +
  guides()

# test for effect 3 -------------------------------------------------------


xsxcpdp = effect3(mod = model$model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=0.11)

xsxccpdp = effect3(mod = model$model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma=0.11)

cpdp = effect3(mod = model$model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "cpdp",gamma=0.11)
pred <- Predictor$new(model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)

mpl = mplot(data, feature, target = "y", eps = 0.5)

xsxcpdp$y=xsxcpdp$y-mean(xsxcpdp$y)
xsxccpdp$y=xsxccpdp$y-mean(xsxccpdp$y)

ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = xsxccpdp, aes(x = x, y = y, col = "xsxccpdp"), lty = 1, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  geom_line(data = xsxcpdp, aes(x = x, y = y, col = "xsxcpdp"), lty = 1, lwd = 1) +
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 1, lwd = 1) +
  guides()
#par(new=TRUE)
#scatter.smooth(x=xsxcpdp$x,y=xsxcpdp$y,lpars = list(col = "blue", lwd = 1, lty = 1))
# test for effect 4 -------------------------------------------------------


xsxcpdp = effect4(mod = model$model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=20)

xsxccpdp = effect4(mod = model$model, data = data, feature = feature, target = "y",
                   gower.power = 1, predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma=25)

cpdp = effect4(mod = model$model, data = data, feature = feature, target = "y",
               gower.power = 1, predict.fun = predict, h = h, method = "cpdp",gamma=25)
pred <- Predictor$new(model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)

mpl = mplot(data, feature, target = "y", eps = 0.5)
xsxcpdp$y=xsxcpdp$y-mean(xsxcpdp$y)
xsxccpdp$y=xsxccpdp$y-mean(xsxccpdp$y)

ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 2) +
  #geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 2, lwd = 2) +
  geom_line(data = xsxcpdp, aes(x = x, y = y, col = "xsxcpdp"), lty = 3, lwd = 2) +
  #geom_line(data = xsxccpdp, aes(x = x, y = y, col = "xsxccpdp"), lty = 4, lwd = 2) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 7, lwd = 2) + 
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 8, lwd = 2) +
  guides()







