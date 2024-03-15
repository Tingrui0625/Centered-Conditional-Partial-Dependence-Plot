#loading the dataset in R
set.seed(7832)
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")
data("kc_housing", package = "mlr3data")

# data cleaning -----------------------------------------------------------

#transforming the feature data to be numeric
dates = anytime(kc_housing$date)
kc_housing$date = as.numeric(difftime(dates, min(dates), units = "days"))

kc_housing$zipcode = as.factor(kc_housing$zipcode)
kc_housing$renovated = as.numeric(!is.na(kc_housing$yr_renovated))
kc_housing$has_basement = as.numeric(!is.na(kc_housing$sqft_basement))
kc_housing$id = NULL
kc_housing$price = kc_housing$price / 1000
# missing values
kc_housing$yr_renovated = NULL
kc_housing$sqft_basement = NULL
kc_housing[,8]=ifelse(kc_housing[,8]==TRUE,1,0)

#printing the variable types
for(i in 1:ncol(kc_housing)){
  print(paste0(i,"th:",colnames(kc_housing)[i]," is ",class(kc_housing[,i])))
}

#transforming the categorical variables to factors
for(i in c(3,8,9,10,11,13)){
  kc_housing[,i] = as.factor(kc_housing[,i])
}

#transforming the integers to numeric values
for(i in c(5,6,12,17,18)){
  kc_housing[,i] = as.numeric(kc_housing[,i])
}
#kc_housing$id=1:nrow(kc_housing)
kc_housing1=kc_housing

#deleting outlier rows
for(i in c(5,6,12,17,18)){
    quartiles <- quantile(kc_housing1[,i], probs=c(.25, .75), na.rm = FALSE)
    IQR <- IQR(kc_housing1[,i])
    
    Lower <- quartiles[1] - 1.5*IQR
    Upper <- quartiles[2] + 1.5*IQR
    if(length(which(kc_housing1[,i] > Lower & kc_housing1[,i] < Upper))!=0){
      kc_housing1_no_outlier= subset(kc_housing1, kc_housing1[,i] > Lower & kc_housing1[,i] < Upper)
    }
    #print(which())
    kc_housing1=kc_housing1_no_outlier
    }
kc_housing1=na.omit(kc_housing1)

#normalization
for(i in 1:ncol(kc_housing1)){
  if(class(kc_housing1[,i])=="numeric"){
    print(i)
    kc_housing1[1:nrow(kc_housing1),i]=(kc_housing1[,i]-mean(kc_housing1[,i]))/sd(kc_housing1[,i])
  }
}
kc_housing1[,"price"]=kc_housing1_no_outlier[,"price"]

#saving the clean data
save(kc_housing1,file="Case2-kc_housing/kc_housing_data.RData")

# loading clean data ------------------------------------------------------
load(file="Case2-kc_housing/kc_housing_data.RData")
# tuning ------------------------------------------------------------------
measure = msrs("regr.mse")
rsmp_cv3 = rsmp("cv", folds = 3)
msr_mse = msr("regr.mse")
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")
future::plan("multisession", workers = 8)

tnr_random_search = tnr("random_search", batch_size = 8)

lrn_effect = lrn("regr.nnet",
                 decay  = to_tune(1e-5, 0.1),#logscale = TRUE
                 size = to_tune(1, 10),
                 trace=FALSE,
                 MaxNWts=84581
)

s=Sys.time()
instance = tune(
  tuner = tnr_random_search,
  task = task,
  learner = lrn_effect,
  resampling = rsmp_cv3,
  measures = msr_mse,
  term_evals = 800
)
instance$result$x_domain
e=Sys.time()
e-s
#needs 2 hours

test_id=sample(nrow(kc_housing1),0.1*nrow(kc_housing1))
kc_housing_test=kc_housing1[test_id,]
kc_housing_train=kc_housing1[-test_id,]

task = as_task_regr(kc_housing_train,target="price")
#The hyperparameters need to be tuned when the test_id is changed
size= 10
decay=0.06229317
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE,MaxNWts=84581)#,size=size,decay=decay,trace=FALSE
model = lrn$train(task = task)

#model evaluation
task_test=as_task_regr(kc_housing_test,target="price")
prediction=model$predict(task_test)
measure = msrs("regr.mse")#msrs("regr.rmse")
prediction$score(measure)

#baseline model
prediction_base = lrn("regr.featureless")$train(task)$predict(task_test)
prediction_base$score(measure)


# Effects and Plotting ----------------------------------------------------------------
h=50 
feature=colnames(kc_housing_train)[18]

#sample plotting observations
plot_id=sample(nrow(kc_housing_train),0.1*nrow(kc_housing_train))
kc_housing_plot=kc_housing_train[plot_id,]

xswcpdp <-effect2 (mod = model$model, data=kc_housing_plot, feature = feature, target="price",
                   predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 0.46,kernel.width=1.3)

xcwcpdp = effect2(mod = model$model, data=kc_housing_plot, feature = feature, target="price",
                  gower.power = 0.46, predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=1.3)

xcwccpdp = effect2(mod = model$model, data=kc_housing_plot, feature = feature, target="price",
                   gower.power = 0.46, predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=1.3)

xswccpdp <- effect2 (mod = model$model, data=kc_housing_plot, feature = feature, target="price",
                     predict.fun = predict, h = h, method = "xs-wccpdp",gower.power = 0.46,kernel.width=1.3)

cpdp = effect2(mod = model$model, data=kc_housing_plot, feature = feature, target="price",
               gower.power = 0.46, predict.fun = predict, h = h, method = "cpdp",kernel.width=1.7)
pred <- Predictor$new(model=model$model, data=kc_housing_plot, y = "price")
q<-quantile(kc_housing_plot[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)

#mpl = mplot(kc_housing1[plot_id,], feature, target="price", eps = 0.5)

xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
xswccpdp$y=xswccpdp$y-mean(xswccpdp$y,na.rm=TRUE)
xcwcpdp$y=xcwcpdp$y-mean(xcwcpdp$y)
xcwccpdp$y=xcwccpdp$y-mean(xcwccpdp$y,na.rm=TRUE)
cpdp$y=cpdp$y-mean(cpdp$y)

pdp$results$.value=pdp$results$.value-mean(pdp$results$.value)
p1=ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  #geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 2) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 2, lwd = 2) +
  geom_line(data = xswcpdp, aes(x = x, y = y, col = "xswcpdp"), lty = 1, lwd = 1)+
  #geom_line(data = xswccpdp, aes(x = x, y = y, col = "xswccpdp"), lty = 4, lwd = 2) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2)+ 
  #geom_line(data = xcwcpdp, aes(x = x, y = y, col = "xcwcpdp"), lty = 5, lwd = 2) +
  #geom_line(data = xcwccpdp, aes(x = x, y = y, col = "xcwccpdp"), lty = 6, lwd = 2) +
  #geom_line(data = mpl, aes(x, mplot, col = "M-plot"),lty = 8, lwd = 2) +
  # geom_line(data = xsxcpdp, aes(x = x, y = y, col = "xsxcpdp"), lty = 3, lwd = 2) +
  # geom_line(data = xsxccpdp, aes(x = x, y = y, col = "xsxccpdp"), lty = 4, lwd = 2) +
  # geom_line(data = xsxcpdp1, aes(x = x, y = y, col = "xsxcpdp1"), lty = 3, lwd = 2) +
  # geom_line(data = xsxccpdp1, aes(x = x, y = y, col = "xsxccpdp1"), lty = 4, lwd = 2) +
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
  labs(x=feature)

ggsave(p1,
       file=paste0("effect2_case2_",feature,"_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )
# effect1 -----------------------------------------------------------------


cpdp <- effect1 (mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price", 
                 predict.fun = predict, h = h, method = "cpdp",gower.power = 2,kernel.width=2)

xswcpdp <-effect1 (mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                   predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 2,kernel.width=2)

xcwcpdp = effect1(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                  gower.power = 2, predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=2)

xcwccpdp = effect1(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                   gower.power = 2, predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=2)

xswccpdp <- effect1 (mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                     predict.fun = predict, h = h, method = "xs-wccpdp",gower.power = 2,kernel.width=2)
xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
xswccpdp$y=xswccpdp$y-mean(xswccpdp$y)
xcwcpdp$y=xcwcpdp$y-mean(xcwcpdp$y)
xcwccpdp$y=xcwccpdp$y-mean(xcwccpdp$y)
cpdp$y=cpdp$y-mean(cpdp$y)

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
       file=paste0("effect1_case2_",feature,"_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )


# effect4 -----------------------------------------------------------------

xsxcpdp = effect4(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                  predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=14.7)
#gamma=13
xsxccpdp = effect4(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                   predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma=14.7)

cpdp = effect4(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
               predict.fun = predict, h = h, method = "cpdp",gamma=14.7)

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
       file=paste0("effect4_case2_",feature,"_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )

# effect3 -----------------------------------------------------------------


xsxcpdp = effect3(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                  predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=0.05)
#gamma=0.13
xsxccpdp = effect3(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
                   predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma=0.05)

cpdp = effect3(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",
               predict.fun = predict, h = h, method = "cpdp",gamma=0.05)

xsxcpdp$y=xsxcpdp$y-mean(xsxcpdp$y)
xsxccpdp$y=xsxccpdp$y-mean(xsxccpdp$y)

p1=ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  #geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
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
       file=paste0("effect3_case2_",feature,"_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )


# effect5 -----------------------------------------------------------------
xsxcpdp = effect5(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",kernel.width=1.3,gower.power=0.46,
                  predict.fun = predict, h = h, method = "xs-xc-pdp")

xsxccpdp = effect5(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",kernel.width=1.3,gower.power=0.46,
                   predict.fun = predict, h = h, method = "xs-xc-cpdp")

cpdp = effect5(mod = model$model, data=kc_housing1[plot_id,], feature = feature, target = "price",kernel.width=1.3,gower.power=0.46,
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
       file=paste0("effect5_case2_",feature,"_all.png"),
       scale = 2.5,
       width = 7.75,             
       height = 5.88,
       units="cm",
       dpi = 500 )

Loss1=function(a,b){
  a=a$y[!is.na(curve$y)]
  a=a-mean(a,na.rm=TRUE)
  b1=b$results$.value[1:length(a)]
#b=pdp$results$.value[1:length(a)]
#b=b-mean(b)
  return(sum((a-b1)^2)/length(a))
}
Loss1(xsxcpdp,ale)
#effect5-> 1.125876

Loss1(xswcpdp,ale)
#effect2-> 1.125876

# ggplot(xswcpdp,aes(x = x, y = y, col = "xswcpdp"))+
#   geom_point() +
#   geom_rug(col="steelblue",alpha=0.1, linewidth=1.5)+
#   geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 2) +
#   geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 7, lwd = 2)
#   #scale_x_continuous(limits=c(0,2))



