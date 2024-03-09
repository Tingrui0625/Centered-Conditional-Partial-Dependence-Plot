
# functions for repeatedly creating data and computing ccpdps -------------

do_once_effect1<-function(feature,decay,size,pars){
  h = 50
  result=matrix(0,ncol=9,nrow=h+1)
  data = create_xor_corr(n = 5000)
  X = data[,setdiff(names(data),"y")]
  
  #switching model
  if(model_name=="nnet"){
    task = as_task_regr(data,target="y")
    lrn = lrn("regr.nnet",size=size,decay=decay)
    model = lrn$train(task = task)$model
  }else{
    #selecting model according to example
    model=switch(example_name,
                 "example1"=gam(y ~ x1+I(x2^2)+I(x3^3)+I(x2*x4), data = data),
                 "example2"=gam(y ~ I(x1^2)+x2,data = data),
                 "example3"=gam(y ~ x1+x2+I(x1*x2),data=data),
                 "example4"=gam(y ~ x1+I(x2^2)+I(x3^3)+x3+I(x2*x4),data=data),
                 "example5"=gam(y ~ I(x1^2)+x2,data=data),
                 "example6"=gam(y ~ I(x1^2)+x2,data=data)
    )
  }
  #reading hyperparameters from csv files
  kernel.width=switch(feature,
                      "x1"=pars[1,1],
                      "x2"=pars[2,1],
                      "x3"=pars[3,1],
                      "x4"=pars[4,1],
                      "x5"=pars[5,1])
  gower.power=switch(feature,
                     "x1"=pars[1,2],
                     "x2"=pars[2,2],
                     "x3"=pars[3,2],
                     "x4"=pars[4,2],
                     "x5"=pars[5,2])
  
  #computing ccpdps with effect one given parameters and created dataset
  cpdp <- effect1 (mod = model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "cpdp",kernel.width=kernel.width,gower.power=gower.power)

  xswcpdp <-effect1 (mod = model, data = data, feature = feature, target = "y",
                     predict.fun = predict, h = h, method = "xs-wcpdp",kernel.width=kernel.width,gower.power=gower.power)
  
  xswccpdp <- effect1 (mod = model, data = data, feature = feature, target = "y",
                       predict.fun = predict, h = h, method = "xs-wccpdp",kernel.width=kernel.width,gower.power=gower.power)
  xcwcpdp = effect1(mod = model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=kernel.width,gower.power=gower.power)
  
  xcwccpdp = effect1(mod = model, data = data, feature = feature, target = "y",
                     predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=kernel.width,gower.power=gower.power)
  
  pred <- Predictor$new(model, data = data, y = "y")
  q<-quantile(data[[feature]], 0:h/h)
  pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
  ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
  
  #recording the results in a data frame
  result[,1]=cpdp$y
  result[,2]=xswcpdp$y
  result[,3]=xswccpdp$y
  
  result[,4]=xcwcpdp$y
  result[,5]=xcwccpdp$y
  
  result[,6]=pdp$results$.value
  result[,7]=ale$results$.value
  result[,8]=quantile(data[[feature]], 0:h/h)
  if(model_name!="nnet"){
    result[1:length(model$coefficients),9]=model$coefficients
  }
  
  return(result)
  #one iteration(one result)has 50 rows 
}

do_once_effect2<-function(feature,decay,size,pars){
  h = 50
  result=matrix(0,ncol=9,nrow=h+1)
  data = create_xor_corr(n = 5000)
  X = data[,setdiff(names(data),"y")]
  if(model_name=="nnet"){
    task = as_task_regr(data,target="y")
    lrn = lrn("regr.nnet",size=size,decay=decay)
    model = lrn$train(task = task)$model
  }else{
    model=switch(example_name,
                 "example1"=gam(y ~ x1+I(x2^2)+I(x3^3)+I(x2*x4), data = data),
                 "example2"=gam(y ~ I(x1^2)+x2,data = data),
                 "example3"=gam(y ~ x1+x2+I(x1*x2),data=data),
                 "example4"=gam(y ~ x1+I(x2^2)+I(x3^3)+x3+I(x2*x4),data=data),
                 "example5"=gam(y ~ I(x1^2)+x2,data=data),
                 "example6"=gam(y ~ I(x1^2)+x2,data=data)
    )
  }
  kernel.width=switch(feature,
                      "x1"=pars[1,3],
                      "x2"=pars[2,3],
                      "x3"=pars[3,3],
                      "x4"=pars[4,3],
                      "x5"=pars[5,3])
  gower.power=switch(feature,
                     "x1"=pars[1,4],
                     "x2"=pars[2,4],
                     "x3"=pars[3,4],
                     "x4"=pars[4,4],
                     "x5"=pars[5,4])
  cpdp <- effect2 (mod = model, data = data, feature = feature, target = "y", 
                   predict.fun = predict, h = h, method = "cpdp",kernel.width=kernel.width,gower.power=gower.power)

  xswcpdp <-effect2 (mod = model, data = data, feature = feature, target = "y",
                     predict.fun = predict, h = h, method = "xs-wcpdp",kernel.width=kernel.width,gower.power=gower.power)
  
  xcwcpdp = effect2(mod = model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xc-wcpdp",kernel.width=kernel.width,gower.power=gower.power)
  
  xcwccpdp = effect2(mod = model, data = data, feature = feature, target = "y",
                     predict.fun = predict, h = h, method = "xc-wccpdp",kernel.width=kernel.width,gower.power=gower.power)
  
  xswccpdp <- effect2 (mod = model, data = data, feature = feature, target = "y",
                       predict.fun = predict, h = h, method = "xs-wccpdp",kernel.width=kernel.width,gower.power=gower.power)
  
  pred <- Predictor$new(model, data = data, y = "y")
  q<-quantile(data[[feature]], 0:h/h)
  pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
  ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
  
  result[,1]=cpdp$y
  result[,2]=xswcpdp$y
  result[,3]=xswccpdp$y
  
  result[,4]=xcwcpdp$y
  result[,5]=xcwccpdp$y
  
  result[,6]=pdp$results$.value
  result[,7]=ale$results$.value
  result[,8]=quantile(data[[feature]], 0:h/h)
  if(model_name!="nnet"){
    result[1:length(model$coefficients),9]=model$coefficients
  }
  return(result)
  #one iteration(one result)has 50 rows 
}

do_once_effect3<-function(feature,decay,size,pars){
  h = 50
  result=matrix(0,ncol=7,nrow=h+1)
  data = create_xor_corr(n = 5000)
  X = data[,setdiff(names(data),"y")]
  if(model_name=="nnet"){
    task = as_task_regr(data,target="y")
    lrn = lrn("regr.nnet",size=size,decay=decay)
    model = lrn$train(task = task)$model
  }else{
    model=switch(example_name,
                 "example1"=gam(y ~ x1+I(x2^2)+I(x3^3)+I(x2*x4), data = data),
                 "example2"=gam(y ~ I(x1^2)+x2,data = data),
                 "example3"=gam(y ~ x1+x2+I(x1*x2),data=data),
                 "example4"=gam(y ~ x1+I(x2^2)+I(x3^3)+x3+I(x2*x4),data=data),
                 "example5"=gam(y ~ I(x1^2)+x2,data=data),
                 "example6"=gam(y ~ I(x1^2)+x2,data=data)
    )
  }
  gamma=switch(feature,
               "x1"=pars[1,5],
               "x2"=pars[2,5],
               "x3"=pars[3,5],
               "x4"=pars[4,5],
               "x5"=pars[5,5])
  cpdp <- effect3(mod = model, data = data, feature = feature, target = "y",
                  predict.fun = predict, h = h, method = "cpdp",gamma=gamma)

  xsxcpdp <-effect3(mod = model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=gamma)
  
  cxsxcpdp <-effect3(mod = model, data = data, feature = feature, target = "y",
                     predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma=gamma)
  
  pred <- Predictor$new(model, data = data, y = "y")
  q<-quantile(data[[feature]], 0:h/h)
  pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
  ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
  
  result[,1]=cpdp$y
  result[,2]=xsxcpdp$y
  result[,3]=cxsxcpdp$y
  result[,4]=pdp$results$.value
  result[,5]=ale$results$.value
  result[,6]=quantile(data[[feature]], 0:h/h)
  if(model_name!="nnet"){
    result[1:length(model$coefficients),7]=model$coefficients
  }
  #result[1:5,6]=model$coefficients
  return(result)
  #one iteration(one result)has 50 rows 
}

do_once_effect4<-function(feature,decay,size,pars){
  h = 50
  result=matrix(0,ncol=7,nrow=h+1)
  data = create_xor_corr(n = 5000)
  X = data[,setdiff(names(data),"y")]
  if(model_name=="nnet"){
    task = as_task_regr(data,target="y")
    lrn = lrn("regr.nnet",size=size,decay=decay)
    model = lrn$train(task = task)$model
  }else{
    model=switch(example_name,
                 "example1"=gam(y ~ x1+I(x2^2)+I(x3^3)+I(x2*x4), data = data),
                 "example2"=gam(y ~ I(x1^2)+x2,data = data),
                 "example3"=gam(y ~ x1+x2+I(x1*x2),data=data),
                 "example4"=gam(y ~ x1+I(x2^2)+I(x3^3)+x3+I(x2*x4),data=data),
                 "example5"=gam(y ~ I(x1^2)+x2,data=data),
                 "example6"=gam(y ~ I(x1^2)+x2,data=data)
    )
  }
  gamma=switch(feature,
               "x1"=pars[1,6],
               "x2"=pars[2,6],
               "x3"=pars[3,6],
               "x4"=pars[4,6],
               "x5"=pars[5,6])
  
  cpdp <- effect4 (mod = model, data = data, feature = feature, target = "y", 
                   predict.fun = predict, h = h, method = "cpdp",gamma = gamma)

  xsxcpdp = effect4(mod = model, data = data, feature = feature, target = "y", 
                    predict.fun = predict, h = h, method = "xs-xc-pdp",gamma = gamma)
  
  cxsxcpdp = effect4(mod = model, data = data, feature = feature, target = "y", 
                     predict.fun = predict, h = h, method = "xs-xc-cpdp",gamma = gamma)
  pred <- Predictor$new(model, data = data, y = "y")
  q<-quantile(data[[feature]], 0:h/h)
  pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
  ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
  
  result[,1]=cpdp$y
  result[,2]=xsxcpdp$y
  result[,3]=cxsxcpdp$y
  result[,4]=pdp$results$.value
  result[,5]=ale$results$.value
  result[,6]=quantile(data[[feature]], 0:h/h)
  if(model_name!="nnet"){
    result[1:length(model$coefficients),7]=model$coefficients
  }
  return(result)
  #one iteration(one result)has 50 rows 
}

do_once_effect5<-function(feature,decay,size,pars){
  h = 50
  result=matrix(0,ncol=7,nrow=h+1)
  data = create_xor_corr(n = 5000)
  X = data[,setdiff(names(data),"y")]
  if(model_name=="nnet"){
    task = as_task_regr(data,target="y")
    lrn = lrn("regr.nnet",size=size,decay=decay)
    model = lrn$train(task = task)$model
  }else{
    model=switch(example_name,
                 "example1"=gam(y ~ x1+I(x2^2)+I(x3^3)+I(x2*x4), data = data),
                 "example2"=gam(y ~ I(x1^2)+x2,data = data),
                 "example3"=gam(y ~ x1+x2+I(x1*x2),data=data),
                 "example4"=gam(y ~ x1+I(x2^2)+I(x3^3)+x3+I(x2*x4),data=data),
                 "example5"=gam(y ~ I(x1^2)+x2,data=data),
                 "example6"=gam(y ~ I(x1^2)+x2,data=data)
    )
  }
  kernel.width=switch(feature,
                      "x1"=pars[1,3],
                      "x2"=pars[2,3],
                      "x3"=pars[3,3],
                      "x4"=pars[4,3],
                      "x5"=pars[5,3])
  gower.power=switch(feature,
                     "x1"=pars[1,4],
                     "x2"=pars[2,4],
                     "x3"=pars[3,4],
                     "x4"=pars[4,4],
                     "x5"=pars[5,4])
  cpdp <- effect5 (mod = model, data = data, feature = feature, target = "y", 
                   predict.fun = predict, h = h, method = "cpdp",kernel.width=kernel.width,gower.power=gower.power)

  xsxcpdp = effect5(mod = model, data = data, feature = feature, target = "y", 
                    predict.fun = predict, h = h, method = "xs-xc-pdp",kernel.width=kernel.width,gower.power=gower.power)
  
  cxsxcpdp = effect5(mod = model, data = data, feature = feature, target = "y", 
                     predict.fun = predict, h = h, method = "xs-xc-cpdp",kernel.width=kernel.width,gower.power=gower.power)
  pred <- Predictor$new(model, data = data, y = "y")
  q<-quantile(data[[feature]], 0:h/h)
  pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
  ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
  
  result[,1]=cpdp$y
  result[,2]=xsxcpdp$y
  result[,3]=cxsxcpdp$y
  result[,4]=pdp$results$.value
  result[,5]=ale$results$.value
  result[,6]=quantile(data[[feature]], 0:h/h)
  if(model_name!="nnet"){
    result[1:length(model$coefficients),7]=model$coefficients
  }
  return(result)
  #one iteration(one result)has 50 rows 
}
