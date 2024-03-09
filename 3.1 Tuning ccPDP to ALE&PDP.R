#Parameters Configuration
example_name="example4"
core_number=6
set_seed_number=41
# plot_id=sample(nrow(kc_housing1),0.1*nrow(kc_housing1))
# kc_housing2=kc_housing1[-plot_id,]
# kc_housing3=kc_housing1[plot_id,]
h=20
# Generating data ---------------------------------------------------------
switch (example_name,
        "example1" = {create_xor_corr = function(n, seed){
          x1 = runif(n, -1, 1)
          x2 = runif(n, -1, 1)
          x3 = runif(n, -1, 1)
          x4 = runif(n, -1, 1)
          x5 = runif(n, -1, 1)
          
          y = x1+x2^2+x3^3+0.8*x2*x4+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3,x4,x5, y)
        }
        size=20
        decay=0.02192304
        feature = c("x1","x2","x3","x4","x5")
        #kernel.width=0.8
        },
        "example2" = {create_xor_corr = function(n, seed){
          x1 = runif(n, -1, 1)
          x2 = 0.8*x1+rnorm(n, sd = 0.1)
          x3 = -x1+rnorm(n, sd = 0.1)
          
          y = x1^2+x2+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3, y)
        }
        size=4
        decay=0.01495741
        feature = c("x1","x2","x3")
        #kernel.width=2
        },
        "example3"={create_xor_corr = function(n, seed){
          x1 = runif(n, -1, 1)
          x3 = runif(n, -1, 1)
          
          x2 = -x1+rnorm(n, sd = 0.1)
          
          y = x1+x2+x1*x2+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3, y)
        }
        size=5
        decay=0.006384976
        feature = c("x1","x2","x3")
        #kernel.width=2
        },
        "example4"={create_xor_corr = function(n, seed){
          x1 = runif(n, -1, 1)
          x2 = runif(n, -1, 1)
          x3 = runif(n, -1, 1)
          
          #1 low corr 0.2822673:
          #x4 = 0.05*x2+rnorm(n, sd = 0.1)
          #size=18
          #decay=0.005944008
          
          #2 middle corr 0.548669:
          x4 = 0.1*x2+rnorm(n, sd = 0.1)
          #size=8
          #decay=0.004630994
          
          #3 high corr 0.9863376:
          #x4 = x2+rnorm(n, sd = 0.1)
          #size=19
          #decay=0.01166048
          
          x5 = -x3+rnorm(n, sd = 0.1)
          
          y = x1+0.5*(3*(x2^2)-1)+0.5*(4*(x3^3)-3*x3)+0.8*x2*x4+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3,x4,x5,y)
        }
        
        #1 low corr 0.2822673:
        # size=18
        # decay=0.005944008
        
        #2 middle corr 0.548669:
        size=8
        decay=0.004630994
        
        #3 high corr 0.9863376:
        #size=19
        #decay=0.01166048
        
        feature = c("x1","x2","x3","x4","x5")
        #kernel.width=1.4
        },
        "example5" = {create_xor_corr = function(n, seed){
          x1 = rexp(n,rate=0.5)
          x2 = 0.8*x1+rnorm(n, sd = 0.1)
          x3 = -x1+rnorm(n, sd = 0.1)
          
          y = x1^2+x2+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3, y)
        }
        size=4
        decay=0.001328558
        feature = c("x1","x2","x3")
        },
        "example6" = {create_xor_corr = function(n, seed){
          mean=c(1,2,3)
          sigma=matrix(c(1,0.8,0.2,0.8,1,0,0.2,0,1),nrow=3,ncol=3)
          d=mvrnorm(n,mean,sigma)
          x1=d[,1]
          x2=d[,2]
          x3=d[,3]
          y = x1^2+x2+rnorm(n, sd = 0.1)
          #x1 and x2 are highly correlated, x3 is weakly correlated with x1,x2
          data.frame(x1,x2,x3, y)
        }
        size=12
        decay=0.00577085
        feature = c("x1","x2","x3")
        }
        
)

# losses ---------------------------------------------------------

#loss for effect1 w.r.t. kernel.width
custom_loss1=function(xs,model,data,feature_i,h,ale) {
  #calculating effect1-xswcpdp with given parameters
  curve=effect1(mod = model, data = data, feature = feature_i, target ="price",
                predict.fun = predict, h = h, method = "xs-wcpdp",gower.power=1,kernel.width=xs)
  a=curve$y[!is.na(curve$y)]            
  a=a-mean(a,na.rm=TRUE)
  
  #selecting whether to push the curve to ALE or cenetred PDP
  b=ale$results$.value[1:length(a)]
  #b=pdp$results$.value[1:length(a)]
  #b=b-mean(b)
  
  #returning the MSE to the selected curve
  return(sum((a-b)^2)/length(a))
}

#loss for effect1 w.r.t. gower.power
custom_loss2 = function(xs,model,data,feature_i,h,ale) {
  #calculating effect1-xcwcpdp with given parameters
  curve=effect1(mod = model, data = data, feature = feature_i, target ="price",
                predict.fun = predict, h = h, method = "xc-wcpdp",gower.power=xs,kernel.width=0.6)
  a=curve$y[!is.na(curve$y)]
  a=a-mean(a,na.rm=TRUE)
  
  #selecting whether to push the curve to ALE or cenetred PDP
  b=ale$results$.value[1:length(a)]
  #b=pdp$results$.value[1:length(a)]
  #b=b-mean(b)
  
  #returning the MSE to the selected curve
  return(sum((a-b)^2)/length(a))
}

#loss for effect2 w.r.t. kernel.width
custom_loss3 = function(xs,model,data,feature_i,h,ale) {
  #calculating effect2-xswcpdp with given parameters
  curve=effect2(mod = model, data = data, feature = feature_i, target ="price",
                predict.fun = predict, h = h, method = "xs-wcpdp",gower.power=1,kernel.width=xs)
  a=curve$y[!is.na(curve$y)]
  a=a-mean(a,na.rm=TRUE)
  
  #selecting whether to push the curve to ALE or cenetred PDP
  b=ale$results$.value[1:length(a)]
  #b=pdp$results$.value[1:length(a)]
  #b=b-mean(b)
  
  #returning the MSE to the selected curve
  return(sum((a-b)^2)/length(a))
}

#loss for effect2 w.r.t. gower.power
custom_loss4 = function(xs,model,data,feature_i,h,ale) {
  #calculating effect2-xswcpdp with given parameters
  curve=effect2(mod = model, data = data, feature = feature_i, target ="price",
                predict.fun = predict, h = h, method = "xc-wcpdp",gower.power=xs,kernel.width=0.6)
  a=curve$y[!is.na(curve$y)]
  a=a-mean(a,na.rm=TRUE)
  
  #selecting whether to push the curve to ALE or cenetred PDP
  b=ale$results$.value[1:length(a)]
  #b=pdp$results$.value[1:length(a)]
  #b=b-mean(b)
  
  #returning the MSE to the selected curve
  return(sum((a-b)^2)/length(a))
}

#loss for effect3 w.r.t. gamma
custom_loss5 = function(xs,model,data,feature_i,h,ale) {
  #calculating effect3-xsxcpdp with given parameters
  curve=effect3(mod = model, data = data, feature = feature_i, target ="price",predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=xs)
  a=curve$y[!is.na(curve$y)]
  a=a-mean(a,na.rm=TRUE)
  
  #selecting whether to push the curve to ALE or cenetred PDP
  b=ale$results$.value[1:length(a)]
  #b=pdp$results$.value[1:length(a)]
  #b=b-mean(b)
  
  #returning the MSE to the selected curve
  return(sum((a-b)^2)/length(a))
}

#loss for effect4 w.r.t. gamma
custom_loss6 = function(xs,model,data,feature_i,h,ale) {
  #calculating effect4-xsxcpdp with given parameters
  curve=effect4(mod = model, data = data, feature = feature_i, target ="price",predict.fun = predict, h = h, method = "xs-xc-pdp",gamma=xs)
  a=curve$y[!is.na(curve$y)]
  a=a-mean(a,na.rm=TRUE)
  
  #selecting whether to push the curve to ALE or cenetred PDP
  b=ale$results$.value[1:length(a)]
  #b=pdp$results$.value[1:length(a)]
  #b=b-mean(b)
  
  #returning the MSE to the selected curve
  return(sum((a-b)^2)/length(a))
}

# tuning function ---------------------------------------------------------
tune=function(combi,feature_i){
  if(combi=="eff3_gamma"|combi=="eff4_gamma"){
    #upper bound for gamma
    upper1=20
  }else{
    #upper bound for kernel.width and gower.power
    upper1=2
  }
  
  #selecting corresponding loss function defined above
  custom_loss=switch(combi,
         "eff1_kw"=custom_loss1,
         "eff1_gp"=custom_loss2,
         "eff2_kw"=custom_loss3,
         "eff2_gp"=custom_loss4,
         "eff3_gamma"=custom_loss5,
         "eff4_gamma"=custom_loss6)
  
  #computing ale or centered pdp
  pred <- Predictor$new(model, data = data, y = "price")
  q<-quantile(data[[feature_i]], 0:h/h)
  ale <- FeatureEffect$new(pred, feature = feature_i, method = "ale", grid.points =q)
  #pdp <- FeatureEffect$new(pred, feature = feature_i, method = "pdp", grid.points =q)
  
  #minimizing the loss function
  lgr::get_logger("bbotk")$set_threshold("warn")
  result = bb_optimize(custom_loss, method = "random_search",lower = 0, upper = upper1,  max_evals = 200,
              data =data,model=model,feature_i=feature_i,h=h,ale=ale)
  return(result$par$x1)
         }

#tuning the 5 parameters for all effects in parallel
data_final=data.frame()
s=Sys.time()
for(i in feature){
  {
    cl<- makeCluster(core_number)      
    registerDoParallel(cl)       
    mydata<- foreach(
      j=c("eff1_kw","eff1_gp","eff2_kw","eff2_gp","eff3_gamma","eff4_gamma"),   
      
      #combining the results column by column
      .combine=cbind,  
      .packages = c("iml","R.utils","mvtnorm","patchwork","data.table",
                    "StatMatch","dplyr","mlr3","mlr3verse","mlr3learners",
                    "nnet","mgcv","MASS","bbotk")
    ) %dopar% {
      set.seed(set_seed_number)
      tune(j,i)
      #do_once(i,decay,size,pars) with tuned parameters in pars.csv
    }
    stopImplicitCluster()
    stopCluster(cl)
  }
  mydata=as.data.frame(mydata)
  names(mydata)=c("eff1_kw","eff1_gp","eff2_kw","eff2_gp","eff3_gamma","eff4_gamma")
  rownames(mydata)=i
  write.csv(mydata,file=paste0(i,"_parameters.csv"))
  data_final=rbind(data_final,mydata)
}
write.csv(data_final,file="parameters.csv")
e=Sys.time()
e-s

