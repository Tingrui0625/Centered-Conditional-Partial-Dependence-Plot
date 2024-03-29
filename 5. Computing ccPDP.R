#This file needs to be sourced multiple times for different examples and effects. 
#The general setting should be set and this file should be called after that.

# general setting for single run---------------------------------------------------------
# core_number = 8
# example_name = "example4"
# effect_name = "effect1"
# model_name="nnet"#"gam"
# kernel.width=0.6
# gamma=15
# set_seed_number = 41
# iteration_number = 8

# data generation function---------------------------------------------------------------
switch (example_name,
        "example1" = {create_xor_corr = function(n){
          x1 = runif(n, -1, 1)
          x2 = runif(n, -1, 1)
          x3 = runif(n, -1, 1)
          x4 = runif(n, -1, 1)
          x5 = runif(n, -1, 1)
          
          y = x1+x2^2+x3^3+0.8*x2*x4+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3,x4,x5, y)
        }
        #tuned parameters with example 1 for NN
        size=20 
        decay=0.02192304
        feature = c("x1","x2","x3","x4","x5")
        },
        "example2" = {create_xor_corr = function(n){
          x1 = runif(n, -1, 1)
          x2 = 0.8*x1+rnorm(n, sd = 0.1)
          x3 = -x1+rnorm(n, sd = 0.1)
          y = x1^2+x2+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3, y)
        }
        size=4
        decay=0.01495741
        feature = c("x1","x2","x3")
        },
        "example3"={create_xor_corr = function(n){
          x1 = runif(n, -1, 1)
          x3 = runif(n, -1, 1)
          x2 = -x1+rnorm(n, sd = 0.1)
          y = x1+x2+x1*x2+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3, y)
        }
        size=5
        decay=0.006384976
        feature = c("x1","x2","x3")
        },
        "example4"={create_xor_corr = function(n){
          x1 = runif(n, -1, 1)
          x2 = runif(n, -1, 1)
          x3 = runif(n, -1, 1)
          #Three levels of correlation
          #1 low corr 0.2822673:
          #x4 = 0.05*x2+rnorm(n, sd = 0.1)
          #size=18
          #decay=0.005944008
          
          #2 middle corr 0.548669:
          #x4 = 0.1*x2+rnorm(n, sd = 0.1)
          #size=8
          #decay=0.004630994
          
          #3 high corr 0.9863376:
          x4 = x2+rnorm(n, sd = 0.1)
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
        #size=8
        #decay=0.004630994
        
        #3 high corr 0.9863376:
        size=19
        decay=0.01166048
        feature = c("x1","x2","x3","x4","x5")
        },
        "example5" = {create_xor_corr = function(n){
          x1 = rexp(n,rate=0.5)
          x3 = rexp(n,rate=0.5)
          x2 = -x1+rnorm(n, sd = 0.1)
          y = x1+x2+x1*x2+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3, y)
        }
        size=7
        decay=0.005851096
        feature = c("x1","x2","x3")
        },
        "example6" = {create_xor_corr = function(n){
          mean=c(1,2,3)
          sigma=matrix(c(1,0.9,0.2,0.9,1,0,0.2,0,1),nrow=3,ncol=3)
          d=mvrnorm(n,mean,sigma)
          x1=d[,1]
          x2=d[,2]
          x3=d[,3]
          y = x1+x2+x1*x2+rnorm(n, sd = 0.1)
          data.frame(x1,x2,x3, y)
        }
        size=20
        decay=0.002771715
        feature = c("x1","x2","x3")
        }
        
)

#changing the multisession result names according to effect name
switch (effect_name,
        "effect1" = {do_once=do_once_effect1
        mydata_name = c("cpdp","xswcpdp","xswccpdp","xcwcpdp","xcwccpdp","pdp","ale","grid")#"coeff"
        },
        "effect2" = {
          do_once=do_once_effect2
          mydata_name = c("cpdp","xswcpdp","xcwcpdp","xcwccpdp","xswccpdp","pdp","ale","grid")
        },
        
        "effect3" = {do_once=do_once_effect3
        mydata_name = c("cpdp","xsxcpdp","cxsxcpdp","pdp","ale","grid")},

        "effect4" = {
          do_once=do_once_effect4
          mydata_name = c("cpdp","xsxcpdp","cxsxcpdp","pdp","ale","grid")
        },
        "effect5" = {
          do_once=do_once_effect5
          mydata_name = c("cpdp","xsxcpdp","cxsxcpdp","pdp","ale","grid")
        }
)

#running in parallel: repeating data generation and ccPDP-calcualtion for multiple times
s=Sys.time()
for(i in feature){
  {
    cl<- makeCluster(core_number)      
    registerDoParallel(cl)       
    mydata<- foreach(
      #repeating iteration_number times
      j=1:iteration_number,  
      
      #combinding the results row by row
      .combine=rbind,  
      
      #loading the needed packages into the interval enviroment
      .packages = c("iml","R.utils","mvtnorm","patchwork","data.table",
                    "StatMatch","dplyr","mlr3","mlr3verse","mlr3learners","nnet","mgcv","MASS")
    ) %dopar% {
      
      #setting the starting seed and making the seed of multisessions different
      set.seed(set_seed_number+j)
      
      #running the function calculating ccPDP
      do_once(i,decay,size)
      #do_once(i,decay,size,pars) #with tuned parameters in pars.csv
      }
    stopImplicitCluster()
    stopCluster(cl)
  }
  mydata=as.data.frame(mydata)
  
  #changing the result names to the earlier initialized "mydata_name"
  names(mydata)=mydata_name
  
  #saving the ccPDP results
  write.csv(mydata,file=paste0(effect_name,"_",example_name,"_",i,".csv"))
}
e=Sys.time()
e-s




