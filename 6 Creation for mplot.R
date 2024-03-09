core_number = 8
set_seed_number = 41
iteration_number = 112
example_name="example5"

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
        #kernel.width=1.4
        },
        "example5" = {create_xor_corr = function(n, seed){
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
        "example6" = {create_xor_corr = function(n, seed){
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

mplot = function(data, feature, target, eps) {
  x = data[, feature]
  y = data[, target]
  x.lower = x - eps
  x.upper = x + eps
  m = n = numeric(length(x))
  for(i in 1:length(x)) {
    ind = x > x.lower[i] & x < x.upper[i]
    m[i] = mean(y[ind])
    n[i] = sum(ind)
  }
  return(data.frame(x = x, mplot = m, n = n))
}

do_once=function(feature){
  h = 50
  data = create_xor_corr(n = 5000)
  mpl = mplot(data=data, feature=feature, target = "y", eps = 0.5)
  return(mpl)
}

for(i in feature){
  {
    cl<- makeCluster(core_number)      
    registerDoParallel(cl)       
    mydata<- foreach(
      j=1:iteration_number,          
      .combine=rbind,  
      .packages = c("iml","R.utils","mvtnorm","patchwork","data.table",
                    "StatMatch","dplyr","mlr3","mlr3verse","mlr3learners","nnet","MASS")
    ) %dopar% {
      set.seed(set_seed_number+j)
      do_once(i)
    }
    stopImplicitCluster()
    stopCluster(cl)
  }
  mydata=as.data.frame(mydata)
  write.csv(mydata,file=paste0(example_name,"_",i,"mplot.csv"))
}

ite=112
example_name="example5"
feature="x1"
mplot=read.csv(paste0(example_name,"_",feature,"mplot.csv"))

