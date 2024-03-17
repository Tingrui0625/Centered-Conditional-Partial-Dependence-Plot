# 1. specialized effect ------------------------------------------------------
#rewriting the ccPDP methods to put all the data points in the calculation
#changing quantile values to the original xs feature values
effect1 = function(mod, data, feature, target, kernel.width, gower.power = 1, predict.fun = predict, h = 999,method){
  gower.power=as.numeric(gower.power)
  kernel.width=as.numeric(kernel.width)
  #replacing the xs with all xs values instead of only quantiles
  ice_vals = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), data[[feature]]))#n row h columns
  
  #xs-proximity-based weight
  xs.weight = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {#lapply: each element as a list, sapply as numeric values
    f.val = data[[feature]]
    xs.dist = exp(-0.5 * ((grid - f.val)^2) / (kernel.width^2))
  }), data[[feature]]))
  
  #xc-proximity-based weight
  xc = data[, setdiff(colnames(data), c(feature, target)), drop = FALSE]
  xc.weight = as.data.frame(1 - as.matrix(cluster::daisy(xc, metric = "gower"))^gower.power)
  
  #index for additional conditioning
  if(method=="xs-wccpdp"|method=="xc-wccpdp"){
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    #the first interval
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    #the last interval
    m2=length(which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1]))
    if(m2>0){
      ind_c[1:m2,h+1]=which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1])
    }else if(q[h]==q[h+1]){
      ind_c[1:length(which(data[[feature]] == q[h+1])),h+1]=which(data[[feature]] == q[h+1])
    }
    
    #the middle intervals
    for(i in 2:h){
      mi=length(which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2))
      if(mi>0){
        ind_c[1:mi,i]=which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2)
      }else if((q[i-1]+q[i])/2==(q[i+1]+q[i])/2){
        ind_c[1:which(data[[feature]] == (q[i-1]+q[i])/2),i]=which(data[[feature]] == (q[i-1]+q[i])/2)
      }
    }
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  weighted_ice_vals=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
  switch(method,
         "xs-wcpdp"={for(i in 1:(h+1)){weighted_ice_vals[,i]=ice_vals[,i] * xs.weight[,i]/sum(xs.weight[,i],na.rm=TRUE)}},
         
         "xs-wccpdp"={for(i in 1:(h+1)){
           if(sum(xs.weight[index[,i], i],na.rm=TRUE)==0){
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=0
           }else{
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=ice_vals[index[,i], i] * xs.weight[index[,i], i]/sum(xs.weight[index[,i], i],na.rm=TRUE)   
           }}},
         
         "xc-wcpdp"={for(i in 1:(h+1)){weighted_ice_vals[,i]=ice_vals[,i] * xc.weight[,i]/sum(xc.weight[,i],na.rm=TRUE)}},#colMeans(xc.weight)/sum(colMeans(xc.weight))
         
         "xc-wccpdp"={for(i in 1:(h+1)){
           if(length(which(index[,i]==TRUE))==1){
             weighted_ice_vals[1,i]=ice_vals[index[,i], i]
             
           }else if(sum(xc.weight[index[,i], i],na.rm=TRUE)==0){
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=0
             
           }else{weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=ice_vals[index[,i], i] * xc.weight[index[,i], i]/sum(xc.weight[index[,i], i],na.rm=TRUE)
           }
           #old: colMeans(xc.weight[index[,i], index[,i]])/sum(colMeans(xc.weight[index[,i], index[,i]]))
         }},
         
         "cpdp"={weighted_ice_vals=ice_vals/nrow(ice_vals)})
  
  weighted_ice_vals=weighted_ice_vals-rowMeans(weighted_ice_vals,na.rm=TRUE)%*%matrix(1,nrow=1,ncol=h+1)
  d = data.frame(x = data[[feature]], y = colSums(weighted_ice_vals,na.rm=TRUE))
  return(d)
}
effect2 = function(mod, data, feature, target, kernel.width, gower.power = 1, predict.fun = predict, h = 999, method){
  gower.power=as.numeric(gower.power)
  kernel.width=as.numeric(kernel.width)
  # Interval bounds
  #q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs
  # compute standard ICE values at quantile grid points
  
  ice_vals = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), data[[feature]]))#n row h columns
  
  # center each ICE curve by subtracting its mean
  #exp_ice_vals = rowMeans(ice_vals)
  ice_vals=ice_vals-rowMeans(ice_vals)%*%matrix(1,nrow=1,ncol=h+1) #exp_ice_vals[i]
  #center over xs, xs taking different values
  
  # compute weights based on similarity of observed Xs-value to grid point
  xs.weight = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {#lapply: each element as a list, sapply as numeric values
    f.val = data[[feature]]
    # scaled L1-distance as done in the gower distance, kernel.width
    #xs.dist = 1 - (abs(grid - f.val)/diff(range(f.val)))^kernel.width
    # L2-distance turned into gaussian kernel
    xs.dist = exp(-0.5 * ((grid - f.val)^2) / (kernel.width^2))#previous:exp(-(grid - f.val)^2/(2*kernel.width))
  }), data[[feature]]))##the distance from all values xs to 20 grid points xs values
  
  
  xc = data[, setdiff(colnames(data), c(feature, target)), drop = FALSE]#setdiff to find different elements
  #xc.weight = 1 - as.matrix(cluster::daisy(xc, metric = "gower"))^(as.numeric(gower.power))
  xc.weight =1-(as.matrix(cluster::daisy(xc, metric = "gower")))**gower.power
  #daisy() is to calculate distance for xc
  if(method=="xs-wccpdp"|method=="xc-wccpdp"){
    
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    m2=length(which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1]))
    if(m2>0){
      ind_c[1:m2,h+1]=which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1])
    }else if(q[h]==q[h+1]){
      ind_c[1:length(which(data[[feature]] == q[h+1])),h+1]=which(data[[feature]] == q[h+1])
    }
    
    for(i in 2:h){
      mi=length(which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2))
      if(mi>0){
        ind_c[1:mi,i]=which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2)
      }else if((q[i-1]+q[i])/2==(q[i+1]+q[i])/2){
        ind_c[1:which(data[[feature]] == (q[i-1]+q[i])/2),i]=which(data[[feature]] == (q[i-1]+q[i])/2)
      }
    }
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  
  # ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
  # ind_c[1:length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)),1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
  # ind_c[1:length(which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1])),h+1]=which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1])
  # for(i in 2:h){
  #   ind_c[1:length(which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2)),i]=which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2)
  # }
  # indicator2origin<-rep(0,h+1)
  # for(k in 1:(h+1)){
  #   indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  # }
  # xc.weight1=matrix(0,nrow=nrow(data),ncol=(h+1))
  # for(i in 1:(h+1)){
  #   xc.weight1[,i]=xc.weight[,indicator2origin[i]]
  # }
  # xc.weight1=as.data.frame(xc.weight1)
  y = vapply(1:ncol(ice_vals), function(i) {
    # index = data[[feature]] >= q[i] & data[[feature]] <= q[i + 1]#for each observation,feature is within quantile interval
    # if (!any(index, na.rm = TRUE)){index = NULL}
    switch(method, 
           # weighted (based on xs) centered pdp
           "xs-wcpdp" = weighted.mean(ice_vals[, i], xs.weight[, i]), ##use all the weights for all points
           
           # weighted (based on xs) conditional (based on obs within quantile interval) centered pdp
           "xs-wccpdp" = {if(sum(xs.weight[index[,i], i],na.rm=TRUE)==0){
             0
           }else{weighted.mean(ice_vals[index[,i], i], xs.weight[index[,i], i],na.rm=TRUE)}},
           ##If xs(i) is within one interval,then use weights for this interval
           
           # weighted (based on xc by averaging weights within quantile interval) conditional (based on obs within quantile interval) centered pdp
           "xc-wcpdp" =weighted.mean(ice_vals[, i], xc.weight[,i]),
           #old: colMeans(xc.weight)),##mean of distance from xc to other xc
           
           # weighted (based on xc by averaging weights within quantile interval) conditional (based on obs within quantile interval) centered pdp
           "xc-wccpdp" = {if(sum(xc.weight[index[,i], i],na.rm=TRUE)==0){
             0
           }else{weighted.mean(ice_vals[index[,i], i], xc.weight[index[,i],i],na.rm=TRUE)}},
           
           #old: weighting all column with the same series of values
           #weighted.mean(ice_vals[index, i], colMeans(xc.weight[index, index]))
           
           # centered pdp
           "cpdp" = mean(ice_vals[[i]]) # this would result in a centered PDP over xc without weights
           
    )
    
  }, FUN.VALUE = NA_real_)
  d = data.frame(x = data[[feature]], y = y)##-mean(y,na.rm=TRUE)
  return(d)
}
effect3 = function(mod, data, feature, target, gamma, predict.fun = predict, h = 999,method){
  gamma=as.numeric(gamma)
  # Interval bounds
  #q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs
  # compute standard ICE values at quantile grid points
  ice_vals = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), data[[feature]]))#n row h columns
  # ICE curve
  
  # indicator2origin<-rep(0,(h+1))
  # for(k in 1:(h+1)){
  #   # if(length(which(data[,feature]==q[k]))!=0){
  #   #   indicator2origin[k]<-which(data[,feature]==q[k])
  #   # }
  #   # else{
  #   indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  #   #   
  #   #   (which(abs(data[,feature]-q[k])==min(abs(data[,feature]-q[k])))[1])}
  # }
  # distance<-mahalanobis.dist(data[, setdiff(colnames(data), c(target)), drop = FALSE])
  # W_kernel_full <- exp(-0.5 * (distance^2) / fixed_sigma^2)
  
  X<-data[, setdiff(colnames(data), target), drop = FALSE]
  W_kernel_full=exp(-gamma*(as.matrix(cluster::daisy(X, metric = "gower")))^2)
  
  if(method=="xs-xc-cpdp"){
    
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    m2=length(which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1]))
    if(m2>0){
      ind_c[1:m2,h+1]=which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1])
    }else if(q[h]==q[h+1]){
      ind_c[1:length(which(data[[feature]] == q[h+1])),h+1]=which(data[[feature]] == q[h+1])
    }
    
    for(i in 2:h){
      mi=length(which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2))
      if(mi>0){
        ind_c[1:mi,i]=which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2)
      }else if((q[i-1]+q[i])/2==(q[i+1]+q[i])/2){
        ind_c[1:which(data[[feature]] == (q[i-1]+q[i])/2),i]=which(data[[feature]] == (q[i-1]+q[i])/2)
      }
    }
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  
  weighted_ice_vals=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
  switch(method,
         "xs-xc-pdp"={for(i in 1:(h+1)){weighted_ice_vals[,i]=ice_vals[,i] * W_kernel_full[,i]/sum(W_kernel_full[,i],na.rm=TRUE)}},
         
         "xs-xc-cpdp"={for(i in 1:(h+1)){
           if(sum(W_kernel_full[index[,i],i],na.rm=TRUE)==0){
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=0
           }else{
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=ice_vals[index[,i], i] * W_kernel_full[index[,i],i]/sum(W_kernel_full[index[,i],i],na.rm=TRUE)}
         }},
         
         "cpdp"={weighted_ice_vals=ice_vals/nrow(ice_vals)})
  
  weighted_ice_vals=weighted_ice_vals-rowMeans(weighted_ice_vals,na.rm=TRUE)%*%matrix(1,nrow=1,ncol=h+1)
  
  d = data.frame(x = data[[feature]], y = colSums(weighted_ice_vals,na.rm=TRUE))
  
  return(d)
}
effect4 = function(mod, data, feature, target, gamma, predict.fun = predict, h = 999,method){
  gamma=as.numeric(gamma)
  # Interval bounds
  #q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs
  #q=ale$results[,feature]
  # compute standard ICE values at quantile grid points
  ice_vals = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), data[[feature]]))#n row h columns
  
  ice_vals=ice_vals-rowMeans(ice_vals)%*%matrix(1,nrow=1,ncol=h+1)
  
  # indicator2origin<-rep(0,h+1)
  # for(k in 1:(h+1)){
  #   indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  # }
  # distance<-mahalanobis.dist(data[, setdiff(colnames(data), c(target)), drop = FALSE])
  # W_kernel_full <- exp(-0.5 * (distance^2) / fixed_sigma^2)
  
  #W_kernel_full=1-gower.dist(data)
  X<-data[, setdiff(colnames(data), target), drop = FALSE]
  W_kernel_full=exp(-gamma*(as.matrix(cluster::daisy(X, metric = "gower")))^2)
  #when gamma=0.2, almost all the weights are 0.9!!!wrong!!!(that is to take all points within an interval)
  if(method=="xs-xc-cpdp"){
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    m2=length(which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1]))
    if(m2>0){
      ind_c[1:m2,h+1]=which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1])
    }else if(q[h]==q[h+1]){
      ind_c[1:length(which(data[[feature]] == q[h+1])),h+1]=which(data[[feature]] == q[h+1])
    }
    
    for(i in 2:h){
      mi=length(which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2))
      if(mi>0){
        ind_c[1:mi,i]=which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2)
      }else if((q[i-1]+q[i])/2==(q[i+1]+q[i])/2){
        ind_c[1:which(data[[feature]] == (q[i-1]+q[i])/2),i]=which(data[[feature]] == (q[i-1]+q[i])/2)
      }
    }
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  
  y = vapply(1:ncol(ice_vals), function(i) {
    # index = data[[feature]] >= q[i] & data[[feature]] <= q[i + 1]#for each observation,feature is within quantile interval
    # if (!any(index, na.rm = TRUE)) 
    #   index = NULL
    switch(method, 
           # weighted (based on xs) conditional (based on obs within quantile interval) centered pdp
           "xs-xc-cpdp"= {
             if(sum(W_kernel_full[index[,i], i],na.rm=TRUE)==0){
               0
             }else{weighted.mean(ice_vals[index[,i], i], W_kernel_full[index[,i], i],na.rm=TRUE)}
           },##If xs(i) is within one interval,then use weights for this interval
           #good
           # weighted (based on xs) centered pdp
           "xs-xc-pdp" = weighted.mean(ice_vals[, i], W_kernel_full[, i]), ##use all the weights for all points
           # centered pdp
           "cpdp" = mean(ice_vals[[i]]), # this would result in a centered PDP over xc without weights
    )
    
  }, FUN.VALUE = NA_real_)
  d = data.frame(x = data[[feature]], y = y)
  
  return(d)
}
effect5 = function(mod, data, feature, target, kernel.width, gower.power = 1, predict.fun = predict, h = 999, method){
  gower.power=as.numeric(gower.power)
  kernel.width=as.numeric(kernel.width)
  q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs
  
  ice_vals = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), data[[feature]]))
  ice_vals=ice_vals-rowMeans(ice_vals)%*%matrix(1,nrow=1,ncol=h+1)
 
  xs.weight = data.table::setDF(setNames(lapply(data[[feature]], function(grid) {#lapply: each element as a list, sapply as numeric values
    f.val = data[[feature]]
    # scaled L1-distance as done in the gower distance, kernel.width
    #xs.dist = 1 - (abs(grid - f.val)/diff(range(f.val)))^kernel.width
    # L2-distance turned into gaussian kernel
    xs.dist = exp(-0.5 * ((grid - f.val)^2) / (kernel.width^2))#previous:exp(-(grid - f.val)^2/(2*kernel.width))
  }), data[[feature]]))##the distance from all values xs to 20 grid points xs values
  
  
  xc = data[, setdiff(colnames(data), c(feature, target)), drop = FALSE]#setdiff to find different elements
  #xc.weight = 1 - as.matrix(cluster::daisy(xc, metric = "gower"))^(as.numeric(gower.power))
  xc.weight =1-(as.matrix(cluster::daisy(xc, metric = "gower")))**gower.power
  
  
  if(method=="xs-xc-cpdp"){
    
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    m2=length(which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1]))
    if(m2>0){
      ind_c[1:m2,h+1]=which(data[[feature]] >= (q[h]+q[h+1])/2 & data[[feature]] <= q[h+1])
    }else if(q[h]==q[h+1]){
      ind_c[1:length(which(data[[feature]] == q[h+1])),h+1]=which(data[[feature]] == q[h+1])
    }
    
    for(i in 2:h){
      mi=length(which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2))
      if(mi>0){
        ind_c[1:mi,i]=which(data[[feature]] >= (q[i-1]+q[i])/2 & data[[feature]] < (q[i+1]+q[i])/2)
      }else if((q[i-1]+q[i])/2==(q[i+1]+q[i])/2){
        ind_c[1:which(data[[feature]] == (q[i-1]+q[i])/2),i]=which(data[[feature]] == (q[i-1]+q[i])/2)
      }
    }
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  y = vapply(1:ncol(ice_vals), function(i) {
    switch(method, 
           "xs-xc-pdp" = weighted.mean(ice_vals[, i], xs.weight[, i]*xc.weight[,i]),##mean of distance from xc to other xc
           
           "xs-xc-cpdp" = {
             if(sum(xs.weight[index[,i], i]*xc.weight[index[,i],i],na.rm=TRUE)==0){
               0
             }else{
               weighted.mean(ice_vals[index[,i], i], xs.weight[index[,i], i]*xc.weight[index[,i],i],na.rm=TRUE)
             }},
           
           "cpdp" = mean(ice_vals[[i]]) # this would result in a centered PDP over xc without weights
           
    )
    
  }, FUN.VALUE = NA_real_)
  d = data.frame(x = q, y = y)##-mean(y,na.rm=TRUE)
  return(d)
}

# 2. losses(model fidelity) ------------------------------------------------------------------
#loss functions to calculate the discrepancies to the prediction values of all data points
#n=1000 observations are generated 
feature_n=switch(example_name,
                 "example3"=3,
                 "example4"=5,
                 "case1"=11)
custom_loss1=function(xs,model,data){
  curve = data.table::setDF(setNames(lapply(1:feature_n, function(i) {
    effect1 (mod = model, data = data, feature = names(data)[i], target = "quality",
             predict.fun = predict, h = 999, method = "xs-wcpdp",gower.power = 10,kernel.width=xs[[i]])$y
  }), 1:feature_n))
  pre=as.numeric(model$fitted.values)
  pre=pre-mean(pre)
  sum((rowSums(curve)-as.numeric(pre))^2)/1000
}
custom_loss2=function(xs,model,data){
  curve = data.table::setDF(setNames(lapply(1:feature_n, function(i) {
    effect1 (mod = model, data = data, feature = names(data)[i], target = "quality",
             predict.fun = predict, h = 999, method = "xc-wcpdp",kernel.width=0.6,gower.power=xs[[i]])$y
  }), 1:feature_n))
  pre=as.numeric(model$fitted.values)
  pre=pre-mean(pre)
  sum((rowSums(curve)-as.numeric(pre))^2)/1000
}
custom_loss3=function(xs,model,data){
  curve = data.table::setDF(setNames(lapply(1:feature_n, function(i) {
    effect2 (mod = model, data = data, feature = names(data)[i], target = "quality",
             predict.fun = predict, h = 999, method = "xs-wcpdp",gower.power = 10,kernel.width=xs[[i]])$y
  }), 1:feature_n))
  pre=as.numeric(model$fitted.values)
  pre=pre-mean(pre)
  sum((rowSums(curve)-as.numeric(pre))^2)/1000
}
custom_loss4=function(xs,model,data){
  curve = data.table::setDF(setNames(lapply(1:feature_n, function(i) {
    effect2 (mod = model, data = data, feature = names(data)[i], target = "quality",
             predict.fun = predict, h = 999, method = "xc-wcpdp",kernel.width=0.6,gower.power=xs[[i]])$y
  }), 1:feature_n))
  pre=as.numeric(model$fitted.values)
  pre=pre-mean(pre)
  sum((rowSums(curve)-as.numeric(pre))^2)/1000
}
custom_loss5=function(xs,model,data){
  curve = data.table::setDF(setNames(lapply(1:feature_n, function(i) {
    effect3 (mod = model, data = data, feature = names(data)[i], target = "quality",
             predict.fun = predict, h = 999, method = "xs-xc-pdp",gamma=xs[[i]])$y
  }), 1:feature_n))
  pre=as.numeric(model$fitted.values)
  pre=pre-mean(pre)
  sum((rowSums(curve)-as.numeric(pre))^2)/1000
}
custom_loss6=function(xs,model,data){
  curve = data.table::setDF(setNames(lapply(1:feature_n, function(i) {
    effect4 (mod = model, data = data, feature = names(data)[i], target = "quality",
             predict.fun = predict, h = 999, method = "xs-xc-pdp",gamma=xs[[i]])$y
  }), 1:feature_n))
  pre=as.numeric(model$fitted.values)
  pre=pre-mean(pre)
  sum((rowSums(curve)-as.numeric(pre))^2)/1000
}

# 3. tuning function -----------------------------------------------------------
tune=function(combi){
  #data generation function is run from 4. Computing ccPDP.R
  data = create_xor_corr(n = 1000)# seed = 123
  X = data[,setdiff(names(data),"y")]
  task = as_task_regr(data,target="y")
  lrn = lrn("regr.nnet",size=size,decay=decay)
  model = lrn$train(task = task)$model

  #functions for case1
  #X = data_train[,setdiff(names(data_train),"quality")]
  #size=11
  #decay=0.0911569
  #task = as_task_regr(data_train,target="quality")
  #lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
  #model = lrn$train(task = task)$model

    if(combi=="eff3_gamma"|combi=="eff4_gamma"){
    #upper bound for gamma
    upper1=rep(20,feature_n)
  }else{
    #upper bound for kernel.width and gower.power
    upper1=rep(2,feature_n)
    }
  
  custom_loss=switch(combi,
                     "eff1_kw"=custom_loss1,
                     "eff1_gp"=custom_loss2,
                     "eff2_kw"=custom_loss3,
                     "eff2_gp"=custom_loss4,
                     "eff3_gamma"=custom_loss5,
                     "eff4_gamma"=custom_loss6)
  
  #loss minimization
  lgr::get_logger("bbotk")$set_threshold("warn")
  result=bb_optimize(custom_loss,lower=rep(0,feature_n),upper=upper1,max_evals=200,data =data_train,model=model$model)
  return(result$par)
}

# 4. parallelization -----------------------------------------------------------

set_seed_number=41
s=Sys.time()
cl<- makeCluster(core_number)      
registerDoParallel(cl)       
mydata<- foreach(
  j=c("eff1_kw","eff1_gp","eff2_kw","eff2_gp","eff3_gamma","eff4_gamma"),
  .combine=rbind,  
  .packages = c("iml","R.utils","mvtnorm","patchwork","data.table",
                    "StatMatch","dplyr","mlr3","mlr3verse","mlr3learners",
                    "nnet","mgcv","MASS","bbotk")
  ) %dopar% {
    set.seed(set_seed_number)
    tune(j)
    }
stopImplicitCluster()
stopCluster(cl)

mydata=as.data.frame(mydata)
rownames(mydata)=c("eff1_kw","eff1_gp","eff2_kw","eff2_gp","eff3_gamma","eff4_gamma")

write.csv(mydata,file=paste0(example_name,"parameters_fidelity.csv"))
e=Sys.time()
e-s


# evaluation --------------------------------------------------------------
set.seed(41)
#data generation function from 4. Computing ccPDP.R file
data = create_xor_corr(n = 1000)
n_feature=switch(example_name,
"example3"=3,
"example4"=5,
"case1"=11)
X = data[,setdiff(names(data),"y")]
task = as_task_regr(data,target="y")
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
model = lrn$train(task = task)
pred <- Predictor$new(model=model$model, data = data, y = "y")
#sum of xswcpdp for all variables of all points
#5 is taken here due to we use example 4
pars=read.csv("Tuned pars/pars_fidelity/example3parameters_fidelity.csv")
xswcpdp_sum = data.table::setDF(setNames(lapply(1:n_feature, function(i) {
  feature_i=paste0("x",i)
  effect2 (mod = model$model, data = data, feature = feature_i, target = "y",
           predict.fun = predict, h = 999, method = "xs-wcpdp",gower.power = 10,kernel.width=pars[3,i+1])$y
}), 1:n_feature))
for(i in 1:n_feature){
  xswcpdp_sum[,i]=xswcpdp_sum[,i]-mean(xswcpdp_sum[,i])
}

#sum of ale for all variables of all points
ale_sum = data.table::setDF(setNames(lapply(1:n_feature, function(i) {
  #h=999
  feature_i=paste0("x",i)
  q<-sort(data[[feature_i]])
  FeatureEffect$new(pred, feature = feature_i, method = "ale", grid.points =q)$results$.value
}), 1:n_feature))
#for(i in 1:n_feature){
#  ale_sum[,i]=ale_sum[,i]-mean(ale_sum[,i])
#}

#sum of ale for all variables of all points
pdp_sum = data.table::setDF(setNames(lapply(1:n_feature, function(i) {
  #h=999
  feature_i=paste0("x",i)
  q<-sort(data[[feature_i]])
  FeatureEffect$new(pred, feature = feature_i, method = "pdp", grid.points =q)$results$.value
}), 1:n_feature))
for(i in 1:n_feature){
  pdp_sum[,i]=pdp_sum[,i]-mean(pdp_sum[,i])
}

#R^2=0.7210609 for PDP decomposition
var(rowSums(pdp_sum))/var(model$model$fitted.values)

#R^2=0.8082447 for ALE decomposition
var(rowSums(ale_sum))/var(model$model$fitted.values)

#R^2=0.8237816 for ccPDP decomposition
var(rowSums(xswcpdp_sum))/var(model$model$fitted.values)

#plotting the dicrepancy
y=rowSums(xswcpdp_sum)
#ccpdp data frame for plotting
data_plot=data.frame(x=data[["x4"]],y=y)

#real prediction data frame for plotting
y_hat=as.numeric(model$model$fitted.values)
data_true=data.frame(x=data[["x4"]],y=y_hat)
ggplot() +
  geom_line(data = data_plot, aes(x = x, y = y, col = "'effect2-xswcpdp'"), lty = 1, lwd = 2) +
  geom_line(data = data_true, aes(x = x, y = y-mean(y), col = "true"), lty = 2, lwd = 2) +
  labs(x="x4")+
  guides()

#plot for 50 sample points
plot_id=sample(1:1000,50)
ggplot() +
  geom_line(data = data_plot[plot_id,], aes(x = x, y = y, col = "'own'"), lty = 1, lwd = 2) +
  geom_line(data = data_true[plot_id,], aes(x = x, y = y, col = "true"), lty = 2, lwd = 2) +
  guides()

#Loss value of xswcpdp
sum((data_plot$y-data_true$y)^2)/1000
#0.04982344 at kernel.width=0.6
