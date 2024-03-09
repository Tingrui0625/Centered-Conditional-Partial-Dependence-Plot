#loss function to calculate the discrepancy between ale and ccPDP
Loss<-function(curve,ale){
  a=curve$y[!is.na(curve$y)]
  # if(length(a)!=length(b)){
  #   b=b[-length(b)]
  # }
  b=ale$y[1:length(a)]
  sum((a-b)^2)/length(a)
}

#function to calcualte Mplot
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

#1. Path1-Restricting first(restricted ICE -> ccICE -> ccPDP)
effect1 = function(mod, data, feature, target, kernel.width, gamma, gower.power = 1, predict.fun = predict, h = 50,method){
  gower.power=as.numeric(gower.power)
  kernel.width=as.numeric(kernel.width)
  
  # interval bounds
  q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs
  
  # compute standard ICE values at quantile grid points
  ice_vals = data.table::setDF(setNames(lapply(q, function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), q))
  
  #calculate xs-proximity-based weights
  xs.weight = data.table::setDF(setNames(lapply(q, function(grid) {
    f.val = data[[feature]]
    xs.dist =exp(-0.5 * ((grid - f.val)^2) / (kernel.width^2))
  }), q))
  
  #calculate xc-proximity-based weights
  xc = data[, setdiff(colnames(data), c(feature, target)), drop = FALSE]
  xc.weight = 1 - as.matrix(cluster::daisy(xc, metric = "gower"))^gower.power
  
  #connecting the xs quantile values and corresponding xc.weight
  indicator2origin<-rep(0,h+1)
  #finding the id of the observation in original dataset having xs quantile values
  for(k in 1:(h+1)){
    indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  }
  xc.weight1=matrix(0,nrow=nrow(data),ncol=(h+1))
  #recording the found observation's xc.weight in xc.weight1
  for(i in 1:(h+1)){
    xc.weight1[,i]=xc.weight[,indicator2origin[i]]
  }
  xc.weight1=as.data.frame(xc.weight1)
  
  #creating index for additional conditioning methods with "-cc"
  if(method=="xs-wccpdp"|method=="xc-wccpdp"){
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    #computing the index indicating the observations within two quantile values next to each other
    #the first interval
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }

    #the second interval
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
    
    #mapping the index to a n*(h+1) matrix for use
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  
  #weighting the ICE curves with the weight matrix and the index mmatrix
  weighted_ice_vals=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
  switch(method,
         
         #xs-proximity-based
         "xs-wcpdp"={for(i in 1:(h+1)){weighted_ice_vals[,i]=ice_vals[,i] * xs.weight[,i]/sum(xs.weight[,i],na.rm=TRUE)}},
         
         #xs-proximity-based + additional conditioning
         "xs-wccpdp"={for(i in 1:(h+1)){
           if(sum(xs.weight[index[,i], i],na.rm=TRUE)==0){
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=0
           }else{
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=ice_vals[index[,i], i] * xs.weight[index[,i], i]/sum(xs.weight[index[,i], i],na.rm=TRUE)   
           }}},
         
         #xc-proximity-based
         "xc-wcpdp"={for(i in 1:(h+1)){weighted_ice_vals[,i]=ice_vals[,i] * xc.weight1[,i]/sum(xc.weight1[,i],na.rm=TRUE)}},
         
         #xc-proximity-based + additional conditioning
         "xc-wccpdp"={for(i in 1:(h+1)){
           #patch for extreme situations
           if(length(which(index[,i]==TRUE))==1){
             weighted_ice_vals[1,i]=ice_vals[index[,i], i]
           }else if(sum(xc.weight1[index[,i], i],na.rm=TRUE)==0){
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=0
           }else{
             #normal situations
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=ice_vals[index[,i], i] * xc.weight1[index[,i], i]/sum(xc.weight1[index[,i], i],na.rm=TRUE)
           }
         }},
         
         #centered PDP
         "cpdp"={weighted_ice_vals=ice_vals/nrow(ice_vals)})
  
  #centering the conditional ICE curves
  weighted_ice_vals=weighted_ice_vals-rowMeans(weighted_ice_vals,na.rm=TRUE)%*%matrix(1,nrow=1,ncol=h+1)
  
  #preparing the data frame for further plotting
  d = data.frame(x = q, y = colSums(weighted_ice_vals,na.rm=TRUE))
  return(d)
}

#2. centering first
effect2 = function(mod, data, feature, target, kernel.width, gower.power = 1, predict.fun = predict, h = 20, method = "pdp") {
  gower.power=as.numeric(gower.power)
  kernel.width=as.numeric(kernel.width)
  # Interval bounds
  q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs
  # compute standard ICE values at quantile grid points
  ice_vals = data.table::setDF(setNames(lapply(q, function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), q))
  
  # centering each ICE curve over xs by subtracting its mean
  ice_vals=ice_vals-rowMeans(ice_vals)%*%matrix(1,nrow=1,ncol=h+1) 
  
  # calculating xs-proximity-based weights
  xs.weight = data.table::setDF(setNames(lapply(q, function(grid) {
    f.val = data[[feature]]
    xs.dist = exp(-0.5 * ((grid - f.val)^2) / (kernel.width^2))
  }), q))
  
  #calculate xc-proximity-based weights
  xc = data[, setdiff(colnames(data), c(feature, target)), drop = FALSE]
  xc.weight =1-(as.matrix(cluster::daisy(xc, metric = "gower")))**as.numeric(gower.power)
  
  #connecting the xs quantile values and corresponding xc.weight
  indicator2origin<-rep(0,h+1)
  #finding the id of the observation in original dataset having xs quantile values
  for(k in 1:(h+1)){
    indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  }
  xc.weight1=matrix(0,nrow=nrow(data),ncol=(h+1))
  #recording the found observation's xc.weight in xc.weight1
  for(i in 1:(h+1)){
    xc.weight1[,i]=xc.weight[,indicator2origin[i]]
  }
  xc.weight1=as.data.frame(xc.weight1)
  
  #creating index for additional conditioning methods with "-cc"
  if(method=="xs-wccpdp"|method=="xc-wccpdp"){
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    #computing the index indicating the observations within two quantile values next to each other
    #the first interval
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    #the second interval
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
    
    #mapping the index to a n*(h+1) matrix for use
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }

  #weighting the centered ICE curves with the weight matrix and the index matrix
  y = vapply(1:ncol(ice_vals), function(i) {
    switch(method, 
           #xs-proximity-based
           "xs-wcpdp" = weighted.mean(ice_vals[, i], xs.weight[, i]), ##use all the weights for all points
           
           #xs-proximity-based + additional conditioning
           "xs-wccpdp" = {if(sum(xs.weight[index[,i], i],na.rm=TRUE)==0){
             0
           }else{weighted.mean(ice_vals[index[,i], i], xs.weight[index[,i], i],na.rm=TRUE)}},
           
           #xc-proximity-based
           "xc-wcpdp" =weighted.mean(ice_vals[, i], xc.weight1[,i]),
           
           #xc-proximity-based + additional conditioning
           "xc-wccpdp" = {if(sum(xs.weight[index[,i], i],na.rm=TRUE)==0){
             0
           }else{
             weighted.mean(ice_vals[index[,i], i], xc.weight1[index[,i],i],na.rm=TRUE)}},

           #centered pdp
           "cpdp" = mean(ice_vals[[i]]) # this would result in a centered PDP over xc without weights
    )
  }, FUN.VALUE = NA_real_)
  
  #preparing the data frame for further plotting
  d = data.frame(x = q, y = y)
  return(d)
}
#3. Restricting first and considering xs-xc-weights
effect3 = function(mod, data, feature, target, gamma, predict.fun = predict, h = 50,method) {
  gamma=as.numeric(gamma)
  
  # Interval bounds
  q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs
  
  # compute standard ICE values at quantile grid points
  ice_vals = data.table::setDF(setNames(lapply(q, function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), q))

  #connecting the xs quantile values and corresponding x
  indicator2origin<-rep(0,(h+1))
  for(k in 1:(h+1)){
    indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  }
  
  #calculating x-proximity-based weights
  X<-data[, setdiff(colnames(data), target), drop = FALSE]
  W_kernel_full=exp(-gamma*(as.matrix(cluster::daisy(X, metric = "gower")))^2)
  
  #creating index for additional conditioning methods with "-cc"
  if(method=="xs-xc-cpdp"){
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    #the first interval
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    #the second interval
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
    
    #mapping the index to a n*(h+1) matrix for use
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  
  #weighting the ICE curves with the weight matrix and the index matrix
  weighted_ice_vals=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
  switch(method,
         
         #x-proximity-based
         "xs-xc-pdp"={for(i in 1:(h+1)){
           weighted_ice_vals[,i]=ice_vals[,i] * W_kernel_full[,indicator2origin[i]]/sum(W_kernel_full[,indicator2origin[i]],na.rm=TRUE)}},
         
         #x-proximity-based + additional conditioning
         "xs-xc-cpdp"={for(i in 1:(h+1)){
           if(sum(W_kernel_full[index[,i],indicator2origin[i]],na.rm=TRUE)==0){
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=0
           }else{
             weighted_ice_vals[1:length(ice_vals[index[,i], i]),i]=ice_vals[index[,i], i] * W_kernel_full[index[,i],indicator2origin[i]]/sum(W_kernel_full[index[,i],indicator2origin[i]],na.rm=TRUE)}
         }},
         
         #centered PDP
         "cpdp"={weighted_ice_vals=ice_vals/nrow(ice_vals)})
  
  #centering the conditional ICE curves
  weighted_ice_vals=weighted_ice_vals-rowMeans(weighted_ice_vals,na.rm=TRUE)%*%matrix(1,nrow=1,ncol=h+1)
  
  #preparing the data frame for further plotting
  d = data.frame(x = q, y = colSums(weighted_ice_vals,na.rm=TRUE))
  return(d)
}
#4. Centering first and considering xs-xc-weights
effect4 = function(mod, data, feature, target, gamma, predict.fun = predict, h = 20, method = "pdp") {
  gamma=as.numeric(gamma)
  
  # Interval bounds
  q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs

  # compute standard ICE values at quantile grid points
  ice_vals = data.table::setDF(setNames(lapply(q, function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), q))
  
  # centering each ICE curve over xs by subtracting its mean
  ice_vals=ice_vals-rowMeans(ice_vals)%*%matrix(1,nrow=1,ncol=h+1)
  
  #connecting the xs quantile values and corresponding x
  indicator2origin<-rep(0,h+1)
  for(k in 1:(h+1)){
    indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  }
  
  #calculating x-proximity-based weights
  X<-data[, setdiff(colnames(data), target), drop = FALSE]
  W_kernel_full=exp(-gamma*(as.matrix(cluster::daisy(X, metric = "gower")))^2)
  #when gamma=0.2, almost all the weights are 0.9! -> close to centered PDP
  
  #creating index for additional conditioning methods with "-cc"
  if(method=="xs-xc-cpdp"){
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    #computing the index indicating the observations within two quantile values next to each other
    #the first interval
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    #the second interval
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
    
    #mapping the index to a n*(h+1) matrix for use
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  
  #weighting the centered ICE curves with the weight matrix and the index matrix
  y = vapply(1:ncol(ice_vals), function(i) {
    switch(method, 
           
           #x-proximity-based + additional conditioning 
           "xs-xc-cpdp"= {
             if(sum(W_kernel_full[index[,i], indicator2origin[i]],na.rm=TRUE)==0){
               0
             }else{weighted.mean(ice_vals[index[,i], i], W_kernel_full[index[,i], indicator2origin[i]],na.rm=TRUE)}
           },
           
           #x-proximity-based
           "xs-xc-pdp" = weighted.mean(ice_vals[, i], W_kernel_full[, indicator2origin[i]]),
           
           # centered pdp
           "cpdp" = mean(ice_vals[[i]]), 
    )
  }, FUN.VALUE = NA_real_)
  #preparing the data frame for further plotting
  d = data.frame(x = q, y = y)
  return(d)
}

#method with xs.weight multiplying xc.weight
effect5 = function(mod, data, feature, target, kernel.width, gower.power = 1, predict.fun = predict, h = 20, method = "pdp") {
  gower.power=as.numeric(gower.power)
  kernel.width=as.numeric(kernel.width)
  
  #Interval bounds
  q = quantile(data[[feature]], 0:h/h)#h quantile values of feature xs

  # compute standard ICE values at quantile grid points
  ice_vals = data.table::setDF(setNames(lapply(q, function(grid) {
    newdata = replace(data, list = which(colnames(data) == feature), values = grid)
    predict.fun(mod, newdata = newdata)
  }), q))
  
  # centering each ICE curve over xs by subtracting its mean
  ice_vals=ice_vals-rowMeans(ice_vals)%*%matrix(1,nrow=1,ncol=h+1)
  
  #connecting the xs quantile values and corresponding xc.weight
  indicator2origin<-rep(0,h+1)
  for(k in 1:(h+1)){
    indicator2origin[k]<-sort(abs(data[,feature]-q[k]),index.return=TRUE)$ix[1]
  }

  #xs-proximity-based
  xs.weight = data.table::setDF(setNames(lapply(q, function(grid) {
    f.val = data[[feature]]
    xs.dist = exp(-0.5 * ((grid - f.val)^2) / (kernel.width^2))
  }), q))
  
  #xc-proximity-based
  xc = data[, setdiff(colnames(data), c(feature, target)), drop = FALSE]#setdiff to find different elements
  xc.weight = 1 - as.matrix(cluster::daisy(xc, metric = "gower"))^gower.power#daisy() is to calculate distance for xc
  
  #creating index for additional conditioning methods with "-cc"
  if(method=="xs-wccpdp"|method=="xc-wccpdp"){
    ind_c=matrix(data=NA,nrow=nrow(data),ncol=h+1)
    
    #computing the index indicating the observations within two quantile values next to each other
    #the first interval
    m1=length(which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2))
    if(m1>0){
      ind_c[1:m1,1]=which(data[[feature]] >= q[1] & data[[feature]] < (q[2]+q[1])/2)
    }else if(q[1]==q[2]){
      ind_c[1:length(which(data[[feature]] == q[1])),1]=which(data[[feature]] == q[1])
    }
    
    #the second interval
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
    
    #mapping the index to a n*(h+1) matrix for use
    index=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
    for(i in 1:(h+1)){
      index[,i] = c(1:nrow(data))%in%ind_c[,i]
    }
  }
  
  #weighting the centered ICE curves with the weight matrix and the index matrix
  y = vapply(1:ncol(ice_vals), function(i) {
    switch(method, 
           #x-proximity-based
           "xs-xc-pdp" = weighted.mean(ice_vals[, i], xs.weight[, i]*xc.weight[,indicator2origin[i]]),
           
           #x-proximity-based + additional conditioning
           "xs-xc-cpdp" = {
             if(sum(xs.weight[index[,i], i]*xc.weight[index[,i],indicator2origin[i]],na.rm=TRUE)==0){
               0
             }else{
               weighted.mean(ice_vals[index[,i], i], xs.weight[index[,i], i]*xc.weight[index[,i],indicator2origin[i]],na.rm=TRUE)
             }},
           
           #centered pdp
           "cpdp" = mean(ice_vals[[i]])
    )
    
  }, FUN.VALUE = NA_real_)
  
  #preparing the data frame for further plotting
  d = data.frame(x = q, y = y)
  return(d)
}


