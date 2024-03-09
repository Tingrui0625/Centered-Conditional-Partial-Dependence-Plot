#parameter configuration
h=50
ite=112

# 1. Plotting ccPDP and recording the losses in for loop -----------------

for(effect_name in c("effect1","effect2","effect3","effect4","effect5")){
  for(example_name in c("example1","example2","example3","example4","example5","example6")){
  a=switch(example_name,
           "example1"=5,
           "example2"=3,
           "example3"=3,
           "example4"=5,
           "example5"=3,
           "example6"=3)
  LOSS=matrix(0,nrow=a,ncol=7)
  for(i in 1:a){
    feature=paste0("x",i)
    
    #loading mplots for further comparison
    mplot=read.csv(paste0("mplot result/",example_name,"_",feature,"mplot.csv"))
    
    #loading ccPDP curves
    eff=read.csv(paste0("Result_tuned_pars/Result_ale_seed41/",effect_name,"_",example_name,"_",feature,".csv"))
    
    #running the plots
    SavePlots(effect_name=effect_name,example_name=example_name,feature=feature,eff=eff,mplot=mplot,centering=TRUE,
              plotting=TRUE,Path="Result_tuned_pars/Result_ale_seed41/Plots/")
    
    #recording the losses in matrix LOSS
    LOSS[i,]=SavePlots(effect_name=effect_name,example_name=example_name,feature=feature,eff=eff,mplot=mplot,centering=TRUE,
                       plotting=FALSE,Path="Result_tuned_pars/Result_ale_seed41/Plots/")
  }
  LOSS=rbind(LOSS,colMeans(LOSS))
  LOSS=as.data.frame(LOSS)
  colnames(LOSS)=c("cpdp","xswcpdp","xswccpdp","xcwcpdp","xcwccpdp","xsxcpdp","cxsxspdp")
  
  #saving the loss data frame
  write.csv(LOSS,file=paste0("Result_tuned_pars/Result_ale_seed41/LOSS/",effect_name,"_",example_name,".csv"))
 }
}

#plotting and loss functions
SavePlots=function(effect_name,example_name,feature,eff,mplot,centering,plotting,h=50,ite=112,Path="Result/Plots/"){
  Path=Path
  
  #computing average mplot
  if(plotting==TRUE){
   mplot=mplot
  colnames(mplot)=c("id","x","mplot","n")
  grid=matrix(mplot[,"x"],ncol=ite,nrow=5000)
  mv=matrix(mplot[,"mplot"],ncol=ite,nrow=5000)
  mv1=matrix(0,ncol=ite,nrow=5000)
  for(i in 1:ite){
    index=sort(grid[,i],index.return=TRUE)$ix
    mv1[,i]=mv[index,i]
    grid[,i]=sort(grid[,i])
  }
  mpl=data.frame(x=rowMeans(grid),y=rowMeans(mv1)) 
  }
  
  #computing average quantiles
  quantile=matrix(eff[,"grid"],ncol=ite,nrow=h+1)
  
  #computing average ale curves
  ale=matrix(eff[,"ale"],ncol=ite,nrow=h+1)
  ale_plot=data.frame(x=rowMeans(quantile),y=rowMeans(ale))
  
  #computing average pdp curves
  pdp=matrix(eff[,"pdp"],ncol=ite,nrow=h+1)
  pdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(pdp))
  
  # -------------------------------------------------------------------------
  
  #preparing cpdp plot
  cpdp=matrix(eff[,"cpdp"],ncol=ite,nrow=h+1)
  if(centering==TRUE){
    cpdp=cpdp-matrix(1,nrow=h+1,ncol=1)%*%colMeans(cpdp)
  }
  cpdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(cpdp))
  
  #creating confidence intervals
  if(plotting==TRUE){
    lower=rep(0,(h+1))
    upper=rep(0,(h+1))
    for(i in 1:(h+1)){
      lower[i]=mean(cpdp[i,])-1.96*sd(cpdp[i,])
      upper[i]=mean(cpdp[i,])+1.96*sd(cpdp[i,])
    }
    cpdp_plot$lower=lower
    cpdp_plot$upper=upper
    
    #plotting and saving
    p1=ggplot(data = cpdp_plot,aes(x, y))+
      geom_ribbon(aes(x=x, ymax=upper, ymin=lower), fill="grey", alpha=.5) +
      geom_line(aes(y = upper), colour = 'grey',lty=3,alpha=.5) +
      geom_line(aes(y = lower), colour = 'grey',lty=3,alpha=.5)+
      geom_line(data = ale_plot, aes(x = x, y = y, color = "ale"), lty = 1,lwd = 2) + 
      geom_line(data = pdp_plot, aes(x = x, y = y, color = "pdp"), lty = 1,lwd = 1)+
      ##geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
      geom_line(data = cpdp_plot,aes(x = x, y = y, color = "cpdp"), lty = 1, lwd = 1)+
      #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
      scale_color_manual(name = "Effects", 
                         values = c("ale" = "lightgreen", 
                                    "pdp" = "orange1", 
                                    "cpdp" = "grey27",
                                    "xswcpdp" = "maroon3",
                                    "xswccpdp"="palevioletred",
                                    "xcwcpdp"="dodgerblue1",
                                    "xcwccpdp"="cornflowerblue",
                                    "mplot"="mistyrose3"
                         ))+
      labs(x=feature)
    #p1=p1+ylim(-0.2,0.5)+xlim(-1,1)
    ggsave(p1,
           file=paste0(Path,effect_name,"_",example_name,"_",feature,"_cpdp.png"),
           scale = 2.5,
           width = 7.75,             
           height = 5.88,
           units="cm",
           dpi = 500 )
  }else{
    Loss_empirical=rep(0,7)
    Loss_empirical[1]=Loss(cpdp_plot,ale_plot)
  }
  
  if(effect_name=="effect1"|effect_name=="effect2"){
    #preparing xswcpdp curve for plotting
    xswcpdp=matrix(eff[,"xswcpdp"],ncol=ite,nrow=h+1)
    if(centering==TRUE){
      xswcpdp=xswcpdp-matrix(1,nrow=h+1,ncol=1)%*%colMeans(xswcpdp)
    }
    xswcpdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(xswcpdp))

    #creating confidence intervals
    if(plotting==TRUE){
      lower=rep(0,(h+1))
      upper=rep(0,(h+1))
      for(i in 1:(h+1)){
        lower[i]=mean(xswcpdp[i,])-1.96*sd(xswcpdp[i,])
        upper[i]=mean(xswcpdp[i,])+1.96*sd(xswcpdp[i,])
      }
      xswcpdp_plot$lower=lower
      xswcpdp_plot$upper=upper
      
      #plotting and saving
      p1=ggplot(data = xswcpdp_plot,aes(x, y))+
        geom_ribbon(aes(x=x, ymax=upper, ymin=lower), fill="grey", alpha=.5) +
        geom_line(aes(y = upper), colour = 'grey',lty=3,alpha=.5) +
        geom_line(aes(y = lower), colour = 'grey',lty=3,alpha=.5)+
        geom_line(data = ale_plot, aes(x = x, y = y, color = "ale"), lty = 1,lwd = 2) + 
        geom_line(data = pdp_plot, aes(x = x, y = y, color = "pdp"), lty = 1,lwd = 1) +
        ##geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = xswcpdp_plot,aes(x = x, y = y, color = "xswcpdp"), lty = 1, lwd = 1)+
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xswcpdp" = "maroon3",
                                      "xswccpdp"="palevioletred",
                                      "xcwcpdp"="dodgerblue1",
                                      "xcwccpdp"="cornflowerblue",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_xswcpdp.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
    }else{
      Loss_empirical[2]=Loss(xswcpdp_plot,ale_plot)
    }
    
    # -------------------------------------------------------------------------
    
    
    #preparing xswccpdp curve for plotting
    xswccpdp=matrix(eff[,"xswccpdp"],ncol=ite,nrow=h+1)

    if(centering==TRUE){
      xswccpdp=xswccpdp-matrix(1,nrow=h+1,ncol=1)%*%colMeans(xswccpdp)
    }
    xswccpdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(xswccpdp))

    #creating confidence intervals
    if(plotting==TRUE){
      lower=rep(0,(h+1))
      upper=rep(0,(h+1))
      for(i in 1:(h+1)){
        lower[i]=mean(xswccpdp[i,])-1.96*sd(xswccpdp[i,])
        upper[i]=mean(xswccpdp[i,])+1.96*sd(xswccpdp[i,])
      }
      xswccpdp_plot$lower=lower
      xswccpdp_plot$upper=upper
      
      #plotting and saving
      p1=ggplot(data = xswccpdp_plot,aes(x, y))+
        geom_ribbon(aes(x=x, ymax=upper, ymin=lower), fill="grey", alpha=.5) +
        geom_line(aes(y = upper), colour = 'grey',lty=3,alpha=.5) +
        geom_line(aes(y = lower), colour = 'grey',lty=3,alpha=.5)+
        geom_line(data = ale_plot, aes(x = x, y = y, color = "ale"), lty = 1,lwd = 2) + 
        geom_line(data = pdp_plot, aes(x = x, y = y, color = "pdp"), lty = 1,lwd = 1)+ 
        ##geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = xswccpdp_plot,aes(x = x, y = y, color = "xswccpdp"), lty = 1, lwd = 1)+
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xswcpdp" = "maroon3",
                                      "xswccpdp"="palevioletred",
                                      "xcwcpdp"="dodgerblue1",
                                      "xcwccpdp"="cornflowerblue",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_xswccpdp.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
    }else{
      Loss_empirical[3]=Loss(xswccpdp_plot,ale_plot)
    }
 
    # -------------------------------------------------------------------------
    
    
    #preparing xcwcpdp curve for plotting
    xcwcpdp=matrix(eff[,"xcwcpdp"],ncol=ite,nrow=h+1)

    if(centering==TRUE){
      xcwcpdp=xcwcpdp-matrix(1,nrow=h+1,ncol=1)%*%colMeans(xcwcpdp)
    }
    xcwcpdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(xcwcpdp))

    #confdence interval
    if(plotting==TRUE){
      
      lower=rep(0,(h+1))
      upper=rep(0,(h+1))
      for(i in 1:(h+1)){
        lower[i]=mean(xcwcpdp[i,])-1.96*sd(xcwcpdp[i,])
        upper[i]=mean(xcwcpdp[i,])+1.96*sd(xcwcpdp[i,])
      }
      xcwcpdp_plot$lower=lower
      xcwcpdp_plot$upper=upper
      
      #plotting and saving
      p1=ggplot(data = xcwcpdp_plot,aes(x, y))+
        geom_ribbon(aes(x=x, ymax=upper, ymin=lower), fill="grey", alpha=.5) +
        geom_line(aes(y = upper), colour = 'grey',lty=3,alpha=.5) +
        geom_line(aes(y = lower), colour = 'grey',lty=3,alpha=.5)+
        geom_line(data = ale_plot, aes(x = x, y = y, color = "ale"), lty = 1,lwd = 2) + 
        geom_line(data = pdp_plot, aes(x = x, y = y, color = "pdp"), lty = 1,lwd = 1)+ 
        ##geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = xcwcpdp_plot,aes(x = x, y = y, color = "xcwcpdp"), lty = 1, lwd = 1)+
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xswcpdp" = "maroon3",
                                      "xswccpdp"="palevioletred",
                                      "xcwcpdp"="dodgerblue1",
                                      "xcwccpdp"="cornflowerblue",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_xcwcpdp.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
    }else{
      Loss_empirical[4]=Loss(xcwcpdp_plot,ale_plot)
    }
    
    # -------------------------------------------------------------------------
    
    
    #preparing xcwccpdp curve for plotting
    xcwccpdp=matrix(eff[,"xcwccpdp"],ncol=ite,nrow=h+1)

    if(centering==TRUE){
      xcwccpdp=xcwccpdp-matrix(1,nrow=h+1,ncol=1)%*%colMeans(xcwccpdp)
    }
    xcwccpdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(xcwccpdp))
    
    #confidence interval
    if(plotting==TRUE){
      lower=rep(0,(h+1))
      upper=rep(0,(h+1))
      for(i in 1:(h+1)){
        lower[i]=mean(xcwccpdp[i,])-1.96*sd(xcwccpdp[i,])
        upper[i]=mean(xcwccpdp[i,])+1.96*sd(xcwccpdp[i,])
      }
      xcwccpdp_plot$lower=lower
      xcwccpdp_plot$upper=upper
      
      p1=ggplot(data = xcwccpdp_plot,aes(x, y))+
        geom_ribbon(aes(x=x, ymax=upper, ymin=lower), fill="grey", alpha=.5) +
        geom_line(aes(y = upper), colour = 'grey',lty=3,alpha=.5) +
        geom_line(aes(y = lower), colour = 'grey',lty=3,alpha=.5)+
        geom_line(data = ale_plot, aes(x = x, y = y, col = "ale"), lty = 1,lwd = 2) + 
        geom_line(data = pdp_plot, aes(x = x, y = y, col = "pdp"), lty = 1,lwd = 1)+ 
        #geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = xcwccpdp_plot,aes(x = x, y = y, col = "xcwccpdp"), lty = 1, lwd = 1)+
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xswcpdp" = "maroon3",
                                      "xswccpdp"="palevioletred",
                                      "xcwcpdp"="dodgerblue1",
                                      "xcwccpdp"="cornflowerblue",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_xcwccpdp.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
      
      # all ---------------------------------------------------------------------
      p1=ggplot() +
        geom_line(data = ale_plot, aes(x = x, y = y, color="ale"), lty = 1,lwd = 2) +
        geom_line(data = pdp_plot, aes(x = x, y = y, color="pdp"), lty = 1, lwd = 1) +
        geom_line(data = cpdp_plot, aes(x = x, y = y, color="cpdp"), lty = 1, lwd = 1,) +
        #geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = xswccpdp_plot, aes(x = x, y = y, color="xswccpdp"), lty = 1, lwd = 1) +
        geom_line(data = xcwcpdp_plot, aes(x = x, y = y,color="xcwcpdp"), lty = 1, lwd = 1) +
        geom_line(data = xcwccpdp_plot, aes(x = x, y = y,color="xcwccpdp"), lty = 1, lwd = 1) +
        geom_line(data = xswcpdp_plot, aes(x = x, y = y, color="xswcpdp"), lty = 1, lwd = 1) +
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xswcpdp" = "maroon3",
                                      "xswccpdp"="palevioletred",
                                      "xcwcpdp"="dodgerblue1",
                                      "xcwccpdp"="cornflowerblue",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_all.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
      
      #plotting shaded ccpdps
      cpdp_plot$id=rep(1,h+1)
      xswcpdp_plot$id=rep(2,h+1)
      xswccpdp_plot$id=rep(3,h+1)
      xcwcpdp_plot$id=rep(4,h+1)
      xcwccpdp_plot$id=rep(5,h+1)
      
      data=rbind(cpdp_plot,xswcpdp_plot,xswccpdp_plot,xcwcpdp_plot,xcwccpdp_plot)
      
      p1=ggplot(data,aes(x=x, y=y,color=id))+
        geom_ribbon(data=cpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="grey27", alpha=.1)+
        geom_ribbon(data=xswcpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="maroon3", alpha=.1)+
        geom_ribbon(data=xswccpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="palevioletred", alpha=.1)+
        geom_ribbon(data=xcwcpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="dodgerblue1", alpha=.1)+
        geom_ribbon(data=xcwccpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="cornflowerblue", alpha=.1)+
        geom_line(data = ale_plot, aes(x = x, y = y, color="ale"), lty = 1,lwd = 2) +
        geom_line(data = pdp_plot, aes(x = x, y = y, color="pdp"), lty = 1, lwd = 1) +
        #geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = cpdp_plot, aes(x = x, y = y, color="cpdp"), lty = 1, lwd = 1,) +
        geom_line(data = xswccpdp_plot, aes(x = x, y = y, color="xswccpdp"), lty = 1, lwd = 1) +
        geom_line(data = xcwcpdp_plot, aes(x = x, y = y,color="xcwcpdp"), lty = 1, lwd = 1) +
        geom_line(data = xcwccpdp_plot, aes(x = x, y = y,color="xcwccpdp"), lty = 1, lwd = 1) +
        geom_line(data = xswcpdp_plot, aes(x = x, y = y, color="xswcpdp"), lty = 1, lwd = 1) +
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xswcpdp" = "maroon3",
                                      "xswccpdp"="palevioletred",
                                      "xcwcpdp"="dodgerblue1",
                                      "xcwccpdp"="cornflowerblue",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_all_shaded.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
    }else{
      Loss_empirical[5]=Loss(xcwccpdp_plot,ale_plot)
    }
    
  }
  
  
  # Effect34 --------------------------------------------------------------------------------------------------------------------
  
  
  if(effect_name=="effect3"|effect_name=="effect4"|effect_name=="effect5"){
    #preparing xsxcpdp curve for plotting
    xsxcpdp=matrix(eff[,"xsxcpdp"],ncol=ite,nrow=h+1)
    #centering
    if(centering==TRUE){
      xsxcpdp=xsxcpdp-matrix(1,nrow=h+1,ncol=1)%*%colMeans(xsxcpdp)
    }
    xsxcpdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(xsxcpdp))
    
    if(plotting==TRUE){
      
      lower=rep(0,(h+1))
      upper=rep(0,(h+1))
      for(i in 1:(h+1)){
        lower[i]=mean(xsxcpdp[i,])-1.96*sd(xsxcpdp[i,])
        upper[i]=mean(xsxcpdp[i,])+1.96*sd(xsxcpdp[i,])
      }
      xsxcpdp_plot$lower=lower
      xsxcpdp_plot$upper=upper
      
      p1=ggplot(data = xsxcpdp_plot,aes(x, y))+
        geom_ribbon(aes(x=x, ymax=upper, ymin=lower), fill="grey", alpha=.5) +
        geom_line(aes(y = upper), colour = 'grey',lty=3,alpha=.5) +
        geom_line(aes(y = lower), colour = 'grey',lty=3,alpha=.5)+
        geom_line(data = ale_plot, aes(x = x, y = y, color = "ale"), lty = 1,lwd = 2) + 
        geom_line(data = pdp_plot, aes(x = x, y = y, color = "pdp"), lty = 1,lwd = 1) +
        #geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = xsxcpdp_plot,aes(x = x, y = y, color = "xsxcpdp"), lty = 1, lwd = 1)+
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xsxcpdp"="blueviolet",
                                      "cxsxcpdp"="mediumpurple4",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_xsxcpdp.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
    }else{
      Loss_empirical[6]=Loss(xsxcpdp_plot,ale_plot)
    }
    
    
    # -------------------------------------------------------------------------
    
    #preparing cxsxcpdp curve for plotting
    cxsxcpdp=matrix(eff[,"cxsxcpdp"],ncol=ite,nrow=h+1)

    if(centering==TRUE){
      cxsxcpdp=cxsxcpdp-matrix(1,nrow=h+1,ncol=1)%*%colMeans(cxsxcpdp)
    }
    cxsxcpdp_plot=data.frame(x=rowMeans(quantile),y=rowMeans(cxsxcpdp))

    if(plotting==TRUE){
      lower=rep(0,(h+1))
      upper=rep(0,(h+1))
      for(i in 1:(h+1)){
        lower[i]=mean(cxsxcpdp[i,])-1.96*sd(cxsxcpdp[i,])
        upper[i]=mean(cxsxcpdp[i,])+1.96*sd(cxsxcpdp[i,])
      }
      cxsxcpdp_plot$lower=lower
      cxsxcpdp_plot$upper=upper
      
      p1=ggplot(data = cxsxcpdp_plot,aes(x, y))+
        geom_ribbon(aes(x=x, ymax=upper, ymin=lower), fill="grey", alpha=.5) +
        geom_line(aes(y = upper), colour = 'grey',lty=3,alpha=.5) +
        geom_line(aes(y = lower), colour = 'grey',lty=3,alpha=.5)+
        geom_line(data = ale_plot, aes(x = x, y = y, color = "ale"), lty = 1,lwd = 2) + 
        geom_line(data = pdp_plot, aes(x = x, y = y, color = "pdp"), lty = 1,lwd = 1) +
        #geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = cxsxcpdp_plot,aes(x = x, y = y, color = "cxsxcpdp"), lty = 1, lwd = 1)+
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xsxcpdp"="blueviolet",
                                      "cxsxcpdp"="mediumpurple4",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_cxsxcpdp.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 ) 
      
      p1=ggplot() +
        geom_line(data = ale_plot, aes(x = x, y = y, color="ale"), lty = 1,lwd = 2) +
        geom_line(data = pdp_plot, aes(x = x, y = y, color="pdp"), lty = 1, lwd = 1) +
        #geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = cpdp_plot, aes(x = x, y = y, color="cpdp"), lty = 1, lwd = 1,) +
        geom_line(data = xsxcpdp_plot,aes(x = x, y = y, color = "xsxcpdp"), lty = 1, lwd = 1)+
        geom_line(data = cxsxcpdp_plot,aes(x = x, y = y, color = "cxsxcpdp"), lty = 1, lwd = 1)+
        #stat_function(data = NULL, fun = function(x) 0, mapping = aes(col = "function f(x5) = 0"))+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xsxcpdp"="blueviolet",
                                      "cxsxcpdp"="mediumpurple4",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_all.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
      
      #plotting shaded ccPDPs
      cpdp_plot$id=rep(1,h+1)
      xsxcpdp_plot$id=rep(2,h+1)
      cxsxcpdp_plot$id=rep(3,h+1)
      data=rbind(cpdp_plot,xsxcpdp_plot,cxsxcpdp_plot)
      p1=ggplot(data,aes(x=x, y=y,color=id))+
        geom_ribbon(data=cpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="grey27", alpha=.1)+
        geom_ribbon(data=xsxcpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="blueviolet", alpha=.1)+
        geom_ribbon(data=cxsxcpdp_plot,aes(x = x,ymin = lower, ymax = upper), inherit.aes = FALSE,fill="mediumpurple4", alpha=.1)+
        geom_line(data = ale_plot, aes(x = x, y = y, color="ale"), lty = 1,lwd = 2) +
        geom_line(data = pdp_plot, aes(x = x, y = y, color="pdp"), lty = 1, lwd = 1) +
        #geom_line(data = mpl, aes(x, y, color = "mplot"),lty = 1, lwd = 1)+
        geom_line(data = cpdp_plot, aes(x = x, y = y, color="cpdp"), lty = 1, lwd = 1,) +
        geom_line(data = xsxcpdp_plot,aes(x = x, y = y, color = "xsxcpdp"), lty = 1, lwd = 1)+
        geom_line(data = cxsxcpdp_plot,aes(x = x, y = y, color = "cxsxcpdp"), lty = 1, lwd = 1)+
        scale_color_manual(name = "Effects", 
                           values = c("ale" = "lightgreen", 
                                      "pdp" = "orange1", 
                                      "cpdp" = "grey27",
                                      "xsxcpdp"="blueviolet",
                                      "cxsxcpdp"="mediumpurple4",
                                      "mplot"="mistyrose3"
                           ))+
        labs(x=feature)
      ggsave(p1,
             file=paste0(Path,effect_name,"_",example_name,"_",feature,"_all_shaded.png"),
             scale = 2.5,
             width = 7.75,             
             height = 5.88,
             units="cm",
             dpi = 500 )
    }else{
      Loss_empirical[7]=Loss(cxsxcpdp_plot,ale_plot)
    }
    }
    if(plotting==FALSE){
      return(Loss_empirical)
    }
  
}


# 2. Varying kernel.width plots-------------------------------------
example_name="example4"
feature="x4"
h = 20
set.seed(41)
data = create_xor_corr(n = 1000)#after running Console.R
X = data[,setdiff(names(data),"y")]
task = as_task_regr(data,target="y")
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
model = lrn$train(task = task)

xswcpdp1 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 2,kernel.width=0.6)
xswcpdp2 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 0.49,kernel.width=0.2)
xswcpdp3 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 0.49,kernel.width=0.4)
xswcpdp4 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 0.49,kernel.width=0.8)
xswcpdp5 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 0.49,kernel.width=1)
xswcpdp6 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 0.49,kernel.width=10)
xswcpdp1$y=xswcpdp1$y-mean(xswcpdp1$y)
xswcpdp2$y=xswcpdp2$y-mean(xswcpdp2$y)
xswcpdp3$y=xswcpdp3$y-mean(xswcpdp3$y)
xswcpdp4$y=xswcpdp4$y-mean(xswcpdp4$y)
xswcpdp5$y=xswcpdp5$y-mean(xswcpdp5$y)
xswcpdp6$y=xswcpdp6$y-mean(xswcpdp6$y)

cpdp = effect2(mod = model$model, data = data, feature = feature, target = "y",
               gower.power = 2, predict.fun = predict, h = h, method = "cpdp",kernel.width=0.6)
pred <- Predictor$new(model=model$model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
ice<-FeatureEffect$new(pred, feature = feature, method = "ice", grid.points =q)


mpl = mplot(data, feature, target = "y", eps = 0.5)

xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
cpdp$y=cpdp$y-mean(cpdp$y)
pdp$results$.value=pdp$results$.value-mean(pdp$results$.value)

ggplot() +
  #stat_function(data = NULL, fun = function(x) (x), mapping = aes(col = "function f(x)")) +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = xswcpdp2, aes(x = x, y = y, col = "xswcpdp kw=0.2"), lty = 3, lwd = 1,alpha=0.5) +
  geom_line(data = xswcpdp1, aes(x = x, y = y, col = "xswcpdp kw=0.6"), lty = 1, lwd = 2) +
  geom_line(data = xswcpdp3, aes(x = x, y = y, col = "xswcpdp kw=0.4"), lty = 3, lwd = 1,) +
  geom_line(data = xswcpdp4, aes(x = x, y = y, col = "xswcpdp kw=0.8"), lty = 3, lwd = 1) +
  geom_line(data = xswcpdp5, aes(x = x, y = y, col = "xswcpdp kw=1"), lty = 3, lwd = 1) +
  geom_line(data = xswcpdp6, aes(x = x, y = y, col = "xswcpdp kw=10"), lty = 3, lwd = 2) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  scale_color_manual(name = "Effects",
                     values = c("ale" = "lightgreen",
                                "pdp" = "orange1",
                                "cpdp" = "grey27",
                                "xswcpdp kw=0.2" = "maroon",
                                "xswcpdp kw=0.4" = "maroon1",
                                "xswcpdp kw=0.6" = "maroon2",
                                "xswcpdp kw=0.8" = "maroon3",
                                "xswcpdp kw=1" = "maroon4",
                                "xswcpdp kw=10" = "maroon4"
                                
                                
                     ))+
  guides()

# 3. Varying gower.power plots--------------------------------------
xcwcpdp1 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xc-wcpdp",gower.power = 2,kernel.width=0.6)
xcwcpdp2 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xc-wcpdp",gower.power = 0.00093,kernel.width=0.6)
xcwcpdp3 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xc-wcpdp",gower.power = 0.01,kernel.width=0.6)
xcwcpdp4 <-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                    predict.fun = predict, h = h, method = "xc-wcpdp",gower.power = 1,kernel.width=0.6)
xcwcpdp1$y=xcwcpdp1$y-mean(xcwcpdp1$y)
xcwcpdp2$y=xcwcpdp2$y-mean(xcwcpdp2$y)
xcwcpdp3$y=xcwcpdp3$y-mean(xcwcpdp3$y)
xcwcpdp4$y=xcwcpdp4$y-mean(xcwcpdp4$y)

cpdp = effect2(mod = model$model, data = data, feature = feature, target = "y",
               gower.power = 2, predict.fun = predict, h = h, method = "cpdp",kernel.width=0.6)
pred <- Predictor$new(model=model$model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)
ice<-FeatureEffect$new(pred, feature = feature, method = "ice", grid.points =q)

mpl = mplot(data, feature, target = "y", eps = 0.5)

cpdp$y=cpdp$y-mean(cpdp$y)
pdp$results$.value=pdp$results$.value-mean(pdp$results$.value)

ggplot() +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = xcwcpdp1, aes(x = x, y = y, col = "xcwcpdp GP=2"), lty = 3, lwd = 1) +
  geom_line(data = xcwcpdp3, aes(x = x, y = y, col = "xcwcpdp GP=0.01"), lty = 3, lwd = 1) +
  geom_line(data = xcwcpdp4, aes(x = x, y = y, col = "xcwcpdp GP=1"), lty = 3, lwd = 1) +
  geom_line(data = xcwcpdp2, aes(x = x, y = y, col = "xcwcpdp GP=0.00093"), lty = 3, lwd = 2) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 2) +
  scale_color_manual(name = "Effects",
                     values = c("ale" = "lightgreen",
                                "pdp" = "orange1",
                                "cpdp" = "grey27",
                                "xcwcpdp GP=0.00093" = "maroon1",
                                "xcwcpdp GP=0.01" = "maroon2",
                                "xcwcpdp GP=1" = "maroon3",
                                "xcwcpdp GP=2" = "maroon4"
                     ))+
  guides()
# 4. GAM v.s. NN under different Corr ----------------------------------------
h = 20
set.seed(41)
data = create_xor_corr(n = 1000)
X = data[,setdiff(names(data),"y")]
#GAM model
model=switch(example_name,
             "example1"=gam(y ~ x1+I(x2^2)+I(x3^3)+I(x2*x4), data = data),
             "example2"=gam(y ~ I(x1^2)+x2,data = data),
             "example3"=gam(y ~ x1+x2+I(x1*x2),data=data),
             "example4"=gam(y ~ x1+I(x2^2)+I(x3^3)+x3+I(x2*x4),data=data),
             "example5"=gam(y ~ I(x1^2)+x2,data=data),
             "example6"=gam(y ~ I(x1^2)+x2,data=data)
)

xswcpdp <-effect2 (mod = model, data = data, feature = feature, target = "y",
                   predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 2,kernel.width=0.55)
cpdp = effect2(mod = model, data = data, feature = feature, target = "y",
               gower.power = 2, predict.fun = predict, h = h, method = "cpdp",kernel.width=0.55)

xswcpdp$y=xswcpdp$y-mean(xswcpdp$y)
cpdp$y=cpdp$y-mean(cpdp$y)

pred <- Predictor$new(model=model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp", grid.points =q)
ale <- FeatureEffect$new(pred, feature = feature, method = "ale", grid.points =q)

#FeatureEffect$predict(data, extrapolate = TRUE)


#nnet
task = as_task_regr(data,target="y")
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
model = lrn$train(task = task)

xswcpdp_nnet<-effect2 (mod = model$model, data = data, feature = feature, target = "y",
                       predict.fun = predict, h = h, method = "xs-wcpdp",gower.power = 2,kernel.width=0.55)
cpdp_nnet= effect2(mod = model$model, data = data, feature = feature, target = "y",
                   gower.power = 2, predict.fun = predict, h = h, method = "cpdp",kernel.width=0.55)
pred_nnet<- Predictor$new(model=model$model, data = data, y = "y")
q<-quantile(data[[feature]], 0:h/h)
pdp_nnet  <- FeatureEffect$new(pred_nnet, feature = feature, method = "pdp", grid.points =q)
ale_nnet  <- FeatureEffect$new(pred_nnet, feature = feature, method = "ale", grid.points =q)
ice_nnet <-FeatureEffect$new(pred_nnet, feature = feature, method = "ice", grid.points =q)

xswcpdp_nnet$y=xswcpdp_nnet$y-mean(xswcpdp_nnet$y)
cpdp_nnet$y=cpdp_nnet$y-mean(cpdp_nnet$y)

#shifting constant
min(ale$results$.value)
#0.4*(x^2)-0.15 for high Corr
#x^2-0.015 for middle Corr
#0.6*(x^2)-0.0078 for low Corr
ggplot() +
  geom_line(data = pdp$results, aes_string(x = feature, y = ".value", col = "'pdp'"), lty = 1, lwd = 1)+
  geom_line(data = pdp_nnet$results, aes_string(x = feature, y = ".value", col = "'pdp_nnet'"), lty = 3, lwd = 1) +
  geom_line(data = cpdp, aes(x = x, y = y, col = "cpdp"), lty = 1, lwd = 1) +
  geom_line(data = cpdp_nnet, aes(x = x, y = y, col = "cpdp_nnet"), lty = 3, lwd = 1) +
  geom_line(data = xswcpdp, aes(x = x, y = y, col = "xswcpdp"), lty = 1, lwd = 2) +
  geom_line(data = xswcpdp_nnet, aes(x = x, y = y, col = "xswcpdp_nnet"), lty = 3, lwd = 1) +
  geom_line(data = ale$results, aes_string(x = feature, y = ".value", col = "'ale'"),lty = 1, lwd = 1) +
  geom_line(data = ale_nnet$results, aes_string(x = feature, y = ".value", col = "'ale_nnet'"),lty = 3, lwd = 1) +
  #stat_function(data = NULL, fun = function(x) 0.6*(x^2)-0.0078, mapping = aes(col = "function f(x)")) +
  scale_color_manual(name = "Effects", 
                     values = c("ale" = "lightgreen", 
                                "pdp" = "orange1", 
                                "cpdp" = "grey27",
                                "xswcpdp" = "maroon3",
                                "ale_nnet" = "lightgreen", 
                                "pdp_nnet" = "orange1", 
                                "cpdp_nnet" = "grey27",
                                "xswcpdp_nnet" = "maroon3",
                                "function f(x)"="blueviolet"
                     ))+
  guides()

# 5. Contour for x2 and x4 ---------------------------------------------------
n=1000
create_xor_corr = function(n){
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
  y=x1+0.5*(3*(x2^2)-1)+0.5*(4*(x3^3)-3*x3)+0.8*x2*x4+rnorm(n, sd = 0.1)
  data.frame(x1,x2,x3,x4,x5,y)
}
h = 20
set.seed(41)
data = create_xor_corr(n = 1000)
input_1 <- data$x2
input_2 <- data$x4

z=data[,c("x2","x4","y")]
ggplot(z, aes(x = x2, y = x4)) +
  geom_density_2d_filled()

# install.packages("plotly")
# library(plotly)
# z=as.matrix(z)
# plot_ly(x = x2, y = x4, z = z) %>% add_surface()

# 6. pdp+ice -----------------------------------------------------------------
#plotting for centeredICE + ccPDP

ice_vals = data.table::setDF(setNames(lapply(q, function(grid) {
  newdata = replace(data, list = which(colnames(data) == feature), values = grid)
  predict(model$model, newdata = newdata)
}), q))
ice_vals=ice_vals-rowMeans(ice_vals)%*%matrix(1,nrow=1,ncol=h+1)
ice_centered=as.vector(t(ice_vals))
pdp <- FeatureEffect$new(pred, feature = feature, method = "pdp+ice", grid.points =q)
pdp$results$.value[1:(h+1)]=xswcpdp$y
pdp$results$.value[(h+2):length(pdp$results$.value)]=ice_centered
pdp$plot()

#plotting for conditional ICE+ccPDP
ice_vals = data.table::setDF(setNames(lapply(q, function(grid) {
  newdata = replace(data, list = which(colnames(data) == feature), values = grid)
  predict(model$model, newdata = newdata)
}), q))
kernel.width=0.6
xs.weight = data.table::setDF(setNames(lapply(q, function(grid) {
  f.val = data[[feature]]
  xs.dist =exp(-0.5 * ((grid - f.val)^2) / (kernel.width^2))
}), q))
weighted_ice_vals=matrix(NA,nrow=nrow(ice_vals),ncol=(h+1))
for(i in 1:(h+1)){
  weighted_ice_vals[,i]=ice_vals[,i] * xs.weight[,i]/sum(xs.weight[,i],na.rm=TRUE)
}
ice_conditional=as.vector(t(weighted_ice_vals))
pdp$results$.value[(h+2):length(pdp$results$.value)]=ice_conditional
pdp$plot()

ice$results$.value=ice_conditional
ice$plot()
# library(yaImpute)
# install.packages("yaImpute")
# eff <- FeatureEffect$new(pred, feature = c("x2", "x4"),method="ale")
# eff$plot()
# eff$plot(show.data = TRUE)