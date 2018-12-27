unicox<-function(dat,var.list,mysurv,filename="coxph"){
  library(survival)
  rownum<-0
  index<-c()
  
  if(!is.Surv(mysurv)){
    if(length(mysurv)==2){
      t1<-dat[,grep(mysurv[1],colnames(dat))]
      d<-dat[,grep(mysurv[2],colnames(dat))]
      mysurv<-Surv(t1,d)
    }else{
      t1<-dat[,grep(mysurv[1],colnames(dat))]
      t2<-dat[,grep(mysurv[2],colnames(dat))]
      d<-dat[,grep(mysurv[3],colnames(dat))]
      mysurv<-Surv(t1,t2,d)
    }
  }
  
  if(length(var.list)==1){
    index <- seq(grep(var.list,colnames(dat)),ncol(dat))
  }else if(length(var.list)==2){
    index<- seq(grep(var.list[1],colnames(dat)),grep(var.list[2],colnames(dat)))
  }else{
    for(i in 1:length(var.list)){
      index<-c(index,grep(var.list[i],colnames(dat)))
    }
  }
  for(i in index){
    if(is.factor(dat[,i])){
      rownum<-rownum+length(levels(dat[,i]))-1
    }else{
      rownum<-rownum+1
    }
  }
  
  table <- matrix(NA,ncol = 4,nrow= rownum)
  colnames(table)<-c("variable","relative risk","p-value","remarks")
  j<-1
  
  for(i in index){ 
    
    fit<-coxph(mysurv~dat[,i])
    zph<-cox.zph(fit)
     
    if(zph$table[nrow(zph$table),3]>=0.05){
      for(k in 1 : length(fit$coefficients)){
        table[j,"variable"]<-paste(names(dat)[i]," ",levels(dat[,i])[k])
        table[j,"p-value"]<-round(coef(summary(fit))[,5][k],3)
        table[j,"relative risk"]<-paste(round(exp(coef(fit))[k],3)," (",round(exp(confint(fit))[k,1],3),", ",round(exp(confint(fit))[k,2],3),")")
         j<-j+1
      }   
    }else{
      table[j,"variable"]<-names(dat)[i]
      table[j,"remarks"]<-"non-proportional"
      j<-j+1
    }
    
  }
  write.csv(table,paste(filename,".csv",sep=""),row.names=F)
  return(table)
}

library(survival)

