#####################
#logistic regression#
####################
unilog<-function(dat,pnum,y,filename="result"){
  library(doBy)
  library(logistf)
  
  rownum<-0
  for(i in 1 : ncol(dat)){
    if(identical(dat[,i],y)|identical(dat[,i],pnum)){next}
    if(is.factor(dat[,i])){
      rownum<-rownum+length(levels(dat[,i]))-1
    }else{
      rownum<-rownum+1
    }
  }
  table <- matrix(NA,ncol = 4,nrow= rownum)
  colnames(table)<-c("variable","odds ratio","p-value","remarks")
  j<-1

  for(i in 1:ncol(dat)){ 
      
    if(identical(dat[,i],y)|identical(dat[,i],pnum)){
       next
    }
    
    table[j,"variable"]<-names(dat)[i]
    t.temp<-table(dat[,i],y)
    
    if(sum(t.temp==0)>=1&length(levels(as.factor(dat[,i])))<=5){
      fit<-logistf(y~dat[,i])
      table[j,"p-value"]<-round(fit$prob[2],3)
      table[j,"odds ratio"]<-paste(round(exp(coef(fit))[2],3)," (",round(exp(confint(fit))[2,1],3),", ",round(exp(confint(fit))[2,2],3),")")
      table[j,"remarks"]<-"firth logistic"
    }else{
      fit<-glm(y~dat[,i],family = binomial)
      table[j,"p-value"]<-round(coef(summary(fit))[,4][2],3)
      table[j,"odds ratio"]<-paste(round(exp(coef(fit))[2],3)," (",round(exp(confint(fit))[2,1],3),", ",round(exp(confint(fit))[2,2],3),")")      
    }
     
    j<-j+1
    
  }
  write.csv(table,paste(filename,".csv",sep=""),row.names=F)
  return(table)
}
