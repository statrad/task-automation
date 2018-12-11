#table 1
setwd("c:/users/ysbios/desktop")
#########

tb1<- function (dat,group,filename="table1"){
  group<-as.factor(group)
  glev<-levels(group)
  table <- matrix(NA,ncol = 4,nrow= 100)
  colnames(table)<-c("variable",levels(as.factor(group)),"p-value")
  j<-1 #row indicator.
  
  for(i in 2 : ncol(dat)){
    dat[,i]<-as.factor(dat[,i])
    lev<-levels(dat[,i])
    if(identical(dat[,i],group)){
      next
    }
    if(length(levels(as.factor(dat[,i])))<=5){ # when the variable is a categorical var.
      for(l in j : (j+length(lev)-1)){
        table[l,1]<-paste(colnames(dat)[i]," ",lev[l-j+1])
      }
      contingency<-table(dat[,i],group)
      print(contingency)
      if(sum(suppressWarnings(chisq.test(contingency)$expected<5))>=1){ #expected value <5 , exact test.
        table[j,"p-value"]<-round(fisher.test(contingency)$p.value,3)
      }else{ #expected value >5, chisq test.
        table[j,"p-value"]<-round(chisq.test(contingency)$p.value,3)
      }
      
      prop<-round(prop.table(contingency,margin=2)*100,1)
      
      for(k in j:(j+length(lev)-1)){
         table[k,2]<-paste(contingency[(k-j+1),1]," (",prop[(k-j+1),1],"%)")
         table[k,3]<-paste(contingency[(k-j+1),2]," (",prop[(k-j+1),2],"%)")
      }          
      
      j<-j+length(lev)
      
    }
    
    if(length(levels(as.factor(dat[,i])))>5){
      dat[,i]<-as.numeric(dat[,i])
      nom.p<-shapiro.test(dat[,i])$p.value#chcek the normality
      tgr1<-as.numeric(subset(dat[,i],group==as.numeric(glev[1])))
      tgr2<-as.numeric(subset(dat[,i],group==as.numeric(glev[2])))
      if(nom.p>0.05){#normality assumption is valid
        table[j,1]<-colnames(dat)[i]
        table[j,"p-value"]<-round(t.test(tgr1,tgr2)$p.value,3)
        table[j,2]<-paste(round(mean(tgr1,na.rm=T))," ¡¾ ",round(sd(tgr1,na.rm=T),3))
        table[j,3]<-paste(round(mean(tgr2,na.rm=T))," ¡¾ ",round(sd(tgr2,na.rm=T),3))
      }else{#not valid
        table[j,1]<-colnames(dat)[i]
        table[j,"p-value"]<-round(wilcox.test(tgr1,tgr2)$p.value,3)
        table[j,2]<-paste(round(median(tgr1,na.rm=T),3)," (",round(quantile(tgr1,na.rm=T)[2],3),", ",round(quantile(tgr1,na.rm=T)[4],3),")")
        table[j,3]<-paste(round(median(tgr2,na.rm=T),3)," (",round(quantile(tgr2,na.rm=T)[2],3),", ",round(quantile(tgr2,na.rm=T)[4],3),")")
      }
      j<-j+1
    }
    
  }
  write.csv(table,paste(filename,".csv",sep=""),row.names = F)
  return(table)
}



#####################
#logistic regression#
####################
unilog<-function(dat,y,filename="result"){
  library(doBy)
  library(logistf)
  table <- matrix(NA,ncol = 4,nrow= 100)
  colnames(table)<-c("variable","odds ratio","p-value","remarks")
  j<-1

  for(i in 2:ncol(dat)){ 
      
    if(identical(dat[,i],y)){
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

data
unilog(data,data$TICI2bor3_1)
fit<-glm(data$TICI2bor3_1~as.factor(data$BALLOON_distal1),family = binomial)
coef(fit)

#################
#cox regression#
################
