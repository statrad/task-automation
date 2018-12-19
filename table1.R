#table 1
#########

tb1<- function (dat,pnum,group,filename="table1"){
  group<-as.factor(group)
  glev<-levels(group)
  rownum<-0
  for(i in 1 : ncol(dat)){
    if(identical(dat[,i],group)|identical(dat[,i],pnum)){next}
    if(is.factor(dat[,i])){
      rownum<-rownum+length(levels(dat[,i]))
    }else{
      rownum<-rownum+1
    }
  }
  table <- matrix(NA,ncol = 4,nrow= rownum)
  colnames(table)<-c("variable",levels(as.factor(group)),"p-value")
  j<-1 #row indicator.
  
  for(i in 1 : ncol(dat)){
    if(identical(dat[,i],group)|identical(dat[,i],pnum)){
      next
    }
    
    lev<-levels(dat[,i])
    
    if(is.factor(dat[,i])){ # when the variable is a categorical var.
      for(l in j : (j+length(lev)-1)){
        table[l,1]<-paste(colnames(dat)[i]," ",lev[l-j+1])
      }
      contingency<-table(dat[,i],group)
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
    
    if(!is.factor(dat[,i])){
      dat[,i]<-as.numeric(dat[,i])
      nom.p<-shapiro.test(dat[,i])$p.value#chcek the normality
      tgr1<-as.numeric(subset(dat[,i],group==as.numeric(glev[1])))
      tgr2<-as.numeric(subset(dat[,i],group==as.numeric(glev[2])))
      if(nom.p>0.05){#normality assumption is valid
        table[j,1]<-colnames(dat)[i]
        table[j,"p-value"]<-round(t.test(tgr1,tgr2)$p.value,3)
        table[j,2]<-paste(round(mean(tgr1,na.rm=T))," ± ",round(sd(tgr1,na.rm=T),3))
        table[j,3]<-paste(round(mean(tgr2,na.rm=T))," ± ",round(sd(tgr2,na.rm=T),3))
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


