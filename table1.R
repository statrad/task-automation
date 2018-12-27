###################################################
#Making table 1(paper like descriptive table form)#
###################################################

# This is the function for making Table1 
# It returns a neat-summarized table resemble reporting format 
# which is used for reporting summary statstics in conventional papers.
# It ignores ID variable and group variable so that 
# able to handle whole dataset at once.

# To get proper result from it, there are few things that 
# you should follow and it is written below.

# 1. The format for input data must be a data frame
#    - It can be done by simple command, "as.data.frame"

# 2. The type of variable must be clarified.
#    - If you want to treat the variable as "categorical" variable
#      you should transform it as factor and it accomplished by R command, "as.factor(x)".
#    - In case of continuous, use "as.numeric(x)" for trasnformation.

# 3. It returns file via "csv" format so that the working directory must be designated
#    or state Full path of it's result to "filename" parameter.

# 5. Parameter description : 
#    -dat : data (data.frame format)
#    -var.list : There are three opions are available.
#                1) var.list = "var1" : summary up var1 to end column of dataset
#                2) var.list = c("var1","var2) : summary up var1 to var2
#                3) var.list = c("var1","var2","var3"...) - more than 2var : summary up to only desginated var.
#    -group : group variable 
#    -filename : "filename", in case of specfied working directory  or full path for the file.

#########################NOTE##########################

# There are two errors you may got while you use.
# 1.  Error: Can't use matrix or array for column indexing  
#     - It may occurs if you put the data as a "tibble" type or something. Please convert it as a data.frame.
# 2.  Error in summary(a)[, 4] : incorrect number of dimensions
#     - This is more typical error. It returns this message when you give it wrong parameter(variable name).
#     - Please confirm your column names by command "colnames(data)" and input right column name.

# It corresponds ONLY "2 Groups".
######################################################

### Example codes will be offered.###

tb1<- function (dat,var.list,group,filename="table1"){
  rownum<-0
  index<-c()
  g.name<-group
  g.index<-grep(group,colnames(dat))
  group<-as.factor(dat[,g.index])
  glev<-levels(group)
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
    if(colnames(dat)[i]==g.name){next}
    if(is.factor(dat[,i])){
      rownum<-rownum+length(levels(dat[,i]))
    }else{
      rownum<-rownum+1
    }
  }
  table <- matrix(NA,ncol = 4,nrow= rownum)
  colnames(table)<-c("variable",levels(as.factor(group)),"p-value")
  j<-1 #row indicator.
  
  for(i in index){
    if(colnames(dat)[i]==g.name){next}
    
    lev<-levels(dat[,i])
    
    if(is.factor(dat[,i])){ # when the variable is a categorical var.
      for(l in j : (j+length(lev)-1)){
        table[l,1]<-paste(colnames(dat)[i]," ",lev[l-j+1])
      }
      contingency<-table(dat[,i],dat[,g.index])
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
