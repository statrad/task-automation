################################
#Univarable logistic regression#
################################

# This is the function for doing repeated univariable logistic regression
# It returns a neat-summarized Odds ratio table similar to reporting format 
# which is used for reporting summary statstics in conventional papers.
# It ignores ID variable and response variable so that 
# able to handle whole dataset at once.

# To get proper result from it, there are few things that 
# you should follow and it is written below.

# 1. The format for input data must be a data frame
#    - It can be done by simple command, "as.data.frame"

# 2. The type of variable must be clarified.
#    - If you want to treat the variable as "categorical" variable
#      you should transform it as factor and it accomplished by R command, "as.factor(x)".
#    - In case of Continuous, use "as.numeric(x)" for trasnformation.

# 3. It returns file via "csv" format so that the working directory must be designated
#    or state Full path of it's result to "filename" parameter.

# 4. "doBy" , "logistf" must be installed in advance.

# 5. Parameter description : 
#    -dat : data (data.frame format) which contains independent variables
#    -pnum : patient(subject) number
#    -y : response variable 
#    -filename : "filename", in case of specfied working directory  or full path for the file.

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
