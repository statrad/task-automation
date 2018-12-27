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
#    - In case of continuous, use "as.numeric(x)" for trasnformation.

# 3. It returns file via "csv" format so that the working directory must be designated
#    or state Full path of it's result to "filename" parameter.

# 4. "doBy" , "logistf" must be installed in advance.

# 5. Parameter description : 
#    -dat : data (data.frame format) which contains independent variables
#    -var.list : There are three opions are available.
#                1) var.list = "var1" : summary up var1 to end column of dataset
#                2) var.list = c("var1","var2) : summary up var1 to var2
#                3) var.list = c("var1","var2","var3"...) - more than 2var : summary up to only desginated var.
#    -y : response variable 
#    -filename : "filename", in case of specfied working directory  or full path for the file.

#########################NOTE##########################

# There are two errors you may got while you use.
# 1.  Error: Can't use matrix or array for column indexing  
#     - It may occurs if you put the data as a "tibble" type or something. Please convert it as a data.frame.
# 2.  Error in summary(a)[, 4] : incorrect number of dimensions
#     - This is more typical error. It returns this message when you give it wrong parameter(variable name).
#     - Please confirm your column names by command "colnames(data)" and input right column name.

######################################################

### Example codes will be offered.###

unilog<-function(dat,var.list,y,filename="result"){
  library(doBy)
  library(logistf)
  
  rownum<-0
  index<-c()
  y.name<-y
  y<-dat[,grep(y,colnames(data))]
  
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
    if(colnames(dat)[i]==y.name){next}
    if(is.factor(dat[,i])){
      rownum<-rownum+length(levels(dat[,i]))-1
    }else{
      rownum<-rownum+1
    }
  }
  
  table <- matrix(NA,ncol = 4,nrow= rownum)
  colnames(table)<-c("variable","odds ratio","p-value","remarks")
  j<-1
  
  for(i in index){ 
    if(colnames(dat)[i]==y.name){next}
    
    
    t.temp<-table(dat[,i],y)
    
    if(sum(t.temp==0)>=1&length(levels(as.factor(dat[,i])))<=5){
      fit<-logistf(y~dat[,i])
 
      for(k in 2 : length(fit$coefficients)){
         table[j,"variable"]<-paste(names(dat)[i]," ",levels(dat[,i])[k])
         table[j,"p-value"]<-round(fit$prob[k],3)
         table[j,"odds ratio"]<-paste(round(exp(coef(fit))[k],3)," (",round(exp(confint(fit))[k,1],3),", ",round(exp(confint(fit))[k,2],3),")")
         table[j,"remarks"]<-"firth logistic"
         j<-j+1
      }   
      
    }else{
      fit<-glm(y~dat[,i],family = binomial)
      
       for(k in 2 : length(fit$coefficients)){
         table[j,"variable"]<-paste(names(dat)[i]," ",levels(dat[,i])[k])
        table[j,"p-value"]<-round(coef(summary(fit))[,4][k],3)
        table[j,"odds ratio"]<-paste(round(exp(coef(fit))[k],3)," (",round(exp(confint(fit))[k,1],3),", ",round(exp(confint(fit))[k,2],3),")")
        j<-j+1
      }   
      
    }
    
    
  }
  write.csv(table,paste(filename,".csv",sep=""),row.names=F)
  return(table)
}



