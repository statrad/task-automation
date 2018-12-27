###############################
###Univarable Cox regression###
###############################

# This is the function for doing repeated univariable Cox regression
# It returns a neat-summarized relative risk table similar to reporting format 
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

# 4. "survival" package must be installed in advance.

# 5. Parameter description : 
#    -dat : data (data.frame format) which contains independent variables
#    -var.list : There are three opions are available.
#                1) var.list = "var1" : summary up var1 to end column of dataset
#                2) var.list = c("var1","var2) : summary up var1 to var2
#                3) var.list = c("var1","var2","var3"...) - more than 2var : summary up to only desginated var.
#    -mysurv : Survival object, (time,indicator) , (time1,time2,indicator) 
#              those of 3 cases are available.
#    -filename : "filename", in case of specfied working directory  or full path for the file.

#########################NOTE##########################

# There are two errors you may got while you use.
# 1.  Error: Can't use matrix or array for column indexing  
#     - It may occurs if you put the data as a "tibble" type or something. Please convert it as a data.frame.
# 2.  Error in summary(a)[, 4] : incorrect number of dimensions
#     - This is more typical error. It returns this message when you give it wrong parameter(variable name).
#     - Please confirm your column names by command "colnames(data)" and input right column name.
# 3. It carries out porportionalty test just before write values on the table.
#    When the proportionality assumption is not valid, the values are not be inputted the output table and
#    "non-proportionality" is wrtten at the remarks tab.
######################################################

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
  for(i in index){a
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

#### Example Code ####

#data(kidney) #load the data
#kidney<-kidney[,-6]

#kidney$sex<-as.factor(kidney$sex) #convert the variable as proper form.(num->factor)

#unicox(kidney,"age",c("time","status"))#do the uni-var cox regression

#mysurv<-Surv(kidney$time,kidney$status)#make the surv object.
#unicox(kidney,c("age","sex"),mysurv)#do the uni-var cox regression with surv obj.

#######################
