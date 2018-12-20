
#################################################################
########## Setting ##############################################
#################################################################
#check and set up the variable type first
#NA's are not allowed in the feature variables
#################################################################

setwd("C:/Users/gradbios/Google 드라이브/Automation")

pMiss <- function(x) {sum(is.na(x))/length(x)*100}



## survival outcome
dat <- read.csv("survival.csv")
str(dat)
dat$EOR_T2 <- as.factor(dat$EOR_T2)
dat$WHO <- as.factor(dat$WHO)
which(apply(dat,2,pMiss)!=0)
dat <- dat[,-which(colnames(dat)=="Diffusion")]

featureSel(dat, "survival", "OS", "Death", "major_axis", "Calvarial_remodeling", seed=20181120) 



## binary outcome
dat <- read.csv("binary.csv")
str(dat)
dat$SR <- as.factor(dat$SR)
which(apply(dat,2,pMiss)!=0)
dat <- dat[,-which(colnames(dat)=="cytology_code")]

featureSel(dat, "binary", "SR", clinical1="size", clinical2="age", clinical3="cytology_code", "ene_1_0", "LL_lrhgle_55_135", seed=20181120) 




#################################################################
data=dd
out.type="binary"
out1="SR"
feature1="ene_1_0"
feature2="LL_lrhgle_55_135"



#################################################################
########## Automation ###########################################
#################################################################
#out.type: "binary", "survival"
#out1: column name of outcome variable (time variable for survival outcome)
#out2: column name of event variable (only for survival outcome)
#clinical: clinical feature name
#feature1: first feature name
#feature2: last feature name
#seed: seed number, current date YYYYMMDD(default)
#################################################################

featureSel <- function(data, out.type, out1, out2=NA, clinical1=NA, clinical2=NA, clinical3=NA, clinical4=NA, feature1, feature2, 
                       seed=NA, filename="Output") 
{
  
  ## Install packages
  list.of.packages <- c("survival", 
                        "caret",
                        "glmnet",
                        "glmnetUtils",
                        "stabs",
                        "rms",
                        "randomForestSRC",
                        "ggRandomForests",
                        "randomForest",
                        "progress",
                        "xlsx",
                        "doParallel")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {install.packages(new.packages)}
  
  library(survival)
  library(caret)
  library(glmnet)
  library(glmnetUtils)
  library(stabs)
  library(rms)
  library(randomForestSRC)
  library(ggRandomForests)
  library(randomForest)
  library(progress)
  library(xlsx)
  library(doParallel)
  registerDoParallel(detectCores())
  detectCores() #the number of logical CPU cores on Window
  options(rf.cores=detectCores(), mc.cores=detectCores()) #Cores for parallel processing
  
  
  ## Dataset for feature selection
  dd <- data[,c(which(names(data) %in% c(out1,out2,clinical1,clinical2,clinical3,clinical4)), which(names(dat)==feature1):which(names(dat)==feature2))]

  
  ## Set seed
  if(is.na(seed)) {seed=format(Sys.Date(), "%Y%b%d"); set.seed(seed)}
  

  
  #################################################################
  # sample code
  #################################################################
  tmp <- function(out.type,out1=NA,s=NA){
    if(s==1) print("imp")
    x=11
    return(x)
  }
  tmp2 <- function(out.type,out1=NA,s=NA) {
    set.seed(s)
    return(c(outcome.type,out1,s))
  }
  tmp("survival",s=1)
  tmp2("survival","sv")
  #################################################################
  
  ## survival
  if(outcome.type=="survival") {
    f.LASSO <- featureLSE(dd, out.type, out1, out2, seed, B=10, "Lasso")
    f.STAB <- "Stability selection is not available for the survival outcome"
    f.ELNET <- featureLSE(dd, out.type, out1, out2, seed, B=10, "Elastic net")
    f.RF <- featureRF(dd, out.type, out1, out2, seed)
  }
  
  
  ## binary
  if(outcome.type=="binary") {
    f.LASSO <- featureLSE(dd, out.type, out1, seed=as.integer(seed), B=10, "Lasso")
    f.STAB <- featureLSE(dd, out.type, out1, seed, B=10, "Stability")
    f.ELNET <- featureLSE(dd, out.type, out1, seed, B=10, "Elastic net")
    f.RF <- featureRF(dd, out.type, out1, seed)
  }
  
  write.xlsx(f.LASSO, file=paste(filename,"_",out1,".xlsx",sep=""), sheetName="Lasso")
  write.xlsx(f.STAB, file=paste(filename,"_",out1,".xlsx",sep=""), sheetName="Stability", append=TRUE)
  write.xlsx(f.ELNET, file=paste(filename,"_",out1,".xlsx",sep=""), sheetName="Elastic net", append=TRUE)
  write.xlsx(f.RF, file=paste(filename,"_",out1,".xlsx",sep=""), sheetName="Random survival forest", append=TRUE)
}





#################################################################
########## featureLSE ###########################################
#################################################################
#B: number of repeated permutation, 100(default)
#method: selection mehtod, "Lasso"(default),"Stability","Elastic net"

#rs.id: resampled id
#data.b: repeated data
#features.b: repeated features 

#type.measure: "deviance" for Cox, "auc" for two-class logistic
#################################################################
featureLSE <- function(data, out.type, out1, out2=NA, seed, B=100, method="LASSO") {
  
  data.b <- list() 
  features.b <- list()
  rs.id <- matrix(0, nrow=nrow(data), ncol=B)
  fit.cv.b <- list()
  min.coef.b <- list()
  
  set.seed(seed)
  features <- data[,-which(names(data) %in% c(out1,out2))]
  for(i in 1:B) {
    rs.id[,i] <- sample(1:nrow(features), nrow(features), replace=T)
    data.b[[i]] <- data[rs.id[,i],]
    features.b[[i]] <- as.matrix(features[rs.id[,i],])
    
    ## Lasso
    if(method=="Lasso") {
      if(out.type=="survival") {
        fit.cv.b[[i]] <- cv.glmnet(features.b[[i]], Surv(data.b[[i]][,out1], data.b[[i]][,out2]), family="cox", type.measure="deviance")
      }
      if(out.type=="binary") {
        fit.cv.b[[i]] <- cv.glmnet(features.b[[i]], data.b[[i]][,out1], family="binomial", type.measure="auc")
      }
      min.coef.b[[i]] <- coef(fit.cv.b[[i]], s="lambda.min")
    }
    
    
    ## Stability
    if(method=="Stability") {
      if(out.type=="survival") {
        return("Stability selection is not available for the survival outcome")
      }
      if(out.type=="binary") {
        x <- model.matrix(~.,data=db[,5:ncol(db)])[,-1] #imaging features
        y <- db$cancertype
        set.seed(180905)
        stab.glmnet <- stabsel(x,y,fitfun=glmnet.lasso,PFER=1,cutoff=0.75)
        stab.glmnet$selected
        plot(stab.glmnet,main="LASSO: glmnet",type="maxsel",cex=0.2,n=50,cex.axis=0.8,ymargin=21)
      }
      min.coef.b[[i]] <- fit.cv.b[[i]]$selected
      or
      min.coef.b[[i]] <- coef(fit.cv.b[[i]], s="lambda.min")
    }
    
    
    ## Elastic net
    if(method=="Elastic net") {
      if(out.type=="survival") {
        elnet.cv <- list()
        alpha <- seq(0, 1, len = 11)^3
        elnet.dev <- rep(0, length(alpha))
        
        elnet.cv[[i]] <- cva.glmnet(features.b[[i]], Surv(data.b[[i]][,out1], data.b[[i]][,out2]), family="cox", type.measure="deviance")
        for(k in 1:11) { elnet.dev[k] <- min(elnet.cv[[i]]$modlist[[k]]$cvm) }
        fit.cv.b[[i]] <- glmnetUtils::glmnet(features.b[[i]], Surv(data.b[[i]][,out1], data.b[[i]][,out2]), family="cox",
                                             alpha=alpha[which.min(elnet.dev)], lambda=elnet.cv[[i]]$modlist[[which.min(elnet.dev)]]$lambda.min)
        min.coef.b[[i]] <- fit.cv.b[[i]]$beta
      }
      if(out.type=="binary") {
        
      }
    }
  }
  
  min.coef <- matrix(NA, ncol=B, nrow=ncol(features), dimnames=list(rownames(min.coef.b[[1]]), 1:B))
  min.coef <- as.data.frame(min.coef)
  
  for(i in 1:B) {min.coef[,i] <- min.coef.b[[i]]}
  min.coef <- as.matrix(min.coef)
  
  m <- c(0); abs_m <- c(0); se <- c(0); rsd <- c(0); abs_rsd <- c(0)
  for(i in 1:nrow(min.coef)) {
    m[i] <- mean(as.numeric(min.coef[i,]))
    abs_m[i] <- abs(m[i])
    se[i] <- sd(as.numeric(min.coef[i,]))
    rsd[i] <- (se[i]/m[i])
    abs_rsd[i] <- (se[i]/abs_m[i])
  }
  
  return(cbind(min.coef, mean=m, abs_m=abs_m, se=se, rsd=rsd, abs_rsd=abs_rsd))
}





#################################################################
#### Random Survival Forests
#################################################################
#ntree: total number of trees
#mtry: number of variables entering in each division as candidates
#      (sqrt(p) in classification, survival / p/3 in regression)
#nodesize: minimum number of observations in a terminal node (3 in survival)
#nsplit: number of random points to explore in continous predictors
#importance=T: prints out variable importance ranking
#proximity=T: compute this metric
#################################################################
featureRF <- function(data, out.type, out1, out2=NA, seed=seed) {

  if(out.type=="survival") {
    form <- as.formula(paste0("Surv(",out1,",",out2,")~."))
    o <- tune(form, data, seed=as.numeric(seed))
    return(var.select(rfsrc(form, data, seed=as.numeric(seed), ntree=o$rf$ntree, mtry=o$rf$mtry, nodesize=o$rf$nodesize, importance=TRUE))$varselect)
  }
  
  if(out.type=="binary") {
    o <- tune()
    selected.var
  }
}




