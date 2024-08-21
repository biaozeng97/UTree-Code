library(GGUM)
library(SimDesign)
library(truncnorm)
library(MASS)
library(mirt)
library(RPushbullet)
library(ggplot2)
library(hrbrthemes)
library(matrixStats)


setwd("C:/Unfolding Tree Model/ORDUTree_1")


############################################################################################
###################### Part 1: data are generated from UIRTree #############################
############################################################################################

## 1. Study Design

Design<-expand.grid(Items = c(5,10,20),
                    Samplesize = c(500,1000,2000),
                    Categories = c(4),
                    Cor = c(0,0.3,0.6))         ## Correlation between focal factor and ERS
## discrimination on the ERS (square root of the Variance of ERS)            
##### Is the Alpha setting needed? #####




## 2. Generate Data
Generate.sample<-function(condition, fixed_objects = NULL) {
  
  N<-condition$Samplesize
  I<-condition$Items
  C<-condition$Categories
  R<-condition$Cor
  A<-condition$Alpha

  
  #### parameter values based on a GGUM model for measuring the theta in Node 1
  alpha1 <- runif(n = I, min = 0.5, max = 2)
  delta1 <- rtruncnorm(n = I,mean=0, sd=1,a=-2, b=2)
  tau1.2 <- runif(n = I, min = 0.5, max = 2)
  tau1.1 <- rep(0,I)
  tau1 <-data.frame(tau1.1,tau1.2)
  
  #### parameter values based on a GGUM model for measuring the theta in Node 2 and 3
  alpha2 <- runif(n = 2*I, min = 0.5, max = 2) 
  delta2 <- rtruncnorm(n = 2*I,mean=0, sd=1,a=-2, b=2)
  tau2.2 <- runif(n = 2*I, min = 0.5, max = 2)
  tau2.1 <- rep(0,2*I)
  tau2 <-data.frame(tau2.1,tau2.2)  
  
  theta <- rnorm(n = N, mean=0, sd=1)
 
    
  ## Data generation function 
  data.gen <- function(alpha1,delta1,tau1,alpha2,delta2,tau2,
                       N,I,theta){

    resp <- matrix(NA,nrow = N,ncol = I)
    
    for (n in 1:N){
      
      for (i in 1:I) {
        
        sum1.0 <- 1+exp(alpha1[i]*3*(theta[n]-delta1[i]))
        sum1.1 <- exp(alpha1[i]*(1*(theta[n]-delta1[i])+sum(tau1[i,1:2])))+exp(alpha1[i]*(2*(theta[n]-delta1[i])+sum(tau1[i,1:2])))
        
        sum2.0 <- 1+exp(alpha2[i]*3*(theta[n]-delta2[i]))
        sum2.1 <- exp(alpha2[i]*(1*(theta[n]-delta2[i])+sum(tau2[i,1:2])))+exp(alpha2[i]*(2*(theta[n]-delta2[i])+sum(tau2[i,1:2])))
        
        sum3.0 <- 1+exp(alpha2[I+i]*3*(theta[n]-delta2[I+i]))
        sum3.1 <- exp(alpha2[I+i]*(1*(theta[n]-delta2[I+i])+sum(tau2[I+i,1:2])))+exp(alpha2[I+i]*(2*(theta[n]-delta2[I+i])+sum(tau2[I+i,1:2])))
        
        ## Probability for the node 1 (ideal point model)
        p1.0 <- sum1.0/(sum1.0 + sum1.1)
        p1.1 <- sum1.1/(sum1.0 + sum1.1)
        
        ## Probability for the node 2 (ideal point model)
        p2.0 <- sum2.0/(sum2.0 + sum2.1)
        p2.1 <- sum2.1/(sum2.0 + sum2.1)
        
        ## Probability for the node 2 (ideal point model)
        p3.0 <- sum3.0/(sum3.0 + sum3.1)
        p3.1 <- sum3.1/(sum3.0 + sum3.1)
        
        ## Probability for observed categories
        p <- c(p1.0*p2.0, # strongly disagree
               p1.0*p2.1, # disagree, UTree set the disagree as 1,  strongly disagree as 0
               p1.1*p3.0, # agree
               p1.1*p3.1) # strongly agree
        
        resp[n, i] <- sample(0:3,size = 1,replace = F,prob = p) # from 0 to 3, the categories are: strongly disagree, disagree, agree, and strongly agree
      }
    }
    
    return(resp)
}
    
    ## Generate data 
    data.4 <- data.gen(alpha1,delta1,tau1,alpha2,delta2,tau2,
                       N,I,theta) # dataset for the GGUM model in node 1, 2 and 3
    mapping.matrix1 <- matrix(c(0,0,1,1,0,1,NA,NA,NA,NA,0,1),nrow=4) # In node 2, the GGUM set the disagree as 1
    data.2 <- matrix(as.vector(mapping.matrix1[(as.vector(data.4)+1),]),nrow = nrow(data.4)) # dataset for the Unfolding models
    # the first I items are to measure whether the responses are in positive(1) or negative(0) direction
    # the second I items are to measure if the responses are in the negative direction, whether they are higher trait (1) or not(0)
    # the last I items are to measure if the responses are in the positive direction, whether they are higher trait (1) or not(0)

    mapping.matrix2 <- matrix(c(0,0,1,1,1,0,NA,NA,NA,NA,0,1),nrow=4) # In node 2, the IRTree set the disagree as 1
    data.6 <- matrix(as.vector(mapping.matrix2[(as.vector(data.4)+1),]),nrow = nrow(data.4)) # dataset for the IRTree models 
    # the first I items are to measure whether the responses are in positive(1) or negative(0) direction
    # the second I items are to measure if the responses are in the negative direction, whether they are extreme(1) or not(0)
    # the last I items are to measure if the responses are in the positive direction, whether they are extreme(1) or not(0)
    
    
    colnames(data.4) <-paste0("V",1:I)
    colnames(data.2) <-paste0("V",1:(3*I))
    colnames(data.6) <-paste0("V",1:(3*I))
    
    
    alpha1.raw <- alpha1
    delta1.raw <- delta1
    tau1.raw <- tau1.2
    alpha2.raw <- alpha2
    delta2.raw <- delta2
    tau2.raw <- tau2.2
    theta1.raw <- theta
    
    dat <- list(data.4,data.2,data.6, alpha1.raw,delta1.raw,tau1.raw,
                alpha2.raw,delta2.raw,tau2.raw,theta1.raw)
    dat
  }
    
    

## 3. Analyze data 

Analyse<-function(condition,dat,fixed_objects=NULL){
  
  N<-condition$Samplesize
  I<-condition$Items
  C<-condition$Categories
  R<-condition$Cor
  A<-condition$Alpha
  
  data.4 <- as.matrix(dat[[1]])
  data.2 <- as.matrix(dat[[2]])
  data.6 <- as.matrix(dat[[3]])
  alpha1.raw <- as.vector(dat[[4]])
  delta1.raw <- as.vector(dat[[5]])
  tau1.raw <- as.vector(dat[[6]])
  alpha2.raw <- as.vector(dat[[7]])
  delta2.raw <- as.vector(dat[[8]])
  tau2.raw <- as.vector(dat[[9]])
  theta1.raw <- as.vector(dat[[10]])

  
  ## Fit models
  if (isTRUE(I==5)) {
    
    mirt.mod.ORDUTree_1 <- "F1=1-15
    PRIOR = (1-15,b1,norm,0,3),(1-15,a1,lnorm,0.2,0.5),(1-15,t1,norm,1,0.5)"
    
    mirt.mod.ERSUTree <- "F1=1-5
    F2=6-15
    COV=F1*F2
    PRIOR = (1-5,b1,norm,0,3),(6-15,d,norm,0,3),(1-5,a1,lnorm,0.2,0.5),(6-15,a2,lnorm,0.2,0.5),(1-5,t1,norm,1,0.5)"
    
    mirt.mod.ORD_1 <- "F1=1-15
    PRIOR = (1-15,d,norm,0,3),(1-15,a1,lnorm,0.2,0.5)"
    
    mirt.mod.ORD_2 <- "F1=1-5
    F2=6-15
    COV=F1*F2
    PRIOR = (1-15,d,norm,0,3),(1-5,a1,lnorm,0.2,0.5),(6-15,a2,lnorm,0.2,0.5)"
    
    mirt.mod.ERS <- "F1=1-5
    F2=6-15
    COV=F1*F2
    PRIOR = (1-15,d,norm,0,3),(1-5,a1,lnorm,0.2,0.5),(6-15,a2,lnorm,0.2,0.5)"
  }
  
  
  if (isTRUE(I==10)) {
    
    mirt.mod.ORDUTree_1 <- "F1=1-30
    PRIOR = (1-30,b1,norm,0,3),(1-30,a1,lnorm,0.2,0.5),(1-30,t1,norm,1,0.5)"
    
    mirt.mod.ERSUTree <- "F1=1-10
    F2=11-30
    COV=F1*F2
    PRIOR = (1-10,b1,norm,0,3),(11-30,d,norm,0,3),(1-10,a1,lnorm,0.2,0.5),(11-30,a2,lnorm,0.2,0.5),(1-10,t1,norm,1,0.5)"
    
    mirt.mod.ORD_1 <- "F1=1-30
    PRIOR = (1-30,d,norm,0,3),(1-30,a1,lnorm,0.2,0.5)"
    
    mirt.mod.ORD_2 <- "F1=1-10
    F2=11-30
    COV=F1*F2
    PRIOR = (1-30,d,norm,0,3),(1-10,a1,lnorm,0.2,0.5),(11-30,a2,lnorm,0.2,0.5)"
    
    mirt.mod.ERS <- "F1=1-10
    F2=11-30
    COV=F1*F2
    PRIOR = (1-30,d,norm,0,3),(1-10,a1,lnorm,0.2,0.5),(11-30,a2,lnorm,0.2,0.5)"
  }
  
  
  if (isTRUE(I==20)) {
    
    mirt.mod.ORDUTree_1 <- "F1=1-60
    PRIOR = (1-60,b1,norm,0,3),(1-60,a1,lnorm,0.2,0.5),(1-60,t1,norm,1,0.5)"
    
    mirt.mod.ERSUTree <- "F1=1-20
    F2=21-60
    COV=F1*F2
    PRIOR = (1-20,b1,norm,0,3),(21-60,d,norm,0,3),(1-20,a1,lnorm,0.2,0.5),(21-60,a2,lnorm,0.2,0.5),(1-20,t1,norm,1,0.5)"
    
    mirt.mod.ORD_1 <- "F1=1-60
    PRIOR = (1-60,d,norm,0,3),(1-60,a1,lnorm,0.2,0.5)"
    
    mirt.mod.ORD_2 <- "F1=1-20
    F2=21-60
    COV=F1*F2
    PRIOR = (1-60,d,norm,0,3),(1-20,a1,lnorm,0.2,0.5),(21-60,a2,lnorm,0.2,0.5)"
    
    mirt.mod.ERS <- "F1=1-20
    F2=21-60
    COV=F1*F2
    PRIOR = (1-60,d,norm,0,3),(1-20,a1,lnorm,0.2,0.5),(21-60,a2,lnorm,0.2,0.5)"
  }
  
  mod.ORDUTree_1 <- mirt(data = data.2,model = mirt.mod.ORDUTree_1, itemtype = "ggum")
  mod.ORD_1 <- mirt(data = data.2,model = mirt.mod.ORD_1, itemtype = "2PL")
  mod.ORD_2 <- mirt(data = data.2,model = mirt.mod.ORD_2, itemtype = "2PL")
  mod.ERS <- mirt(data = data.6,model = mirt.mod.ERS, itemtype = "2PL")
  
  
  # report item paramter and transform d to beta in 2PL models
  para.ERS_origin <- coef(mod.ERS,simplify=T)$item
  ERS_beta1 <- para.ERS_origin[,3]/(-para.ERS_origin[,1])
  ERS_beta2 <- para.ERS_origin[,3]/(-para.ERS_origin[,2])
  para.ERS <- cbind(para.ERS_origin,ERS_beta1,ERS_beta2)
  
  para.ORD_1_origin <- coef(mod.ORD_1,simplify=T)$item
  ORD_1_beta <- para.ORD_1_origin[,2]/(-para.ORD_1_origin[,1]) 
  para.ORD_1 <- cbind(para.ORD_1_origin,ORD_1_beta)
  para.ORD_2_origin <- coef(mod.ORD_2,simplify=T)$item
  ORD_2_beta1 <- para.ORD_2_origin[,3]/(-para.ORD_2_origin[,1])
  ORD_2_beta2 <- para.ORD_2_origin[,3]/(-para.ORD_2_origin[,2])
  para.ORD_2 <- cbind(para.ORD_2_origin,ORD_2_beta1,ORD_2_beta2)
  
  para.ORDUTree_1 <- coef(mod.ORDUTree_1,simplify=T)$items
  
  ##     theta and eta    ##
  theta.ERS <- fscores(mod.ERS, full.scores.SE = T)
  theta.ORD_1 <- fscores(mod.ORD_1, full.scores.SE = T)
  theta.ORD_2 <- fscores(mod.ORD_2, full.scores.SE = T)
  theta.ORDUTree_1 <- fscores(mod.ORDUTree_1, full.scores.SE = T)

  
  if(isTRUE(max(abs(para.ORDUTree_1[,2:3]))>=10)) {stop("Delta estimates reached the limit")}
  
  AIC.mod.ERS <- extract.mirt(mod.ERS,"AIC")
  AIC.mod.ORD_1 <- extract.mirt(mod.ORD_1,"AIC")
  AIC.mod.ORD_2 <- extract.mirt(mod.ORD_2,"AIC")
  AIC.mod.ORDUTree_1 <- extract.mirt(mod.ORDUTree_1,"AIC")
  BIC.mod.ERS <- extract.mirt(mod.ERS,"BIC")
  BIC.mod.ORD_1 <- extract.mirt(mod.ORD_1,"BIC")
  BIC.mod.ORD_2 <- extract.mirt(mod.ORD_2,"BIC")
  BIC.mod.ORDUTree_1 <- extract.mirt(mod.ORDUTree_1,"BIC")
  SABIC.mod.ERS <- extract.mirt(mod.ERS,"SABIC")
  SABIC.mod.ORD_1 <- extract.mirt(mod.ORD_1,"SABIC")
  SABIC.mod.ORD_2 <- extract.mirt(mod.ORD_2,"SABIC")
  SABIC.mod.ORDUTree_1 <- extract.mirt(mod.ORDUTree_1,"SABIC")
  
  
  sign.ERS.beta1 <- sign(cor(para.ERS[1:I,6],delta1.raw))
  sign.ERS.beta2 <- sign(cor(para.ERS[(I+1):(3*I),7],delta2.raw))
  sign.ORD_1.beta1 <- sign(cor(para.ORD_1[1:I,5],delta1.raw))
  sign.ORD_1.beta2 <- sign(cor(para.ORD_1[(I+1):(3*I),5],delta2.raw))
  sign.ORD_2.beta1 <- sign(cor(para.ORD_2[1:I,6],delta1.raw))
  sign.ORD_2.beta2 <- sign(cor(para.ORD_2[(I+1):(3*I),7],delta2.raw))
  sign.ORDUTree_1.delta1 <- sign(cor(para.ORDUTree_1[1:I,2], delta1.raw))
  sign.ORDUTree_1.delta2 <- sign(cor(para.ORDUTree_1[(I+1):(3*I),2], delta2.raw))

  
  # transfer theta and eta from wrong direction to right direction with true values
  sign.ERS.theta <- sign(cor(theta.ERS[,1], theta1.raw))
  theta.esti.ERS <- theta.ERS[,1] *sign.ERS.theta
  theta.combine.ERS  <- data.frame(cbind(theta.esti.ERS ,theta1.raw))
  sign.ERS.eta <- sign(cor(theta.ERS[,2], theta1.raw))
  eta.esti.ERS <- theta.ERS[,2]*sign.ERS.eta
  eta.combine.ERS <- data.frame(cbind(eta.esti.ERS,theta1.raw))
  
  sign.ORD_1.theta <- sign(cor(theta.ORD_1[,1], theta1.raw))
  theta.esti.ORD_1 <- theta.ORD_1[,1] *sign.ORD_1.theta
  theta.combine.ORD_1  <- data.frame(cbind(theta.esti.ORD_1 ,theta1.raw))
  
  sign.ORD_2.theta1 <- sign(cor(theta.ORD_2[,1], theta1.raw))
  theta1.esti.ORD_2 <- theta.ORD_2[,1] *sign.ORD_2.theta1
  theta1.combine.ORD_2  <- data.frame(cbind(theta1.esti.ORD_2 ,theta1.raw))
  sign.ORD_2.theta2 <- sign(cor(theta.ORD_2[,2], theta1.raw))
  theta2.esti.ORD_2 <- theta.ORD_2[,2]*sign.ORD_2.theta2
  theta2.combine.ORD_2 <- data.frame(cbind(theta2.esti.ORD_2,theta1.raw))
  
  sign.ORDUTree_1.theta <- sign(cor(theta.ORDUTree_1[,1], theta1.raw))
  theta.esti.ORDUTree_1 <- theta.ORDUTree_1[,1] *sign.ORDUTree_1.theta
  theta.combine.ORDUTree_1  <- data.frame(cbind(theta.esti.ORDUTree_1 ,theta1.raw))
  
  
  # 95% confidence interval (theta and eta) include true value
  theta_95CI.ERS <- cbind(theta.esti.ERS-1.96*theta.ERS[,3],theta.esti.ERS+1.96*theta.ERS[,3])
  eta_95CI.ERS <- cbind(eta.esti.ERS-1.96*theta.ERS[,4],eta.esti.ERS+1.96*theta.ERS[,4])
  theta_95CI.ORD_1 <- cbind(theta.esti.ORD_1-1.96*theta.ORD_1[,2],theta.esti.ORD_1+1.96*theta.ORD_1[,2])
  theta1_95CI.ORD_2 <- cbind(theta1.esti.ORD_2-1.96*theta.ORD_2[,3],theta1.esti.ORD_2+1.96*theta.ORD_2[,3])
  theta2_95CI.ORD_2 <- cbind(theta2.esti.ORD_2-1.96*theta.ORD_2[,4],theta2.esti.ORD_2+1.96*theta.ORD_2[,4])
  
  theta_95CI.ORDUTree_1 <- cbind(theta.esti.ORDUTree_1-1.96*theta.ORDUTree_1[,2],theta.esti.ORDUTree_1+1.96*theta.ORDUTree_1[,2])
 
  
  
  #### Summarize results
  ret<-c(
    Abis.alpha1.ERS <- mean(abs(para.ERS[1:I,1]-alpha1.raw)),
    Abis.alpha2.ERS <- mean(abs(para.ERS[(I+1):(3*I),2]-alpha2.raw)),
    Abis.beta1.ERS <- mean(abs(para.ERS[1:I,6]*sign.ERS.beta1-delta1.raw)),
    Abis.beta2.ERS <- mean(abs(para.ERS[(I+1):(3*I),7]*sign.ERS.beta2-delta2.raw)),
    Abis.cor.ERS <- abs(abs(coef(mod.ERS,simplify=T)$cov[1,2])-R),
    
    Abis.alpha1.ORD_1 <- mean(abs(para.ORD_1[1:I,1]-alpha1.raw)),
    Abis.alpha2.ORD_1 <- mean(abs(para.ORD_1[(I+1):(3*I),1]-alpha2.raw)),
    Abis.beta1.ORD_1 <- mean(abs(para.ORD_1[1:I,5]*sign.ORD_1.beta1-delta1.raw)),
    Abis.beta2.ORD_1 <- mean(abs(para.ORD_1[(I+1):(3*I),5]*sign.ORD_1.beta2-delta2.raw)),
    
    Abis.alpha1.ORD_2 <- mean(abs(para.ORD_2[1:I,1]-alpha1.raw)),
    Abis.alpha2.ORD_2 <- mean(abs(para.ORD_2[(I+1):(3*I),2]-alpha2.raw)),
    Abis.beta1.ORD_2 <- mean(abs(para.ORD_2[1:I,6]*sign.ORD_2.beta1-delta1.raw)),
    Abis.beta2.ORD_2 <- mean(abs(para.ORD_2[(I+1):(3*I),7]*sign.ORD_2.beta2-delta2.raw)),
    Abis.cor.ORD_2 <- abs(abs(coef(mod.ORD_2,simplify=T)$cov[1,2])-R),
    
    Abis.alpha1.ORDUTree_1 <- mean(abs(para.ORDUTree_1[1:I,1]-alpha1.raw)),
    Abis.alpha2.ORDUTree_1 <- mean(abs(para.ORDUTree_1[(I+1):(3*I),1]-alpha2.raw)),
    Abis.delta1.ORDUTree_1 <- mean(abs(para.ORDUTree_1[1:I,2]*sign.ORDUTree_1.delta1-delta1.raw)),
    Abis.delta2.ORDUTree_1 <- mean(abs(para.ORDUTree_1[(I+1):(3*I),2]*sign.ORDUTree_1.delta2-delta2.raw)),
    Abis.tau1.ORDUTree_1 <- mean(abs(para.ORDUTree_1[1:I,3]-tau1.raw)),
    Abis.tau2.ORDUTree_1 <- mean(abs(para.ORDUTree_1[(I+1):(3*I),3]-tau2.raw)),
    
    Bias.alpha1.ERS <- mean(para.ERS[1:I,1]-alpha1.raw),
    Bias.alpha2.ERS <- mean(para.ERS[(I+1):(3*I),2]-alpha2.raw),
    Bias.beta1.ERS <- mean(para.ERS[1:I,6]*sign.ERS.beta1-delta1.raw),
    Bias.beta2.ERS <- mean(para.ERS[(I+1):(3*I),7]*sign.ERS.beta2-delta2.raw),
    Bias.cor.ERS <- abs(coef(mod.ERS,simplify=T)$cov[1,2])-R,
    
    Bias.alpha1.ORD_1 <- mean(para.ORD_1[1:I,1]-alpha1.raw),
    Bias.alpha2.ORD_1 <- mean(para.ORD_1[(I+1):(3*I),1]-alpha2.raw),
    Bias.beta1.ORD_1 <- mean(para.ORD_1[1:I,5]*sign.ORD_1.beta1-delta1.raw),
    Bias.beta2.ORD_1 <- mean(para.ORD_1[(I+1):(3*I),5]*sign.ORD_1.beta2-delta2.raw),
    
    Bias.alpha1.ORD_2 <- mean(para.ORD_2[1:I,1]-alpha1.raw),
    Bias.alpha2.ORD_2 <- mean(para.ORD_2[(I+1):(3*I),2]-alpha2.raw),
    Bias.beta1.ORD_2 <- mean(para.ORD_2[1:I,6]*sign.ORD_2.beta1-delta1.raw),
    Bias.beta2.ORD_2 <- mean(para.ORD_2[(I+1):(3*I),7]*sign.ORD_2.beta2-delta2.raw),
    Bias.cor.ORD_2 <- abs(coef(mod.ORD_2,simplify=T)$cov[1,2])-R,
    
    Bias.alpha1.ORDUTree_1 <- mean(para.ORDUTree_1[1:I,1]-alpha1.raw),
    Bias.alpha2.ORDUTree_1 <- mean(para.ORDUTree_1[(I+1):(3*I),1]-alpha2.raw),
    Bias.delta1.ORDUTree_1 <- mean(para.ORDUTree_1[1:I,2]*sign.ORDUTree_1.delta1-delta1.raw),
    Bias.delta2.ORDUTree_1 <- mean(para.ORDUTree_1[(I+1):(3*I),2]*sign.ORDUTree_1.delta2-delta2.raw),
    Bias.tau1.ORDUTree_1 <- mean(para.ORDUTree_1[1:I,3]-tau1.raw),
    Bias.tau2.ORDUTree_1 <- mean(para.ORDUTree_1[(I+1):(3*I),3]-tau2.raw),
    
    RMSE.alpha1.ERS <- RMSE(unname(para.ERS[1:I,1]),alpha1.raw),
    RMSE.alpha2.ERS <- RMSE(unname(para.ERS[(I+1):(3*I),2]),alpha2.raw),
    RMSE.beta1.ERS <- RMSE(unname(para.ERS[1:I,6])*sign.ERS.beta1,delta1.raw),
    RMSE.beta2.ERS <- RMSE(unname(para.ERS[(I+1):(3*I),7])*sign.ERS.beta2,delta2.raw),
    RMSE.cor.ERS <- RMSE(abs(coef(mod.ERS,simplify=T)$cov[1,2]),R),
    
    RMSE.alpha1.ORD_1 <- RMSE(unname(para.ORD_1[1:I,1]),alpha1.raw),
    RMSE.alpha2.ORD_1 <- RMSE(unname(para.ORD_1[(I+1):(3*I),1]),alpha2.raw),
    RMSE.beta1.ORD_1 <- RMSE(unname(para.ORD_1[1:I,5])*sign.ORD_1.beta1,delta1.raw),
    RMSE.beta2.ORD_1 <- RMSE(unname(para.ORD_1[(I+1):(3*I),5])*sign.ORD_1.beta2,delta2.raw),
    
    RMSE.alpha1.ORD_2 <- RMSE(unname(para.ORD_2[1:I,1]),alpha1.raw),
    RMSE.alpha2.ORD_2 <- RMSE(unname(para.ORD_2[(I+1):(3*I),2]),alpha2.raw),
    RMSE.beta1.ORD_2 <- RMSE(unname(para.ORD_2[1:I,6])*sign.ORD_2.beta1,delta1.raw),
    RMSE.beta2.ORD_2 <- RMSE(unname(para.ORD_2[(I+1):(3*I),7])*sign.ORD_2.beta2,delta2.raw),
    RMSE.cor.ORD_2 <- RMSE(abs(coef(mod.ORD_2,simplify=T)$cov[1,2]),R),
    
    RMSE.alpha1.ORDUTree_1 <- RMSE(unname(para.ORDUTree_1[1:I,1]),alpha1.raw),
    RMSE.alpha2.ORDUTree_1 <- RMSE(unname(para.ORDUTree_1[(I+1):(3*I),1]),alpha2.raw),
    RMSE.delta1.ORDUTree_1 <- RMSE(unname(para.ORDUTree_1[1:I,2])*sign.ORDUTree_1.delta1,delta1.raw),
    RMSE.delta2.ORDUTree_1 <- RMSE(unname(para.ORDUTree_1[(I+1):(3*I),2])*sign.ORDUTree_1.delta2,delta2.raw),
    RMSE.tau1.ORDUTree_1 <- RMSE(unname(para.ORDUTree_1[1:I,3]),tau1.raw),
    RMSE.tau2.ORDUTree_1 <- RMSE(unname(para.ORDUTree_1[(I+1):(3*I),3]),tau2.raw),    
    
   
    AIC.ERS <- ifelse(AIC.mod.ERS==min(AIC.mod.ORDUTree_1,
                                       AIC.mod.ORD_1,AIC.mod.ORD_2,AIC.mod.ERS),1,0),
    AIC.ORD_1 <- ifelse(AIC.mod.ORD_1==min(AIC.mod.ORDUTree_1,
                                           AIC.mod.ORD_1,AIC.mod.ORD_2,AIC.mod.ERS),1,0),
    AIC.ORD_2 <- ifelse(AIC.mod.ORD_2==min(AIC.mod.ORDUTree_1,
                                           AIC.mod.ORD_1,AIC.mod.ORD_2,AIC.mod.ERS),1,0),
    AIC.ORDUTree_1 <- ifelse(AIC.mod.ORDUTree_1==min(AIC.mod.ORDUTree_1,
                                                     AIC.mod.ORD_1,AIC.mod.ORD_2,AIC.mod.ERS),1,0),

    BIC.ERS <- ifelse(BIC.mod.ERS==min(BIC.mod.ORDUTree_1,
                                       BIC.mod.ORD_1,BIC.mod.ORD_2,BIC.mod.ERS),1,0),
    BIC.ORD_1 <- ifelse(BIC.mod.ORD_1==min(BIC.mod.ORDUTree_1,
                                           BIC.mod.ORD_1,BIC.mod.ORD_2,BIC.mod.ERS),1,0),
    BIC.ORD_2 <- ifelse(BIC.mod.ORD_2==min(BIC.mod.ORDUTree_1,
                                           BIC.mod.ORD_1,BIC.mod.ORD_2,BIC.mod.ERS),1,0),
    BIC.ORDUTree_1 <- ifelse(BIC.mod.ORDUTree_1==min(BIC.mod.ORDUTree_1,
                                                     BIC.mod.ORD_1,BIC.mod.ORD_2,BIC.mod.ERS),1,0),
    
    SABIC.ERS <- ifelse(SABIC.mod.ERS==min(SABIC.mod.ORDUTree_1,
                                           SABIC.mod.ORD_1,SABIC.mod.ORD_2,SABIC.mod.ERS),1,0),
    SABIC.ORD_1 <- ifelse(SABIC.mod.ORD_1==min(SABIC.mod.ORDUTree_1,
                                               SABIC.mod.ORD_1,SABIC.mod.ORD_2,SABIC.mod.ERS),1,0),
    SABIC.ORD_2 <- ifelse(SABIC.mod.ORD_2==min(SABIC.mod.ORDUTree_1,
                                               SABIC.mod.ORD_1,SABIC.mod.ORD_2,SABIC.mod.ERS),1,0),
    SABIC.ORDUTree_1 <- ifelse(SABIC.mod.ORDUTree_1==min(SABIC.mod.ORDUTree_1,
                                                         SABIC.mod.ORD_1,SABIC.mod.ORD_2,SABIC.mod.ERS),1,0),

    AIC.mod.ERS,
    AIC.mod.ORD_1,
    AIC.mod.ORD_2,
    AIC.mod.ORDUTree_1,
    BIC.mod.ERS,
    BIC.mod.ORD_1,
    BIC.mod.ORD_2,
    BIC.mod.ORDUTree_1,
    SABIC.mod.ERS,
    SABIC.mod.ORD_1,
    SABIC.mod.ORD_2,
    SABIC.mod.ORDUTree_1,

    # correlation of estimate values and true values
    theta.cor.ERS <- cor(theta.combine.ERS)[2,1],
    eta.cor.ERS <- cor(eta.combine.ERS)[2,1],
    theta.cor.ORD_1 <- cor(theta.combine.ORD_1)[2,1],
    theta1.cor.ORD_2 <- cor(theta1.combine.ORD_2)[2,1],
    theta2.cor.ORD_2 <- cor(theta2.combine.ORD_2)[2,1],
    
    theta.cor.ORDUTree_1 <- cor(theta.combine.ORDUTree_1)[2,1],

    # 95% confidence interval (theta and eta) include true value
    theta_95CItrue.ERS <- mean(ifelse(data.table::between(theta1.raw,theta_95CI.ERS[,1],theta_95CI.ERS[,2]),1,0)),
    eta_95CItrue.ERS <- mean(ifelse(data.table::between(theta1.raw,eta_95CI.ERS[,1],eta_95CI.ERS[,2]),1,0)),
    theta_95CItrue.ORD_1 <- mean(ifelse(data.table::between(theta1.raw,theta_95CI.ORD_1[,1],theta_95CI.ORD_1[,2]),1,0)),
    theta1_95CItrue.ORD_2 <- mean(ifelse(data.table::between(theta1.raw,theta1_95CI.ORD_2[,1],theta1_95CI.ORD_2[,2]),1,0)),
    theta2_95CItrue.ORD_2 <- mean(ifelse(data.table::between(theta1.raw,theta2_95CI.ORD_2[,1],theta2_95CI.ORD_2[,2]),1,0)),
    theta_95CItrue.ORDUTree_1 <- mean(ifelse(data.table::between(theta1.raw,theta_95CI.ORDUTree_1[,1],theta_95CI.ORDUTree_1[,2]),1,0)),

    
    ####################################
    ###### reliablity #####
    ####################################  
    # reliability 1 (Jin et al., 2021)
    theta.relia_Jin.ERS <- 1-(var(theta.esti.ERS-theta1.raw)/var(theta1.raw)),
    eta.relia_Jin.ERS <-1-(var(eta.esti.ERS-theta1.raw)/var(theta1.raw)),
    theta.relia_Jin.ORD_1 <- 1-(var(theta.esti.ORD_1-theta1.raw)/var(theta1.raw)),
    theta1.relia_Jin.ORD_2 <- 1-(var(theta1.esti.ORD_2-theta1.raw)/var(theta1.raw)),
    theta2.relia_Jin.ORD_2 <-1-(var(theta2.esti.ORD_2-theta1.raw)/var(theta1.raw)),
    theta.relia_Jin.ORDUTree_1 <- 1-(var(theta.esti.ORDUTree_1-theta1.raw)/var(theta1.raw)),
   
    # reliablity 2 (mirt)
    theta.relia_mirt.ERS <- empirical_rxx(theta.ERS)[[1]],
    eta.relia_mirt.ERS <- empirical_rxx(theta.ERS)[[2]],
    theta.relia_mirt.ORD_1 <- empirical_rxx(theta.ORD_1)[[1]],
    theta1.relia_mirt.ORD_2 <- empirical_rxx(theta.ORD_2)[[1]],
    theta2.relia_mirt.ORD_2 <- empirical_rxx(theta.ORD_2)[[2]],
    theta.relia_mirt.ORDUTree_1 <- empirical_rxx(theta.ORDUTree_1)[[1]],
    
    # 计算标准误差的均值
    theta_SE_mean.ERS <- mean(theta.ERS[, 3]),
    eta_SE_mean.ERS <- mean(theta.ERS[, 4]),
    theta_SE_mean.ORD_1 <- mean(theta.ORD_1[, 2]),
    theta1_SE_mean.ORD_2 <- mean(theta.ORD_2[, 3]),
    theta2_SE_mean.ORD_2 <- mean(theta.ORD_2[, 4]),
    theta_SE_mean.ORDUTree_1 <- mean(theta.ORDUTree_1[, 2]),
    
    # 计算绝对偏差（Abs）
    Abs.theta.ERS <- mean(abs(theta.combine.ERS[,1] - theta.combine.ERS[,2])),
    Abs.eta.ERS <- mean(abs(eta.combine.ERS[,1] - eta.combine.ERS[,2])),
    Abs.theta.ORD_1 <- mean(abs(theta.combine.ORD_1[,1] - theta.combine.ORD_1[,2])),
    Abs.theta1.ORD_2 <- mean(abs(theta1.combine.ORD_2[,1] - theta1.combine.ORD_2[,2])),
    Abs.theta2.ORD_2 <- mean(abs(theta2.combine.ORD_2[,1] - theta2.combine.ORD_2[,2])),
    Abs.theta.ORDUTree_1 <- mean(abs(theta.combine.ORDUTree_1[,1] - theta.combine.ORDUTree_1[,2])),
    
    # 计算偏差（Bias）
    Bias.theta.ERS <- mean(theta.combine.ERS[,1] - theta.combine.ERS[,2]),
    Bias.eta.ERS <- mean(eta.combine.ERS[,1] - eta.combine.ERS[,2]),
    Bias.theta.ORD_1 <- mean(theta.combine.ORD_1[,1] - theta.combine.ORD_1[,2]),
    Bias.theta1.ORD_2 <- mean(theta1.combine.ORD_2[,1] - theta1.combine.ORD_2[,2]),
    Bias.theta2.ORD_2 <- mean(theta2.combine.ORD_2[,1] - theta2.combine.ORD_2[,2]),
    Bias.theta.ORDUTree_1 <- mean(theta.combine.ORDUTree_1[,1] - theta.combine.ORDUTree_1[,2]),
    
    # 计算均方根误差（RMSE）
    RMSE.theta.ERS <- RMSE(unname(theta.combine.ERS[,1]), theta.combine.ERS[,2]),
    RMSE.eta.ERS <- RMSE(unname(eta.combine.ERS[,1]), eta.combine.ERS[,2]),
    RMSE.theta.ORD_1 <- RMSE(unname(theta.combine.ORD_1[,1]), theta.combine.ORD_1[,2]),
    RMSE.theta1.ORD_2 <- RMSE(unname(theta1.combine.ORD_2[,1]), theta1.combine.ORD_2[,2]),
    RMSE.theta2.ORD_2 <- RMSE(unname(theta2.combine.ORD_2[,1]), theta2.combine.ORD_2[,2]),
    RMSE.theta.ORDUTree_1 <- RMSE(unname(theta.combine.ORDUTree_1[,1]), theta.combine.ORDUTree_1[,2])
    

  )
  
  
  return(ret)
  
  }


## 4. Summarize results

Summarise <- function(condition, results, fixed_objects = NULL){
  
  ret<-colMeans(results,na.rm = T)
  ret
  
}




## 5. Run simulation 
Sim.ORDUTree_1data <- runSimulation(design = Design,
                                 replications =100,
                                 parallel = T ,
                                 generate = Generate.sample,
                                 analyse = Analyse,
                                 summarise = Summarise,
                                 #notification = "complete", #have revised, the original code is "warnings_as_errors = T"
                                 CI=T,
                                 save_results = T,
                                 progress = T,
                                 ncores = 100,
                                 max_errors = 20,
                                 verbose = T,
                                 packages=c("mirt","MASS","truncnorm","GGUM"))



colnames(Sim.ORDUTree_1data) <- c(colnames(Sim.ORDUTree_1data)[1:4],
                                  "Abis.alpha1.ERS","Abis.alpha2.ERS","Abis.beta1.ERS","Abis.beta2.ERS","Abis.cor.ERS",
                                  "Abis.alpha1.ORD_1","Abis.alpha2.ORD_1","Abis.beta1.ORD_1","Abis.beta2.ORD_1",
                                  "Abis.alpha1.ORD_2","Abis.alpha2.ORD_2","Abis.beta1.ORD_2","Abis.beta2.ORD_2","Abis.cor.ORD_2",
                                  "Abis.alpha1.ORDUTree_1","Abis.alpha2.ORDUTree_1","Abis.delta1.ORDUTree_1","Abis.delta2.ORDUTree_1","Abis.tau1.ORDUTree_1","Abis.tau2.ORDUTree_1",
                                  
                                  "Bias.alpha1.ERS","Bias.alpha2.ERS","Bias.beta1.ERS","Bias.beta2.ERS","Bias.cor.ERS",
                                  "Bias.alpha1.ORD_1","Bias.alpha2.ORD_1","Bias.beta1.ORD_1","Bias.beta2.ORD_1",
                                  "Bias.alpha1.ORD_2","Bias.alpha2.ORD_2","Bias.beta1.ORD_2","Bias.beta2.ORD_2","Bias.cor.ORD_2",
                                  "Bias.alpha1.ORDUTree_1","Bias.alpha2.ORDUTree_1","Bias.delta1.ORDUTree_1","Bias.delta2.ORDUTree_1","Bias.tau1.ORDUTree_1","Bias.tau2.ORDUTree_1",
                                  
                                  "RMSE.alpha1.ERS","RMSE.alpha2.ERS","RMSE.beta1.ERS","RMSE.beta2.ERS","RMSE.cor.ERS",
                                  "RMSE.alpha1.ORD_1","RMSE.alpha2.ORD_1","RMSE.beta1.ORD_1","RMSE.beta2.ORD_1",
                                  "RMSE.alpha1.ORD_2","RMSE.alpha2.ORD_2","RMSE.beta1.ORD_2","RMSE.beta2.ORD_2","RMSE.cor.ORD_2",
                                   "RMSE.alpha1.ORDUTree_1","RMSE.alpha2.ORDUTree_1","RMSE.delta1.ORDUTree_1","RMSE.delta2.ORDUTree_1","RMSE.tau1.ORDUTree_1","RMSE.tau2.ORDUTree_1",
                                  
                                  "AIC.ERS","AIC.ORD_1","AIC.ORD_2","AIC.ORDUTree_1",
                                  "BIC.ERS","BIC.ORD_1","BIC.ORD_2","BIC.ORDUTree_1",
                                  "SABIC.ERS","SABIC.ORD_1","SABIC.ORD_2","SABIC.ORDUTree_1",
                                  
                                  "AIC.ERS_origin","AIC.ORD_1_origin","AIC.ORD_2_origin","AIC.ORDUTree_1_origin",
                                  "BIC.ERS_origin","BIC.ORD_1_origin","BIC.ORD_2_origin","BIC.ORDUTree_1_origin",
                                  "SABIC.ERS_origin","SABIC.ORD_1_origin","SABIC.ORD_2_origin","SABIC.ORDUTree_1_origin",
                                  
                                  
                                  "theta.cor.ERS","eta.cor.ERS","theta.cor.ORD_1","theta1.cor.ORD_2","theta2.cor.ORD_2",
                                  "theta.cor.ORDUTree_1",
                                  "theta_95CItrue.ERS","eta_95CItrue.ERS","theta_95CItrue.ORD_1","theta1_95CItrue.ORD_2","theta2_95CItrue.ORD_2",
                                  "theta_95CItrue.ORDUTree_1",
                                  "theta.relia_Jin.ERS","eta.relia_Jin.ERS","theta.relia_Jin.ORD_1","theta1.relia_Jin.ORD_2","theta2.relia_Jin.ORD_2",
                                  "theta.relia_Jin.ORDUTree_1",
                                  "theta.relia_mirt.ERS","eta.relia_mirt.ERS","theta.relia_mirt.ORD_1","theta1.relia_mirt.ORD_2","theta2.relia_mirt.ORD_2",
                                  "theta.relia_mirt.ORDUTree_1",
                                  
                                  "theta_SE_mean.ERS", "eta_SE_mean.ERS", "theta_SE_mean.ORD_1", "theta1_SE_mean.ORD_2", "theta2_SE_mean.ORD_2", "theta_SE_mean.ORDUTree_1",
                                  "Abs.theta.ERS", "Abs.eta.ERS", "Abs.theta.ORD_1", "Abs.theta1.ORD_2", "Abs.theta2.ORD_2", "Abs.theta.ORDUTree_1",
                                  "Bias.theta.ERS", "Bias.eta.ERS", "Bias.theta.ORD_1", "Bias.theta1.ORD_2", "Bias.theta2.ORD_2", "Bias.theta.ORDUTree_1",
                                  "RMSE.theta.ERS", "RMSE.eta.ERS", "RMSE.theta.ORD_1", "RMSE.theta1.ORD_2", "RMSE.theta2.ORD_2", "RMSE.theta.ORDUTree_1",
                                  
                                  "REPLICATIONS","SIM_TIME","COMPLETED","SEED","ERRORS")



write.table(Sim.ORDUTree_1data,"sim_4models_ORDUTree_1data_final_rep100_20240802(delete ORDUTree_2 ERSUTree).csv",row.names = F,sep = ",")
saveRDS(Sim.ORDUTree_1data,"Sim_4models_ORDUTree_1data_final_rep100_20240802(delete ORDUTree_2 ERSUTree).rds")
save.image("Sim_4models_ORDUTree_1data_final_rep100_20240802(delete ORDUTree_2 ERSUTree).RData")



