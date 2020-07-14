

library(MASS)
library(dplyr)
library(parallel)
library(dfoptim)



##############################################
##============================================
## Generate sample from multivaraite normal
##============================================
##############################################
simu.weibull <- function(n1, 
                         ## number of observation from disease cat 1
                         n2, 
                         ## number of observation from disease cat 2
                         n3, 
                         ## number of observation from disease cat 3
                         k,
                         ## vector of shape parameters for d biomarkers
                         lambda, 
                         # vector of scale parameters for 3 disease categories
                         d,
                         seed){
  
  set.seed(seed)
  
  sim1 <- matrix(0,nrow=n1,ncol=d)
  colname <- paste0("X", 1:d)
  colnames(sim1) <- colname
  sim1 <- as.data.frame(sim1)
  
  sim2 <- matrix(0,nrow=n2,ncol=d)
  colnames(sim2) <- colname
  sim2 <- as.data.frame(sim2)
  
  sim3 <- matrix(0,nrow=n3,ncol=d)
  colnames(sim3) <- colname
  sim3 <- as.data.frame(sim3)
  
  for(i in 1:d){
    var <- paste0("X", i)
    sim1[, var] <- rweibull(n1, shape=k[i], scale=lambda[1])
  }
  
  for(i in 1:d){
    var <- paste0("X", i)
    sim2[, var] <- rweibull(n1, shape=k[i], scale=lambda[2])
  }
  
  for(i in 1:d){
    var <- paste0("X", i)
    sim3[, var] <- rweibull(n1, shape=k[i], scale=lambda[3])
  }
  
  sim1$Y <- rep(0, n1)
  sim2$Y <- rep(1, n2)
  sim3$Y <- rep(2, n3)
  
  sim <- rbind(sim1, sim2, sim3)
  
  sim$Y <- as.factor(sim$Y)
  
  sim$Xmax <- as.numeric(apply(sim[,colname], 1, max))
  sim$Xmin <- as.numeric(apply(sim[,colname], 1, min))
  
  return(sim)
}




#######################################
##=====================================
## 1. EHUM
##=====================================
#######################################

EHUM <- function(theta, sim, d){
  theta <- as.numeric(theta)
  beta <- c(theta, 1)
  # beta <- as.numeric(beta)
  
  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)
  
  n1 <- nrow(sim1)
  n2 <- nrow(sim2)
  n3 <- nrow(sim3)
  
  var <- paste0("X", 1:d)
  
  l1 <- t(beta)%*%t(sim1[, var])
  l2 <- t(beta)%*%t(sim2[, var])
  l3 <- t(beta)%*%t(sim3[, var])
  
  S <- 0
  for(j in 1:n2){
    S <- S+ sum(l3>l2[j])*sum(l2[j]>l1)
  }
  
  return(-S/(n1*n2*n3))
}




SHUM <- function(theta, sim, d){
  theta <- as.numeric(theta)
  beta <- c(theta, 1)
  # beta <- as.numeric(beta)
  
  sn <- function(x){
    return(1/(1+exp(-10*x)))
  }
  
  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)
  
  n1 <- nrow(sim1)
  n2 <- nrow(sim2)
  n3 <- nrow(sim3)
  
  var <- paste0("X", 1:d)
  
  l1 <- t(beta)%*%t(sim1[, var])
  l2 <- t(beta)%*%t(sim2[, var])
  l3 <- t(beta)%*%t(sim3[, var])
  
  S <- 0
  for(j in 1:n2){
    sn_i <- sn(l3-l2[j])
    sn_k <- sn(l2[j]-l1)
    
    # sn_i <- sapply(1:n1, function(i) sn(l1[i], l2[j]))
    # sn_k <- sapply(1:n3, function(k) sn(l2[j], l3[k]))
    
    S <- S+ sum(sn_i)*sum(sn_k)
  }
  
  return(-S/(n1*n2*n3))
}




NHUM <- function(theta, sim, d){
  theta <- as.numeric(theta)
  beta <- c(theta, 1)
  # beta <- as.numeric(beta)
  
  sn <- function(x){
    return(pnorm(10*x))
  }
  
  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)
  
  n1 <- nrow(sim1)
  n2 <- nrow(sim2)
  n3 <- nrow(sim3)
  
  var <- paste0("X", 1:d)
  
  l1 <- t(beta)%*%t(sim1[, var])
  l2 <- t(beta)%*%t(sim2[, var])
  l3 <- t(beta)%*%t(sim3[, var])
  
  S <- 0
  for(j in 1:n2){
    sn_i <- sn(l3-l2[j])
    sn_k <- sn(l2[j]-l1)
    
    # sn_i <- sapply(1:n1, function(i) sn(l1[i], l2[j]))
    # sn_k <- sapply(1:n3, function(k) sn(l2[j], l3[k]))
    
    S <- S+ sum(sn_i)*sum(sn_k)
  }
  
  return(-S/(n1*n2*n3))
}



MM <- function(theta, sim){
  theta <- as.numeric(theta)
  beta <- theta
  # beta <- as.numeric(beta)
  
  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)
  
  n1 <- nrow(sim1)
  n2 <- nrow(sim2)
  n3 <- nrow(sim3)
  
  var <- c("Xmax", "Xmin")
  
  l1 <- t(beta)%*%t(sim1[, var])
  l2 <- t(beta)%*%t(sim2[, var])
  l3 <- t(beta)%*%t(sim3[, var])
  
  S <- 0
  for(j in 1:n2){
    S <- S+ sum(l3>l2[j])*sum(l2[j]>l1)
  }
  
  return(-S/(n1*n2*n3))
}



ULBA <- function(theta, sim, d){
  theta <- as.numeric(theta)
  beta <- c(theta, 1)
  # beta <- as.numeric(beta)
  
  S <- rep(0,2)
  for(i in 1:2){
    sim1 <- sim %>% filter(Y==(i-1))
    sim2 <- sim %>% filter(Y==i)
    
    n1 <- nrow(sim1)
    n2 <- nrow(sim2)
    
    var <- paste0("X", 1:d)
    
    l1 <- t(beta)%*%t(sim1[, var])
    l2 <- t(beta)%*%t(sim2[, var])
    
    for(j in 1:n2){
      S[i] <- S[i] + sum(l2[j]>l1)
    }
    
    S[i] <- S[i]/(n1*n2)
  }
  return(-min(S))
}






##############################################
##============================================
## Specify parameters' values and process the data
##============================================
##############################################

d <- 3

k <- c(0.5, 1, 1.5) ## for 3 biomarkers, c(0.5, 1, 2, 5) for 4 biomarkers
lambda <- c(1, 2, 3) ## 3 disease categories
 
lambda.n <- 5      #sqrt(n1+n2+n3)

## computing true beta
n1 <- 1000
n2 <- 1000
n3 <- 1000
seed <- 123
sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)

theta0 <- rep(5, d-1)
opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
theta.true <- opt$par
hum.true <-  -EHUM(theta.true, sim, d)
beta <-  c(theta.true, 1)
beta

MCsz <- 1000



########################################
sim1 <- sim %>% filter(Y==0)
sim2 <- sim %>% filter(Y==1)
sim3 <- sim %>% filter(Y==2)

# delta <- apply(sim2[,var], 2, mean) - apply(sim1[,var], 2, mean) 
# Sigma <- cov(sim[,var])
# beta <- solve(Sigma)%*%delta
# beta <- beta/beta[d]

var <- paste0("X", 1:d)
l1 <- t(beta)%*%t(sim1[, var])
l2 <- t(beta)%*%t(sim2[, var])
l3 <- t(beta)%*%t(sim3[, var])

postscript("density1.pdf", height=4, width=6)
plot(density(as.numeric(l2)), lwd=2, col="red", main="Density plot", xlab="Tumor size", 
     ylab="Density")
lines(density(as.numeric(l1)), lwd=2, col="blue")
dev.off()



#########################################
#########################################
## (n1, n2, n3) = (60, 60, 60)
#########################################
#########################################
n1 <- 60
n2 <- 60
n3 <- 60




##############################################
##============================================
## 1. EHUM approach
##============================================
##############################################
beta.ehum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.ehum.bf <- rep(0, MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  EHUM(beta[1:2], sim, d)
  ## optimize EHUM w.r.t theta
  
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
  theta.ehum.bf <- opt$par
  
  hum.ehum.bf[t] <-  -EHUM(theta.ehum.bf, sim, d)
  beta.ehum.bf[,t] <-  c(theta.ehum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.ehum.bf)
sd(hum.ehum.bf)

apply(beta.ehum.bf, 1, mean)
apply(beta.ehum.bf, 1, mean) - beta
apply(beta.ehum.bf, 1, sd)





#######################################
##=====================================
## 2. SHUM
##=====================================
#######################################

beta.shum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.shum.bf <- rep(0, MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  
  ## optimize shum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, SHUM, sim=sim, d=d, method="BFGS")
  theta.shum.bf <- opt$par
  
  hum.shum.bf[t] <-  -EHUM(theta.shum.bf, sim, d)
  beta.shum.bf[,t] <-  c(theta.shum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.shum.bf)
sd(hum.shum.bf)

apply(beta.shum.bf, 1, mean)
apply(beta.shum.bf, 1, mean) - beta
apply(beta.shum.bf, 1, sd)





#######################################
##=====================================
## 3. NHUM
##=====================================
#######################################
beta.nhum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.nhum.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  
  ## optimize nhum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, NHUM, sim=sim, d=d, method="BFGS")
  theta.nhum.bf <- opt$par
  
  hum.nhum.bf[t] <-  -EHUM(theta.nhum.bf, sim, d)
  beta.nhum.bf[,t] <-  c(theta.nhum.bf, 1)
 }
print(proc.time()-ptm)

mean(hum.nhum.bf)
sd(hum.nhum.bf)

apply(beta.nhum.bf, 1, mean)
apply(beta.nhum.bf, 1, sd)




#######################################
##=====================================
## 4. Min-Max method
##=====================================
#######################################
beta.mm.bf <- matrix(0, ncol=MCsz, nrow=2)
hum.mm.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
 
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize mm w.r.t theta
  theta0 <- rep(2, 2)
  opt <- optim(par=theta0, MM, sim=sim, method="BFGS")
  theta.mm.bf <- opt$par
  
  hum.mm.bf[t] <-  -MM(theta.mm.bf, sim)
  beta.mm.bf[,t] <-  theta.mm.bf
}
print(proc.time()-ptm)

mean(hum.mm.bf)
sd(hum.mm.bf)

apply(beta.mm.bf, 1, mean)
apply(beta.mm.bf, 1, sd)





# 
# #######################################
# ##=====================================
# ## 5. PMNA
# ##=====================================
# #######################################
# PMNA <- function(theta, sim, d){
#   theta <- as.numeric(theta)
#   beta <- c(theta, 1)
#   
#   sim1 <- sim %>% filter(Y==0)
#   sim2 <- sim %>% filter(Y==1)
#   sim3 <- sim %>% filter(Y==2)
#   
#   n1 <- nrow(sim1)
#   n2 <- nrow(sim2)  
#   n3 <- nrow(sim3)
#   
#   var <- paste0("X", 1:d)
#   ## estimate mu1, mu2, mu3
#   mu1 <- apply(sim1[,var], 2, mean)
#   mu2 <- apply(sim2[,var], 2, mean)
#   mu3 <- apply(sim3[,var], 2, mean)
#   
#   ## estimate Sigma1, Sigma2, Sigma3
#   Sigma1 <- t(as.matrix(sim1[,var]))%*%as.matrix(sim1[,var])/n1 - mu1%*%t(mu1)
#   Sigma2 <- t(as.matrix(sim2[,var]))%*%as.matrix(sim2[,var])/n2 - mu2%*%t(mu2)
#   Sigma3 <- t(as.matrix(sim3[,var]))%*%as.matrix(sim3[,var])/n3 - mu3%*%t(mu3)  
#   
#   b1 <- sqrt(t(beta)%*%Sigma2%*%beta)/ sqrt(t(beta)%*%Sigma1%*%beta)
#   a1 <- t(beta)%*%(mu2-mu1)/sqrt(t(beta)%*%Sigma1%*%beta)
#   
#   b2 <- sqrt(t(beta)%*%Sigma2%*%beta)/ sqrt(t(beta)%*%Sigma3%*%beta)
#   a2 <- t(beta)%*%(mu3-mu2)/sqrt(t(beta)%*%Sigma3%*%beta)
#   
#   z <- rnorm(100000, 0, 1)
#   H <- sapply(1:100000, function(i) pnorm(b1*z[i] + a1) * pnorm(-b2*z[i]+a2))
#   
#   HUM <- mean(H)
#   return(-HUM)
# }
# 
# 
# 
# MCsz <- 10
# beta.pmna.bf <- matrix(0, ncol=MCsz, nrow=d)
# hum.pmna.bf <- rep(0, MCsz)
# 
# 
# ptm <- proc.time()
# for(t in 1:MCsz){
#   seed <- 123+t
#   
#   sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
#   
#   ## optimize pmna w.r.t theta
#   theta0 <- rep(2, d-1)
#   opt <- optim(par=theta0, PMNA, sim=sim, d=d, method="BFGS")
#   theta.pmna.bf <- opt$par
#   
#   hum.pmna.bf[t] <-  -EHUM(theta.pmna.bf, sim, d)
#   beta.pmna.bf[,t] <-  c(theta.pmna.bf, 1)
#   
#   print(t)
# }
# print(proc.time()-ptm)
# 
# mean(hum.pmna.bf)
# sd(hum.pmna.bf)
# 
# apply(beta.pmna.bf, 1, mean)
# apply(beta.pmna.bf, 1, sd)
# 



#######################################
##=====================================
## 6. ULBA
##=====================================
#######################################
beta.ulba.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.ulba.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize ulba w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, ULBA, sim=sim, d=d, method="BFGS")
  theta.ulba.bf <- opt$par
  
  hum.ulba.bf[t] <-  -EHUM(theta.ulba.bf, sim, d)
  beta.ulba.bf[,t] <-  c(theta.ulba.bf, 1)
}
print(proc.time()-ptm)

mean(hum.ulba.bf)
sd(hum.ulba.bf)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)



##====================================
## Results
##====================================
df.hum.60 <- c(mean(hum.ehum.bf), sd(hum.ehum.bf), mean(hum.mm.bf), sd(hum.mm.bf), 
               mean(hum.ulba.bf), sd(hum.ulba.bf), mean(hum.shum.bf), sd(hum.shum.bf),
                 mean(hum.nhum.bf), sd(hum.nhum.bf))

# df.hum.sd.60 <- c(sd(hum.ehum.bf), sd(hum.mm.bf), sd(hum.ulba.bf), 
#                  sd(hum.shum.bf), sd(hum.nhum.bf))



apply(beta.ehum.bf, 1, mean)
apply(beta.ehum.bf, 1, sd)

apply(beta.mm.bf, 1, mean)
apply(beta.mm.bf, 1, sd)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)

apply(beta.shum.bf, 1, mean)
apply(beta.shum.bf, 1, sd)

apply(beta.nhum.bf, 1, mean)
apply(beta.nhum.bf, 1, sd)


df.beta <- data.frame(beta=beta, 
                      EHUM.mean=apply(beta.ehum.bf, 1, mean), 
                      EHUM.bias=apply(beta.ehum.bf, 1, mean)-beta, EHUM.sd=apply(beta.ehum.bf, 1, sd),
                      ULBA.mean=apply(beta.ulba.bf, 1, mean), 
                      ULBA.bias=apply(beta.ulba.bf, 1, mean)-beta, ULBA.sd=apply(beta.ulba.bf, 1, sd),
                      SHUM.mean=apply(beta.shum.bf, 1, mean), 
                      SHUM.bias=apply(beta.shum.bf, 1, mean)-beta, SHUM.sd=apply(beta.shum.bf, 1, sd),
                      NHUM.mean=apply(beta.nhum.bf, 1, mean), 
                      NHUM.bias=apply(beta.nhum.bf, 1, mean)-beta, NHUM.sd=apply(beta.nhum.bf, 1, sd))

df.beta.60 <- round(df.beta, digits = 3)












#############################################
#############################################
## (n1, n2, n3) = (90, 90, 90)
#############################################
#############################################
n1 <- 90
n2 <- 90
n3 <- 90





##############################################
##============================================
## 1. EHUM approach
##============================================
##############################################


beta.ehum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.ehum.bf <- rep(0, MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  EHUM(beta[1:2], sim, d)
  ## optimize EHUM w.r.t theta
  
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
  theta.ehum.bf <- opt$par
  
  hum.ehum.bf[t] <-  -EHUM(theta.ehum.bf, sim, d)
  beta.ehum.bf[,t] <-  c(theta.ehum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.ehum.bf)
sd(hum.ehum.bf)

apply(beta.ehum.bf, 1, mean)
apply(beta.ehum.bf, 1, sd)





#######################################
##=====================================
## 2. SHUM
##=====================================
#######################################

beta.shum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.shum.bf <- rep(0, MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  
  ## optimize shum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, SHUM, sim=sim, d=d, method="BFGS")
  theta.shum.bf <- opt$par
  
  hum.shum.bf[t] <-  -EHUM(theta.shum.bf, sim, d)
  beta.shum.bf[,t] <-  c(theta.shum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.shum.bf)
sd(hum.shum.bf)

apply(beta.shum.bf, 1, mean)
apply(beta.shum.bf, 1, sd)





#######################################
##=====================================
## 3. NHUM
##=====================================
#######################################
beta.nhum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.nhum.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  
  ## optimize nhum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, NHUM, sim=sim, d=d, method="BFGS")
  theta.nhum.bf <- opt$par
  
  hum.nhum.bf[t] <-  -EHUM(theta.nhum.bf, sim, d)
  beta.nhum.bf[,t] <-  c(theta.nhum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.nhum.bf)
sd(hum.nhum.bf)

apply(beta.nhum.bf, 1, mean)
apply(beta.nhum.bf, 1, sd)




#######################################
##=====================================
## 4. Min-Max method
##=====================================
#######################################
beta.mm.bf <- matrix(0, ncol=MCsz, nrow=2)
hum.mm.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize mm w.r.t theta
  theta0 <- rep(2, 2)
  opt <- optim(par=theta0, MM, sim=sim, method="BFGS")
  theta.mm.bf <- opt$par
  
  hum.mm.bf[t] <-  -MM(theta.mm.bf, sim)
  beta.mm.bf[,t] <-  theta.mm.bf
}
print(proc.time()-ptm)

mean(hum.mm.bf)
sd(hum.mm.bf)

apply(beta.mm.bf, 1, mean)
apply(beta.mm.bf, 1, sd)





#######################################
##=====================================
## 6. ULBA
##=====================================
#######################################
beta.ulba.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.ulba.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize ulba w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, ULBA, sim=sim, d=d, method="BFGS")
  theta.ulba.bf <- opt$par
  
  hum.ulba.bf[t] <-  -EHUM(theta.ulba.bf, sim, d)
  beta.ulba.bf[,t] <-  c(theta.ulba.bf, 1)
}
print(proc.time()-ptm)

mean(hum.ulba.bf)
sd(hum.ulba.bf)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)



##====================================
## Results
##====================================
df.hum.90 <- c(mean(hum.ehum.bf), sd(hum.ehum.bf), mean(hum.mm.bf), sd(hum.mm.bf), 
               mean(hum.ulba.bf), sd(hum.ulba.bf), mean(hum.shum.bf), sd(hum.shum.bf),
               mean(hum.nhum.bf), sd(hum.nhum.bf))


apply(beta.ehum.bf, 1, mean)
apply(beta.ehum.bf, 1, sd)

apply(beta.mm.bf, 1, mean)
apply(beta.mm.bf, 1, sd)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)

apply(beta.shum.bf, 1, mean)
apply(beta.shum.bf, 1, sd)

apply(beta.nhum.bf, 1, mean)
apply(beta.nhum.bf, 1, sd)


df.beta <- data.frame(beta=beta, 
                      EHUM.mean=apply(beta.ehum.bf, 1, mean), 
                      EHUM.bias=apply(beta.ehum.bf, 1, mean)-beta, EHUM.sd=apply(beta.ehum.bf, 1, sd),
                      ULBA.mean=apply(beta.ulba.bf, 1, mean), 
                      ULBA.bias=apply(beta.ulba.bf, 1, mean)-beta, ULBA.sd=apply(beta.ulba.bf, 1, sd),
                      SHUM.mean=apply(beta.shum.bf, 1, mean), 
                      SHUM.bias=apply(beta.shum.bf, 1, mean)-beta, SHUM.sd=apply(beta.shum.bf, 1, sd),
                      NHUM.mean=apply(beta.nhum.bf, 1, mean), 
                      NHUM.bias=apply(beta.nhum.bf, 1, mean)-beta, NHUM.sd=apply(beta.nhum.bf, 1, sd))

df.beta.90 <- round(df.beta, digits = 3)







#############################################
#############################################
## (n1, n2, n3) = (120, 120, 120)
#############################################
#############################################
n1 <- 120
n2 <- 120
n3 <- 120



##############################################
##============================================
## 1. EHUM approach
##============================================
##############################################
beta.ehum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.ehum.bf <- rep(0, MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  EHUM(beta[1:2], sim, d)
  ## optimize EHUM w.r.t theta
  
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
  theta.ehum.bf <- opt$par
  
  hum.ehum.bf[t] <-  -EHUM(theta.ehum.bf, sim, d)
  beta.ehum.bf[,t] <-  c(theta.ehum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.ehum.bf)
sd(hum.ehum.bf)

apply(beta.ehum.bf, 1, mean)
apply(beta.ehum.bf, 1, sd)





#######################################
##=====================================
## 2. SHUM
##=====================================
#######################################

beta.shum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.shum.bf <- rep(0, MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  
  ## optimize shum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, SHUM, sim=sim, d=d, method="BFGS")
  theta.shum.bf <- opt$par
  
  hum.shum.bf[t] <-  -EHUM(theta.shum.bf, sim, d)
  beta.shum.bf[,t] <-  c(theta.shum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.shum.bf)
sd(hum.shum.bf)

apply(beta.shum.bf, 1, mean)
apply(beta.shum.bf, 1, sd)





#######################################
##=====================================
## 3. NHUM
##=====================================
#######################################
beta.nhum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.nhum.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  
  ## optimize nhum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, NHUM, sim=sim, d=d, method="BFGS")
  theta.nhum.bf <- opt$par
  
  hum.nhum.bf[t] <-  -EHUM(theta.nhum.bf, sim, d)
  beta.nhum.bf[,t] <-  c(theta.nhum.bf, 1)
}
print(proc.time()-ptm)

mean(hum.nhum.bf)
sd(hum.nhum.bf)

apply(beta.nhum.bf, 1, mean)
apply(beta.nhum.bf, 1, sd)




#######################################
##=====================================
## 4. Min-Max method
##=====================================
#######################################
beta.mm.bf <- matrix(0, ncol=MCsz, nrow=2)
hum.mm.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize mm w.r.t theta
  theta0 <- rep(2, 2)
  opt <- optim(par=theta0, MM, sim=sim, method="BFGS")
  theta.mm.bf <- opt$par
  
  hum.mm.bf[t] <-  -MM(theta.mm.bf, sim)
  beta.mm.bf[,t] <-  theta.mm.bf
}
print(proc.time()-ptm)

mean(hum.mm.bf)
sd(hum.mm.bf)

apply(beta.mm.bf, 1, mean)
apply(beta.mm.bf, 1, sd)





#######################################
##=====================================
## 6. ULBA
##=====================================
#######################################
beta.ulba.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.ulba.bf <- rep(0, MCsz)


ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize ulba w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, ULBA, sim=sim, d=d, method="BFGS")
  theta.ulba.bf <- opt$par
  
  hum.ulba.bf[t] <-  -EHUM(theta.ulba.bf, sim, d)
  beta.ulba.bf[,t] <-  c(theta.ulba.bf, 1)
}
print(proc.time()-ptm)

mean(hum.ulba.bf)
sd(hum.ulba.bf)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)



##====================================
## Results
##====================================
df.hum.120 <- c(mean(hum.ehum.bf), sd(hum.ehum.bf), mean(hum.mm.bf), sd(hum.mm.bf), 
               mean(hum.ulba.bf), sd(hum.ulba.bf), mean(hum.shum.bf), sd(hum.shum.bf),
               mean(hum.nhum.bf), sd(hum.nhum.bf))


apply(beta.ehum.bf, 1, mean)
apply(beta.ehum.bf-beta%*%t(rep(1,MCsz)), 1, sd)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf-beta%*%t(rep(1,MCsz)), 1, sd)

apply(beta.shum.bf, 1, mean)
apply(beta.shum.bf-beta%*%t(rep(1,MCsz)), 1, sd)

apply(beta.nhum.bf, 1, mean)
apply(beta.nhum.bf-beta%*%t(rep(1,MCsz)), 1, sd)

# df.beta <- data.frame(beta=beta, EHUM.mean=apply(beta.ehum.bf, 1, mean), EHUM.sd=apply((beta.ehum.bf-beta%*%t(rep(1,MCsz)))^2, 1, mean),
#                       PMNA.mean=apply(beta.ehum.bf, 1, mean), PMNA.sd=apply((beta.ehum.bf-beta%*%t(rep(1,MCsz)))^2, 1, mean),
#                       ULBA.mean=apply(beta.ulba.bf, 1, mean), ULBA.sd=apply((beta.ulba.bf-beta%*%t(rep(1,MCsz)))^2, 1, mean),
#                       SSHUM.mean=apply(beta.shum.bf, 1, mean), SHUM.sd=apply((beta.shum.bf-beta%*%t(rep(1,MCsz)))^2, 1, mean),
#                       NSHUM.mean=apply(beta.nhum.bf, 1, mean), NHUM.sd=apply((beta.nhum.bf-beta%*%t(rep(1,MCsz)))^2, 1, mean))


df.beta <- data.frame(beta=beta,
                      EHUM.mean=apply(beta.ehum.bf, 1, mean), 
                      EHUM.bias=apply(beta.ehum.bf, 1, mean)-beta, EHUM.sd=apply(beta.ehum.bf, 1, sd),
                      ULBA.mean=apply(beta.ulba.bf, 1, mean), 
                      ULBA.bias=apply(beta.ulba.bf, 1, mean)-beta, ULBA.sd=apply(beta.ulba.bf, 1, sd),
                      SHUM.mean=apply(beta.shum.bf, 1, mean), 
                      SHUM.bias=apply(beta.shum.bf, 1, mean)-beta, SHUM.sd=apply(beta.shum.bf, 1, sd),
                      NHUM.mean=apply(beta.nhum.bf, 1, mean), 
                      NHUM.bias=apply(beta.nhum.bf, 1, mean)-beta, NHUM.sd=apply(beta.nhum.bf, 1, sd))

df.beta.120 <- round(df.beta, digits = 3)
df.beta.120



## combined results
df.hum <- t(round(data.frame(n60=df.hum.60, n90=df.hum.90, n120=df.hum.120), digits=3))
df.hum

df.beta.60
df.beta.90
df.beta.120

filename <- paste0("result_weibull_d3_", 4, ".RData")
save(df.hum.60, df.hum.90, df.hum.120, df.beta.60, df.beta.90, df.beta.120, file=filename)














MCsz <- 100
################################
## PMNA method
################################
###############################
###############################
PMNA <- function(theta, sim, d){
  theta <- as.numeric(theta)
  beta <- c(theta, 1)

  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)

  n1 <- nrow(sim1)
  n2 <- nrow(sim2)
  n3 <- nrow(sim3)

  var <- paste0("X", 1:d)
  ## estimate mu1, mu2, mu3
  mu1 <- apply(sim1[,var], 2, mean)
  mu2 <- apply(sim2[,var], 2, mean)
  mu3 <- apply(sim3[,var], 2, mean)

  ## estimate Sigma1, Sigma2, Sigma3
  Sigma1 <- t(as.matrix(sim1[,var]))%*%as.matrix(sim1[,var])/n1 - mu1%*%t(mu1)
  Sigma2 <- t(as.matrix(sim2[,var]))%*%as.matrix(sim2[,var])/n2 - mu2%*%t(mu2)
  Sigma3 <- t(as.matrix(sim3[,var]))%*%as.matrix(sim3[,var])/n3 - mu3%*%t(mu3)

  b1 <- sqrt(t(beta)%*%Sigma2%*%beta)/ sqrt(t(beta)%*%Sigma1%*%beta)
  a1 <- t(beta)%*%(mu2-mu1)/sqrt(t(beta)%*%Sigma1%*%beta)

  b2 <- sqrt(t(beta)%*%Sigma2%*%beta)/ sqrt(t(beta)%*%Sigma3%*%beta)
  a2 <- t(beta)%*%(mu3-mu2)/sqrt(t(beta)%*%Sigma3%*%beta)

  z <- rnorm(100000, 0, 1)
  H <- sapply(1:100000, function(i) pnorm(b1*z[i] + a1) * pnorm(-b2*z[i]+a2))

  HUM <- mean(H)
  return(-HUM)
}




n1 <- 60
n2 <- 60
n3 <- 60

pmna.hum <- function(n1,n2,n3,k,lambda,d, t){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize pmna w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, PMNA, sim=sim, d=d, method="BFGS")
  theta.pmna.bf <- opt$par
  
  hum.pmna.bf <-  -EHUM(theta.pmna.bf, sim, d)
  beta.pmna.bf <-  c(theta.pmna.bf, 1)
  
  filename <- paste0("pmna_hum_beta_weibull_n60_", t, ".RData")
  save(hum.pmna.bf, beta.pmna.bf, file=filename)
}

ptm <- proc.time()
mclapply(1:MCsz, function(t) pmna.hum(n1,n2,n3,k,lambda,d, t), 
         mc.cores=detectCores())
print(proc.time()-ptm)


## Combine the results
hum.pmna.vec <- rep(0, MCsz)
beta.pmna.mat <- matrix(0, ncol=MCsz, nrow=d)
for(t in 1:MCsz){
  filename <- paste0("pmna_hum_beta_weibull_n60_", t, ".RData")
  load(filename)
  
  hum.pmna.vec[t] <-  hum.pmna.bf
  beta.pmna.mat[,t] <- beta.pmna.bf 
}

mean.hum.60.weib <- mean(hum.pmna.vec)
sd.hum.60.weib <- sd(hum.pmna.vec)

mean.beta.60.weib <- apply(beta.pmna.mat, 1, mean)
bias.beta.60.weib <- apply(beta.pmna.mat, 1, mean)-beta
sd.beta.60.weib <- apply(beta.pmna.mat, 1, sd)



###############################
###############################
n1 <- 90
n2 <- 90
n3 <- 90

pmna.hum <- function(n1,n2,n3,k,lambda,d, t){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize pmna w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, PMNA, sim=sim, d=d, method="BFGS")
  theta.pmna.bf <- opt$par
  
  hum.pmna.bf <-  -EHUM(theta.pmna.bf, sim, d)
  beta.pmna.bf <-  c(theta.pmna.bf, 1)
  
  filename <- paste0("pmna_hum_beta_weibull_n90_", t, ".RData")
  save(hum.pmna.bf, beta.pmna.bf, file=filename)
}

ptm <- proc.time()
mclapply(1:MCsz, function(t) pmna.hum(n1,n2,n3,k,lambda,d,t), 
         mc.cores=detectCores())
print(proc.time()-ptm)


## Combine the results
hum.pmna.vec <- rep(0, MCsz)
beta.pmna.mat <- matrix(0, ncol=MCsz, nrow=d)
for(t in 1:MCsz){
  filename <- paste0("pmna_hum_beta_weibull_n90_", t, ".RData")
  load(filename)
  
  hum.pmna.vec[t] <-  hum.pmna.bf
  beta.pmna.mat[,t] <- beta.pmna.bf 
}
mean.hum.90.weib <- mean(hum.pmna.vec)
sd.hum.90.weib <- sd(hum.pmna.vec)

mean.beta.90.weib <- apply(beta.pmna.mat, 1, mean)
bias.beta.90.weib <- apply(beta.pmna.mat, 1, mean)-beta
sd.beta.90.weib <- apply(beta.pmna.mat, 1, sd)



###############################
###############################
n1 <- 120
n2 <- 120
n3 <- 120

pmna.hum <- function(n1,n2,n3,k,lambda,d, t){
  seed <- 123+t
  
  sim <- simu.weibull(n1,n2,n3,k,lambda,d,seed)
  
  ## optimize pmna w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, PMNA, sim=sim, d=d, method="BFGS")
  theta.pmna.bf <- opt$par
  
  hum.pmna.bf <-  -EHUM(theta.pmna.bf, sim, d)
  beta.pmna.bf <-  c(theta.pmna.bf, 1)
  
  filename <- paste0("pmna_hum_beta_weibull_n120_", t, ".RData")
  save(hum.pmna.bf, beta.pmna.bf, file=filename)
}

ptm <- proc.time()
mclapply(1:MCsz, function(t) pmna.hum(n1,n2,n3,k,lambda,d,t), 
         mc.cores=detectCores())
print(proc.time()-ptm)


## Combine the results
hum.pmna.vec <- rep(0, MCsz)
beta.pmna.mat <- matrix(0, ncol=MCsz, nrow=d)
for(t in 1:MCsz){
  filename <- paste0("pmna_hum_beta_weibull_n120_", t, ".RData")
  load(filename)
  
  hum.pmna.vec[t] <-  hum.pmna.bf
  beta.pmna.mat[,t] <- beta.pmna.bf 
}
mean.hum.120.weib <- mean(hum.pmna.vec)
sd.hum.120.weib <- sd(hum.pmna.vec)

mean.beta.120.weib <- apply(beta.pmna.mat, 1, mean)
bias.beta.120.weib <- apply(beta.pmna.mat, 1, mean)-beta
sd.beta.120.weib <- apply(beta.pmna.mat, 1, sd)

filename <- paste0("result_pmna_weibull_d3", ".RData")
save(mean.hum.60.weib, sd.hum.60.weib, 
     mean.beta.60.weib, bias.beta.60.weib, sd.beta.60.weib, 
     mean.hum.90.weib, sd.hum.90.weib, 
     mean.beta.90.weib, bias.beta.90.weib, sd.beta.90.weib, 
     mean.hum.120.weib, sd.hum.120.weib, 
     mean.beta.120.weib, bias.beta.120.weib, sd.beta.120.weib, 
     file=filename)

df.hum.weib <- round(c(mean.hum.60.weib, sd.hum.60.weib, mean.hum.90.weib, sd.hum.90.weib, 
                mean.hum.120.weib, sd.hum.120.weib), digits=3)
df.hum.weib

df.beta.60.weib <- round(data.frame(mean.beta.60.weib, bias.beta.60.weib, sd.beta.60.weib), digits = 3)
df.beta.90.weib <- round(data.frame(mean.beta.90.weib, bias.beta.90.weib, sd.beta.90.weib), digits = 3)
df.beta.120.weib <- round(data.frame(mean.beta.120.weib, bias.beta.120.weib, sd.beta.120.weib), digits = 3)
df.beta.60.weib
df.beta.90.weib
df.beta.120.weib




