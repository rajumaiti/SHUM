
## Date: 02/11/2018

## Estimate beta by maximizing different approx. of HUM in the restricted space B={beta: beta_d=1}


library(MASS)
library(dplyr)
library(parallel)
library(dfoptim)



##############################################
##============================================
## Generate sample from multivaraite normal
##============================================
##############################################
## function to generate samples from multivariate normal
simu.normal<- function(n1, 
                       ## number of observation from disease cat 1
                       n2, 
                       ## number of observation from disease cat 2
                       n3, 
                       ## number of observation from disease cat 3
                       mu1, 
                       ## mean vector of biomarkers for disease cat 1
                       mu2, 
                       ## mean vector of biomarkers for disease cat 2
                       mu3, 
                       ## mean vector of biomarkers for disease cat 3
                       cor.structure, 
                       ## The correlation structure as a string: "independence","AR-1", or "exchangeable"	 
                       rho, 
                       ## correlation parameter
                       d,
                       ## dimension of the biomarkers
                       seed){
  set.seed(seed)
  
  if(cor.structure=="independence"){
    Sigma <- matrix(1,d,d)
    for(i in 1:d){
      for(j in 1:d){
        if(i!=j)
          Sigma[i, j] <- 0
      }
    }
  }
  else if(cor.structure=="exchangeable"){
    Sigma <- matrix(1,d,d)
    for(i in 1:d){
      for(j in 1:d){
        if(i!=j)
          Sigma[i, j] <- rho
      }
    }
  }
  else if(cor.structure=="AR-1"){
    Sigma <- matrix(1,d,d)
    for(i in 1:d){
      for(j in 1:d){
        if(i!=j)
          Sigma[i, j] <- rho^abs(i-j)
      }
    }
  }
  
  sim1 <- as.data.frame(mvrnorm(n1, mu1, Sigma))
  sim1$Y <- rep(0, n1)
  
  sim2 <- as.data.frame(mvrnorm(n2, mu2, Sigma))
  sim2$Y <- rep(1, n2)
  
  sim3 <- as.data.frame(mvrnorm(n3, mu3, Sigma))
  sim3$Y <- rep(2, n3)
  
  sim <- rbind(sim1, sim2, sim3)
  
  colname <- paste0("X", 1:d)
  colnames(sim) <- c(colname, "Y")
  
  sim$Y <- as.factor(sim$Y)
  
  sim$Xmax <- as.numeric(apply(sim[,colname], 1, max))
  sim$Xmin <- as.numeric(apply(sim[,colname], 1, min))
  
  return(as.data.frame(sim))
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
d <- 4

delta <- seq(1, length.out=d, by=0.1)
delta <- sort(delta, decreasing=TRUE)
mu <- rep(0, d)

mu1 <- mu
mu2 <- mu + delta
mu3 <- mu + 2*delta


cor.structure <- "exchangeable"
rho <- 0.3
lambda.n <- 5      #sqrt(n1+n2+n3)


## compute the true beta
sigma.fn <- function(d, rho, cor.structure){
  if(cor.structure=="independence"){
    Sigma <- matrix(1,d,d)
    for(i in 1:d){
      for(j in 1:d){
        if(i!=j)
          Sigma[i, j] <- 0
      }
    }
  }
  else if(cor.structure=="exchangeable"){
    Sigma <- matrix(1,d,d)
    for(i in 1:d){
      for(j in 1:d){
        if(i!=j)
          Sigma[i, j] <- rho
      }
    }
  }
  else if(cor.structure=="AR-1"){
    Sigma <- matrix(1,d,d)
    for(i in 1:d){
      for(j in 1:d){
        if(i!=j)
          Sigma[i, j] <- rho^abs(i-j)
      }
    }
  }
  return(Sigma)
}

## true beta
beta <- as.numeric(solve(sigma.fn(d, rho, cor.structure))%*%delta)
beta <- beta/beta[d]

sim <- simu.normal(1000, 1000, 1000, mu1, mu2, mu3, cor.structure, rho, d, seed=123)
hum.true.ind <- -EHUM(beta[1:(d-1)], sim, d)
hum.true.ind






MCsz <- 1000
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  EHUM(beta[1:2], sim, d)
  ## optimize EHUM w.r.t theta
  
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
  theta.ehum.bf <- opt$par
  
  hum.ehum.bf[t] <-  -EHUM(theta.ehum.bf, sim, d)
  beta.ehum.bf[,t] <-  c(theta.ehum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  
  ## optimize shum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, SHUM, sim=sim, d=d, method="BFGS")
  theta.shum.bf <- opt$par
  
  hum.shum.bf[t] <-  -EHUM(theta.shum.bf, sim, d)
  beta.shum.bf[,t] <-  c(theta.shum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  
  ## optimize nhum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, NHUM, sim=sim, d=d, method="BFGS")
  theta.nhum.bf <- opt$par
  
  hum.nhum.bf[t] <-  -EHUM(theta.nhum.bf, sim, d)
  beta.nhum.bf[,t] <-  c(theta.nhum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  ## optimize mm w.r.t theta
  theta0 <- rep(2, 2)
  opt <- optim(par=theta0, MM, sim=sim, method="BFGS")
  theta.mm.bf <- opt$par
  
  hum.mm.bf[t] <-  -MM(theta.mm.bf, sim)
  beta.mm.bf[,t] <-  theta.mm.bf
  
  print(t)
}
print(proc.time()-ptm)

mean(hum.mm.bf)
sd(hum.mm.bf)

apply(beta.mm.bf, 1, mean)
apply(beta.mm.bf, 1, sd)






#######################################
##=====================================
## 5. PMNA
##=====================================
#######################################
hum.pmna.vec <- rep(0, MCsz)
beta.pmna.mat <- matrix(0, nrow=d, ncol=MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)
  
  var <- paste0("X", 1:d)
  mu1.hat <- apply(sim1[, var], 2, mean)
  mu2.hat <- apply(sim2[, var], 2, mean)
  mu3.hat <- apply(sim3[, var], 2, mean)
  
  delta1 <- mu2.hat - mu1.hat
  delta2 <- mu3.hat - mu2.hat
  delta.hat <- (delta1+delta2)/2
  
  Sigma.hat <- cov(sim[,var])
  
  theta.pmna.bf <- solve(Sigma.hat)%*%delta.hat
  theta.pmna.bf <- theta.pmna.bf/theta.pmna.bf[d]
  
  hum.pmna.vec[t] <-  -EHUM(theta.pmna.bf[1:(d-1)], sim, d)
  beta.pmna.mat[,t] <-  as.numeric(theta.pmna.bf)
}
print(proc.time()-ptm)

mean(hum.pmna.vec)
sd(hum.pmna.vec)

apply(beta.pmna.mat, 1, mean)
apply(beta.pmna.mat, 1, sd)




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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  ## optimize ulba w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, ULBA, sim=sim, d=d, method="BFGS")
  theta.ulba.bf <- opt$par
  
  hum.ulba.bf[t] <-  -EHUM(theta.ulba.bf, sim, d)
  beta.ulba.bf[,t] <-  c(theta.ulba.bf, 1)
  
  print(t)
}
print(proc.time()-ptm)

mean(hum.ulba.bf)
sd(hum.ulba.bf)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)



##====================================
## Results
##====================================
df.hum.60 <- c(mean(hum.ehum.bf), sd(hum.ehum.bf), mean(hum.pmna.vec), sd(hum.pmna.vec),
               mean(hum.mm.bf), sd(hum.mm.bf), 
               mean(hum.ulba.bf), sd(hum.ulba.bf), mean(hum.shum.bf), sd(hum.shum.bf),
               mean(hum.nhum.bf), sd(hum.nhum.bf))

df.beta <- data.frame(beta=beta, EHUM.mean=apply(beta.ehum.bf, 1, mean), 
                      EHUM.bias=apply(beta.ehum.bf, 1, mean)-beta, EHUM.sd=apply(beta.ehum.bf, 1, sd),
                      PMNA.mean=apply(beta.pmna.mat, 1, mean), 
                      PMNA.bias=apply(beta.pmna.mat, 1, mean)-beta, PMNA.sd=apply(beta.pmna.mat, 1, sd),
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  EHUM(beta[1:2], sim, d)
  ## optimize EHUM w.r.t theta
  
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
  theta.ehum.bf <- opt$par
  
  hum.ehum.bf[t] <-  -EHUM(theta.ehum.bf, sim, d)
  beta.ehum.bf[,t] <-  c(theta.ehum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  
  ## optimize shum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, SHUM, sim=sim, d=d, method="BFGS")
  theta.shum.bf <- opt$par
  
  hum.shum.bf[t] <-  -EHUM(theta.shum.bf, sim, d)
  beta.shum.bf[,t] <-  c(theta.shum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  
  ## optimize nhum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, NHUM, sim=sim, d=d, method="BFGS")
  theta.nhum.bf <- opt$par
  
  hum.nhum.bf[t] <-  -EHUM(theta.nhum.bf, sim, d)
  beta.nhum.bf[,t] <-  c(theta.nhum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  ## optimize mm w.r.t theta
  theta0 <- rep(2, 2)
  opt <- optim(par=theta0, MM, sim=sim, method="BFGS")
  theta.mm.bf <- opt$par
  
  hum.mm.bf[t] <-  -MM(theta.mm.bf, sim)
  beta.mm.bf[,t] <-  theta.mm.bf
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  ## optimize ulba w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, ULBA, sim=sim, d=d, method="BFGS")
  theta.ulba.bf <- opt$par
  
  hum.ulba.bf[t] <-  -EHUM(theta.ulba.bf, sim, d)
  beta.ulba.bf[,t] <-  c(theta.ulba.bf, 1)
  
  print(t)
}
print(proc.time()-ptm)

mean(hum.ulba.bf)
sd(hum.ulba.bf)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)


#######################################
##=====================================
## 5. PMNA
##=====================================
#######################################
hum.pmna.vec <- rep(0, MCsz)
beta.pmna.mat <- matrix(0, nrow=d, ncol=MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)
  
  var <- paste0("X", 1:d)
  mu1.hat <- apply(sim1[, var], 2, mean)
  mu2.hat <- apply(sim2[, var], 2, mean)
  mu3.hat <- apply(sim3[, var], 2, mean)
  
  delta1 <- mu2.hat - mu1.hat
  delta2 <- mu3.hat - mu2.hat
  delta.hat <- (delta1+delta2)/2
  
  Sigma.hat <- cov(sim[,var])
  
  theta.pmna.bf <- solve(Sigma.hat)%*%delta.hat
  theta.pmna.bf <- theta.pmna.bf/theta.pmna.bf[d]
  
  hum.pmna.vec[t] <-  -EHUM(theta.pmna.bf[1:(d-1)], sim, d)
  beta.pmna.mat[,t] <-  as.numeric(theta.pmna.bf)
}
print(proc.time()-ptm)

mean(hum.pmna.vec)
sd(hum.pmna.vec)

apply(beta.pmna.mat, 1, mean)
apply(beta.pmna.mat, 1, sd)



##====================================
## Results
##====================================
df.hum.90 <- c(mean(hum.ehum.bf), sd(hum.ehum.bf), mean(hum.pmna.vec), sd(hum.pmna.vec), 
               mean(hum.mm.bf), sd(hum.mm.bf), 
               mean(hum.ulba.bf), sd(hum.ulba.bf), mean(hum.shum.bf), sd(hum.shum.bf),
               mean(hum.nhum.bf), sd(hum.nhum.bf))




df.beta <- data.frame(beta=beta, EHUM.mean=apply(beta.ehum.bf, 1, mean), 
                      EHUM.bias=apply(beta.ehum.bf, 1, mean)-beta, EHUM.sd=apply(beta.ehum.bf, 1, sd),
                      PMNA.mean=apply(beta.pmna.mat, 1, mean), 
                      PMNA.bias=apply(beta.pmna.mat, 1, mean)-beta, PMNA.sd=apply(beta.pmna.mat, 1, sd),
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  EHUM(beta[1:2], sim, d)
  ## optimize EHUM w.r.t theta
  
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
  theta.ehum.bf <- opt$par
  
  hum.ehum.bf[t] <-  -EHUM(theta.ehum.bf, sim, d)
  beta.ehum.bf[,t] <-  c(theta.ehum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  
  ## optimize shum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, SHUM, sim=sim, d=d, method="BFGS")
  theta.shum.bf <- opt$par
  
  hum.shum.bf[t] <-  -EHUM(theta.shum.bf, sim, d)
  beta.shum.bf[,t] <-  c(theta.shum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  
  ## optimize nhum w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, NHUM, sim=sim, d=d, method="BFGS")
  theta.nhum.bf <- opt$par
  
  hum.nhum.bf[t] <-  -EHUM(theta.nhum.bf, sim, d)
  beta.nhum.bf[,t] <-  c(theta.nhum.bf, 1)
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  ## optimize mm w.r.t theta
  theta0 <- rep(2, 2)
  opt <- optim(par=theta0, MM, sim=sim, method="BFGS")
  theta.mm.bf <- opt$par
  
  hum.mm.bf[t] <-  -MM(theta.mm.bf, sim)
  beta.mm.bf[,t] <-  theta.mm.bf
  
  print(t)
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
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  ## optimize ulba w.r.t theta
  theta0 <- rep(2, d-1)
  opt <- optim(par=theta0, ULBA, sim=sim, d=d, method="BFGS")
  theta.ulba.bf <- opt$par
  
  hum.ulba.bf[t] <-  -EHUM(theta.ulba.bf, sim, d)
  beta.ulba.bf[,t] <-  c(theta.ulba.bf, 1)
  
  print(t)
}
print(proc.time()-ptm)

mean(hum.ulba.bf)
sd(hum.ulba.bf)

apply(beta.ulba.bf, 1, mean)
apply(beta.ulba.bf, 1, sd)


#######################################
##=====================================
## 5. PMNA
##=====================================
#######################################
hum.pmna.vec <- rep(0, MCsz)
beta.pmna.mat <- matrix(0, nrow=d, ncol=MCsz)

ptm <- proc.time()
for(t in 1:MCsz){
  seed <- 123+t
  
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  sim1 <- sim %>% filter(Y==0)
  sim2 <- sim %>% filter(Y==1)
  sim3 <- sim %>% filter(Y==2)
  
  var <- paste0("X", 1:d)
  mu1.hat <- apply(sim1[, var], 2, mean)
  mu2.hat <- apply(sim2[, var], 2, mean)
  mu3.hat <- apply(sim3[, var], 2, mean)
  
  delta1 <- mu2.hat - mu1.hat
  delta2 <- mu3.hat - mu2.hat
  delta.hat <- (delta1+delta2)/2
  
  Sigma.hat <- cov(sim[,var])
  
  theta.pmna.bf <- solve(Sigma.hat)%*%delta.hat
  theta.pmna.bf <- theta.pmna.bf/theta.pmna.bf[d]
  
  hum.pmna.vec[t] <-  -EHUM(theta.pmna.bf[1:(d-1)], sim, d)
  beta.pmna.mat[,t] <-  as.numeric(theta.pmna.bf)
}
print(proc.time()-ptm)

mean(hum.pmna.vec)
sd(hum.pmna.vec)

apply(beta.pmna.mat, 1, mean)
apply(beta.pmna.mat, 1, sd)


##====================================
## Results
##====================================
df.hum.120 <- c(mean(hum.ehum.bf), sd(hum.ehum.bf), mean(hum.pmna.vec), sd(hum.pmna.vec),
                mean(hum.mm.bf), sd(hum.mm.bf), 
                mean(hum.ulba.bf), sd(hum.ulba.bf), mean(hum.shum.bf), sd(hum.shum.bf),
                mean(hum.nhum.bf), sd(hum.nhum.bf))


df.beta <- data.frame(beta=beta, EHUM.mean=apply(beta.ehum.bf, 1, mean), 
                      EHUM.bias=apply(beta.ehum.bf, 1, mean)-beta, EHUM.sd=apply(beta.ehum.bf, 1, sd),
                      PMNA.mean=apply(beta.pmna.mat, 1, mean), 
                      PMNA.bias=apply(beta.pmna.mat, 1, mean)-beta, PMNA.sd=apply(beta.pmna.mat, 1, sd),
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


## save all the above results
filename <- paste0("result_normal_ind_d4", ".RData")
save(df.hum.60, df.hum.90, df.hum.120, df.beta.60, df.beta.90, df.beta.120, file=filename)

