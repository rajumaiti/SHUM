
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
d <- 3
MCsz <- 100

delta <- seq(1, length.out=d, by=0.1)
delta <- sort(delta, decreasing=TRUE)
mu <- rep(0, d)

mu1 <- mu
mu2 <- mu + delta
mu3 <- mu + 2*delta


cor.structure <- "independence"
rho <- 0
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







#########################################
#########################################
## (n1, n2, n3) = (60, 60, 60)
#########################################
#########################################

n1 <- 60
n2 <- 60
n3 <- 60

##============================================
## 1. EHUM approach
##============================================

beta.ehum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.ehum.bf <- rep(0, MCsz)

ehum.method <- function(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, t){
  seed <- 123+t
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  
  ## EHUM
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, EHUM, sim=sim, d=d, method="BFGS")
  theta.ehum.bf <- opt$par
  hum.ehum.bf <-  -EHUM(theta.ehum.bf, sim, d)
  beta.ehum.bf <-  c(theta.ehum.bf, 1)
  
  
  
  
  filename <- paste0("ehum_method_ind_normal_d4_n", n1, "_", t, ".RData")
  save(hum.ehum.bf, beta.ehum.bf, file=filename)
}


mclapply(1:MCsz, function(t) ehum.method(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, t), 
         mc.cores = detectCores()-1)

hum.ehum.vec <- rep(NA, MCsz)
beta.ehum.mat <- matrix(NA, d, MCsz)
for(t in 1:MCsz){
  filename <- paste0("ehum_method_ind_normal_d4_n", n1, "_", t, ".RData")
  load(filename)
  
  hum.ehum.vec[t] <- hum.ehum.bf
  beta.ehum.mat[,t] <- beta.ehum.bf
}



##============================================
## 2. SHUM approach
##============================================

beta.shum.bf <- matrix(0, ncol=MCsz, nrow=d)
hum.shum.bf <- rep(0, MCsz)

shum.method <- function(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, t){
  seed <- 123+t
  sim <- simu.normal(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, seed)
  
  theta0 <- rep(1, d-1)
  opt <- optim(par=theta0, SHUM, sim=sim, d=d, method="BFGS")
  theta.shum.bf <- opt$par
  
  hum.shum.bf <-  -SHUM(theta.shum.bf, sim, d)
  beta.shum.bf <-  c(theta.shum.bf, 1)
  
  filename <- paste0("shum_method_ind_normal_d4_n", n1, "_", t, ".RData")
  save(hum.shum.bf, beta.shum.bf, file=filename)
}


mclapply(1:MCsz, function(t) shum.method(n1, n2, n3, mu1, mu2, mu3, cor.structure, rho, d, t), 
         mc.cores = detectCores()-1)

hum.shum.vec <- rep(NA, MCsz)
beta.shum.mat <- matrix(NA, d, MCsz)
for(t in 1:MCsz){
  filename <- paste0("shum_method_ind_normal_d4_n", n1, "_", t, ".RData")
  load(filename)
  
  hum.shum.vec[t] <- hum.shum.bf
  beta.shum.mat[,t] <- beta.shum.bf
}
# 
# mean.shum.60 <- mean(hum.shum.vec)
# sd.shum.60 <- sd(hum.shum.vec)
# 
# mean.shum.beta.60 <- apply(beta.shum.mat, 1, mean)
# bias.shum.beta.60 <- apply(beta.shum.mat, 1, mean) - beta
# sd.shum.beta.60 <- apply(beta.shum.mat, 1, sd)











#################
## OUTPUT
#################

df.hum <- round(c(mean.hum.60, sd.hum.60, mean.hum.90, sd.hum.90, 
                       mean.hum.120, sd.hum.120), digits=3)
df.hum

df.beta.60 <- round(data.frame(mean.beta.60, bias.beta.60, sd.beta.60), digits = 3)
df.beta.90 <- round(data.frame(mean.beta.90, bias.beta.90, sd.beta.90), digits = 3)
df.beta.120 <- round(data.frame(mean.beta.120, bias.beta.120, sd.beta.120), digits = 3)
df.beta.60
df.beta.90
df.beta.120

