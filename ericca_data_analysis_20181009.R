
## Date: 09 Oct 2018

## ERICCA data analysis






## required R packages
## for data managing
library(dplyr)
library(parallel)


###################################
## Reading data
###################################

## ericca data set
data <- read.csv("errica_ind_ngal_21june2016.csv", header=TRUE, na.strings = "")

## select N0, N6, N12, N24 and AKI
data_tmp <- data %>% select(NGAL.0.hours, NGAL.6.hours, NGAL.12.hours, NGAL.24.hours, AKI)

data_tmp <- data_tmp %>% mutate(X1 = log(NGAL.0.hours), X2 = log(NGAL.6.hours),
                                X3 = log(NGAL.12.hours), X4 = log(NGAL.24.hours), 
                                Y=AKI)

## data set for the analysis
data_mvn <-  data_tmp %>% select(X1, X2, X3, X4, Y)
data_mvn <- data_mvn %>% filter(X1>5.5) 
# data_mvn <- data_mvn %>% filter(X1>=4.75, X1<=5.5) #filter(X1>5.5)

## delete all missing observation
dat <- na.omit(data_mvn)
n <- nrow(dat)

dat_3cat <- dat
dat_3cat$Y <- ifelse(dat_3cat$Y==3, 2, dat_3cat$Y)

dat_3cat$Xmax <- as.numeric(apply(dat_3cat[, c("X1", "X2", "X3", "X4")], 1, max))
dat_3cat$Xmin <- as.numeric(apply(dat_3cat[, c("X1", "X2", "X3", "X4")], 1, min))
dat_3cat$Xmean <- as.numeric(apply(dat_3cat[, c("X1", "X2", "X3", "X4")], 1, mean))

d <- 4


##################################################################
##================================================================
## Section 1: est of HUM (se), and est of beta (se) using different methods
##================================================================
##################################################################


##=========================
## Method 1: Simple average
##=========================
## It is given after Method 7

##================================
## Method 2: EHUM
##================================
EHUM <- function(dat, beta, d){
  beta <- as.numeric(beta)
  
  I <- function(x){
    x <- as.numeric(x)
    return(1*(x[3]>x[2])*(x[2]>x[1]))
  }
  
  dat1 <- dat %>% filter(Y==0)
  dat2 <- dat %>% filter(Y==1)
  dat3 <- dat %>% filter(Y==2)
  
  n1 <- nrow(dat1)
  n2 <- nrow(dat2)  
  n3 <- nrow(dat3)
  
  var <- paste0("X", 1:d)
  
  l1 <- t(beta) %*% t(dat1[, var])
  l2 <- t(beta) %*% t(dat2[, var])
  l3 <- t(beta) %*% t(dat3[, var])
  
  L <- expand.grid(l1=l1, l2=l2, l3=l3)
  
  S <- apply(L, 1, I)
  
  return(-sum(S)/(n1*n2*n3))
}


# ptm <- proc.time()
# beta0 <- rep(0.5, 4)
# opt <- optim(par=beta0, EHUM, dat=dat_3cat, method="L-BFGS-B", lower=-1, upper=1)
# beta.ehum.qn <- opt$par/sqrt(sum(opt$par^2))
# hum.ehum.qn <- -EHUM(dat_3cat, beta.ehum.qn, 4)
# hum.ehum.qn
# print(proc.time()-ptm)


ptm <- proc.time()
beta0 <- rep(0.5, d)
opt <- optim(par=beta0, EHUM, dat=dat_3cat, d=d, method="Nelder-Mead")
beta.ehum.nm <- opt$par/sqrt(sum(opt$par^2))
hum.ehum.nm <- -EHUM(dat_3cat, beta.ehum.nm, d)
hum.ehum.nm
print(proc.time()-ptm)



## compute se of beta and HUM
boot.fn.ehum <- function(dat, d, b){
  set.seed(123+b)
  
  n <- nrow(dat)
  indx <- sample(1:n, n, replace=TRUE)
  boot.dat <- dat[indx,]
  
  beta0 <- rep(0.5, d)
  # opt <- optim(par=beta0, EHUM, dat=boot.dat, method="L-BFGS-B", lower=-1, upper=1)
  opt <- optim(par=beta0, EHUM, dat=boot.dat, d=d, method="Nelder-Mead")
  beta.ehum <- opt$par/sqrt(sum(opt$par^2))
  hum.ehum <- -EHUM(boot.dat, beta.ehum, d)
  
  filename <- paste0("ehum_boot_ericca_", b, ".RData")
  save(hum.ehum, beta.ehum, file=filename)
}


ptm <- proc.time()
B <- 100
mclapply(1:B, function(b) boot.fn.ehum(dat_3cat, d, b), mc.cores = detectCores()-1)
print(proc.time()-ptm)

hum.ehum.vec <- rep(NA, B)
beta.ehum.mat <- matrix(NA, d, B)
for(b in 1:B){
  filename <- paste0("ehum_boot_ericca_", b, ".RData")
  load(filename)
  
  hum.ehum.vec[b] <- hum.ehum
  beta.ehum.mat[,b] <- beta.ehum
}

se.hum.ehum <- sd(hum.ehum.vec)
se.beta.ehum <- apply(beta.ehum.mat, 1, sd)

se.hum.ehum
se.beta.ehum



##=========================
## Method 3: SSHUM
##=========================

SHUM <- function(dat, beta, d){
  
  I_sig <- function(x){
    return(1/(1+exp(-10*(x[2]-x[1]))) * 1/(1+exp(-10*(x[3]-x[2]))) )
  }
  
  beta <- as.numeric(beta)
  
  dat1 <- dat %>% filter(Y==0)
  dat2 <- dat %>% filter(Y==1)
  dat3 <- dat %>% filter(Y==2)
  
  n1 <- nrow(dat1)
  n2 <- nrow(dat2)  
  n3 <- nrow(dat3)
  
  var <- paste0("X", 1:d)
  
  l1 <- t(beta) %*% t(dat1[, var])
  l2 <- t(beta) %*% t(dat2[, var])
  l3 <- t(beta) %*% t(dat3[, var])
  
  L <- expand.grid(l1=l1, l2=l2, l3=l3)
  
  S <- apply(L, 1, I_sig)
  
  return(-sum(S)/(n1*n2*n3))
}


ptm <- proc.time()
beta0 <- rep(0.5, d)
opt <- optim(par=beta0, SHUM, dat=dat_3cat, d=d, method="L-BFGS-B", lower=-1, upper=1)
beta.shum.qn <- opt$par/sqrt(sum(opt$par^2))
hum.shum.qn <- -EHUM(dat_3cat, beta.shum.qn, d)
hum.shum.qn
print(proc.time()-ptm)



boot.fn.shum <- function(dat, d, b){
  set.seed(123+b)
  
  n <- nrow(dat)
  indx <- sample(1:n, n, replace=TRUE)
  boot.dat <- dat[indx,]
  
  beta0 <- rep(0.5, d)
  opt <- optim(par=beta0, SHUM, dat=boot.dat, d=d, method="L-BFGS-B", lower=-1, upper=1)
  beta.shum <- opt$par/sqrt(sum(opt$par^2))
  hum.shum <- -EHUM(boot.dat, beta.shum,d)
  
  filename <- paste0("shum_boot_ericca_", b, ".RData")
  save(hum.shum, beta.shum, file=filename)
}


ptm <- proc.time()
B <- 100
mclapply(1:B, function(b) boot.fn.shum(dat_3cat, d, b), mc.cores = detectCores()-1)
print(proc.time()-ptm)


hum.shum.vec <- rep(NA, B)
beta.shum.mat <- matrix(NA, d, B)
for(b in 1:B){
  filename <- paste0("shum_boot_ericca_", b, ".RData")
  load(filename)
  
  hum.shum.vec[b] <- hum.shum
  beta.shum.mat[,b] <- beta.shum
}

se.hum.shum <- sd(hum.shum.vec)
se.beta.shum <- apply(beta.shum.mat, 1, sd)

se.hum.shum
se.beta.shum





##=========================
## Method 4: NSHUM
##=========================

NHUM <- function(dat, beta, d){
  
  I_norm <- function(x){
    return(pnorm(10*(x[3]-x[2])) * pnorm(10*(x[2]-x[1])))
  }
  
  
  beta <- as.numeric(beta)
  
  dat1 <- dat %>% filter(Y==0)
  dat2 <- dat %>% filter(Y==1)
  dat3 <- dat %>% filter(Y==2)
  
  n1 <- nrow(dat1)
  n2 <- nrow(dat2)  
  n3 <- nrow(dat3)
  
  var <- paste0("X", 1:d)
  
  l1 <- t(beta) %*% t(dat1[, var])
  l2 <- t(beta) %*% t(dat2[, var])
  l3 <- t(beta) %*% t(dat3[, var])
  
  L <- expand.grid(l1=l1, l2=l2, l3=l3)
  
  S <- apply(L, 1, I_norm)
  
  return(-sum(S)/(n1*n2*n3))
}


ptm <- proc.time()
beta0 <- rep(0.5, d)
opt <- optim(par=beta0, NHUM, dat=dat_3cat, d=d, method="L-BFGS-B", lower=-1, upper=1)
beta.nhum.qn <- opt$par/sqrt(sum(opt$par^2))
hum.nhum.qn <- -EHUM(dat_3cat, beta.nhum.qn, d)
hum.nhum.qn
print(proc.time()-ptm)



boot.fn.nhum <- function(dat, d, b){
  set.seed(123+b)
  
  n <- nrow(dat)
  indx <- sample(1:n, n, replace=TRUE)
  boot.dat <- dat[indx,]
  
  beta0 <- rep(0.5, d)
  opt <- optim(par=beta0, NHUM, dat=boot.dat, d=d, method="L-BFGS-B", lower=-1, upper=1)
  beta.nhum <- opt$par/sqrt(sum(opt$par^2))
  hum.nhum <- -EHUM(boot.dat, beta.nhum,d)
  
  filename <- paste0("nhum_boot_ericca_", b, ".RData")
  save(hum.nhum, beta.nhum, file=filename)
}


ptm <- proc.time()
B <- 100
mclapply(1:B, function(b) boot.fn.nhum(dat_3cat, d, b), mc.cores = detectCores()-1)
print(proc.time()-ptm)

hum.nhum.vec <- rep(NA, B)
beta.nhum.mat <- matrix(NA, d, B)
for(b in 1:B){
  filename <- paste0("nhum_boot_ericca_",b,".RData")
  load(filename)
  
  hum.nhum.vec[b] <- hum.nhum
  beta.nhum.mat[,b] <- beta.nhum
}

se.hum.nhum <- sd(hum.nhum.vec)
se.beta.nhum <- apply(beta.nhum.mat, 1, sd)

se.hum.nhum
se.beta.nhum




##=========================
## Method 5: ULBA
##=========================

PM <- function(dat, beta, d) {

I.bin <- function(x){
  x <- as.numeric(x)
  return(1*(x[2]>x[1]))
}

  S <- rep(0,2)
  for(i in 1:2){
    dat1 <- dat %>% filter(Y==(i-1))
    dat2 <- dat %>% filter(Y==i)
    
    n1 <- nrow(dat1)
    n2 <- nrow(dat2)  
    
    var <- paste0("X", 1:d)
  
    l1 <- t(beta) %*% t(dat1[, var])
    l2 <- t(beta) %*% t(dat2[, var])  
     
    L <- expand.grid(l1=l1, l2=l2)
    
    S[i] <- sum(apply(L, 1, I.bin))/(n1*n2)
  }
  
  return(-min(S))
}

# 
# ptm <- proc.time()
# beta0 <- rep(0.5, d)
# opt <- optim(par=beta0, PM, dat=dat_3cat, d=d,  method="L-BFGS-B", lower=-1, upper=1)
# beta.pm.qn <- opt$par/sqrt(sum(opt$par^2))
# hum.pm.qn <- -EHUM(dat_3cat, beta.pm.qn, d)
# hum.pm.qn
# print(proc.time()-ptm)


ptm <- proc.time()
beta0 <- rep(0.5, d)
opt <- optim(par=beta0, PM, dat=dat_3cat, d=d, method="Nelder-Mead")
beta.pm.nm <- opt$par/sqrt(sum(opt$par^2))
hum.pm.nm <- -EHUM(dat_3cat, beta.pm.nm, d)
hum.pm.nm
print(proc.time()-ptm)


boot.fn.pm <- function(dat, d, b){
  set.seed(123+b)
  
  n <- nrow(dat)
  indx <- sample(1:n, n, replace=TRUE)
  boot.dat <- dat[indx,]
  
  beta0 <- rep(0.5, d)
  opt <- optim(par=beta0, PM, dat=boot.dat, d=d, method="L-BFGS-B", lower=-1, upper=1)
  beta.pm <- opt$par/sqrt(sum(opt$par^2))
  hum.pm <- -EHUM(boot.dat, beta.pm, d)
  
  filename <- paste0("pm_boot_ericca_", b, ".RData")
  save(hum.pm, beta.pm, file=filename)
}


ptm <- proc.time()
B <- 100
mclapply(1:B, function(b) boot.fn.pm(dat_3cat, d, b), mc.cores = detectCores()-1)
print(proc.time()-ptm)

hum.pm.vec <- rep(NA, B)
beta.pm.mat <- matrix(NA, d, B)
for(b in 1:B){
  filename <- paste0("pm_boot_ericca_", b, ".RData")
  load(filename)
  
  hum.pm.vec[b] <- hum.pm
  beta.pm.mat[,b] <- beta.pm
}

se.hum.pm <- sd(hum.pm.vec)
se.beta.pm <- apply(beta.pm.mat, 1, sd)

se.hum.pm
se.beta.pm




##=========================
## Method 6: PMNA
##=========================

PMNA <- function(beta, dat, d){
  beta <- as.numeric(beta)
  
  dat1 <- dat %>% filter(Y==0)
  dat2 <- dat %>% filter(Y==1)
  dat3 <- dat %>% filter(Y==2)
  
  n1 <- nrow(dat1)
  n2 <- nrow(dat2)  
  n3 <- nrow(dat3)
  
  var <- paste0("X", 1:d)  
  
  ## estimate mu1, mu2, mu3
  mu1 <- as.numeric(apply(dat1[,var], 2, mean)) #c(mean(dat1$X1), mean(dat1$X2), mean(dat1$X3), mean(dat1$X4))
  mu2 <- as.numeric(apply(dat2[,var], 2, mean)) #c(mean(dat2$X1), mean(dat2$X2), mean(dat2$X3), mean(dat2$X4))
  mu3 <- as.numeric(apply(dat3[,var], 2, mean)) #c(mean(dat3$X1), mean(dat3$X2), mean(dat3$X3), mean(dat3$X4))
  
  ## estimate Sigma1, Sigma2, Sigma3
  Sigma1 <- t(as.matrix(dat1[, var]))%*%as.matrix(dat1[,var])/n1 - mu1%*%t(mu1)
  Sigma2 <- t(as.matrix(dat2[, var]))%*%as.matrix(dat2[,var])/n2 - mu2%*%t(mu2)
  Sigma3 <- t(as.matrix(dat3[, var]))%*%as.matrix(dat3[,var])/n3 - mu3%*%t(mu3)  
  
  b1 <- sqrt(t(beta)%*%Sigma2%*%beta)/ sqrt(t(beta)%*%Sigma1%*%beta)
  a1 <- t(beta)%*%(mu2-mu1)/sqrt(t(beta)%*%Sigma1%*%beta)
  
  b2 <- sqrt(t(beta)%*%Sigma2%*%beta)/ sqrt(t(beta)%*%Sigma3%*%beta)
  a2 <- t(beta)%*%(mu3-mu2)/sqrt(t(beta)%*%Sigma3%*%beta)
  
  z <- rnorm(5000,0,1)
  H <- sapply(1:5000, function(i) pnorm(b1*z[i] + a1) * pnorm(-b2*z[i]+a2))
  
  hum <- mean(H)
  return(-hum)
}



ptm <- proc.time()
beta0 <- rep(0.5, d)
opt <- optim(par=beta0, PMNA, dat=dat_3cat, d=d, method="L-BFGS-B", lower=-1, upper=1)
beta.pmna.qn <- opt$par/sqrt(sum(opt$par^2))
hum.pmna.qn <- -EHUM(dat_3cat, beta.pmna.qn, d)
hum.pmna.qn
print(proc.time()-ptm)

# 
# ptm <- proc.time()
# beta0 <- rep(0.5, d)
# opt <- optim(par=beta0, PMNA, dat=dat_3cat, d=d, method="Nelder-Mead")
# beta.pmna.nm <- opt$par/sqrt(sum(opt$par^2))
# hum.pmna.nm <- -EHUM(dat_3cat, beta.pmna.nm, d)
# hum.pmna.nm
# print(proc.time()-ptm)


boot.fn.pmna <- function(dat, d, b){
  set.seed(123+b)
  
  n <- nrow(dat)
  indx <- sample(1:n, n, replace=TRUE)
  boot.dat <- dat[indx,]
  
  beta0 <- rep(0.5, d)
  opt <- optim(par=beta0, PMNA, dat=boot.dat, d=d, method="L-BFGS-B", lower=-1, upper=1)
  # opt <- optim(par=beta0, PMNA, dat=boot.dat, d=d, method="Nelder-Mead")
  beta.pmna <- opt$par/sqrt(sum(opt$par^2))
  hum.pmna <- -EHUM(boot.dat, beta.pmna, d)
  
  filename <- paste0("pmna_boot_ericca_", b, ".RData")
  save(hum.pmna, beta.pmna, file=filename)
}


ptm <- proc.time()
B <- 100
mclapply(1:B, function(b) boot.fn.pmna(dat_3cat, d, b), mc.cores = detectCores()-1)
print(proc.time()-ptm)

hum.pmna.vec <- rep(NA, B)
beta.pmna.mat <- matrix(NA, d, B)
for(b in 1:B){
  filename <- paste0("pmna_boot_ericca_", b, ".RData")
  load(filename)
  
  hum.pmna.vec[b] <- hum.pmna
  beta.pmna.mat[,b] <- beta.pmna
}

se.hum.pmna <- sd(hum.pmna.vec)
se.beta.pmna <- apply(beta.pmna.mat, 1, sd)

se.hum.pmna
se.beta.pmna



##=========================
## Method 7: Min-Max
##=========================

MM <- function(dat, beta){
beta <- as.numeric(beta)

I <- function(x){
  x <- as.numeric(x)
  return(1*(x[3]>x[2])*(x[2]>x[1]))
}
  
  dat1 <- dat %>% filter(Y==0)
  dat2 <- dat %>% filter(Y==1)
  dat3 <- dat %>% filter(Y==2)
  
  n1 <- nrow(dat1)
  n2 <- nrow(dat2)  
  n3 <- nrow(dat3)
  
  var <- c("Xmax", "Xmin")
  
  l1 <- t(beta) %*% t(dat1[, var])
  l2 <- t(beta) %*% t(dat2[, var])  
  l3 <- t(beta) %*% t(dat3[, var])      
  
  L <- expand.grid(l1=l1, l2=l2, l3=l3)    
  S <- apply(L, 1, I)
  
  return(-sum(S)/(n1*n2*n3))
}


# 
# ptm <- proc.time()
# beta0 <- rep(0.5, 2)
# opt <- optim(par=beta0, MM, dat=dat_3cat, method="L-BFGS-B", lower=-1, upper=1)
# beta.mm.qn <- opt$par/sqrt(sum(opt$par^2))
# hum.mm.qn <- -MM(dat_3cat, beta.mm.qn)
# hum.mm.qn
# print(proc.time()-ptm)


ptm <- proc.time()
beta0 <- rep(0.5, 2)
opt <- optim(par=beta0, MM, dat=dat_3cat, method="Nelder-Mead")
beta.mm.nm <- opt$par/sqrt(sum(opt$par^2))
hum.mm.nm <- -MM(dat_3cat, beta.mm.nm)
hum.mm.nm
print(proc.time()-ptm)




boot.fn.mm <- function(dat, b){
  set.seed(123+b)
  
  n <- nrow(dat)
  indx <- sample(1:n, n, replace=TRUE)
  boot.dat <- dat[indx,]
  
  beta0 <- rep(0.5, 2)
  # opt <- optim(par=beta0, MM, dat=boot.dat, method="L-BFGS-B", lower=-1, upper=1)
  opt <- optim(par=beta0, MM, dat=boot.dat, method="Nelder-Mead")
  beta.mm <- opt$par/sqrt(sum(opt$par^2))
  hum.mm <- -MM(boot.dat, beta.mm)
  
  filename <- paste0("mm_boot_ericca_", b, ".RData")
  save(hum.mm, beta.mm, file=filename)
}


ptm <- proc.time()
B <- 100
mclapply(1:B, function(b) boot.fn.mm(dat_3cat, b), mc.cores = detectCores()-1)
print(proc.time()-ptm)

hum.mm.vec <- rep(NA, B)
beta.mm.mat <- matrix(NA, 2, B)
for(b in 1:B){
  filename <- paste0("mm_boot_ericca_", b, ".RData")
  load(filename)
  
  hum.mm.vec[b] <- hum.mm
  beta.mm.mat[,b] <- beta.mm
}

se.hum.mm <- sd(hum.mm.vec)
se.beta.mm <- apply(beta.mm.mat, 1, sd)

se.hum.mm
se.beta.mm


##=========================
## Method 1: Simple average 
##=========================
beta.avg <- c(1/d, d)
hum.avg <- -EHUM(dat_3cat, beta.avg, d)


data.frame(beta.avg=beta.avg, beta.ehum=beta.ehum.nm, se=se.beta.ehum, 
           beta.shum=beta.shum.qn, se=se.beta.shum, 
           beta.nhum=beta.nhum.qn, se=se.beta.nhum, 
           beta.pm=beta.pm.nm, se=se.beta.pm, 
           beta.pmna=beta.pmna.qn, se=se.beta.pmna)


data.frame(hum.avg=hum.avg, hum.ehum=hum.ehum.nm, se=se.hum.ehum,
           hum.shum=hum.shum.qn, se=se.hum.shum,
           hum.nhum=hum.nhum.qn, se=se.hum.nhum,
           hum.pm=hum.pm.nm, se=se.hum.pm,
           hum.pmna=hum.pmna.qn, se=se.hum.pmna,
           hum.mm=hum.mm.nm, se=se.hum.mm)



##################################################################
##================================================================
## Section 2: HUM value for individual biomarker (se)
##================================================================
##################################################################
indv.ehum.vec <- rep(0, d)
for(i in 1:d){
  beta <- rep(0, d)
  beta[i] <- 1
  
  indv.ehum.vec[i] <- -EHUM(dat_3cat, beta, d)
}

beta <- c(1,0,0,0)
hum_X1 <- -EHUM(dat_3cat, beta, d)

beta <- c(0,1,0,0)
hum_X2 <- -EHUM(dat_3cat, beta, d)

beta <- c(0,0,1,0)
hum_X3 <- -EHUM(dat_3cat, beta, d)

beta <- c(0,0,0,1)
hum_X4 <- -EHUM(dat_3cat, beta, d)

data.frame(hum_X1=hum_X1, hum_X2=hum_X2, hum_X3=hum_X3, hum_X4=hum_X4)


B <- 100
hum.X1.vec <- rep(0, B)
hum.X2.vec <- rep(0, B)
hum.X3.vec <- rep(0, B)
hum.X4.vec <- rep(0, B)

for(b in 1:B){
  set.seed(123+b)
  
  n <- nrow(dat_3cat)
  indx <- sample(1:n, n, replace=TRUE)
  boot.dat <- dat_3cat[indx,]
  
  beta <- c(1,0,0,0)
  hum.X1.vec[b] <- -EHUM(boot.dat, beta, d)
  
  beta <- c(0,1,0,0)
  hum.X2.vec[b] <- -EHUM(boot.dat, beta, d)
  
  beta <- c(0,0,1,0)
  hum.X3.vec[b] <- -EHUM(boot.dat, beta, d)
  
  beta <- c(0,0,0,1)
  hum.X4.vec[b] <- -EHUM(boot.dat, beta, d)
}

mean.hum <- c(mean(hum.X1.vec), mean(hum.X2.vec), mean(hum.X3.vec), mean(hum.X4.vec))
sd.hum <- c(sd(hum.X1.vec), sd(hum.X2.vec), sd(hum.X3.vec), sd(hum.X4.vec))
mean.hum
sd.hum


##################################################################
##================================================================
## Section 3: Box and density plots for individual and combine biomarkers
##================================================================
##################################################################
beta.shum.qn2 <- beta.shum.qn/sum(beta.shum.qn)
dat_3cat$Xshum <- as.numeric(as.matrix(dat_3cat[, c("X1", "X2", "X3", "X4")]) %*% beta.shum.qn2)
dat_3cat$Y <-as.factor(dat_3cat$Y)

## ggplot
library(ggplot2)
library(Rmisc)  # for multiplot function

g1 <- ggplot(dat_3cat, aes(x = AKI, y = X1, fill=AKI)) + geom_boxplot(show.legend=F) + 
  labs(y=expression(log(NGAL))) + ggtitle("NGAL at 0 hours")  

g2 <- ggplot(dat_3cat, aes(x = AKI, y = X2, fill=AKI)) + geom_boxplot(show.legend=F) + 
  labs(y=expression(log(NGAL))) + ggtitle("NGAL at 6 hours")  

g3 <- ggplot(dat_3cat, aes(x = AKI, y = X3, fill=AKI)) + geom_boxplot(show.legend=F) + 
  labs(y=expression(log(NGAL))) + ggtitle("NGAL at 12 hours")   

g4 <- ggplot(dat_3cat, aes(x = Y, y = X4, fill=Y)) + geom_boxplot(show.legend=F) + 
  labs(y=expression(log(NGAL)), x="AKI") + ggtitle("NGAL at 24 hours")  

gE <- ggplot(dat_3cat, aes(x = Y, y = Xshum, fill=Y)) + geom_boxplot(show.legend=F) + 
  labs(y=expression(log(NGAL)), x="AKI") + ggtitle("SSHUM")  

gS <- ggplot(dat_3cat, aes(x = Y, y = Xmean, fill=Y)) + geom_boxplot(show.legend=F) + 
  labs(y=expression(log(NGAL)), x="AKI") + ggtitle("Simple average")  


layout1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow=3, byrow=TRUE)
list1 <- list(g1, g2, g3, g4, gS, gE)


pdf("boxplot_ericca_v2.pdf", width = 10, height = 10)
multiplot(plotlist = list1, layout = layout1)
dev.off()


