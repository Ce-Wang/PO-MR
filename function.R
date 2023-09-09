library(pseudo)
library(survRM2)
library(geepack)
library(locfit)
############################
######true value############
############################
true_race <- function(
  n.sample = n.sample,
  time = time,
  seed = seed,
  lamda = 3){
  
  set.seed(seed)
  n.sample <- n.sample 
  
  rho=sqrt(0.2)
  nrho=sqrt(1-rho^2)
  
  Z0=rnorm(n.sample)
  Z1=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
  Z2=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
  Z3=pmin(4,pmax(-4,rnorm(n.sample)))
  u <- runif(n.sample) 
  
  L1 <- -3-1+0.5*(Z1+Z2)
  L0 <- -3+0+0.5*(Z1+Z2)
  
  T1 <- (-log(u)/(0.005*exp(L1)))^(1/3)  
  T0 <- (-log(u)/(0.005*exp(L0)))^(1/lamda)

  time <- time 
  
  Y1_r <- sum(T1 >= time)/n.sample
  Y0_r <- sum(T0 >= time)/n.sample
  A <- Y1_r - Y0_r #
  return(A)
}


############################################
data_generate <- function(
  n.sample = n.sample, 
  seed = seed,  
  lamda = 3,
  censor_rate = -3.5
){
  
  set.seed(seed)
  n.sample <- n.sample 
  
  ###generate covariate####
  rho=sqrt(0.2)
  nrho=sqrt(1-rho^2)
  
  Z0=rnorm(n.sample)
  Z1=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
  Z2=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
  Z3=pmin(4,pmax(-4,rnorm(n.sample)))
  
  X1 <- Z1^2
  X2 <- Z2^2
  X3 <- Z3^2
  
  ###generate treatment###
  l <- expit((Z1+Z2+Z3)/3)
  A <- rbinom( n.sample, 1,l)
  L <- rep(NA, n.sample)
  
  u <- runif(n.sample)
  
  L[A==1] <- -3-1+0.5*(Z1[A==1]+Z2[A==1])
  L[A==0] <- -3+0+0.5*(Z1[A==0]+Z2[A==0])
  
  TT <- rep(NA, n.sample)
  
  TT[A==1] <- ((-log(u[A==1]))/(0.005*exp(L[A==1])))^(1/3)
  TT[A==0] <- ((-log(u[A==0]))/(0.005*exp(L[A==0])))^(1/lamda)

  C=rexp(n.sample,exp(censor_rate))
  ttime=pmin(TT,C)
  event=1*(TT<C)
  data.tt <- as.data.frame(cbind(Z1, Z2, Z3, X1, X2, X3, A, ttime, event))
  return(data.tt)
}



double_estimator <- function(
  n.sample = n.sample,
  seed = seed,
  pmodel="C",tmodel="C",
  est_time = time, 
  lamda = 3,
  censor_rate = -3.5
){
  dat <- data_generate(n.sample = n.sample, seed = seed,lamda = lamda, censor_rate = censor_rate) 
  pseudo=pseudosurv(time=dat$ttime, event=dat$event,tmax=est_time)$pseudo ##伪观察计算时间
  dat$pseudo <- pseudo
  
  id <- c(1:n.sample)
  dat <- cbind(dat,id)
  data <- dat
  
  n <- NROW(data) 
  m <- sum(data$A) 
  m.c <- sum(data$A == 0) 
  
  ##################
  ####ps model######
  ##################
  p1 <- predict(glm(A~Z1+Z2+Z3,family="binomial", data = data),type="response") 
  p2 <- predict(glm(A~ X1 + X2 + X3,family="binomial", data = data),type="response")
  
  #################
  #####OR model####
  #################
  beta1=geese(pseudo~ Z1 + Z2 + A,scale.fix=T,data=data,family=gaussian, id=id, jack=F, mean.link="cloglog",corstr="independence")$beta
  m11 =  1- exp(-exp(beta1[1] + beta1[2]* data$Z1+beta1[3]*data$Z2 + beta1[4]))
  m10 =  1- exp(-exp(beta1[1] + beta1[2]* data$Z1+beta1[3]*data$Z2))
  
  
  beta2=geese(pseudo~X1+X2+A,scale.fix=T,data=data,family=gaussian, id=id, jack=F, mean.link="cloglog",corstr="independence")$beta
  m21 =  1- exp(-exp(beta2[1] + beta2[2]* data$X1+beta2[3]*data$X2 + beta2[4]))
  m20 =  1- exp(-exp(beta2[1] + beta2[2]* data$X1+beta2[3]*data$X2))
  
  ##########################
  #####ps###################
  ##########################
  ##ps1
  weight_p1 <- data$A/p1 + (1-data$A)/(1-p1)
  ewei_11 <- sum(weight_p1[data$A == 1] * pseudo[data$A == 1])/n
  ewei_10 <- sum(weight_p1[data$A == 0] * pseudo[data$A == 0])/n
  r1 <- ewei_11 - ewei_10
  
  ##ps2
  weight_p2 <- data$A/p2 + (1-data$A)/(1-p2)
  ewei_21 <- sum(weight_p2[data$A == 1] * pseudo[data$A == 1])/n
  ewei_20 <- sum(weight_p2[data$A == 0] * pseudo[data$A == 0])/n
  r2 <- ewei_21 - ewei_20
  
  ##########################
  #####or#################
  ##########################
  ##OR1
  r3 <- mean(m11) - mean(m10)
  ##OR2
  r4 <- mean(m21) - mean(m20)
  
  ##########################
  ######MR##################
  ##########################
  Fn1 <- function(rho, ghat){ -sum(log(1 + ghat %*% rho)) }
  Grd1 <- function(rho, ghat){ -colSums(ghat / c(1 + ghat %*% rho)) }
  
  Fn0 <- function(rho, chat){ -sum(log(1 + chat %*% rho)) }
  Grd0 <- function(rho, chat){ -colSums(chat / c(1 + chat %*% rho)) }
  ##########
  ###1010###
  ##########
  J <- 1 
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1 
  
  c.hat[, 1] <- 1-p1 
  
  g.hat[, 2] <- m11 
  
  c.hat[, 2] <- m10 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1010 <- estimate.t - estimate.c 
  
  ##########
  ###1001###
  ##########
  J <- 1 
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1 
  
  c.hat[, 1] <- 1-p1
  
  g.hat[, 2] <- m21 
  
  c.hat[, 2] <- m20 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1001 <- estimate.t - estimate.c
  
  ##########
  ###0110###
  ##########
  J <- 1 
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p2 
  
  c.hat[, 1] <- 1-p2 
  
  g.hat[, 2] <- m11 
  
  c.hat[, 2] <- m10
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR0110 <- estimate.t - estimate.c 
  
  ##########
  ###0101###
  ##########
  J <- 1 
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p2 
  
  c.hat[, 1] <- 1-p2 
  
  g.hat[, 2] <- m21 
  
  c.hat[, 2] <- m20 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR0101 <- estimate.t - estimate.c
  
  ##########
  ###1110###
  ##########
  J <- 2 
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1
  g.hat[, 2] <- p2 
  
  c.hat[, 1] <- 1-p1
  c.hat[, 2] <- 1-p2 
  
  g.hat[, 3] <- m11 
  
  c.hat[, 3] <- m10
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1110 <- estimate.t - estimate.c 
  
  ##########
  ###1110###
  ##########
  J <- 2
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1 
  g.hat[, 2] <- p2
  
  c.hat[, 1] <- 1-p1
  c.hat[, 2] <- 1-p2 
  
  g.hat[, 3] <- m21
  
  c.hat[, 3] <- m20 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1101 <- estimate.t - estimate.c 
  
  ##########
  ###1011###
  ##########
  J <- 1 
  K <- 2 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1 
  c.hat[, 1] <- 1-p1
  
  g.hat[, 2] <- m11
  c.hat[, 2] <- m10 
  
  g.hat[, 3] <- m21
  c.hat[, 3] <- m20 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1011 <- estimate.t - estimate.c 
  
  
  ##########
  ###0111###
  ##########
  J <- 1 
  K <- 2 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p2 
  c.hat[, 1] <- 1-p2
  
  g.hat[, 2] <- m11
  c.hat[, 2] <- m10 
  
  g.hat[, 3] <- m21
  c.hat[, 3] <- m20 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR0111 <- estimate.t - estimate.c 
  
  ##########
  ###1111###
  ##########
  J <- 2
  K <- 2
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1 
  c.hat[, 1] <- 1-p1
  
  g.hat[, 2] <- p2
  c.hat[, 2] <- 1-p2
  
  g.hat[, 3] <- m11
  c.hat[, 3] <- m10 
  
  g.hat[, 4] <- m21
  c.hat[, 4] <- m20 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1111 <- estimate.t - estimate.c 
  
  
  ##########
  ###1000###
  ##########
  J <- 1 
  K <- 0 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1 
  
  c.hat[, 1] <- 1-p1
 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1000 <- estimate.t - estimate.c 
  
  ##########
  ###0100###
  ##########
  J <- 1 
  K <- 0
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p2 
  
  c.hat[, 1] <- 1-p2
  

  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR0100 <- estimate.t - estimate.c 
  
  
  ##########
  ###0010###
  ##########
  J <- 0 
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  
  g.hat[, 1] <- m11 
  
  c.hat[, 1] <- m10 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR0010 <- estimate.t - estimate.c 
  
  ##########
  ###0001###
  ##########
  J <- 0 
  K <- 1 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  
  g.hat[, 1] <- m21 
  
  c.hat[, 1] <- m20 
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR0001 <- estimate.t - estimate.c
  
  
  
  ##########
  ###1100###
  ##########
  J <- 2 
  K <- 0 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
  g.hat[, 1] <- p1 
  g.hat[, 2] <- p2
  
  c.hat[, 1] <- 1- p1 
  c.hat[, 2] <- 1- p2

  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR1100 <- estimate.t - estimate.c
  
  
  ##########
  ###0011###
  ##########
  J <- 0 
  K <- 2 
  
  g.col <- J + K
  c.col <- J + K
  
  g.hat <- matrix(0, n, g.col)
  c.hat <- matrix(0, n, c.col)
  
 
  g.hat[, 1] <- m11
  c.hat[, 1] <- m10 
  
  g.hat[, 2] <- m21 
  c.hat[, 2] <- m20 
  
  
  g.hat <- scale(g.hat, center = TRUE, scale = F)[data$A == 1, ]
  g.hat <- matrix(data = g.hat, nrow = m)
  
  c.hat <- scale(c.hat, center = TRUE, scale = F)[data$A == 0, ]
  c.hat <- matrix(data = c.hat, nrow = m.c, )
  
  # calculate the weights
  rho.hat1 <- constrOptim(theta = rep(0, g.col), f = Fn1, grad = Grd1, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
  wts1 <- c(1 / m / (1 + g.hat %*% rho.hat1))
  wts1 <- wts1 / sum(wts1)
  estimate.t <- sum(data$pseudo[data$A == 1] * wts1)
  
  rho.hat0 <- constrOptim(theta = rep(0, c.col), f = Fn0, grad = Grd0, ui = c.hat, ci = rep(1 / m.c - 1, m.c), chat = c.hat)$par
  wts0 <- c(1 / m.c / (1 + c.hat %*% rho.hat0))
  wts0 <- wts0 / sum(wts0)
  estimate.c <- sum(data$pseudo[data$A == 0] * wts0)
  MR0011 <- estimate.t - estimate.c 
  
  
  r.all <- cbind(r1,r2,r3,r4,MR1000,MR0100,MR0010,MR0001,MR1100,MR1010,MR1001,MR0110,MR0101,MR0011,MR1110,MR1101,MR1011,MR0111,MR1111)
  return(r.all)
}