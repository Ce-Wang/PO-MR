source("function.R")
boot <- function(n = n, dat = dat, est_time = est_time){
  r.all.all <- c()
  for(i in 1:200){
    cat(paste("####",i,"####\n"))
    set.seed(i)
    # i <- 136
    line <- sample(1:n,size=n,replace=T)
    data <- dat[line,] 
    pseudo=pseudosurv(time=data$ttime, event=data$event,tmax=est_time)$pseudo ##PO
    data$pseudo <- pseudo
    
    id <- c(1:n)
    data <- cbind(data,id)
    
    
    n <- NROW(data) 
    m <- sum(data$A) # number of observed subjects treatment
    m.c <- sum(data$A == 0) 
    
    ##################
    ####PS model######
    ##################
    p1 <- predict(glm(A~Z1+Z2+Z3,family="binomial", data = data),type="response") 
    p2 <- predict(glm(A~ X1 + X2 + X3,family="binomial", data = data),type="response") 
    
    #################
    #####回归模型####
    #################
    beta1=geese(pseudo~ Z1 + Z2 + A,scale.fix=T,data=data,family=gaussian, id=id, jack=F, mean.link="cloglog",corstr="independence")$beta
    m11 =  1- exp(-exp(beta1[1] + beta1[2]* data$Z1+beta1[3]*data$Z2 + beta1[4]))
    m10 =  1- exp(-exp(beta1[1] + beta1[2]* data$Z1+beta1[3]*data$Z2))
    
    
    beta2=geese(pseudo~X1+X2+A,scale.fix=T,data=data,family=gaussian, id=id, jack=F, mean.link="cloglog",corstr="independence")$beta
    m21 =  1- exp(-exp(beta2[1] + beta2[2]* data$X1+beta2[3]*data$X2 + beta2[4]))
    m20 =  1- exp(-exp(beta2[1] + beta2[2]* data$X1+beta2[3]*data$X2))
    
    ##########################
    #####PS###################
    ##########################
    ##PS1
    weight_p1 <- data$A/p1 + (1-data$A)/(1-p1)
    ewei_11 <- sum(weight_p1[data$A == 1] * pseudo[data$A == 1])/n
    ewei_10 <- sum(weight_p1[data$A == 0] * pseudo[data$A == 0])/n
    r1 <- ewei_11 - ewei_10
    
    ##PS2
    weight_p2 <- data$A/p2 + (1-data$A)/(1-p2)
    ewei_21 <- sum(weight_p2[data$A == 1] * pseudo[data$A == 1])/n
    ewei_20 <- sum(weight_p2[data$A == 0] * pseudo[data$A == 0])/n
    r2 <- ewei_21 - ewei_20
    
    ##########################
    #####OR###################
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
    MR1000 <- estimate.t - estimate.c #倾向模型指定正确 回归模型指定错
    
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
    
    r.all <- c(r1,r2,r3,r4,MR1000,MR0100,MR0010,MR0001,MR1100,MR1010,MR1001,MR0110,MR0101,MR0011,MR1110,MR1101,MR1011,MR0111,MR1111)
    r.all.all <- rbind(r.all.all, r.all)
  }
  sd.all <- apply(r.all.all, 2, sd)
  return(sd.all)
}

################################################################################
########################time = 10, lamda = 3, censor_rate=-4, n = 200###########
################################################################################
#true value
A_all <- c()
for (i in 1:1000){
  A <- true_race(n.sample = 1000000, time = 10,seed = i,lamda = 3) #PH
  A_all <- c(A_all, A)
}
mean(A_all)

cr.all <- c()
cr.all_sd <- c()
estimate_all <- c()
for(a in 1:1000){
  cat(paste("####",a,"####\n"))
  dat <- data_generate(n.sample = 200, seed = a,lamda = 3,censor_rate = -4) ##lamda determine PH ;censor rate
  sd.r<-boot(n = 200, dat = dat, est_time = 10) #n sample size
  B <- double_estimator(n.sample = 200,seed = a,est_time = 10,lamda = 3, censor_rate = -4) 
  cr.l <- B - 1.96*sd.r
  cr.u <- B + 1.96*sd.r
  cr.r <- ifelse(cr.l <= A & cr.u >= A,1,0)
  cr.all <- rbind(cr.all, cr.r)
  cr.all_sd <- rbind(cr.all_sd,sd.r)
  estimate_all <-  rbind(estimate_all, B)
  # sd.all<- rbind(sd.all,sd.r)
}
all_results <- list(cr.all, cr.all_sd, estimate_all)
