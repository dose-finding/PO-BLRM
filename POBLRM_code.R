#######################################################################
#######################################################################
##                                                                   ##
##  Partial Ordering Bayesian Logistic Regression Model for Phase I  ##
##  Combination trials and Computationally Efficient Approach to     ##
##  Operational Prior Specification                                  ##
##                                                                   ##
##  R code                                                           ##
##                                                                   ##
##  Weishi Chen and Pavel Mozgunov                                   ##
##                                                                   ##
##  Last edit: 6th May 2024                                          ##
##                                                                   ##
#######################################################################
#######################################################################

#
#  R packeges used
#
library(nnet)
library(dplyr); library(magrittr); library(tidyr); library(knitr)
library(rjags); library(coda)

#
#  Define standardised dose levels
#
StandardDose <- function(p1 = NULL, nu = NULL, k, mu1, mu2, sigma2){
  # Inputs:
  # -------------------------------------------------
  # p1    : prior toxicity at the lost dose
  # nu    : spacing between prior toxicity probs
  # k     : number of dose levels
  # mu1   : prior mean of theta1
  # mu2   : prior mean of log(theta2)
  # sigma2: prior sd of log(theta2)
  #
  # Outputs:
  # -------------------------------------------------
  # d: vector of length k, standardised dose levels 
  p_0 <- seq(from = p1, by = nu, length.out = k)
  theta1 <- mu1
  theta2 <- exp(mu2 + sigma2^2/2)
  d <- (log(p_0/(1-p_0)) - theta1) / theta2
  return(d)
}

#
# KL divergence to match priors
#
KL <- function(pseudo, mu1, mu2, sigma1, sigma2, d, n.sim) {
  # Inputs:
  # -------------------------------------------------------
  # pseudo        : starting point of the optimisation algorithm
  #                 pseudo prior in the order (y_{-1}, n_{-1}, y_0, n_0)
  # mu1, mu2      : prior means of theta1 and log(theta2)
  # sigma1, sigma2: prior var of theta1 and log(theta2)
  # d             : k-vector, standardised dose levels
  # n.sim         : number of samples simulated from the prior
  #
  # Outputs:
  # ------------------------------------------------------
  # KL : matched pseudo prior from the given normal prior
  #      4-vector of the form (y_{-1}, n_{-1}, y_0, n_0)
  KL_int <- function(theta1, theta2, y_minus1, y_0, n_minus1, n_0, mu1, mu2, sigma1, sigma2, d) {
    k <- length(d)
    KL_int <- -log(2*pi) - log(sigma1) - log(sigma2) - log(theta2) - (theta1-mu1)^2/(2*sigma1^2) -
      (log(theta2)-mu2)^2/(2*sigma2^2)
    KL_int <- KL_int - log(d[k]-d[1]) + log(beta(y_minus1,n_minus1-y_minus1)) + log(beta(y_0, n_0-y_0))
    KL_int <- KL_int + y_minus1*log(1+exp(-theta1-theta2*d[1])) + 
      (n_minus1-y_minus1)*log(1+exp(theta1+theta2*d[1]))
    KL_int <- KL_int + y_0*log(1+exp(-theta1-theta2*d[k])) + (n_0-y_0)*log(1+exp(theta1+theta2*d[k]))
    KL_int
  }
  if(pseudo[1]>pseudo[2] | pseudo[3]>pseudo[4] | any(pseudo<0)){
    KL <- 1e+8
  } else{
    y_minus1 <- pseudo[1]; n_minus1 <- pseudo[2]
    y_0 <- pseudo[3]; n_0 <- pseudo[4]
    theta.data <- data.frame(theta1 = rnorm(n.sim, mu1, sigma1),
                             theta2 = rlnorm(n.sim, mu2, sigma2))
    KL <- theta.data %>% as_tibble() %>%
      transmute(KL = KL_int(theta1, theta2, y_minus1, y_0, n_minus1, n_0,
                            mu1, mu2, sigma1, sigma2, d)) %>% 
      filter(!is.infinite(KL)) %>%
      colMeans() %>% as.numeric()
  }
  return(KL)
}
#
#  posterior distribution for BLRM selection
#
posterior <- function(theta1, theta2, mu1, mu2, sigma1, sigma2, d, y, n){
  prior <- exp(dnorm(theta1, mu1, sigma1, log=T) + dlnorm(theta2, mu2, sigma2, log=T))
  p <- exp(theta1+theta2*d)/(1+exp(theta1+theta2*d))
  lik <- exp(sum(dbinom(y, n, p, log=TRUE)))
  prior*lik
}

##########################################################################################
#
#  Fit the Bayesian POCRM (Wages et al., 2011)
#
bpocrm<-function(p0, p.skel, sigma, ttr, cohortsize, ncohort, start.comb){
  if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel))
  nord.tox = nrow(p.skel)                 # M
  mprior.tox = rep(1/nord.tox, nord.tox)  # p(m)
  bcrmh <- function(a,p,y,n){
    s2 <- sigma^2                               # the variance of a
    lik <- exp(-0.5/s2*a^2)                  # the prior g(a)
    for(j in 1:length(p)){
      pj <- p[j]^exp(a)                      # psi_m(d_j, a)
      lik <- lik*pj^y[j]*(1-pj)^(n[j]-y[j])  # binomial mass
    }
    return(lik)
  }
  bcrmht <- function(a, p, y, n){
    lik <- a * bcrmh(a, p, y, n)
    return(lik)
  }
  ncomb <- ncol(p.skel)  # k
  y <- rep(0, ncomb) 
  n <- rep(0, ncomb)
  #  start the trail at a dose level pre-specified by start.comb
  comb.curr <- start.comb
  ptox.hat <- numeric(ncomb)  # \hat{R}(d_j)
  comb.select <- rep(0, ncomb)
  i <- 1      # index for number of cohorts
  #  Now, we start running the trail
  while (i <= ncohort) {
    y[comb.curr] <- y[comb.curr] + rbinom(1, cohortsize, p0[comb.curr])
    n[comb.curr] <- n[comb.curr] + cohortsize
    marginal.tox <- rep(0, nord.tox)
    for (k in 1:nord.tox) {
      marginal.tox[k] <- integrate(bcrmh, lower=-Inf, upper=Inf, p=p.skel[k,], y=y, n=n, abs.tol=0)$value
    }
    if(nord.tox>1){
      mtox.sel <- which.is.max(marginal.tox)   # h
    } else{
      mtox.sel <- 1
    }
    # estimate a under selected model
    est.tox <- integrate(bcrmht, lower=-5, upper=5, p.skel[mtox.sel,], y, n, abs.tol=0)$value / marginal.tox[mtox.sel]
    #  Under the selected model *mtox.sel*, update toxicity prob.s
    ptox.hat <- p.skel[mtox.sel,]^(exp(est.tox)) 
    distance <- abs(ptox.hat-ttr)
    comb.best <- which.is.max(-distance)  # since we want to minimise
    comb.curr <- comb.best
    i <- i+1 # recruit the next cohort
  }
  comb.select[comb.curr] <- comb.select[comb.curr] + 1
  return(list(comb.select=comb.select, tox.data=y, pt.allocation=n, ptox.hat=ptox.hat))
}

###########################################################################################
#
#  Fit the proposed POBLRM model
#
POBLRM <- function(d, orders, pm, mu1, mu2, sigma1, sigma2, pseudo, p_true, delta,
                   p0=NULL, c2=1, S=1e+4, burnin=2500, chains=1,
                   start.d=1, ncohort=12, cohortsize=3, n.stop=30, n.sim=1e+4, 
                   print=FALSE, randomised=FALSE){
  # Simulate from posterior
  jags.script <- "
  model {
    # Likelihood
    for(j in 1:length(y)) {
      p[j] <- 1/(1+exp(-theta1-theta2*x[j]))
      y[j] ~ dbin(p[j], n[j])
    }
    # Prior
    theta[1:2] ~ dmnorm(priorMean, priorPrec)
    theta1 <- theta[1]
    theta2 <- exp(theta[2])
  }
  "
  if (randomised) {
    # set up
    M <- length(pm)
    k <- length(p_true)
    y <- y0 <- n <- decisions <- numeric()
    tox.data <- pt.allocation <- d.select <- rep(0, k)
    tox.data[2] <- tox.data[2] + pseudo[1]
    pt.allocation[2] <- pt.allocation[2] + pseudo[2]
    tox.data[k] <- tox.data[k] + pseudo[3]
    pt.allocation[k] <- pt.allocation[k] + pseudo[4]
    d0 <- (log(p0/(1-p0)) - mu1) / exp(mu2 + sigma2^2/2)
    # decisions: n-vector, decisions made at each step
    d.curr <- start.d
    stop <- 0   # indicate if trial stops early
    j <- 1      # index for number of cohorts
    #  Now, we start running the trail
    while (j <= ncohort) {
      #  generate data for the new cohort Y ~ Bin(n, \psi_m(d_j,a))
      y[j] <- rbinom(1, cohortsize, p_true[d.curr])
      y0[j] <- rbinom(1, c2, p_true[1])
      tox.data[d.curr] <- tox.data[d.curr] + y[j]
      tox.data[1] <- tox.data[1] + y0[j]
      n[j] <- cohortsize
      pt.allocation[d.curr] <- pt.allocation[d.curr] + cohortsize
      pt.allocation[1] <- pt.allocation[1] + c2
      decisions[j] <- d.curr
      #  stop if the number of patients treated at the current level exceed *n.stop*.
      if(any(pt.allocation>n.stop)){
        stop <- 0
        break
      }
      # AIC model selection
      mod.select <- function(m, y, n){
        d_m <- c(d0, d)[order(orders[m,])]
        logistic.fit <- glm(cbind(y, n-y) ~ d_m, family = 'binomial')
        logistic.fit$aic
      }
      AIC <- sapply(1:M, mod.select, y=tox.data, n=pt.allocation)
      resample <- function(x, ...) x[sample.int(length(x), ...)]
      if(M>1){
        mtox.sel <- resample(which(AIC==min(AIC)), 1)
      } else{
        mtox.sel <- 1
      }
      # update theta under selected order
      order.sel <- orders[mtox.sel,]
      d.sel <- c(d0,d)[order(orders[mtox.sel,])]
      model.fit <- jags.model(textConnection(jags.script), quiet = TRUE,
                              data=list(x=c(rep(d0,j), d.sel[decisions]), 
                                        y=c(y0, y), priorMean=c(mu1, mu2), 
                                        priorPrec=diag(c(sigma1^(-2), sigma2^(-2))), 
                                        n=c(rep(c2,j), n)),
                              n.chains=chains, n.adapt=S)
      update(model.fit, S, progress.bar="none")
      tt<-jags.samples(model.fit, c('theta1','theta2'), S, progress.bar="none")
      theta1 <- mean(tt$theta1[1,,])
      theta2 <- mean(tt$theta2[1,,])
      ptox.hat <- 1/(1+exp(-theta1-theta2*d.sel))
      d.curr <- which.is.max(-abs(ptox.hat-delta))
      if(print){
        cat(paste('cohort', j, '\n'))
        cat(paste('Selected order: ', mtox.sel, '\n'))
        print(theta.hat)
        print(ptox.hat)
        cat(paste('Selected dose: ', d.curr, '\n'))
        cat('-------------------------------------------\n')
      }
      j <- j+1 # recruit the next cohort
    }
  } else {
    # set up
    M <- length(pm)
    k <- length(p_true)
    y <- n <- decisions <- numeric()
    tox.data <- pt.allocation <- d.select <- rep(0, k)
    tox.data[1] <- tox.data[1] + pseudo[1]
    pt.allocation[1] <- pt.allocation[1] + pseudo[2]
    tox.data[k] <- tox.data[k] + pseudo[3]
    pt.allocation[k] <- pt.allocation[k] + pseudo[4]
    # decisions: n-vector, decisions made at each step
    d.curr <- start.d
    stop <- 0   # indicate if trial stops early
    j <- 1      # index for number of cohorts
    #  Now, we start running the trail
    while (j <= ncohort) {
      #  generate data for the new cohort Y ~ Bin(n, \psi_m(d_j,a))
      y[j] <- rbinom(1, cohortsize, p_true[d.curr])
      tox.data[d.curr] <- tox.data[d.curr] + y[j]
      n[j] <- cohortsize
      pt.allocation[d.curr] <- pt.allocation[d.curr] + cohortsize
      decisions[j] <- d.curr
      #  stop if the number of patients treated at the current level exceed *n.stop*.
      if(any(pt.allocation>n.stop)){
        stop <- 0
        break
      }
      # AIC model selection
      mod.select <- function(m, y, n){
        d_m <- d[order(orders[m,])]
        logistic.fit <- glm(cbind(y, n-y) ~ d_m, family = 'binomial')
        logistic.fit$aic
      }
      AIC <- sapply(1:M, mod.select, y=tox.data, n=pt.allocation)
      resample <- function(x, ...) x[sample.int(length(x), ...)]
      if(M>1){
        mtox.sel <- resample(which(AIC==min(AIC)), 1)
      } else{
        mtox.sel <- 1
      }
      # update theta under selected order
      order.sel <- orders[mtox.sel,]
      d.sel <- d[order(orders[mtox.sel,])]
      model.fit <- jags.model(textConnection(jags.script), quiet = TRUE,
                              data=list(k=k, x=d.sel[decisions], y=y, n=n,
                                        priorMean=c(mu1, mu2), 
                                        priorPrec=diag(c(sigma1^(-2), sigma2^(-2)))),
                              n.chains=chains, n.adapt=S)
      update(model.fit, S, progress.bar="none")
      tt<-jags.samples(model.fit, c('theta1','theta2'), S, progress.bar="none")
      theta1 <- mean(tt$theta1[1,,])
      theta2 <- mean(tt$theta2[1,,])
      ptox.hat <- 1/(1+exp(-theta1-theta2*d.sel))
      d.curr <- which.is.max(-abs(ptox.hat-delta))
      if(print){
        cat(paste('cohort', j, '\n'))
        cat(paste('Selected order: ', mtox.sel, '\n'))
        print(theta.hat)
        print(ptox.hat)
        cat(paste('Selected dose: ', d.curr, '\n'))
        cat('-------------------------------------------\n')
      }
      j <- j+1 # recruit the next cohort
    }
  }
  
  # final selected dose
  d.select[d.curr] <- d.select[d.curr] + 1
  return(list(d.select=d.select, tox.data=tox.data-c(0, pseudo[1], rep(0, k-3), pseudo[3]),
              pt.allocation=pt.allocation-c(0, pseudo[2], rep(0, k-3), pseudo[4]),
              ptox.hat=ptox.hat))
}

############################################################################################################
#
#  Fit the two-dimensional BLRM model (Neuenschwander et al.,2015)
#
# rjags code for MCMC
model1.string <- "
model {
  for (i in 1:N1){
    lin.pred.p.1[i] <- alpha0[1] + alpha1[1] * d[i]
    lin.pred.p.2[i] <- ifelse(lin.pred.p.1[i] < -10, -10, ifelse(lin.pred.p.1[i]> 10, 10, lin.pred.p.1[i]))
    odds1[i]<- exp(lin.pred.p.2[i])
  }

  for (j in 1:N2){
    lin.pred.q.1[j] <- alpha0[2] + alpha1[2] * f[j]
    lin.pred.q.2[j] <- ifelse(lin.pred.q.1[j] < -10, -10, ifelse(lin.pred.q.1[j]> 10, 10, lin.pred.q.1[j]))
    odds2[j]<- exp(lin.pred.q.2[j])
  }

  for (i in 1:N1){
    for (j in 1:N2){
      odds0[i,j]<-odds1[i] + odds2[j] + odds1[i] * odds2[j]
      odds[i,j]<-odds0[i,j]*exp(eta * pow(d[i],power) * pow(f[j],power) )
      tox[i,j]<-odds[i,j]/(1+odds[i,j])
      s[i,j] ~ dbin(tox[i,j], n[i,j])		
    }
  }
  eta~dnorm(meanEta,priorEta)
  for(t in 1:2){
    theta[1:2, t] ~ dmnorm(priorMean[1:2, t], priorPrec[1:2, 1:2, t])
    ## extract actual coefficients
    alpha0[t] <- theta[1, t]
    alpha1[t] <- exp(theta[2, t])
  }
}
"
model1.spec<-textConnection(model1.string)

BLRM2 <- function(ptrue, a, b, mu_a1, mu_a2, mu_b1, mu_b2, mu_zeta, 
                  sigma_a1, sigma_a2, sigma_b1, sigma_b2, sigma_zeta,
                  delta, start.d=c(1, 1), ncohort=12, cohortsize=3, n.stop=36, 
                  iter=1e+4, print=FALSE) {
  # Transforming the prior variance to prior precision
  priorPrec <- array(0, dim = c(2,2,2))
  priorPrec[,,1] <- matrix(c(sigma_a1^(-2), 0, 0, sigma_a2^(-2)), nrow=2)
  priorPrec[,,2] <- matrix(c(sigma_b1^(-2), 0, 0, sigma_b2^(-2)), nrow=2)
  # prior mean
  priorMean <- matrix(c(mu_a1, mu_a2, mu_b1, mu_b2), nrow=2)
  # Creating matrices to store the estimates
  # Rows are the dose levels of a, Columns are the dose level of b
  prob.mtd <- toxicity <- d.select <- pt.allocation <- mat.or.vec(length(a),length(b))
  n <- 1
  d.curr <- start.d
  while (n<=ncohort) {
    toxicity[d.curr[1], d.curr[2]] <- toxicity[d.curr[1], d.curr[2]] + 
      rbinom(1, cohortsize, ptrue[d.curr[1], d.curr[2]])
    pt.allocation[d.curr[1], d.curr[2]] <- pt.allocation[d.curr[1], d.curr[2]] + cohortsize
    # Uploading MCMC Model
    model1.spec<-textConnection(model1.string)
    # Running MCMC Code
    mydata <- list(n=pt.allocation, s=toxicity, N1=length(a), N2=length(b),
                   d=a, f=b, priorMean=priorMean, priorPrec=priorPrec, 
                   priorEta=sigma_zeta^(-2), meanEta=mu_zeta, power=1)
    jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iter,quiet=TRUE)
    update(jags, iter,progress.bar="none")
    # Storing Samples of the Posterior Distribution
    tt<-jags.samples(jags,c('alpha0','alpha1','eta'),iter,progress.bar="none")
    a01<-tt$alpha0[1,,]
    a02<-tt$alpha0[2,,]
    a11<-(tt$alpha1[1,,])
    a12<-(tt$alpha1[2,,])
    e<-(tt$eta[1,,])
    # Computing the estimates using the Posteriors of the parameters
    posterior.tox.samples<-array(9,dim=c(length(a),length(b),iter))
    for (j in 1:length(a)){
      for (k in 1:length(b)){
        help1 <- a01 + a11 * a[j]
        help1[which(help1<(-10))] <- (-10); help1[which(help1>( 10))] <- ( 10)
        help2 <- a02 + a12 * b[k]
        help2[which(help2<(-10))] <- (-10); help2[which(help2>( 10))] <- ( 10)
        odds1 <- exp(help1); odds2 <- exp(help2)
        odds0 <- odds1 + odds2 + odds1*odds2
        odds <- odds0 * exp(e * a[j] * b[k])
        distr <- odds/(1+odds)
        posterior.tox.samples[j,k,] <- distr
        prob.mtd[j,k]<-mean(distr)
      }
    }
    d.curr <- which.is.max(-abs(prob.mtd-delta))
    row <- ceiling(d.curr/length(a))
    column <- ifelse((d.curr%%length(b))!=0, d.curr%%length(a), length(b))
    d.curr <- c(row, column)
    if(print){
      cat(paste('cohort', n, '\n y: '))
      print(toxicity)
      print(prob.mtd)
      cat(paste('Selected dose: ', d.curr, '\n'))
      cat('-------------------------------------------\n')
    }
    n <- n + 1
  }
  # final selected dose
  d.select[d.curr[1],d.curr[2]] <- d.select[d.curr[1], d.curr[2]] + 1
  return(list(d.select=d.select, tox.data=toxicity, pt.allocation=pt.allocation, ptox.hat=prob.mtd))
}

##############################################################################################
#
# Define scenarios
#
k <- 9
C <- 20
p_true_mat <- matrix(ncol=k, nrow=C)
p_true_mat[1,] <- c(0.40, 0.45, 0.60, 0.50, 0.60, 0.70, 0.60, 0.70, 0.80)
p_true_mat[2,] <- c(0.30, 0.40, 0.50, 0.40, 0.50, 0.65, 0.50, 0.60, 0.70)
p_true_mat[3,] <- c(0.20, 0.30, 0.50, 0.40, 0.50, 0.60, 0.50, 0.60, 0.70)
p_true_mat[4,] <- c(0.10, 0.20, 0.30, 0.20, 0.40, 0.50, 0.40, 0.50, 0.60)
p_true_mat[5,] <- c(0.20, 0.40, 0.50, 0.30, 0.50, 0.60, 0.50, 0.60, 0.70)
p_true_mat[6,] <- c(0.10, 0.15, 0.20, 0.20, 0.30, 0.50, 0.40, 0.50, 0.60)
p_true_mat[7,] <- c(0.10, 0.20, 0.25, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60)
p_true_mat[8,] <- c(0.05, 0.10, 0.25, 0.10, 0.20, 0.50, 0.30, 0.50, 0.60)
p_true_mat[9,] <- c(0.05, 0.10, 0.40, 0.10, 0.20, 0.50, 0.20, 0.30, 0.60)
p_true_mat[10,] <- c(0.05, 0.15, 0.20, 0.15, 0.20, 0.25, 0.20, 0.25, 0.30)
p_true_mat[11,] <- c(0.20, 0.30, 0.50, 0.30, 0.40, 0.60, 0.40, 0.50, 0.70)
p_true_mat[12,] <- c(0.10, 0.30, 0.40, 0.20, 0.40, 0.50, 0.30, 0.50, 0.60)
p_true_mat[13,] <- c(0.10, 0.20, 0.30, 0.30, 0.40, 0.50, 0.40, 0.50, 0.60)
p_true_mat[14,] <- c(0.10, 0.20, 0.30, 0.20, 0.30, 0.40, 0.40, 0.50, 0.60)
p_true_mat[15,] <- c(0.05, 0.15, 0.30, 0.15, 0.20, 0.50, 0.30, 0.50, 0.60)
p_true_mat[16,] <- c(0.05, 0.10, 0.30, 0.10, 0.20, 0.40, 0.20, 0.30, 0.50)
p_true_mat[17,] <- c(0.10, 0.20, 0.40, 0.20, 0.30, 0.50, 0.30, 0.50, 0.60)
p_true_mat[18,] <- c(0.05, 0.15, 0.25, 0.15, 0.20, 0.30, 0.30, 0.50, 0.60)
p_true_mat[19,] <- c(0.05, 0.20, 0.25, 0.15, 0.25, 0.30, 0.25, 0.30, 0.40)
p_true_mat[20,] <- c(0.10, 0.20, 0.30, 0.20, 0.30, 0.50, 0.30, 0.40, 0.60)

#
# True MTC
#
comb_true <- list()
comb_true[[1]] <- 1
comb_true[[2]] <- 1
comb_true[[3]] <- 2
comb_true[[4]] <- 3
comb_true[[5]] <- 4
comb_true[[6]] <- 5
comb_true[[7]] <- 6
comb_true[[8]] <- 7
comb_true[[9]] <- 8
comb_true[[10]] <- 9
comb_true[[11]] <- c(2, 4)
comb_true[[12]] <- c(2, 7)
comb_true[[13]] <- c(3, 4)
comb_true[[14]] <- c(3, 5)
comb_true[[15]] <- c(3, 7)
comb_true[[16]] <- c(3, 8)
comb_true[[17]] <- c(5, 7)
comb_true[[18]] <- c(6, 7)
comb_true[[19]] <- c(6, 8)
comb_true[[20]] <- c(3, 5, 7)

# Define orderings
orders <- matrix(c(1:9,
                   1, 4, 7, 2, 5, 8, 3, 6, 9,
                   1, 2, 4, 3, 5, 7, 6, 8, 9,
                   1, 4, 2, 7, 5, 3, 8, 6, 9,
                   1, 2, 4, 7, 5, 3, 6, 8, 9,
                   1, 4, 2, 3, 5, 7, 8, 6, 9), ncol = k, byrow = T)
M <- nrow(orders)
pm <- rep(1/M, M)

#######################################################################################
#
#  POBLRM non-randomised setting
#
# setup the trail design
delta <- 0.3
cohortsize <- 3
ncohort <- 15
start.comb <- 1
start.d <- 1

# design parameters
p1 <- 0.15
nu <- 0.01
mu1 <- 1
mu2 <- -1
sigma1 <- 1
sigma2 <- 1
d <- StandardDose(p1=p1, nu=nu, k=k, mu1=mu1, mu2=mu2, sigma2=sigma2)
pseudo_KL <- c(0.4529323, 1.4953184, 0.5722332, 1.6481884)
skeleton <- seq(from=p1, by=nu, len=k)
p.skel <- matrix(nrow = M, ncol = k)
for (m in 1:M) {
  p.skel[m,] <- skeleton[order(orders[m,])]
}

B <- 10
PCS <- tox <- pt <- ptox <- matrix(nrow=C, ncol=k)
for (c in 1:C) {
  p_true <- p_true_mat[c,]
  PCS.temp <- tox.temp <- pt.temp <- ptox.temp <- matrix(nrow=B, ncol=k)
  for (b in 1:B) {
    sim <- POBLRM(d, orders, pm, mu1, mu2, sigma1, sigma2, pseudo_KL, p_true, delta, start.d=1, ncohort)
    PCS.temp[b,] <- sim$d.select
    tox.temp[b,] <- sim$tox.data
    pt.temp[b,] <- sim$pt.allocation
    ptox.temp[b,] <- sim$ptox.hat
  }
  PCS[c,] <- colSums(PCS.temp)
  tox[c,] <- colSums(tox.temp)
  pt[c,] <- colSums(pt.temp)
  ptox[c,] <- colSums(ptox.temp)
}

#############################################################################################
#
#  POBLRM simulation - randomised setting
#

# setup the trail design
delta <- 0.3
cohortsize <- 3
ncohort <- 15
start.d <- 1
start.comb <- 1
S <- 1e+4
burnin <- 2500
chains <- 1
# design parameters
p1 <- 0.10
nu <- 0.07
mu1 <- 1
mu2 <- -1
sigma1 <- 1
sigma2 <- 1
c2 <- 1
p0 <- 0.05
start.d <- 2

PCS <- tox <- pt <- ptox <- matrix(nrow=C, ncol=k)
for (c in 1:C) {
  p_true <- p_true_mat[c,]
  PCS.temp <- tox.temp <- pt.temp <- ptox.temp <- matrix(nrow=B, ncol=k)
  for (b in 1:B) {
    sim <- POBLRM(d, orders, pm, mu1, mu2, sigma1, sigma2, pseudo_KL, p_true, delta,
                  p0, c2, S, burnin, chains, start.d, ncohort, randomised=TRUE)
    PCS.temp[b,] <- sim$d.select
    tox.temp[b,] <- sim$tox.data
    pt.temp[b,] <- sim$pt.allocation
    ptox.temp[b,] <- sim$ptox.hat
  }
  PCS[c,] <- colSums(PCS.temp)
  tox[c,] <- colSums(tox.temp)
  pt[c,] <- colSums(pt.temp)
  ptox[c,] <- colSums(ptox.temp)
}

############################################################################################################
#
# Simulation under 2BLRM
#
# setup the trail 
#
delta <- 0.3                        # TTL
cohortsize <- 3                     # cohort size m
ncohort <- 15                       # number of cohorts N
n.stop <- cohortsize * ncohort      # no early stop
start.d <- c(1, 1)                  # start at the lowest combination
iter <- 2e+4                        # MCMC iterations

#
# Design parameters 
#
sigma_a1 <- sigma_b1 <- 2
sigma_a2 <- sigma_b2 <- 1
mu_zeta <- 0
sigma_zeta <- 1
p_a1 <- 0.1
p_b1 <- 0.1
nu_a <- 0.1
nu_b <- 0.1
mu_a1 <- log(p_a1/(1-p_a1))
mu_b1 <- log(p_b1/(1-p_b1))
mu_a2 <- 0
mu_b2 <- 0
a <- StandardDose(p1=p_a1, nu=nu_a, k=3, mu1=mu_a1, mu2=mu_a2, sigma2=sigma_a2)
b <- StandardDose(p1=p_b1, nu=nu_b, k=3, mu1=mu_b1, mu2=mu_b2, sigma2=sigma_b2)

B <- 10
PCS <- tox <- pt <- ptox <- matrix(nrow=C, ncol=k)
for (c in 1:C) {
  p_true <- p_true_mat[c,]
  p_true <- matrix(p_true, byrow = T, nrow=3)
  PCS.temp <- tox.temp <- pt.temp <- ptox.temp <- matrix(nrow=B, ncol=k)
  for (b.sim in 1:B) {
    sim <- BLRM2(p_true, a, b, mu_a1, mu_a2, mu_b1, mu_b2, mu_zeta,
                 sigma_a1, sigma_a2, sigma_b1, sigma_b2, sigma_zeta,
                 delta, start.d, ncohort, cohortsize, n.stop, iter, print=FALSE)
    PCS.temp[b.sim,] <- as.vector(t(sim$d.select))
    tox.temp[b.sim,] <- as.vector(t(sim$tox.data))
    pt.temp[b.sim,] <- as.vector(t(sim$pt.allocation))
    ptox.temp[b.sim,] <- as.vector(t(sim$ptox.hat))
  }
  PCS[c,] <- colSums(PCS.temp)
  tox[c,] <- colSums(tox.temp)
  pt[c,] <- colSums(pt.temp)
  ptox[c,] <- colSums(ptox.temp)
}
