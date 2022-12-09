library(raster)
library(fields)
library(RColorBrewer)
library(rasterVis)
library(truncnorm)
library(mvtnorm)
library(gridExtra)
library(fBasics)
library(coda)



## Cross entropy loss calculation ##

cross.entropy <- function(p, phat){  # from https://rpubs.com/juanhklopper/cross_entropy
  x <- 0
  for (i in 1:length(p)){
    x <- x + (p[i] * log(phat[i]))
  }
  return(-x)
}

## link function g(.) ##

link <- function(nu){ptruncnorm(nu, a=0, b=Inf, mean = 0, sd = 1)} # Link function g(.)


## obtains first order neighborhood of a RasterLayer object ##

neighborhood <- function(raster){  
  nn <- matrix(NA, length(raster[]), 4)
  for(i in 1:dim(nn)[1]){
    loc <- adjacent(raster,i)[,2]
    ln <- loc[which((loc+1)==i)]
    rn <- loc[which((loc-1)==i)]
    bn <- loc[which((loc-dim(raster)[2])==i)]
    tn <- loc[which((loc+dim(raster)[2])==i)]
    nn[i,1] <- if(length(ln)>0){ln}else{0}
    nn[i,2] <- if(length(rn)>0){rn}else{0}
    nn[i,3] <- if(length(bn)>0){bn}else{0}
    nn[i,4] <- if(length(tn)>0){tn}else{0}
  }
  nn
}



## Propagator matrix for diffusion PDE ##
propagator.plain <- function(NN, mu, lambda, dx, dy, dt){  
  
  H <- matrix(0,dim(NN)[1],dim(NN)[1])
  for(i in 1:dim(H)[1]){
    if(length(which(NN[i,]>0))==4){
      H[i,i] <- 1-2*mu[i]*(dt/dx^2 + dt/dy^2) + dt*lambda[i]
      H[i,NN[i,1]] <- dt/dx^2*mu[i]
      H[i,NN[i,2]] <- dt/dx^2*mu[i]
      H[i,NN[i,3]] <- dt/dy^2*mu[i]
      H[i,NN[i,4]] <- dt/dy^2*mu[i]}
  }
  H
}


## RasterStack of c(s,t) ##
calc.c <- function(H,c0,t.steps,t.keep){
  c.all <- c0
  c.all[] <- H%*%c0[]
  c.all <- stack(mget(rep("c.all",t.steps)))
  for(t in 2:t.steps){
    c.all[[t]][] <- H%*%c.all[[t-1]][]
  }
  c.all[[t.keep]]
}


## Bayesian parameter update ##

fit.model.mcmc <- function(n.iter,print.iter,thin.u,save.u,
                           alpha.start,beta.start,gamma.start,theta.start,phi.start,
                           alpha.prior.var,beta.prior.var,gamma.prior.var,theta.prior,phi.prior,
                           alpha.tune,beta.tune,gamma.tune,theta.tune,phi.tune,
                           y,X,K,t.stop,t.keep,dt,spatial.covariates,diffusion.coef,growth.coef,d,
                           us.fact){
  
  # Link function
  link <- function(nu){ptruncnorm(nu, a=0, b=Inf, mean = 0, sd = 1)}
  # Broad-scale grid cells
  us.cells <- aggregate(spatial.covariates$cell,fact=us.fact,FUN=mean)
  us.cells[] <- 1:length(us.cells[])
  us.res <- res(us.cells)
  # Fine-scale grid cells
  ds.cells <- spatial.covariates$cell
  # First-order neighborhood matrix
  NN <- neighborhood(us.cells)
  # Distance matrix for initial state
  D <- rdist(data.frame(SpatialPoints(spatial.covariates)),data.frame(d))
  # Spatial covariates z(s) and w(s)
  Z <- as.matrix(spatial.covariates[[which(diffusion.coef==1)]])
  W <- as.matrix(spatial.covariates[[which(growth.coef==1)]])
  # Variables that will be saved
  theta <- matrix(NA,n.iter+1,1)
  phi <- matrix(NA,n.iter+1,1)
  alpha <- matrix(NA,n.iter+1,sum(diffusion.coef))
  beta <- matrix(NA,n.iter+1,dim(X)[2])
  gamma <- matrix(NA,n.iter+1,sum(growth.coef))
  accept <- matrix(NA,n.iter+1,4)
  
  
  if(save.u==TRUE){
    u.save <- matrix(NA,n.iter/thin.u+1, max(spatial.covariates$cell[])*length(t.keep));
    keep.u <- seq(thin.u,n.iter,by=thin.u)
  }
  
  #### Initial values setting
  theta[1,] <- theta.start
  phi[1,] <- phi.start
  alpha[1,] <- alpha.start
  beta[1,] <- beta.start
  gamma[1,] <- gamma.start
  
  colnames(theta) <- c("theta")
  colnames(phi) <- c("phi")
  colnames(alpha) <- names(spatial.covariates)[which(diffusion.coef==1)]
  colnames(beta) <- colnames(X)
  colnames(gamma) <- names(spatial.covariates)[which(growth.coef==1)]
  colnames(accept) <- c("accept.rate.theta.phi","accept.rate.alpha", "accept.rate.beta","accept.rate.gamma")
  
  
  
  
  ### Begin MCMC loop ###
  
  for(i in 1:n.iter){
    
    if(i %% 50 == 0){cat('iter = ',i,'\n')}
    
    #### Sample u(s,t) on fine-scale grid
    
    nu <- exp(X%*%beta[i,])
    mu <- exp(Z%*%alpha[i,])
    lambda <- W%*%gamma[i,]
    ds.cells$mu <- mu
    ds.cells$lambda <- lambda/mu
    mu.bar <- aggregate(ds.cells$mu,fact=us.fact, fun=function(x,na.rm){(1/mean(1/x,na.rm=TRUE))})
    lambda.bar <- mu.bar*aggregate(ds.cells$lambda,fact=us.fact,
                                   fun=function(x,na.rm){(mean(x,na.rm=TRUE))})
    H <- propagator.plain(NN,mu.bar[],lambda.bar[],us.res[1],us.res[2],dt=dt)
    ds.cells$u0 <- (exp(-D^2/phi[i,]^2)/sum(exp(-D^2/phi[i,]^2))*theta[i,])
    us.cells$c0 <- extract(ds.cells$mu*ds.cells$u0, SpatialPoints(us.cells))
    c.all <- calc.c(H,us.cells$c0,t.stop,t.keep)
    u.all <- disaggregate(c.all,us.fact)/ds.cells$mu
    u <- vec(u.all[])
    if(length(which(keep.u==i))==1 & save.u==TRUE){u.save[which(keep.u==i)+1,] <- u}
    u <- u[K]
    
    
    
    
    #### Sample phi and theta
    
    # print('phi and theta')
    phi.star <- rnorm(1,phi[i,1],phi.tune) # propose phi_star
    theta.star <- rnorm(1,theta[i,],theta.tune) # propose theta_star
    
    if(min(phi.star,theta.star)>0){
      
      ds.cells$u0.star <- (exp(-D^2/phi.star^2)/sum(exp(-D^2/phi.star^2))*theta.star)
      us.cells$c0.star <- extract(ds.cells$mu*ds.cells$u0.star, SpatialPoints(us.cells))
      c.all.star <- calc.c(H,us.cells$c0.star,t.stop,t.keep)
      u.all.star <- disaggregate(c.all.star,us.fact)/ds.cells$mu
      u.star <- vec(u.all.star[])[K]
      
      mh1 <- sum(dbinom(y,1,link(nu*u.star),log=TRUE)) +
        log(dtruncnorm(theta.star, mean=0, sd=theta.prior^0.5, a=0, b=Inf)) +
        log(dtruncnorm(phi.star, mean=0, sd=phi.prior^0.5, a=0, b=Inf))
      mh2 <- sum(dbinom(y,1,link(nu*u),log=TRUE)) +
        log(dtruncnorm(theta[i,], mean=0, sd=theta.prior^0.5, a=0, b=Inf)) +
        log(dtruncnorm(phi[i,], mean=0, sd=phi.prior^0.5, a=0, b=Inf))
      mh <- exp(mh1-mh2)
    }
    else{
      mh=0
    }
    # acceptance probability calculation
    # print(mh)
    if(runif(1) < mh){ # accept
      phi[i+1,] <- phi.star;
      theta[i+1,] <- theta.star;
      us.cells$c0 <- us.cells$c0.star;
      accept[i+1,1] <- 1; 
      u <- u.star
    } else{ # rejection
      phi[i+1,] <- phi[i,];
      theta[i+1,] <- theta[i,];
      accept[i+1,1] <- 0
    }
    
    
    
    #### Sample alpha
    # print('alpha')
    alpha.star <- rmvnorm(1,alpha[i,],diag(alpha.tune,dim(Z)[2])) # propose alpha from Multivariate Normal
    mu.star <- exp(Z%*%t(alpha.star))
    lambda <- W%*%gamma[i,]
    ds.cells$mu <- mu.star
    ds.cells$lambda <- lambda/mu.star
    mu.bar <- aggregate(ds.cells$mu,fact=us.fact,
                        fun=function(x,na.rm){(1/mean(1/x,na.rm=TRUE))})
    lambda.bar <- mu.bar*aggregate(ds.cells$lambda,fact=us.fact,
                                   fun=function(x,na.rm){(mean(x,na.rm=TRUE))})
    H.star <- propagator.plain(NN,mu.bar[],lambda.bar[],us.res[1],us.res[2],dt=dt)
    c.all.star <- calc.c(H.star,us.cells$c0,t.stop,t.keep)
    u.all.star <- disaggregate(c.all.star,us.fact)/ds.cells$mu
    u.star <- vec(u.all.star[])[K]
    # acceptance probability calculation
    mh1 <- sum(dbinom(y,1,link(nu*u.star),log=TRUE)) + sum(dnorm(alpha.star,0,alpha.prior.var^0.5,log=TRUE))
    mh2 <- sum(dbinom(y,1,link(nu*u),log=TRUE)) +  sum(dnorm(alpha[i,],0,alpha.prior.var^0.5,log=TRUE))
    mh <- exp(mh1-mh2)
    
    if(runif(1) < mh){ # acceptance
      alpha[i+1,] <- alpha.star;
      accept[i+1,2] <- 1;
      u <- u.star
    } else{ # rejection
      alpha[i+1,] <- alpha[i,];
      accept[i+1,2] <- 0
    }
    
    
    # Sample beta
    # print('beta')
    beta.star <- rmvnorm(1,beta[i,],diag(beta.tune,dim(X)[2]))
    nu.star <- exp(X%*%t(beta.star))
    mh1 <- sum(dbinom(y,1,link(nu.star*u),log=TRUE)) +
      sum(dnorm(beta.star,0,beta.prior.var^0.5,log=TRUE))
    mh2 <- sum(dbinom(y,1,link(nu*u),log=TRUE)) +
      sum(dnorm(beta[i,],0,beta.prior.var^0.5,log=TRUE))
    mh <- exp(mh1-mh2)
    if(runif(1) < mh){
      beta[i+1,] <- beta.star;
      accept[i+1,3] <- 1;
      nu <- nu.star} 
    else{
      beta[i+1,] <- beta[i,];
      accept[i+1,3] <- 0
    }
    
    
    #### Sample gamma
    # print('gamma')
    gamma.star <- rmvnorm(1,gamma[i,], diag(gamma.tune,dim(W)[2]))
    lambda.star <- W%*%t(gamma.star)
    mu <- exp(Z%*%alpha[i+1,])
    ds.cells$mu <- mu
    ds.cells$lambda <- lambda.star/mu
    mu.bar <- aggregate(ds.cells$mu,fact=us.fact, fun=function(x,na.rm){(1/mean(1/x,na.rm=TRUE))})
    lambda.bar <- mu.bar*aggregate(ds.cells$lambda,fact=us.fact, fun=function(x,na.rm){(mean(x,na.rm=TRUE))})
    H.star <- propagator.plain(NN,mu.bar[],lambda.bar[],us.res[1],us.res[2],dt=dt)
    c.all.star <- calc.c(H.star,us.cells$c0,t.stop,t.keep)
    u.all.star <- disaggregate(c.all.star,us.fact)/ds.cells$mu
    u.star <- vec(u.all.star[])[K]
    #acceptance probability calculation
    mh1 <- sum(dbinom(y,1,link(nu*u.star),log=TRUE)) + sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
    mh2 <- sum(dbinom(y,1,link(nu*u),log=TRUE)) + sum(dnorm(gamma[i,],0,gamma.prior.var^0.5,log=TRUE))
    mh <- exp(mh1-mh2)
    
    if(runif(1) < mh){ # acceptance
      gamma[i+1,] <- gamma.star;
      accept[i+1,4] <- 1
    }else{ # rejection
      gamma[i+1,] <- gamma[i,];
      accept[i+1,4] <- 0
    }
    if(print.iter==TRUE){if(i %% 1==0) print(i)}
  }
  if(save.u==FALSE){u.save="NA"}
  list(accept = accept, alpha = alpha, beta = beta, gamma = gamma, phi = phi, theta = theta, u = u.save)
}



