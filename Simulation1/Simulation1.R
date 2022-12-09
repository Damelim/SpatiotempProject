source('G:/My Drive/spatiotemp_project/functions.R')

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


#### Step 1) Set model parameters and setting for simulations ####

theta <- 100
alpha <- c(-8, 0.1)
gamma <- c(0.03, 0.02)
beta <- 0.6
phi <- 0.1
d <- c(0.5, 0.5) ; # Location of Gaussian kernel for initial conditions
q <- 100^2 # Number of spatial grid cells


grid = expand.grid(seq(0,1,length = 100), seq(0,1,length = 100)) # generate covariates
x_grid = grid$Var1 ; y_grid = grid$Var2

w <- exp(-0.5*((9*(x_grid-0.2))^2 + (5*(y_grid-0.2))^2)) + exp(-0.5*((4*(x_grid-0.8))^2 + (8*(y_grid-0.8))^2)) # length q
z <- 0.5*(x_grid-0.5)^3 + (y_grid-0.5)^3 # runif(q, -0.5, 0.5) # length q

xmn = 0 ; xmx = 1 ; ymn = 0 ; ymx = 1

mu     <- raster( nrows = q^0.5, ncols = q^0.5, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
mu[]     <- exp(model.matrix(~z) %*% alpha) # model.matrix(~z) : q by 2 matrix

lambda <- raster(nrows = q^0.5, ncols = q^0.5, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
lambda[] <- model.matrix(~w) %*% gamma # model.matrix(~w) : q by 2 matrix


png('mu.png')
plot(mu, main = expression(paste('Diffusion rate ', mu(s) > 0)))
dev.off()

png('lambda.png')
plot(lambda, main = expression(paste('Growth rate ',lambda(s))))
dev.off()



us.fact <- 5 # Scaling factor

mu.bar <- aggregate(mu, fact = us.fact, na.rm = TRUE, fun = function(x, na.rm) { # Diffusion rate for homogenized pde
  1/mean(1/x)
})
mu.bar #; plot(mu.bar)

lambda.bar <- mu.bar * aggregate(lambda/mu, fact = us.fact, na.rm = TRUE, FUN = mean) # Growth rate for homogenized pde
lambda.bar #; plot(lambda.bar)

NN <- neighborhood(mu.bar) ; dim(NN) ; NN[1:3,] # First-order neighborhood matrix

 
H <- propagator.plain(NN = NN, 
                      mu = mu.bar[], 
                      lambda = lambda.bar[], 
                      dx = us.fact/q^0.5,
                      dy = us.fact/q^0.5, 
                      dt = 1)  # Propagator matrix

dim(H) ; c(min(H), max(H))



u0 <- raster(  nrows = q^0.5, ncols = q^0.5, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
D <- rdist(data.frame(SpatialPoints(u0)), matrix(d, 1, 2))
u0[] <- (exp(-D^2/phi^2)/sum(exp(-D^2/phi^2))) * theta
c0 <- raster(  nrows = q^0.5/us.fact, ncols = q^0.5/us.fact, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
c0[] <- extract(mu * u0, SpatialPoints(mu.bar))



#Make raster layers
c1 <- c0; c2 <- c0; c3 <- c0; c4 <- c0; c5 <- c0; c6 <- c0; c7 <- c0 ; c8 <- c0; c9 <- c0; c10 <- c0
c11 <- c0; c12 <- c0; c13 <- c0; c14 <- c0 ; c15 <- c0; c16 <- c0; c17 <- c0; c18 <- c0; c19 <- c0; c20 <- c0
c21 <- c0; c22 <- c0; c23 <- c0; c24 <- c0 ; c25 <- c0; c26 <- c0; c27 <- c0; c28 <- c0; c29 <- c0; c30 <- c0


#Fill raster layers
c1[] <- H%*%c0[]
c2[] <- H%*%c1[]
c3[] <- H%*%c2[]
c4[] <- H%*%c3[]
c5[] <- H%*%c4[]
c6[] <- H%*%c5[]
c7[] <- H%*%c6[]
c8[] <- H%*%c7[]
c9[] <- H%*%c8[]
c10[] <- H%*%c9[]
c11[] <- H%*%c10[]
c12[] <- H%*%c11[]
c13[] <- H%*%c12[]
c14[] <- H%*%c13[]
c15[] <- H%*%c14[]
c16[] <- H%*%c15[]
c17[] <- H%*%c16[]
c18[] <- H%*%c17[]
c19[] <- H%*%c18[]
c20[] <- H%*%c19[]  
c21[] <- H%*%c20[]
c22[] <- H%*%c21[]
c23[] <- H%*%c22[]
c24[] <- H%*%c23[]
c25[] <- H%*%c24[]
c26[] <- H%*%c25[]
c27[] <- H%*%c26[]
c28[] <- H%*%c27[]
c29[] <- H%*%c28[]
c30[] <- H%*%c29[] ; length(c30[]) ; class(c30[])

c.all <- stack(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
               c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
               c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) #Stack raster layers

length_time = 30

c.all <- calc.c(H,c0,length_time,1:length_time)



color <- colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1)
png('simulated_plain_diffusion_process_c(s,t).png')
levelplot(c.all, cuts = 254, margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:length_time), 
          col.regions = colorRampPalette(rev(brewer.pal(11,"Spectral")),bias = 1),
          main = "Simulated plain diffusion process c(s,t)")
dev.off()


# Calculate u(s,t)
u.all <- disaggregate(c.all, us.fact)/mu
png('simulated_homogenized_ecological_diffusion_process_u(s,t).png')
levelplot(u.all, cuts = 254, margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:length_time), col.regions = colorRampPalette(rev(brewer.pal(11,"Spectral")), bias = 1),
          main = "Simulated homogenized ecological diffusion process u(s,t)")
dev.off()


x <- runif(length(vec(u.all[])), -0.5, 0.5) ; class(x) ; summary(x) ; length(x)



png('True_probability_of_binaryevent.png')
p <- link(exp(x * beta) * vec(u.all[]))
p.all <- u.all
p.all[] <- p
levelplot(p.all, cuts = 254, margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:length_time), col.regions = colorRampPalette(rev(brewer.pal(11,"Spectral")), bias = 1),
          main = "True probabilities p(s,t)")
dev.off()


png('True_binaryevent.png')
y <- rbinom(length(p), 1, p)
y.all <- u.all
y.all[] <- y

totalpoints_vec = rep(NA, length_time)
for(i in 1:length(totalpoints_vec)){
  mn = 10000 * (i-1) + 1
  mx = 10000 * i
  totalpoints_vec[i] = sum(y[mn:mx])
}


levelplot(y.all, margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:length_time, '  :  ', totalpoints_vec ,sep = ""), col.regions = c("black","white"), cuts = 1, colorkey = FALSE,
          main = "True event (white = 1, black = 0)")
dev.off()



df1 <- data.frame(y = vec(y.all[]), 
                  x = x, 
                  z = rep(z, length_time), 
                  w = rep(w, length_time),
                  t = rep(1:length_time, each = q), 
                  s1 = data.frame(SpatialPoints(y.all))[, 1],
                  s2 = data.frame(SpatialPoints(y.all))[, 2], 
                  cell = rep(1:q, length_time)
                  ) 
names(df1) 

# df1$t    : time indicator
# df1$cell : 1~10000, 1~10000, .... : cell indicator

df1$s1[1:10]
df1$s1[9990:10000]

df1$s1[10001:10010]
















#### 2. MCMC setting and run MCMC ####

set.seed(4223)
K <- sample(1:dim(df1)[1], replace = FALSE, size = 50000) ; class(K) ; length(K)
df2 <- df1[K, ] ; dim(df2)

spatial.covariates <- raster(  nrows = sqrt(q), ncols = sqrt(q), xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA)
spatial.covariates$cell <- 1:q
spatial.covariates$intercept <- 1
spatial.covariates$z <- z # numeric
spatial.covariates$w <- w # numeric


set.seed(1) ; n.iter = 20000 ; thin.u = 10
ptm <- proc.time()
chain1 <-fit.model.mcmc(
  n.iter = n.iter, print.iter=FALSE,
  thin.u = thin.u, save.u = TRUE,
  alpha.start = alpha + rnorm(2, 0, sd = 0.2),
  beta.start = beta + rnorm(1,0,0.2),
  gamma.start = gamma + rnorm(2, 0, sd = 0.2),  
  theta.start = theta + rnorm(1,0,0.2),  
  phi.start = phi + rnorm(1, 0, 0.01), 
  alpha.prior.var = 10^6, beta.prior.var = 10^6, gamma.prior.var = 10^6,
  theta.prior = 10^6, phi.prior = 10^6,
  alpha.tune = c(1/500,1/1500), beta.tune = 1/40, gamma.tune = c(1/10^5,1/600), theta.tune = 2, phi.tune = 1/500,
  # alpha.tune = c(1/600,1/3000), beta.tune = 1/60, gamma.tune = c(1.4/10^5,1/440), theta.tune = 2, phi.tune = 1/500,
  y = df2$y, X = model.matrix(~x-1,data=df2), K = K,
  t.stop = length_time, t.keep = seq(1,length_time,by=1),
  dt = 1, # time step
  spatial.covariates = spatial.covariates,
  diffusion.coef = c(0,1,1,0), growth.coef = c(0,1,0,1),
  d = SpatialPoints(cbind(0.5,0.5)),
  us.fact = 5
)
runtime = proc.time() - ptm


burn.in <- n.iter/2
mcmc.parms <- as.mcmc(do.call(cbind, chain1[-7])[-c(1:burn.in), ])
summary(mcmc.parms) ; effectiveSize(mcmc.parms)


parms <- 11
keep.sample <- seq(thin.u, n.iter, by = thin.u)
temp <- cbind(do.call(cbind, chain1[-7])[keep.sample, ], do.call(cbind, chain1[7])[-1, ])
temp <- temp[which(burn.in < keep.sample), ]


u.est <- raster(  nrows = q^0.5, ncols = q^0.5, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
u.est[] <- 1
u.est <- stack(u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est,
               u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est,
               u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est, u.est)
names(u.est) <- paste("t =", 1:length_time)
u.est[] <- colMeans((temp[, -c(1:parms)]))

u.est[["layer.1"]] ; u.est[["layer.30"]]

estimated <- levelplot(u.est, main = list(expression("E(u(" * bold(s) * "," * italic(t) * ")|" * bold(y) * ")"), cex = 3), 
                       cuts = 254, margin = FALSE, scales = list(draw = FALSE),
                       names.attr = paste("t =", 1:length_time), col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))


truth <- levelplot(u.all, main = list(expression("Truth (u(" * bold(s) * "," * italic(t) * "))"), cex = 3), 
                   cuts = 254, margin = FALSE, scales = list(draw = FALSE), names.attr = paste("t =", 1:length_time),
                   col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))   







#### 3. Estimated probabilities ####

p.est <- raster(nrows = q^0.5, ncols = q^0.5, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
p.est[] <- 1

p.est <- stack(p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est,
               p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est,
               p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est, p.est)
names(p.est) <- paste("t =", 1:length_time)
probability_mat = apply((t(temp[, -c(1:parms)]) * exp(df1$x * temp[, 7])), 1, link)
postmean = colMeans(probability_mat)
p.est[] <- postmean
p.est[["layer.1"]] ; p.est[["layer.30"]]


estimated <- levelplot(p.est, main = list(expression("E(" * italic(p) * "|" * bold(y) * ")"), cex = 1.2), cuts = 254, margin = FALSE, scales = list(draw = FALSE), 
                       names.attr = paste("t =", 1:length_time), col.regions = colorRampPalette(rev(brewer.pal(11,"Spectral")), bias = 1))

png('Estimated_probability_of_binaryevent.png')
estimated
dev.off()




#### 4. LOWER and UPPER QUANTILES of estimated probabilities ####

lower_quantile = apply(probability_mat, MARGIN = 2, FUN = function(x){quantile(x, 0.025)})
upper_quantile = apply(probability_mat, MARGIN = 2, FUN = function(x){quantile(x, 0.975)})

# sanity check
all( lower_quantile[299000:300000] <=  upper_quantile[299000:300000])
all( lower_quantile[299000:300000] <= postmean[299000:300000])
all( upper_quantile[299000:300000] >= postmean[299000:300000])



png('lower.png')

lower <- raster(nrows = q^0.5, ncols = q^0.5, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
lower[] <- 1

lower <- stack(lower, lower, lower, lower, lower, lower, lower, lower, lower, lower,
               lower, lower, lower, lower, lower, lower, lower, lower, lower, lower,
               lower, lower, lower, lower, lower, lower, lower, lower, lower, lower)
names(lower) <- paste("t =", 1:length_time)

lower[] <- lower_quantile
lower[["layer.1"]] ; lower[["layer.30"]]

lower_plot <- levelplot(lower, main = expression(paste("0.025 quantile : ",pi,"(p|y)")), cuts = 254, margin = FALSE, scales = list(draw = FALSE), 
                       names.attr = paste("t =", 1:length_time), col.regions = colorRampPalette(rev(brewer.pal(11,"Spectral")), bias = 1))


lower_plot
dev.off()


png('upper.png')

upper <- raster(nrows = q^0.5, ncols = q^0.5, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = NA)
upper[] <- 1

upper <- stack(upper, upper, upper, upper, upper, upper, upper, upper, upper, upper,
               upper, upper, upper, upper, upper, upper, upper, upper, upper, upper,
               upper, upper, upper, upper, upper, upper, upper, upper, upper, upper)
names(upper) <- paste("t =", 1:length_time)

upper[] <- upper_quantile
upper[["layer.1"]] ; upper[["layer.30"]]

upper_plot <- levelplot(upper, main = expression(paste("0.975 quantile : ",pi,"(p|y)")), cuts = 254, margin = FALSE, scales = list(draw = FALSE), 
                        names.attr = paste("t =", 1:length_time), col.regions = colorRampPalette(rev(brewer.pal(11,"Spectral")), bias = 1))

upper_plot
dev.off()








#### 5. posterior predictive (by plug in for the sake of time) ####

### originally, posterior predictive check is done by iterating 1) sampling from posterior dist'n and 2) generating Y replicate given the sampled paramter sample.
### I instead plug in the posterior mean as the parameter update and generate new data for the sake of time


png('Predicted_binaryevent.png')
p_est = colMeans(apply((t(temp[, -c(1:parms)]) * exp(df1$x * temp[, 7])), 1, link))
y_pred <- rbinom(n = length(p_est), size = 1, prob = p_est)
y_pred.all <- u.all
y_pred.all[] <- y_pred

totalpoints_pred_vec = rep(NA, length_time)
for(i in 1:length(totalpoints_pred_vec)){
  mn = 10000 * (i-1) + 1
  mx = 10000 * i
  totalpoints_pred_vec[i] = sum(y_pred[mn:mx])
}

levelplot(y_pred.all, margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:length_time, '  :  ', totalpoints_pred_vec ,sep = ""), col.regions = c("black","white"), cuts = 1, colorkey = FALSE,
          main = "Posterior predictive Y (white = 1, black = 0)")
dev.off()





#### 6. acceptance probability of each parameter ####

acceptance_info = mcmc.parms[,1:4] 
acceptance_rate = colMeans(acceptance_info) ; acceptance_rate


#### 7. traceplot of each parameter ####
mcmc.parms = mcmc.parms[,5:11] 

png('traceplot.png')
# matrix(c(1,2,3,3,4,5,6,7),4,2, byrow =T)
layout(matrix(c(1,2,3,3,4,5,6,7),4,2, byrow =T)) # ; layout.show(7)
ts.plot(mcmc.parms[,1], main = "intercept_z", ylab = "value") #; abline(v = burn.in, col = "red")
ts.plot(mcmc.parms[,2], main = "z", ylab = "value") # ; abline(v = burn.in, col = "red")
ts.plot(mcmc.parms[,3], main = "x", ylab = "value") # ; abline(v = burn.in, col = "red")
ts.plot(mcmc.parms[,4], main = "intercept_w", ylab = "value")  #  ; abline(v = burn.in, col = "red")
ts.plot(mcmc.parms[,5], main = "w", ylab = "value")   #; abline(v = burn.in, col = "red")
ts.plot(mcmc.parms[,6], main = "phi", ylab = "value")   #; abline(v = burn.in, col = "red")
ts.plot(mcmc.parms[,7], main = "theta", ylab = "value")   #; abline(v = burn.in, col = "red")
dev.off()


#### 8. Goodness of fit ####

p_est[p_est == 0] = 1e-10 # to avoid NaN
crossentropy = cross.entropy(p, p_est) ; paste('Cross entropy loss = ', crossentropy)

MAE = mean(abs(p_est - p))
accuracy = mean(y == y_pred) ; paste('Accuracy = ',accuracy)



#### 9. save all the results ####

save(p, p_est, y, y_pred, totalpoints_vec, mcmc.parms, acceptance_rate, MAE, crossentropy, accuracy, runtime, file = "result.RData")



load('result.rdata')
MAE
crossentropy
accuracy


colMeans(as.matrix(mcmc.parms))
c(alpha, beta, gamma, phi, theta)
