
library(glmmPen)

setupLambda_pglmm <- function(X, y, family, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)
  
  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    fit <- glm(y~X[, -ind], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }
  if (family=="gaussian") {
    zmax <- .Call("maxprod", X, fit$residuals, ind, penalty.factor) / n
  } else {
    zmax <- .Call("maxprod", X, residuals(fit, "working") * fit$weights, ind, penalty.factor) / n
  }
  lambda.max <- zmax/alpha
  
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  } else {
    # lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
    lambda <- exp(seq(log(lambda.max+0.05),log(lambda.min*lambda.max),len=nlambda))
  }
  
  if (length(ind)!=p) lambda[1] <- lambda[1] * 1.000001
  
  return(lambda)
}

LambdaRange = function(X, y, family, alpha = 1, lambda.min = NULL, nlambda = 40, 
                       penalty.factor = NULL){
  n = nrow(X)
  p = ncol(X)
  
  # X_std = std(X)
  
  if(family == "gaussian"){
    yy = y = mean(y)
  }else{
    yy = y
  }
  
  if(is.null(lambda.min)){
    lambda.min = ifelse(n>p, 0.001, 0.05)
  }
  
  if(is.null(penalty.factor)){
    penalty.factor = rep(1, p)
  }
  
  lambda = setupLambda_pglmm(X, yy, family, alpha, lambda.min, nlambda, penalty.factor)
  
  return(lambda)
  
}


N = 500
K = 5
seeds = 1:100
sd_options = c(1.0, 2.0)

s_l = length(seeds)
sd_l = length(sd_options)

data_sim_Tab4 = list()

for(stddev in 1:length(sd_options)){
  for(s in seeds){
    data_init =
      sim.data2(n = N, ptot = 10, pnonzero = 2, nstudies = K,
                sd_raneff = sd_options[stddev], family = 'binomial',
                slopes = T,
                seed = s, imbalance = 1, pnonzerovar = 0, beta = c(0, 1, 1))
    data_sim_Tab4[[s+(stddev-1)*s_l]] = c(data_init, seed = s)
  }
}


select_sim = function(dat, const = 1, M = 10^5, nlambda = 20){
  # Calculate lambda0 and lambda1 ranges
  range_fixed = LambdaRange(dat$X[,-1], dat$y, family = "binomial", nlambda = nlambda)
  range_random = range_fixed
  
  opt_res = matrix(0, nrow = 2, ncol = 9)
  opt_fef = matrix(0, nrow = 2, ncol = ncol(dat$X))

  set.seed(dat$seed)
  out_fit = select_tune(dat, nMC = 100, lambda0_range = range_fixed,
                        lambda1_range = range_random, family = "binomial",
                        penalty = "grMCP", returnMC = T,
                        conv = 0.001, nMC_max = 3000, trace = 0, ufull = NULL, coeffull = NULL,
                        gibbs = T, maxitEM = 50, alpha = 1,
                        c = const, M = M)


  opt_res[1,] = out_fit$result[which.min(out_fit$result[,9]),]
  opt_res[2,] = out_fit$result[which.min(out_fit$result[,3]),]

  opt_fef[1,] = out_fit$coef[which.min(out_fit$result[,9]),1:ncol(dat$X)]
  opt_fef[2,] = out_fit$coef[which.min(out_fit$result[,3]),1:ncol(dat$X)]

  rownames(opt_res) = c("BIC","BICh")
  rownames(opt_fef) = rownames(opt_res)

  return(list(opt_res = opt_res, opt_fef = opt_fef))

}

library(parallel)
library(doParallel)
numCores = detectCores()

output = mclapply(data_sim_Tab4, FUN = select_sim, mc.cores = numCores)

save(output, file = "sim_table4_output.RData")

# Calculated desired results: optimized by BIC
B1_sum = numeric(sd_l)
B2_sum = numeric(sd_l)
TP_sum = numeric(sd_l)
FP_sum = numeric(sd_l)

for(j in 1:sd_l){
  for(i in (1+(j-1)*s_l):(j*s_l)){
    B1_sum[j] = B1_sum[j] + output[[i]]$opt_fef[1,2]
    B2_sum[j] = B2_sum[j] + output[[i]]$opt_fef[1,3]
    TP_sum[j] = TP_sum[j] + sum(output[[i]]$opt_fef[1,c(2,3)] != 0)
    FP_sum[j] = FP_sum[j] + sum(output[[i]]$opt_fef[1,-c(2,3)] != 0)
  }
}

BIC_output = list(B1_est = B1_sum/s_l, B2_est = B2_sum/s_l, TP = TP_sum/s_l, FP = FP_sum/s_l)
save(BIC_output, file = "sim_table4_BIC_output.RData")

# Calculated desired results: optimized by BICh (hybrid)
B1_sum = numeric(sd_l)
B2_sum = numeric(sd_l)
TP_sum = numeric(sd_l)
FP_sum = numeric(sd_l)

for(j in 1:sd_l){
  for(i in (1+(j-1)*s_l):(j*s_l)){
    B1_sum[j] = B1_sum[j] + output[[i]]$opt_fef[2,2]
    B2_sum[j] = B2_sum[j] + output[[i]]$opt_fef[2,3]
    TP_sum[j] = TP_sum[j] + sum(output[[i]]$opt_fef[2,c(2,3)] != 0)
    FP_sum[j] = FP_sum[j] + sum(output[[i]]$opt_fef[2,-c(2,3)] != 0)
  }
}

BICh_output = list(B1_est = B1_sum/s_l, B2_est = B2_sum/s_l, TP = TP_sum/s_l, FP = FP_sum/s_l)
save(BICh_output, file = "sim_table4_BICh_output.RData")