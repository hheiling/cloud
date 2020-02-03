# Test whether convergence criteria change in fit_dat improves performance of glmmPen process
# Also use adaptive random walk in Metropolis-within-Gibbs during E steps in fit_dat

library(glmmPen)
library(stringr)

source("test_fit_dat.R")

N = 500
K = 5
seeds = 101:125
sd_options = c(0.5, 1.0)

# seeds = 1:2
# sd_options = c(0.5, 1.0)

M = 5000

s_l = length(seeds)
sd_l = length(sd_options)

data_sim = list()

for(sd_ in 1:length(sd_options)){
  for(s in 1:length(seeds)){
    dat =
      sim.data2(n = N, ptot = 2, pnonzero = 2, nstudies = K,
                sd_raneff = sd_options[sd_], family = 'binomial',
                slopes = T,
                seed = seeds[s], imbalance = 1, pnonzerovar = 0, beta = c(0, 2, 2))
    
    data_sim[[s+(sd_-1)*s_l]] = c(dat, seed = seeds[s], sd_ranef = sd_options[sd_])
  }
}

# for(i in 1:length(data_sim)){
#   print(data_sim[[i]]$seed)
#   print(data_sim[[i]]$sd_ranef)
# } 

fit_sim = function(dat, M, adapt_RW_options = adaptControl(batch_length = 100, 
                                                           burnin_batchnum = 1000,
                                                           offset = 8500)){
  
  set.seed(dat$seed)
  # Set gibbs = T
  fit_glmmPen = fit_dat_test(dat, lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 100, 
                        family = "binomial", trace = 0, penalty = "grMCP",
                        alpha = 1, nMC_max = 5000, t = 5,
                        returnMC = T, ufull = NULL, coeffull = NULL, gibbs = T, maxitEM = 75, 
                        ufullinit = NULL, MwG_sampler = "random_walk",
                        adapt_RW_options = adapt_RW_options) 
  
  coef = fit_glmmPen$coef
  sigma = fit_glmmPen$sigma
  u = fit_glmmPen$u
  # gibbs_accept_rate = fit_glmmPen$gibbs_accept_rate
  proposal_SD = fit_glmmPen$proposal_SD
  
  coef_record_all = fit_glmmPen$coef_record_all
  rownames(coef_record_all) = 1:nrow(coef_record_all)
  used_iter = apply(coef_record_all, 1, function(x) !any(is.na(x)))
  coef_record = coef_record_all[used_iter,]
  coef_record_tail = coef_record[(nrow(coef_record)-9):nrow(coef_record),]
  
  fit0_record = fit_glmmPen$fit0_record
  fit0_intercepts = fit0_record[used_iter,1]
  
  Znew2 = fit_glmmPen$extra$Znew2
  okindex = fit_glmmPen$extra$ok_index
  J = fit_glmmPen$J
  covgroup = fit_glmmPen$covgroup
  
  # Posterior draws taken using Adaptive Random Walk with Random Scan
  post_list = sample_mc_adapt(fit = fit_glmmPen$fit, cov = sigma, y = dat$y, X = dat$X,
                              Z = Znew2, nMC = M, family = "binomial",
                              group = dat$group, d = nlevels(dat$group), okindex = okindex,
                              nZ = Znew2, gibbs = T, uold = u, trace = 2,
                              proposal_SD = proposal_SD, batch = adapt_RW_options$burnin_batchnum,
                              batch_length = adapt_RW_options$batch_length, 
                              offset = adapt_RW_options$offset,
                              burnin_batchnum = adapt_RW_options$burnin_batchnum)
  
  post_U = post_list$u0
  gibbs_accept_rate = post_list$gibbs_accept_rate
  
  # Log-likelihood - glmer
  data = data.frame(y = dat$y, dat$X[,-1], group = dat$group)
  colnames(data) = c("y","X1","X2","group")
  
  fit_glmer = glmer(formula = y ~ X1 + X2 + (X1 + X2 | group), data = data, family = "binomial")
  
  return(list(dat = dat, coef = coef, sigma = sigma, coef_record_tail = coef_record_tail,
              fit0_intercepts = fit0_intercepts, 
              fit_extras = list(Znew2 = Znew2, covgroup = covgroup, J = J, okindex = okindex),
              post_U = post_U, 
              gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD,
              fit_glmer = fit_glmer)) 
}


output = list()

# for(i in 1:length(data_sim)){
#   print("Start dataset i")
#   output[[i]] = fit_sim(dat = data_sim[[i]], M = M)
# }

for(i in 1){
  cat("Start dataset ", i, "\n")
  output[[i]] = fit_sim(dat = data_sim[[i]], M = M)
}

sim_conv_crit = output

save(sim_conv_crit, file = "sim_conv_critera_output.RData")

