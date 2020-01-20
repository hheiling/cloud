# Test whether convergence criteria change in fit_dat improves performance of glmmPen process
# Also use adaptive random walk in Metropolis-within-Gibbs during E steps in fit_dat

library(glmmPen)
library(parallel)
library(foreach)
library(doParallel)
library(stringr)

N = 500
K = 5
seeds = 1:7
sd_options = 1.0

M = 10^4
M_str = c("10^4")

s_l = length(seeds)
sd_l = length(sd_options)

data_sim = list()

for(sd_ in 1:length(sd_options)){
  for(s in seeds){
    dat =
      sim.data2(n = N, ptot = 2, pnonzero = 2, nstudies = K,
                sd_raneff = sd_options[sd_], family = 'binomial',
                slopes = T,
                seed = s, imbalance = 1, pnonzerovar = 0, beta = c(0, 1, 1))
    data_sim[[s+(sd_-1)*s_l]] = c(dat, seed = s)
  }
}

fit_sim = function(dat, M){
  
  set.seed(dat$seed)
  # Set gibbs = T
  fit_glmmPen = fit_dat(dat, lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 100, 
                        family = "binomial", trace = 0, penalty = "grMCP",
                        alpha = 1, nMC_max = 5000, t = 7,
                        returnMC = T, ufull = NULL, coeffull = NULL, gibbs = T, maxitEM = 100, 
                        ufullinit = NULL, 
                        adapt_RW_options = adaptControl(batch_length = 100,
                                                        burnin_batchnum = 500,
                                                        offset = 9000)) 
  
  # Log-likelihood - Pajor Method
  # Posterior draws taken using Adaptive Random Walk with Random Scan
  post_list = sample_mc_adapt(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                              Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial",
                              group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                              nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2,
                              proposal_SD = matrix(1.0, nrow = K, ncol = q), batch = 0.0,
                              batch_length = batch_length, offset = offset, burnin = burnin)
  
  post_U = post_list$u0
  gibbs_accept_rate = post_list$gibbs_accept_rate
  proposal_SD = post_list$proposal_SD
  
  ll_pajor = numeric(length(M))
  
  for(m in 1:length(M)){
    ll_pajor[m] = CAME_IS(posterior = post_U, y = dat$y, X = dat$X, Z = dat$Z, 
                          group = dat$group, coef = fit_glmmPen$coef, sigma = fit_glmmPen$sigma, 
                          family = "binomial", M = M[m])
  }
  
  # Log-likelihood - glmer
  data = data.frame(y = dat$y, dat$X[,-1], group = dat$group)
  colnames(data) = c("y","X1","X2","group")
  
  fit_glmer = glmer(formula = y ~ X1 + X2 + (X1 + X2 | group), data = data, family = "binomial")
  
  ll_glmer = logLik(fit_glmer)
  
  ll_pdif = (ll_pajor - ll_glmer)/ll_glmer
  
  return(list(dat = dat, ll_pajor = ll_pajor, ll_glmer = ll_glmer, ll_pdif = ll_pdif, 
              post_U = post_U, gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD,
              fit_glmer = fit_glmer)) 
}


output =
  foreach(i = 1:length(data_sim)) %do% {
    library(glmmPen)
    fit_sim(data_sim[[i]], M = M)
  }

save(output, file = "sim_conv_critera_output.RData")

