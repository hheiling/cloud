# Diagnostics for sim_table2_Pajor.R

# Question of interest:
# Why it there such a discrepancy between the glmer and glmmPen results in some cases?

# Comparison of Pajor IS CAME logLik calculation method with glmer results

library(glmmPen)
library(parallel)
library(foreach)
library(doParallel)
library(stringr)

N = 500
K = 5
seeds = 1:100
sd_options = c(0.5, 1.0)

M = c(10^4,10^5,3*10^5)
M_str = c("10^4","10^5","3*10^5")
M_label = str_c("M_", M_str)

s_l = length(seeds)
sd_l = length(sd_options)

data_diagnostics = list()

for(sd_ in 1:length(sd_options)){
  for(s in seeds){
    dat =
      sim.data2(n = N, ptot = 2, pnonzero = 2, nstudies = K,
                sd_raneff = sd_options[sd_], family = 'binomial',
                slopes = T,
                seed = s, imbalance = 1, pnonzerovar = 0, beta = c(0, 1, 1))
    data_diagnostics[[s+(sd_-1)*s_l]] = c(dat, seed = s, sd_ranef = sd_)
  }
}

fit_sim = function(dat, M){
  
  set.seed(dat$seed)
  # Set gibbs = T
  fit_glmmPen = fit_dat(dat, lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 100, 
                        family = "binomial", trace = 0, penalty = "grMCP",
                        alpha = 1, nMC_max = 4000, 
                        returnMC = T, ufull = NULL, coeffull = NULL, gibbs = T, maxitEM = 150, 
                        ufullinit = NULL) 
  
  # Log-likelihood - Pajor Method
  ## Set gibbs = T to test thinning of posterior
  post_list = sample.mc2(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                         Z = fit_glmmPen$extra$Znew2, nMC = 10^4, family = "binomial", 
                         group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                         nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2)
  
  post_U = post_list$u0
  
  ll_pajor = numeric(length(M))
  
  for(m in 1:length(M)){
    ll_pajor[m] = CAME_IS(posterior = post_U, y = dat$y, X = dat$X, Z = dat$Z, 
                          group = dat$group, coef = fit_glmmPen$coef, sigma = fit_glmmPen$sigma, 
                          family = "binomial", M = M[m])
  }
  
  names(ll_pajor) = M_label
  
  # fit_dat diagnostic output:
  p = ncol(dat$X)
  beta_glmmPen = fit_glmmPen$coef[1:p]
  modes_glmmPen = matrix(colMeans(post_U), nrow = K, byrow = T)
  sigma_glmmPen = fit_glmmPen$sigma
  gibbs_accept_rate = post_list$gibbs_accept_rate
    colnames(gibbs_accept_rate) = c("(Intercept)","X1","X2")
    rownames(gibbs_accept_rate) = c("grp1","grp2","grp3","grp4","grp5")
  
  # Log-likelihood - glmer
  data = data.frame(y = dat$y, dat$X[,-1], group = dat$group)
  colnames(data) = c("y","X1","X2","group")
  
  fit_glmer = glmer(formula = y ~ X1 + X2 + (X1 + X2 | group), data = data, family = "binomial")
  
  ll_glmer = logLik(fit_glmer)
  
  ll_pdif = (ll_pajor - ll_glmer)/ll_glmer
  names(ll_pdif) = M_label
  
  # glmer diagnostic output:
  
  beta_glmer = fixef(fit_glmer)
    names(beta_glmmPen) = names(beta_glmer)
  modes_glmer = ranef(fit_glmer)[[1]]
    colnames(modes_glmmPen) = colnames(modes_glmer)
    rownames(modes_glmmPen) = rownames(modes_glmer)
  sigma_glmer = cov(modes_glmer)
  
  return(list(dat = dat, beta_glmmPen = beta_glmmPen, beta_glmer = beta_glmer,
              modes_glmmPen = modes_glmmPen, modes_glmer = modes_glmer,
              sigma_glmer = sigma_glmer, sigma_glmmPen = sigma_glmmPen,
              gibbs_accept_rate = gibbs_accept_rate,
              ll_glmer = ll_glmer, ll_pajor = ll_pajor, ll_pdif = ll_pdif))
}

# In parallel:
numCores = 4
registerDoParallel(numCores)
output =
  foreach(i = 1:length(data_diagnostics)) %dopar% {
    library(glmmPen)
    fit_sim(data_diagnostics[[i]], M = M)
  }

# In serial:
# output =
#   foreach(i = 1:length(data_diagnostics)) %do% {
#     library(glmmPen)
#     fit_sim(data_diagnostics[[i]], M = M)
#   }

save(output, file = "glmer_Pajor_diagnostics.RData")


