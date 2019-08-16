N = 500
K = 5
seeds = 1:100
sd_options = c(1.0, 2.0)
const = seq(from = 1, to = 3, by = 0.1)
M = c(10^4,10^5,10^6)
M_str = c("10^4","10^5","10^6")

s_l = length(seeds)
sd_l = length(sd_options)

data_sim_Tab2 = list()


for(stddev in 1:length(sd_options)){
  for(s in seeds){
    dat =
      sim.data2(n = N, ptot = 2, pnonzero = 2, nstudies = K,
                sd_raneff = sd_options[stddev], family = 'binomial',
                slopes = T,
                seed = s, imbalance = 1, pnonzerovar = 0, beta = c(0, 1, 1))
    data_sim_Tab2[[s+(stddev-1)*s_l]] = c(dat, seed = s)
  }
}

fit_sim = function(dat, const, M){
  
  ll_glmmPen = matrix(0, nrow = length(const), ncol = length(M))
  
  set.seed(dat$seed)
  fit_glmmPen = fit_dat(dat, lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 100, 
                        family = "binomial", trace = 0, penalty = "grMCP",
                        alpha = 1, nMC_max = 3000, 
                        returnMC = T, ufull = NULL, coeffull = NULL, gibbs = F, maxitEM = 100, 
                        ufullinit = NULL,
                        c = const[1], M = M[1])
  # Log-likelihood
  ll_glmmPen[1,1] = fit_glmmPen$ll
  
  set.seed(dat$seed)
  for(m in 1:length(M)){
    for(c in 1:length(const)){
      if(c == 1 && m == 1) next
      
      # Log-likelihood
      ll_glmmPen[c,m] = logLik_imp(dat$y, dat$X, dat$Z, U = fit_glmmPen$u, 
                                   sigma = fit_glmmPen$sigma, 
                                   dat$group, fit_glmmPen$coef, fit_glmmPen$J, family = "binomial", 
                                   df = 10, c, m)
      
    }
  }
  
  data = data.frame(y = dat$y, dat$X[,-1], group = dat$group)
  fit_glmer = glmer(formula = y ~ X1 + X2 + (X1 + X2 | group), data = data, family = "binomial")
  
  ll_glmer = logLik(fit_glmer)
  
  ll_pdif = (ll_glmmPen - ll_glmer)/ll_glmer
  
  if(ncol(ll_pdif) == 1){
    ll_pdif_vec = ll_pdif
  }else{
    ll_pdif_vec = ll_pdif[,1]
    for(j in 2:ncol(ll_pdif)){
      ll_pdif_vec = rbind(ll_pdif_vec, ll_pdif[,j])
    }
  }
  
  return(matrix(ll_pdif_vec, nrow = 1))
}

library(glmmPen)
library(parallel)
library(foreach)
library(doParallel)

numCores = detectCores()
registerDoParallel(numCores)
output = 
  foreach(i = 1:length(data_sim_Tab2), .combine = rbind) %dopar% {
    library(glmmPen)
    fit_sim(data_sim_Tab2[[i]], const = const, M = M)
  }

library(stringr)
M_label = rep(str_c("M:", M_str), each = length(const))
colnames(output) = str_c("c:", const, ",", M_label)

# save(output, file = "sim_table2_output.RData")

# save(output)

