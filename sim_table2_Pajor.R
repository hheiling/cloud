
# Comparison of Pajor IS CAME logLik calculation method with glmer results

library(glmmPen)
library(parallel)
library(foreach)
library(doParallel)
library(stringr)

N = 500
K = 5
seeds = 1:100
sd_options = c(0.5, 1.0, 2.0)

M = c(10^4,10^5,3*10^5)
M_str = c("10^4","10^5","3*10^5")

s_l = length(seeds)
sd_l = length(sd_options)

data_sim_Tab2 = list()

for(sd_ in 1:length(sd_options)){
  for(s in seeds){
    dat =
      sim.data2(n = N, ptot = 2, pnonzero = 2, nstudies = K,
                sd_raneff = sd_options[sd_], family = 'binomial',
                slopes = T,
                seed = s, imbalance = 1, pnonzerovar = 0, beta = c(0, 1, 1))
    data_sim_Tab2[[s+(sd_-1)*s_l]] = c(dat, seed = s)
  }
}

fit_sim = function(dat, M){
  
  set.seed(dat$seed)
  # Set gibbs = T
  fit_glmmPen = fit_dat(dat, lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 100, 
                        family = "binomial", trace = 0, penalty = "grMCP",
                        alpha = 1, nMC_max = 4000, 
                        returnMC = T, ufull = NULL, coeffull = NULL, gibbs = T, maxitEM = 150, 
                        ufullinit = NULL,
                        c = 1, M = 10^3) 
                # These c and M  are unrelated to the Pajor method under consideration
  
  # Log-likelihood - Pajor Method
  ## Set gibbs = T to test thinning of posterior
  post_list = sample.mc2(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                         Z = fit_glmmPen$extra$Znew2, nMC = 10^4, family = "binomial", 
                         group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                         nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u)
  
  post_U = post_list$u0
  
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
  
  # Log-likelihood - glmm ?
  
  # ll_pdif_vec = rbind(ll_pdif_glmer, ll_pdif_glmm)
  
  return(matrix(ll_pdif, nrow = 1))
}


numCores = detectCores()
registerDoParallel(numCores)
output = 
  foreach(i = 1:length(data_sim_Tab2), .combine = rbind) %dopar% {
    library(glmmPen)
    fit_sim(data_sim_Tab2[[i]], M = M)
  }

M_label = str_c("M_", M_str)
colnames(output) = M_label

row_labels = c(rep(str_c("sd_", sd_options[1]), times = s_l),
               rep(str_c("sd_", sd_options[2]), times = s_l),
               rep(str_c("sd_", sd_options[3]), times = s_l))
rownames(output) = row_labels

save(output, file = "sim_table2_Pajor_output.RData")

# save(output)

