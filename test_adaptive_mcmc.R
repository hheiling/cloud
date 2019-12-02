

# Examining adaptive MCMC results

# Objectives: 
## Identify datasets from glmer_vs_Pajor_diagnostics.R (glmer_Pajor_diagnostics.RData) that
## have generally low acceptance rates from the Metropolis-within-Gibbs sampling approach initially
## implemented in sample.mc2 from the subset of datasets with sd_ranef = 1.0

## Using these data-sets, re-run fit_dat and now perform the adaptive Metropolis-within-Gibbs 
## algorithm (sample.mc3)

## Output the following: fit_dat object, autocorrelation and sample.path plots from the u output
## by the fit_dat function, the sample.mc3 object (adaptive Metropolis-within-Gibbs object),
## and autocorrelation and sample.path plots from this sample.mc3 object

# glmmmPen version 1.4.3, before changing proposal_var name to proposal_SD name

library(glmmPen)

library(glmmPen)
library(stringr)
library(ggplot2)
library(reshape2)

# Number of studies
K = 5 
d = 5
# Number of variables (Intercept, X1, X2)
q = 3

## Load simulation results
load("glmer_Pajor_diagnostics.RData")
output_1 = output[c(101:200)]

# Acceptance rates for sd_ranef = 1.0
acc_rate2 = matrix(0, nrow = 100, ncol = d*q)

vars_q = rep(c("(Intercept)","X1","X2"), each = d)
grp = rep(1:5, times = q)
colnames(acc_rate2) = str_c(vars_q, ":", grp)

# Convert matrix of acceptance rates into a row vector for each dataset
for(i in 1:100){
  gibbs_acc_rate = output_1[[i]]$gibbs_accept_rate
  rate_vec = c(gibbs_acc_rate[,1], gibbs_acc_rate[,2], gibbs_acc_rate[,3])
  acc_rate2[i,] = matrix(rate_vec, nrow = 1)
}

full_accept_rates_sd1 = acc_rate2

# Pick datasets where the original gibbs acceptance rate is generally low
low_acc_rate = function(vector){
  (sum(vector < 0.25) >= 6)
}

low_rate_results = apply(full_accept_rates_sd1, 1, FUN = low_acc_rate)

sum(low_rate_results) # Number of datasets with generally low gibbs acceptance rates

output_low = output_1[low_rate_results]
length(output_low)

## Code to get fit_dat results
fit_dat_out = function(dat){
  
  set.seed(dat$seed)
  # Set gibbs = T
  fit_glmmPen = fit_dat(dat, lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 100,
                        family = "binomial", trace = 0, penalty = "grMCP",
                        alpha = 1, nMC_max = 4000,
                        returnMC = T, ufull = NULL, coeffull = NULL, gibbs = T, maxitEM = 150,
                        ufullinit = NULL)
  
  return(fit_glmmPen)
}

## Testing Adaptive Metropolis-within-Gibbs functions
AMCMC_test = function(fit_glmmPen, dat, M){
  
  # Set variables
  K = 5
  q = 3
  
  # Test acceptance rates for new Adaptive Metropolis-within-gibbs
  post_list = sample.mc3(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                         Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial", 
                         group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                         nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2,
                         proposal_var = matrix(1.0, nrow = K, ncol = q), batch = 0.0)
  
  post_U = post_list$u0
  
  gibbs_accept_rate = post_list$gibbs_accept_rate
  
  proposal_var = post_list$proposal_var
  
  return(list(post_U = post_U, gibbs_accept_rate = gibbs_accept_rate, proposal_var = proposal_var))
}

## Testing manual setting of proposal_var in Metropolis-within-Gibbs functions
manual_MCMC_test = function(fit_glmmPen, dat, proposal_var, M){
  
  # Test acceptance rates for new Adaptive Metropolis-within-gibbs
  post_list = sample.mc_test(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                             Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial", 
                             group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                             nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2,
                             proposal_var)
  
  post_U = post_list$u0
  
  gibbs_accept_rate = post_list$gibbs_accept_rate
  
  proposal_var = post_list$proposal_var
  
  return(list(post_U = post_U, gibbs_accept_rate = gibbs_accept_rate, proposal_var = proposal_var))
}


## Evaluate the performance of the chain
mcmc_diagnostics = function(post_U){ # modified version of plot_mcmc.pglmmObj
  
  d = 5
  var_num = 3
  
  type = c("sample.path","histogram","cumsum","autocorr")
  
  U_keep = post_U
  var_names = c("Intercept","X1","X2")
  grp_names = str_c("grp", 1:5)
  var_str = rep(var_names, each = d)
  grp_str = rep(grp_names, times = q)
  U_cols = str_c(var_str, ":", grp_str)
  
  vars = "all"
  grps = "all"
  
  U_t = data.frame(U_keep, t = 1:nrow(U_keep))
  colnames(U_t) = c(U_cols, "t")
  U_long = melt(U_t, id = "t")
  U_plot = data.frame(U_long, var_names = rep(var_names, each = d*nrow(U_keep)),
                      grp_names = rep(rep(grp_names, each = nrow(U_keep)), times = var_num))
  
  plots_return = list()
  
  if("sample.path" %in% type){
    plot_sp = ggplot(U_plot, mapping = aes(x = t, y = value)) + geom_path() +
      facet_grid(var_names ~ grp_names) + xlab("iteration t") + ylab("draws")
    
    plots_return$sample_path = plot_sp
  }
  if("histogram" %in% type){
    hist_U = ggplot(U_plot) + geom_histogram(mapping = aes(x = value)) + 
      facet_grid(var_names ~ grp_names) + xlab("draws")
    
    plots_return$histogram = hist_U
  }
  if("cumsum" %in% type){
    U_means = colMeans(U_keep)
    U_means = data.frame(rbind(U_means))[rep.int(1L, nrow(U_keep)), , drop = FALSE]
    U_tmeans = apply(U_keep, 2, cumsum) / 1:nrow(U_keep)
    U_tdiff = U_tmeans - U_means
    U_cumsum = apply(U_tdiff, 2, cumsum)
    U_t = data.frame(U_cumsum, t = 1:nrow(U_cumsum))
    colnames(U_t) = c(colnames(U_keep), "t")
    U_long = melt(U_t, id = "t")
    U_plot = data.frame(U_long, var_names = rep(var_names, each = d*nrow(U_keep)),
                        grp_names = rep(rep(grp_names, each = nrow(U_keep)), times = var_num)) 
    plot_cumsum = ggplot(U_plot) + geom_smooth(mapping = aes(x = t, y = value), color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(var_names ~ grp_names) + xlab("iteration t") + ylab("Cumulative Sum")
    
    plots_return$cumsum = plot_cumsum
  }
  if("autocorr" %in% type){
    grp_index = rep(grp_names, times = var_num)
    var_index = rep(var_names, each = d)
    for(j in 1:ncol(U_keep)){
      ACF = acf(U_keep[,j], plot=F, lag.max = 40)
      ACF_df = with(ACF, data.frame(lag,acf))
      ACF_df$grp_names = grp_index[j]
      ACF_df$var_names = var_index[j]
      if(j == 1){
        ACF_all = ACF_df
      }else{
        ACF_all = rbind(ACF_all, ACF_df)
      }
    }
    
    plot_acf = ggplot(data = ACF_all, mapping = aes(x = lag, y = acf)) +
      geom_hline(mapping = aes(yintercept = 0)) + 
      geom_segment(mapping = aes(xend = lag, yend = 0)) +
      facet_grid(var_names ~ grp_names)
    
    plots_return$autocorr = plot_acf
  }
  
  return(plots_return)
  
}

## Plots: sample_path, histogram, cumsum, autocorr
# test_plots = mcmc_diagnostics(post_U = test_obj1$post_U)

fit_sim = function(output, M = 2*10^4){
  
  K = 5
  q = 3
  
  dat = output$dat
  fit_glmmPen = fit_dat_out(dat)
  
  original = sample.mc2(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                        Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial", 
                        group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                        nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2)
  
  post_U_original = original$u0
  accept_rate_original = original$gibbs_accept_rate
  # original_plots = mcmc_diagnostics(post_U_original)
  
  original_list = list(post_U = post_U_original, gibbs_accept_rate = accept_rate_original)
  
  acc_rate_vec = c(accept_rate_original)
  proposal_var_vec = numeric(length(acc_rate_vec))
  for(i in 1:length(proposal_var_vec)){
    if(acc_rate_vec[i] < 0.35){
      proposal_var_vec[i] = 1.5
    }else if(acc_rate_vec[i] > 0.55){
      proposal_var_vec[i] = 1/1.5
    }else{
      proposal_var_vec[i] = 1.0
    }
  }
  
  proposal_var_manual = matrix(proposal_var_vec, nrow = K, ncol = q, byrow = F)
  
  # Output: post_U, gibbs_accept_rate, proposal_var
  manual = manual_MCMC_test(fit_glmmPen, dat, proposal_var_manual, M)
  # manual_plots = mcmc_diagnostics(manual$post_U)
  
  # Ouput: post_U, gibbs_accept_rate, proposal_var
  adaptive = AMCMC_test(fit_glmmPen, dat, M)
  # adaptive_plots = mcmc_diagnostics(adaptive$post_U)
  
  # plots = list(original = original_plots, manual = manual_plots,
  #              adaptive = adaptive_plots)
  
  # keep = list(dat = dat, fit_glmmPen = fit_glmmPen, original = original_list, 
  #             manual = manual, adaptive = adaptive, plots = plots)
  keep = list(dat = dat, fit_glmmPen = fit_glmmPen, original = original_list, 
              manual = manual, adaptive = adaptive)
  
  return(keep)
}

# test_sim = fit_sim(output_low[[1]])

test_results = list()
for(i in 1:length(output_low)){
  test_results[[i]] = fit_sim(output_low[[i]], M = 2*10^4)
}

test_results_noplots = test_results

save(test_results_noplots, file = "adaptive_vs_set_mcmc_noplots.RData")


test_results_abbrev = list()
for(i in 1:7){
  original = test_results[[i]]$original
  manual = test_results[[i]]$manual
  adaptive = test_results[[i]]$adaptive

  test_results_abbrev[[i]] = list(original = original, manual = manual, adaptive = adaptive)
  
}

save(test_results_abbrev, file = "adaptive_vs_set_mcmc_abbrev.RData")

# load("adaptive_vs_set_mcmc_abbrev.RData")
# test_results_abbrev2 = list()
# for(i in 1:7){
#   original = list(post_U = test_results_abbrev[[i]]$original$post_U,
#                   gibbs_accept_rate = test_results_abbrev[[i]]$original$gibbs_accept_rate)
#   # original$post_U = original$post_U[1:10^4,]
#   manual = list(post_U = test_results_abbrev[[i]]$manual$post_U,
#                 gibbs_accept_rate = test_results_abbrev[[i]]$manual$gibbs_accept_rate,
#                 proposal_var = test_results_abbrev[[i]]$manual$proposal_var)
#   # manual$post_U = manual$post_U[1:10^4,]
#   adaptive = list(post_U = test_results_abbrev[[i]]$adaptive$post_U,
#                   gibbs_accept_rate = test_results_abbrev[[i]]$adaptive$gibbs_accept_rate,
#                   proposal_var = test_results_abbrev[[i]]$adaptive$proposal_var)
#   # adaptive$post_U = adaptive$post_U[1:10^4,]
#   
#   test_results_abbrev2[[i]] = list(original = original, manual = manual, adaptive = adaptive)
# }
# 
# save(test_results_abbrev2, file = "adaptive_vs_set_mcmc_abbrev2.RData")

