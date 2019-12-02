# Extension of test_adaptive_mcmc.R

# Updates to code: 
## proposal_var to proposal_SD 
##    (below: updated argument and output names in functions from test_adaptive_mcmc.R)
## Testing different increment specifications

library(remotes)

# Package version 1.4.4, from adaptive branch (created 11/26/2019)
install_github("hheiling/glmmPen", ref = "adaptive", force = TRUE)

library(glmmPen)
library(stringr)
library(ggplot2)
library(reshape2)

# Number of studies
K = 5 
d = 5
# Number of variables (Intercept, X1, X2)
q = 3

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
                         proposal_SD = matrix(1.0, nrow = K, ncol = q), batch = 0.0)
  
  post_U = post_list$u0
  
  gibbs_accept_rate = post_list$gibbs_accept_rate
  
  proposal_SD = post_list$proposal_SD
  
  return(list(post_U = post_U, gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD))
}

## Testing manual setting of proposal_SD in Metropolis-within-Gibbs functions
manual_MCMC_test = function(fit_glmmPen, dat, proposal_SD, M){
  
  # Test acceptance rates for new Adaptive Metropolis-within-gibbs
  post_list = sample.mc_test(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                             Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial", 
                             group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                             nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2,
                             proposal_SD)
  
  post_U = post_list$u0
  
  gibbs_accept_rate = post_list$gibbs_accept_rate
  
  proposal_SD = post_list$proposal_SD
  
  return(list(post_U = post_U, gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD))
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
  proposal_SD_vec = numeric(length(acc_rate_vec))
  for(i in 1:length(proposal_SD_vec)){
    if(acc_rate_vec[i] < 0.35){
      proposal_SD_vec[i] = 1.5
    }else if(acc_rate_vec[i] > 0.55){
      proposal_SD_vec[i] = 1/1.5
    }else{
      proposal_SD_vec[i] = 1.0
    }
  }
  
  proposal_SD_manual = matrix(proposal_SD_vec, nrow = K, ncol = q, byrow = F)
  
  # Output: post_U, gibbs_accept_rate, proposal_SD
  manual = manual_MCMC_test(fit_glmmPen, dat, proposal_SD_manual, M)
  # manual_plots = mcmc_diagnostics(manual$post_U)
  
  # Ouput: post_U, gibbs_accept_rate, proposal_SD
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

# Re-do test using updates to code in version 1.4.4

fit_sim2 = function(output, M = 10^4){
  
  K = 5
  q = 3
  
  dat = output$dat
  fit_glmmPen = output$fit_glmmPen
  
  original = sample.mc2(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                        Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial", 
                        group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                        nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2)
  
  post_U_original = original$u0
  accept_rate_original = original$gibbs_accept_rate
  # original_plots = mcmc_diagnostics(post_U_original)
  
  original_list = list(post_U = post_U_original, gibbs_accept_rate = accept_rate_original)
  
  # acc_rate_vec = c(accept_rate_original)
  # proposal_SD_vec = numeric(length(acc_rate_vec))
  # for(i in 1:length(proposal_SD_vec)){
  #   if(acc_rate_vec[i] < 0.35){
  #     proposal_SD_vec[i] = 1.5
  #   }else if(acc_rate_vec[i] > 0.55){
  #     proposal_SD_vec[i] = 1/1.5
  #   }else{
  #     proposal_SD_vec[i] = 1.0
  #   }
  # }
  # 
  # proposal_SD_manual = matrix(proposal_SD_vec, nrow = K, ncol = q, byrow = F)
  
  # Output: post_U, gibbs_accept_rate, proposal_SD
  # manual = manual_MCMC_test(fit_glmmPen, dat, proposal_SD_manual, M)
  # manual_plots = mcmc_diagnostics(manual$post_U)
  
  # Ouput: post_U, gibbs_accept_rate, proposal_SD
  adaptive = AMCMC_test(fit_glmmPen, dat, M)
  # adaptive_plots = mcmc_diagnostics(adaptive$post_U)
  
  # plots = list(original = original_plots, manual = manual_plots,
  #              adaptive = adaptive_plots)
  
  # keep = list(dat = dat, fit_glmmPen = fit_glmmPen, original = original_list, 
  #             manual = manual, adaptive = adaptive, plots = plots)
  keep = list(dat = dat, fit_glmmPen = fit_glmmPen, original = original_list, 
              adaptive = adaptive)
  
  return(keep)
}

# load test_results_noplots object
load("adaptive_vs_set_mcmc_noplots.RData")

test_results_2 = list()
for(i in 3){
  test_results_2[[i]] = fit_sim2(test_results_noplots[[i]], M = 12500)
}

x = c(test_results_2[[3]]$original$gibbs_accept_rate)
y = c(test_results_2[[3]]$adaptive$gibbs_accept_rate)
plot(x, y, ylim = c(0,1), xlim = c(0,1))
abline(a = 0, b = 1)

# save(test_results_2, file = "adaptive_mcmc_test2.RData")

plots_test = list()
for(i in 3){
  plots_test[[i]] = mcmc_diagnostics(test_results_2[[i]]$adaptive$post_U)
}

plots_test[[3]]$sample_path
plots_test[[3]]$autocorr

######################
# Random code experimentation / exploration
batch_full = seq(from = 500, to = 5*10^4, by = 500)
increment = 1 / sqrt(8000 + batch_full)
delta = pmin(0.01, increment)

# Worst case scenario:
sd_cum = 1*exp(cumsum(delta))
log_sd = log(sd_cum)


######################
