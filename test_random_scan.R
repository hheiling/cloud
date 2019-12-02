# Testing of adaptive random walk Metropolis-within-Gibbs code with random scan incorporated
# Extention of test_random_walk.R

library(remotes)

# Package version 1.4.4.1, from adaptive branch (created 11/27/2019)
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
AMCMC_test = function(fit_glmmPen, dat, M, batch_length, offset){
  
  # Set variables
  K = 5
  q = 3
  
  # U = matrix(rnorm(K*q), nrow = 1)
  
  # Test acceptance rates for new Adaptive Metropolis-within-gibbs
  post_list = sample_mc_rw_rs(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                           Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial",
                           group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                           nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2,
                           proposal_SD = matrix(1.0, nrow = K, ncol = q), batch = 0.0,
                           batch_length = batch_length, offset = offset)
  
  # post_list = sample_mc_rw_rs(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
  #                             Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial", 
  #                             group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
  #                             nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = U, trace = 2,
  #                             proposal_SD = matrix(1.0, nrow = K, ncol = q), batch = 0.0,
  #                             batch_length = batch_length, offset = offset)
  
  post_U = post_list$u0
  
  gibbs_accept_rate = post_list$gibbs_accept_rate
  
  proposal_SD = post_list$proposal_SD
  
  return(list(post_U = post_U, gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD))
}


## Evaluate the performance of the chain
mcmc_diagnostics = function(post_U, dat){ # modified version of plot_mcmc.pglmmObj
  
  # glmer modes:
  data = data.frame(y = dat$y, X1 = dat$X[,2], X2 = dat$X[,3], group = dat$group)
  fit_glmer = glmer(formula = y ~ X1 + X2 + (1 + X1 + X2 | group), data = data, family = "binomial")
  
  mode_df = ranef(fit_glmer)[[1]]
  modes = c(as.matrix(mode_df))
  
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
  
  modes_data = data.frame(modes = modes, var_names = var_str, grp_names = grp_str)
  print(modes_data)
  
  plots_return = list()
  
  if("sample.path" %in% type){
    plot_sp = ggplot(U_plot, mapping = aes(x = t, y = value)) + geom_path() +
      facet_grid(var_names ~ grp_names) + xlab("iteration t") + ylab("draws") +
      geom_hline(data = modes_data, aes(yintercept = modes), color = "red")
    
    plots_return$sample_path = plot_sp
  }
  if("histogram" %in% type){
    hist_U = ggplot(U_plot) + geom_histogram(mapping = aes(x = value)) + 
      facet_grid(var_names ~ grp_names) + xlab("draws") +
      geom_vline(data = modes_data, aes(xintercept = modes), color = "red")
    
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


# Test of random walk code

fit_sim_rw_rs = function(output, M = 10^4, batch_length = 500, offset = 5000){
  
  K = 5
  q = 3
  
  dat = output$dat
  fit_glmmPen = output$fit_glmmPen
  
  # original = sample.mc2(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
  #                       Z = fit_glmmPen$extra$Znew2, nMC = M, family = "binomial", 
  #                       group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
  #                       nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2)
  # 
  # post_U_original = original$u0
  # accept_rate_original = original$gibbs_accept_rate
  # 
  # original_list = list(post_U = post_U_original, gibbs_accept_rate = accept_rate_original)
  
  original_list = output$original
  
  
  # Ouput: post_U, gibbs_accept_rate, proposal_SD
  random_walk = AMCMC_test(fit_glmmPen, dat, M, batch_length, offset)
  
  keep = list(dat = dat, fit_glmmPen = fit_glmmPen, original = original_list, 
              random_walk = random_walk)
  
  return(keep)
}

# load test_results_noplots object
load("adaptive_vs_set_mcmc_noplots.RData")

test_randscan = list()

set.seed(146)
for(i in 4:6){
  test_randscan[[i]] = fit_sim_rw_rs(test_results_noplots[[i]], M = 2000,
                            batch_length = 250, offset = 0)
}

par(mfrow = c(1,3))
for(element in 4:6){
  x = c(test_randscan[[element]]$original$gibbs_accept_rate)
  y = c(test_randscan[[element]]$random_walk$gibbs_accept_rate)
  plot(x, y, ylim = c(0,1), xlim = c(0,1))
  abline(a = 0, b = 1)
  abline(a = 0.44, b = 0, col = "red")
}


plot_test = mcmc_diagnostics(post_U = test_randscan[[5]]$random_walk$post_U,
                             dat = test_randscan[[5]]$dat)

plot_test$sample_path

# The End