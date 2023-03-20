library(tidyverse)
library(readr)
library(glmnet)
library(cowplot)

set.seed(1)

###################################################################
# functions to get first, second, and third order regression fits # 
# each takes a dataframe of presence absence data, fits the model #
# to the specified order                                          #
###################################################################

#helper function
get_folds <- function(df, unique_communities, n_folds = 10){
  
  n_points <- length(unique_communities) 
  
  base_fold_size <- floor(n_points/n_folds)
  rem_fold_size <- floor(n_points %% n_folds)
  diff <- n_folds - rem_fold_size
  
  if (diff != n_folds){
    full_sample <- c(rep(seq(1:diff), base_fold_size), rep((diff + 1): n_folds, base_fold_size + 1))
  } else {
    full_sample <- rep(seq(1:n_folds), base_fold_size)
  }
  
  fold_sample <- sample(full_sample, size = length(full_sample), replace = FALSE)
  fold_data <- data.frame(community = unique_communities, fold_id = fold_sample)
  
  fold_id_data <- left_join(df, fold_data, by = 'community')
  fold_ids <- fold_id_data$fold_id
  
  return(fold_ids)
}

###### First order regression and leave-one-out ######
get_first_order_fit <- function(df){
  
  #first pull out replicates 
  communities <- sapply(1:nrow(df), function(i) paste0(df[i,] %>% dplyr::select(-fitness), collapse = '-'))
  df <- df %>% mutate(community = communities)
  unique_communities <- unique(communities) 
  
  #set up regression
  y <- as.matrix(df$fitness)   
  f <- as.formula(y ~ .)
  x <- model.matrix(f, df %>% dplyr::select(all_of(all_species)))
  
  fold_ids <- get_folds(df, unique_communities, n_folds = 10)
  
  cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
  
  predicted_y <- predict(cv_fit, x, s = "lambda.1se") 
  
  #plot(cv_fit, sign.lambda = 1)
  
  res_cv <- data.frame(exp = df$community, observed = y, predicted = array(predicted_y)) %>% 
    group_by(exp) %>% mutate(mean_observed = mean(observed))
  res_cv_unique <- res_cv %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  r2 <- cor(res_cv_unique$mean_observed, res_cv_unique$predicted)^2
  
  maxlim <- max(c(res_cv$mean_observed, res_cv$predicted))*1.01
  minlim <- min(c(res_cv$mean_observed, res_cv$predicted))
  minlim <- ifelse(minlim > 0, minlim*.99, minlim*1.01) #check for 0 condition
  
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  p1 <- res_cv %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + 
    #geom_point(aes(x = observed, y = predicted), alpha = .2) + #this line for replicates
    geom_abline() + 
    theme_bw() + xlab('Observed') + ylab('Predicted')  + #ggtitle('Regression to 2nd order') + 
    xlim(minlim, maxlim) + ylim(minlim, maxlim) +
    annotate('text', x = .9*max(res_cv$observed), y = 1.1*min(res_cv$predicted), label = mylabel)
  
  #save infit RData 
  save(res_cv, file = paste0("../../Results/model_fit_plots/first_order_infit_", name, ".RData")) 
  
  ##### basic linear model with loo ######
  loo_cv_res <- tibble()
  for (i in 1:length(unique_communities)){
    #get experiments to drop
    exp <- unique_communities[i]
    exp_index <- which(df$community == exp)
    
    y <- as.matrix(df[-exp_index,]$fitness)
    f <- as.formula(y ~ .)
    x <- model.matrix(f, df[-exp_index,] %>% dplyr::select(all_of(all_species)))
    
    exp_index <- which(df$community == exp )
    fold_ids <- get_folds(df[-exp_index,], unique_communities[-i], n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_index[1],]$fitness)
    f_test <- as.formula(y_test ~ .)
    x_test <- model.matrix(f_test, df[exp_index[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se")
    
    for (j in 1:length(exp_index)){
      tmp <- list(exp = exp, observed = df[exp_index[j],]$fitness, predicted = array(y_pred_cv))
      loo_cv_res <- rbind(loo_cv_res, tmp)
    }
  }
  
  loo_cv_res <- data.frame(loo_cv_res) %>% group_by(exp) %>% mutate(mean_observed = mean(observed))
  loo_res_unique <- loo_cv_res %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  #r2 on mean observed value across replicates
  r2_loo <- cor(loo_res_unique$mean_observed, loo_res_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2_loo, digits = 3)))
  p2 <- loo_cv_res %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + geom_abline() + 
    #geom_point(aes(x = observed, y = predicted), alpha = .2) +
    theme_bw() + xlab('Observed') + ylab('Predicted') + #ggtitle('Regression to 1st order, loo') + 
    xlim(minlim, maxlim) + ylim(minlim, maxlim) +
    annotate('text', x = .9*max(loo_cv_res$observed), y = 1.1*min(loo_cv_res$predicted), label = mylabel) 
  
  
  save(loo_cv_res, file = paste0("../../Results/model_fit_plots/first_order_oof_", name, ".RData")) 
  
  return(list(p1 = p1, p2 = p2, r2 = r2, r2_loo = r2_loo)) 
}


###################################################
### Regularized regression to second order ###
###################################################
get_second_order_fit <- function(df, all_species){
  #first pull out replicates 
  communities <- sapply(1:nrow(df), function(i) paste0(df[i,] %>% dplyr::select(-fitness), collapse = '-'))
  df <- df %>% mutate(community = communities)
  unique_communities <- unique(communities) 
  
  #set up regression
  y <- as.matrix(df$fitness)   
  f <- as.formula(y ~ .*.)
  x <- model.matrix(f, df %>% dplyr::select(all_of(all_species)))
  
  fold_ids <- get_folds(df, unique_communities, n_folds = 10)
  
  cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
  
  predicted_y <- predict(cv_fit, x, s = "lambda.1se") 
  
  #plot(cv_fit, sign.lambda = 1)
  
  res_cv <- data.frame(exp = df$community, observed = y, predicted = array(predicted_y)) %>% 
    group_by(exp) %>% mutate(mean_observed = mean(observed))
  res_cv_unique <- res_cv %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  r2 <- cor(res_cv_unique$mean_observed, res_cv_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  p3 <- res_cv %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + 
    #geom_point(aes(x = observed, y = predicted), alpha = .2) + 
    geom_abline() + 
    theme_bw() + xlab('Observed') + ylab('Predicted')  + #ggtitle('Regression to 2nd order') + 
    annotate('text', x = .9*max(res_cv$observed), y = 1.1*min(res_cv$predicted), label = mylabel)
  
  print(p3)
  
  save(res_cv, file = paste0("../../Results/model_fit_plots/second_order_infit_", name, ".RData")) 
  ### leave one out ### 
  loo_cv_res <- tibble()
  for (i in 1:length(unique_communities)){
    
    print(paste0(i, 'out of ', length(unique_communities)))
    
    exp <- unique_communities[i]
    exp_index <- which(df$community == exp )
    
    #pairwise regression
    y <- as.matrix(df[-exp_index,]$fitness)
    f <- as.formula(y ~ .*.)
    x <- model.matrix(f, df[-exp_index,] %>% dplyr::select(all_of(all_species)))
    
    fold_ids <- get_folds(df[-exp_index,], unique_communities[-i], n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_index[1],]$fitness)
    f_test <- as.formula(y_test ~ .*.)
    x_test <- model.matrix(f_test, df[exp_index[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se")
    
    for (j in 1:length(exp_index)){
      tmp <- list(exp = exp, observed = df[exp_index[j],]$fitness, predicted = array(y_pred_cv))
      loo_cv_res <- rbind(loo_cv_res, tmp)
    }
  }
  
  loo_cv_res <- data.frame(loo_cv_res) %>% group_by(exp) %>% mutate(mean_observed = mean(observed))
  loo_cv_res_unique <- loo_cv_res %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  #get r2 for mean value across replicates
  r2_loo <- cor(loo_cv_res_unique$mean_observed, loo_cv_res_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2_loo, digits = 3)))
  p4 <- loo_cv_res %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + 
    geom_point(aes(x = observed, y = predicted), alpha = .2) + geom_abline() + 
    theme_bw() + xlab('Observed') + ylab('Predicted') + #ggtitle('Regression to 2nd order, loo') +
    annotate('text', x = .9*max(loo_cv_res$observed), y = 1.1*min(loo_cv_res$predicted), label = mylabel)
  
  
  save(loo_cv_res, file = paste0("../../Results/model_fit_plots/second_order_oof_", name, ".RData")) 
  
  return(list(p3 = p3, p4 = p4, r2 = r2, r2_loo = r2_loo))
  
}

###################################################
### Regularized regression to third order ###
###################################################
get_third_order_fit <- function(df, all_species){
  #first pull out replicates 
  communities <- sapply(1:nrow(df), function(i) paste0(df[i,] %>% dplyr::select(-fitness), collapse = '-'))
  df <- df %>% mutate(community = communities)
  unique_communities <- unique(communities) 
  
  #set up regression
  y <- as.matrix(df$fitness)   
  f <- as.formula(y ~ .*.*.)
  x <- model.matrix(f, df %>% dplyr::select(all_of(all_species)))
  
  fold_ids <- get_folds(df, unique_communities, n_folds = 10)
  
  cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
  
  predicted_y <- predict(cv_fit, x, s = "lambda.1se") 
  
  #plot(cv_fit, sign.lambda = 1)
  
  res_cv <- data.frame(exp = df$community, observed = y, predicted = array(predicted_y)) %>% 
    group_by(exp) %>% mutate(mean_observed = mean(observed))
  res_cv_unique <- res_cv %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  r2 <- cor(res_cv_unique$mean_observed, res_cv_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  p3 <- res_cv %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + 
    geom_point(aes(x = observed, y = predicted), alpha = .2) + geom_abline() + 
    theme_bw() + xlab('Observed') + ylab('Predicted')  + #ggtitle('Regression to 3rd order') + 
    theme(panel.border = element_rect(fill=NA, colour = "black", size = 1)) +
    annotate('text', x = .9*max(res_cv$observed), y = 1.1*min(res_cv$predicted), label = mylabel)
  
  
  ### leave one out ### 
  loo_cv_res <- tibble()
  count <- 1
  for (i in 1:length(unique_communities)){
    
    print(paste0(i, 'out of ', length(unique_communities)))
    count <- count + 1
    
    exp <- unique_communities[i]
    exp_index <- which(df$community == exp )
    
    #third order regression
    y <- as.matrix(df[-exp_index,]$fitness)
    f <- as.formula(y ~ .*.*.)
    x <- model.matrix(f, df[-exp_index,] %>% dplyr::select(all_of(all_species)))
    
    fold_ids <- get_folds(df[-exp_index,], unique_communities[-i], n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_index[1],]$fitness)
    f_test <- as.formula(y_test ~ .*.*.)
    x_test <- model.matrix(f_test, df[exp_index[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se")
    
    for (j in 1:length(exp_index)){
      tmp <- list(exp = exp, observed = df[exp_index[j],]$fitness, predicted = array(y_pred_cv))
      loo_cv_res <- rbind(loo_cv_res, tmp)
    }
  }
  
  loo_cv_res <- data.frame(loo_cv_res) %>% group_by(exp) %>% mutate(mean_observed = mean(observed))
  loo_cv_res_unique <- loo_cv_res %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  #get r2 for mean value across replicates
  r2_loo <- cor(loo_cv_res_unique$mean_observed, loo_cv_res_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2_loo, digits = 3)))
  p4 <- loo_cv_res %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + 
    geom_point(aes(x = observed, y = predicted), alpha = .2) + geom_abline() + 
    theme_bw() + xlab('Observed') + ylab('Predicted') + #ggtitle('Regression to 3rd order, loo') +
    theme(panel.border = element_rect(fill=NA, colour = "black", size = 1)) + 
    annotate('text', x = .9*max(loo_cv_res$observed), y = 1.1*min(loo_cv_res$predicted), label = mylabel)
  
  #list(res_cv = res_cv, loo_cv_res = loo_cv_res)
  
  save(res_cv, file = paste0("../../Results/model_fit_plots/third_order_infit_", name, ".RData")) 
  save(loo_cv_res, file = paste0("../../Results/model_fit_plots/third_order_oof_", name, ".RData")) 
  
  return(list(p3 = p3, p4 = p4, r2 = r2, r2_loo = r2_loo))
  
}

data_path <- "../../Data/"
all_datasets <- list.files(data_path, pattern = '.csv') 

r2_res <- tibble()
for (dataset in all_datasets){
  
  df <- read_csv(paste0(data_path, dataset))
  all_species <- colnames(df %>% dplyr::select(-fitness))
  name <- sub(".csv", "", dataset)
  print(name)
  
  ### make fourier basis ###
  fourier_basis <- df %>% dplyr::select(all_of(all_species))
  fourier_basis[fourier_basis == 0] <- -1 
  fourier_df <- data.frame(fourier_basis, df %>% dplyr::select(-all_of(all_species))) 
  
  ### get first order fit ### 
  first_order_res <- get_first_order_fit(fourier_df)
  
  p_first_order <- plot_grid(first_order_res$p1, first_order_res$p2)
  
  ggsave(paste0("../../Results/model_fit_plots/first_order_plot_", name, ".png"), plot = p_first_order, width = 6.67, height = 5)
  
  ### get second order fit ###
  second_order_res <- get_second_order_fit(fourier_df, all_species)
  
  p_second_order <- plot_grid(second_order_res$p3, second_order_res$p4)
  show(p_second_order)
  
  ggsave(paste0("../../Results/model_fit_plots/second_order_plot_", name, ".png"), plot = p_second_order, width = 6.67, height = 5)
  
  #save(file = paste0("../../Results/model_fit_plots/second_order_loo_", name, ".RData"), second_order_res$p4)#, width = 6.67, height = 5)
  
  ### get third order fit ###
  third_order_res <- get_third_order_fit(fourier_df, all_species)
  
  p_third_order <- plot_grid(third_order_res$p3, third_order_res$p4)
  show(p_third_order)
  
  ggsave(paste0("../../Results/model_fit_plots/third_order_plot_", name, ".png"), plot = p_third_order, width = 6.67, height = 5)
  
  tmp <- list(name = name, r2_1o = first_order_res$r2, r2_loo_1o = first_order_res$r2_loo, 
              r2_2o = second_order_res$r2, r2_loo_2o = second_order_res$r2_loo, 
              r2_3o = third_order_res$r2, r2_loo_3o = third_order_res$r2_loo)
  r2_res <- rbind(r2_res, tmp)
  print(r2_res)
}

saveRDS(r2_res, file = '../../Results/model_fit_data/r2_res.RDS')


###############################
# just for first order reruns #
###############################

data_path <- "../../Data/"
all_datasets <- list.files(data_path, pattern = '.csv') 

new_r2_res <- tibble()
for (dataset in all_datasets){
  
  df <- read_csv(paste0(data_path, dataset))
  all_species <- colnames(df %>% dplyr::select(-fitness))
  name <- sub(".csv", "", dataset)
  print(name)
  
  ### make fourier basis ###
  fourier_basis <- df %>% dplyr::select(all_of(all_species))
  fourier_basis[fourier_basis == 0] <- -1 
  fourier_df <- data.frame(fourier_basis, df %>% dplyr::select(-all_of(all_species))) 
  
  ### get first order fit ### 
  first_order_res <- get_first_order_fit(fourier_df)
  
  p_first_order <- plot_grid(first_order_res$p1, first_order_res$p2)
  
  ggsave(paste0("../../Results/model_fit_plots/first_order_plot_", name, ".png"), plot = p_first_order, width = 6.67, height = 5)
  
  ### get second order fit ###
  #second_order_res <- get_second_order_fit(fourier_df, all_species)
  
  #p_second_order <- plot_grid(second_order_res$p3, second_order_res$p4)
  #show(p_second_order)
  
  #ggsave(paste0("../../Results/model_fit_plots/second_order_plot_", name, ".png"), plot = p_second_order, width = 6.67, height = 5)
  
  #save( second_order_res$p3,file = paste0("../../Results/model_fit_plots/second_order_infit_", name, ".RData"))#, width = 6.67, height = 5)
  #save(file = paste0("../../Results/model_fit_plots/second_order_loo_", name, ".RData"), second_order_res$p4)#, width = 6.67, height = 5)
  
  ### get third order fit ###
  third_order_res <- get_third_order_fit(fourier_df, all_species)
  
  #p_third_order <- plot_grid(third_order_res$p3, third_order_res$p4)
  #show(p_third_order)
  
  #ggsave(paste0("../Results/model_fit_plots/third_order_plot_", name, ".png"), plot = p_third_order, width = 6.67, height = 5)
  
  tmp <- list(name = name, r2_1o = first_order_res$r2, r2_loo_1o = first_order_res$r2_loo) #, 
  #r2_2o = second_order_res$r2, r2_loo_2o = second_order_res$r2_loo) #, 
  #r2_3o = third_order_res$r2, r2_loo_3o = third_order_res$r2_loo)
  new_r2_res <- rbind(new_r2_res, tmp)
  print(new_r2_res)
}

#####################################################
# just for jared data (don't run this all the time) #
#####################################################

df <- read_csv("~/linear_approach/Original_Data/output_chip1_sucrose.csv")
name <- 'jared_data'
all_species <- c('A10','A3','A7','B5','B9','C11','C4','C6','D11','E10','F10','F5','G4','H7')
N <- length(all_species)

#subset
df <- df[c(3:16,25)]

df <- df[-which(is.na(df$t3)),] #some fitness values not measured
colnames(df)[which(colnames(df) == 't3')] <- 'fitness'

#make presence/absence
df[,1:N][df[,1:N] > 0] <- 1

### make fourier basis ###
fourier_basis <- df %>% dplyr::select(all_of(all_species))
fourier_basis[fourier_basis == 0] <- -1 
fourier_df <- data.frame(fourier_basis, df %>% dplyr::select(-all_of(all_species))) 

### get first order fit ### 
first_order_res <- get_first_order_fit(fourier_df)

p_first_order <- plot_grid(first_order_res$p1, first_order_res$p2)
show(p_first_order)

#ggsave(paste0("../.../Results/model_fit_plots/first_order_plot_", name, ".png"), plot = p_first_order, width = 6.67, height = 5)

### get second order fit ###
second_order_res <- get_second_order_fit(fourier_df, all_species)

p_second_order <- plot_grid(second_order_res$p3, second_order_res$p4)
show(p_second_order)

#ggsave(paste0("../../Results/model_fit_plots/second_order_plot_", name, ".png"), plot = p_second_order, width = 6.67, height = 5)

### get third order fit ###
third_order_res <- get_third_order_fit(fourier_df, all_species)

p_third_order <- plot_grid(third_order_res$p3, third_order_res$p4)
show(p_third_order)

#ggsave(paste0("../Results/model_fit_plots/third_order_plot_", name, ".png"), plot = p_third_order, width = 6.67, height = 5)


###########
get_first_order_fit_unregularized <- function(df){
  
  #first pull out replicates 
  communities <- sapply(1:nrow(df), function(i) paste0(df[i,] %>% dplyr::select(-fitness), collapse = '-'))
  df <- df %>% mutate(community = communities)
  unique_communities <- unique(communities) 
  
  #first order linear regression, all data infit
  lin_fit <- lm(fitness ~ . , df %>% dplyr::select(-community))
  
  res_linear <- data.frame(exp = df$community, observed = df$fitness, predicted = lin_fit$fitted.values) %>% 
    group_by(exp) %>% mutate(mean_observed = mean(observed)) %>% ungroup()
  res_linear_unique <- res_linear %>% distinct(exp, .keep_all = TRUE)
  
  r2 <- cor(res_linear_unique$mean_observed, res_linear_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  p1 <- res_linear %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + geom_abline() + 
    geom_point(aes(x = observed, y = predicted), alpha = .2) +
    theme_bw() + xlab('Observed') + ylab('Predicted') + #ggtitle('Regression to 1st order') + 
    annotate('text', x = .9*max(res_linear$observed), y = 1.1*min(res_linear$predicted), label = mylabel) 
  
  
  ##### basic linear model with loo ######
  loo_res <- tibble()
  for (exp in unique_communities){
    exp_index <- which(df$community == exp)
    lin_fit <- lm(fitness ~ . , df[-exp_index,] %>% dplyr::select(-community))
    predicted_vals <- predict(lin_fit, df[exp_index[1],] %>% dplyr::select(-community)) 
    
    for (i in 1:length(exp_index)){
      tmp <- list(exp = exp, observed = df[exp_index[i],]$fitness, predicted = predicted_vals)
      loo_res <- rbind(loo_res, tmp)
    }
  }
  
  loo_res <- data.frame(loo_res) %>% group_by(exp) %>% mutate(mean_observed = mean(observed))
  loo_res_unique <- loo_res %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  #r2 on mean observed value across replicates
  r2_loo <- cor(loo_res_unique$mean_observed, loo_res_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2_loo, digits = 3)))
  p2 <- loo_res %>% ggplot(aes(x = mean_observed, y = predicted)) + geom_point() + geom_abline() + 
    geom_point(aes(x = observed, y = predicted), alpha = .2) +
    theme_bw() + xlab('Observed') + ylab('Predicted') + #ggtitle('Regression to 1st order, loo') + 
    annotate('text', x = .9*max(loo_res$observed), y = 1.1*min(loo_res$predicted), label = mylabel) 
  
  save(res_linear, file = paste0("../../Results/model_fit_plots/first_order_infit_", name, ".RData")) 
  save(loo_res, file = paste0("../../Results/model_fit_plots/first_order_oof_", name, ".RData")) 
  
  return(list(p1 = p1, p2 = p2, r2 = r2, r2_loo = r2_loo)) 
}

