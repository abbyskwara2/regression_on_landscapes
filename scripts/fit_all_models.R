library(tidyverse)
library(readr)
library(glmnet)

set.seed(1)

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


#####################################
## function to fit regularized     ##
## regression to a specified order ##
#####################################

get_model_fit <- function(df, all_species, order){
  #first pull out replicates 
  communities <- sapply(1:nrow(df), function(i) paste0(df[i,] %>% dplyr::select(-fitness), collapse = '-'))
  df <- df %>% mutate(community = communities)
  unique_communities <- unique(communities) 
  
  #set up regression
  y <- as.matrix(df$fitness)   
  if (order == 1){
    f <- as.formula(y ~ .)
  } else if (order == 2){
    f <- as.formula(y ~ .*.)
  } else if (order == 3){
    f <- as.formula(y ~ .*.*.)
  }
  
  x <- model.matrix(f, df %>% dplyr::select(all_of(all_species)))
  
  #fit model
  fold_ids <- get_folds(df, unique_communities, n_folds = 10)
  cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
  
  #make predictions
  predicted_y <- predict(cv_fit, x, s = "lambda.1se") 
  
  #gather
  res_cv <- data.frame(exp = df$community, observed = y, predicted = array(predicted_y)) %>% 
    group_by(exp) %>% mutate(mean_observed = mean(observed))
  res_cv_unique <- res_cv %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  r2 <- cor(res_cv_unique$mean_observed, res_cv_unique$predicted)^2
  
  if (order == 1){
    save(res_cv, file = paste0("../../Results/model_fit_plots/first_order_infit_", name, ".RData")) 
  } else if (order == 2){
    save(res_cv, file = paste0("../../Results/model_fit_plots/second_order_infit_", name, ".RData")) 
  } else if (order == 3){
    save(res_cv, file = paste0("../../Results/model_fit_plots/third_order_infit_", name, ".RData")) 
  }
  
  ### leave one out ### 
  loo_cv_res <- tibble()
  for (i in 1:length(unique_communities)){
    
    print(paste0(i, 'out of ', length(unique_communities)))
    
    exp <- unique_communities[i]
    exp_index <- which(df$community == exp )
    
    #pairwise regression
    y <- as.matrix(df[-exp_index,]$fitness)
    
    if (order == 1){
      f <- as.formula(y ~ .)
    } else if (order == 2){
      f <- as.formula(y ~ .*.)
    } else if (order == 3){
      f <- as.formula(y ~ .*.*.)
    }
  
    x <- model.matrix(f, df[-exp_index,] %>% dplyr::select(all_of(all_species)))
    
    #fit model
    fold_ids <- get_folds(df[-exp_index,], unique_communities[-i], n_folds = 10)
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_index[1],]$fitness)
    
    if (order == 1){
      f_test <- as.formula(y_test ~ .)
    } else if (order == 2){
      f_test <- as.formula(y_test ~ .*.)
    } else if (order == 3){
      f_test <- as.formula(y_test ~ .*.*.)
    }
    
    #make predictions
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
  
  if (order == 1){
    save(loo_cv_res, file = paste0("../../Results/model_fit_plots/first_order_oof_", name, ".RData")) 
  } else if (order == 2){
    save(loo_cv_res, file = paste0("../../Results/model_fit_plots/second_order_oof_", name, ".RData")) 
  } else if (order == 3){
    save(loo_cv_res, file = paste0("../../Results/model_fit_plots/third_order_oof_", name, ".RData")) 
  }
  
  
  return(list(r2 = r2, r2_loo = r2_loo))
}


##### fit models to all datasets #### 
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
  first_order_res <- get_model_fit(fourier_df, all_species, order = 1)
  
  ### get second order fit ###
  second_order_res <- get_model_fit(fourier_df, all_species, order = 2)
  
  ### get third order fit ###
  third_order_res <- get_model_fit(fourier_df, all_species, order = 3)
  
 
  tmp <- list(name = name, r2_1o = first_order_res$r2, r2_loo_1o = first_order_res$r2_loo, 
              r2_2o = second_order_res$r2, r2_loo_2o = second_order_res$r2_loo, 
              r2_3o = third_order_res$r2, r2_loo_3o = third_order_res$r2_loo)
  r2_res <- rbind(r2_res, tmp)
  print(r2_res)
}

saveRDS(r2_res, file = '../../Results/model_fit_data/r2_res.RDS')
