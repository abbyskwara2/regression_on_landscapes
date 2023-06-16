library(tidyverse)
library(extrafont)

#####################################################
# Code to make residuals plots from fitted models   #
# (here, second order regressions are used,         #
# although this can easily be changed).             #
# Standard residuals are plotted, as well as        #
# residuals as a function of community richness     #
#####################################################

make_residuals_plot <- function(res_cv, my_title){
  #processing for fourier basis
  coms <- sapply(1:nrow(res_cv), function(i) gsub('--1', '-0', res_cv$exp[i]))
  coms <- sapply(1:length(coms), function(i) ifelse(substr(coms[i],1,2) == '-1', 
                                                    paste(0, substr(coms[i], 4, 100), sep = '-'), 
                                                    coms[i] ))
  coms <- sapply(1:length(coms), function(i) strsplit(coms[i], '-'))
  
  
  pres_abs <- matrix(0, nrow = nrow(res_cv), ncol = length(coms[[1]]))
  for (i in 1:length(coms)){
    pres_abs[i,] <- as.numeric(coms[[i]])
  }
  
  richness <- rowSums(pres_abs)
  
  #mean observed value for distinct communities
  resids <- res_cv$mean_observed - res_cv$predicted
  
  resids_df <- data.frame(resids = resids, richness = richness, 
                          mean_obs = res_cv$mean_observed, exp = res_cv$exp) %>% 
    distinct(exp, .keep_all = TRUE)
  
  ### residuals vs community richness ###
  
  scale <- max(abs(resids_df$resids))
  p <- resids_df %>% 
    ggplot(aes(x = richness, y = resids)) + 
    geom_point(size = 1) + 
    geom_hline(yintercept = 0) +
    theme_classic() +  
    ggtitle(my_title) +
    scale_y_continuous(limits = c(-scale, scale)) +
    xlab('Community richness')  + ylab(expression(y - hat(y))) + 
    theme(text = element_text(family = "Arial"),
      plot.title = element_text(hjust = .5, vjust = 1),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()) 
  
  ### residuals ###
  p2 <- resids_df %>%
    ggplot(aes(x = mean_obs, y = resids)) + 
    geom_point(size = 1) + 
    theme_classic() +  
    xlab('Observed value (y)')  + ylab(expression(y - hat(y))) +
    geom_hline(yintercept = 0) + 
    scale_y_continuous(limits = c(-scale, scale)) +
    theme(text = element_text(family = "Arial"),
      plot.title = element_text(hjust = .5, vjust = 1),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()) +
    ggtitle(my_title) 
  
  return(list(p = p, p2 = p2))
}

#where are model fits saved? 
data_path <- '../results/' 

#where should these plots be saved?
save_path <- '../results/figures/residual_plots/'

#data
dataset_names <- c('Langenheder', 'Butyrate', 'starch_data', 'pyoverdine', 'GE_biomass' , 'kehe_data') 

my_titles <- list('Langenheder' = 'Growth on xylose', 
                  'Butyrate' = 'Butyrate production',
                  'starch_data' = 'Amylase secretion', 
                  'pyoverdine' = 'Siderophore\n production', 
                  'GE_biomass' = 'Total biomass', 
                  'kehe_data' = 'Abundance of a\nfocal strain')

p_list <- list()
p_res_list <- list()
for (name in dataset_names){
  load(file = paste0(data_path, "second_order_infit_", name, ".RData"))
  
  out <- make_residuals_plot(res_cv, my_titles[[name]]) 
  
  #community richness vs residuals
  p_list[[name]] <- out$p2
  
  #residuals
  p_res_list[[name]] <- out$p
  
}

### community richness vs residuals plot ###
plot_grid(p_list[[1]], NULL,
          p_list[[2]],NULL,
          p_list[[3]],NULL,
          NULL,NULL,NULL,
          NULL,NULL,NULL,
          p_list[[4]], NULL,
          p_list[[5]], NULL,
          p_list[[6]], NULL,
          nrow = 3,
          #arrange plots
          rel_widths = c(1, .05, 
                         1, .05, 
                         1, .05, 
                         1, .05,
                         1, .05,
                         1, .05), 
          rel_heights = c(1, .05, 1),
          #labels
          labels = c('A', '', 
                     'B', '',
                     'C', '', 
                     '','','',
                     '','','',
                     'D', '', 
                     'E', '', 
                     'F', ''))


ggsave(filename = paste0(save_path, 'residuals_vs_richness.pdf'),
       device = 'pdf', dpi = 750, height = 5, width = 7.5)

### actual residuals plot ###
plot_grid(p_res_list[[1]], NULL,
          p_res_list[[2]], NULL, 
          p_res_list[[3]], NULL, 
          NULL, NULL, NULL, 
          NULL, NULL, NULL, 
          p_res_list[[4]], NULL,
          p_res_list[[5]], NULL,
          p_res_list[[6]], NULL, 
          nrow = 3, 
          #arrange plots
          rel_widths = c(1, .05, 
                        1, .05, 
                        1, .05, 
                        1, .05, 
                        1, .05, 
                        1, .05), 
          rel_heights = c(1, .05, 1),
          #labels
          labels = c('A', '',
                     'B', '',
                     'C', '', 
                     '','','',
                     '','','',
                     'D', '',
                     'E', '', 
                     'F', ''))


ggsave(filename = paste0(save_path, 'second_order_residuals.pdf'),
       device = 'pdf', dpi = 600, height = 5, width = 7.5)


