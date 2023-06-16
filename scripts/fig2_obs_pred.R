library(tidyverse)
library(readr)
library(extrafont)
library(cowplot)

make_infit_plot <- function(res_cv, xlabelval, ylabelval, my_fontsize = 17, replicates = FALSE){
  #filter experimental replicates
  res_cv_unique <- res_cv %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  #get r2 for mean value across replicates
  r2 <- cor(res_cv_unique$mean_observed, res_cv_unique$predicted)^2
  
  #generate plot
  if (replicates == FALSE){
    maxlim <- max(c(res_cv$mean_observed, res_cv$predicted))
    minlim <- min(c(res_cv$mean_observed, res_cv$predicted))
    mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
    
    p <- res_cv %>% ggplot(aes(x = mean_observed, y = predicted)) + 
      geom_point() + 
      geom_abline() + 
      theme_bw() + xlab('Observed') + ylab('Predicted') +
      scale_x_continuous(limits = c(minlim, maxlim)) + 
      scale_y_continuous(limits = c(minlim, maxlim)) + 
      theme(text = element_text(size = my_fontsize, family = "Arial"), 
            panel.border = element_rect(fill=NA, colour = "black", size = 1), 
            plot.background = element_rect(fill = "transparent"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      annotate('text', x = xlabelval, y = ylabelval, label = mylabel) 
    
  } else {
    maxlim <- max(c(res_cv$observed, res_cv$predicted))
    minlim <- min(c(res_cv$observed, res_cv$predicted))
    mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
    
    p <- res_cv %>% ggplot(aes(x = mean_observed, y = predicted)) + 
      geom_point() + 
      geom_point(aes(x = observed, y = predicted), alpha = .25) + 
      geom_abline() + 
      theme_bw() + xlab('Observed') + ylab('Predicted') + 
      scale_x_continuous(limits = c(minlim, maxlim)) + 
      scale_y_continuous(limits = c(minlim, maxlim)) + 
      theme(text = element_text(size = my_fontsize, family = "Arial"), 
            panel.border = element_rect(fill=NA, colour = "black", size = 1), 
            plot.background = element_rect(fill = "transparent"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      annotate('text', x = xlabelval, y = ylabelval, label = mylabel) 
  }
  
  return(p)
}

make_oof_plot <- function(loo_cv_res, xlabelval, ylabelval){
  loo_cv_res_unique <- loo_cv_res %>% ungroup() %>% distinct(exp, .keep_all = TRUE) 
  
  maxlim <- max(c(loo_cv_res$mean_observed, loo_cv_res$predicted))#*1.025
  minlim <- min(c(loo_cv_res$mean_observed, loo_cv_res$predicted))
  
  posneg <- ifelse(minlim > 0, 1, -1)
  
  #get r2 for mean value across replicates
  r2_loo <- cor(loo_cv_res_unique$mean_observed, loo_cv_res_unique$predicted)^2
  
  mylabel = bquote(italic(R)^2 == .(format(r2_loo, digits = 3)))
  p <- loo_cv_res %>% ggplot(aes(x = mean_observed, y = predicted)) + 
    geom_point() + 
    geom_abline() + 
    theme_bw() + xlab('Observed') + 
    ylab('Predicted') +  
    theme(text = element_text(size = 15.75, family = "Arial")) +
    scale_x_continuous(limits = c(minlim, maxlim), expand = c(0.05, 0)) + 
    scale_y_continuous(limits = c(minlim, maxlim)) +  
    theme(panel.border = element_rect(fill= NA, colour = "black"),
          panel.background = element_rect(fill = 'transparent', color = NA),
          plot.background = element_rect(fill = 'transparent', color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_cartesian(clip = 'off') + 
    annotate('text', xlabelval, y = ylabelval, label = mylabel)
  
  return(p)
}

#where are model fits?
data_path <- '../results/'

#where should plots be saved? 
save_path <- '../results/figures/'

## generate plots for all datasets ## 

# Langenheder data
load(file = paste0(data_path, "second_order_infit_Langenheder.RData"))
load(file = paste0(data_path, "second_order_oof_Langenheder.RData"))

p2_Langenheder <- make_infit_plot(res_cv, 1.26, .47)
p2_oof_Langenheder <- make_oof_plot(loo_cv_res, 1.26, .45)

# Butyrate data

load(file = paste0(data_path, "second_order_infit_Butyrate.RData"))
load(file = paste0(data_path, "second_order_oof_Butyrate.RData"))

p2_butyrate <- make_infit_plot(res_cv, 51, 4)
p2_butyrate_replicates <- make_infit_plot(res_cv, 59, 3, replicates = TRUE)
p2_oof_butyrate <- make_oof_plot(loo_cv_res,  51, 2.5)

ggsave(filename = paste0(save_path, 'so_fit_butyrate.pdf'), plot = p2_butyrate, 
       device = 'pdf', dpi = 700, height = 4, width = 4 )
ggsave(filename = paste0(save_path,'so_fit_replicates_butyrate.pdf'), plot = p2_butyrate_replicates, 
       device = 'pdf', dpi = 700, height = 4, width = 4 )
ggsave(filename = paste0(save_path,'so_fit_loo_butyrate.pdf'), plot = p2_oof_butyrate,
       device = 'pdf', dpi = 750, height = 4, width = 4)

# starch 
load(file = paste0(data_path, "second_order_infit_starch_data.RData"))
load(file = paste0(data_path, "second_order_oof_starch_data.RData"))

p2_starch_data <- make_infit_plot(res_cv, 31, 2.75)
p2_oof_starch_data <- make_oof_plot(loo_cv_res, 33.5,2)

# Pyoverdine secretion data 
load(file = paste0(data_path, "second_order_infit_pyoverdine.RData"))
load(file = paste0(data_path, "second_order_oof_pyoverdine.RData"))

p2_pyoverdine <- make_infit_plot(res_cv, .41, .04)
p2_pyoverdine_replicates <- make_infit_plot(res_cv, .42, .042, replicates = TRUE)
p2_oof_pyoverdine <- make_oof_plot(loo_cv_res, .42, .04)

# Total biomass 
load(file = paste0(data_path, "second_order_infit_GE_biomass.RData"))
load(file = paste0(data_path, "second_order_oof_GE_biomass.RData"))

p2_GE_biomass <- make_infit_plot(res_cv, .084, .022)
p2_GE_biomass_replicates <- make_infit_plot(res_cv, .084, .022, replicates = TRUE)
p2_oof_GE_biomass <- make_oof_plot(loo_cv_res, .084, .022)

# kehe data
load(file = paste0(data_path, "second_order_infit_kehe_data.RData"))
load(file = paste0(data_path, "second_order_oof_kehe_data.RData"))

p2_kehe_data <- make_infit_plot(res_cv, 7900, 1900) 
p2_kehe_data_replicates <- make_infit_plot(res_cv, 7900, 1500, replicates = TRUE)
p2_oof_kehe_data <- make_oof_plot(loo_cv_res, 7900, 1900)


#############################
##### alignment all infit ### 
##############################

palign <- align_plots(p2_butyrate, 
                      p2_Langenheder,
                      p2_starch_data, 
                      p2_pyoverdine,
                      p2_GE_biomass, 
                      p2_kehe_data)

ggsave(filename = paste0(save_path,'so_fit_Butyrate.pdf'), plot = palign[[1]], 
       device = 'pdf', dpi = 600, height = 4, width = 4) 

ggsave(filename = paste0(save_path, "so_fit_Langenheder.pdf"), plot = palign[[2]], 
       device = 'pdf', dpi = 600, height = 4, width = 4) 

ggsave(filename = paste0(save_path,'so_fit_starch_data.pdf'), plot = palign[[3]], 
       device = 'pdf', dpi = 600, height = 4, width = 4) 

ggsave(filename = paste0(save_path,'so_fit_pyoverdine.pdf'), plot = palign[[4]], 
       device = 'pdf', dpi = 600, height = 4, width = 4) 

ggsave(filename = paste0(save_path,'so_fit_GE_biomass.pdf'), plot = palign[[5]], 
       device = 'pdf', dpi = 600, height = 4, width = 4) 

ggsave(filename = paste0(save_path,'so_fit_kehe_data.pdf'), plot = palign[[6]], 
       device = 'pdf', dpi = 600, height = 4, width = 4) 


###############################
## alignment leave-one-out  ###
###############################

palign <- align_plots(p2_oof_butyrate, 
                      p2_oof_Langenheder, 
                      p2_oof_starch_data, 
                      p2_oof_pyoverdine, 
                      p2_oof_GE_biomass, 
                      p2_oof_kehe_data)

plot_grid(p2_oof_butyrate, 
          p2_oof_Langenheder, 
          p2_oof_starch_data, 
          p2_oof_pyoverdine, 
          p2_oof_GE_biomass, 
          p2_oof_kehe_data)

h <- 3
w <- 3

ggsave(filename = paste0(save_path,'so_fit_loo_Butyrate.pdf'), plot = palign[[1]], 
       device = 'pdf', dpi = 600, bg = 'transparent', height = h, width = w) 

ggsave(filename = paste0(save_path,'so_fit_loo_Langenheder.pdf'), plot = palign[[2]], 
       device = 'pdf', dpi = 600, bg = 'transparent',  height = h, width = w) 

ggsave(filename = paste0(save_path,'so_fit_loo_starch_data.pdf'), plot = palign[[3]], 
       device = 'pdf', dpi = 600, bg = 'transparent',  height = h, width = w) 

ggsave(filename = paste0(save_path,'so_fit_loo_pyoverdine.pdf'), plot = palign[[4]], 
       device = 'pdf', dpi = 600, bg = 'transparent', height = h, width = w) 

ggsave(filename = paste0(save_path,'so_fit_loo_GE_biomass.pdf'), plot = palign[[5]], 
       device = 'pdf', dpi = 600, bg = 'transparent',  height = h, width = w) 

ggsave(filename = paste0(save_path,'so_fit_loo_kehe_data.pdf'), plot = palign[[6]], 
       device = 'pdf', dpi = 600, bg = 'transparent', height = h, width = w) 



################################
## alignment with replicates ### 
###############################

palign <- align_plots(p2_butyrate_replicates, 
                      p2_pyoverdine_replicates, 
                      p2_GE_biomass_replicates, 
                      p2_kehe_data_replicates)

ggsave(filename = paste0(save_path,'so_fit_replicates_Butyrate.pdf'), plot = palign[[1]], 
       device = 'pdf', dpi = 600, height = 4, width = 4, bg= 'transparent') 

ggsave(filename = paste0(save_path,'so_fit_replicates_pyoverdine.pdf'), plot = palign[[2]], 
       device = 'pdf', dpi = 600, height = 4, width = 4, bg= 'transparent') 

ggsave(filename = paste0(save_path,'so_fit_replicates_GE_biomass.pdf'), plot = palign[[3]], 
       device = 'pdf', dpi = 600, height = 4, width = 4, bg= 'transparent') 

ggsave(filename = paste0(save_path,'so_fit_replicates_kehe_data.pdf'), plot = palign[[4]], 
       device = 'pdf', dpi = 600, height = 4, width = 4, bg= 'transparent') 

