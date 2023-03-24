library(tidyverse)
library(readr)
library(glmnet)
library(stringr)

rmse <- function(observed, predicted){
  return(sqrt(mean((observed - predicted)^2)))
}

data_path <- "../../Data/" 
all_datasets <- list.files(data_path, pattern = '.csv')

###############################
# calculate for all datasets ##
###############################
rs_res <- tibble()
for (dataset in all_datasets){
  
  df <- read_csv(paste0(data_path, dataset))
  name <- sub(".csv", "", dataset)
  
  ### make fourier basis ###
  fourier_basis <- df %>% dplyr::select(-fitness)
  fourier_basis[fourier_basis == 0] <- -1 
  fourier_df <- data.frame(fourier_basis, df %>% dplyr::select(fitness)) 
  
  #first order fit 
  lin_fit <- lm(fitness ~ . , fourier_df)
  add_coefs <- coef(lin_fit)[2:ncol(fourier_df)]
  
  res_linear <- data.frame(fitness = df$fitness, 
                           predicted = lin_fit$fitted.values)
  
  r <- rmse(res_linear$fitness, res_linear$predicted)
  s <- mean(abs(add_coefs))
  
  rs <- r/s 
  
  #add empirical value
  tmp <- list(name = name, rs = rs, rep = 'empirical') 
  rs_res <- rbind(rs_res, tmp)
  
  #randomized landscape comparisons
  scramble_res <- tibble()
  for (rep in 1:100){
    scramble_df <- fourier_df
    scramble_df$fitness <- sample(fourier_df$fitness, length(fourier_df$fitness))
    
    scramble_fit <- lm(fitness ~ . , scramble_df)
    add_coefs <- coef(scramble_fit)[2:ncol(scramble_df)]
    
    res_linear <- data.frame(fitness = scramble_df$fitness, 
                             predicted = scramble_fit$fitted.values)
    
    #calculate r,s
    r <-rmse(res_linear$fitness, res_linear$predicted)
    s <- mean(abs(add_coefs))
    
    rs <- r/s 
    
    tmp <- list(name = name, rs = rs, rep = rep)
    scramble_res <- rbind(scramble_res, tmp)
    rs_res <- rbind(rs_res, tmp)
  }
}


##### plot relative/empirical 
rs_plot_res <- left_join(rs_res %>% filter(rep != 'empirical'), 
                         rs_res %>% filter(rep == 'empirical'), by = 'name') %>% 
               mutate(rel_rs = rs.y/rs.x)

rs_plot_res$name <- factor(rs_plot_res$name, levels = c('Butyrate', 
                                                        'Langenheder', 
                                                        'starch_data', 
                                                        'pyoverdine', 
                                                        'GE_biomass', 
                                                        'jared_data'))

my_labs <- c('Butyrate production', 
             'Xylose oxidation', 
             'Starch hydrolysis', 
             'Siderophore production', 
             'Total biomass', 
             'Abundance of a single strain')

#name aesthetics
my_labs <- str_wrap(my_labs, width = 12)

rs_plot_res %>% 
  ggplot(aes(x = name, y = rel_rs)) +
  geom_boxplot() +
  xlab('Dataset') + 
  ylab('Normalized Landscape Ruggedness (r/s)') + 
  theme_bw() +
  ylim(0,1) + 
  scale_x_discrete(labels = my_labs, guide = guide_axis(n.dodge = 1)) + 
  theme(text = element_text(size = 13.5), 
        panel.border = element_rect(fill=NA, colour = "black", size = 1.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

ggsave(filename = '../../Figures/rs_plot_norm.pdf', height = 4.55, width = 6.7, device = 'pdf', dpi = 600)
