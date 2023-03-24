library(tidyverse)
library(readr)
library(glmnet)
library(ggpubr)

rmse <- function(observed, predicted){
  return(sqrt(mean((observed - predicted)^2)))
}

data_path <- "../../Data/" #data_path <- "~/linear_approach/Data/" 
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
  
  #print(name)
  
  ###### basic linear model #########
  lin_fit <- lm(fitness ~ . , fourier_df)
  add_coefs <- coef(lin_fit)[2:ncol(fourier_df)]
  
  res_linear <- data.frame(fitness = df$fitness, predicted = lin_fit$fitted.values)
  
  r <- rmse(res_linear$fitness, res_linear$predicted)
  s <- mean(abs(add_coefs))
  
  rs <- r/s 
  print(rs)
  
  tmp <- list(name = name, rs = rs, rep = 'empirical') 
  rs_res <- rbind(rs_res, tmp)
  
  scramble_res <- tibble()
  for (rep in 1:100){
    scramble_df <- fourier_df
    scramble_df$fitness <- sample(fourier_df$fitness, length(fourier_df$fitness))
    
    scramble_fit <- lm(fitness ~ . , scramble_df)
    add_coefs <- coef(scramble_fit)[2:ncol(scramble_df)]
    
    res_linear <- data.frame(fitness = scramble_df$fitness, predicted = scramble_fit$fitted.values)
    
    #calculate
    r <-rmse(res_linear$fitness, res_linear$predicted)
    s <- mean(abs(add_coefs))
    
    rs <- r/s 
    
    tmp <- list(name = name, rs = rs, rep = rep)
    scramble_res <- rbind(scramble_res, tmp)
    rs_res <- rbind(rs_res, tmp)
  }
}


##### plot relative/empirical #
rs_plot_res <- left_join(rs_res %>% filter(rep != 'empirical'), rs_res %>% filter(rep == 'empirical'), by = 'name') %>% mutate(rel_rs = rs.y/rs.x)

rs_plot_res$name <- factor(rs_plot_res$name, levels = c('Butyrate', 'Langenheder', 'starch_data', 'pyoverdine', 'GE_biomass', 'jared_data'))

my_labs <- c('Butyrate production', 'Xylose oxidation', 'Starch hydrolysis', 'Total biomass', 'Siderophore production', 'Abundance of a single strain')

#gives points as in old fig
#rs_plot_res %>% 
#ggplot(aes(x = name, y = rel_rs)) +
#geom_point() +
#xlab('Dataset') + ylab('Normalized Landscape Ruggedness (r/s)') + 
#theme_bw() +
#ylim(0,1) + 
#theme(text = element_text(size = 16), panel.border = element_rect(fill=NA, colour = "black", size = 1.5)) 

library(stringr)
my_labs <- str_wrap(my_labs, width = 12)

rs_plot_res %>% 
  ggplot(aes(x = name, y = rel_rs)) +
  geom_boxplot() +
  xlab('Dataset') + ylab('Normalized Landscape Ruggedness (r/s)') + 
  theme_bw() +
  ylim(0,1) + 
  scale_x_discrete(labels = my_labs, guide = guide_axis(n.dodge = 1)) + #angle = 5)) +
  theme(text = element_text(size = 13.5), 
        panel.border = element_rect(fill=NA, colour = "black", size = 1.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

ggsave(filename = '../../Figures/rs_plot_norm.png', height = 4.55, width = 6.7, device = 'png', dpi = 700)


####
rs_res <- rs_res %>% mutate(norm_rs = rs/scramble_mean)

rs_res %>% filter(rep != 'empirical') %>% 
  ggplot(aes(x = name, y = norm_rs)) + geom_violin() +#fill = '#619CFF', alpha = .8) + 
  geom_point(aes(x = name, y = norm_rs), 
             color = '#F8766D', shape = 9, size = 3, 
             data = rs_res %>% filter(rep == 'empirical')) +
  xlab('Dataset') + ylab('Landscape Ruggedness (r/s)') + 
  theme_bw() +
  theme(text = element_text(size = 16), panel.border = element_rect(fill=NA, colour = "black", size = 2)) 

ggsave(filename = '../../Figures/rs_plot_norm_violin.png', height = 5, width = 8, device = 'png', dpi = 700)

#normalize by mean of scramble
scramble_mean <- rs_res %>% group_by(name) %>% filter(rep!= 'empirical') %>%
  mutate(scramble_mean = mean(rs)) %>% distinct(name, scramble_mean)

scramble_sd <- rs_res %>% group_by(name) %>% filter(rep!= 'empirical') %>%
  mutate(scramble_sd = sd(rs)) %>% distinct(name, scramble_sd)

rs_res <- left_join(rs_res, scramble_mean, by = 'name')
rs_res <- left_join(rs_res, scramble_sd, by = 'name')

rs_res <- rs_res %>% mutate(norm_rs = rs/scramble_mean)

rs_res %>% filter(rep != 'empirical') %>% 
  ggplot(aes(x = name, y = norm_rs)) + geom_violin() +#fill = '#619CFF', alpha = .8) + 
  geom_point(aes(x = name, y = norm_rs), 
             color = 'red', shape = 18, size = 4, 
             data = rs_res %>% filter(rep == 'empirical')) +
  xlab('Dataset') + ylab('Normalized Landscape Ruggedness (r/s)') + 
  theme_bw() +
  theme(text = element_text(size = 16), panel.border = element_rect(fill=NA, colour = "black", size = 1)) 

ggsave(filename = '../../Figures/rs_plot_norm_violin.png', height = 5, width = 8, device = 'png', dpi = 750)

rs_res %>% filter(rep != 'empirical') %>% 
  ggplot(aes(x = name, y = norm_rs)) + geom_boxplot() +#fill = '#619CFF', alpha = .8) + 
  geom_point(aes(x = name, y = norm_rs), 
             color = 'red', shape = 18, size = 4, 
             data = rs_res %>% filter(rep == 'empirical')) +
  xlab('Dataset') + ylab('Normalized Landscape Ruggedness (r/s)') + 
  theme_bw() +
  theme(text = element_text(size = 16), panel.border = element_rect(fill=NA, colour = "black", size = 1)) 


ggsave(filename = '../../Figures/rs_plot_norm_boxplot.png', height = 5, width = 8, device = 'png', dpi = 750)

######################
# unnormalized below #
######################

rs_res %>% filter(rep != 'empirical') %>% 
  ggplot(aes(x = name, y = rs)) + geom_boxplot(alpha = .8) + #fill = '#619CFF',
  geom_point(aes(x = name, y = rs), color = 'red', shape = 18, size = 4, data = rs_res %>% filter(rep == 'empirical')) +
  theme_bw() + xlab('Dataset') + ylab('r/s')

rs_res %>% filter(rep != 'empirical') %>% 
  ggplot(aes(x = name, y = scramble_mean)) + geom_point() + #fill = '#619CFF',
  geom_linerange(aes(ymin = scramble_mean - scramble_sd, ymax = scramble_mean + scramble_sd)) +
  geom_point(aes(x = name, y = rs), color = 'red', shape = 18, size = 3, data = rs_res %>% filter(rep == 'empirical')) +
  theme_bw() + xlab('Dataset') + ylab('Landscape Ruggedness (r/s)') +
  theme(text = element_text(size = 16), panel.border = element_rect(fill=NA, colour = "black", size = 1))

ggsave(filename = '../../Figures/rs_plot_all.png', height = 5, width = 8, device = 'png', dpi = 750)


rs_res %>% 
  ggplot(aes(x = as.factor(rep == 'empirical'), y = rs)) + geom_boxplot() + 
  facet_wrap(~name) + 
  geom_point(aes(x = name, y = rs), color = 'purple3', shape = 9, size = 2, data = rs_res %>% filter(rep == 'empirical')) +
  theme_bw() + xlab('Dataset') 

rs_res %>% filter(rep!= 'empirical') %>%
  ggplot(aes(y = rs)) + geom_boxplot() + 
  facet_wrap(~name, scales = 'free') + 
  geom_point(aes(x = 0, y = rs), color = 'purple3', shape = 9, size = 2, data = rs_res %>% filter(rep == 'empirical')) +
  theme_bw() + xlab('Dataset') 

