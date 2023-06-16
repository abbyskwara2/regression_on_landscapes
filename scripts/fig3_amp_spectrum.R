library(tidyverse)
library(cowplot)
library(extrafont)

set.seed(1)
choose_font('Arial')

#######################################
### fig 3A, full amplitude spectrum  ##
### for langenheder et al. data      ##
#######################################

df <- read_csv("../data/Langenheder.csv", show_col_types = FALSE)
df <- rbind(df, rep(0, 7))
all_species <- colnames(df %>% select(-fitness))

### make fourier basis 
fourier_basis <- df %>% select(all_of(all_species))
fourier_basis[fourier_basis == 0] <- -1 
fourier_df <- data.frame(fourier_basis, df %>% select(-all_of(all_species))) 

### find fourier coefs 
y <- as.matrix(fourier_df$fitness)   
f <- as.formula(y ~ .^6)
x <- model.matrix(f, fourier_df %>% dplyr::select(all_of(all_species)))

fourier_coefs <- MASS::ginv(x) %*% fourier_df$fitness 

### get coef names 
names <- strsplit(colnames(x), ':') 
orders <- sapply(1:length(names), function(i) length(names[[i]]))
orders[1] <- 0 #intercept 

fourier_coef_data <- data.frame(fourier_coefs) %>% 
  mutate(names = colnames(x), order = orders) %>%  
  group_by(order) %>% 
  mutate(amplitude = sum(fourier_coefs^2)) %>% 
  distinct(order, .keep_all = TRUE) %>% 
  mutate(norm_amplitude = 
           amplitude/var(fourier_df$fitness)*(nrow(fourier_df))/(nrow(fourier_df) - 1)) %>%  #sample variance correction
  ungroup()


#name empirical data
empirical_coef_data <- fourier_coef_data

#check amplitudes sum to 1 
#empirical_coef_data %>% filter(order !=0) %>% summarize(sum(norm_amplitude))

# randomized landscape comparison
amplitude_res <- tibble() 
n_reps <- 100

for (rep in 1:n_reps){
  rand_df <- fourier_df
  ### get randomized landscape
  rand_df$fitness <- sample(rand_df$fitness, length(rand_df$fitness))
  
  ### get fourier coefs
  y <- as.matrix(rand_df$fitness)   
  f <- as.formula(y ~ .^6)
  x <- model.matrix(f, rand_df %>% dplyr::select(all_of(all_species)))
  
  fourier_coefs <- MASS::ginv(x) %*% rand_df$fitness
  
  ### massage
  names <- strsplit(colnames(x), ':')
  orders <- sapply(1:length(names), function(i) length(names[[i]]))
  orders[1] <- 0 #intercept 
  
  fourier_coef_data <- data.frame(fourier_coefs) %>% 
    mutate(names = colnames(x), order = orders) %>% 
    group_by(order) %>% 
    mutate(amplitude = sum(fourier_coefs^2), 
           norm_amplitude = amplitude/var(fourier_df$fitness)*(nrow(fourier_df))/(nrow(fourier_df) - 1)) %>% 
    ungroup()
  
  coef_distinct <- fourier_coef_data %>% 
                   filter(order!= 0) %>% 
                   distinct(order, amplitude, norm_amplitude) %>% 
                   mutate(rep = rep)
  
  amplitude_res <- rbind(amplitude_res, coef_distinct)
}

#######################################
# plot empirical and randomized data  #
#######################################

#normalized
amplitude_res_toplot <- amplitude_res %>% 
                        group_by(order) %>%
                        mutate(mean_amp = mean(amplitude), 
                        sd_amp = sd(amplitude), 
                        mean_na = mean(norm_amplitude), 
                        sd_na = sd(norm_amplitude)) %>% 
                        distinct(order, mean_amp, sd_amp, mean_na, sd_na) %>% 
                        ungroup()

empirical_coef_data %>% filter(order > 0) %>%
  ggplot(aes(x = as.factor(order), y = norm_amplitude)) + 
  geom_point(aes(color = 'red'), color = 'red', size = 1.6) + 
  geom_path(group = 1, color = 'red', size = 1.2) + 
  geom_point(aes(x = order, y = mean_na), alpha = .5, size = 1.6,
             data = amplitude_res_toplot) + 
  geom_path(aes(x = order, y = mean_na), alpha = .5, size = 1.2,
            data = amplitude_res_toplot)+
  geom_linerange(aes(x = order, ymin = mean_na - sd_na, ymax = mean_na + sd_na),
                 alpha = .5, color = 'black',
                 data = amplitude_res_toplot, inherit.aes = FALSE) + 
  ylim(0,1) + 
  theme_bw() + 
  xlab('Order') + 
  ylab('Fraction Variance Explained') +
  theme(text = element_text(size = 18, family = "Arial"),
        panel.border = element_rect(fill=NA, colour = "black", size = 2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave(filename = '../results/figures/amplitude_spectrum_langenheder.pdf', 
       device = 'pdf', dpi = 600, height = 5, width = 6.2)


#########################################
## Fig 3C, inferred 3rd order spectrum ##
#########################################

#read in fits from fit_all_models.R
my_path <- "../results/"
my_pattern <- 'third_order_infit_coef_data'
coef_datasets <- list.files(path = my_path, pattern = my_pattern)

for (dataset in coef_datasets){
  name <- gsub(my_pattern, dataset, replace = '')
  name <- gsub(".RData", name, replace = '')
  
  load(file = paste0(my_path, dataset)) 
  
  colnames(coef_data) <- c('value', 'order', 'amplitude', 'norm_amplitude')
  
  total_fit_amplitude <- coef_data %>% filter(order > 0) %>% 
                         distinct(order, amplitude) %>% 
                         summarize(total = sum(amplitude))
  total_fit_amplitude <- total_fit_amplitude$total 
  
  coef_data <- coef_data %>% mutate(norm_amplitude = amplitude/total_fit_amplitude)
  
  p <- coef_data %>% filter(order > 0) %>% 
    distinct(order, amplitude, norm_amplitude) %>% 
    ggplot(aes(x = as.factor(order), y = norm_amplitude)) +  
    geom_point(color = 'red') + 
    geom_path(group = 1, size = 1, color = 'red') + 
    theme_bw() +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0,1)) +
    xlab('Order') + 
    ylab('Percent Variance Explained') +
    theme(text = element_text(size = 12.5, family = 'Arial'),
          panel.border = element_rect(fill=NA, colour = "black", size = 1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  print(p)
  
  ggsave(filename = paste0("../results/figures/amp_spectrum_reg_", name, ".pdf"), 
         height = 68, width = 65, units = 'mm', 
         plot = p, device = 'pdf', dpi = 600, limitsize = F)
}
