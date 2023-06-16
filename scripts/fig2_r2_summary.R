library(tidyverse)
library(readr)
library(extrafont)
library(RColorBrewer)

#read in r2 summary from all datasets (saved from fit_all_models.R)
r2_res <- readRDS("../results/r2_res.RDS")

#collect names 
all_datasets <- r2_res$name

#long data format
infit_res <- r2_res %>% select(c(name, r2_1o, r2_2o, r2_3o)) 
colnames(infit_res) <- c('name', '1', '2','3')
infit_res <- infit_res %>% 
  pivot_longer(-c(name), names_to = 'order') %>% 
  mutate(infit = 'full dataset')


loo_res <- r2_res %>% select(c(name, r2_loo_1o, r2_loo_2o, r2_loo_3o))
colnames(loo_res) <- c('name', '1', '2', '3')
loo_res <- loo_res %>% 
  pivot_longer(-c(name), names_to = 'order') %>% 
  mutate(infit = 'LOO-CV') 


r2_res_all <- rbind(infit_res, loo_res)
r2_res_all$value <- as.numeric(r2_res_all$value)
plot_list <- list()

# make r2 summary plot, all datasets
color_1 <- brewer.pal(8, 'Greys')[6]
for (ind in 1:length(all_datasets)){
  
  dataset <- all_datasets[ind]
  p <- r2_res_all %>% filter(name == dataset) %>% 
    ggplot(aes(x = order, y = value, fill = infit, color = infit)) + 
    geom_col(width = .5, position = 'dodge') + theme_bw() +
    xlab('Model Order') +
    ylab(bquote(italic(R)^2)) + 
    coord_cartesian(ylim = c(0, 1)) + 
    scale_fill_manual(values = alpha(c(color_1, color_1), c(1, .01))) + 
    scale_color_manual(values = alpha(c(color_1, color_1), c(1, 1))) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0)) +
    guides(fill = guide_legend(''), color = guide_legend('')) + 
    theme_classic() + 
    theme(text = element_text(size = 18, family = 'Arial'),
          legend.position="bottom", 
          axis.title.x = element_text(size = 18, color = "black"), 
          axis.title.y = element_text(size = 18, face="bold", color = "black"),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = 'transparent', color = NA),
          legend.background = element_rect(fill = 'transparent'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  print(p)
  
  scale_val <- 1.75
  ggsave(filename = paste0('../results/figures/r2_summary_', dataset, '.pdf'), 
         plot = p, device = 'pdf', height = 37*scale_val, width = 50*scale_val, units = 'mm', 
         limitsize = FALSE,
         dpi = 500)
}

