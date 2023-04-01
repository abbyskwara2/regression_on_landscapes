library(tidyverse)
library(readr)

data_path <- "../../original_data/"
all_datasets <- list.files(data_path) 

#################
## Clark et al. #
#################
butyrate_df <- read_csv(paste0(data_path, "2020_02_28_MasterDF.csv"))
butyrate_df <- as.data.frame(butyrate_df)

#remove data with HB
butyrate_df <- butyrate_df[-c(which(butyrate_df$HB == 1)),]
butyrate_df <- butyrate_df %>% dplyr::select(-HB)

# remove data with contamination flag
butyrate_df <- butyrate_df[c(which(butyrate_df$`Contamination?` == 'No')),]

#make presence/absence df
species_cols <- 9:33
fitness <- butyrate_df$Butyrate 
df <- cbind(butyrate_df[,species_cols], fitness)

write.csv(x = df, file = "../../data/Butyrate.csv", row.names = FALSE)

###############################
## Sanchez-gorostiaga et al. ##
###############################

starch_df <- read_csv(paste0(data_path, "measured_starch_dataframe.csv"))

#split by community, get species
communities <- strsplit(starch_df$foc, split = '-')
all_species <- unique(unlist(communities))
all_species <- all_species[all_species != ""] 

#create presence absence data
species_presence <- data.frame(matrix(nrow = length(communities), ncol = length(all_species)))
for (species in 1:length(all_species)){
  species_presence[,species] <- sapply(1:length(communities), 
                                       function(i) ifelse(all_species[species] %in% communities[[i]],1,0))
}

colnames(species_presence) <- all_species

#add function col 
df <- species_presence %>% mutate(fitness = starch_df$F)

write.csv(x = df, file = "../../data/starch_data.csv", row.names = FALSE)


##################################################
## Diaz-Colunga et al. (Siderophore production) ##
##################################################

pyoverdine_df <- read.table(paste0(data_path, "ge_pyoverdine_final.txt"), header = TRUE)

pyoverdine_df <- pyoverdine_df %>% pivot_longer(-community) %>% filter(!is.na(value)) 

#background communities 
communities <- strsplit(pyoverdine_df$community, ',') 
fitness <- pyoverdine_df$value

all_species <- 1:8
sp_names <- paste0('SP_', all_species)
N <- length(all_species)

#make presence/absence data
species_presence <- data.frame(matrix(nrow = length(communities), ncol = N))
for (species in 1:N){
  species_presence[,species] <- sapply(1:length(communities), 
                                       function(i) ifelse(all_species[species] %in% communities[[i]],1,0))
}

colnames(species_presence) <- sp_names

df <- cbind(species_presence, fitness)

write.csv(x = df, file = "../../data/pyoverdine.csv", row.names = FALSE)

########################
## Langenheder et al. ##
########################
df <- read_csv(paste0(data_path, "Langenheder.etal.2010.xylose.48h.csv"))

#rename function col
colnames(df)[7] <- 'fitness'

write.csv(x = df, file = "../../data/oxidation.csv", row.names = FALSE)

#########################################
## Diaz-Colunga et al. (total biomass) ##
#########################################

GE_df <- read.table(paste0(data_path, 'GE_biomass.txt'), header = TRUE)

#background communities 
communities <- sapply(1:nrow(GE_df), 
                      function(i) as.numeric(strsplit(as.character(GE_df$background_community[i]), "")[[1]]))
fitness <- GE_df$background_fun
communities_knock_in <- lapply(1:nrow(GE_df), 
                               function(i) c(communities[[i]], as.numeric(GE_df$knock_in[i]))) 
fitness_knock_in <- GE_df$background_fun + GE_df$d_fun

communities <- c(communities,communities_knock_in)
fitness <- c(fitness, fitness_knock_in)

all_species <- 1:8
sp_names <- paste0('SP_', all_species)
N <- length(all_species)

#create presence absence data
species_presence <- data.frame(matrix(nrow = length(communities), ncol = N))
for (species in 1:N){
  species_presence[,species] <- sapply(1:length(communities),
                                       function(i) ifelse(all_species[species] %in% communities[[i]],1,0))
}

colnames(species_presence) <- sp_names

df <- species_presence %>% mutate(fitness = fitness) %>% distinct()

write.csv(x = df, file = "../../data/GE_biomass.csv", row.names = FALSE)

#################
## Kehe et al. ##
#################

kehe_df <- read_csv(paste0(data_path,"output_chip1_sucrose.csv"))
all_species <- c('A10','A3','A7','B5','B9','C11','C4','C6','D11','E10','F10','F5','G4','H7')
N <- length(all_species)

#subset
kehe_df <- kehe_df[c(3:16,25)]

kehe_df <- kehe_df[-which(is.na(kehe_df$t3)),] #some t3 values not measured
colnames(kehe_df)[which(colnames(kehe_df) == 't3')] <- 'fitness'

#make binary presence/absence
kehe_df[,1:N][kehe_df[,1:N] > 0] <- 1

write.csv(x = kehe_df, file = "../../data/kehe_data.csv", row.names = FALSE)
