#Alpha diversity multivariate analysis
library(tidyverse)
library(vegan)
library(cowplot)

setwd('')

#join with site metadata
sites <- read_csv('DATA/site_metadata.csv')
metadata <- read_csv('DATA/metadata.csv') %>%
  inner_join(sites, by = c('locality' = 'Location'))

otu_data <- read_tsv('RESULTS/OTUs/OTUsotu_table.txt') %>%
  slice(-1)

otu_data %>%
  pivot_longer(ACMB2019_10_A01:ACMB2019_9_H11) %>% 
  pivot_wider(names_from = OTU_ID, values_from = value) -> OTU_long

OTU_long %>%
  column_to_rownames(var = "name") %>%
  mutate(total = rowSums(.,)) %>%
  ggplot(aes(x = total)) + 
  geom_histogram(binwidth = 50) +
  xlim(0,500)
  
#see where bulk of data sits
#log
OTU_long %>%
  column_to_rownames(var = "name") %>%
  mutate(total = rowSums(.,)) %>%
  ggplot(aes(x=1,y=total)) +
  geom_jitter() +
  scale_y_log10()

#rarefy to 10 seqs
rarefied_richness <- OTU_long %>%
  column_to_rownames(var = "name") %>%
  rarefy(., 10) %>%
  as_tibble(rownames ='name') %>%
  select(name, rarefied_richness = value)

#non-rarefied richness
richness <- OTU_long %>%
  column_to_rownames(var = "name") %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  mutate(total = rowSums(.,)) %>%
  select(total) %>%
  as_tibble(rownames ='name') %>%
  select(name, richness = total) %>%
  inner_join(rarefied_richness)


#join sample and alpha data
alpha_stats_metadata <- inner_join(richness, metadata, by = c('name' = 'project_sample_id'))
write_csv(alpha_stats_metadata[,1:3],'OTUs/alpha_div_stats.csv')

#plot country
a <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(country, richness, median, .desc =TRUE),, x = richness)) +
  geom_boxplot() +
  ylab('Country') +
  theme(axis.title.x=element_blank())
b <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(country, rarefied_richness, median, .desc =TRUE),, x = rarefied_richness)) +
  geom_boxplot() + 
  theme(axis.title=element_blank())
#host genus
c <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(Insect_genus, richness, median, .desc =TRUE), x = richness))+ 
  geom_boxplot() +
  ylab('Beetle genus') +
  theme(axis.title.x=element_blank())
d <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(Insect_genus, rarefied_richness, median, .desc =TRUE), x = rarefied_richness)) + 
  geom_boxplot() +
  theme(axis.title=element_blank())
#site
e <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(locality, richness, median, .desc =TRUE), 
                                 x = richness)) + 
  geom_boxplot() +
  ylab('Site') +
  theme(axis.title.x=element_blank())
f <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(locality, rarefied_richness, median, .desc =TRUE), x = rarefied_richness)) + 
  geom_boxplot() +
  theme(axis.title=element_blank())
#forest composition
g <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(Forest.comp, richness, median, .desc =TRUE), 
                                 x = richness)) + 
  geom_boxplot() +
  ylab('Forest Composition') +
  xlab('Fungal OTU richness') +
  scale_y_discrete(breaks=c("MixBroad","MixCon","Monosp", "Mixed.CB"),
                   labels=c("Mixed Broadleaf","Mixed Conifer","Monospecific","Conifer + Broadleaf"))
h <- ggplot(alpha_stats_metadata, aes(y = fct_reorder(Forest.comp, rarefied_richness, median, .desc =TRUE), x = rarefied_richness)) + 
  geom_boxplot() +
  theme(axis.title.y=element_blank()) +
  xlab('Fungal OTU rarefied richness') +
  scale_y_discrete(breaks=c("MixBroad","Mixed.CB","MixCon","Monosp"),
                   labels=c("Mixed Broadleaf","Conifer + Broadleaf","Mixed Conifer","Monospecific"))

png("FIGURES/Richness_Summary_UKSurvey.png",  res = 400, width = 2700, height = 3800)
plot_grid(a,b,c,d,e,f,g,h, ncol = 2, align = 'hv', axis = 'lb', labels = 'auto')
dev.off()
