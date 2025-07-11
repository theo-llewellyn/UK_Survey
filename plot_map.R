library(ggplot2)
library(maps)
library(tidyverse)
library(ggrepel)

setwd(dir = "")

#Data
worldmap = map_data('world')
sites <- read_csv('DATA/site_metadata.csv')
#Option 2
(map <- ggplot() + 
    geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),
                        fill = 'gray90', 
                        color = 'grey70') +
  geom_sf() + coord_sf(xlim = c(-10,3), 
                       ylim = c(49, 62)) +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, colour = Forest.comp)) + 
  geom_text_repel(data = sites, 
                  aes(x = Longitude, y = Latitude, label = Location),
                  size = 3) +
  xlab('longitude') +
  ylab('latitude') +
  scale_colour_discrete(name="Forest\nComposition",
                      breaks=c("Monosp", "MixBroad", "MixCon","Mixed.CB"),
                      labels=c("Monospecific", "Mixed Broadleaf", "Mixed Conifer","Conifer + Broadleaf")) +
  theme_bw())

png("FIGURES/map_UKSurvey.png",  res = 400, width = 2480, height = 2480)
map
dev.off()
