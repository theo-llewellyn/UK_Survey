library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(stringi)

setwd('')

#read in TBAS and RDP taxonomy
taxonomy <- read_csv('assignments_reportUA4EAZEB.csv')
#read in dnaconsensus results
#BLAST using dnabarcoder cut-offs
dnabarcoder_tab <- read_tsv('OTUs/OTUstaxonomy.txt')
dnabarcoder_tab$OTU_ID <- gsub(':','',dnabarcoder_tab$OTU_ID)
#combine
taxonomy_UKSurvey <- left_join(taxonomy,dnabarcoder_tab, by = c('Query sequence' = 'OTU_ID'))

#Remove underscores from species
taxonomy_UKSurvey$species <- gsub("_", " ", taxonomy_UKSurvey$species)
taxonomy_UKSurvey$`Taxon assignment` <- gsub("_", " ", taxonomy_UKSurvey$`Taxon assignment`)
taxonomy_UKSurvey$rdp_species <- gsub("_", " ", taxonomy_UKSurvey$rdp_species)
#remove CBS codes from TBAS
taxonomy_UKSurvey$`Taxon assignment` <- gsub("^(([^ ]+ ){1}[^ ]+).*", "\\1", taxonomy_UKSurvey$`Taxon assignment`)
#remove species hypothesis codes from rdp
taxonomy_UKSurvey$rdp_species <- gsub("\\|.*", "", taxonomy_UKSurvey$rdp_species)


#Add column for whether dnabarcoder and T-BAS agree
taxonomy_UKSurvey[,c('phylum.match','class.match','order.match','family.match','genus.match','species.match')] <- NA


#dnabarcoder vs TBAS
for (i in 1:nrow(taxonomy_UKSurvey)) {
    #if blast species matches TBAS species then add Y(es) to match column else add N(o)
    if(grepl(taxonomy_UKSurvey$species[i], taxonomy_UKSurvey$`Taxon assignment`[i], fixed = TRUE) == TRUE){
      taxonomy_UKSurvey$species.match[i] <- 'Y'
    }else{
      taxonomy_UKSurvey$species.match[i] <- 'N'
    }
    if(grepl(taxonomy_UKSurvey$genus[i], taxonomy_UKSurvey$`Most common Genus-level assignment`[i], fixed = TRUE) == TRUE){
      taxonomy_UKSurvey$genus.match[i] <- 'Y'
    }else{
      taxonomy_UKSurvey$genus.match[i] <- 'N'
    }
    if(grepl(taxonomy_UKSurvey$family[i], taxonomy_UKSurvey$`Most common family-level assignment`[i], fixed = TRUE) == TRUE){
      taxonomy_UKSurvey$family.match[i] <- 'Y'
    }else{
      taxonomy_UKSurvey$family.match[i] <- 'N'
    }
    if(grepl(taxonomy_UKSurvey$order[i], taxonomy_UKSurvey$`Most common order-level assignment`[i], fixed = TRUE) == TRUE){
      taxonomy_UKSurvey$order.match[i] <- 'Y'
    }else{
      taxonomy_UKSurvey$order.match[i] <- 'N'
    }
    if(grepl(taxonomy_UKSurvey$class[i], taxonomy_UKSurvey$`Most common class-level assignment`[i], fixed = TRUE) == TRUE){
      taxonomy_UKSurvey$class.match[i] <- 'Y'
    }else{
      taxonomy_UKSurvey$class.match[i] <- 'N'
    }
    if(grepl(taxonomy_UKSurvey$phylum[i], taxonomy_UKSurvey$`Most common phylum-level assignment`[i], fixed = TRUE) == TRUE){
      taxonomy_UKSurvey$phylum.match[i] <- 'Y'
    }else{
      taxonomy_UKSurvey$phylum.match[i] <- 'N'
    }
}


#Add column for whether TBAS and RDP
taxonomy_UKSurvey[,c('phylum.match.1','class.match.1','order.match.1','family.match.1','genus.match.1','species.match.1')] <- NA

#TBAS vs RDP
for (i in 1:nrow(taxonomy_UKSurvey)) {
  #if blast species matches TBAS species then add Y(es) to match column else add N(o)
  if(grepl(taxonomy_UKSurvey$rdp_species[i], taxonomy_UKSurvey$`Taxon assignment`[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$species.match.1[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$species.match.1[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$rdp_genus[i], taxonomy_UKSurvey$`Most common Genus-level assignment`[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$genus.match.1[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$genus.match.1[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$rdp_family[i], taxonomy_UKSurvey$`Most common family-level assignment`[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$family.match.1[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$family.match.1[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$rdp_order[i], taxonomy_UKSurvey$`Most common order-level assignment`[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$order.match.1[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$order.match.1[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$rdp_class[i], taxonomy_UKSurvey$`Most common class-level assignment`[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$class.match.1[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$class.match.1[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$rdp_phylum[i], taxonomy_UKSurvey$`Most common phylum-level assignment`[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$phylum.match.1[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$phylum.match.1[i] <- 'N'
  }
}


#dnabarcoder vs RDP
taxonomy_UKSurvey[,c('phylum.match.2','class.match.2','order.match.2','family.match.2','genus.match.2','species.match.2')] <- NA
for (i in 1:nrow(taxonomy_UKSurvey)) {
  #if blast species matches TBAS species then add Y(es) to match column else add N(o)
  if(grepl(taxonomy_UKSurvey$species[i], taxonomy_UKSurvey$rdp_species[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$species.match.2[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$species.match.2[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$genus[i], taxonomy_UKSurvey$rdp_genus[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$genus.match.2[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$genus.match.2[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$family[i], taxonomy_UKSurvey$rdp_family[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$family.match.2[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$family.match.2[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$order[i], taxonomy_UKSurvey$rdp_order[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$order.match.2[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$order.match.2[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$class[i], taxonomy_UKSurvey$rdp_class[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$class.match.2[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$class.match.2[i] <- 'N'
  }
  if(grepl(taxonomy_UKSurvey$phylum[i], taxonomy_UKSurvey$rdp_phylum[i], fixed = TRUE) == TRUE){
    taxonomy_UKSurvey$phylum.match.2[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$phylum.match.2[i] <- 'N'
  }
}


#BLAST vs RDP vs TBAS
taxonomy_UKSurvey[,c('phylum.match.3','class.match.3','order.match.3','family.match.3','genus.match.3','species.match.3')] <- NA
for (i in 1:nrow(taxonomy_UKSurvey)) {
  #if blast species matches TBAS species then add Y(es) to match column else add N(o)
  if((taxonomy_UKSurvey$species[i] == taxonomy_UKSurvey$rdp_species[i] & taxonomy_UKSurvey$species[i] == taxonomy_UKSurvey$`Taxon assignment`[i]) == TRUE){
    taxonomy_UKSurvey$species.match.3[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$species.match.3[i] <- 'N'
  }
  if((taxonomy_UKSurvey$genus[i] == taxonomy_UKSurvey$rdp_genus[i] & taxonomy_UKSurvey$genus[i] == taxonomy_UKSurvey$`Most common Genus-level assignment`[i]) == TRUE){
    taxonomy_UKSurvey$genus.match.3[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$genus.match.3[i] <- 'N'
  }
  if((taxonomy_UKSurvey$family[i] == taxonomy_UKSurvey$rdp_family[i] & taxonomy_UKSurvey$family[i] == taxonomy_UKSurvey$`Most common family-level assignment`[i]) == TRUE){
    taxonomy_UKSurvey$family.match.3[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$family.match.3[i] <- 'N'
  }
  if((taxonomy_UKSurvey$order[i] == taxonomy_UKSurvey$rdp_order[i] & taxonomy_UKSurvey$order[i] == taxonomy_UKSurvey$`Most common order-level assignment`[i]) == TRUE){
    taxonomy_UKSurvey$order.match.3[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$order.match.3[i] <- 'N'
  }
  if((taxonomy_UKSurvey$class[i] == taxonomy_UKSurvey$rdp_class[i] & taxonomy_UKSurvey$class[i] == taxonomy_UKSurvey$`Most common class-level assignment`[i]) == TRUE){
    taxonomy_UKSurvey$class.match.3[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$class.match.3[i] <- 'N'
  }
  if((taxonomy_UKSurvey$phylum[i] == taxonomy_UKSurvey$rdp_phylum[i] & taxonomy_UKSurvey$phylum[i] == taxonomy_UKSurvey$`Most common phylum-level assignment`[i]) == TRUE){
    taxonomy_UKSurvey$phylum.match.3[i] <- 'Y'
  }else{
    taxonomy_UKSurvey$phylum.match.3[i] <- 'N'
  }
}

#count how ASVs ID agree at each taxon rank
dnabarcoder_vs_TBAS_agree_columns <- colnames(taxonomy_UKSurvey)[c(38:43)]
dnabarcoder_vs_TBAS <- c()
for(col in dnabarcoder_vs_TBAS_agree_columns){
  dnabarcoder_vs_TBAS <- c(dnabarcoder_vs_TBAS,length(grep('Y',taxonomy_UKSurvey[[col]]))/length(taxonomy_UKSurvey[[col]]))
}

TBAS_vs_RDP_agree_columns <- colnames(taxonomy_UKSurvey)[c(44:49)]
TBAS_vs_RDP <- c()
for(col in TBAS_vs_RDP_agree_columns){
  TBAS_vs_RDP <- c(TBAS_vs_RDP,length(grep('Y',taxonomy_UKSurvey[[col]]))/length(taxonomy_UKSurvey[[col]]))
}

dnabarcoder_vs_RDP_agree_columns <- colnames(taxonomy_UKSurvey)[c(50:55)]
dnabarcoder_vs_RDP <- c()
for(col in dnabarcoder_vs_RDP_agree_columns){
  dnabarcoder_vs_RDP <- c(dnabarcoder_vs_RDP,length(grep('Y',taxonomy_UKSurvey[[col]]))/length(taxonomy_UKSurvey[[col]]))
}

dnabarcoder_vs_RDP_vs_TBAS_agree_columns <- colnames(taxonomy_UKSurvey)[c(56:61)]
dnabarcoder_vs_RDP_vs_TBAS <- c()
for(col in dnabarcoder_vs_RDP_vs_TBAS_agree_columns){
  dnabarcoder_vs_RDP_vs_TBAS <- c(dnabarcoder_vs_RDP_vs_TBAS,length(grep('Y',taxonomy_UKSurvey[[col]]))/length(taxonomy_UKSurvey[[col]]))
}

#combine into table
props_agree_table <- tibble(dnabarcoder_vs_TBAS,TBAS_vs_RDP,dnabarcoder_vs_RDP,dnabarcoder_vs_RDP_vs_TBAS, c('phylum','class','order','family','genus','species')) %>%
  pivot_longer(1:4)

(agree_plot <- ggplot(data = props_agree_table, aes(x = `c("phylum", "class", "order", "family", "genus", "species")`, y = value*100, col = name, group = name)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(limits = c('phylum','class','order','family','genus','species')) +
  ylab('Taxonomy ID agreement (%)') +
  xlab('Taxon rank') +
  theme(aspect.ratio=1) +
  ylim(c(0,100)) +
  labs(col = "Comparison"))


#plot percentage of identified taxa at each rank
#TBAS
TBAS_columns <- colnames(taxonomy_UKSurvey)[c(2,4,6,8,10,12)]
TBAS <- c()
for(col in TBAS_columns){
  TBAS <- c(TBAS,length(grep('\\.|\\,',taxonomy_UKSurvey[[col]]))/length(taxonomy_UKSurvey[[col]]))
}
#RDP
rdp_columns <- colnames(taxonomy_UKSurvey)[c(20:25)]
RDP <- c()
for(col in rdp_columns){
  RDP <- c(RDP,length(grep('unidentified| sp',taxonomy_UKSurvey[[col]]))/length(taxonomy_UKSurvey[[col]]))
}
#dnabarcoder
rank_order <- c('phylum','class','order','family','genus','species')
dnabarcoder <- taxonomy_UKSurvey %>% count(rank) %>%
  mutate(rank = factor(rank, levels = rank_order)) %>%
  arrange(rank) %>%
  mutate(cumulative = rev(cumsum(rev(n))), prop = cumulative/nrow(taxonomy_UKSurvey)) 

#combine into table
props_table <- tibble(TBAS = 1-TBAS,RDP = 1-RDP, dnabarcoder = dnabarcoder$prop, rank_order) %>%
  pivot_longer(1:3)

(id_plot <- ggplot(data = props_table, aes(x = rank_order, y = value*100, col = name, group = name)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(limits = c('phylum','class','order','family','genus','species')) +
  ylab('OTUs identified (%)') +
  xlab('Taxon rank') +
  theme(aspect.ratio=1) +
  ylim(c(0,100)) +
  labs(col = "Tool"))

library(cowplot)
png("../FIGURES/prop_OTUs_ID_taxonomy_and_agreement_UKSurvey.png",  res = 400, width = 2480, height = 2480)
plot_grid(id_plot,agree_plot, labels = c('(b)', '(c)'), ncol = 1, align = T)
dev.off()
