library(tidyverse)
library(RColorBrewer)

data <- readRDS('data/tcga_pancan_cna.rds')

#print(data %>% filter(total_samp==253) %>% select(cancer_type) %>% unique())
data <- data %>% 
	select(cancer_type, total_samp, hugoGeneSymbol, n, alteration) %>% 
	group_by(cancer_type, hugoGeneSymbol, alteration) %>% 
	summarize(n=sum(n), total_samp=unique(total_samp)) %>% 
	ungroup() %>% 
	group_by(cancer_type, alteration) %>% 
	summarize(n=length(unique(hugoGeneSymbol)), total_samp=unique(total_samp)) %>%
	pivot_wider(names_from = 'alteration', values_from = 'n')
print(data)
colourCount <- nrow(data)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
data <- data %>% filter(total_samp!=253)
gg<-ggplot(data)+
 geom_bar(mapping=aes(x=cancer_type,y=total_samp,fill=cancer_type), stat='identity')+
 scale_fill_manual(values = getPalette(colourCount))+
 geom_text(mapping=aes(x=cancer_type,y=total_samp, label=del), color='#b5265f', size=2.5, vjust=-0.5)+
 geom_text(mapping=aes(x=cancer_type,y=total_samp, label=dup), color='#449e26',size=2.5, vjust =-2)+
 theme_bw()+
 xlab('')+
 ylab('N Samples')+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
 guides(fill='none')
ggsave(filename='figs/tcgastats.png', width=3000,height=2500,units='px')


