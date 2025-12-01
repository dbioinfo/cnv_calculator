library(tidyverse)
library(RColorBrewer)

options(repr.plot.width=25, repr.plot.height=12)

data <- read_csv('data/seerstats.csv')

colourCount <- nrow(data)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

gg<-ggplot(data)+
 geom_bar(mapping=aes(x=cancer_type,y=total_incidence,fill=cancer_type), stat='identity')+
 scale_fill_manual(values = getPalette(colourCount))+
 geom_text(mapping=aes(x=cancer_type,y=total_incidence, label=male), color='#3b91ed',vjust=-0.5, size=3)+
 geom_text(mapping=aes(x=cancer_type,y=total_incidence, label=female), color='#f542d7',vjust=-2, size=3)+
 theme_bw()+
 scale_y_continuous(trans='log2')+
 xlab('')+
 ylab('Incidence')+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
 guides(fill='none')
ggsave(filename='figs/seerstats.png', width=3000,height=2500,units='px')
