library(ggplot2)
library(reshape)
library(wesanderson)
library(ggforce)
library(dplyr)
require("scales")
library(gridExtra)

theme_Publication_simple <- function(base_size=14, base_family="Helvetica", xangle=0, 
                                     position="bottom", direction="horizontal", 
                                     legend.title=element_text(face="italic"),
                                     legend.text=element_text()) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.text.x = element_text(angle=xangle,hjust=0.5),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = position,
            legend.direction = direction,
            legend.key.size= unit(1, "cm"),
            legend.key.width = unit(1.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = legend.title,
            legend.text = legend.text,
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

setwd('/Users/edimopoulos/Marek/second_round_captures')

sfs <- read.csv("SFS_res.txt", header=TRUE, sep='\t')


ggplot(sfs, aes(x=Number_of_SNPs, y=Number_of_sites)) + geom_bar(stat="identity") + 
  ylab('Number of sites') + xlab('Number of SNPs') +
  theme_Publication_simple(base_size = 20, base_family = "Helvetica") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     breaks =c(0, 5, 50, 500, 5000, 15000, 90000)) 












