## PLOT MDV-POSITIVE CHICKENS ##


library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("readxl")
library("dplyr")
library("tidyr")
library("ggrepel")
library("viridis")


setwd("~/Desktop/MDV/locations_figure/")


locations <- read_excel(path = "20220308_from_antony_thesis_samples_reads_updated_for_rebuttal.xlsx")
locationsrm_rpts <- subset(locations, Sample != "OL1389" & Sample != "OL1986" & Sample != "OL1987" & Sample != "OL1999" & Sample != "OL2268" & Sample != "OL2272")
locationsrm_rpts$Short_loc_cov_final <- gsub('_', '\n', locationsrm_rpts$Short_loc_cov_final)
#locations <- subset(locations, Sample != "OL2178")
#Split up the comma-limited coordinates into two separate columns of lat/long
#locations <- separate(locations,col="Safe_coordinates", into=c("lat","long"), sep=",")
#Get rid of NA values
#locations <- subset(locations,locations$lat != "N/A" | locations$long != "N/A")
#Set the populations as factors so they remain in the same order
#locations$Population <- factor(locations$Population, levels = c("Africa", "EastAsia", "NorthAmerica", "SouthAsia", "WestAsia", "bankiva", "gallus", "murghi", "jabouillei", "spadiceus", "lafayetii", "sonnerattii", "varius"  ))
#Get the world data to plot the map
world <- ne_countries(scale = "medium", returnclass = "sf")

#Get a list of 13 colours to be used to colour the points. 
#library(RColorBrewer)
#col <- colorRampPalette(c("black", "blue", "green", "red", "purple", "orange"))(12)


plot <- ggplot(data = world)+ 
  geom_sf(color = "grey10", fill= "antiquewhite", size = 0.1) + 
  #Set the map limits
  coord_sf(xlim = c(-10, 80), ylim = c(30,80), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "slategray1")) +
  #Use geom_count - this factors in the overplotting and increases the size of the point based on the number of samples
  #geom_count(shape = 20, data = locations,  aes(x=as.numeric(Long) , y=as.numeric(Lat), colour = Average_sequencing_depth), position = "jitter") +
  #scale_colour_manual(values = col)+
  theme(text = element_text(size=12), 
        axis.text = element_text(colour = "black")) +
  scale_colour_viridis_b()+
  geom_text_repel(data = locationsrm_rpts, aes(x=Long, y=Lat, label=Short_location), segment.size = 0.5, size = 4, force = 50, arrow = arrow(length = unit(0.01, "npc")), box.padding = 1)
  

plot

#ggsave(filename = "20220308_plot_MDV_locations.svg", plot = plot, width=8, height=8)
#ggsave(filename = "20220308_plot_MDV_locations.png", plot = plot, width=8, height=8)

######SECOND PLOT######
########DATES##########

dates_known <- filter(locations, Numeric_date > 0)
dates_unknown <- filter(locations, Numeric_date < 0)



plot2 <- ggplot(data = world) +
  geom_sf(color = "grey90", fill = "grey90", size = 0.1) +
  #Set the map limits
  coord_sf(xlim = c(-10, 80), ylim = c(30,70), expand = FALSE) + 
  #xlab("Longitude") + 
  #ylab("Latitude") + 
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "white")) +
  geom_point(shape = 21, size = 3, alpha = 0.7, colour = "black", data = dates_known, aes(x=as.numeric(Long), y=as.numeric(Lat), fill = Numeric_date)) +
  #geom_point(shape = 21, size = 3, alpha = 0.7, colour  = "black", fill = "gray45", data = dates_unknown, aes(x=as.numeric(Long), y=as.numeric(Lat))) +
  #theme(text = element_text(size=8), 
  #      axis.text = element_text(colour = "black")) +
  #scale_fill_viridis(option = "C")+
  scale_fill_gradient(low = "red", high = "white")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  geom_text_repel(data = locationsrm_rpts, aes(x=Long, y=Lat, label=Short_loc_cov_final), lineheight = 0.8, segment.size = 0.2, size = 2.5, force = 20, arrow = arrow(length = unit(0.01, "npc")), box.padding = 1, max.overlaps = 12)+
  labs(fill="Sample Date")


plot2

ggsave(filename = "20230310_plot_MDV.svg", plot = plot2, width=6.5, height=6.5)
ggsave(filename = "20230310_plot_MDV.png", plot = plot2, width=6.5, height=6.5)




basemap <- ggplot(data = world) +
  geom_sf(color = "grey90", fill = "grey90", size = 0.1) +
  #Set the map limits
  coord_sf(xlim = c(-10, 80), ylim = c(30,70), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "white"))
  #geom_point(shape = 21, size = 3, alpha = 0.7, colour = "black", data = dates_known, aes(x=as.numeric(Long), y=as.numeric(Lat), fill = Numeric_date)) +
  #geom_point(shape = 21, size = 3, alpha = 0.7, colour  = "black", fill = "gray45", data = dates_unknown, aes(x=as.numeric(Long), y=as.numeric(Lat))) +
  #theme(text = element_text(size=8), 
  #      axis.text = element_text(colour = "black")) +
  #scale_fill_viridis(option = "C")+
  #scale_fill_gradient(low = "red", high = "white")+
  #geom_text_repel(data = locationsrm_rpts, aes(x=Long, y=Lat, label=Short_loc_cov_final), lineheight = 0.8, segment.size = 0.2, size = 2.5, force = 20, arrow = arrow(length = unit(0.01, "npc")), box.padding = 1, max.overlaps = 12)+
  #labs(fill="Sample Date")


basemap
ggsave(filename = "20221117_basemap.svg", plot = basemap, width=6.5, height=6.5)


  
