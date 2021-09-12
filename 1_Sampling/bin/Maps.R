#https://www.diegovalle.net/mxmaps/
# Verónica Reyes, febrero 2020
# Create CDMX´s map

# Install
#if (!require("devtools")) {
#  install.packages("devtools")
#}
#devtools::install_github("diegovalle/mxmaps")

library("mxmaps")
library(ggplot2)


library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering, height = unit(1.25, "cm")) +
  coord_sf(xlim = c(-116.15, -75.12), ylim = c(10.65, 33.97))+
  theme_bw()+
  annotate("point", x = -99.301, y = 19.285, colour = "#4cacdb", size = 2)+
  labs(x="", y="")
ggsave("../outputs/1_mexico_map.png")

data("df_mxmunicipio")
df_mxmunicipio$value <- as.numeric(df_mxmunicipio$municipio_name =="La Magdalena Contreras")


mxmunicipio_choropleth(df_mxmunicipio, num_colors = 2,
                       zoom = subset(df_mxmunicipio, state_abbr_official %in% "DF")$region) 

ggsave("../outputs/1_CDMX_map.png")




