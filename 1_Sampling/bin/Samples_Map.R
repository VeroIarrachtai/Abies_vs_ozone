library(maps)
library(maptools)

# library(mapdata);tambien puedo utilizar archivos del portal de bioinformación de la CONABIO(http://www.conabio.gob.mx/informacion/gis/)


#Cargar archivo CONABIO y coordenadas dentro del mapa CONABIO
rbiogeo<-readShapePoly("../metadata/rbiog4mgw_p1/rbiog4mgw.shp") #Archivo CONABIO ".shp"

#Plotear regiones biogeograficas
plot(rbiogeo, ) #Solo plotea las regiones b

levels(rbiogeo$PROVINCIA) #Enlista las regiones biogeograficas

mycols<-c(rep("white", 7), "mediumpurple1", rep("white", 5), "gold", "darkolivegreen1", rep("white", 4)) #Elijo el color de cada region
palette(mycols)
palette() #Me muestra el patron del coloreado


plot(rbiogeo, col=rbiogeo$PROVINCIA)

#Agregar puntos de coordenadas

# Crear una tabla con todos los datos. Sigo en "/Users/geyev15/Documents/Homeworks_and_Labworks/Bitacora/SEPTIEMBRE/Mapas_con_R/rbiog4mgw_p1"
# caracteristicas<-read.csv("../Coordenadas_ju_myb.csv") #NO PUDE ASI
caracteristicas<-read.delim("../Coordenadas_ju_myb") #Caracteristicas de poblaciones de Juniperus

caracteristicas$ID #Checar algun dato de la tabla(ID,Longitud, Latitud, etc)


#Agregar puntos

# points(caracteristicas$Longitud, caracteristicas$Latitud,pch=17, col="cyan")
# triangulo= J. monticola y circulo J. blancoi; azul=baja altitud y rosa=alta altitud

palette(c("deeppink3", "navy"))
points(caracteristicas$Longitud, caracteristicas$Latitud, pch=c(rep(17,13), rep(19, 8)), col=caracteristicas$Habitat)

#Escala del mapa
map.scale(ratio = FALSE)