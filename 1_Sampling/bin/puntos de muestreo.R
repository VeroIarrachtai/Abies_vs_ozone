## 14 Enero
## Verónica Reyes
## Crear mapas con R

#Estoy en HOME ~, debo elegir una carpeta para ordenar lo generado.
#Trabajaré en la carpeta de Bitacora

setwd("..//Mapas_con_R/rbiog4mgw_p1") #Direccion completa para cambiar el wd
getwd() #Confirmar el wd

#Cargar paqueterias
library(maps)
library(maptools)
# library(mapdata);tambien puedo utilizar archivos del portal de bioinformación de la CONABIO(http://www.conabio.gob.mx/informacion/gis/)


#Cargar archivo CONABIO y coordenadas dentro del mapa CONABIO
rbiogeo<-readShapePoly("rbiog4mgw") #Archivo CONABIO ".shp"

#Plotear regiones biogeograficas
plot(rbiogeo, ) #Solo plotea las regiones b.
levels(rbiogeo$PROVINCIA) #Enlista las regiones biogeograficas

#Colocar colores, mediante un patron de repeticion.

mycols<-c(rep("white", 7), "mediumorchid3", rep("white", 5), "orange", "yellowgreen", rep("white", 4)) #Elijo el color de cada region
palette(mycols)
palette() #Me muestra el patron del coloreado
plot(rbiogeo, col=rbiogeo$PROVINCIA)


# Crear una tabla con todos los datos. Sigo en "/Users/geyev15/Documents/Homeworks_and_Labworks/Bitacora/SEPTIEMBRE/Mapas_con_R/rbiog4mgw_p1"
# caracteristicas<-read.csv("../Coordenadas_ju_myb.csv") #NO PUDE ASI
caracteristicas<-read.delim("../Coordenadas_ju_myb") #Caracteristicas de poblaciones de Juniperus

caracteristicas$ID #Checar algun dato de la tabla(ID,Longitud, Latitud, etc)

#Agregar puntos

# points(caracteristicas$Longitud, caracteristicas$Latitud,pch=17, col="cyan")
# triangulo= J. monticola y circulo J. blancoi; azul=baja altitud y rosa=alta altitud

palette(c("red", "blue"))
points(caracteristicas$Longitud, caracteristicas$Latitud, pch=c(rep(20,18), rep(17, 18)), col=caracteristicas$Habitat)

#Escala del mapa
map.scale(ratio = FALSE)
