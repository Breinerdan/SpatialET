#SpatialET
#Código simple para calcular ETP y ETR por diversos métodos

#Se implementa de manera espacial los métodos empíricos de Hargreaves, 
#Thornthwaite, Turc, Cenicafé y Countagne, además de la conversión a ETR usando Budyko

#___________________________________________________________________________________
#Primero cargamos las librerías necesarias

library(raster) #cargar la libreria raster
library(rgdal) #cargar la libreria rgdal
library(sf) #cargar la libreria sf
library(sp) #cargar la libreria sp
library(spData) #cargar la libreria sp

path<-"C:/RScript/SpatialET"
setwd(path) #localizar el directorio de trabajo
#___________________________________________________________________________________
# Definimos dos funciones que se usarán repetidamente a lo largo del código
# una función para graficar las entradas y salidas y una función para convertir la
# Evapotranspiración potencial en real, usando la expresión de Budyko

# Primero definimos la función de Budyko

ETR_Budyko <- function(ETPf,Precipf) {
  ETRf <- sqrt(ETPf*Precipf*tanh(Precipf/ETPf)*(1-cosh(ETPf/Precipf)+sinh(ETPf/Precipf)))
  return (ETRf)
}

# Ahora definimos la función de graficar
# Con esta función graficamos la variable anual con la zona de estudio
graficar <-function(variable_anual,nombre_variable) {
  plot(variable_anual, xlab="Este (m)", ylab="Norte (m)",
       main =nombre_variable)
  plot(st_geometry(Studyzone), lwd =2, axes=TRUE, add=TRUE)
  legend('bottomleft',legend = "Zona de estudio", col='black',lty=1,lwd=2 )
}

#___________________________________________________________________________________

#Asignar las variables de trabajo
#Leemos el shape de la zona de estudio
Studyzone <- st_read("00_Zona_estudio/Zona_Estudio.shp")
#Leemos el Modelo Digital de Elevación - MDT
High<-raster("01_Topografía/mdt.tif")
#Podemos ir visualizando nuestros datos para asegurar un adecuado cargue de los mismos
graficar(High, "Elevación del terreno (m.s.n.m)")
#___________________________________________________________________________________
#Leemos las precipitaciones mensuales multianuales
setwd("02_Precip/")
files_prec <- list.files(pattern = '*.tif')
Monthly_Prec <- list(length(files_prec)) # Este será nuestra lista raster que almacena la precipitación mensual multianual
for (i in 1:length(files_prec)) {
  prec <- raster(files_prec[i])
  Monthly_Prec[i] <- prec 
}
Monthly_Prec <- stack(Monthly_Prec) # La lista raster se convierte en un stack para facilitar su manejo
plot(Monthly_Prec) # Graficamos el ciclo anual de precipitación para verificar un adecuado cargue de los datos
# Calculamos la precipitación multianual
Annual_Prec <- calc(Monthly_Prec,sum)

# Graficamos este resultado:
graficar(Annual_Prec, "Precipitación anual (mm/año)")
#___________________________________________________________________________________

#Leemos las temperaturas medias mensuales multianuales
setwd("C:/RScript/SpatialET/03_Temp")
files_temp <- list.files(pattern = '*.tif')
Monthly_temp <- list(length(files_temp)) # Este será nuestra lista raster que almacena la temperatura media mensual multianual
for (i in 1:length(files_temp)) {
  temp <- raster(files_temp[i])
  Monthly_temp[i] <- temp 
}
Monthly_temp <- stack(Monthly_temp) # La lista raster se convierte en un stack para facilitar su manejo
plot(Monthly_temp) # Graficamos el ciclo anual de temperatura para verificar un adecuado cargue de los datos
# Calculamos la temperatura media multianual
Annual_temp <- calc(Monthly_temp,mean)

# Graficamos este resultado:
graficar(Annual_temp,"Temperatura media anual (°C)")

#___________________________________________________________________________________

#Leemos la radiación solar incidente mensual multianual
setwd("C:/RScript/SpatialET/04_RadS_Inc")
files_rad <- list.files(pattern = '*.tif')
Monthly_rad <- list(length(files_rad)) # Este será nuestra lista raster que almacena la radiación solar incidente mensual multianual
for (i in 1:length(files_rad)) {
  rad <- raster(files_rad[i])
  Monthly_rad[i] <- rad
}
Monthly_rad <- stack(Monthly_rad) # La lista raster se convierte en un stack para facilitar su manejo
plot(Monthly_rad) # Graficamos el ciclo anual de temperatura para verificar un adecuado cargue de los datos
# Calculamos la temperatura media multianual
Annual_rad <- calc(Monthly_rad,mean)

# Graficamos este resultado:
graficar(Annual_rad, "Radiación solar incidente anual (kJ m-2 day-1)")

#___________________________________________________________________________________


#Aqui se empieza a aplicar los métodos

# 1) Método Cenicafé - Budyko

ETP_cc<-(1700.17*exp(-0.0002*High)) # Esta es la ETP por Cenicafe
#Aquí graficamosn el resultado con fines de verificación
graficar(ETP_cc, "Evapotranspiración potencial - Cenicafé (mm/año)")
ETP_cc<-resample(ETP_cc,Annual_Prec) # consitencia espacial entre ETP y precipitación
#aqui escribimos la ETP por Cenicafe
writeRaster(ETP_cc,"C:/RScript/SpatialET/06_resultados/ETP/anuales/ETP_cc.tif",overwrite=TRUE) 

#Aquí calculamos la ETR por Cenicafé - Budyko:
ETR_cb<-ETR_Budyko(ETP_cc,Annual_Prec)
#Aquí graficamosn el resultado con fines de verificación
graficar(ETR_cb,"Evapotranspiración real - Cenicafé - Budyko (mm/año)")
#aqui escribimos la ETR por Cenicafe-Budyko
writeRaster(ETR_cb,"C:/RScript/SpatialET/06_resultados/ETR/anuales/ETR_cb.tif",overwrite=TRUE) 

#___________________________________________________________________________________

# 2) Método de Turc
L<-300+25*Annual_temp+0.05*(Annual_temp^3) # Calculo del parametro heliotermico de Turc
L<-resample(L,Annual_Prec) # consistencia espacial con la precipitación
ETR_T<-(Annual_Prec/sqrt(0.9+(Annual_Prec/L)^2))
ETR_T[Annual_Prec/L<0.316]<-Annual_Prec # Calculo condicional de ETR por Turc
# Aqui se grafica la ETR por Turc
graficar(ETR_T, "Evapotranspiración real - Turc (mm/año)")
writeRaster(ETR_T,"C:/RScript/SpatialET/06_resultados/ETR/anuales/ETR_T.tif",overwrite=TRUE) #aqui escribimos la ETR por Turc

#___________________________________________________________________________________

# 3) Método de Thornthwaite

# Primero se calcula el factor 
coef_I<-12*((Annual_temp/5)^1.514) # calculo coeficiente I
# calculo coeficiente a
coef_a<-0.000000675*(coef_I^3)-0.0000771*(coef_I^2)+0.0179*(coef_I)+0.492 

# Ahora aplicamos la ecuación de Thornthwaite mensual:
Monthly_ETP_twt <- list(length(files_temp)) # Este será nuestra lista raster que almacena la ETP mensual multianual calculada
for (i in 1:length(files_temp)){
 ETP<-16*((10*Monthly_temp[[i]]/coef_I)^coef_a)
 Monthly_ETP_twt[i]<-ETP
}
Monthly_ETP_twt <- stack(Monthly_ETP_twt) # La lista raster se convierte en un stack para facilitar su manejo
plot(Monthly_ETP_twt) # Graficamos el ciclo anual de ETP para verificar los cálculos
# Calculamos la temperatura media multianual
ETP_twt <- calc(Monthly_ETP_twt,sum)

#Aquí graficamosn el resultado con fines de verificación
graficar(ETP_twt, "Evapotranspiración potencial - Thornthwaite (mm/año)")
ETP_twt<-resample(ETP_twt,Annual_Prec) # consitencia espacial entre ETP y precipitación
#aqui escribimos la ETP por Thornthwaite
writeRaster(ETP_twt,"C:/RScript/SpatialET/06_resultados/ETP/anuales/ETP_twt.tif",overwrite=TRUE) 

#Aquí calculamos la ETR por Thornthwaite - Budyko:
  ETR_twt<-ETR_Budyko(ETP_twt,Annual_Prec) 
#Aquí graficamosn el resultado con fines de verificación
graficar(ETR_twt,"Evapotranspiración real - Thornthwaite - Budyko (mm/año)")
#aqui escribimos la ETR por Thornthwaite-Budyko
writeRaster(ETR_cb,"C:/RScript/SpatialET/06_resultados/ETR/anuales/ETR_twt.tif",overwrite=TRUE) 

#___________________________________________________________________________________

# 4) Método de Hargreaves

# Ahora aplicamos la ecuación de Hargreaves mensual:
Monthly_ETP_hg <- list(length(files_rad)) # Este será nuestra lista raster que almacena la ETP mensual multianual calculada
for (i in 1:length(files_rad)){
  ETP<-0.0135*31*(Monthly_temp[[i]]+17.78)*Monthly_rad[[i]]/2450
  Monthly_ETP_hg[i]<-ETP
}
Monthly_ETP_hg <- stack(Monthly_ETP_hg) # La lista raster se convierte en un stack para facilitar su manejo
plot(Monthly_ETP_hg) # Graficamos el ciclo anual de ETP para verificar los cálculos
# Calculamos la temperatura media multianual
ETP_hg <- calc(Monthly_ETP_hg,sum)

#Aquí graficamosn el resultado con fines de verificación
graficar(ETP_hg, "Evapotranspiración potencial - Hargreaves (mm/año)")
ETP_hg<-resample(ETP_hg,Annual_Prec) # consitencia espacial entre ETP y precipitación
#aqui escribimos la ETP por Thornthwaite
writeRaster(ETP_hg,"C:/RScript/SpatialET/06_resultados/ETP/anuales/ETP_hg.tif",overwrite=TRUE) 

#Aquí calculamos la ETR por Hargreaves - Budyko:
ETR_hg<-ETR_Budyko(ETP_hg,Annual_Prec) 
#Aquí graficamosn el resultado con fines de verificación
graficar(ETR_hg,"Evapotranspiración real - Hargreaves - Budyko (mm/año)")
#aqui escribimos la ETR por Hargreaves-Budyko
writeRaster(ETR_hg,"C:/RScript/SpatialET/06_resultados/ETR/anuales/ETR_hg.tif",overwrite=TRUE) 

#___________________________________________________________________________________

# 5) Método de Countagne

lambda<-1/(0.8+0.14*Annual_temp)
P<-0.001*Annual_Prec
P<-resample(P,Annual_temp)
P<-crop(P,Annual_temp)
Einf<-1/(8*lambda)
Esup<-1/(2*lambda)
Emin<-1/(4*lambda)
ETR_cou<-P-lambda*((P)^2)
ETR_cou[P<Einf]<-P
ETR_cou[P>Esup]<-Emin[P>Esup]
ETR_cou <-1000*ETR_cou
#Aquí graficamosn el resultado con fines de verificación
graficar(ETR_cou,"Evapotranspiración real - countagne (mm/año)")
#aqui escribimos la ETR por Countagne
writeRaster(ETR_hg,"C:/RScript/SpatialET/06_resultados/ETR/anuales/ETR_cou.tif",overwrite=TRUE) 

#___________________________________________________________________________________

# 6) Método externo: en este caso se carga la ETP mensual multianual proveniente de otras fuentes (p.ej renaálisis climáticos)
# y en este código se calcula la ETR resultante, considerando la precipitación utilizada en los otros métodos

setwd("C:/RScript/SpatialET/05_ETP_Alt")
files_EtpAlt <- list.files(pattern = '*.tif')
Monthly_EtpAlt <- list(length(files_rad)) # Este será nuestra lista raster que almacena la ETP mensual multianual ingresada
for (i in 1:length(files_EtpAlt)) {
  EtpAlt <- raster(files_EtpAlt[i])
  Monthly_EtpAlt[i] <- EtpAlt
}
Monthly_EtpAlt <- stack(Monthly_EtpAlt) # La lista raster se convierte en un stack para facilitar su manejo
plot(Monthly_EtpAlt) # Graficamos el ciclo anual de ETP para verificar un adecuado cargue de los datos
# Calculamos la ETP multianual del método alterno
ETP_alt<- calc(Monthly_EtpAlt,sum)
# Graficamos este resultado:
graficar(ETP_alt, "Evapotranspiración potencial - alterna (mm/año)")
ETP_alt<-resample(ETP_alt,Annual_Prec) # consitencia espacial entre ETP y precipitación
#aqui escribimos la ETP multianual
writeRaster(ETP_alt,"C:/RScript/SpatialET/06_resultados/ETP/anuales/ETP_alt.tif",overwrite=TRUE) 

#Aquí convertimos la ETP alterna en ETR por Budyko
  ETR_alt<-ETR_Budyko(ETP_alt,Annual_Prec) 
#Aquí graficamosn el resultado con fines de verificación
graficar(ETR_alt,"Evapotranspiración real - alterna - Budyko (mm/año)")
#aqui escribimos la ETR alterna por Budyko
writeRaster(ETR_hg,"C:/RScript/SpatialET/06_resultados/ETR/anuales/ETR_alt.tif",overwrite=TRUE) 

#___________________________________________________________________________________

# Finalmente, podemos calcular el índice de Aridez definido por el IDEAM (2010), usando los diversos métodos
# Creamos una función para tal fin

#Las opciones son:
#Hargreaves
#Thornthwaite
#Cenicafe

Indice_aridez <- function(metodo) {
  if (metodo == "Hargreaves"){
    IA <-(ETP_hg-ETR_hg)/ETP_hg
    graficar(IA,"Índice Aridez - Hargreaves")
    ruta <- paste("C:/RScript/SpatialET/06_resultados/Indice_Aridez/IA",metodo)
    writeRaster(IA,ruta,format="GTiff",overwrite=TRUE)
  } else if (metodo == "Cenicafe"){
    IA <-(ETP_cc-ETR_cb)/ETP_cc
    graficar(IA,"Índice Aridez - Cenicafé")
    ruta <- paste("C:/RScript/SpatialET/06_resultados/Indice_Aridez/IA",metodo)
    writeRaster(IA,ruta,format="GTiff",overwrite=TRUE)
  }else if
  (metodo == "Thornthwaite"){
    IA <-(ETP_twt-ETR_twt)/ETP_twt
    graficar(IA,"Índice Aridez - Thornthwaite")
    ruta <- paste("C:/RScript/SpatialET/06_resultados/Indice_Aridez/IA",metodo)
    writeRaster(IA,ruta,format="GTiff",overwrite=TRUE)
  }else 
  {print("especifique un método valido:Hargreaves,Thornthwaite o Cenicafe")}}

Indice_aridez("Cenicafe") # aquí calculamos el IA de Cenicafé, visualizamos y exportamos
Indice_aridez("Hargreaves") # aquí calculamos el IA de Hargreaves, visualizamos y exportamos
Indice_aridez("Thornthwaite") # aquí calculamos el IA de Hargreaves, visualizamos y exportamos

#___________________________________________________________________________________
#Una opción importante es poder exportar los mapas de ETP mensuales de interés, para ello creamos una función
#De forma tal que solo se deba especificar que ETP exportar. Es importante aclarar que esto solo para la ETP, ya que la ETR
#mensual no puede calcularse por Budyko, debería clacularse por otros métodos como modelación hidrológica o balance hídrico

nombres<-c("01ETP", "02ETP","03ETP", "04ETP", "05ETP", "06ETP",
           "07ETP", "08ETP", "09ETP", "10ETP", "11ETP", "12ETP")

Exportar_mensual <- function(ETPm) {
 for (i in 1:length(ETPm)) {
   mesi <-ETPm[[i]]
  ruta <- paste("C:/RScript/SpatialET/06_resultados/ETP/mensuales/",nombres[i])
  writeRaster(mesi,ruta,format="GTiff",overwrite=TRUE)
 }}

# Las opciones para la función son: 
#Monthly_ETP_hg # para exportar las ETP mensuales del método de Hargreaves
#Monthly_ETP_twt # para exportar las ETP mensuales del método de Thorthwaite

Exportar_mensual(Monthly_ETP_twt) # en este ejemplo exportamos las de Hargreaves

#FIN DEL CÓDIGO


