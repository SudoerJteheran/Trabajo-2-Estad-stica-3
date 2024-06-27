library(forecast)
library(TSA)
library(car)
library(lmtest)
library(fANCOVA)

rm(list=ls(all=TRUE))
#----------------------------------------------------------------------------------------------------------------------------------------------------
# Cargando funciones de usuario necesarias desde repositorio en github
#----------------------------------------------------------------------------------------------------------------------------------------------------
#Cargando archivo con las funciones exp.crit.inf.resid() para calcular AIC y BIC, y amplitud.cobertura() para calcular amplitud media y coberturas de 
#los I.P
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")

#Cargando archivo con la funcion regexponencialv02() para hacer la regresion exponencial polinomial
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencialv02.R")

#Cargando archivo con la funcion Mipoly() para generar las variables de los polinomios
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Mipoly.R")

#Cargando archivo con la funcion Mytrigon() para generar las variables trigonometricas para frecuencias dadas
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Mytrigon.R")

#Cargando archivo con las funciones exp.crit.inf.resid() para calcular AIC y BIC, y amplitud.cobertura() para calcular amplitud media y coberturas de los I.P
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")

#Cargando archivo con funciones de usuario BP.LB.test() y pruebaDW1() 
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")

#Cargando archivo con la función de usuario SelectModel() para identificar modelos AR con criterios AIC y BIC
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SelectModel.R")



#----------------------------------------------------------------------------------------------------------------------------------------------------
#Lectura de los datos y graficos descriptivos
#----------------------------------------------------------------------------------------------------------------------------------------------------

setwd("C:/Users/admin/Desktop/Séptimo semestre/Estadística III/Trabajos/Trabajo II/Trabajo 2 Estadística III/")
Datos31=read.table("anex-EMMET-TotalNacional-oct2023-Fabricacion de Otros Productos Quimicos.csv",header=T,sep=";",skip=14,dec=",",colClasses=c(rep("NULL",4),"numeric",rep("NULL",6)))

#Utilizando los datos como dataframe, para que sea posible que se indiquen
#Los índices de forma "correcta" en la ACF

#Presentando en una sola ventana las serie
win.graph(width=6,height=6)
plot(Datos31$Produccion.nominal,xaxp=c(2001,2022,7),ylab="Yt",xlab="Tiempo", lwd=1, type='l')
#abline(v=seq(2001,2022,1), col="red", lty=2)
#Grafica de la serie transformada 
win.graph(width=6,height=6)
plot(log(Datos31$Produccion.nominal),ylab="log(Yt)",xlab="Tiempo",xaxp=c(2001,2022,7), lwd=1, type='l')





#Para ajustar la ACF sobre el log de la serie
win.graph(width=6,height=6)
acf(log(Datos31$Produccion.nominal),
    lag.max=36, 
    ci.type="ma",
    col=4,
    ci.col=2,
    main="",
    xlab = "Orden de Rezago",
    ylab = "Autocorrelación Estimada")

#Ahora sí, tornamos todo una serie de tiempo

Datos31=ts(Datos31,freq=12,start=c(2001,1))

