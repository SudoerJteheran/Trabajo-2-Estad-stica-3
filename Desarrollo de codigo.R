setwd("C:/Users/admin/Desktop/Séptimo semestre/Estadística III/Trabajos/Trabajo II/Git Trabajo 2/Trabajo-2-Estad-stica-3")
library(forecast)
library(car)
library(TSA)
library(lmtest)
library(fANCOVA)

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

#Cargando archivo con la funcion regexponencialv02() para hacer la regresion exponencial polinomial
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencialv02.R")

#Cargando archivo con funciones de usuario BP.LB.test() y pruebaDW1() 
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")

#Cargando archivo con la funcion Mipoly() para generar las variables de los polinomios
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Mipoly.R")

source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SelectModel.R")

Datos31=read.table("anex-EMMET-TotalNacional-oct2023-Fabricacion de Otros Productos Quimicos.csv",header=T,sep=";",skip=14,dec=",",colClasses=c(rep("NULL",4),"numeric",rep("NULL",6)))
Datos31=ts(Datos31,freq=12,start=c(2001,1))

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ACF de los datos
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

acf(as.numeric(log(Datos31)),lag.max=36,ci.type="ma",
    col='black',ci.col=2,lwd=2,main='')

#----------------------------------------------------------------------------------------------------------------------------------------------------
#Definiendo variables necesarias para el ajuste y pronostico
#----------------------------------------------------------------------------------------------------------------------------------------------------
#definir variables para el ajuste

n=length(Datos31)-12

#indice de tiempo en los periodos de ajuste 
t=1:n
yt=ts(Datos31[t],freq=12,start=c(2001,1)) #valores de la serie en muestra de ajuste

#variables para el ajuste, en el polinomio de grado 4
poli=Mipoly(tiempo=t,grado=4)

mes=seasonaldummy(yt) #Matriz con las 11 primeras variables Indicadoras 


#Matriz de diseno con las potencias de t en polinomio de grado 3 e indicadoras para la estacionalidad, en el ajuste
X1=data.frame(poli,mes)


#head(X1,n=3)  #para ver las primeras tres filas de X1


#Valores de las variables en la validacion cruzada
tnuevo=(n+1):length(Datos31) #indice de tiempo en los pronosticos
ytnuevo=ts(Datos31[tnuevo],freq=12,start=c(2022,11)) #valores reales de la serie en la muestra de validacion cruzada

#variables para el pronostico, en el polinomio de grado 4
polinuevo=Mipoly(tiempo=tnuevo,grado=4)

mesnuevo=seasonaldummy(yt,h=12) #matriz con las 11 indicadoras para los tiempos de pronostico

#Matriz de diseno con las potencias de t en polinomio de grado 3 e indicadoras para la estacionalidad, en los pronosticos
X1nuevo=data.frame(polinuevo,mesnuevo)


#---------------------------------------------------------------------------------------------------------------------------------------------
#Ajustes y pronosticos Modelo 1: log-polinomial grado 4 estacional con indicadoras; usa matriz X1 en el ajuste y X1nuevo en el pronostico
#---------------------------------------------------------------------------------------------------------------------------------------------

mod1=lm(log(yt)~.,data=X1)
summary(mod1)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Graficos de supuestos
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

layout(rbind(c(1,2)))
plot.ts(residuals(mod1),ylab='Residual mod1',xlab='Tiempo')
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),col=2)
plot(fitted(mod1),residuals(mod1))
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),col=2)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prueba de Lungbox
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BP.LB.test(residuals(mod1),maxlag=36,type="Ljung")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prueba de Durbin Watson
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pruebaDW1(mod1)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Grafica del ACF y la PACF
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

layout(rbind(c(1,2)))
acf(residuals(mod1),ci.type="ma",lag.max=36)
pacf(residuals(mod1),lag.max=36)

## Serie de tiempo de los residuales 

 
# Residuales del modelo en serie de tiempo 
resi <- ts(mod1$residuals,freq=12,start=c(2001,1))

eacf(resi,24,24)

auto.arima(resi,ic='aic')
auto.arima(resi,ic='bic')

auto.arima(residuals(mod1), ic= "aic")
auto.arima(residuals(mod1), ic= "bic")

SelectModel(resi,24,Criterion='AIC',ARModel='AR')
SelectModel(resi,24,Criterion='BIC',ARModel='AR')

win.graph()
plot(armasubsets(resi, ar.method = "ml", nar = 18, nma = 18))


#---------------------------------------
#Ajuste para el modelo 4
#---------------------------------------
#Creamos el vector donde vamos a fijar los parámetros a utilizar

fijados=c(NA, NA, rep(0,4), NA, 0, NA, rep(0,5), NA, #Parámetros para parte AR
        rep(0,6),NA,rep(0,8), NA, # -- -- -- MA
        rep(NA,16)) # -- -- tendencia y estacionalidad

modelo4 <- Arima(log(yt),
                 order = c(15, 0, 16),
                 xreg = as.matrix(X1),
                 fixed = fijados,
                 method = "ML")

#Calcule grados de libertad del MSE del modelo4
k4=length(coef(modelo4)[coef(modelo4)!=0]);k4 #número de parámetros del modelo4
dfmodelo4=n-k4

#Valor de la varianza para las innovaciones
modelo4$sigma2

#Construya tabla de parámetros estimados con estadísticos T0 y valores P para cada parámetro
coeftest(modelo4,df=dfmodelo4) 

#Modelo estimado para la escala original
ythat4=exp(modelo4$fitted)*exp(modelo4$sigma2/2) #este objeto ya 
                                                 #queda con las 
                                                #fechas de la serie

#Gráfico de la serie y su ajuste
win.graph()
plot(Datos31)
lines(ythat4,col=2)
legend("topleft",legend=c("datos","ajuste modelo4"),lty=1,col=1:2)

#Gráficos de residuales
win.graph()
plot(residuals(modelo4))
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),lty=2)
win.graph()
plot(as.numeric(modelo4$fitted),residuals(modelo4))
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),lty=2)

#Calculo de AIC y BIC usando Exp(Cn*(p))
Res.orig4=yt-ythat4 #Cálculo de pseudo residuales
Criteriosmodelo4=exp.crit.inf.resid(residuales=Res.orig4,n.par=k4)
Criteriosmodelo4

#------------------------------------------
#Validación de supuestos sobre innovaciones
#------------------------------------------
#ACF sobre residuales de ajuste en el modelo4. Use valor para m el que se indica en la guía del trabajo
win.graph()
m = 36
acf(as.numeric(residuals(modelo4)),ci.type="ma",lag.max=m,ci.col=2)
#PACF sobre residuales de ajuste en el modelo4. Use valor para m el que se indica en la guía del trabajo
win.graph()
pacf(as.numeric(residuals(modelo4)),lag.max=m,ci.col=2)
BP.LB.test(residuals(modelo4),maxlag=m,type="Ljung") #test Ljung-Box use máximo m igual al de ACF y PACF
#Normalidad sobre residuales de ajuste en el modelo4. Sólo si no se detectó correlaciones
shapiro.test(residuals(modelo4))
win.graph()
qqnorm(residuals(modelo4))
legend("topleft",legend="modelo4")
qqline(residuals(modelo4),col=2)

#---------------------------------------
# Pronósticos para la validación cruzada
#---------------------------------------
#Cálculo del pronóstico con I.P del 95%, en escala original
predmodelo4=exp(as.data.frame(forecast(modelo4,xreg=as.matrix(X1nuevo),
                                       level=95)))*exp(modelo4$sigma2/2)
predmodelo4=ts(predmodelo4,freq=12,start=start(ytnuevo))
predmodelo4
ytpronmodelo4=predmodelo4[,1] #Tomando el pronóstico puntual
#Medidas precisión pronósticos
accuracy(ytpronmodelo4,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmodelo4[,2],LSP=predmodelo4[,3])

#Guardando en excels las estimaciones

#Exportando al directorio de trabajo archivo excel con tabla de par´ametros estimados
write.csv2(coeftest(modelo4,df=dfmodelo4),file="tablamodelo4trabajo2.csv",row.names = TRUE)
#Exportando al directorio de trabajo archivo excel con pronosticos
write.csv2(predmodelo4,file="pronosticosmodelo4trabajo2.csv",row.names = paste(trunc(time(predmodelo4)),cycle(predmodelo4),sep="/"))
