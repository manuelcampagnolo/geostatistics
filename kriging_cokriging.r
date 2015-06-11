###############################################################################
#
# interpolação espacial de observações pontuais: Introdução 
# Kriging, Co-kriging
#

temps<-readOGR(dsn=getwd(),layer="stations_iberia_wgs84",encoding="ISO8859-1")
temps@proj4string
iberia<-readOGR(dsn=getwd(),layer="IberianPeninsula",encoding="ISO8859-1")
iberia@proj4string
tempib<-spTransform(temps,iberia@proj4string)
colnames(tempib@data)

# há observações que não são "continentais"
plot(tempib)
plot(iberia,add=TRUE)
# vamos seleccionar apenas as "continentais"
its<-gIntersects(tempib,iberia,byid=TRUE) # devolve matriz 2*56
# calcular vector lógico dos pontos que estão em 1 dos 2 polígonos "continentais" (Portugal, Espanha)
ctl<-as.logical(apply(its,2,max)) 
# seleccionar de tempib os pontos "continentais":
tempib<-tempib[ctl,] # tem 45 pontos
xy<-coordinates(tempib)
plot(iberia)
text(xy[,1],xy[,2],tempib@data$name,cex=.6)

# preparar dados
# criar data.frame com x, y, e variável resposta
xyt<-data.frame(coordinates(tempib))
colnames(xyt)<-c("x","y")

#Create sampling grid
x.range <- as.integer(range(xyt[, "x"]))
y.range <- as.integer(range(xyt[, "y"]))
grd <- expand.grid(x = seq(x.range[1],x.range[2], length.out=20), 
                   y = seq(y.range[1],y.range[2], length.out=20))
grd.sp <- SpatialPoints(grd, proj4string=iberia@proj4string)

#################################################################################
#

library(gstat) # krige, idw
# ver http://www.css.cornell.edu/faculty/dgr2/teach/R/R_ck.pdf

# organizar dados: guardar só os que vão ser necessários:
#
elev<-as.numeric(as.character(tempib@data$grelev))
jan<-as.numeric(as.character(tempib@data$jan))

# criar SpatialPoints com essa informação
co<-data.frame(x=xyt[,"x"], y=xyt[,"y"], jan=jan,elev=elev) # 45*2 
# eliminar NAs
co<-co[!is.na(co$elev),]
# o comando abaixo converte a data.frame em SpatialPointDataFrame
coordinates(co) <- ~ x+y
# associar um CRS a codf
proj4string(co) <- iberia@proj4string


##################################################################################
#
# modelar uma só variável: temperatura em janeiro
#
# A - simples regressão sobre x e y 
k<-krige(jan~x+y,locations=co,newdata=grd.sp)

# seleccionar os pontos de predição no interior da península
its<-gIntersects(k,iberia,byid=TRUE) # devolve matriz 2*56
# calcular vector lógico dos pontos que estão em 1 dos 2 polígonos "continentais" (Portugal, Espanha)
itspi<-as.logical(apply(its,2,max)) 
# seleccionar de tempib os pontos "continentais":
k<-k[itspi,] # 

xyk<-coordinates(k)
# imagem
if (export)  png(paste(aulas,"mapa_regressao_multipla_sobre_x_y.png",sep="\\"), width=800, height=800, res=120)
plot(iberia)
text(xyt[,"x"],xyt[,"y"],round(co@data[,"jan"]),col="blue",cex=.5)
text(xyk[,1],xyk[,2],round(k@data[,"var1.pred"]),cex=.5)
if (export) graphics.off()

# comparar com regressão linear para ver que é igual
mod<-lm(jan~x+y, data=cbind(coordinates(co),co@data))
z<-predict(mod, newdata=as.data.frame(coordinates(grd.sp)))
text(coordinates(grd.sp)[,1],coordinates(grd.sp)[,2],round(z),cex=.5,col="red")

#
# B - usar variograma para a temperatura de Janeiro
# 

# visualizar variograma
# variograma, com pressuposto de anisotropia, limitado a 1/3 da distância maior e 15 "lags"
plot(variogram( jan~x+y,locations=co))
# variograma, com pressuposto de anisotropia, limitado a 1/3 da distância maior e com lag de 20km
plot(variogram( jan~x+y,locations=co,width=20000))

# modelos disponíveis para o variograma
vgm() #gstat
?vgm

# ajustar parâmetros do variograma depois de escolher o modelo
if (export)  png(paste(aulas,"variaograma_temperatura.png",sep="\\"), width=800, height=800, res=120)
v<-variogram(jan~x+y,locations=co)
plot(v) # para ajudar a escolher ponto inicial da pesquisa com vgm
vf<-fit.variogram(v,vgm(psill=100000,model="Exp",range=300000,nugget=0))
plot(v,vf)
if (export) graphics.off()

# nova predição usando o modelo para o variograma
kt<-krige(jan~x+y,locations=co,newdata=grd.sp,model=vf)
# seleccionar os pontos de predição no interior da península
its<-gIntersects(kt,iberia,byid=TRUE) # devolve matriz 2*56
# calcular vector lógico dos pontos que estão em 1 dos 2 polígonos "continentais" (Portugal, Espanha)
itspi<-as.logical(apply(its,2,max)) 
# seleccionar de tempib os pontos "continentais":
kt<-kt[itspi,] # tem 59 pontos dos 100 iniciais
xykt<-coordinates(kt)

# imagem
if (export)  png(paste(aulas,"mapa_sem_co_regionalizacao.png",sep="\\"), width=800, height=800, res=120)
plot(iberia)
text(xyt[,"x"],xyt[,"y"],round(co@data[,"jan"]),col="blue",cex=.5)
text(xykt[,1],xykt[,2],round(kt@data[,"var1.pred"]),cex=.5)
if (export) graphics.off()

# imagem com "bubbles"
bubble(kt,"var1.pred",maxsize=1.5)

# validação de resultados com validação cruzada
if (export)  png(paste(aulas,"bubbles_residuos_sem_co_regionalizacao.png",sep="\\"), width=800, height=800, res=120)
cv<-krige.cv(jan~x+y,locations=co,model=vf,nfold=10) # cv ainda é SpatialPointsDataFrame 
bubble(cv,"residual")
if (export) graphics.off()

##################################################################################
#
# usar uma co-variável: elevação
#

# relação entre elevação (variável grelev) e temperatura de Janeiro:

plot(jan~elev)

# ajustar variograma para variável temperatura de Janeiro (como atrás)
plot(variogram(jan ~x+y, locations=co))
vft<-fit.variogram(v,vgm(psill=100000,model="Exp",range=400000,nugget=0))
plot(v,vft)

# observar variagrama para variável elevação e ajustar variograma
plot(variogram(elev ~x+y, locations=co))
vfe<-fit.variogram(v,vgm(psill=200000,model="Exp",range=400000,nugget=0))
plot(v,vfe)

# criar objecto gstat 
g <- gstat(NULL, id = "jan", form = jan ~ x+y, data=co)
g <- gstat(g, id = "elev", form = elev ~ x+y, data=co)

#  variogramas cruzados (há 3 "frames")
v.cross <- variogram(g)
plot(v.cross, pl=T)

# ajustar um modelo linear de "co-regionalização"
# usar o mesmo "range" e estrutura para todas as 3 "frames"
g <- gstat(g, id = "jan", model = vft, fill.all=T)
# ajustar um modelo linear de co-regionalização
g <- fit.lmc(v.cross, g)
if (export)  png(paste(aulas,"variograma_v_cross.png",sep="\\"), width=800, height=800, res=120)
plot(variogram(g), model=g$model)
if (export) graphics.off()

# fazer predição (c de co-regionalização)
kc <- predict.gstat(g, newdata=grd.sp)

# mostrar resultado
its<-gIntersects(kc,iberia,byid=TRUE) # devolve matriz 
# calcular vector lógico dos pontos que estão em 1 dos 2 polígonos "continentais" (Portugal, Espanha)
itspi<-as.logical(apply(its,2,max)) 
# seleccionar de tempib os pontos "continentais":
kc<-kc[itspi,] # 
xykc<-coordinates(kc)
# imagem
if (export)  png(paste(aulas,"mapa_com_co_regionalizacao.png",sep="\\"), width=800, height=800, res=120)
plot(iberia)
text(xyt[,"x"],xyt[,"y"],round(co@data[,"jan"]),col="blue",cex=.5)
text(xykc[,1],xykc[,2],round(kc@data[,"jan.pred"]),cex=.5)
if (export) graphics.off()


# validação de resultados com validação cruzada
cvc <- gstat.cv(g,nfold=10)
if (export)  png(paste(aulas,"bubbles_residuos_com_co_regionalizacao.png",sep="\\"), width=800, height=800, res=120)
bubble(cvc,"residual")
if (export) graphics.off()

# comparar valores observados e previstos pelos modelos kt e kc
if (export)  png(paste(aulas,"comparacao_com_sem_co_regionalizacao.png",sep="\\"), width=800, height=800, res=120)
par(mfrow=c(1,1),mar=c(3,3,.2,1))
auxt<-cv@data
auxc<-cvc@data
plot(auxt$observed,auxt$var1.pred,xlab="temperatura observada",ylab="temperatura estimada",pch=16,asp=1)
points(auxc$observed,auxc$jan.pred,col="red",pch=16)
abline(0,1)
text(min(auxt$observed), max(auxt$var1.pred), "sem co-regionalização",pos=4)
text(min(auxt$observed), quantile(auxt$var1.pred,0.8), "com co-regionalização", col="red",pos=4)
if (export) graphics.off()
