# Manuel Campagnolo
# 9 de fevereiro de 2023

library(sp)
library(spatialreg)
library(spdep)
library(gstat)

# Analyze effect of interaction
# Compare response witout spatial autocorrelation (Y0) and response when spatial autocorrelation occurs (YA)

# lista para guardar resultados:
RESULTS=list()

# ciclo para passar pelos sites e pelo número de vizinhos
for (SITE in LETTERS[1:4]) # "D") #
for (K in 8:8) # 
{
print(paste(K,SITE))

results=list()


DF=read.csv2("Dados_SpatialCorkOak.csv",header=TRUE,sep=";",dec = ",")
head(DF)
dim(DF)
YA=DF$du #du_annual_growth #
mort=DF$Morta
x=DF$coordsX
y=DF$coordsY
mais_velha=DF$mais_velha
sites=DF$Site
# "probmort"
vars=c("Cea_1m_1px","Cea_0.5m_1px","Altimetria_1px","slope_1px","TPI_tree_1px","TWI_tree_1px","probmort")
results$vars=vars

# select data depending on cond
cond=(sites==SITE & mort==0 & mais_velha==0 & YA>=0)
n=sum(cond)
YA=YA[cond] # has columns X1,..., and Y
x=x[cond] # has columns X1,..., and Y
y=y[cond] # has columns X1,..., and Y
DF=DF[cond,]
dfxy=cbind(data.frame(x,y),DF[,vars])
df=dfxy
df$x=NULL
df$y=NULL
coordinates(dfxy)=c("x","y") # dfxy is a sp object

# scale response variable
#df$Y=scale(YA)
head(df)

# nlist
# to avoid empty neighborhoods
set.ZeroPolicyOption(TRUE)
get.ZeroPolicyOption()
#nlist=spdep::dnearneigh(dfxy,d1=0,d2=d2) # problems with empty neighborhoods
nlist=spdep::knn2nb(knearneigh(dfxy,k=K))

############################################################# analize data set
if (FALSE)
{
  # estimate rho (lento)
  lagsarlm(Y~.,data=df,listw=nb2listw(nlist,style="W",zero.policy = TRUE))$rho
  # estimate lambda
  errorsarlm(Y~.,data=df,listw=nb2listw(nlist,style="W",zero.policy = TRUE))$lambda
}

# plot values
plot(dfxy$x,dfxy$y,pch='.') # locations (x,y)
text(dfxy$x,dfxy$y,label=round(YA,2),cex=0.5) # locations (x,y)

# distance based (1 connection is 1 unit) Moran's I
plot(sp.correlogram(nlist,as.vector(YA),order=10,method="I",style="W"))

# variogram
dfxy$Y=YA
myvgm=variogram(YA~1, dfxy)
plot(myvgm)
mymodel <- gstat::fit.variogram(myvgm,vgm("Sph"),fit.kappa=TRUE)
plot(myvgm, model=mymodel)

############################################################## fit models

# use a model that doesn't "model" autocorrelation (regular linear model)
df$Y=YA
model=lm(Y~.,data=df)
summary(model)
Ypred=fitted(model)
plot(Ypred~YA)
results$lrr2=paste('LR R2=',1-sum((Ypred-YA)^2)/(var(YA)*(n-1)))

# lagmodel 
df$Y=YA
model=lagsarlm(Y~.,data=df,listw=nb2listw(nlist,style="W",zero.policy = TRUE))
summary(model)
Ypred=fitted(model)
plot(Ypred~YA)
results$lagr2=paste('lagmodel R2=',1-sum((Ypred-YA)^2)/(var(YA)*(n-1)))

# spatial error model 
df$Y=YA
model=errorsarlm(Y~.,data=df,listw=nb2listw(nlist,style="W",zero.policy = TRUE))
summary(model)
Ypred=fitted(model)
plot(Ypred~YA)
results$errr2=paste('spatial error model R2=',1-sum((Ypred-YA)^2)/(var(YA)*(n-1)))

# both
df$Y=YA
model=sacsarlm(Y~.,data=df,listw=nb2listw(nlist, style="W",zero.policy = TRUE))
summary(model)
Ypred=fitted(model)
plot(Ypred~YA)
results$errr2=paste('spatial error model R2=',1-sum((Ypred-YA)^2)/(var(YA)*(n-1)))

# compare models with spdep::lm.LMtests
# Lagrange multiplier diagnostics for spatial dependence (pg 435)
model=lm(YA~.,data=df)
res=spdep::lm.LMtests(model,nb2listw(nlist,style="W",zero.policy = TRUE),test="all")
for (mod in c("LMerr",  "LMlag",  "RLMerr", "RLMlag" ,"SARMA")) results[[mod]]=res[[mod]]$p.value

# atualizar lista de resultados
RESULTS[[paste(SITE,K)]]=results
}

