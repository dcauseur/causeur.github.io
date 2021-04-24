

setwd("C:/Documents and Settings/David Causeur/Mes documents/HOME/david/ens/Master Tase/TD")
library(gstat)
library(akima)
library(lattice)

### Importation des données

bt = read.table("données/barre_t.txt",sep="\t",dec=",",header=TRUE)
names(bt)[4] = "y"
coordinates(bt) = ~x+y

### Cartes

plot(coordinates(bt),bty="n",pch=16,col="blue",xlab="x",ylab="y",cex.lab=1.25)
title("Points de mesure \n Parcelle Barre Thomas")

bubble(bt,"OC",col="blue",key.entries=seq(10,20,1),
   main="Teneur en Carbone Organique (g/kg) \n Parcelle barre Thomas",
   maxsize=2)

### Interpolation linéaire

x = coordinates(bt)[,1]
y = coordinates(bt)[,2]

grid = interp(x,y,bt$OC,xo=seq(min(x), max(x), length=50), yo=seq(min(y), max(y), length=50))
egrid = expand.grid(grid$x,grid$y)

par(bg = "slategray")
persp(grid, phi = 30 , theta = 45 , col = "green3", 
scale = TRUE,shade = 0.75, border = NA, box = TRUE,xlab="X",ylab="Y",
zlab="Carbone Organique (g/kg)") 
par(bg = "white")

tcol <- terrain.colors(12)
opar <- par(pty = "s", bg = "lightcyan")
contour(grid, lty = "solid", 
add = FALSE,vfont = c("sans serif", "plain"),col=tcol[2])
title("Teneur en Carbone Organique (g/kg) \n Parcelle Barre-Thomas", font = 4)
     
egrid = expand.grid(grid$x,grid$y)
bt.grid = data.frame(x=egrid$Var1,y=egrid$Var2,OC=as.vector(grid$z))
contourplot(OC~x+y, bt.grid, aspect = "iso",region=TRUE,
main="Teneur en Carbone Organique (g/kg) \n Parcelle Barre-Thomas")

### Interpolation IDW

x = coordinates(bt)[,1]
y = coordinates(bt)[,2]

grid = interp(x,y,bt$OC,xo=seq(min(x), max(x), length=10), yo=seq(min(y), max(y), length=10))
egrid = expand.grid(grid$x,grid$y)
names(egrid) = c("x","y")
coordinates(egrid) =~x+y

g = gstat(id="OC",formula=OC~1,data=bt)
oc.idw = predict(g,newdata=egrid)
contourplot(OC~x+y, data=data.frame(coordinates(oc.idw),OC=oc.idw$OC.pred), aspect = "iso",region=TRUE,
main="Teneur en Carbone Organique (g/kg) \n Parcelle Barre-Thomas")

### Loess

x = coordinates(bt)[,1]
y = coordinates(bt)[,2]

alpha=0.3
rdt.lo = loess(OC~x*y,data=data.frame(x=x,y=y,OC=bt$OC),span=alpha)

grid = interp(x,y,bt$OC,xo=seq(min(x), max(x), length=50), yo=seq(min(y), max(y), length=50))
egrid = expand.grid(grid$x,grid$y)
names(egrid) = c("x","y")

fit = predict(rdt.lo,newdata=egrid)

contourplot(OC~x+y, data=data.frame(egrid,OC=as.vector(fit)), aspect = "iso",region=TRUE,
main="Teneur en Carbone Organique (g/kg) \n Parcelle Barre-Thomas")

seqalpha = seq(0.15,0.9,0.01)
aic = rep(0,length(seqalpha))
n = nrow(bt)

for (i in 1:length(seqalpha)) {
rdt.lo = loess(OC~x*y,data=data.frame(x=x,y=y,OC=bt$OC),span=seqalpha[i])
aic[i] = n*log(sum(rdt.lo$residuals^2)/n)+2*rdt.lo$enp
print(rdt.lo$enp)
}


plot(seqalpha,aic,xlab=expression(alpha),bty="l",lwd=2,col="blue",ylab="AIC",
type="l",cex.lab=1.25,main=expression(Choix~de~alpha))
text(0.6,120,expression(AIC(alpha)==n*log(SCER(alpha)/n)+2*p(alpha)),cex=1.25)

alpha=seqalpha[which.min(aic)]
rdt.lo = loess(OC~x*y,data=data.frame(x=x,y=y,OC=bt$OC),span=alpha)

grid = interp(x,y,bt$OC,xo=seq(min(x), max(x), length=50), yo=seq(min(y), max(y), length=50))
egrid = expand.grid(grid$x,grid$y)
names(egrid) = c("x","y")

fit = predict(rdt.lo,newdata=egrid)

contourplot(OC~x+y, data=data.frame(egrid,OC=as.vector(fit)), aspect = "iso",region=TRUE,
main="Teneur en Carbone Organique (g/kg) \n Parcelle Barre-Thomas")

### Krigeage

g=gstat(id="OC", formula=OC~1, data=bt, beta=mean(bt$OC))
variog=variogram(g,cutoff=150,alpha=c(0,90))
plot(variog,pch=16,cex=1.25,col="blue",cex.lab=1.5,xlab="h",main="Variogrammes directionnels")

variog=variogram(g,cutoff=150)
plot(variog,pch=16,cex=1.25,col="blue",cex.lab=1.5,xlab="h",main="Variogramme")
g.vgm=fit.variogram(variog, model=vgm(nugget=0, "Sph", range=100, psill=6), fit.sills=c(FALSE,TRUE))

plot(variog, model=g.vgm, xlab="distance (m)", ylab="semi-variance ((g/kg)^2)",ylim=c(0,8),pch=16,col="blue",
main="Ajustement d'un variogramme sphérique \n Nugget = 0, Psill = 7.22, Range = 149.24",cex=1.25)

x = coordinates(bt)[,1]
y = coordinates(bt)[,2]

grid = interp(x,y,bt$OC,xo=seq(min(x), max(x), length=30), yo=seq(min(y), max(y), length=30))
egrid = expand.grid(grid$x,grid$y)
names(egrid) = c("x","y")
coordinates(egrid) =~x+y
g=gstat(id="OC", formula=OC~1, data=bt, beta=mean(bt$OC),model=g.vgm)
oc.krige = predict(g,newdata=egrid)
contourplot(OC~x+y, data=data.frame(coordinates(oc.krige),OC=oc.krige$OC.pred), aspect = "iso",region=TRUE,
main="Teneur en Carbone Organique (g/kg) \n Parcelle Barre-Thomas")

### Co-krigeage

bt = read.table("données/barre_t.txt",sep="\t",dec=",",header=TRUE)
names(bt)[4] = "y"
bt = na.omit(data.frame(bt,resOC=residuals(rdt.lo)))
coordinates(bt) = ~x+y

x = coordinates(bt)[,1]
y = coordinates(bt)[,2] 
grid = interp(x,y,bt$OC,xo=seq(min(x), max(x), length=30), yo=seq(min(y), max(y), length=30))
egrid = expand.grid(grid$x,grid$y)
names(egrid) = c("x","y")
coordinates(egrid) =~x+y


round(cor(data.frame(bt)[,c("OC","CEC.metson")],use="complete"),2)

g = gstat(id = "OC",formula = OC~1, data = na.omit(bt), beta=mean(bt$OC))
g = gstat(g, id = "CEC",formula = CEC.metson~1, data = na.omit(bt), beta=mean(bt$CEC.metson,na.rm=TRUE))

variog=variogram(g,cutoff=150)

plot(variog, xlab="distance (m)", ylab="semi-variance",pch=16,col="blue",
main="Ajustement d'un covariogramme sphérique",cex=1.25)

savePlot("pic/covario2dajust1",type="pdf")

g = fit.lmc(variog, g, model=vgm(psill=0.15, range=125, nugget=0, model="Sph"),
fit.sills=c(FALSE,TRUE,FALSE,TRUE,FALSE,TRUE))
plot(variog, model=g,xlab="distance (m)", ylab="semi-variance",pch=16,col="blue",
main="Ajustement d'un covariogramme sphérique",cex=1.25)

oc.cokrige = predict(g,newdata=egrid)

contourplot(OC~x+y, data=data.frame(coordinates(oc.cokrige),OC=oc.cokrige$OC.pred), aspect = "iso",region=TRUE,
main="Teneur en Carbone Organique (g/kg) \n Parcelle Barre-Thomas")

contourplot(CEC~x+y, data=data.frame(coordinates(oc.cokrige),CEC=oc.cokrige$CEC.pred), aspect = "iso",region=TRUE,
main="CEC (cmol/kg) \n Parcelle Barre-Thomas")




