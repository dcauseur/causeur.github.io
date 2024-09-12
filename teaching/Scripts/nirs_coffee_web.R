
# Import data

dta <- read.table("nirscafe.txt",sep=";",
                   header=TRUE,
                   stringsAsFactors=TRUE,
                   fileEncoding = "latin1")
dta$Localisation <- factor(dta$Localisation)

x <- as.matrix(dta[,-(1:8)])   # NIRS
snv_x = t(scale(t(x)))         # SNV-transformed NIRS 

wl = seq(402,2500,2)           # Wavelengths
matplot(wl,t(snv_x),type="l",lwd=2,col="orange",lty=1,
        bty="l",xlab=expression(Wavelength~(nm)),
        ylab="SNV-transformed NIRS",main="NIRS data of coffee samples",
        cex.main=1.25,cex.lab=1.25,cex.axis=1.25)
grid()

# Correlation between NIRS and biochemical profile

cor_xy <- cor(snv_x,dta[,c(3,5)])

matplot(wl,cor_xy,type="l",lty=1,
   col=c("blue","orange"),
   lwd=2,xlab=expression(Wavelength~(nm)),
   ylab="Corrrelation NIRS-Biochemical profile",
   main="Correlation between NIRS and biochemical profile",
   cex.main=1.25,cex.lab=1.25,cex.axis=1.25)
grid()
abline(h=0,lwd=2,col="darkgray")
legend("bottomleft",bty="n",col=c("blue","orange"),
   lty=1,lwd=2,legend=c("Fat","DM"),cex=1.25)

# Comparison of mean NIRS regarding their origins

m_loc <- apply(snv_x,2,function(x,g) tapply(x,INDEX=g,FUN=mean,na.rm=TRUE),g=dta$Localisation)
m_loc <- m_loc[c(1,6),]

matplot(wl,t(m_loc),type="l",lty=1,
   col=c("blue","orange"),
   lwd=2,xlab=expression(Wavelength~(nm)),
   ylab="Mean NIRS",
   main="Mean NIRS in all sites",
   cex.main=1.25,cex.lab=1.25,cex.axis=1.25)
grid()
abline(h=0,lwd=2,col="darkgray")
legend("topleft",bty="n",col=c("blue","orange"),
   lty=1,lwd=2,legend=c(1,6),cex=1.25,
   title="Localisation")

matplot(wl[(wl>=402)&(wl<=700)],t(m_loc[,(wl>=402)&(wl<=700)]),type="l",lty=1,
   col=c("blue","orange"),
   lwd=2,xlab=expression(Wavelength~(nm)),
   ylab="Mean NIRS",
   main="Mean NIRS in all sites",
   cex.main=1.25,cex.lab=1.25,cex.axis=1.25)
grid()
abline(h=0,lwd=2,col="darkgray")
legend("topright",bty="n",col=c("blue","orange"),
   lty=1,lwd=2,legend=c(1,6),cex=1.25,
   title="Localisation")








