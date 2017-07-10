
data = read.table("../data/freq_2months.csv", header=TRUE, sep=',')
freqK = data$K
freqT = data$T
year = data$year
xat = seq(from=2012, to=2017, by=1)
yat = seq(from=0.0, to=1.0, by=0.5)

pdf("../figure/Figure1A.pdf", width=3.5, height=2.0)
par(xpd=T, mar=c(3,3,0.5,4.5))
plot(year, freqT, type='o', pch=16, col='red', xlab='', ylab='', bty="l", yaxt = 'n', xaxt = 'n', lwd=1.5)
lines(year, freqK, type='o', col=rgb(0.1,0.1,0.8), lwd=1.5, pch=16)
xis(side=1, at=xat, cex.axis=0.7, mgp=c(3, 0.6, 0), las=2, tck=-0.05)
axis(side=1, at=xat, cex.axis=0.7, mgp=c(3, 0.6, 0), las=2, tck=-0.05)
axis(side=2, at=yat, cex.axis=0.7, mgp=c(3, 0.6, 0), las=2, tck=-0.05)
mtext(side=1, text="year", line=1.8)
mtext(side=2, text="frequency", line=2.0)
legend("topright",legend=c('K160', 'T160'), col=c(rgb(0.1, 0.1, 0.8), 'red'), lty=c(1,1), pch=16, inset=c(-0.43, -0.0), bty="n", lwd=1.5, cex=0.8)
dev.off()


