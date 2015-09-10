#Note that before you read in the allevents.dat file, you need to make sure the header line does not start with a "#" or else it will be ignored. 

dat=read.table("allevents.dat", header=TRUE)
#########################################################
# HISTOGRAMS of plain ol' T vs F events. 

fpdat=dat[dat$true==0,]
tpdat=dat[dat$true==1 | dat$true==3,]  # if the 'true' value is 3, that means that the event doesn't perfectly match a true event, but is equivalent (is part of a linear combination) of a true event. 

mylwd=2
myylim=c(0,1000)
fph=hist(fpdat$Lscore, plot=FALSE)
tph=hist(tpdat$Lscore, plot=FALSE)
plot(tph$mids, tph$counts, type="S", col="black", lwd=mylwd, ylim=myylim, xlab="", ylab="", cex.lab=1.3)
lines(fph$mids, fph$counts, type="S", col="grey", lwd=mylwd, ylim=myylim)
legend("topleft", legend=c(
paste("True (", prettyNum(nrow(tpdat), big.mark=","), ")", sep=""), 
paste("False (", prettyNum(nrow(fpdat), big.mark=","),")", sep="")), 
lwd=mylwd, col=c("black", "grey"), bty="n")
mtext("likelihood score", side=1, line=2.5)
mtext("number of events", side=2, line=2.5)

tptot= sum(tpdat$Lscore>0.5)
fptot= sum(fpdat$Lscore>0.5)
mtext(paste("TP: ", tptot, ", FP: ", fptot, ", acc: ", round(tptot/(tptot+fptot), digits=4), sep=""), side=3)
