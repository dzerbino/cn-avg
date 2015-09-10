# Need to make sure that the header line of the file does not begin with a '#' before reading it in. 

dat=read.table("simevntsp50.txt", header=TRUE)
edgedat=read.table("simedgesp00.txt", header=TRUE)

totaledges=edgedat$TP+edgedat$FP+edgedat$TN+edgedat$FN
trueedges=edgedat$TP+edgedat$FN+edgedat$TN
connectivity=trueedges/edgedat$truenodes

TP=dat$TP
FP=dat$FP
TN=dat$TN
FN=dat$FN

total=TP+FP+TN+FN
totaltrue=dat$TP+dat$FN
precision=TP/(FP+TP)
recall=TP/totaltrue
f1score=2*(precision*recall)/(precision+recall)
f1score[is.nan(f1score)]=0
	
colors=c(rgb(27,158,119, maxColorValue=255), 
rgb(217,95,2,maxColorValue=255),
rgb(117,112,179, maxColorValue=255))

ptcolors=dat$blocks
ptcolors[dat$blocks==100]=colors[1]
ptcolors[dat$blocks==150]=colors[2]
ptcolors[dat$blocks==200]=colors[3]

shapes=c(15,16,17)
ptshapes=dat$blocks
ptshapes[dat$blocks==100]=shapes[1]
ptshapes[dat$blocks==150]=shapes[2]
ptshapes[dat$blocks==200]=shapes[3]

#################----------------
# Function for getting the pvalues out of a lm
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}


#########################################################################
##  Plotting begins here. 

legendtext=c(100,150,200)
x=100
legendtext[1]=paste(x, " (", min(totaltrue[dat$blocks==x]), "-", max(totaltrue[dat$blocks==x]), ")", sep="")
x=150
legendtext[2]=paste(x, " (", min(totaltrue[dat$blocks==x]), "-", max(totaltrue[dat$blocks==x]), ")", sep="")
x=200
legendtext[3]=paste(x, " (", min(totaltrue[dat$blocks==x]), "-", max(totaltrue[dat$blocks==x]), ")", sep="")

myplotfunction <- function(connectivity, y, ylabel, ptcolors, ptshapes, legendtext, xaxtval="s", xlabel=""){
	mylabsize=1.5
	plot(connectivity, y, col=ptcolors, pch=ptshapes, xaxt=xaxtval, ylab=ylabel,cex.lab=mylabsize, xlab=xlabel)
	# The lines below can be substituted for the 2 lines following for an exponential decay fit rather than a linear fit.  The exponential decay is a better fit for histories of multi-clonal cell populations (not using the --int option in the cn-avg pipeline). 
	#points(connectivity[betterpreds], y[betterpreds], pch=1)
	#mylm=lm(y~exp(-connectivity))
	#xrange=seq(min(connectivity), max(connectivity), by=0.01)
	#ylm=mylm$coefficients[1]+mylm$coefficients[2]*exp(-xrange)
	#lines(xrange, ylm, col="gray")
	mylm=lm(y~connectivity)
	abline(mylm$coefficients, col="gray")
	mtext(paste("R-squared:",signif(summary(mylm)$r.squared, digits=3), ", p-value: ", signif(lmp(mylm), digits=3)), outer=FALSE, side=3, cex=0.7) 
	legend("topright", legend=legendtext, title="number of blocks (#events)", pch=shapes, col=colors, bty="n")
}

par(mfrow=c(3,1))
margins=c(2,4,0,2)+0.1
par(oma=c(4,4,2,2)+0.1)
par(mar=margins)
betterpreds=dat$truescores>dat$minscores
y=recall
myplotfunction(connectivity, y, ylabel="sensitivity", ptcolors, ptshapes, legendtext, xaxtval="n")
y=precision
myplotfunction(connectivity, y, ylabel="specificity", ptcolors, ptshapes, legendtext, xaxtval="n")
y=f1score
myplotfunction(connectivity, y, ylabel="f1score", ptcolors, ptshapes, legendtext, xaxtval="s", xlabel="connectivity")
mtext("connectivity", side=1, line=2, cex=1)
