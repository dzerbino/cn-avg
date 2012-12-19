#!/usr/bin/env Rscript

### This code uses DNAcopy to find smooth CN for a matrix of raw CN's from
### several samples. Takes as input the file with the dChip output.
### Returns a file whose rows are segments and whose
### columns are: sample number, chromosome, beginning bp, end bp,
### number SNPs in segment, and mean copy numbers for
### segment. 


gpCBS<-function(infile, outfile, nPerm=10000, alpha=0.01, seed=12345678, libdir='', Rdir='') {
	if(libdir!='') {
		source(paste(libdir, "common.R", sep=''))
		
		if(!is.package.installed(libdir, "DNAcopy")) {
			install.package(libdir, "DNAcopy_1.12.0.zip", "DNAcopy_1.12.0.tgz", "DNAcopy_1.12.0.tar.gz")
		}
	}

    seed <- as.numeric(seed)
	set.seed(seed)
	nPerm <- as.numeric(nPerm)
	alpha <- as.numeric(alpha)
	splitJob <-F
	
	array.names <- get.array.names.cn(infile)

	cols <- length(array.names)
	
	iter.cols <- 1:cols
#	if(array.list.file!='') {
#		iter.cols <- get.column.indices(array.names, array.list.file)
#	}
	
	
	if(splitJob) {		
	   lsf <- suppressWarnings(library('Rlsf', logical.return=T))
		if(!lsf) {
		   splitJob <- F
		} else {
		   library(Rlsf)
	   }
	}
	jobs <- NULL
	cn <- read.cn(infile, iter.cols)
	if(!splitJob) {
	   library(DNAcopy)
	   result <- run.cbs(cn, -1, array.names, nPerm, alpha, NULL)
	   write.table(result, file=outfile, sep="\t", quote=F, row.names=F)
	} else {
      for(i in iter.cols) {
         ctrl <- lsf.ctrl()
         ctrl$R <- paste(R.home("bin"), "/R", sep='')
         savelist <- c("cn", "i", "array.names", "nPerm", "alpha", "libdir")
         ctrl$savelist <- savelist
         j <- lsf.submit2(func=run.cbs, cn, i, array.names, nPerm, alpha, libdir, ctrl=ctrl)
         jobs <- c(jobs, j)
         # ctrl$$tmpPath <- 
      } 
      interval <- 10
	   maxSleep <- 60 * 60
		wait.for.lsf.jobs(jobs, interval, maxSleep, F)
		for(job in jobs) {
		   result <- lsf.get.result(job)
		   write.table(result, file=outfile, sep="\t", quote=F, row.names=F, append=T)
		}
	}

}


run.cbs<-function(cn, column.index=-1, array.names, nPerm, alpha, libdir=NULL){
	if(!is.null(libdir)) {
		source(paste(libdir, "common.R", sep=''))
	   setLibPath(libdir)	
	   library(DNAcopy)
	}

	genomdat <- NULL
	chrom <- NULL
	maploc <- NULL
	
	if(column.index != -1) {
	   array.names <- array.names[[column.index]]	
	   genomdat <- cn$genomdat[,column.index]
	   chrom <- cn$chrom[,column.index]
	   maploc <- cn$maploc[,column.index]
	} else {
	   genomdat <- cn$genomdat
	   chrom <- cn$chrom
	   maploc <- cn$maploc
	}
	
	CNA.object <- CNA(genomdat, chrom, maploc, data.type = "logratio", sampleid = array.names)
	smoothed.CNA.object <- smooth.CNA(CNA.object)
	
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose = F, nperm=nPerm, alpha=alpha, undo.splits="sdundo")
	
	return(segment.smoothed.CNA.object$output)
}

main<-function(){
	args <- commandArgs(trailingOnly=TRUE)
	gpCBS(args[1], args[2], libdir='~/.local/lib/R/')	
}

main()
