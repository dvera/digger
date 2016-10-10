# function to take a feature counts table and a list of genes, make a log counts table (rpkm) and plot a heatmap for those genes
# depends on conifur and travis library and gplots

rpkmHeatmap <- function(countsTable, genefiles, samplenames , marg=c(10,10), rpkmout=FALSE, pdfOut=NULL,threads=getOption("threads",1L) ) {
	library(edgeR)

	cnts=tsvRead(countsTable,col_names=TRUE)		        # read in feature counts table
	numsamples <- ncol(cnts)-6				                # first 6 cols are extra info
	if(missing(samplenames))
		samplenames <- colnames(cnts)[7:ncol(cnts)]	        # default, samplenames are file paths to samples
	else
		stopifnot(length(samplenames)==numsamples)	        # exit if samplenames vector not the correct length
	colnames(cnts)[7:ncol(cnts)] <- samplenames		        # set column names of table

	rownames(cnts)=toupper(cnts[,1])			            # make all genes uppercase for searching
	genelengths=cnts[,6]				        	        # extract gene lengths from table for rpkm calc
	names(genelengths)=row.names(cnts)
	cnts=cnts[-(1:6)]					                    # remove geneid,chr,start,end,length cols from cnts table, gives simple counts table
	cnts=data.matrix(cnts)
	rpkms=rpkm.default(cnts,genelengths)			        # make rpkm counts table

	if(rpkmout) {
		fo=paste0( removeext( basename(countsTable) ), "_rpkm.tsv")
		tsvWrite(as.data.frame(rpkms),fo,col_names=T,row_names=T)	# print rpkm table
	}

	e=which(filelines(genefiles)==0)			# check for empty files
	if(length(e) != 0){
		cat("Empty files not considered:\n")
		for( i in 1:length(e) )
			print(genefiles[e[i]])
		genefiles=genefiles[-e]				                # remove empty files from list
	}
	genes=tsvRead(genefiles)				                # read in list of genes to use in map
	genes=toupper( unique( unlist(genes) ) )	    	    # unlist and remove gene duplicates

	w=which( rownames(rpkms) %in% genes )			        # subset rpkm table to just genes from genelist
	rpkmsSub=rpkms[w,]

	if(!is.null(pdfOut)) {
		pdf(pdfOut)
	}

	heatmap.2(log2(1+rpkmsSub),trace="none",margins=marg)	# make heatMap of log2 rpkms genes subset

	if(!is.null(pdfOut)) {
		dev.off()
	}

}
