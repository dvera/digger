edger <- function( counts , grouping=NULL , samples=NULL , tabletype="featureCounts", dispersion=NULL, threads=getOption("threads",1L)  ){

  library(edgeR)
  library(readr)
  # library(biomaRt)
  # library(GO.db)
  library(plyr)
  # library(pathview)
  # library(KEGGREST)
  # library(gage)
  # library(gageData)
  # data(go.sets.hs)
  # data(go.subs.hs)
  # data(kegg.gs)


  if(tabletype=="cuffdiff"){
    rawcnts=read.tsv(counts,header=T)
    numsamples=(ncol(rawcnts)-1)/5
    cnts=rawcnts[,(2+(0:(numsamples-1))*5)]
    row.names(cnts)<-rawcnts[,1]
    rm(rawcnts)
  } else if ( tabletype=="counts" ){
    cnts=read.tsv(counts,header=TRUE,row.names=1)
  } else if ( tabletype=="featureCounts" ){
    cnts=read.tsv(counts,header=TRUE,row.names=1)
    genelengths=cnts[,5]
    names(genelengths)=row.names(cnts)
    cnts=cnts[,-(1:5)]
  }
  if(is.null(grouping)){
    grp <- colnames(cnts)
  } else{
    grp <- grouping
  }
  if(!is.null(samples)){
    cnts <- cnts[,samples]
    numsamples <- length(samples)
    grp <- grp[samples]
  }

  cnts <- cnts[order(rownames(cnts)),]

  # make expression tables
  #cpms<-cpm.default(cnts)
  #rpkms<-rpkm.default(cnts,genelengths)
  #cnts <- round(cnts)
  #cnts=cnts[order(row.names(cnts)),]
  #rpkms=rpkms[order(row.names(rpkms)),]
  #log1cnts=log2(1+cnts)
  #logrpkms=log2(1+rpkms)

  # define comparisons
  groups=unique(grp)
  eg=expand.grid(groups,groups,stringsAsFactors=FALSE)
  eg=eg[-which(eg[,1]==eg[,2]),]
  row.names(eg)<-1:nrow(eg)


  if(length(unique(grp))==length(grp)){
    replicated=FALSE
    cat("no replicates found\n")
    if(is.null(dispersion)){stop("must have dispersion value if no replicates, try 0.4")}
  } else{
    replicated=TRUE
  }

  # below assumes no replicates
  dge=DGEList(counts=cnts,group=grp)

  compstrings<-paste0(eg[,2],"_over_",eg[,1],".edger")

  print(compstrings)

  dump<-mclapply(seq_len(nrow(eg)), function(g){


    #compstring=compstring[g]
    #dir.create(compstring)
    s1a<-rowMeans(cnts[,which(grouping==eg[g,1]),drop=F])
    s2a<-rowMeans(cnts[,which(grouping==eg[g,2]),drop=F])
    avg<-data.frame(s1a,s2a)
    colnames(avg)=c(eg[g,1],eg[g,2])

    if(replicated){
      dge <- calcNormFactors     (dge)
      dge <- estimateCommonDisp  (dge)
      dge <- estimateTagwiseDisp (dge)
      et  <- exactTest(dge,pair=c(as.character(eg[g,1]),as.character(eg[g,2])))
    } else{
      et=exactTest(dge,dispersion=dispersion,pair=c(as.character(eg[g,1]),as.character(eg[g,2])))
    }

    ett=et$table
    ett=ett[order(row.names(ett)),]
    ett$QValue=p.adjust(ett$PValue,method="fdr")

    ettt=cbind(row.names(ett),genelengths,cnts,ett)
    colnames(ettt)[1] <- "gene"
    colnames(ettt)[2] <- "length"


    write.tsv(ettt,file=compstrings[g],colnames=TRUE)
  },mc.cores=threads, mc.preschedule=F)

  return(paste0(compstrings))

}
