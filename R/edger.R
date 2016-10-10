edger <- function( counts , grouping=NULL , samples=NULL , tabletype="featureCounts", dispersion=NULL, rpkmout=F, regOut=F, pva;l=FALSE, threads=getOption("threads",1L)  ){
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
#    rownames(cnts)=cnts[,1]
#    colnames(cnts)=cnts[1,]
#    cnts=cnts[-1,-1]
    genelengths=cnts[,5]
    names(genelengths)=row.names(cnts)
    cnts=cnts[,-(1:5)]
  }
  if(is.null(grouping)){
    grp <- colnames(cnts)
  } else{
    if(length(grouping) != length(cnts) ){stop("all samples in counts table must have grouping")}
    grp <- grouping
  }
  if(!is.null(samples)){
    cnts <- cnts[,samples]
    numsamples <- length(samples)
    grp <- grp[samples]
  }
  if(rpkmout) {
    cnts=data.matrix(cnts)
    rpkms=rpkm.default(cnts,as.numeric(genelengths))			# make rpkm counts table
    fo=paste0( removeext( basename(counts) ), "_rpkm.tsv")
    cat(paste0("rpkm table: ",fo,"\n"))
    t=cbind(rownames(rpkms),rpkms)
    colnames(t)[1]="gene"
    tsvWrite(as.data.frame(t),fo,col_names=T)	# print rpkm table
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
  cnts=data.matrix(cnts)
  dge=DGEList(counts=cnts,group=grp)

  if(replicated) {
    compstrings<-paste0(eg[,2],"_over_",eg[,1],".edger")
  }
  else {
    compstrings<-paste0(eg[,2],"_over_",eg[,1],"_disp",dispersion,".edger")
  }

  print(compstrings)

  dump<-mclapply(seq_len(nrow(eg)), function(g){


    #compstring=compstring[g]
    #dir.create(compstring)
    s1a<-rowMeans(cnts[,which(grp==eg[g,1]),drop=F])
    s2a<-rowMeans(cnts[,which(grp==eg[g,2]),drop=F])
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

  if(regOut) {
    dump<-mclapply(1:length(compstrings), function(g) {
    x=tsvRead(compstrings[g],col_names=T)
    if(pval) {
      upreg=which(x$PValue<=0.05 & x$logFC>0)
      downreg=which(x$PValue<=0.05 & x$logFC<0)
      s="_Pval_"
    }
    else {
      upreg=which(x$QValue<=0.05 & x$logFC>0)
      downreg=which(x$QValue<=0.05 & x$logFC<0)
      s="_Qval_"
    }
    xup=x[upreg,]
    xdown=x[downreg,]
    tsvWrite(xup,file=paste0(remove.suffix(compstrings[g],".edger"),s,"upreg.txt"))
    tsvWrite(xdown,file=paste0(remove.suffix(compstrings[g],".edger"),s,"downreg.txt"))
    tsvWrite(as.data.frame(xup$gene),file=paste0(remove.suffix(compstrings[g],".edger"),s,"upreg_genes.txt"))
    tsvWrite(as.data.frame(xdown$gene),file=paste0(remove.suffix(compstrings[g],".edger"),s,"downreg_genes.txt"))
    },mc.cores=threads, mc.preschedule=F)
  }

  return(paste0(compstrings))

}
