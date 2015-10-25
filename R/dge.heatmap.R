dge.heatmap <- function ( counts , genelists , savetables=FALSE , listnames=NULL, samples=NULL, color=colorRampPalette(c("black","blue","yellow","red"))  ){

  library(edgeR)

  cnts=read.tsv(counts,header=TRUE,row.names=1)
  rownames(cnts) <- toupper(rownames(cnts))
  genelengths=cnts[,5]
  names(genelengths)=row.names(cnts)
  cnts=cnts[,-(1:5)]
  cnts=data.matrix(cnts)
  if(!is.null(samples)){ cnts=cnts[,samples] }

  cpms<-cpm.default(cnts)
  rpkms<-rpkm.default(cnts,genelengths)
  log1cnts=log2(1+cnts)
  log1rpkms=log2(1+rpkms)

  if(is.null(listnames)){
    listnames = basename((genelists))
  }

  numlists <- length(genelists)
  lists <- lapply(lapply( genelists, read_tsv, col_names=F ),unlist)
  genes <- lapply(1:numlists,function(x){
    unique(toupper(lists[[x]]))
  })
  generows <- lapply(1:numlists, function(x){
    which(toupper(row.names(cnts)) %in% toupper(gsub(" ","",gsub("Α","A",genes[[x]]))))
  })



  pdf(file=paste0(basename(removeext(counts)),"_log1pRPKM_heatmaps.pdf"),width=12,height=12)

  for(i in seq_len(numlists)){

    sublog1rpkms <- log1rpkms[generows[[i]],]

    subrpkms <- rpkms[generows[[i]],]
    if(savetables){ write.tsv(subrpkms, file=paste0(basename(removeext(counts)),"_",listnames[i],".tsv") ,colnames=T, rownames=T)}


    heatmap.2(sublog1rpkms,trace="none",Colv=F,main=paste(listnames[i]), col=color)
    #heatmap.2(log1rpkms[generows[[i]],],trace="none",Rowv=F,Colv=F,main=paste(listnames[i],"log( 1 + CPM )"))
  }

  dev.off()


}
