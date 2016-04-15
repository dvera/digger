dge.heatmap <- function ( counts , genelists , suffix="" , savetables=FALSE , sortby=1, hcluster=FALSE , PDF=TRUE, listnames=NULL, samples=NULL, color=colorRampPalette(c("black","blue","yellow","red"))  ){

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
    gene <- toupper(lists[[x]])
    gene <- gsub("Α","A",gene)
    gene <- gsub(" ","",gene)
    gene <- gsub("-","",gene)
    gene <- unique(gene)
    gene
  })


  generows <- lapply(1:numlists, function(x){
    which(row.names(cnts) %in% genes[[x]])
  })
  badgenes <- lapply(1:numlists, function(x){
    genes[[x]][which(genes[[x]] %ni% row.names(cnts) ) ]

  })



  if(PDF){
    pdf(file=paste0(basename(removeext(counts)),"_log1pRPKM_heatmaps",suffix,".pdf"),width=12,height=12)
  }

  for(i in seq_len(numlists)){

    sublog1rpkms <- log1rpkms[generows[[i]],]

    subrpkms <- rpkms[generows[[i]],]
    if(savetables){ write.tsv(subrpkms, file=paste0(basename(removeext(counts)),"_",listnames[i],".tsv") ,colnames=T, rownames=T)}

    sublog1rpkms <- sublog1rpkms[order(sublog1rpkms[,sortby],decreasing=F),]
    #sublog1rpkms <- sublog1rpkms[order(round(sublog1rpkms[,1]),round(sublog1rpkms[,2]),decreasing=F),]
    if(!PDF){
      tiff(file=paste0(basename(removeext(counts)),"_log1pRPKM_heatmaps",suffix,"_",listnames[i],".tiff"),width=400,height=1000)
    }
    heatmap.2(sublog1rpkms,trace="none",Colv=F,Rowv=if(hcluster){TRUE} else{FALSE},main=paste(listnames[i],"(",nrow(sublog1rpkms),")"), col=color)
    if(!PDF){
      dev.off()
    }
    xx=as.data.frame(sublog1rpkms)
    print(xx)
    write.tsv(xx,file=paste0(basename(removeext(counts)),"_log1pRPKM",suffix,"_",listnames[i],".tsv"),rownames=T,colnames=T)
    #heatmap.2(log1rpkms[generows[[i]],],trace="none",Rowv=F,Colv=F,main=paste(listnames[i],"log( 1 + CPM )"))
  }

  dev.off()

  print(badgenes)

}
