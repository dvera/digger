# significantOnly option
### 0 = false
### 1 = dark and light color bars (Qval <= 0.1)
### 2 = only light color bars (most significant, Qval <= 0.05)
# dges = edger or cuffdiff (cuffdiff=T) files

dge.bar <- function( dges , genelists=NULL , listnames=NULL , sortByChange=TRUE , debug=getOption("verbose"), cuffdiff=FALSE, significantOnly=0, savetables=FALSE, logTwo=TRUE ){

  #library(edgeR)
  #library(DESeq2)
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
  if(significantOnly < 0) {
    cat(">> Invalid significantOnly Value! Defaulting to 0\n")
    significantOnly = 0
  } else if(significantOnly > 2) {
    cat(">> Invalid significantOnly Value! Defaulting to 2\n")
    significantOnly = 2
  }

  numdes <- length(dges)

  if(is.null(listnames)){
    listnames = basename((genelists))
  }


  for(d in seq_len(numdes)){
    if(logTwo) {
      ext="_log2bar.pdf"
    } else {
      ext="_bar.pdf"
    }
    pdf(file=paste0(basename(removeext(dges[d])),"_sig",significantOnly,ext))
    de <- read.tsv(dges[d],header=T)

    if(cuffdiff){
      colnames(de)[13] <- "QValue"
      colnames(de)[10] <- "logFC"
      de$logFC[is.infinite(de$logFC)]<-10
    }

    # find differentially-expressed genes
    #deg=which(de$QValue<=0.05)
    upreg=which(de$QValue<=0.05 & de$logFC>0)
    doreg=which(de$QValue<=0.05 & de$logFC<0)
    upreg2=which(de$QValue>0.05 & de$QValue <= 0.1 & de$logFC>0)
    doreg2=which(de$QValue>0.05 & de$QValue <= 0.1 & de$logFC<0)
    numupreg<-length(upreg)
    numdoreg<-length(doreg)


    # assign colors based on status
    de$color="black"
    de$color[upreg]="green"
    de$color[doreg]="red"
    de$color[upreg2]="darkgreen"
    de$color[doreg2]="darkred"




    if(!is.null(genelists)){

      numlists <- length(genelists)
      lists <- lapply( genelists, read_tsv, col_names=F )
      genes <- lapply(1:numlists,function(x){

        if(ncol(lists[[x]])>2){
          unique(sort(lists[[x]][ which(grepl("HUMAN|",lists[[x]][,11] ) ),3]))
        } else{
          unique(unlist(as.vector(lists[[x]])))
        }

      })
      generows <- lapply(1:numlists, function(x){
        which(toupper(de$gene) %in% toupper(gsub(" ","",gsub("Α","A",genes[[x]]))))
      })
    }


    par(mfrow=c(2,1), mar=c(4,4,2,1), xpd=TRUE)

    for(i in seq_len(numlists)){

      subde <- de[generows[[i]],]

      if(debug){
        cat("i=",i,"\n")
        print(subde)
      }



      if(sortByChange){subde<-subde[order(subde$logFC, decreasing=T),]}
      if(significantOnly > 0) {
        w=which(subde$color=="black")
        if(length(w) != 0)
          subde<-subde[-w,]
        if(significantOnly == 2) {
          w=which(subde$color=="darkgreen")
          if(length(w) != 0)
            subde<-subde[-w,]
          w=which(subde$color=="darkred")
          if(length(w) != 0)
            subde<-subde[-w,]
        }
      }

      if(savetables){ write.tsv(subde, file=paste0(basename(removeext(dges[d])),"_",removeext(listnames[i]),".edger"), colnames=T )}   # eh....


      if(nrow(subde)>0){
        if(logTwo) {
          barplot(
            subde$logFC,
            #log1cnts[,which(colnames(cnts)==eg[g,2])],
            main=removeext(listnames[i]),
            col=subde$color,
            names=subde$gene,
            las=3,
            #pch=20,
            #xlab="log2(mean of counts)",
            #ylab=paste0("log2(",colnames(de)[3],"/",colnames(de)[2],")")
            ylab="log2(Fold Change)"
          )
          abline(h=0,v=(1:nrow(subde)*1.2) , col=rgb(0,0,0,20,maxColorValue=255) )
          #par(mar=c(2,2,2,1), xpd=FALSE)
        } else {
          barplot(
            2^subde$logFC,
            #log1cnts[,which(colnames(cnts)==eg[g,2])],
            main=removeext(listnames[i]),
            col=subde$color,
            names=subde$gene,
            las=3,
            #pch=20,
            #xlab="log2(mean of counts)",
            #ylab=paste0(colnames(de)[3],"/",colnames(de)[2],")")
            ylab="Fold Change"
          )
          abline(h=0,v=(1:nrow(subde)*1.2) , col=rgb(0,0,0,20,maxColorValue=255) )
          #par(mar=c(2,2,2,1), xpd=FALSE)
        }
        barplot(
          log10(subde$QValue),
          #log1cnts[,which(colnames(cnts)==eg[g,2])],
          #main=basename(genelists[[i]]),
          #col=ett$color,
          #names=subde$gene,
          col=subde$color,
          las=3,
          #pch=20,
          #xlab="log2(mean of counts)",
          ylab="log10( q-value )"
        )

        abline(h=c(-1,-1.3), v=(1:nrow(subde)*1.2) , col=rgb(0,0,0,20,maxColorValue=255) )
      }
      # subcnts <- cnts[generows[[i]],]
      # subrpkms <- rpkms[generows[[i]],]
      # sublogrpkms <- logrpkms[generows[[i]],]
      # subdelog1cnts <- log1cnts[generows[[i]],]


    }
    dev.off()
  }

}
