limma_xys <- function( xysFiles , outname , contrasts=NULL , grouping=NULL , annos=NULL , pval=0.05 , adjust="BH" ){
  #colnames of xysFiles should be FileName and Target

  #design <- model.matrix(~0+factor(xysFiles$Target))
  u=unique(xysFiles$Target)
  eg <- expand.grid(u,u)
  eg <- eg[-which(eg[,1]==eg[,2]),]
  cont <- paste0(eg[,1],"-",eg[,2])
  f=factor(u)
  design <- model.matrix(~0+xysFiles$Target)
  #colnames(design) <- gsub("\\$","!",colnames(design))
  colnames(design) <- remove.prefix(colnames(design),"Target")
  cont.matrix <- makeContrasts(contrasts=cont,levels=u)

  if(is.data.frame(xysFiles)){
    grp=xysFiles[,2]
    xys=xysFiles[,1]
  } else{
    grp=grouping
  }

  numfiles = length(xys)


  ###
  # ug <- unique(grp)
  #
  # design <- matrix(rep(0,length(grp)*2),ncol=2)
  # design[which(grp==ug[1]),1]<-1
  # design[which(grp==ug[2]),2]<-1
  # colnames(design)<-ug

  ###

  maqc=read.xysfiles(xys)
  eset=rma(maqc)
  fit <- lmFit(eset, design)
  #cont.matrix <- makeContrasts(test=paste0(ug[2],"-",ug[1]), levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  eb <- eBayes(fit2)
  tt <- topTable(eb,n=5000000, adjust=adjust )
  e <- exprs(eset)

  tt <- tt[order(rownames(tt)),]
  e <- e[order(rownames(e)),]



  newt=cbind(probe=rownames(e),e,tt)
  if(!identical(colnames(newt)[2:(numfiles+1)],basename(xys))){stop("somethings wrong")}

  colnames(newt)[2:(numfiles+1)] <- paste0(grp,"__",basename(xys))


  if(!is.null(annos)){


    anno=tsvRead(annos)

    if(is.data.frame(anno)){
      anno <- list(anno)
    }

    for(i in 1:length(annos)){
      annos2=anno[[i]][,c(4,7)]
      annos3=unique(annos2)
      annos4=aggregate(annos3[,2],by=list(annos3[,1]),FUN=paste0)
      m=match(newt$probe,annos4[,1])
      newt[[basename(removeext(annos[i]))]] <- annos4[m,2]
      #newt2=merge(newt,annos4,by.x="probe",by.y="Group.1",all.x=T,all.y=F)
      newt[[basename(removeext(annos[i]))]] <- unlist(lapply(newt[[basename(removeext(annos[i]))]],paste0,collapse=","))
      #colnames(newt)[which(colnames(newt)=="x")] <- "xenoRef"
    }
  }

  pos24=newt[which(newt$P.Value<=pval & newt$logFC >=2 & newt$logFC < 4 ),]
  neg24=newt[which(newt$P.Value<=pval & newt$logFC <=(-2) & newt$logFC > -4 ),]
  neg4=newt[which(newt$P.Value<=pval & newt$logFC < -4 ),]
  pos4=newt[which(newt$P.Value<=pval & newt$logFC > 4 ),]


  #outname=paste0(ug[1],"_vs_",ug[2],".tsv")
  pos24name=paste0(ug[1],"_vs_",ug[2],"_upreg_logFc2to4.tsv")
  neg24name=paste0(ug[1],"_vs_",ug[2],"_downreg_logFc2to4.tsv")
  pos4name=paste0(ug[1],"_vs_",ug[2],"_upreg_logFcGreaterThan4.tsv")
  neg4name=paste0(ug[1],"_vs_",ug[2],"_downreg_logFcGreaterThan4.tsv")

  tsvWrite(newt,file=outname,col_names=T)
  if(nrow(pos24)>0){ tsvWrite(pos24,file=pos24name,col_names=T) }
  if(nrow(neg24)>0){ tsvWrite(neg24,file=neg24name,col_names=T) }
  if(nrow(pos4)>0){ tsvWrite(pos4,file=pos4name,col_names=T) }
  if(nrow(neg4)>0){ tsvWrite(neg4,file=neg4name,col_names=T) }

}





  # par(mfrow=c(2,2))
  #
  # x1= rowMeans(e[,which(design[,1]==1)])
  # x2= rowMeans(e[,which(design[,2]==1)])
  #
  # plot(x1,x2,col=tt$color,pch=20,xlab=ug[1],ylab=ug[2])
  # plot(tt$AveExpr,tt$logFC,col=tt$color,xlab="mean expression",ylab="log2 ratio")
  # plot(density(as.vector(e)))
  # plot(density(as.vector(e)))



  # tt$color="black"
  # tt$color[which(tt$logFC<0 & tt$P.Value <=0.05)]<-"darkred"
  # tt$color[which(tt$logFC<0 & tt$adj.P.Val <=0.05)]<-"darkred"
  # tt$color[which(tt$logFC>0 & tt$P.Value <=0.05)]<-"darkgreen"
  # tt$color[which(tt$logFC>0 & tt$adj.P.Val <=0.05)]<-"green"
