# upfiles = plain text gene list of upreg genes, 1 gene per row in file
# downfiles = plain text gene list of downreg genes, 1 gene per row in file
# outValues = if true, creates files of genes in each category of the venn diagram
# valString = suffixes for files created if outValues true
# catNames is passed to venn.diagram to label the categories of the diagram
# ... passed to venn.diagram() function
# depends on conifur library
deVenn <- function(upfiles,downfiles,pdfName=NULL, outValues=T,valStrings=letters[1:length(upfiles)],catNames=letters[1:length(upfiles)],threads=getOption("threads",1L), ...) {
  library(conifur)
  library(VennDiagram)

  if(is.null(pdfName))
    pdfName="deVenn.pdf"

  up=tsvRead(upfiles,threads=3,drop=F)
  down=tsvRead(downfiles,threads=3,drop=F)
#  x=list(ba=upva[[1]],ca=upva[[2]],da=upva[[3]])
  uup=unlist(up, recursive=F)
  udown=unlist(down,recursive=F)
#  p=get.venn.partitions(tt,keep.elements=F)
  uupp=get.venn.partitions(uup)
  udownp=get.venn.partitions(udown)
#  v=get.venn.partitions(tt)$..values..
  pdf(pdfName)
  venn.plot <- venn.diagram(uup,NULL,,main="Upregulated",category.names=catNames,...)
  grid.draw(venn.plot)
  grid.newpage()
  venn.plot <- venn.diagram(udown,NULL,,main="Downregulated",category.names=catNames,...)
  grid.draw(venn.plot)
  dev.off()

  # writing intersection lists to files
  if(outValues) {
    upvals=paste0("up_",valStrings,".tsv")
    dwnvals=paste0("down_",valStrings,".tsv")
    invisible( mclapply(1:length(upvals),function(i) {
      tsvWrite(as.data.frame(uupp$..values..[i]),upvals[i])
      tsvWrite(as.data.frame(udownp$..values..[i]),dwnvals[i])
    },mc.cores=threads,mc.preschedule=F) )
  }

} # end function
