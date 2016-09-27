autoGOrilla <- function (targetFiles, filePrefixes, backgroundFiles=NULL, outputdir=NULL, runmode=NULL,ontology=NULL,pvalue=NULL,name=NULL,email=NULL,includeDups=NULL,fast=NULL, threads=getOption("threads",1L) ){
  if(is.null(outputdir))
    outputdir="."
  if(!is.null(backgroundFiles))
    cmdString=paste0("perl automateGOrilla.pl -targets ", targetFiles, " -background ", backgroundFiles, " -outputdir ", outputdir, " -fileprefix ", filePrefixes )
  else
    cmdString=paste0("perl automateGOrilla.pl -targets ", targetFiles, " -outputdir ", outputdir, " -fileprefix ", filePrefixes )
  if(!is.null(runmode))
    cmdString=paste0(cmdString," -runmode ",runmode)
  if(!is.null(ontology))
    cmdString=paste0(cmdString," -ontology ",ontology)
  if(!is.null(pvalue))
    cmdString=paste0(cmdString," -pvalue ",pvalue)
  if(!is.null(name))
    cmdString=paste0(cmdString," -name ",name)
  if(!is.null(email))
    cmdString=paste0(cmdString," -email ",email)
  if(!is.null(includeDups))
    cmdString=paste0(cmdString," -includedups! ",includeDups)
  if(!is.null(fast))
    cmdString=paste0(cmdString," -fast! ",fast)

  cmdRun(cmdString,threads)
}
