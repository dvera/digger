# organism input possible:
###       ARABIDOPSIS_THALIANA
###       SACCHAROMYCES_CEREVISIAE
###       CAENORHABDITIS_ELEGANS
###       DROSOPHILA_MELANOGASTER
###       DANIO_RERIO
###       HOMO_SAPIENS  (default)
###       MUS_MUSCULUS
###       RATTUS_NORVEGICUS
# runmodes input possible:
###       mhg (single ranked list of genes)
###       hg  (target and background gene lists)
# ontologies input possible:
###       proc (process)
###       func (function)
###       comp (component)
###       all
# includeDups, revigo, fast input possible:
###       0   (false)
###       1   (true)

autoGOrilla <- function (targetFiles, backgroundFiles=NULL, outputdir=NULL, organism="HOMO_SAPIENS",runMode="mhg",ontology="all",pValue=0.01,name=NULL,email=NULL,includeDups=0,revigo=1,fast=1, scpPath=NULL, threads=getOption("threads",1L) ){

  scriptURL="http://www.veralab.org/kkyle/automateGOrilla.pl"
   l=list.files("~",pattern="automateGOrilla.pl",recursive=T)
   if(length(l)==0){
     cat("\nScript not found! Downloading...\n\n")
     l="./automateGOrilla.pl"
     download.file(scriptURL,destfile=l,method="wget")
   }
   else
     l=paste0("~/",l)

  if(is.null(outputdir))
    outputdir="."

  cmdString=paste("perl", l[1], "-targets", targetFiles, "-outputdir", outputDir)

  if(!is.null(backgroundFiles))
    cmdString=paste(cmdString, "-background", backgroundFiles)

  cmdString=paste(cmdString,"-organism",organism,"-runmode",runMode,"-ontology",ontology,"-pvalue",pValue,"-includedups",includeDups,"-revigo",revigo,"-fast",fast)


  if(!is.null(name))
    cmdString=paste(cmdString,"-name",name)
  if(!is.null(email))
    cmdString=paste(cmdString,"-email",email)

  print(cmdString)

  cmdRun(cmdString,threads)

  if(!is.null(scpPath))
    cmdString=paste("scp -r",outputDir,scpPath)

  cmdRun(cmdString,threads)

}
