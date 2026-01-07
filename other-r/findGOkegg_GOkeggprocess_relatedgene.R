library(openxlsx)
source("F:\\getGoTerm.R")
GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
# save(GO_DATA, file = "GO_DATA.RData")
findGO <- function(pattern, method = "key"){
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method == "key"){
    pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
  } else if(method == "gene"){
    pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  colnames(pathways) = "pathway"
  if(length(pathways) == 0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
}

getGO <- function(ID){
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)     
  } else{
    cat("No results!\n")
  }
}


# 寻找含有指定关键字的 pathway name 的 pathway
findGO("insulin")
# 寻找含有指定基因名的 pathway
findGO("INS", method = "gene") 
# 获取指定 GO ID 的 gene set
data <- read.xlsx("E:\\temp\\zhangxu\\miRNA_relatedGOprocess.xlsx",colNames = F)


for (i in 1:length(data[,1])){
  datalist <- getGO(data$X1[i])
  datalist <- as.data.frame(datalist)
  write.xlsx(datalist,file=paste(colnames(datalist),".xlsx"),rowNames = FALSE)
}
###kegg_process
getKEGG <- function(ID){
  library("KEGGREST")
  gsList = list()
  for(xID in ID){
    gsInfo = keggGet(xID)[[1]]
    if(!is.null(gsInfo$GENE)){
      geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])   
      xgeneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])          
      NAME = sapply(strsplit(gsInfo$NAME, " - "), function(x) x[1])
      names(xgeneSet) = NAME
      gsList[NAME] = xgeneSet 
    } else{
      cat(" ", xID, "No corresponding gene set in specific database.\n")
    }
  }
  return(gsList)
}




