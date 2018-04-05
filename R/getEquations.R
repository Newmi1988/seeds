getEquations <- function(model){

  extractModel <- function(model) {
    dpTestModel <- deparse(model,width.cutoff = 500)
    matches <- grepl("\\b[dx,dy][1-9]*",dpTestModel)
    eq <- dpTestModel[matches] 
    listBoo <- grepl("list",eq)
    eq <- eq[!listBoo]
    listBoo = grepl("function", eq)
    eq <- eq[!listBoo]
    return(eq) #return equation without a list statement if present
  }



  formatModelEq <- function(modelEq) {
      trim <- function (x)  sub("^\\s+", "", x)   # get rid of the delimiter
      trimedModelEq <- tolower(gsub("\\[|\\]|[,]", "",trim(modelEq)))
      trimedModelEq <- gsub("([x,y]+)\\s*([0-9]+)","\\1\\2",trimedModelEq)
  
      if(grepl("d+[x,y]",x = trimedModelEq[1])){
        trimDxy <- gsub("(d+[x,y])","",trimedModelEq) # trim the dx/dy
      }
      else {
        trimDxy <- gsub("([y]+)\\s*([0-9]+)","\\2",trimedModelEq)
      }
      if(length(trimDxy)>1) {
        order <- as.integer(substr(x = trimDxy, start = 1, stop = 2)) #convert character of all first elements of the string vector into numbers
        newModelEq <- trimedModelEq[order]
      }
      else {
        newModelEq <- trimedModelEq
      }
      
      
      return(newModelEq)
  }
  
trimSpace <- function (x) gsub("\\s", "",x)
  
getPara <- function(model){
  dpTestModel <- deparse(model,width.cutoff = 500)
  functionHead <- strsplit(dpTestModel[1], split = ",")[[1]]
  paraInd <- trimSpace(functionHead[grepl("para",functionHead)])
  paraInd = gsub(pattern = "[^a-zA-Z]", replacement = "", x = paraInd)
  
  paras <- dpTestModel[grepl(paste0(paraInd,"\\[[0-9]*"),dpTestModel)]
  
  paras = strsplit(x = trimSpace(paras), split = "=|<-")
  paras = unlist(lapply(paras, '[[',1))
  
  return(paras)
}

paras <- as.character(getPara(model))

if(class(model)=="function"){
  modelStr <- formatModelEq(extractModel(model))
  return(list(strM = modelStr, strP = paras))
}
else {
  modelStr <- formatModelEq(model)
  return(list(strM = modelStr, strP = paras))
}


}



