getEquations <- function(model){

extractModel <- function(model) {
  dpTestModel <- deparse(model,width.cutoff = 500)
  #matches <- grepl("d*\\b[x,y]\\[*[1-9]*\\]*",dpTestModel) #finds the equations of the ode
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

if(class(model)=="function"){
  return(formatModelEq(extractModel(model)))
}
else {
  return(formatModelEq(model))
}


}



