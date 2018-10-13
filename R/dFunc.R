# write a dummy function to parse string to odeModel class

writeDummy <- function(eqList){
  
  funcStr = paste(paste0('dx',1:length(eqList$reac)),eqList$reac, sep = ' = ')
  measStr = paste(paste0('y',1:length(eqList$meas)),eqList$meas, sep = ' = ')
  
  writeDummyFile <- function(str){
    formatStr <- gsub(pattern = "[d]{1}||\\d",'',strsplit(split = ' = ',str)[[1]][1])
    
    if(formatStr == 'x'){
      funcHead <- 'modelFunc <- function(t,x,para){ \n with (as.list(parameters),{ \n'
      funcBody <- paste0('\t',str)
      funcEnd <- paste('\n\tlist(c(',paste0('x',1:length(str),collapse = ','),'))\n})\n}')
    } else {
      funcHead <- 'measFunc <- function(x) {\n'
      funcBody <- paste0('\t',gsub('(x{1})([0-9]+)', replacement = '\\1[,\\2]',str))
      funcEnd <- paste('\n\treturn(list(',paste0('y',1:length(str),collapse = ','),')) \n}')
    }
     
    str <- c(funcHead,funcBody,funcEnd)
    
    file.create('tmp.R')
    tmp <- file('tmp.R')
    writeLines(str,tmp)
    close(tmp)
    
    source('tmp.R')
    if(formatStr == 'x'){
      return(modelFunc)
    } else {
      return(measFunc) 
    }
    
  }
  
  reacFunc <- writeDummyFile(funcStr)
  measFunc <- writeDummyFile(measStr)
   
  return(list("reac"= reacFunc,'meas'= measFunc ))
}
