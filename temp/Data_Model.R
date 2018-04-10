Data_Model <- function(NAME = 'JAKSTAT'){
  if (NAME == 'JAKSTAT') {
    #PARAMETERS <- 10^c(k1=0.31, k2=-1.0, k3=-0.49, k4=0.42, s1=-0.21, s2=-0.34);
    PARAMETERS  <- c(k1=2.4290, k2=975.4280, k3=0.1157, k4=0, s1=10^-0.21, s2=10^-0.34);
    X_0   <- c(x1 =10^0.31, x2 = 0.0, x3 = 0.0, x4 = 0.0)
    dataDir <- "./temp/Data/"
    INPUTDATA <- paste(dataDir,"DATA1_hall_inp.txt", sep="")
    INPUTDATA <- read.table(file=INPUTDATA)
    colnames(INPUTDATA) <- c("time", "EpoRp")
    INPUTDATA$EpoRp[INPUTDATA$time>50] <- 0.009
    
    OUTPUTDATA <- paste(dataDir,"DATA1_hall.txt",sep="")
    OUTPUTDATA <- read.table(file=OUTPUTDATA)
    colnames(OUTPUTDATA) <- c("time", "STAT5p_cyt" ,"sd_STAT5p_cyt","STAT5ptot_cyt","sd_STAT5ptot_cyt")
    
    VARIANCE <- cbind(OUTPUTDATA['time'],OUTPUTDATA['sd_STAT5p_cyt'],OUTPUTDATA['sd_STAT5ptot_cyt'])
    
    OBSERVATIONS <- cbind(OUTPUTDATA['time'],((OUTPUTDATA['STAT5ptot_cyt']/PARAMETERS['s2'])-(OUTPUTDATA['STAT5p_cyt']/PARAMETERS['s1'])),OUTPUTDATA['STAT5ptot_cyt'],OUTPUTDATA['STAT5p_cyt'],(X_0['x1_0']-(OUTPUTDATA['STAT5ptot_cyt']/PARAMETERS['s2']))/2/(1400/450))
    OBSERVATIONS[OBSERVATIONS<0] <- 0
    colnames(OBSERVATIONS) <- c("time", "STAT5" ,"STAT5p_cyt","STAT5ptot_cyt","STAT5_n")
    N <- length(OBSERVATIONS)-1
    
    LIST <- list(OBSERVATIONS,VARIANCE,N,INPUTDATA,PARAMETERS,X_0,'JAKSTAT')
    names(LIST) <- c("observations","variance","N","inputData","parameters","X_0","name")
    return(LIST)}
}
