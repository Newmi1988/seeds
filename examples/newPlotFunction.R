load('plotExample.RData')

graphics.off()

# normales Plotten
plot(res)

# plotten mit Annotationen
annotationStates <- c("CS", "CD", "CDCS", "UVR8M", "UCS", "UVR8D", "RUP", "UR","","", "HY5", "DWD", "CDW")
annotationMeas <- c("UM total","COP1 total", "UVR8D obs.", "HY5 obs", "UVR8M obs")

plotAnno(res, stateAnno = annotationStates, measAnno = annotationMeas)
