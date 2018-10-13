modelStr <- 'BIOMD0000000545_url.xml'
uvb <- importSBML(modelStr)

y <- uvbData[,1:3]

# error in parameters. correction with new vector
uvbParameter = c(  ks1=0.23,
                   ks2=4.0526,
                   kdr1=0.1,
                   kdr2=0.2118,
                   k1=0.0043,
                   k2=161.62,
                   ka1=0.0372,
                   ka2=0.0611,
                   ka3=4.7207,
                   kd1=94.3524,
                   kd2=50.6973,
                   kd3=0.5508,
                   ks3=0.4397,
                   kdr3=1.246,
                   uv=1,
                   ka4=10.1285,
                   kd4=1.1999,
                   n1=3,
                   n2=2,
                   n3=3.5,
                   kdr3a=0.9735,
                   kdr3b=0.406,
                   ksr=0.7537,
                   fhy3=5)

uvb <- setMeas(theObject = uvb, y)
uvb <- setParms(theObject = uvb, uvbParameter)

plot(nominalSol(odeModel = uvb))
