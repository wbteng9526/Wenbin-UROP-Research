outform <- read.csv(file.choose(),stringsAsFactors = F)
attach(outform)
AIClist <- NULL;BIClist <- NULL;deltaR2list <- NULL;Cplist <- NULL;Adjr2list <- NULL;PRESSlist <- NULL;PRESSnestlist <- NULL

for (i in 1:420) {
  
  # build selecting vectors
  AICq <- outform$AIC[(1000*i-999):(1000*i)]
  BICq <- outform$BIC[(1000*i-999):(1000*i)]
  deltaR2q <- outform$deltaR2[(1000*i-999):(1000*i)]
  Cpq <- outform$Cp[(1000*i-999):(1000*i)]
  Adjr2q <- outform$Adjr2[(1000*i-999):(1000*i)]
  PRESSq <- outform$PRESS[(1000*i-999):(1000*i)]
  PRESSnestq <- outform$PRESSnest[(1000*i-999):(1000*i)]
  
  # Answer q5
  AICq5 <- NULL;BICq5 <- NULL;deltaR2q5 <- NULL;Cpq5 <- NULL;Adjr2q5 <- NULL;PRESSq5 <- NULL;PRESSnestq5 <- NULL
  for (j in 1:4) {
    newAICq5 <- assign(paste("AICv",toString(j),sep=''),length(grep(toString(j),AICq)))
    newBICq5 <- assign(paste("BICv",toString(j),sep=''),length(grep(toString(j),BICq)))
    newdeltaR2q5 <- assign(paste("deltaR2v",toString(j),sep=''),length(grep(toString(j),deltaR2q)))
    newCpq5 <- assign(paste("Cpv",toString(j),sep=''),length(grep(toString(j),Cpq)))
    newAdjr2q5 <- assign(paste("Adjr2v",toString(j),sep=''),length(grep(toString(j),Adjr2q)))
    newPRESSq5 <- assign(paste("PRESSv",toString(j),sep=''),length(grep(toString(j),PRESSq)))
    newPRESSnestq5 <- assign(paste("PRESSnestv",toString(j),sep=''),length(grep(toString(j),PRESSnestq)))
    
    AICq5 <- c(AICq5,newAICq5)
    BICq5 <- c(BICq5,newBICq5)
    deltaR2q5 <- c(deltaR2q5,newdeltaR2q5)
    Cpq5 <- c(Cpq5,newCpq5)
    Adjr2q5 <- c(Adjr2q5,newAdjr2q5)
    PRESSq5 <- c(PRESSq5,newPRESSq5)
    PRESSnestq5 <- c(PRESSnestq5,newPRESSnestq5)
  }
  
  # Answer q1
  AICq1 <- length(which(AICq == paste(c(1:2),collapse = "")))
  BICq1 <- length(which(BICq == paste(c(1:2),collapse = "")))
  deltaR2q1 <- length(which(deltaR2q == paste(c(1:2),collapse = "")))
  Cpq1 <- length(which(Cpq == paste(c(1:2),collapse = "")))
  Adjr2q1 <- length(which(Adjr2q == paste(c(1:2),collapse = "")))
  PRESSq1 <- length(which(PRESSq == paste(c(1:2),collapse = "")))
  PRESSnestq1 <- length(which(PRESSnestq == paste(c(1:2),collapse = "")))
  # Answer q2-q4 & q6
  # Set counter for each protocol
  ctrAICq2 <- 0;ctrAICq3 <- 0;ctrAICq4 <- 0;ctrAICq6 <- 0
  ctrBICq2 <- 0;ctrBICq3 <- 0;ctrBICq4 <- 0;ctrBICq6 <- 0
  ctrdeltaR2q2 <- 0;ctrdeltaR2q3 <- 0;ctrdeltaR2q4 <- 0;ctrdeltaR2q6 <- 0
  ctrCpq2 <- 0;ctrCpq3 <- 0;ctrCpq4 <- 0;ctrCpq6 <- 0
  ctradjr2q2 <- 0;ctradjr2q3 <- 0;ctradjr2q4 <- 0;ctradjr2q6 <- 0
  ctrPRESSq2 <- 0;ctrPRESSq3 <- 0;ctrPRESSq4 <- 0;ctrPRESSq6 <- 0
  ctrPRESSnestq2 <- 0;ctrPRESSnestq3 <- 0;ctrPRESSnestq4 <- 0;ctrPRESSnestq6 <- 0
  
  for (t in 1:1000) {
    
    # split each numeric number
    AICv <- unlist(as.vector(strsplit(as.character(AICq[t]),"")))
    BICv <- unlist(as.vector(strsplit(as.character(BICq[t]),"")))
    deltaR2v <- unlist(as.vector(strsplit(as.character(deltaR2q[t]),"")))
    Cpv <- unlist(as.vector(strsplit(as.character(Cpq[t]),"")))
    Adjr2v <- unlist(as.vector(strsplit(as.character(Adjr2q[t]),"")))
    PRESSv <- unlist(as.vector(strsplit(as.character(PRESSq[t]),"")))
    PRESSnestv <- unlist(as.vector(strsplit(as.character(PRESSnestq[t]),"")))
    
    # Answer q2
    if (all(c(1:2) %in% AICv)==T){ctrAICq2<-ctrAICq2+1} 
    if (all(c(1:2) %in% BICv)==T){ctrBICq2<-ctrBICq2+1} 
    if (all(c(1:2) %in% deltaR2v)==T){ctrdeltaR2q2<-ctrdeltaR2q2+1} 
    if (all(c(1:2) %in% Cpv)==T){ctrCpq2<-ctrCpq2+1}
    if (all(c(1:2) %in% Adjr2v)==T){ctradjr2q2<-ctradjr2q2+1} 
    if (all(c(1:2) %in% PRESSv)==T){ctrPRESSq2<-ctrPRESSq2+1}
    if (all(c(1:2) %in% PRESSnestv)==T){ctrPRESSnestq2<-ctrPRESSnestq2+1}
    
    # Answer q3,q4
    if (any(c(1:2) %in% AICv)==T){
      ctrAICq3<-ctrAICq3+1} else {ctrAICq4 <- ctrAICq4 + 1}
    if (any(c(1:2) %in% BICv)==T){
      ctrBICq3<-ctrBICq3+1} else {ctrBICq4 <- ctrBICq4 + 1}
    if (any(c(1:2) %in% deltaR2v)==T){
      ctrdeltaR2q3<-ctrdeltaR2q3+1} else {ctrdeltaR2q4 <- ctrdeltaR2q4 + 1}
    if (any(c(1:2) %in% Cpv)==T){
      ctrCpq3<-ctrCpq3+1} else {ctrCpq4 <- ctrCpq4 + 1}
    if (any(c(1:2) %in% Adjr2v)==T){
      ctradjr2q3<-ctradjr2q3+1} else {ctradjr2q4 <- ctradjr2q4 + 1}
    if (any(c(1:2) %in% PRESSv)==T){
      ctrPRESSq3<-ctrPRESSq3+1} else {ctrPRESSq4 <- ctrPRESSq4 + 1}
    if (any(c(1:2) %in% PRESSnestv)==T){
      ctrPRESSnestq3<-ctrPRESSnestq3+1} else {ctrPRESSnestq4 <- ctrPRESSnestq4 + 1}
    
    # Answer q3,q6
    if (all((AICv) %in% c(1:2))==T) {ctrAICq6 <- ctrAICq6 + 1}
    if (all((BICv) %in% c(1:2))==T) {ctrBICq6 <- ctrBICq6 + 1}
    if (all((deltaR2v) %in% c(1:2))==T) {ctrdeltaR2q6 <- ctrdeltaR2q6 + 1}
    if (all((Cpv) %in% c(1:2))==T) {ctrCpq6 <- ctrCpq6 + 1}
    if (all((Adjr2v) %in% c(1:2))==T) {ctradjr2q6 <- ctradjr2q6 + 1}
    if (all((PRESSv) %in% c(1:2))==T) {ctrPRESSq6 <- ctrPRESSq6 + 1}
    if (all((PRESSnestv) %in% c(1:2))==T) {ctrPRESSnestq6 <- ctrPRESSnestq6 + 1}
  }
  AICq2 <- ctrAICq2-AICq1;AICq3 <- ctrAICq6-AICq1;AICq6 <- ctrAICq3-AICq1-AICq2-AICq3;AICq4 <- ctrAICq4;
  AIClist <- c(AIClist,AICq1,AICq2,AICq3,AICq6,AICq4,AICq5)
  
  BICq2 <- ctrBICq2-BICq1;BICq3 <- ctrBICq6-BICq1;BICq6 <- ctrBICq3-BICq1-BICq2-BICq3;BICq4 <- ctrBICq4;
  BIClist <- c(BIClist,BICq1,BICq2,BICq3,BICq6,BICq4,BICq5)
  
  deltaR2q2 <- ctrdeltaR2q2-deltaR2q1;deltaR2q3 <- ctrdeltaR2q6-deltaR2q1;deltaR2q6 <- ctrdeltaR2q3-deltaR2q1-deltaR2q2-deltaR2q3;deltaR2q4 <- ctrdeltaR2q4;
  deltaR2list <- c(deltaR2list,deltaR2q1,deltaR2q2,deltaR2q3,deltaR2q6,deltaR2q4,deltaR2q5)
  
  Cpq2 <- ctrCpq2-Cpq1;Cpq3 <- ctrCpq6-Cpq1;Cpq6 <- ctrCpq3-Cpq1-Cpq2-Cpq3;Cpq4 <- ctrCpq4;
  Cplist <- c(Cplist,Cpq1,Cpq2,Cpq3,Cpq6,Cpq4,Cpq5)
  
  Adjr2q2 <- ctradjr2q2-Adjr2q1;Adjr2q3 <- ctradjr2q6-Adjr2q1;Adjr2q6 <- ctradjr2q3-Adjr2q1-Adjr2q2-Adjr2q3;Adjr2q4 <- ctradjr2q4;
  Adjr2list <- c(Adjr2list,Adjr2q1,Adjr2q2,Adjr2q3,Adjr2q6,Adjr2q4,Adjr2q5)
  
  PRESSq2 <- ctrPRESSq2-PRESSq1;PRESSq3 <- ctrPRESSq6-PRESSq1;PRESSq6 <- ctrPRESSq3-PRESSq1-PRESSq2-PRESSq3;PRESSq4 <- ctrPRESSq4;
  PRESSlist <- c(PRESSlist,PRESSq1,PRESSq2,PRESSq3,PRESSq6,PRESSq4,PRESSq5)
  
  PRESSnestq2 <- ctrPRESSnestq2-PRESSnestq1;PRESSnestq3 <- ctrPRESSnestq6-PRESSnestq1;PRESSnestq6 <- ctrPRESSnestq3-PRESSnestq1-PRESSnestq2-PRESSnestq3;PRESSnestq4 <- ctrPRESSnestq4;
  PRESSnestlist <- c(PRESSnestlist,PRESSnestq1,PRESSnestq2,PRESSnestq3,PRESSnestq6,PRESSnestq4,PRESSnestq5)
}



outanalysis <- matrix(NA,nrow=3780,ncol=12);ctra <- 0
for (n in c(50,100,200,400)) {
  for (s in seq(0,1,0.05)) {
    for (delta in c(0.5,0.8,1.0,1.2,1.5)) {
      for (i in 1:9) {
        ctra <- ctra + 1
        if (i < 6) {
          outanalysis[ctra,] <- c(ctra,n,s,delta,paste('Q',toString(i),sep=''),AIClist[ctra],BIClist[ctra],deltaR2list[ctra],Cplist[ctra],Adjr2list[ctra],PRESSlist[ctra],PRESSnestlist[ctra])
        } else {
          outanalysis[ctra,] <- c(ctra,n,s,delta,paste('V',toString(i-5),sep=''),AIClist[ctra],BIClist[ctra],deltaR2list[ctra],Cplist[ctra],Adjr2list[ctra],PRESSlist[ctra],PRESSnestlist[ctra])
        }
      }
    }
  }
}
outanalysis <- as.data.frame(outanalysis,stringsAsFactors = FALSE)
colnames(outanalysis) <- c("sim.ctr","obs","discounter","delta","Qst","AICQ","BICQ","deltaR2Q","CpQ","Adjr2Q","PRESSQ","PRESSnestQ")

write.csv(outanalysis,"C:/Users/thinkpad/Desktop/BU/Research/Prof. Harbaugh/wenbin/output/outanalysis4var05.csv")