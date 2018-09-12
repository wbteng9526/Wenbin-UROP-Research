library(dplyr)
library(MPV)
sdxvec <- function(sdx) {
  sdxvector <- NULL 
  for (i in sdx){
    for (j in sdx){
      newsdx <- (i*j)^0.5
      sdxvector <- c(sdxvector,newsdx)
    }
  }
  return(sdxvector)
}
#chisquare test functions
chitestAIC <- function(outmat,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12) {
  lm.cur <- lm(y ~ 1);ctr <- 0; flag <- TRUE;
  while (flag == TRUE & ctr < length(outmat)) {
    ctr <- ctr + 1
    tmp.var <- names(select(outmat[ctr,],which((outmat[ctr,]==1)==T)))
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~.",tmp.var),collapse="+")))
    chiAIC <- AIC(lm.cur)-AIC(lm.nxt)
    if (ctr==1) {
      flag <- TRUE
    } else {
      if (chiAIC < 0) {
        flag <- FALSE
      } else {
        flag <- pchisq(chiAIC, df = 1,lower.tail = F) < 0.05
      }
    }
    lm.old <- lm.cur
    lm.cur <- lm.nxt
  }
  if (flag == TRUE){
    lm.final <- lm.cur
  } else {lm.final <- lm.old}
  return(lm.final)
}

chitestBIC <- function(outmat,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12) {
  lm.cur <- lm(y ~ 1);ctr <- 0; flag <- TRUE;
  while (flag == TRUE & ctr < length(outmat)) {
    ctr <- ctr + 1
    tmp.var <- names(select(outmat[ctr,],which((outmat[ctr,]==1)==T)))
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~",tmp.var),collapse="+")))
    chiBIC <- BIC(lm.cur)-BIC(lm.nxt)
    if (ctr==1) {
      flag <- TRUE
    } else {
      if (chiBIC < 0) {
        flag <- FALSE
      } else {
        flag <- pchisq(chiBIC, df = 1,lower.tail = F) < 0.05
      }
    }
    lm.old <- lm.cur
    lm.cur <- lm.nxt
  }
  if (flag == TRUE){
    lm.final <- lm.cur
  } else {lm.final <- lm.old}
  return(lm.final)
}


# get the ideal outmat table
getoutmat <- function(outmat){
  outmat <- data.frame(data.table(data.matrix(data.frame(outmat))-1))
  for (i in 1:ncol(outmat)){
    if (prod(outmat[,i] == 0)==1){
      outmat[,i] <- 1
    }
  }
  return(outmat)
}

# nested models determination
nest <- function(outmat){
  tmp.varcheck <- sort(apply(outmat,2,sum),decreasing = T)
  if (prod(sort(apply(outmat,2,sum),decreasing = T) == length(outmat):1) == 1){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
 

# F-test for nested models
FdeltaR2 <- function(outmat,n,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12){
  lm.cur <- lm(y ~ 1);ctr <- 0; flag <- TRUE
  while (flag == TRUE & ctr < length(outmat)) {
    ctr <- ctr + 1
    tmp.var <- names(select(outmat[ctr,],which((outmat[ctr,]==1)==T)))
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~.",tmp.var),collapse="+")))
    Frr <- (summary(lm.nxt)$r.squared-summary(lm.cur)$r.squared)/((1-summary(lm.nxt)$r.squared)/(n-ctr-1))
    if (ctr==1) {
      flag <- TRUE
    } else {
      if (Frr < 0){
        flag <- FALSE
      } else {
        flag <- pf(Frr, df1=1,df2=n-ctr-1,lower.tail = F) < 0.05
      }
    }
    
    lm.old <- lm.cur
    lm.cur <- lm.nxt
  }
  if (flag == TRUE){
    lm.final <- lm.cur
  } else {lm.final <- lm.old}
  return(lm.final)
}
# nested cp-p model selection
Fcp_pnest <- function(outmat,outsummary,n,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12){
  lm.cur <- lm(y ~ 1);ctr <- 0; flag <- TRUE;
  while (flag == TRUE & ctr < length(outmat)) {
    ctr <- ctr + 1
    tmp.var <- names(select(outmat[ctr,],which((outmat[ctr,]==1)==T)))
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~",tmp.var),collapse="+")))
    
    if (ctr==1) {
      flag <- TRUE
    } else{
      Fcpnest <- ((outsummary$cp[ctr]-ctr)-(outsummary$cp[ctr-1]-(ctr-1)))/((outsummary$cp[ctr]-ctr)/(n-ctr-1))
      if (Fcpnest > 0){
        flag <- FALSE
      } else {
        flag <- pf(-Fcpnest, df1=1,df2=n-ctr-1,lower.tail = F) < 0.05
      }
    }
    
    lm.old <- lm.cur
    lm.cur <- lm.nxt
  }
  if (flag == TRUE){
    lm.final <- lm.cur
  } else {lm.final <- lm.old}
  return(lm.final)
}

# for nested models, adjusted R-squared
Fadjr2nest <- function(outmat,outsummary,n,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12){
  lm.cur <- lm(y ~ 1);ctr <- 0; flag <- TRUE;
  while (flag == TRUE & ctr < length(outmat)) {
    ctr <- ctr + 1
    tmp.var <- names(select(outmat[ctr,],which((outmat[ctr,]==1)==T)))
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~",tmp.var),collapse="+")))
    
    if (ctr==1) {
      flag <- TRUE
    } else{
      F_adjr2nest <- (outsummary$adjr2[ctr]-outsummary$adjr2[ctr-1])/(outsummary$adjr2[ctr]/(n-ctr-1))
      if (F_adjr2nest < 0){
        flag <- FALSE
      } else {
        flag <- pf(F_adjr2nest, df1=1,df2=n-ctr-1,lower.tail = F) < 0.05
      }
    }
    
    lm.old <- lm.cur
    lm.cur <- lm.nxt
  }
  if (flag == TRUE){
    lm.final <- lm.cur
  } else {lm.final <- lm.old}
  return(lm.final)
}

# for nested models, PRESS
Fpressnest <- function(outmat,n,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12){
  lm.cur <- lm(y ~ 1);ctr <- 0; flag <- TRUE;
  while (flag == TRUE & ctr < length(outmat)) {
    ctr <- ctr + 1
    tmp.var <- names(outmat[ctr,][which((outmat[ctr,]==1)==T)])
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~",tmp.var),collapse="+")))
    
    if (ctr==1) {
      flag <- TRUE
    } else{
      F_pressnest <- (PRESS(lm.nxt)-PRESS(lm.cur))/(PRESS(lm.nxt)/(n-ctr-1))
      if (F_pressnest > 0){
        flag <- FALSE
      } else {
        flag <- pf(-F_pressnest, df1=1,df2=n-ctr-1,lower.tail = F) < 0.05
      }
    }
    
    lm.old <- lm.cur
    lm.cur <- lm.nxt
  }
  if (flag == TRUE){
    lm.final <- lm.cur
  } else {lm.final <- lm.old}
  return(lm.final)
}
# for non-nested models, choose the smallest PRESS
Fpress <- function(outmat,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12) {
  lm.cur <- lm(y ~ 1);ctr <- 0
  
  while (ctr < length(outmat)) {
    ctr <- ctr + 1
    tmp.var <- names(outmat[ctr,][which((outmat[ctr,]==1)==T)])
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~",tmp.var),collapse="+")))
    print(PRESS(lm.nxt))
    if (ctr==1) {
      lm.min <- lm.nxt
      lm.cur <- lm.nxt
    } else {
      F_press <- PRESS(lm.nxt)-PRESS(lm.min)
      if (F_press <= 0) {
        lm.min <- lm.nxt
        lm.cur <- lm.nxt
      } else {
        lm.min <- lm.cur
        lm.cur <- lm.nxt
      }
    }
  }
  return(lm.min)
}
# for non-nested models, choose the smallest cp-p value
Fcp_p <- function(outmat,outsummary,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12) {
  lm.cur <- lm(y ~ 1);ctr <- 0
  ctrend <- which.min(outsummary$cp-c(1:length(outmat)))
  while (ctr < ctrend) {
    ctr <- ctr + 1
    tmp.var <- names(select(outmat[ctr,],which((outmat[ctr,]==1)==T)))
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~",tmp.var),collapse="+")))
    lm.cur <- lm.nxt
    lm.final <- lm.cur
  }
  return(lm.final)
}

# for non-nested models, choose the smallest cp-p value
Fadjr2 <- function(outmat,outsummary,y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12) {
  lm.cur <- lm(y ~ 1);ctr <- 0
  ctrend <- which.max(outsummary$adjr2)
  while (ctr < ctrend) {
    ctr <- ctr + 1
    tmp.var <- names(select(outmat[ctr,],which((outmat[ctr,]==1)==T)))
    lm.nxt <- update(lm.cur,as.formula(paste(c(".~",tmp.var),collapse="+")))
    lm.cur <- lm.nxt
    lm.final <- lm.cur
  }
  return(lm.final)
}

modelcombineAIC <- function(outmat,newAIC){
  index10 <- match(10,newAIC)
  index11 <- match(11,newAIC)
  index12 <- match(12,newAIC)
  newAICindex <- replace(newAIC,index10,"A")
  newAICindex <- replace(newAICindex,index11,"B")
  newAICindex <- replace(newAICindex,index12,"C")
  newpaste <- paste(newAICindex,collapse = "")
  return(newpaste)
}

modelcombineBIC <- function(outmat,newBIC){
  index10 <- match(10,newBIC)
  index11 <- match(11,newBIC)
  index12 <- match(12,newBIC)
  newBICindex <- replace(newBIC,index10,"A")
  newBICindex <- replace(newBICindex,index11,"B")
  newBICindex <- replace(newBICindex,index12,"C")
  newpaste <- paste(newBICindex,collapse = "")
  return(newpaste)
}

modelcombinenewdeltaR2 <- function(outmat,newdeltaR2){
  index10 <- match(10,newdeltaR2)
  index11 <- match(11,newdeltaR2)
  index12 <- match(12,newdeltaR2)
  newdeltaR2index <- replace(newdeltaR2,index10,"A")
  newdeltaR2index <- replace(newdeltaR2index,index11,"B")
  newdeltaR2index <- replace(newdeltaR2index,index12,"C")
  newpaste <- paste(newdeltaR2index,collapse = "")
  return(newpaste)
}

modelcombinecp_p <- function(outmat,newcp_p){
  index10 <- match(10,newcp_p)
  index11 <- match(11,newcp_p)
  index12 <- match(12,newcp_p)
  newcp_pindex <- replace(newcp_p,index10,"A")
  newcp_pindex <- replace(newcp_pindex,index11,"B")
  newcp_pindex <- replace(newcp_pindex,index12,"C")
  newpaste <- paste(newcp_pindex,collapse = "")
  return(newpaste)
}

modelcombineadjr2 <- function(outmat,newadjr2){
  index10 <- match(10,newadjr2)
  index11 <- match(11,newadjr2)
  index12 <- match(12,newadjr2)
  newadjr2index <- replace(newadjr2,index10,"A")
  newadjr2index <- replace(newadjr2index,index11,"B")
  newadjr2index <- replace(newadjr2index,index12,"C")
  newpaste <- paste(newadjr2index,collapse = "")
  return(newpaste)
}

