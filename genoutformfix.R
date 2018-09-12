library(MASS)
library(leaps)
library(data.table)
library(base)
library(locfit)
library(Matrix)
library(tibble)
library(MPV)
library(parallel)
suppressMessages(library(foreach))
library(iterators)
library(doMC)



# determine the number of cores to use in this device
num_cores <- detectCores()
# num_cores <- num_cores - 1
max_cores <- 20
num_cores <- min(num_cores,max_cores)
registerDoMC(num_cores)

# Set seeds for each cores that would be using in parallel loop
set.seed(1234)
my.new.seeds <- NULL
for (i in 1:num_cores) {
  tmp.ind <-  floor(runif(1,2,627))
  next.rand <- .Random.seed[tmp.ind]
  my.new.seeds <- c(my.new.seeds,next.rand)
}


source('functions.R')

intb <- 5
mxlist <- rep(0,20)
sdxlist <- rep(1,20)

parallel.4var.list <- foreach(core.ctr = 1:num_cores) %dopar% {
  
  set.seed(my.new.seeds[core.ctr])
  
  ctr <- 0
  out4var <- matrix(NA,nrow = 420000,ncol = 10)
  for (n in c(50,100,200,400)) {
    for (s in seq(0,1,0.05)) {
      for (delta in c(0.5,0.8,1.0,1.2,1.5)) {
        for (i in 1:50) {
          # data simulation
          
          mx <- mxlist[1:4]
          sdx <- sdxlist[1:4]
          sdxvector <- sdxvec(sdx) # a new function created in source script
          
          # generate 3-type vars correlation matrix
          corr <- c(1,0,0.5,0,0,1,0.5,0,0.5,0.5,1,0,0,0,0,1)
          covmatrix <- matrix(corr*sdxvector,4,4)
          x <- mvrnorm(n,mx,covmatrix)
          sde <- delta*(sqrt(1-0.1024*(1+s^2)))
          e <- rnorm(n,0,sde)
          newdata <- data.frame(x)
          for (h in 1:4){
            assign(paste("x",toString(h),sep=''),x[,h])
            names(newdata[,h]) <- paste("x",toString(h),sep='')
          }
          y <- intb + s*0.32*x[,1] + 0.32*x[,2] + e
          newdata <- cbind(y,newdata)
          
          # multiple regression
          mod <- y ~ 1
          for (p in 1:4){
            mod <-update(mod, as.formula(paste(".~.+",paste("x",toString(p),sep=''),sep='')))
          }
          out <- regsubsets(mod,nvmax=4,data=newdata)
          outsummary <- summary(out)
          
          
          outmat <- outsummary$outmat
          outmat <- getoutmat(outmat)
          
          # generate AIC comparison
          newAIC <- as.numeric(sub("x","",sort(names(chitestAIC(outmat,y,x1,x2,x3,x4)$coeff[-1]))))
          storeAIC <- modelcombineAIC(outmat,newAIC)
          
          
          # generate BIC comparison
          newBIC <- as.numeric(sub("x","",sort(names(chitestBIC(outmat,y,x1,x2,x3,x4)$coeff[-1]))))
          storeBIC <- modelcombineBIC(outmat,newBIC)
          
          
          # Test for nested models
          
          if (nest(outmat)==TRUE){
            
            # Generate delta R-square F test
            newdeltaR2 <- as.numeric(sub("x","",sort(names(FdeltaR2(outmat,n,y,x1,x2,x3,x4)$coeff[-1]))))
            storedeltaR2 <- modelcombinenewdeltaR2(outmat,newdeltaR2)
            
            # choose the model with the smallesr cp-p
            newcp_p <- as.numeric(sub("x","",sort(names(Fcp_pnest(outmat,outsummary,n,y,x1,x2,x3,x4)$coeff[-1]))))
            storecp_p <- modelcombinecp_p(outmat,newcp_p)
            
            
            # choose the model with the largest adjusted R square
            newadjr2 <- as.numeric(sub("x","",sort(names(Fadjr2nest(outmat,outsummary,n,y,x1,x2,x3,x4)$coeff[-1]))))
            storeadjr2 <- modelcombineadjr2(outmat,newadjr2)
            
            # choose the model with the smallest PRESS
            newpress <- as.numeric(sub("x","",sort(names(Fpressnest(outmat,n,y,x1,x2,x3,x4)$coeff[-1]))))
            storepress <- modelcombinepress(outmat,newpress)
            
          } else {
            
            # choose the model with the smallest cp-p
            newcp_p <- as.numeric(sub("x","",sort(names(Fcp_p(outmat,outsummary,y,x1,x2,x3,x4)$coeff[-1]))))
            storecp_p <- modelcombinecp_p(outmat,newcp_p)
            
            # choose the model with the largest adjusted R square
            newadjr2 <- as.numeric(sub("x","",sort(names(Fadjr2(outmat,outsummary,y,x1,x2,x3,x4)$coeff[-1]))))
            storeadjr2 <- modelcombineadjr2(outmat,newadjr2)
            
            # choose the model with the smallest PRESS
            newpress <- as.numeric(sub("x","",sort(names(Fpress(outmat,y,x1,x2,x3,x4)$coeff[-1]))))
            storepress <- modelcombinepress(outmat,newpress)
            
            storedeltaR2 <- c("NA")
          }
          ctr <- ctr + 1
          out4var[ctr,] <- c(n,s,delta,storeAIC,storeBIC,storedeltaR2,storecp_p,storeadjr2,storepress,nest(outmat))
        }
      }
    }
  }
  return(list(out4var))
}
outwhole <- NULL
for (j in 1:420) {
  for (i in 1:num_cores) {
    newout <- parallel.4var.list[[i]][[1]][(50*j-(50-1)):(50*j),]
    outwhole <- rbind(outwhole,newout)
  }
}

outwhole <- as.data.frame(outwhole,stringsAsFactors=FALSE)
colnames(outwhole) <- c("obs","discounter","delta","AIC","BIC","deltaR2","Cp","Adjr2","PRESS","nest")
outwhole <- add_column(outwhole,int.ctr = rep(c(1:1000),420),.after = "delta")

write.csv(outwhole,"output/out4varfix.csv")

