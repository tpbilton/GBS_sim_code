#####################################################
## Simulation code for Section 4.4.1
## - Simulation depends on a number of external packages
## - File paths will need to be changed
## - File is a copy of File S3 in the supplementary material for Bilton et al. (2018) Genetics, 209(1).
#####################################################

### Code for generating the datasets for the simulations
### Author: Timothy Bilton
### Date:   20/02/17
### Edited: 20/11/17

## set the current working directory
setwd("./")

#### Create directory for simulation data to go into
system("mkdir simData")

#### Note: It is assumed that there are the following subfolders
## LepMap2 - Program code for LepMap2 as downloaded from https://sourceforge.net/projects/lepmap2
##   Note: Some modifications to the code where made, these were
##     - Creating another function to compute map distance in morgans 
##     - Fixed a bug when generating the output file
## CriMap  - Program code for CriMap complied in this folder, where the code was downloaded from https://www.animalgenome.org/tools/share/crimap/
## JoinMap - Folder containing the results from the JoinMap simulations.
##           These simulations were performed separately on a Windows OS using the required software and could not
##           be coded in this script.

############ Load packages that are required:
############ Note: assumed to be pre-installed
library(onemap)
library(GUSMap)    # Needs to be v0.1.1 (https://github.com/tpbilton/GUSMap/releases/tag/v0.1.1)
library(extrafont) # for Arial fonts
library(foreach)
library(doSNOW)
#############################################

############# Additional functions used in the simulations ############################################################

## Function for defining the morgan mapping distance in onemap
## Note: Code is taken from the haldane mapping function in onemap and changed for the margon mapping function
# input variable:
#  - rcmb: A recombination fraction value
morgan <- function(rcmb){
  if (is.numeric(rcmb)) {
    if (rcmb >= 0 && rcmb <= 0.5) 
      return(rcmb)
    else stop("the recombination fraction must be between 0 and 0.5")
  }
  else stop("the argument to 'morgan' function must be ", 
            dQuote("numeric"))
}

## Function for setting the mapping function in onemap
## Note: this is an adjustment to the original function in onemap but is change 
## to allow for the morgan mapping function defined above.
set.map.fun <- function (type = c("kosambi", "haldane", "morgan")) 
{
  type <- match.arg(type, c("kosambi", "haldane", "morgan"))
  if (type == "haldane") 
    assign(".map.fun", "haldane", envir = .onemapEnv)
  else if (type == "morgan")
    assign(".map.fun", "morgan", envir = .onemapEnv)
  else assign(".map.fun", "kosambi", envir = .onemapEnv)
}

## Function for converting the phase from onemap to OPGP format
phaseToOPGP_OM <- function(x){
  ## code from here taken from the onemap function print.sequence()
  link.phases <- matrix(NA, length(x$seq.num), 2)
  link.phases[1, ] <- rep(1, 2)
  for (i in 1:length(x$seq.phases)) {
    switch(EXPR = x$seq.phases[i],
           link.phases[i + 1, ] <- link.phases[i, ] * c(1, 1),
           link.phases[i +  1, ] <- link.phases[i, ] * c(1, -1),
           link.phases[i + 1, ] <- link.phases[i, ] * c(-1, 1),
           link.phases[i + 1, ] <- link.phases[i, ] * c(-1, -1))
  }
  if (class(get(x$data.name, pos = 1)) == "outcross") {
    link.phases <- apply(link.phases, 1, function(x) paste(as.character(x), collapse = "."))
    parents <- matrix("", length(x$seq.num), 4)
    for (i in 1:length(x$seq.num)) 
      parents[i, ] <- return.geno(get(x$data.name, pos = 1)$segr.type[x$seq.num[i]], link.phases[i])
    ## Our code below
    #transpose the the parents and set to baseline
    parents[which(parents == 'a')] <-'A'
    parents[which(parents == 'b')] <- 'B'
    
    parents = t(parents)
    if(parents[1,which(apply(parents[1:2,],2,function(x) !(all(x=='A'))))[1]] == 'B')
      parents[1:2,] <- parents[2:1,]
    if(parents[3,which(apply(parents[3:4,],2,function(x) !(all(x=='A'))))[1]] == 'B')
      parents[3:4,] <- parents[4:3,]
    
    ## Now from the parental haplotypes, determine the OPGPs
    return(GUSMap:::parHapToOPGP(parents))
  }
}

## Function for converting phase information in the form of inheritance vectors
## into to OPGPs
phaseToOPGP <- function(phase, nSnps, config, swap=FALSE){
  parHap <- matrix("A",nrow=4,ncol=nSnps)
  if(swap){
    mat <- substr(phase,1,1)
    pat <- substr(phase,2,2)
  }
  else{
    pat <- substr(phase,1,1)
    mat <- substr(phase,2,2)
  }
  
  ## Sort out any duplicates
  if(any(pat=='d')|any(mat=='u')){
    pat[which(pat == 'd')] <- pat[which(pat == 'd')-1]
    mat[which(mat == 'd')] <- pat[which(mat == 'd')-1]
  }
  
  parHap[cbind(suppressWarnings(as.numeric(pat))+1,1:nSnps)] <- "B"
  parHap[cbind(suppressWarnings(as.numeric(mat))+3,1:nSnps)] <- "B"
  parHap[cbind(c(rep(1,sum(config==5)),rep(2,sum(config==5))),which(config==5))] <- "B"
  parHap[cbind(c(rep(3,sum(config==3)),rep(4,sum(config==3))),which(config==3))] <- "B"
  
  ## make sure the baseline phase is correct
  if(parHap[1,which(apply(parHap[1:2,],2,function(x) !(all(x=='A'))))[1]] == 'B')
    parHap[1:2,] <- parHap[2:1,]
  if(parHap[3,which(apply(parHap[3:4,],2,function(x) !(all(x=='A'))))[1]] == 'B')
    parHap[3:4,] <- parHap[4:3,]
  
  ## Convert the parental haplotypes into phase numbers
  OPGP <- GUSMap:::parHapToOPGP(parHap)
}

####### Functions for computing the r.f. estimates for the simulated data

trim_fn <- function(x) gsub("^\\/+|\\/+$", "" , trimws(x))

## GusMap
simEst_GM <- function(sim,direct,ss_mp=F,sim2=NULL){
  
  simPath = trim_fn(paste0(direct,"/",sim))
  
  if(!ss_mp){
    simPara = dget(paste0(simPath,"_info.txt"))
    nSnps = simPara$nSnps
    NoDS = simPara$NoDS
    config = simPara$config
    
    rf.est1 <- matrix(0,nrow=NoDS,ncol=nSnps-1)
    OPGP.est <- matrix(0,nrow=NoDS,ncol=nSnps) 
    time.est1 <- numeric(NoDS)
    epsilon_est <- numeric(NoDS)
    
    for(i in 1:NoDS){
      pstart <- proc.time()[3]
      ## Read in the genon data and depth data
      depth_Ref <- as.matrix(read.table(file=paste0(simPath,i,"_depth_Ref_SEQ.txt")))
      depth_Alt <- as.matrix(read.table(file=paste0(simPath,i,"_depth_Alt_SEQ.txt")))
      
      OPGP.est[i,] <- infer_OPGP_FS(depth_Ref,depth_Alt,config, epsilon=0.001)
      
      MLE <- rf_est_FS(0.01,epsilon = 0.0001,depth_Ref = list(depth_Ref), depth_Alt = list(depth_Alt),
                       OPGP = list(OPGP.est[i,]), reltol=1e-20)
      
      rf.est1[i,] <- MLE$rf
      epsilon_est[i] <- MLE$epsilon
      
      time.est1[i] <- proc.time()[3] - pstart
    }
    return(list(rf1=rf.est1,time1=time.est1,OPGP=OPGP.est,epsilon_est=epsilon_est))
  }
  else{
    simPath2 = trim_fn(paste0(direct,"/",sim2))
    
    simPara_1 = dget(paste0(simPath,"_info.txt"))
    nSnps = simPara_1$nSnps
    NoDS = simPara_1$NoDS
    config_1 = simPara_1$config
    simPara_2 = dget(paste0(simPath2,"_info.txt"))
    config_2 = simPara_2$config
    OPGP <- list(simPara_1$OPGP, simPara_2$OPGP)
    ps = sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:8)))))[-1] - 1
    ms = sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,9:12))))))[-1] - 1
    
    rf.est1_f1 <- rf.est1_m1 <- matrix(0,nrow=NoDS,ncol=nSnps-1)   
    OPGP.est_1 <- OPGP.est_2 <- matrix(13,nrow=NoDS,ncol=nSnps) 
    time.est1 <- numeric(NoDS)
    epsilon_est <- numeric(NoDS)
    
    for(i in 1:NoDS){
      pstart <- proc.time()[3]
      ## Read in the genon data and depth data
      depth_Ref <- list(as.matrix(read.table(file=paste0(simPath,i,"_depth_Ref_SEQ.txt"))),
                        as.matrix(read.table(file=paste0(simPath2,i,"_depth_Ref_SEQ.txt"))))
      depth_Alt <- list(as.matrix(read.table(file=paste0(simPath,i,"_depth_Alt_SEQ.txt"))),
                        as.matrix(read.table(file=paste0(simPath2,i,"_depth_Alt_SEQ.txt"))))
      
      ### Infer the phase
      OPGP.est_1[i,] <- infer_OPGP_FS(depth_Ref[[1]][,which(!(config_1 %in% 6:9))], depth_Alt[[1]][,which(!(config_1 %in% 6:9))], epsilon=0.0001)
      OPGP.est_2[i,] <- infer_OPGP_FS(depth_Ref[[2]][,which(!(config_2 %in% 6:9))], depth_Alt[[2]][,which(!(config_2 %in% 6:9))], epsilon=0.0001)
      
      ## Estimate the parameters        
      rf.est1 <- rf_est_FS(0.01,epsilon = 0.0001,depth_Ref = depth_Ref, depth_Alt = depth_Alt,
                           OPGP = list(OPGP.est_1[i,],OPGP.est_2[i,]), noFam=2, sexSpec = T, reltol=1e-20)
      rf.est1_f1[i,ps] <- rf.est1$rf_p
      rf.est1_m1[i,ms] <- rf.est1$rf_m
      epsilon_est[i] <- rf.est1$epsilon
      
      time.est1[i] <- proc.time()[3] - pstart
    }
  }
  return(list(rf_f1=rf.est1_f1,rf_m1=rf.est1_m1,
              time1=time.est1,OPGP1=OPGP.est_1,OPGP2=OPGP.est_2,epsilon_est=epsilon_est))
}

## Onemap
simEst_OM <- function(sim,direct){
  
  simPath = trim_fn(paste0(direct,"/",sim))
  simPara = dget(paste0(simPath,"_info.txt"))
  nSnps = simPara$nSnps
  NoDS = simPara$NoDS
  config = simPara$config
  
  rf.est.OM <- matrix(nrow=NoDS,ncol=nSnps-1)   
  OPGP.est.OM <- matrix(nrow=NoDS,ncol=nSnps) 
  time.est.OM <- numeric(NoDS)
  
  set.map.fun('morgan')
  for(i in 1:NoDS){
    
    tstart <- proc.time()[3]
    
    ## load in the data
    OM_simData <<- read.outcross(file=paste0(simPath,i,"_OneMap.txt"))
    
    ## Compute the 2 point estimates
    OM_2p <<- rf.2pts(OM_simData,LOD=0)
    
    ## Create the map
    OM_temp <- make.seq(OM_2p,'all')
    OM_temp <- map(OM_temp,1:12)
    rf.est.OM[i,] <- OM_temp$seq.rf
    OPGP.est.OM[i,] <- phaseToOPGP_OM(OM_temp)
    
    time.est.OM[i] <- proc.time()[3] - tstart
  }
  ## Clean up
  rm(OM_simData, OM_2p, pos=".GlobalEnv")
  
  return(list(rf=rf.est.OM,time=time.est.OM,OPGP=OPGP.est.OM))
}

## LepMap
simEst_LM <- function(sim,direct,simPara, ss_mp=F, sim2=NULL){
  
  simPath = trim_fn(paste0(direct,"/",sim))
  
  if(!ss_mp){
    simPara = dget(paste0(simPath,"_info.txt"))
    nSnps = simPara$nSnps
    NoDS = simPara$NoDS
    config = simPara$config
    
    rf.est.LM <- matrix(nrow=NoDS,ncol=nSnps-1)   
    OPGP.est.LM <- matrix(nrow=NoDS,ncol=nSnps)   
    rf.est.LM_err <- matrix(nrow=NoDS,ncol=nSnps-1)   
    OPGP.est.LM_err <- matrix(nrow=NoDS,ncol=nSnps)   
    time.est.LM <- time.est.LM_err <- numeric(NoDS)
    
    for(i in 1:NoDS){
      
      #### approach 1: No error parameter
      tstart <- proc.time()[3]
      
      ## Excute the relevant LepMap2 java command in linux
      system(paste0("java -cp LepMap2/bin/ OrderMarkers evaluateOrder=fixOrder.txt data=",simPath,
                    i,"_LepMap.txt improveOrder=0 sexAveraged=1 useKosambi=2 removeDuplicates=0 learnErrorParameters=0 initError=0 >temp_LepMap.txt"))
      
      ## Read in the data and extract the r.f's and phase
      LP.res <- read.table(file="temp_LepMap.txt",sep="\t",quote="#",skip=3,row.names = 1,stringsAsFactors=F)
      rf.est.LM[i,] <- diff(LP.res[,1])
      OPGP.est.LM[i,] <- phaseToOPGP(LP.res[,4],nSnps,config)
      
      time.est.LM[i] <- proc.time()[3] - tstart
      
      #### approach 2: With error parameter
      tstart <- proc.time()[3]
      
      ## Excute the relevant LepMap2 java command in linux
      system(paste0("java -cp LepMap2/bin/ OrderMarkers evaluateOrder=fixOrder.txt data=",simPath,
                    i,"_LepMap.txt improveOrder=0 sexAveraged=1 useKosambi=2 removeDuplicates=0 >temp_LepMap.txt"))
      
      ## Read in the data and extract the r.f's and phase
      LP.res <- read.table(file="temp_LepMap.txt",sep="\t",quote="#",skip=3,row.names = 1,stringsAsFactors=F)
      rf.est.LM_err[i,] <- diff(LP.res[,1])
      OPGP.est.LM_err[i,] <- phaseToOPGP(LP.res[,4],nSnps,config)
      
      time.est.LM_err[i] <- proc.time()[3] - tstart
    } 
    return(list(list(rf=rf.est.LM,time=time.est.LM,OPGP=OPGP.est.LM),
                list(rf=rf.est.LM_err,time=time.est.LM_err,OPGP=OPGP.est.LM_err)))
  }
  else{
    simPara = dget(paste0(simPath,"_info.txt"))
    nSnps = simPara$nSnps
    NoDS = simPara$NoDS
    
    rf.est.LM_f <- rf.est.LM_m <- rf.est.LM_err_f <- rf.est.LM_err_m <- matrix(nrow=NoDS,ncol=nSnps-1)  
    OPGP.est.LM_f1 <- OPGP.est.LM_f2 <-  matrix(nrow=NoDS,ncol=nSnps)   
    OPGP.est.LM_err_f1 <- OPGP.est.LM_err_f2 <-  matrix(9,nrow=NoDS,ncol=nSnps)   
    time.est.LM <- time.est.LM_err <- numeric(NoDS)
    
    simPath2 = trim_fn(paste0(direct,"/",sim2))
    config_1 = simPara$config
    config_2 = dget(paste0(simPath2,"_info.txt"))$config
    
    for(i in 1:NoDS){
      
      #### approach 1: No error parameter
      tstart <- proc.time()[3]
      
      ## Excute the relevant LepMap2 java command in linux
      system(paste0("java -cp LepMap2/bin/ OrderMarkers evaluateOrder=fixOrder.txt data=",simPath,
                    i,"_LepMap.txt ",simPath2,i,"_LepMap.txt improveOrder=0 useKosambi=2 removeDuplicates=0 learnErrorParameters=0 initError=0 >temp_LepMap.txt"))
      
      ## Read in the data and extract the r.f's and phase
      LP.res <- read.table(file="temp_LepMap.txt",sep="\t",quote="#",skip=3,row.names = 1,stringsAsFactors=F)
      rf.est.LM_f[i,] <- diff(LP.res[,1])
      rf.est.LM_m[i,] <- diff(LP.res[,2])
      OPGP.est.LM_f1[i,] <- phaseToOPGP(LP.res[,4],nSnps,config_1)
      OPGP.est.LM_f2[i,] <- phaseToOPGP(LP.res[,5],nSnps,config_2)
      
      time.est.LM[i] <- proc.time()[3] - tstart
      
      #### approach 2: With error parameter
      tstart <- proc.time()[3]
      
      ## Excute the relevant LepMap2 java command in linux
      system(paste0("java -cp LepMap2/bin/ OrderMarkers evaluateOrder=fixOrder.txt data=",simPath,
                    i,"_LepMap.txt ",simPath2,i,"_LepMap.txt improveOrder=0 useKosambi=2 removeDuplicates=0 >temp_LepMap.txt"))
      
      ## Read in the data and extract the r.f's and phase
      LP.res <- read.table(file="temp_LepMap.txt",sep="\t",quote="#",skip=3,row.names = 1,stringsAsFactors=F)
      rf.est.LM_err_f[i,] <- diff(LP.res[,1])
      rf.est.LM_err_m[i,] <- diff(LP.res[,2])
      OPGP.est.LM_err_f1[i,] <- phaseToOPGP(LP.res[,4],nSnps,config_1)
      OPGP.est.LM_err_f2[i,] <- phaseToOPGP(LP.res[,5],nSnps,config_2)
      
      time.est.LM_err[i] <- proc.time()[3] - tstart
    }
    return(list(list(rf_f=rf.est.LM_f,rf_m=rf.est.LM_m,time=time.est.LM,OPGP1=OPGP.est.LM_f1,OPGP2=OPGP.est.LM_f2),
                list(rf_f=rf.est.LM_err_f,rf_m=rf.est.LM_err_m,time=time.est.LM_err,OPGP1=OPGP.est.LM_err_f1,OPGP2=OPGP.est.LM_err_f2)))  
  }
}

#### CriMap
simEst_CM <- function(sim,direct,ss_mp=F, simPath2=NULL){
  
  if(!ss_mp){
    
    simPath = trim_fn(paste0(direct,"/",sim))
    simPara = dget(paste0(simPath,"_info.txt"))
    nSnps = simPara$nSnps
    NoDS = simPara$NoDS
    config = simPara$config
    
    rf.est.CM <- matrix(nrow=NoDS,ncol=nSnps-1)
    OPGP.est.CM <- matrix(NA,nrow=NoDS,ncol=nSnps)
    time.est.CM <- numeric(NoDS)
    
    setwd("./simData")
    for(i in 1:NoDS){
      
      tstart <- proc.time()[3]
      
      ## Run crimap and produce the output required
      system(paste0("../CriMap/bin/crimap 1_",trim_fn(sim),i," prepare <cm_run.txt"))
      system(paste0("../CriMap/bin/crimap 1_",trim_fn(sim),i," fixed >cm_out.txt"))
      
      ## load output in R and extract the results
      out <- readLines("cm_out.txt")
      findx <- which(substr(out,1,6) == "Sex_av") + 3
      rf.est.CM[i,] <- as.numeric(sapply(seq(findx,findx+(nSnps-1)*2-1,2),function(x){
        return(strsplit(out[x],"\\s+")[[1]][2])
      }))
      
      time.est.CM[i] <- proc.time()[3] - tstart
      ## Delete the created files
      system(paste0("rm chr1_",trim_fn(sim),i,".dat"))
      system(paste0("rm chr1_",trim_fn(sim),i,".loc"))
      system(paste0("rm chr1_",trim_fn(sim),i,".par"))
    } 
    setwd("../")
    return(list(rf=rf.est.CM,time=time.est.CM,OPGP=OPGP.est.CM))
  }
  else{
    
    simPara = dget(paste0(trim_fn(paste0(direct,"/",simPath2)),"_info.txt"))
    nSnps = simPara$nSnps
    NoDS = simPara$NoDS
    
    rf.est.CM_f <- rf.est.CM_m <- matrix(nrow=NoDS,ncol=nSnps-1)
    OPGP.est.CM_f1 <- OPGP.est.CM_f2 <- matrix(NA,nrow=NoDS,ncol=nSnps)
    time.est.CM <- numeric(NoDS)
    
    setwd("./simData")
    for(i in 1:NoDS){
      
      tstart <- proc.time()[3]
      
      ## Run crimap and produce the output required
      system(paste0("../CriMap/bin/crimap 1_",trim_fn(sim),i," prepare <cm_run2.txt"))
      system(paste0("../CriMap/bin/crimap 1_",trim_fn(sim),i," fixed >cm_out.txt"))
      
      out <- readLines("cm_out.txt")
      findx <- which(substr(out,1,8) == "Sex-spec") + 3
      rf.est.CM_f[i,] <- as.numeric(sapply(seq(findx,findx+(nSnps-1)*2-1,2),function(x){
        return(strsplit(out[x],"\\s+")[[1]][4])
      }))
      rf.est.CM_m[i,] <- as.numeric(sapply(seq(findx,findx+(nSnps-1)*2-1,2),function(x){
        return(strsplit(out[x],"\\s+")[[1]][2])
      }))
      
      time.est.CM[i] <- proc.time()[3] - tstart
      ## Delete the created files
      system(paste0("rm chr1_",trim_fn(sim),i,".dat"))
      system(paste0("rm chr1_",trim_fn(sim),i,".loc"))
      system(paste0("rm chr1_",trim_fn(sim),i,".par"))
    }
    setwd("../")
    return(list(rf_f=rf.est.CM_f,rf_m=rf.est.CM_m,time=time.est.CM,OPGP1=OPGP.est.CM_f1,OPGP1=OPGP.est.CM_f2))
  }
}

#### JoinMap
## It is assumed that the results have already been generated using some custom scripts via 
## Visual studio and coded UI tests and the results were written to files.
simEst_JM <- function(sim, direct){
  
  simPath = paste0(trim_fn(direct),"/",trim_fn(sim))
  simPara = dget(paste0(simPath,"_ds_info.txt"))
  nSnps = simPara$nSnps
  NoDS = simPara$NoDS
  config = simPara$config
  
  rf.est.JM <- matrix(nrow=NoDS,ncol=nSnps-1)
  OPGP.est.JM <- matrix(nrow=NoDS,ncol=nSnps)
  time.est.JM <- scan(paste0(simPath,"_time_short.loc"))/1000
  
  for(i in 1:NoDS){
    
    ## Specify file names
    resultsFile <- paste0(simPath,"_results",i,"_JoinMap.txt")
    phaseFile <- paste0(simPath,"_phase",i,"_JoinMap.loc")
    
    ## Extract the r.f.estimates
    rf.est.JM[i,] <- inv.haldane(sapply(which(readLines(resultsFile)=="best map order:")+5+seq(0,10,1),function(x){
      tempV <- strsplit(readLines(resultsFile)[x],"\\s+")[[1]]
      return(as.numeric(tempV[length(tempV)-1]))
    }))
    
    ## Extract the phase and compute the OPGPs
    phase <- sapply(8+0:11,function(x){
      return(substr(strsplit(readLines(phaseFile)[x],"\\s+")[[1]][3],2,3))
    })
    OPGP.est.JM[i,] <- phaseToOPGP(phase,nSnps=nSnps,config=config, swap=TRUE)
  }
  return(list(rf=rf.est.JM,time=time.est.JM,OPGP=OPGP.est.JM))
}

### Function for plotting the recombination fractions
plotRFs <- function(results, nSnps, myCol, names, rf, plotName, parent=0, ext="pdf",
                    leg.insert=-0.2){
  
  nMethods <- length(results)
  
  # Compute the expectation of the parameters (each method is represented by a column)
  exp.est <- matrix(unlist(lapply(results,colMeans)),ncol=nSnps-1, nrow=nMethods, byrow=T)
  
  quan.est <- lapply(results, function(x) {
    sapply(1:(nSnps-1), function(y) quantile(x[,y],prob=c(0.025,0.25,0.5,0.75,0.975)) ) })
  quan.05 <- matrix(unlist(lapply(quan.est,function(x) x[1,]),use.names=F),ncol=nSnps-1,byrow=T)
  quan.25 <- matrix(unlist(lapply(quan.est,function(x) x[2,]),use.names=F),ncol=nSnps-1,byrow=T)
  quan.50 <- matrix(unlist(lapply(quan.est,function(x) x[3,]),use.names=F),ncol=nSnps-1,byrow=T)
  quan.75 <- matrix(unlist(lapply(quan.est,function(x) x[4,]),use.names=F),ncol=nSnps-1,byrow=T)
  quan.95 <- matrix(unlist(lapply(quan.est,function(x) x[5,]),use.names=F),ncol=nSnps-1,byrow=T)
  
  tiff(paste0(plotName,".tiff"),width=3375,height=2475, res=450, pointsize=12, family="Arial",type='cairo',compression='lzw')
  par(mfrow=c(1,1),mar=c(5.8,4.6,1.1,1.1))
  
  xCoor <- 1:nMethods + as.vector(sapply(seq(0,(nMethods+2)*10,nMethods+2),function(x) rep(x,nMethods)))
  xCoorlab <- mean(c(1,nMethods)) + seq(0,(nMethods+2)*10,nMethods+2)
  plot(NULL,NULL,xaxt='n',xlim=c(1,max(xCoor)), xlab="",
       ylab = switch(parent+1, expression(italic(hat(r)[j])), expression(italic(hat(r)[j]^p)), expression(italic(hat(r)[j]^m))),
       ylim = c(min(unlist(lapply(quan.est,function(z) min(z[,1:(nSnps-1)])))),
                max(unlist(lapply(quan.est,function(z) max(z[,1:(nSnps-1)]))))))
  abline(h=rf,lty=2)
  if(parent %in% c(1,2))
    abline(h=2*rf-2*rf^2,lty=3)
  
  points(xCoor,as.vector(exp.est[,1:(nSnps-1)]),col=myCol,pch=19)
  # end lines
  segments(xCoor-0.4, as.vector(quan.05[,1:(nSnps-1)]), xCoor+0.4, as.vector(quan.05[,1:(nSnps-1)]), col=myCol)
  segments(xCoor-0.4, as.vector(quan.25[,1:(nSnps-1)]), xCoor+0.4, as.vector(quan.25[,1:(nSnps-1)]), col=myCol)
  segments(xCoor-0.4, as.vector(quan.50[,1:(nSnps-1)]), xCoor+0.4, as.vector(quan.50[,1:(nSnps-1)]), col=myCol)
  segments(xCoor-0.4, as.vector(quan.75[,1:(nSnps-1)]), xCoor+0.4, as.vector(quan.75[,1:(nSnps-1)]), col=myCol)
  segments(xCoor-0.4, as.vector(quan.95[,1:(nSnps-1)]), xCoor+0.4, as.vector(quan.95[,1:(nSnps-1)]), col=myCol)
  # segment lines
  segments(xCoor, as.vector(quan.05[,1:(nSnps-1)]), xCoor, as.vector(quan.25[,1:(nSnps-1)]), col=myCol, lty=2)
  segments(xCoor, as.vector(quan.25[,1:(nSnps-1)]), xCoor, as.vector(quan.75[,1:(nSnps-1)]), col=myCol)
  segments(xCoor, as.vector(quan.75[,1:(nSnps-1)]), xCoor, as.vector(quan.95[,1:(nSnps-1)]), col=myCol, lty=2)
  axis(side=1,at=c(xCoorlab),lwd.ticks=0,
       labels=c(1:(nSnps-1)))
  mtext(expression(italic(j)), side=1, line=2.5)
  legend('bottom', legend=levels(names),pch=19,col=myCol,horiz=TRUE,lty=1,inset = c(0,leg.insert), xpd=T, bty='o')
  
  dev.off()
}

# Function for plotting the Haldane map distances of the overall map
plotMAPDIST <- function(results, myCol, names, rf, parent=0, depth, epsilon,ctex=1){
  
  
  nMethods <- length(myCol)
  if(length(levels(names)) != nMethods || any(unlist(lapply(results,function(x) length(x) != nMethods))))
    stop("Number of methods do not match up for all the different objects")
  
  ncols = 3;
  nrows = ceiling(length(results) %/% ncols)
  
  ## Compute the map distances
  sim_mapDist <- lapply(results,function(List){
    lapply(List,function(x) apply(x,1,function(y) sum(haldane(y))))
  })
  
  ## recode infinity estimates of map distance
  for(i in 1:length(sim_mapDist)){
    for(j in 1:nMethods){
      sim_mapDist[[i]][[j]][which(is.infinite(sim_mapDist[[i]][[j]]))] <- 1836.84
    }
  }  
  
  exp.est <- matrix(unlist(lapply(sim_mapDist, function(x) lapply(x,mean))),nrow=9,ncol=nMethods,byrow=T)
  
  quan.dist <- lapply(sim_mapDist, function(x){
    matrix(unlist(lapply(x, function(y) quantile(y,prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=T))), nrow=5,ncol=nMethods,byrow=F)
  })
  
  par(mfrow=c(nrows,ncols),mar=c(0,0,0,0), oma=c(4.5,6,4.5,5))
  
  #ylim_vec <- sapply(0:(nrows-1),function(x) c(min(unlist(quan.dist[1:ncols+ncols*x])),max(unlist(quan.dist[1:ncols+ncols*x]))))
  ylim_vec <- sapply(0:(nrows-1),function(x) c(0,max(unlist(quan.dist[1:ncols+ncols*x]))))
  
  for(sim in 1:(ncols*nrows)){
    if(sim > length(results))
      plot(names, rep(0,nMethods), ylim=ylim_vec[2,ceiling(sim/ncols)], las=2, ylab="",
           xaxt='n', yaxt=ifelse((sim %% ncols) == 1, 's', 'n'),xpd=F)
    else{
      x = 1:nMethods; y = exp.est[sim,]; z = quan.dist[[sim]]
      plot(names, rep(0,nMethods), ylim=ylim_vec[,ceiling(sim/ncols)], las=2, ylab="",
           xaxt='n', yaxt=ifelse((sim %% ncols) == 1, 's', 'n'),xpd=F,cex.axis=ctex)
      abline(h=sum(haldane(rep(rf,11))),lty=2)
      points(x,y,pch=19,col=myCol)
      # end lines
      segments(x-0.15,z[1,],x+0.15,z[1,], col=myCol)
      segments(x-0.15,z[2,],x+0.15,z[2,], col=myCol)
      segments(x-0.15,z[3,],x+0.15,z[3,], col=myCol)
      segments(x-0.15,z[4,],x+0.15,z[4,], col=myCol)
      segments(x-0.15,z[5,],x+0.15,z[5,], col=myCol)
      # segment lines
      segments(x,z[1,],x,z[2,], col=myCol, lty=2)
      segments(x,z[2,],x,z[4,], col=myCol)
      segments(x,z[4,],x,z[5,], col=myCol, lty=2)
    }
    if(sim %in% (ncols*nrows - 0:2))
      axis(1,at=c(1:nMethods), labels=levels(names), las=2, cex.axis=ctex)
    if(sim %in% c(1:ncols)){
      mtext(bquote(epsilon==.(epsilon[sim])),side=3,padj=-0.5)
    }
    if(sim==median(1:ncols))
      mtext("Sequencing Error", side=3,padj=-3)
    if(sim==median(seq(0,nrows*ncols-1,ncols)+1))
      mtext("Map Distance (cM)",side=2,padj=-4)
    if((sim %% ncols)==0)
      mtext(bquote(italic(mu[d[j]])==.(depth[sim %/% ncols])),side=4,padj=1)
    if(sim==median(seq(ncols,nrows*ncols,nrows)))
      mtext("Mean Read Depth", side=4,padj=4)
  }
}

## Function for plotting the run times
plotTIME <- function(results, myCol, names, NoDS){
  
  nMethods <- length(results)
  
  exp.est <- log(unlist(lapply(results, mean)))
  quan.time <- t(log(matrix(unlist(lapply(results, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))), byrow=T, ncol =5)))
  
  x = 1:nMethods
  plot(names, rep(0,nMethods), ylim=c(min(quan.time),max(quan.time)), las=2,
       ylab="Run time (log(s))")
  points(x,exp.est,pch=19,col=myCol)
  # end lines
  segments(x-0.15,quan.time[1,],x+0.15,quan.time[1,], col=myCol)
  segments(x-0.15,quan.time[2,],x+0.15,quan.time[2,], col=myCol)
  segments(x-0.15,quan.time[3,],x+0.15,quan.time[3,], col=myCol)
  segments(x-0.15,quan.time[4,],x+0.15,quan.time[4,], col=myCol)
  segments(x-0.15,quan.time[5,],x+0.15,quan.time[5,], col=myCol)
  # segment lines
  segments(x,quan.time[1,],x,quan.time[2,], col=myCol, lty=2)
  segments(x,quan.time[2,],x,quan.time[4,], col=myCol)
  segments(x,quan.time[4,],x,quan.time[5,], col=myCol, lty=2)
  
}

plotERROR <- function(results, trueEp, meanDepth, myCol=AgRgreen,cex){
  
  nMethods = length(results)
  
  exp.est <- unlist(lapply(results, mean))
  quan.time <- t(matrix(unlist(lapply(results, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))), byrow=T, ncol =5))
  
  par(mar=c(5,5,2,2),cex=cex)
  
  x = 1:nMethods
  plot(NULL, xlim=c(1,nMethods), ylim=c(min(quan.time),max(quan.time)), las=2,
       ylab="",xaxt='n', xlab="")
  title(ylab=expression("Sequencing Error Estimate, " ~ hat(epsilon)), line=3.5)
  axis(side=1,at=c(1:9),labels=rep(meanDepth,3), las=1)
  title(xlab=expression("Mean Depth," ~ italic(mu[d[j]])), line=3)
  mtext(side=3, at=c(2,5,8), paste0("\u03b5=",trueEp), cex=cex)
  abline(v=3.5)
  abline(v=6.5)
  lines(x=c(-0.5,3.5),y=rep(trueEp[1],2),lty=3)
  lines(x=c(3.5,6.5),y=rep(trueEp[2],2),lty=3)
  lines(x=c(6.5,9.5),y=rep(trueEp[3],2),lty=3)
  #lines(x = c(3.5,3.5),y=c(-0.011,max(quan.time)+0.0006), lty=2,xpd=NA)
  #lines(x = c(6.5,6.5),y=c(-0.011,max(quan.time)+0.0006), lty=2,xpd=NA)
  points(x,exp.est,pch=19,col=myCol)
  # end lines
  segments(x-0.15,quan.time[1,],x+0.15,quan.time[1,], col=myCol)
  segments(x-0.15,quan.time[2,],x+0.15,quan.time[2,], col=myCol)
  segments(x-0.15,quan.time[3,],x+0.15,quan.time[3,], col=myCol)
  segments(x-0.15,quan.time[4,],x+0.15,quan.time[4,], col=myCol)
  segments(x-0.15,quan.time[5,],x+0.15,quan.time[5,], col=myCol)
  # segment lines
  segments(x,quan.time[1,],x,quan.time[2,], col=myCol, lty=2)
  segments(x,quan.time[2,],x,quan.time[4,], col=myCol)
  segments(x,quan.time[4,],x,quan.time[5,], col=myCol, lty=2)
}

### Function for extracting the proportion of the data sets which correctly inferred the OPGP.
correctOPGP <- function(results,OPGP){
  nMethods <- length(results)
  return( sapply(1:nMethods, function(x) 100 * mean(apply(results[[x]], 1, function(y) all(y == OPGP)))) )
}

## Inverse Haldane mapping function
inv.haldane <- function(d){
  return((1-exp(-2*d/100))/2)
}

stop("end source here")


###################################################################################################################
#### Code for running the simulations

##########################################
############ Simulation 1: ###############

## Specify the simulation parameters
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=20,filename="Sim1_ds", direct="simData", thres=11,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)
nSnps <- dget("simData/Sim1_ds_info.txt")$nSnps
### Create the file for LepMAP to specific a fixed order
write(1:nSnps,file = "simData/fixOrder.txt", ncolumns = 1)
### Create the txt file required to automate CRI-MAP
cat(paste0(rep("n",4),sep="\n"),"4\n",paste(0:(nSnps-1),collapse=" ")," *\n","y", file="simData/cm_run.txt", sep="")

## Estimate the rf using the various methods
sim1_GM <- simEst_GM("Sim1_ds","simData")
sim1_OM <- simEst_OM("Sim1_ds","simData")
LPres <- simEst_LM("Sim1_ds","simData")
sim1_LM <- LPres[[1]]
sim1_LM_err <- LPres[[2]]
sim1_CM <- simEst_CM("Sim1_ds","simData")
sim1_JM <- simEst_JM("Sim1","JoinMap/Results")

#save(sim1_GM,sim1_LM,sim1_LM_err,sim1_OM,sim1_CM,sim1_JM,file="Sim1_res_Final.RData")

##########################################
############ Simulation 2: ###############
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0.002
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=20,filename="Sim2_ds", direct="simData", thres=11,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)

## estimate the r.f.'s
sim2_GM <- simEst_GM("Sim2_ds","simData")
sim2_OM <- simEst_OM("Sim2_ds","simData")
LPres <- simEst_LM("Sim2_ds","simData")
sim2_LM <- LPres[[1]]
sim2_LM_err <- LPres[[2]]
sim2_CM <- simEst_CM("Sim2_ds","simData")
sim2_JM <- simEst_JM("Sim2","JoinMap/Results")

#save(sim2_GM,sim2_LM,sim2_LM_err,sim2_OM,sim2_CM,sim2_JM,file="Sim2_res_Final.RData")

##########################################
############ Simulation 3: ###############
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0.01
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=20,filename="Sim3_ds", direct="simData", thres=11,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)

## estimate the r.f.'s
sim3_GM <- simEst_GM("Sim3_ds","simData")
sim3_OM <- simEst_OM("Sim3_ds","simData")
LPres <- simEst_LM("Sim3_ds","simData")
sim3_LM <- LPres[[1]]
sim3_LM_err <- LPres[[2]]
sim3_CM <- simEst_CM("Sim3_ds","simData")
sim3_JM <- simEst_JM("Sim3","JoinMap/Results")

#save(sim3_GM,sim3_LM,sim3_LM_err,sim3_OM,sim3_CM,sim3_JM,file="Sim3_res_Final.RData")


##########################################
############ Simulation 4: ###############

## Specify the simulation parameters
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=2,filename="Sim4_ds", direct="simData", thres=NULL,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)
nSnps = dget("simData/Sim1_ds_info.txt")$nSnps
### Create the file for LepMAP to specific a fixed order
write(1:nSnps,file = "simData/fixOrder.txt", ncolumns = 1)
### Create the txt file required to automate CRI-MAP
cat(paste0(rep("n",4),sep="\n"),"4\n",paste(0:(nSnps-1),collapse=" ")," *\n","y", file="simData/cm_run.txt", sep="")

## Estimate the rf using the various methods
sim4_GM <- simEst_GM("Sim4_ds","simData")
sim4_OM <- simEst_OM("Sim4_ds","simData")
LPres <- simEst_LM("Sim4_ds","simData")
sim4_LM <- LPres[[1]]
sim4_LM_err <- LPres[[2]]
sim4_CM <- simEst_CM("Sim4_ds","simData")
sim4_JM <- simEst_JM("Sim4","JoinMap/Results")

#save(sim4_GM,sim4_LM,sim4_LM_err,sim4_OM,sim4_CM,sim4_JM,file="Sim4_res_Final.RData")

##########################################
############ Simulation 5: ###############

## Specify the simulation parameters
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0.002
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=2,filename="Sim5_ds", direct="simData", thres=NULL,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)
nSnps = dget("simData/Sim5_ds_info.txt")$nSnps
### Create the file for LepMAP to specific a fixed order
write(1:nSnps,file = "simData/fixOrder.txt", ncolumns = 1)
### Create the txt file required to automate CRI-MAP
cat(paste0(rep("n",4),sep="\n"),"4\n",paste(0:(nSnps-1),collapse=" ")," *\n","y", file="simData/cm_run.txt", sep="")

## Estimate the rf using the various methods
sim5_GM <- simEst_GM("Sim5_ds","simData")
sim5_OM <- simEst_OM("Sim5_ds","simData")
LPres <- simEst_LM("Sim5_ds","simData")
sim5_LM <- LPres[[1]]
sim5_LM_err <- LPres[[2]]
sim5_CM <- simEst_CM("Sim5_ds","simData")
sim5_JM <- simEst_JM("Sim5","JoinMap/Results")

#save(sim5_GM,sim5_LM,sim5_LM_err,sim5_OM,sim5_CM,sim5_JM,file="Sim5_res_Final.RData")


##########################################
############ Simulation 6: ###############

## Specify the simulation parameters
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0.01
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=2,filename="Sim6_ds", direct="simData", thres=NULL,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)
nSnps = dget("simData/Sim6_ds_info.txt")$nSnps
### Create the file for LepMAP to specific a fixed order
write(1:nSnps,file = "simData/fixOrder.txt", ncolumns = 1)
### Create the txt file required to automate CRI-MAP
cat(paste0(rep("n",4),sep="\n"),"4\n",paste(0:(nSnps-1),collapse=" ")," *\n","y", file="simData/cm_run.txt", sep="")

## Estimate the rf using the various methods
sim6_GM <- simEst_GM("Sim6_ds","simData")
sim6_OM <- simEst_OM("Sim6_ds","simData")
LPres <- simEst_LM("Sim6_ds","simData")
sim6_LM <- LPres[[1]]
sim6_LM_err <- LPres[[2]]
sim6_CM <- simEst_CM("Sim6_ds","simData")
sim6_JM <- simEst_JM("Sim6","JoinMap/Results")

#save(sim6_GM,sim6_LM,sim6_LM_err,sim6_OM,sim6_CM,sim6_JM,file="Sim6_res_Final.RData")

##########################################
############ Simulation 7: ###############

## Specify the simulation parameters
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=10,filename="Sim7_ds", direct="simData", thres=6,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)
nSnps = dget("simData/Sim7_ds_info.txt")$nSnps
### Create the file for LepMAP to specific a fixed order
write(1:nSnps,file = "simData/fixOrder.txt", ncolumns = 1)
### Create the txt file required to automate CRI-MAP
cat(paste0(rep("n",4),sep="\n"),"4\n",paste(0:(nSnps-1),collapse=" ")," *\n","y", file="simData/cm_run.txt", sep="")

## Estimate the rf using the various methods
sim7_GM <- simEst_GM("Sim7_ds","simData")
sim7_OM <- simEst_OM("Sim7_ds","simData")
LPres <- simEst_LM("Sim7_ds","simData")
sim7_LM <- LPres[[1]]
sim7_LM_err <- LPres[[2]]
sim7_CM <- simEst_CM("Sim7_ds","simData")
sim7_JM <- simEst_JM("Sim7","JoinMap/Results")

#save(sim7_GM,sim7_LM,sim7_LM_err,sim7_OM,sim7_CM,sim7_JM,file="Sim7_res_Final.RData")

##########################################
############ Simulation 8: ###############

## Specify the simulation parameters
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0.002
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=10,filename="Sim8_ds", direct="simData", thres=6,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)
nSnps = dget("simData/Sim8_ds_info.txt")$nSnps
### Create the file for LepMAP to specific a fixed order
write(1:nSnps,file = "simData/fixOrder.txt", ncolumns = 1)
### Create the txt file required to automate CRI-MAP
cat(paste0(rep("n",4),sep="\n"),"4\n",paste(0:(nSnps-1),collapse=" ")," *\n","y", file="simData/cm_run.txt", sep="")

## Estimate the rf using the various methods
sim8_GM <- simEst_GM("Sim8_ds","simData")
sim8_OM <- simEst_OM("Sim8_ds","simData")
LPres <- simEst_LM("Sim8_ds","simData")
sim8_LM <- LPres[[1]]
sim8_LM_err <- LPres[[2]]
sim8_CM <- simEst_CM("Sim8_ds","simData")
sim8_JM <- simEst_JM("Sim8","JoinMap/Results")

#save(sim8_GM,sim8_LM,sim8_LM_err,sim8_OM,sim8_CM,sim8_JM,file="Sim8_res_Final.RData")


##########################################
############ Simulation 9: ###############

## Specify the simulation parameters
nInd <- 100   
NoDS <- 1000
rf <- 0.01; ep = 0.01
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## simulate the data
simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=10,filename="Sim9_ds", direct="simData", thres=6,
      formats=list(gusmap=T,onemap=T,crimap=T,lepmap=T,joinmap=F), NoDS=NoDS, rd_dist = rd_dist,
      seed1=11984321, seed2=91544322)
nSnps = dget("simData/Sim9_ds_info.txt")$nSnps
### Create the file for LepMAP to specific a fixed order
write(1:nSnps,file = "simData/fixOrder.txt", ncolumns = 1)
### Create the txt file required to automate CRI-MAP
cat(paste0(rep("n",4),sep="\n"),"4\n",paste(0:(nSnps-1),collapse=" ")," *\n","y", file="simData/cm_run.txt", sep="")

## Estimate the rf using the various methods
sim9_GM <- simEst_GM("Sim9_ds","simData")
sim9_OM <- simEst_OM("Sim9_ds","simData")
LPres <- simEst_LM("Sim9_ds","simData")
sim9_LM <- LPres[[1]]
sim9_LM_err <- LPres[[2]]
sim9_CM <- simEst_CM("Sim9_ds","simData")
sim9_JM <- simEst_JM("Sim9","JoinMap/Results")

#save(sim9_GM,sim9_LM,sim9_LM_err,sim9_OM,sim9_CM,sim9_JM,file="Sim9_res_Final.RData")

######################################
#### Plots of the simulation results

AgRgreen <- rgb(141,199,63, maxColorValue = 255)
lblue <- rgb(0,150,255, maxColorValue = 255)
gold <- rgb(180,180,0, maxColorValue = 255)

names <- factor(rep(0,6),levels=1:6,label=c("GM", "LM2\u03b5", "LM2", "OM", "CM","JM"))
myCol <- c(AgRgreen,"orange","red",lblue,'purple',gold)

results <- list(
  list(sim1_GM[[1]],sim1_LM_err[[1]],sim1_LM[[1]],sim1_OM[[1]],sim1_CM[[1]],sim1_JM[[1]]),
  list(sim2_GM[[1]],sim2_LM_err[[1]],sim2_LM[[1]],sim2_OM[[1]],sim2_CM[[1]],sim2_JM[[1]]),
  list(sim3_GM[[1]],sim3_LM_err[[1]],sim3_LM[[1]],sim3_OM[[1]],sim3_CM[[1]],sim3_JM[[1]]),
  list(sim7_GM[[1]],sim7_LM_err[[1]],sim7_LM[[1]],sim7_OM[[1]],sim7_CM[[1]],sim7_JM[[1]]),
  list(sim8_GM[[1]],sim8_LM_err[[1]],sim8_LM[[1]],sim8_OM[[1]],sim8_CM[[1]],sim8_JM[[1]]),
  list(sim9_GM[[1]],sim9_LM_err[[1]],sim9_LM[[1]],sim9_OM[[1]],sim9_CM[[1]],sim9_JM[[1]]),
  list(sim4_GM[[1]],sim4_LM_err[[1]],sim4_LM[[1]],sim4_OM[[1]],sim4_CM[[1]],sim4_JM[[1]]),
  list(sim5_GM[[1]],sim5_LM_err[[1]],sim5_LM[[1]],sim5_OM[[1]],sim5_CM[[1]],sim5_JM[[1]]),
  list(sim6_GM[[1]],sim6_LM_err[[1]],sim6_LM[[1]],sim6_OM[[1]],sim6_CM[[1]],sim6_JM[[1]])
)
epsilon = c(0,0.002,0.01)
depths=c(20,10,2)

## plots of Map distance
tiff(paste0("Figure1.tiff"),width=7.5*450,height=8*450, type='cairo', res=450, family="Arial",compression='lzw')
plotMAPDIST(results, myCol,names, 0.01,depth = depths, epsilon=epsilon, ctex=1.5)
dev.off()         

# Plot of error parameter estimates from GM
results_err <- list(sim1_GM[[4]],sim7_GM[[4]],sim4_GM[[4]],sim2_GM[[4]],sim8_GM[[4]],sim5_GM[[4]],sim3_GM[[4]],sim9_GM[[4]],sim6_GM[[4]])
tiff(paste0("Figure2.tiff"),width=3.6*450,height=3.6*450, type='cairo', res=450, family="Arial",compression='lzw')
plotERROR(results_err,epsilon,depths,AgRgreen,cex=0.8)
dev.off()

# Plot the computation time
results.time <- list(
  c(sim1_GM[[2]],sim2_GM[[2]],sim3_GM[[2]],sim4_GM[[2]],sim5_GM[[2]],sim6_GM[[2]],sim7_GM[[2]],sim8_GM[[2]],sim9_GM[[2]]),
  c(sim1_LM_err[[2]],sim2_LM_err[[2]],sim3_LM_err[[2]],sim4_LM_err[[2]],sim5_LM_err[[2]],sim6_LM_err[[2]],sim7_LM_err[[2]],sim8_LM_err[[2]],sim9_LM_err[[2]]),
  c(sim1_LM[[2]],sim2_LM[[2]],sim3_LM[[2]],sim4_LM[[2]],sim5_LM[[2]],sim6_LM[[2]],sim7_LM[[2]],sim8_LM[[2]],sim9_LM[[2]]),
  c(sim1_OM[[2]],sim2_OM[[2]],sim3_OM[[2]],sim4_OM[[2]],sim5_OM[[2]],sim6_OM[[2]],sim7_OM[[2]],sim8_OM[[2]],sim9_OM[[2]]),
  c(sim1_CM[[2]],sim2_CM[[2]],sim3_CM[[2]],sim4_CM[[2]],sim5_CM[[2]],sim6_CM[[2]],sim7_CM[[2]],sim8_CM[[2]],sim9_CM[[2]]),
  c(sim1_JM[[2]],sim2_JM[[2]],sim3_JM[[2]],sim4_JM[[2]],sim5_JM[[2]],sim6_JM[[2]],sim7_JM[[2]],sim8_JM[[2]],sim9_JM[[2]])
)

tiff(paste0("Figure3.tiff"),width=3.6*450,height=3.6*450, type='cairo', res=450, family="Arial",compression='lzw')
par(mar=c(4.1,4.1,1.1,1.1), cex=0.8)
plotTIME(results = results.time,myCol,names,NoDS = 1000*9)
dev.off()

## plot the distribution of the recombination fractions
for(sim in 1:9)
  plotRFs(results[[sim]], 12, myCol, names, 0.01, paste0("FigureS",sim), parent=0,
          leg.insert=-0.28)

## Subset of the map distances for low depth
results_subset <- list(
  list(sim4_GM[[1]],sim4_LM_err[[1]]),
  list(sim5_GM[[1]],sim5_LM_err[[1]]),
  list(sim6_GM[[1]],sim6_LM_err[[1]])
)

tiff(paste0("FigureS10.tiff"),width=7.5*450,height=8*450/3, type='cairo', res=450, family="Arial",compression='lzw')
plotMAPDIST(results_subset, myCol[1:2],factor(rep(0,2),levels=1:2,label=c("GM", "LM2\u03b5")), 0.01,depth = depths[3], epsilon=epsilon, ctex=1.5)
dev.off()


### Output the results of the OPGPs inference
results2 <- list(
  list(sim1_GM[[3]],sim1_LM_err[[3]],sim1_LM[[3]],sim1_OM[[3]],sim1_CM[[3]],sim1_JM[[3]]),
  list(sim2_GM[[3]],sim2_LM_err[[3]],sim2_LM[[3]],sim2_OM[[3]],sim2_CM[[3]],sim2_JM[[3]]),
  list(sim3_GM[[3]],sim3_LM_err[[3]],sim3_LM[[3]],sim3_OM[[3]],sim3_CM[[3]],sim3_JM[[3]]),
  list(sim7_GM[[3]],sim7_LM_err[[3]],sim7_LM[[3]],sim7_OM[[3]],sim7_CM[[3]],sim7_JM[[3]]),
  list(sim8_GM[[3]],sim8_LM_err[[3]],sim8_LM[[3]],sim8_OM[[3]],sim8_CM[[3]],sim8_JM[[3]]),
  list(sim9_GM[[3]],sim9_LM_err[[3]],sim9_LM[[3]],sim9_OM[[3]],sim9_CM[[3]],sim9_JM[[3]]),
  list(sim4_GM[[3]],sim4_LM_err[[3]],sim4_LM[[3]],sim4_OM[[3]],sim4_CM[[3]],sim4_JM[[3]]),
  list(sim5_GM[[3]],sim5_LM_err[[3]],sim5_LM[[3]],sim5_OM[[3]],sim5_CM[[3]],sim5_JM[[3]]),
  list(sim6_GM[[3]],sim6_LM_err[[3]],sim6_LM[[3]],sim6_OM[[3]],sim6_CM[[3]],sim6_JM[[3]])
) 
OPGPres <- matrix(unlist(lapply(results2,function(x) correctOPGP(x,dget("simData/Sim9_ds_info.txt")$OPGP))), nrow=9,ncol=6, byrow=T)



##########################################
############ Second set of simulations
## Look at the difference in depth in the rf estimates for a fixed sequencing effort

## Set the simulation parameters
NoDS <- 1000
nSnps <- 12
rf <- 0.01; ep=0.002
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

## To speed things up, use the foreach package
## Set up the clusters
cl <- makeCluster(4)
registerDoSNOW(cl)

## Undertake the simulation
nInd_V <- c(833,416,277,208,166,138,119,104,92,83,75,69,64,59,55)
depths <- 10000/(12*nInd_V)
rf.est <- foreach(i = iter(1:length(nInd_V)), .packages="GUSMap", .combine='c') %dopar% {
  rf.v <- matrix(nrow=NoDS,ncol=nSnps)
  for(run in 1:NoDS){
    seed = 106935 + 5*run
    newData <- simFS(0.01,config=config,nInd=nInd_V[i],epsilon=ep,meanDepth=depths[i],
                     formats=NULL, NoDS=1, rd_dist = rd_dist, seed1=11984321, seed2=seed)
    orig_est <- rf_est_FS(0.01,epsilon=0.0001,depth_Ref = list(newData$depth_Ref), depth_Alt = list(newData$depth_Alt),OPGP = list(newData$OPGP), maxit=1000,reltol=1e-20)
    rf.v[run,] <- c(orig_est$rf,orig_est$epsilon)
  }
  return(list(rf.v))
}

stopCluster(cl)

## Compute the mean square error
MSE <- unlist(lapply(rf.est,function(x) sum(apply(x,2,function(y)  mean((y-0.01)^2)+(mean(y-0.01))^2 ))))

## check for badestimates
sapply(1:15, function(x) which(rf.est[[x]] > 0.2))

## plot the results:
tiff("Figure4.tiff",width=3.6*450,height=3.6*450, type='cairo', res=450, family="Arial",compression='lzw')
par(mar=c(5.1,4.1,1.1,1.1), cex=0.8)
plot(MSE~depths, xlab=expression(paste("Mean Read Depth, ",italic(mu[d[j]]))), ylab = "Sum of Mean Square Errors", pch=19,ylim=c(0,max(MSE)))
dev.off()

#save.image("sim_fixSeq.RData")

##########################################
## Look at the relationship between the number of individuals and rf and the ability to infer phase.

## Set up the simulation parameters
NoDS <- 1000
nSnps <- 12
rf <- 0.01; ep=0.002
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
rd_dist="NegBinom"

depths <- c(1,2,5,10)
rf_V <- c(0.01,0.05,0.1,0.2)
nInd_V <- c(25,50,75,100,150,200)

## Undertake the simulation
OPGP.est <- replicate(length(depths),list(replicate(length(rf_V),list(replicate(length(nInd_V),matrix(nrow=NoDS,ncol=nSnps),simplify=FALSE)))))
for(meanDepth in depths){
  for(rf in rf_V){
    print(rf)
    for(nInd in nInd_V){
      print(nInd)
      for(i in 1:NoDS){
        print(i)
        seed = 74832 + 13*i
        newData <- simFS(rf,epsilon=ep,config=config,nInd=nInd,meanDepth=meanDepth,
                         formats=NULL, NoDS=1, rd_dist = rd_dist, seed1=11984321, seed2=seed)
        for(j in (10)^(-(2:11))){
          tryOPGP <- tryCatch( infer_OPGP_FS(newData$depth_Ref, newData$depth_Alt, newData$config, 0.001, ndeps=rep(j,sum(config==1)*2+sum(config!=1)-1)), error =function(x) NULL )
          if(!is.null(tryOPGP)){
            OPGP.est[[which(meanDepth == depths)]][[which(rf == rf_V)]][[which(nInd == nInd_V)]][i,] <- tryOPGP
            break
          }
        }
      }
    }
  }
}
## find the ture OPGP value
OPGP <- newData$OPGP

#### Plot the results
tiff("FigureS14.tiff",width=7.5*450,height=7.5*450, type='cairo', res=450, family="Arial",compression='lzw')
par(mfrow=c(2,2), mar=c(6.1,4.1,1.1,1.1))
for(i in 1:length(depths)){
  OPGPc <- lapply(OPGP.est[[i]][[1]], function(x) mean(apply(x,1,function(y) all(y == OPGP)))*100)
  plot(nInd_V, unlist(OPGPc), ylim=c(0,100), xlim=c(min(nInd_V), max(nInd_V)), type='l',
       ylab="Correct OPGP (%)", xlab="Number of Progeny")
  points(nInd_V, unlist(OPGPc), col=1, pch=15)
  #mtext(switch(i,"(a)","(b)","(c)","(d)"), line = 4.5, side=1)
  for(j in 2:length(rf_V)){
    OPGPc <- lapply(OPGP.est[[i]][[j]], function(x) mean(apply(x,1,function(y) all(y == OPGP)))*100)
    lines(nInd_V, unlist(OPGPc), col=j)
    points(nInd_V, unlist(OPGPc), col=j, pch=15+j)
  }
  d = depths[i]
  legend('bottomright', legend=bquote(italic(mu[d[j]])==.(d)), bty='n')
}
legend('right', legend=rf_V, title=expression(italic(r[j])), col=1:length(rf_V),
       pch=c(15,15+2:length(depths)), lty=1, bty = 'n')
dev.off()

#save.image("sim_OPGP.RData")