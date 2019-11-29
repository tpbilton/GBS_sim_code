#####################################################
## Simulation code for Section 6.4.1
## - Generation of simulated data and constructing of GRMs
## - Was designed to run on HPC (namely NeSI) with slurm job scheduler (with two arugments supplied for ploid and meanDepth)
#####################################################

## Set parameters
alpha =  1/7   # Double reduction parameter
epsilon = 0.01 # sequencing error parameter
nSires <- 20   # Number of sires
nDams <- 5     # Number of dams per sire
nProg <- 2     # Number of progeny for each pair of parents
nSelf <- 1     # Number of progeny from mating of siblings
nInd <- nSires+(nDams*nSires)*(1+nProg)+ nSelf*choose(nProg,2)*nDams*nSires # Number of individuals in total
nSnps <- 5000  # number of SNPs
ploid <- as.integer(commandArgs(trailingOnly=T))[2]
meanDepth <- as.integer(commandArgs(trailingOnly=T))[3]
set = paste0("ploid",ploid*2,"_DR")

# sort out structure of the population
#matrix(sapply(1:nSires,function(x) rbind(rep(x,nDams), 20+(1:nDams)+(x-1)*5)), ncol=2,byrow=T)
sirepos <- as.vector(sapply(1:nSires,rep,times=nDams*nProg))
dampos <- as.vector(sapply(1:(nDams*nSires) + nSires,rep,times=nProg))
progpos <- 1:(nProg*nDams*nSires) + nSires + nDams*nSires
# index the siblings in the 
uSib <- matrix(dampos,nrow=nProg*nDams*nSires,ncol=nProg*nDams*nSires) ==  matrix(dampos,nrow=nProg*nDams*nSires,ncol=nProg*nDams*nSires,byrow=T)
uSib[lower.tri(uSib)] <- F
diag(uSib) <- F
uSib <- which(uSib,arr.ind=T) + nSires + nDams*nSires
sibpos <- max(uSib)+1:nrow(uSib)

## index which individuals should have the same self-relatedness value
sr1 <- 1:((nDams*nSires) + nSires)
sr2 <- max(sr1) + 1:(nProg*nDams*nSires)
sr3 <- max(sr2)+1:(length(sr2)/2)
## index which individuals should have the same relatedness value
library(AGHmatrix)
ped <- cbind(1:(nSires*nDams+nSires),rep(0,nSires*nDams+nSires),rep(0,nSires*nDams+nSires))
ped <- rbind(ped, cbind(progpos,sirepos,dampos))
ped <- rbind(ped, cbind(sibpos,uSib))
colnames(ped) <- c("Ind","Par1","Par2")
ped <- as.data.frame(ped)
for(i in 1:ncol(ped))
  ped[,i] <- factor(ped[,i])

Amat <- Amatrix(ped, ploidy = 2*ploid, w=alpha, slater = TRUE)

tab <- table(Amat)[1:4]
rr1 <- which(Amat == names(tab)[1], arr.ind=T)
rr1 <- rr1[-which(rr1[,2] < rr1[,1]),]
rr2 <- which(Amat == names(tab)[2], arr.ind=T)
rr2 <- rr2[-which(rr2[,2] < rr2[,1]),]
rr3 <- which(Amat == names(tab)[3], arr.ind=T)
rr3 <- rr3[-which(rr3[,2] < rr3[,1]),]
rr4 <- which(Amat == names(tab)[4], arr.ind=T)
rr4 <- rr4[-which(rr4[,2] < rr4[,1]),]
tab <- as.numeric(c(names(tab)[1:4], Amat[sr1[1],sr1[1]], Amat[sr2[1],sr2[1]], Amat[sr3[1],sr3[1]] ))

run <- as.integer(commandArgs(trailingOnly=T))[1]

if(run == 1) save.image(paste0("/nesi/nobackup/uoo02457/polysim/Inbred/",set,"/",set,"_depth",meanDepth,"simPara.Rdata"))

### VanRaden GRM
calcGRM_vanRaden <- function(geno, depth = NULL, phat, ep=0.01, thres=0.01){
  if(is.null(depth)){
    snpsubset <- which(phat < 1-thres & phat > thres)
    nsnpsub <- length(snpsubset)
    phat <- phat[snpsubset]
    genon0 <- geno[,snpsubset] - (2*ploid)*rep.int(phat, rep(nrow(geno), nsnpsub))
    genon0[is.na(geno[,snpsubset])] <- 0
    P0 <- matrix(phat,nrow=nInd,ncol=nsnpsub,byrow=T)
    P1 <- 1-P0
    P0[is.na(geno[,snpsubset])] <- 0
    P1[is.na(geno[,snpsubset])] <- 0
    div0 <- (2*ploid)*tcrossprod(P0,P1)
    GRM <- tcrossprod(genon0)/div0
    return(GRM)
  } else{
    snpsubset <- which(phat < 1-thres & phat > thres)
    depth <- depth[,snpsubset]
    nsnpsub <- length(snpsubset)
    phat <- phat[snpsubset]
    if(length(ep) > 1) ep = matrix(ep[snpsubset], nrow=nrow(geno), ncol=nsnpsub, byrow=T)
    else               ep = matrix(rep(ep,nsnpsub), nrow=nrow(geno), ncol=nsnpsub, byrow=T)
    
    ## Compute the adjusted GRM
    genon0 <- geno[,snpsubset] - (2*ploid)*rep.int(phat, rep(nrow(geno), nsnpsub))
    genon0[depth<1] <- 0
    P0 <- matrix(phat,nrow=nInd,ncol=nsnpsub,byrow=T)
    P1 <- 1-P0
    P0[depth<1] <- 0
    P1[depth<1] <- 0
    ep[depth<1] <- 0
    div0 <- (2*ploid)*tcrossprod(P0,P1)
    GRM <- (tcrossprod(genon0/sqrt(1-4*ep*(1-ep))) - tcrossprod(sqrt(((2*ploid*ep)^2*(1-4*P0*P1)/(1-4*ep*(1-ep))))) )/div0
    genon0[depth<2] <- 0
    depth.temp <- depth
    depth.temp[which(depth < 2)] <- 0
    depth.temp <- 1/depth.temp
    depth.temp2 <- depth.temp
    depth.temp[is.infinite(depth.temp)] <- 1
    depth.temp2[is.infinite(depth.temp2)] <- 0
    P0[depth<2] <- 0
    P1[depth<2] <- 0
    ep[depth<2] <- 0
    adj <- (2*ploid)^2*(P0*P1*(depth.temp + 4*ep*(1-ep)*(1-depth.temp)) +
                                      ep*(ep+(1-ep)*depth.temp - 4*P0*P1))
    diag(GRM) <- rowSums(  (genon0^2 - adj)/((1-depth.temp2)*(1-4*ep*(1-ep))))/(2*ploid*rowSums(P0*P1))
    return(GRM)
  }
}

### Weir GRM
calcGRM_weir <- function(geno, depth=NULL, ep=NULL){
  if(!is.null(ep) & !is.null(depth)){
    geno[which(depth < 1)] <- NA
    snpsubset <- which(!(colMeans(geno, na.rm=T) %in% c(0, 2*ploid)))
    geno <- geno[, snpsubset]
    depth <- depth[, snpsubset]
    na_indx <- !is.na(geno)
    na_mat <- tcrossprod(na_indx,na_indx)
    geno[which(is.na(geno))] <- ploid

    epMat <- matrix(ep, nrow=nrow(geno), ncol=ncol(geno))
    epMat[which(depth < 1)] <- 0
    
    mat <- 1/2 + 2/(2*ploid)^2*(tcrossprod((geno-ploid)/sqrt(1-4*epMat*(1-epMat))) - 
                                  tcrossprod(sqrt((((2*ploid)^2/4)*na_indx)/(1-4*epMat*(1-epMat)))) + tcrossprod(sqrt((2*ploid)^2/4*na_indx)) +
                                  tcrossprod(sqrt((2*ploid)^2*epMat*(1-epMat)/(1-4*epMat*(1-epMat)))))/na_mat
    
    geno[which(depth < 2)] <- NA
    drat <- 1/depth
    drat[which(depth < 2)] <- 0
    epMat <- matrix(ep, nrow=nrow(geno), ncol=ncol(geno))
    epMat[which(depth < 2)] <- 0
    na_indx <- !is.na(geno)
    na_mat <- tcrossprod(na_indx,na_indx)
    geno[which(is.na(geno))] <- ploid
    diag(mat) <- 1/2 + 2/(2*ploid)^2*rowSums( (((geno-ploid)^2 - (2*ploid)^2/4)/(1-drat) + (2*ploid)^2*epMat*(1-epMat))/(1-4*epMat*(1-epMat)) + (2*ploid)^2/4)/diag(na_mat)
  }
  else if(!is.null(depth)){
    geno[which(depth < 2)] <- NA
    snpsubset <- which(!(colMeans(geno, na.rm=T) == 2*ploid))
    geno <- geno[, snpsubset]
    na_mat <- !is.na(geno)
    na_mat <- tcrossprod(na_mat,na_mat)
    geno[which(is.na(geno))] <- ploid
    
    drat <- 1/depth[,snpsubset]
    drat[which(depth[,snpsubset] < 2)] <- 0
    
    mat <- 1/2 + 2/(2*ploid)^2*tcrossprod(geno - ploid)/na_mat
    diag(mat) <- 1/2 + 2/(2*ploid)^2*rowSums(((geno-ploid)^2 - (2*ploid)^2/4)/(1-drat) + (2*ploid)^2/4)/diag(na_mat)
  } else{
    snpsubset <- which(!(colMeans(geno, na.rm=T) == 2*ploid))
    geno <- geno[, snpsubset]
    nSnps = length(snpsubset)
    na_mat <- !is.na(geno)
    na_mat <- tcrossprod(na_mat,na_mat)
    geno[which(is.na(geno))] <- ploid
    
    mat <- 1/2 + 2/(2*ploid)^2*tcrossprod(geno - ploid)/na_mat
  }
  mat_sum <- (sum(mat) - sum(diag(mat)))/(nrow(mat)*(nrow(mat)-1))
  GRM_weir <- 2*ploid*(mat - mat_sum)/(1-mat_sum)
  return(GRM_weir)
}

IBDgmat <- function(IBD, ploid){
  nInd <- nrow(IBD)/ploid
  gmat <- matrix(NA, nrow=nInd, ncol=nInd)
  for(i in 1:nInd){
    mat1 <- IBD[rep(1:ploid + (i-1)*ploid,ploid),]
    for(j in i:nInd){
      mat2 <- IBD[rep(1:ploid + (j-1)*ploid,rep(ploid,ploid)),]
      theta <- sum(mat1==mat2)/(ncol(IBD)*ploid^2)
      gmat[i,j] <- theta*ploid
    }
  }
  gmat[lower.tri(gmat)] <- t(gmat)[lower.tri(gmat)]
  return(gmat)
}

### Run the simulation

set.seed(run*39+meanDepth*5+ploid*111)
## simulate the minor allele frequencies for each SNP
p <- runif(n = nSnps, min=0, max = 1)

## List of Computed G matrices
Gmat <- vector(mode = "list", length=0)

DR_sam_IBD <- function(){
  sam_allele <- rep(NA,ploid)
  countV <- numeric(2*ploid)
  sam_allele[1] <- sample(1:(2*ploid), size=1)
  countV[sam_allele[1]] <- countV[sam_allele[1]] + 1
  for(i in 2:ploid){
    ind_sam <- which(countV == 1)
    ind_unsam <- which(countV == 0)
    sam_allele[i] <- sample(c(ind_unsam,ind_sam), size=1, 
                            prob = c(rep((1-alpha)/length(ind_unsam),length(ind_unsam)),
                                      rep(alpha/length(ind_sam),length(ind_sam))))
    countV[sam_allele[i]] <- countV[sam_allele[i]] + 1
  }
  return(sam_allele)
}

### Simulate the parents 
sim_par <- sapply(p, function(p){
  nFound <- nSires+nDams*nSires
  ibd <- nFound*(2*ploid)
  allele <- numeric(nFound*2*ploid)
  
  ibd <- as.vector(sapply(1:(nFound*2),function(x) sort(as.numeric(factor(DR_sam_IBD()))) + (x-1)*ploid))
  allele <- sample(0:1,size=nFound*2*ploid,prob=c(1-p,p),replace=T)[ibd]
  geno <- colSums(matrix(allele, nrow=2*ploid))
  return(list(paste0("a",ibd),allele,geno))
}, simplify = F)

IBD_par <- do.call("cbind" ,lapply(sim_par, function(x) x[[1]]))
allele_par <- do.call("cbind" ,lapply(sim_par, function(x) x[[2]]))
geno_par <- do.call("cbind" ,lapply(sim_par, function(x) x[[3]]))

rm(sim_par)
gc()

### simulate the progeny
sim_prog <- sapply(1:nSnps, function(x){
  indx_prog_mat <- as.vector(replicate(length(progpos), DR_sam_IBD())) + rep((dampos-1)*2*ploid,rep(ploid,length(dampos)))
  indx_prog_pat <- as.vector(replicate(length(progpos), DR_sam_IBD())) + rep((sirepos-1)*2*ploid,rep(ploid,length(sirepos)))
  indx_prog <- as.vector(rbind(matrix(indx_prog_mat,nrow=ploid),matrix(indx_prog_pat,nrow=ploid)))
  ibd <- IBD_par[indx_prog,x]
  allele <- allele_par[indx_prog,x]
  geno <- colSums(matrix(allele, nrow=2*ploid))
  return(list(ibd,allele,geno))
}, simplify = F)

IBD_prog <- do.call("cbind" ,lapply(sim_prog, function(x) x[[1]]))
allele_prog <- do.call("cbind" ,lapply(sim_prog, function(x) x[[2]]))
geno_prog <- do.call("cbind" ,lapply(sim_prog, function(x) x[[3]]))

rm(sim_prog)
gc()

# simulate genotypes from progeny of sibling mating
uSib1 <- uSib[,1] - min(uSib)
uSib2 <- uSib[,2] - min(uSib)

sim_sib <- sapply(1:nSnps, function(x){
  indx_mat <- as.vector(replicate(length(uSib1), DR_sam_IBD())) + rep(uSib1*2*ploid,rep(ploid,length(uSib1)))
  indx_pat <- as.vector(replicate(length(uSib2), DR_sam_IBD())) + rep(uSib2*2*ploid,rep(ploid,length(uSib2)))
  indx <- as.vector(rbind(matrix(indx_mat,nrow=ploid),matrix(indx_pat,nrow=ploid)))
  ibd <- IBD_prog[indx,x]
  allele <- allele_prog[indx,x]
  geno <- colSums(matrix(allele, nrow=2*ploid))
  return(list(ibd,allele,geno))
}, simplify = F)

IBD_sib <- do.call("cbind" ,lapply(sim_sib, function(x) x[[1]]))
allele_sib <- do.call("cbind" ,lapply(sim_sib, function(x) x[[2]]))
geno_sib <- do.call("cbind" ,lapply(sim_sib, function(x) x[[3]]))

rm(sim_sib)
gc()

## Combine information
IBD <- rbind(IBD_par,IBD_prog,IBD_sib)
geno <- rbind(geno_par,geno_prog,geno_sib)

## simulation the GBS data:
depth <- matrix(rnbinom(nInd*nSnps,mu=meanDepth, size=2*ploid),ncol=nSnps)
Acounts <- matrix(rbinom(nInd*nSnps,depth,geno/(2*ploid)),ncol=nSnps)
genoGBS <- (2*ploid)*Acounts/depth

## simulate the errors
Bcounts <- depth - Acounts
aCountsFinal <- matrix(rbinom(nInd*nSnps,Acounts,prob=1-epsilon),ncol=nSnps) + matrix(rbinom(nInd*nSnps,Bcounts,prob=epsilon),ncol=nSnps)
genoGBS_ep <- (2*ploid)*aCountsFinal/depth


##########################################
############# Now compute the GRMs

###### Case 1: VanRaden (known p) ######
# Genotypes
snpsubset <- which(p < 1 & p > 0)
nsnpsub <- length(snpsubset)
psub <- p[snpsubset] 
Gmat$Geno <- tcrossprod(geno[,snpsubset]-(2*ploid)*rep.int(psub,rep(nrow(geno),nsnpsub)))/(2*ploid*sum(psub*(1-psub)))

# GBS data
Gmat$GBS <- calcGRM_vanRaden(genoGBS, depth=depth, p, 0)

# GBS data with no adjustment for depth 
Gmat$GBS_noAdj <- calcGRM_vanRaden(genoGBS, depth=NULL, p, 0)

# GBSep data with (ep=ep)
Gmat$GBSep <- calcGRM_vanRaden(genoGBS_ep, depth=depth, p, epsilon)

# GBSep data with (ep=0)
Gmat$GBSep_0 <- calcGRM_vanRaden(genoGBS_ep, depth=depth, p, 0)

###### Case 2: VanRaden (Estimated p from whole pedigree) ######
# Genotypes
phat <- colMeans(geno, na.rm=T)/(2*ploid)
snpsubset <- which(phat < 1-0.01 & phat > 0.01)
nsnpsub <- length(snpsubset)
phat_sub <- phat[snpsubset]
Gmat$Geno_p <- tcrossprod(geno[,snpsubset]-(2*ploid)*rep.int(phat_sub,rep(nrow(geno),nsnpsub)))/(2*ploid*sum(phat_sub*(1-phat_sub)))

# GBS data
phat <- colMeans(genoGBS, na.rm=T)/(2*ploid)
Gmat$GBS_p <- calcGRM_vanRaden(genoGBS, depth=depth, phat, 0)

# GBS data with no adjustment for depth 
phat <- colMeans(genoGBS, na.rm=T)/(2*ploid)
Gmat$GBS_p_noAdj <- calcGRM_vanRaden(genoGBS, depth=NULL, phat, 0)

# GBSep data with (ep=ep)
phat <- colMeans(genoGBS_ep, na.rm=T)/(2*ploid)
Gmat$GBSep_p <- calcGRM_vanRaden(genoGBS_ep, depth=depth, phat, epsilon)

# GBSep data with (ep=0)
phat <- colMeans(genoGBS_ep, na.rm=T)/(2*ploid)
Gmat$GBSep_p_0 <- calcGRM_vanRaden(genoGBS_ep, depth=depth, phat, 0)

###### Case 3: VanRaden (Estimated p from founders) ######
# Genotypes
phat <- colMeans(geno[sr1,], na.rm=T)/(2*ploid)
snpsubset <- which(phat < 1-0.01 & phat > 0.01)
nsnpsub <- length(snpsubset)
phat_sub <- phat[snpsubset]
Gmat$Geno_p_foun <- tcrossprod(geno[,snpsubset]-(2*ploid)*rep.int(phat_sub,rep(nrow(geno),nsnpsub)))/(2*ploid*sum(phat_sub*(1-phat_sub)))

# GBS data
phat <- colMeans(genoGBS[sr1,], na.rm=T)/(2*ploid)
Gmat$GBS_p_foun <- calcGRM_vanRaden(genoGBS, depth=depth, phat, 0)

# GBS data with no adjustment for depth 
phat <- colMeans(genoGBS[sr1,], na.rm=T)/(2*ploid)
Gmat$GBS_p_foun_noAdj <- calcGRM_vanRaden(genoGBS, depth=NULL, phat, 0)

# GBSep data with (ep=ep)
phat <- colMeans(genoGBS_ep[sr1,], na.rm=T)/(2*ploid)
Gmat$GBSep_p_foun <- calcGRM_vanRaden(genoGBS_ep, depth=depth, phat, epsilon)

# GBSep data with (ep=0)
phat <- colMeans(genoGBS_ep[sr1,], na.rm=T)/(2*ploid)
Gmat$GBSep_p_foun_0 <- calcGRM_vanRaden(genoGBS_ep, depth=depth, phat, 0)

###### Case 4: Weir method ######
# Genotypes
Gmat$Weir <- calcGRM_weir(geno)

# GBS data
Gmat$GBS_Weir <- calcGRM_weir(genoGBS, depth, ep=0)

# GBS data with no adjustment for depth 
Gmat$GBS_Weir_noAdj <- calcGRM_weir(genoGBS)

# GBSep data with (ep=ep)
Gmat$GBSep_Weir <- calcGRM_weir(genoGBS_ep, depth=depth, ep=epsilon)

# GBSep data with (ep=0)
Gmat$GBSep_Weir_0 <- calcGRM_weir(genoGBS_ep, depth=depth, ep=0)

# True IBD sharing
Gmat$true <- IBDgmat(IBD, 2*ploid)

# Pedigree estimator
#Gmat$Ped <- Amat

save(Gmat, file=paste0("/nesi/nobackup/uoo02457/polysim/Inbred/",set,"/",set,"_depth",meanDepth,"_pop_Gmat_run",run,".Rdata"))

