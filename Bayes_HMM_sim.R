#####################################################
## Simulation code for Section 5.3.4
## - Generation of simulated data and fitting the Bayesian hierarchical HMM
## - Was designed to run on HPC (namely NeSI) with slurm job scheduler
## - Requires using the external program PedigreeSim (which is java code)
## - Note: Path names will need to be changed
#####################################################

.libPaths("/scale_wlg_persistent/filesets/home/bilti119/R/3.5")

## Simulation parameters
sim = as.integer(commandArgs(trailingOnly=T))[1] 
#nInd <- 177     # Number of progeny for each pair of parents
nInd = as.integer(commandArgs(trailingOnly=T))[2]
set.seed(58932+5*sim+13*nInd)

## Depths from the data 
load("/nesi/project/uoo02457/GUSMap_Bayes/Simulations/alleleCount.Rdata")
meanDepth <- c(2.8475,3.096,2.435,42.8644,4.8362,103.4689,2.6045,44.7684,3.5537,4.7175,43.887,
               5.9661,43.1186,64.096,54.9379,4.3446,43.1977,4.4633,2.3955,4.6667,91.1073,52.8192,4.6497,
               4.3333,44,4.9718,59.1751,3.4576,5.8023,50.8814,5.7345,56.6893,48.1073,2.8814,4.548,5.904,3.1638,
               50.6723,3.6836,82.5141,3.7062,5.0056,5.8588,4.9492,44.7175,40.3898,5.4181,5.1356,3.1243,4.0282,53.9492,
               2.1243,53.1695,2.5819,59.4915,3.3729,5.661,3.1864,3.3277,41.3051,5.4576,5.5254,46.2203,5.2542,39.6723,
               51.9944,4.6271,51.5367,2.6045,53.6102,4.9605,4.9661,5.6328,46.7797,5.661,44.3164,4.1073,3.2768,2.7571,
               4.9153,42.9831,3.5085,3.7458,45.5367,4.1356,4.5706,45.8644,5.9096,43.3559,45.1977,5.5311,4.2825,48.4859,
               39.8531,43.5932,3.9379,5.1582,55.9661,5.5876,4.4463,43.0226,49.4181,3.339,40.2938,5.9605,39.5763,59.5085,
               3.8927,4.0339,4.0395,5.7684,3.9266,4.1186,2.6723,3.1299,4.7062,4.1525,5.0282,3,39.3898,4.5537,2.4633,
               43.774,4.2316,5.2938,5.2994,5.9096,42.1695,4.1751,3.3333,46.0169,38.2147,5.0508,43.7684,4.7175,40.3333,
               5.3672,2.8927,4.7401,4.5028,4.8927,5.2542,44.661,49.8927,3.435,2.3729,4.8983,4.4802,43.5989)
samDepth <- c(10.3624,28.8322,26.0537,52.1745,32.1074,47.7383,10.9128,25.5906,20.4564,11.9463,35.5302,37.9597,28.7047,
              13.1812,9.4228,34.2416,25.1611,8.094,45.5101,20.5034,42.7919,26.4765,13.2013,50.4631,14.7047,25.0537,
              16.5503,21.6913,19.6644,15.094,2.4295,21.3893,17.3691,22.8322,16.2752,15.8255,19.443,18.0067,17.7315,
              11.4497,18.7248,21.7584,21.1208,26.2685,22.7852,7.3893,28.8322,12.604,32.3289,29.1074,20.5034,31.1141,
              5.2752,8.8859,16.3893,26.7718,35.2953,15.8121,6.1477,25.8725,10.6309,10.3221,25.7047,7.8255,11.8725,
              9.8054,42.2886,3.1477,14.4027,34.7315,20.094,9.8523,25.7919,35.2081,11.1477,25.9463,18.6644,13.6913,
              9.1007,27.5839,26.4027,20.7718,26.3221,28.6443,31.1879,28.4161,32.5101,6.8389,13.0201,32.2483,24.4698,
              47.5772,15.5034,26.8591,43.9732,12.1007,31.5436,25.3691,11.6644,40.0738,33.3423,25.1678,14.7785,7.8591,
              30.906,20.1946,8.3557,38.4228,19.6242,25.9195,24.7047,12.3758,31.8121,14.3893,21.5839,12.7987,28.1007,
              13.3087,22.745,21.8054,18.8859,2.8188,25.0403,24.1275,22.7919,16.2685,17.9866,19.953,21.0805,19.8121,
              13.9463,16.3893,21.5168,25.2349,21.2819,7.8725,25.557,14.5235,14.9262,13.8322,19.1074,5.5772,8.255,
              14.9195,27.3691,38.5302,14.6577,3.953,27.5168,11.5503,9.2148,22.9396,7.7047,9.4027,11.5034,20.4161,
              3.5705,18.4027,12.4094,36.4631,19.6711,8.2282,29.5235,25.8255,10.1678,21.1007,11.604,8.4094,8.5503,
              28.604,17.1611,23.5034,23.9396,23.6309,14.2886,30.5772,7.1946)
OPGP <- c(9,3,9,3,10,9,10,5,5,9,3,3,9,10,10,5,5,6,9,9,5,6,9,6,5,10,9,5,5,3,3,9,10,10,5,9,10,10,4,10,10,6,10,6,5,5,9,9,
          10,9,4,6,3,3,3,10,5,10,5,6,10,9,9,10,10,6,6,4,10,5,6,5,10,10,9,10,10,10,6,5,6,6,10,10,5,6,6,5,9,9,5,6,10,10,
          9,5,3,1,5,10,5,9,10,9,5,10,9,10,10,1,9,6,3,2,9,6,6,5,3,10,6,9,10,5,10,9,9,10,10,5,5,5,5,9,10,5,6,3,10,9,6,6,
          1,5,1,5,6,9,9)
rf_est <- c(0.5506,0.4151,0.3101,0.5462,0.6664,1.1655,1.2783,0.325,7.3798,0.6405,0.1909,0.4007,0.9475,0.2244,0.3669,
            0.7294,0.7525,0.9417,0.4879,1.3078,0.2656,0.7674,1.4355,0.428,0.4111,0.3091,1.0243,0.5922,0.6056,0.192,
            0.6542,0.6631,0.354,2.9248,0.478,2.6444,0.2774,0.7247,0.2755,0.2628,0.4301,0.8339,1.4252,0.7763,0.271,
            0.6145,1.8296,0.3172,0.2645,0.2731,0.3027,0.2067,0.177,0.1788,0.5642,0.3415,0.364,0.4523,0.246,0.3309,
            0.4444,0.9488,0.2327,0.3567,0.4168,0.2604,0.849,1.4836,0.4472,0.3317,0.3106,0.4252,0.2765,0.7334,0.241,
            0.2458,0.265,0.3308,1.2717,0.4098,0.2356,0.3251,0.3931,0.2617,0.2059,0.4726,0.1906,0.1911,0.2436,0.3628,
            0.3788,1.0417,0.761,0.4597,0.3114,0.2585,0.9245,0.4719,0.3728,0.324,0.2332,0.667,0.224,0.3334,0.3781,0.3601,
            0.2189,0.2385,1.6609,0.3288,0.2223,0.4038,0.2721,0.4431,0.5336,0.3267,0.3225,0.3067,0.4075,0.338,0.2608,
            0.3465,0.188,0.1803,0.1926,0.1826,0.3273,0.504,0.2748,0.283,0.2731,0.2656,0.2703,0.3037,0.1834,0.1858,0.3136,
            0.1918,0.3585,0.3545,0.3463,0.3451,0.1825,0.1876,0.1993,0.2487,0.2496,0.2499)
nSnps <- length(OPGP)
#epsilon = rbeta(nSnps, shape1=1, shape2=499)
epsilon = rbeta(nSnps, shape1=3, shape2=1497)
#epsilon <- scan(paste0("/nesi/project/uoo02457/GUSMap_Bayes/Ep/MK_both/Post_samples/ep_",sim,".txt"))
low <- which(meanDepth < 6)
high <- which(meanDepth > 6)

## filenames
genfile = paste0("sim",sim,".gen")
chromfile = paste0("sim",sim,".chrom")
parfile = paste0("sim",sim,".par")
mapfile = paste0("sim",sim,".map")
outfile = paste0("out_sim",sim)

###############################
#### Generate the parents genotypes file:

OPGPtoParHap <- function(OPGP){
  parHap <- matrix("A", nrow=4, ncol=length(OPGP))
  parHap[1,which(OPGP %in% c(2,4,6,8,11,12,15,16))] <- "B"
  parHap[2,which(OPGP %in% c(1,3,5,7,11,12,15,16))] <- "B"
  parHap[3,which(OPGP %in% c(3,4,7,8,10,12,14,16))] <- "B"
  parHap[4,which(OPGP %in% c(1,2,7,8,9,11,14,16))] <- "B"
  return(parHap)
}
geno_par <- cbind(paste0("Mar",1:nSnps),t(OPGPtoParHap(OPGP)))
colnames(geno_par) <- c("marker", paste0(c(rep("P1_",2),rep("P2_",2)),1:2))
write.table(geno_par, file=genfile, quote=FALSE, row.names=F, sep="\t")

###############################
#### Generate the map file:
snpdist <- rgamma(nSnps-1, shape=rf_est/sum(rf_est)*(nSnps-1), rate=1.85)
pos <- c(0,round(cumsum(snpdist),6))
linkfile <- cbind(marker=paste0("Mar",1:nSnps), chromosome="chr1", position=pos)
write.table(linkfile, file=mapfile, quote=FALSE, row.names=F, sep="\t")

###############################
#### Generate the chromosome file:
cat(c("chromosome",	"length",	"centromere"),"\n",sep="\t", file=chromfile)
cat(c(paste0("chr",1),max(pos), round(5/8*max(pos),6)),"\n", sep="\t", file=chromfile, append=T)

###############################
#### Generate the parameter file
parinfo <- paste(paste0("PLOIDY = 2"),
                 paste0("MAPFUNCTION = HALDANE"),
                 "MISSING = NA",
                 paste0("CHROMFILE = ",chromfile),
                 paste0("MAPFILE = ", mapfile),
                 paste0("FOUNDERfile = ", genfile),
                 paste0("OUTPUT = ", outfile,"\n"),
                 "POPTYPE = F1",
                 paste0("POPSIZE = ",nInd),
		 paste0("SEED = ",2984 + 13*sim + 17*nInd),"\n",
                 sep="\n")
cat(parinfo, file=parfile)

Sys.sleep(1)

## Run PedigreeSim
system(paste0("java -Xmx1g -jar /home/bilti119/PedigreeSim/PedigreeSim.jar sim",sim,".par")) # Path will need to be changes

Sys.sleep(1)

## Now simulate the GBS data
# determine genotypes
geno <- read.table(paste0("out_sim",sim,"_alleledose.dat"), header=T)
geno <- t(as.matrix(geno[,-c(1:3)]))
for(snp in 1:nSnps){
  if((sum(geno[,snp]==0) > sum(geno[,snp]==2)) & (sum(geno[,snp]==2) == 0)){
    geno[,snp] <- 2-geno[,snp]
  }
}

# specify individual read depths
depth_indx <- sample(1:177,size=nInd)
depth <- ref + alt
depth <- depth[depth_indx,]

# sample alleles from true genotype
Acounts <- matrix(rbinom(nInd*nSnps,depth,geno/2),nrow=nInd, ncol=nSnps)
Bcounts <- depth - Acounts

# simulate sequencing errors
ref <- matrix(as.integer(rbinom(nInd*nSnps,Acounts,prob=1-epsilon)),ncol=nSnps) + 
  matrix(as.integer(rbinom(nInd*nSnps,Bcounts,prob=epsilon)),ncol=nSnps)
alt <- depth - ref

# Determine the OPGP
seg <- (OPGP < 5)*1 + (OPGP > 4 & OPGP < 7)*2 + (OPGP > 8 & OPGP < 11)*4
OPGP_temp <- as.integer(GUSMap:::infer_OPGP_FS(ref, alt, seg, method="EM"))

noFam=as.integer(1)

library(ggplot2, lib.loc="/scale_wlg_persistent/filesets/home/bilti119/R/3.5")
library(htmlwidgets, lib.loc="/scale_wlg_persistent/filesets/home/bilti119/R/3.5")
library(GUSMap)

## find the MLE estimates for Frequentist appraoch
esttime <- list()
timestart <- proc.time()[3]
MLE_both <- GUSMap:::rf_est_FS(init_r = 0.001, ep = 0.001, ref = list(ref), 
                 alt = list(alt), OPGP = list(OPGP_temp), 
                 sexSpec = F, seqErr = T, method = "EM", 
                 nThreads = as.integer(1), multiErr = TRUE)
timestop <- proc.time()[3]
esttime$MLE_both <- timestop - timestart

output <- list()
output$ref <- ref
output$alt <- alt
output$OPGP_infer <- OPGP_temp
output$pos <- pos
output$epsilon <- epsilon
output$nInd <- nInd
output$low <- low
output$high <- high
output$MLE_both <- MLE_both


## find the MLE using the Hirarchical Bayesian approach
library(doParallel)
library(BUSMap)

## MCMC parameters
nchain <- 3
simPar <- c(5000,25000)

## starting values
startVal <- replicate(nchain, c(rnorm(nSnps-1,-5,1), rnorm(nSnps,-5,0.8),
                        rnorm(1,-5,1),rnorm(1,-5,0.8),rlnorm(1, sdlog=0.5),rlnorm(1, sdlog=0.5)), simplify=F)

## Fit the model
registerDoParallel(nchain)
out  <- foreach(i=1:nchain) %dopar% {
	timestart <- proc.time()[3]
        res <- MH_Bayes_Hir(ref,alt,OPGP_temp,nInd,nSnps,startVal[[i]],simPar,11093+13*i+sim*5)
	timet <-  proc.time()[3] - timestart
	return(list(time=timet,samples=res))
}
esttime$Bayes_both <- unname(unlist(lapply(out, function(x) x$time)))
output$Bayes_both <- lapply(out, function(x) x$samples)

output$time     <- esttime
output$MH_para <- list(nchain = nchain, nadapt = simPar[1], nmh = simPar[2])

save(output, file=paste0("/nesi/nobackup/uoo02457/GUSMap_Bayes/Simulations/Data/sim",sim,".Rdata"))

