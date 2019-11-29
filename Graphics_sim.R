#####################################################
## Code for generating the simulated data and plots 
## for Section 3.4
#####################################################

library(GUSbase)

## Simulation variables
ploid=1    
nInd = 200
nSnps = 50000
meanDepth <- rgamma(nSnps, shape = 1, rate=0.1)
set.seed(5743)
p <- rbeta(n = nSnps, 0.5, 0.5)
## function for simulating genotypes under HWE
samGeno <- function(x) sum(sample(x=c(rep(1,x),rep(0,2*ploid-x)),size=ploid))
samGeno <- Vectorize(samGeno)

## Simulate data
geno <- matrix(rbinom(nInd*nSnps,2*ploid,p),ncol=nSnps,byrow=T)
depth <- matrix(sapply(meanDepth, function(x) rnbinom(nInd,mu=x, size=2*ploid)),ncol=nSnps)
Acounts <- matrix(rbinom(nInd*nSnps,depth,geno/(2*ploid)),ncol=nSnps)
genoGBS <- (2*ploid)*Acounts/depth
Bcounts <- depth - Acounts
## Check for SNPs with no data and remove
badSNPs <- apply(depth, 2, function(x) all(x == 0))
ref <- Acounts[,-which(badSNPs)]
alt <- Bcounts[,-which(badSNPs)]

## Comet plot
cometPlot(ref, alt, ploid = 2, cex = 1.3, maxdepth = 200)

## RocketPlot
rocketPlot(ref, alt, ploid = 2, cex = 1.3, maxdepth = 200)
