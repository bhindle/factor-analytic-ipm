### This code provides a framework for running IPMs using a factor analytic approach. Here we use
### demographic data simulated using an IBM. We run the FA model using JAGs and use these parameter 
### estimates to parameterise an IPM. 

##################################################################################################
## Load required packages etc
##################################################################################################
require(runjags); require(dplyr); require(tidyr);
require(purrr); require(MASS)

set.seed(11011989)

##################################################################################################
## Simulate demographic data using an IBM 
##################################################################################################

## NB: demographic parameters are set within SimulateIBM.R

## Set parameters
tot.yrs <- 5000 ## Total length of IBM simulation
run.in <- 1000 ## Run in period
init.pop.size <- 500 ## Initial population size
n.yrs <- 12 ## Length of data set to extract
n.ind <- 20:150 ## Range of no. of individuals per year

## Run IBM
source("SimulateIBM.R")

##################################################################################################
## Run JAGs model 
##################################################################################################

## Format data 
grow <- filter(sim.data, Surv==1)
recs <- filter(sim.data, !is.na(Rcsz))

dataList <- list(## Survival
                 Nsurv = nrow(sim.data),
                 zsu = sim.data$z,
                 surv = sim.data$Surv,
                 tsu = sim.data$Year,
                 ## Growth
                 Ngrow = nrow(grow),
                 zgr = grow$z,
                 grow = grow$z1,
                 tgr = grow$Year,
                 ## Reproduction
                 rep = grow$Repr,
                 ## Recruit size
                 Nrs = nrow(recs),
                 zrs = recs$z,
                 rs = recs$Rcsz,
                 trs = recs$Year, 
                 Yrs = length(unique(sim.data$Year)))

## Parameters to monitor
submod <- c("su", "gr", "rep", "rs")
parKeep <- c(paste0("b0", submod), paste0("bz", submod), paste0("bq", submod), 
             "Q", "ts", "tg", "trp", "tr", 
             "sigma.ts", "sigma.tg", "sigma.trp", "sigma.tr", 
             "sigma.rs", "sigma.gr")

## Set initial parameters here

## Run model
postsamp <- run.jags("FAmodel.jags", monitor=parKeep, data=dataList, n.chains= 2, 
                     burnin=10000, adapt=1000, thin=200, 
                     sample=10000)
# It is simple to parallise this step using e.g. method = parallel or rjparallel

## NB: various checks should be carried out here e.g. convergence of chains and  
## that a single latent parameter is sufficient to account for the covariation among the vital rates
## See Appendix A1c for more information

##################################################################################################
## Run IPM
##################################################################################################

source("IPMfunctions.R")

## Extract parameter estimates from jags output
samp <- do.call(rbind, postsamp$mcmc)
outmns <- data.frame(t(apply(samp, 2, mean)))
## NB: for simplicity here we just use the parameter means as estimates - would be better to sample 
## from posterior - see effect of parameter uncertainty on model output

## Iterate model and calculate population growth rate
stlamb <- stochlamb(outmns, num.yrs = 5000, nrunin = 1000)

## NB: various checks should also be carried out at this point e.g. that no eviction is occuring 
## from the IPM and comparing model predictions to observed population dynamics



