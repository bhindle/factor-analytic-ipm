
################# Code to simulate demographic data #########################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Define 'true' parameter vector
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m.par <- c(## survival
  surv.int  = 3.44,
  surv.z = -0.62,
  yr.surv = 1.05,
  ## growth 
  grow.int  =  2.38,
  grow.z    =  0.24,
  grow.sd   =  0.14,
  yr.grow = 0.04,
  ## reproduce or not
  rep.int  =  -7.37,
  rep.z    =  2.88, 
  yr.rep = 0.46,
  ## recruit or not
  recr.int  =  0.80,
  ## recruit size
  rcsz.int  =  1.22,
  rcsz.z    =  0.41,
  rcsz.sd   =  0.20,
  yr.recsz = 0.07)

## Correlation matrix
corrmat <- matrix(c(1, 0.1, 0.3, 0.05, 
                    0.1, 1, 0.05, 0.6,
                    0.3, 0.05, 1, 0.1, 
                    0.05, 0.6, 0.1, 1), nrow = 4)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run simulation
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## initial size distribution (assuming everyone is a new recuit)

z <- rnorm(init.pop.size, mean = m.par["rcsz.int"] +  m.par["rcsz.z"] * 3.2, sd = m.par["rcsz.sd"])

## vector to store pop size 

stoch.lamb <- pop.size.t <- rep(NA, tot.yrs)

# Estimate year effects for length simulation

## Covariance matrix
covmat <- matrix(NA, nrow=4, ncol=4)
yrvar <- c("yr.surv", "yr.grow", "yr.rep", "yr.recsz")
for (i in 1:4){
  for (j in 1:4){
    covmat[i,j] <- corrmat[i,j] * m.par[yrvar[i]] * m.par[yrvar[j]]    
  }
}

# Sample from a multivariate distribution with means of 0 and covariance matrix covmat for the number of 
# years in the simulation 

yref <- mvrnorm(tot.yrs, mu=rep(0, 4), covmat)
colnames(yref) <- c("surv", "grow", "rep", "recsz") 

## Iterate the model using the "true" parameters and store data in a data.frame
yr <- 1
sim.data <- list()

for (i in 1:tot.yrs){
  
  ## calculate current population size
  pop.size <- length(z)
  
  ## generate binomial random number for survival, where survival depends on your size z,
  ## this is a vector of 0's and 1's, you get a 1 if you survive
  mn <- m.par["surv.int"] + m.par["surv.z"] * z + yref[i,"surv"]
  surv <- rbinom(n=pop.size, prob=1/(1+exp(-mn)), size=1)
  
  ## generate the size of surviving individuals next year (subset z so that only surviving 
  ## individuals are included)
  i.subset <- which(surv == 1)
  z1 <- rep(NA, pop.size)
  E.z1 <- m.par["grow.int"] + m.par["grow.z"] * z[i.subset] + yref[i, "grow"]
  z1[i.subset] <- rnorm(n = length(i.subset), mean = E.z1, sd = m.par["grow.sd"])
  
  ## generate a binomial random number for reproduction from surviving individuals
  repr <- rep(NA, pop.size)
  mn <- m.par["rep.int"] + m.par["rep.z"] * z[i.subset] + yref[i,"rep"]
  repr[i.subset] <- rbinom(n = length(i.subset), prob = 1/(1+exp(-mn)), size=1)
  
  ## generate a binomial random number for offspring sex (female==1)
  i.subset <- which(surv == 1 & repr == 1)
  osex <- rep(NA, pop.size)
  osex[i.subset] <- rbinom(n = length(i.subset), prob = 1/2, size=1)
  
  ## generate a binomial random number for offspring recruitment from surviving / reproducing individuals
  i.subset <- which(surv == 1 & repr == 1 & osex == 1)
  recr <- rep(NA, pop.size)
  recr[i.subset] <- rbinom(n = length(i.subset), prob=m.par["recr.int"], size=1)
  
  ## generate the size of new recruits
  i.subset <- which(surv == 1 & repr == 1 & osex==1 & recr == 1)
  z1.rec <- rep(NA, pop.size)
  E.rec.z1 <- m.par["rcsz.int"] + m.par["rcsz.z"] * z[i.subset] + yref[i, "recsz"]
  z1.rec[i.subset] <- rnorm(n = length(i.subset), mean = E.rec.z1, sd = m.par["rcsz.sd"])
  
  ## create new population body size vector, store population size 
  zall <- c(z1.rec[which(surv == 1 & repr == 1 & osex==1 & recr == 1)], z1[which(surv == 1)])
  if(yr>1) stoch.lamb[yr] <- log(length(zall)/500)
  pop.size.t[yr] <- length(zall)
  
  ## store the simulation data for the current year
  sim.data[[i]] <- data.frame(z, surv, z1, repr, osex, recr, z1.rec, year=i, popsize=length(zall))
  
  ## take sample of 500 from current pop size
  z <- zall[sample.int(length(zall), 500, replace=T)]
  
  ## iterate the year
  yr <- yr+1
  
}

## bind list of populations
sim.data <- bind_rows(sim.data)

## reassign column names to match jags models
names(sim.data) <- c("z","Surv","z1","Repr","Sex","Recr","Rcsz", "Year", "popsize")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Extract data set of realistic length for demographic study
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Extract years
start.yr <- sample(run.in:(tot.yrs-n.yrs), 1)
sim.data <- sim.data %>% filter(Year %in% start.yr:(start.yr+(n.yrs-1))) %>%
  mutate(Year = Year-start.yr+1)

## Extract individuals
samp <- data.frame(Year = unique(sim.data$Year), Ind = sample(n.ind, n.yrs))

sim.data <- sim.data %>% nest(-Year) %>% 
  left_join(samp, by = "Year") %>%
  mutate(Sample = map2(data, Ind, sample_n)) %>%
  unnest(Sample)
