 
############ Functions for running IPM ###################

## Extract parameters from jags output ##
createpvec <- function(jagsmns, num.years){
  ## Set up matrix - with a column for each variable + row for each year
  nyrparam <- 15
  p.vec <- matrix(NA, nrow = num.years, ncol = nyrparam)
  colnames(p.vec) <- c("tLat", "tSu", "tGr", "tRep", "tRs", 
                       "IntSu", "IntGr", "IntRep", "IntRs",
                       "SzSu", "SzGr", "SzRep", "SzRs", 
                       "SdGr", "SdRs")
  # intercept variation due stochastic environment
  with(jagsmns, {
    p.vec[, "tLat"] <- rnorm(num.years, 0, 1)  # latent parameter (NB in JAGs model st dev of latent parameter is set to 1)
    p.vec[, "tSu"] <- rnorm(num.years, 0, sigma.ts)  # survival year effect
    p.vec[, "tGr"] <- rnorm(num.years, 0, sigma.tg)  # growth year effect
    p.vec[, "tRep"] <- rnorm(num.years, 0, sigma.trp)  # reproduction year effect
    p.vec[, "tRs"] <- rnorm(num.years, 0, sigma.tr)  # recruit size year effect
    # parameter estimates (intercepts and slopes)
    p.vec[, "IntSu"] <- b0su + bqsu * p.vec[,"tLat"] + p.vec[,"tSu"] # survival intercept
    p.vec[, "SzSu"] <- bzsu # survival size slope
    p.vec[, "IntGr"] <- b0gr + bqgr * p.vec[,"tLat"] + p.vec[,"tGr"] # growth intercept
    p.vec[, "SzGr"] <- bzgr # growth size slope
    p.vec[, "SdGr"] <- sigma.gr # growth sd
    p.vec[, "IntRep"] <- b0rep + bqrep * p.vec[,"tLat"] + p.vec[,"tRep"]
    p.vec[, "SzRep"] <- bzrep
    p.vec[, "IntRs"] <- b0rs + bqrs * p.vec[,"tLat"] + p.vec[,"tRs"]
    p.vec[, "SzRs"] <- bzrs
    p.vec[, "SdRs"] <- sigma.rs
    return(p.vec)
  })
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Vital rate functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Survival function, logistic regression
s_z <- function(z, m.par)
{
  linear.p <- m.par["IntSu"] + m.par["SzSu"] * z     # linear predictor
  p <- 1/(1+exp(-linear.p))                          # logistic transformation to probability
  return(p)
}

## Growth function, given you are size z at time t returns the probability density function of size z1 at t+1
g_z1z <- function(z1, z, m.par)
{
  mean <- m.par["IntGr"] + m.par["SzGr"] * z           # mean size next year
  sd <- m.par["SdGr"]                                  # sd about mean
  p.den.grow <- dnorm(z1, mean = mean, sd = sd)        # pdf that you are size z1 given you were size z
  return(p.den.grow)
}

## Reproduction function, logistic regression
rp_z <- function(z, m.par)
{
  linear.p <- m.par["IntRep"] + m.par["SzRep"] * z     # linear predictor
  p <- 1/(1+exp(-linear.p))                            # logistic transformation to probability
  return(p)
}

## Recruit size function, gaussian
rs_z1z <- function(z1, z, m.par)
{
  mean <- m.par["IntRs"] + m.par["SzRs"] * z           # linear predictor
  sd <- m.par["SdRs"]                                  # sd about mean
  p.den.grow <- dnorm(z1, mean = mean, sd = sd)        # pdf that offspring is size z1 given you were size z
  return(p.den.grow)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z, m.par) {
  return( s_z(z, m.par) * g_z1z(z1, z, m.par) )
}

## Define the reproduction kernel
F_z1z <- function (z1, z, m.par) {
  return(s_z(z, m.par) * rp_z(z, m.par) * (1/2) * rs_z1z(z1, z, m.par))
}
# NB: multiplied by 0.5 as single sex model i.e. assuming equal sex ratio

## Build the discretized kernel
mk_K <- function(m, m.par, L, U) {
  ## Calculate the 'implementation' parameters - m is the number of meshpoints and L and U are the lower 
  ## and upper limits respectively
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  ## Here P is the survival-growth kernel and F is the reproduction kernel
  P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
  F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
  K <- P + F
  return(K)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Functions iterate IPM and calculate population growth rate
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stochlamb <- function(jagsout, num.yrs, nrunin, m=100, L=0.89, U=4.70, init=c(rep(0,8), 500, rep(0,91))){
  # Extract parameters from jagsoutput
  m.par <- createpvec(jagsout, num.yrs)
  # Set up 
  nt <- init
  Rt <- numeric()
  # Iterate model
  for (t in 1:num.yrs){
    mat <- mk_K(m, m.par[t,], L, U)
    nt1 <- mat %*% nt
    Rt[t] <- log(sum(nt1))
    nt <- nt1/sum(nt1)
  }  
  # Calculate lambda
  stgr <- mean(Rt[nrunin:num.yrs], na.rm=T)
  return(stgr)
}
