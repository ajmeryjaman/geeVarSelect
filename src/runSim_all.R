##These codes run the simulations and store the results
rm(list = ls())

source("src/geeFUNCTIONS.R")
source("src/penalizedGEE.R")
source("src/geeBOOST.R")
source("src/simulation.R")

load("data/seeds.Rdata")

lam_max <- 1
lam_min <- 0.01*lam_max
lambda.seq <- sort(seq(lam_min,lam_max,(lam_max-lam_min)/99), decreasing=TRUE)

nsim <- 100

out <- vector(nsim, mode="list")
for(sim in 1:nsim){
  out[[sim]] <- simFUNnew(n=50, n.i=4, sigma2.e=1, alpha.true=0, maxitr.penGEE=100, maxitr.boosting=1000,
                          lambda.seq=lambda.seq , working.strs= c("independence", "exchangeable_fixed", "exchangeable"),
                          seed=seeds[sim], ntest = 500)
  # cat("Simulation:",sim,"\n") #,"Errors:",errors,"\n")
  save(out, file="results/out_rho0.Rdata")
}

out <- vector(nsim, mode="list")
for(sim in 1:nsim){
  out[[sim]] <- simFUNnew(n=50, n.i=4, sigma2.e=1, alpha.true=0.3, maxitr.penGEE=100, maxitr.boosting=1000,
                          lambda.seq=lambda.seq , working.strs= c("independence", "exchangeable_fixed", "exchangeable"),
                          seed=seeds[100+sim], ntest = 500)
  # cat("Simulation:",sim,"\n") #,"Errors:",errors,"\n")
  save(out, file="results/out_rho03.Rdata")
}

out <- vector(nsim, mode="list")
for(sim in 1:nsim){
  out[[sim]] <- simFUNnew(n=50, n.i=4, sigma2.e=1, alpha.true=0.7, maxitr.penGEE=100, maxitr.boosting=1000,
                          lambda.seq=lambda.seq , working.strs= c("independence", "exchangeable_fixed", "exchangeable"),
                          seed=seeds[200+sim], ntest = 500)
  # cat("Simulation:",sim,"\n") #,"Errors:",errors,"\n")
  save(out, file="results/out_rho07.Rdata")
}

out <- vector(nsim, mode="list")
for(sim in 1:nsim){
  out[[sim]] <- simFUNnew(n=50, n.i=4, sigma2.e=1, alpha.true=0.9, maxitr.penGEE=100, maxitr.boosting=1000,
                          lambda.seq=lambda.seq, working.strs= c("independence", "exchangeable_fixed", "exchangeable"),
                          seed=seeds[300+sim], ntest = 500)
  # cat("Simulation:",sim,"\n") #,"Errors:",errors,"\n")
  save(out, file="results/out_rho09.Rdata")
}
