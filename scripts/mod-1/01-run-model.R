## Testing sandbox
a = 1.5
b = 8
c = 0.3
x = seq(0, 1, 0.05)
y1 = a*exp(-b*x) + c


d = -0.2
cp = 0.6
y2 = c + d*(x-cp)
plot(x, y1, ylim = c(0, 1.3))
points(x, y2, col = "red")
abline(h = 0.3)
abline(v = 0.6)

### Control script to run model
library(tidyverse)
library(rjags)
library(mcmcplots)

pv <- read_csv("data_clean/pv_comb.csv")

dat_list <- list(invP = 1/pv$P.MPa,
                 x = pv$mass_lost,
                 N = nrow(pv),
                 id = pv$ID,
                 Nid = length(unique(pv$ID)))

inits <- function(){
  list(mu.a = rnorm(1, 0, 10),
       mu.b = rnorm(1, 0, 10),
       mu.c = rnorm(1, 0, 10),
       mu.d = rnorm(1, 0, 10),
       tau = runif(1, 0, 1),
       sig.a = runif(1, 0, 1),
       sig.b = runif(1, 0, 1),
       sig.c = runif(1, 0, 1),
       sig.d = runif(1, 0, 1))
}
inits_list <- list(inits(), inits(), inits())

# Run model
jm <- jags.model("scripts/mod-1/simple_v1.jags",
                 data = dat_list,
                 inits = inits_list,
                 n.chains = 3)

update(jm, n.iter = 100000)

params <- c("a", "b", "c", "d",
            "mu.a", "mu.b", "mu.c", "mu.d",
            "tau", "sig", 
            "sig.a", "sig.b", "sig.c", "sig.d")
jm_coda <- coda.samples(jm,
                        variable.names = params,
                        n.iter = 50000,
                        thin = 5)

mcmcplot(jm_coda, parms = c("a", "b", "c", "d",
                            "mu.a", "mu.b", "mu.c", "mu.d",
                            "tau", "sig", 
                            "sig.a", "sig.b", "sig.c", "sig.d"))

caterplot(jm_coda, parms = "a", reorder = FALSE)
caterplot(jm_coda, parms = "b", reorder = FALSE)
caterplot(jm_coda, parms = "c", reorder = FALSE)
caterplot(jm_coda, parms = "d", reorder = FALSE)
caterplot(jm_coda, regex = "mu", reorder = FALSE)
