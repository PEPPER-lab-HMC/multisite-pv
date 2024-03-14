### Control script to run model
# Changepoint exponential to linear, hierarchical with sample as RE
# Constrain mean parameters mu.a, mu.b, and mu.c to be positive 

library(tidyverse)
library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(broom.mixed)

# Load data
pv <- read_csv("data_clean/pmass_comb_20240308.csv") |>
  filter(keep == TRUE)

# Quick plot
# Log y
pv |> 
  ggplot(aes(x = mass_lost,
             y = log(1/P.MPa))) +
  geom_point(aes(color = factor(ID))) +
  geom_line(aes(color = factor(ID)))
# Untransformed y
pv |> 
  ggplot(aes(x = mass_lost,
             y = 1/P.MPa)) +
  geom_vline(xintercept = c(0.075, 0.4),
             lty = "dotted") +
  geom_point(aes(color = factor(ID))) +
  geom_line(aes(color = factor(ID)))

# Create list of data
dat_list <- list(y = 1/pv$P.MPa,
                 x = pv$mass_lost,
                 N = nrow(pv),
                 id = factor(pv$ID),
                 Nid = length(unique(pv$ID)))

# Create list of initials
inits <- function() {
  list(mu.log.a = rnorm(1, 0, 10),
       mu.log.b = rnorm(1, 0, 5),
       mu.log.c = runif(1, 0, 5),
       mu.cp = runif(1, 0.1, 0.35),
       tau = runif(1, 0, 1),
       sig.log.a = runif(1, 0, 1),
       sig.log.b = runif(1, 0, 1),
       sig.log.c = runif(1, 0, 1),
       sig.cp = runif(1, 0, 1))
}
inits_list <- list(inits(), inits(), inits())

# Or, load saved state
# load(saved_state, file = "scripts/mod-2/inits/inits_mod2.Rdata")

# Compile model
jm <- jags.model("scripts/mod-2/mod2.JAGS",
                 data = dat_list,
                 # inits = saved_state[[2]],
                 # inits = inits_list,
                 inits = list(saved_state[[2]][[3]],
                              saved_state[[2]][[3]],
                              saved_state[[2]][[2]]),
                 n.chains = 3)
update(jm, 10000)
dic.samples(jm, 10000)

# Monitor posterior samples
params <- c("deviance", "Dsum", "R2", # model fit parameters
            "a", "mu.a", "mu.log.a", "sig.log.a",
            "b", "mu.b", "mu.log.b", "sig.log.b",
            "c", "mu.c", "mu.log.c", "sig.log.c",
            "cp", "mu.cp", "sig.cp",
            "tau", "sig", # observation precision terms
            "tlp", "mu.tlp", "mean.tlp" # calculated turgor loss points
)
jm_coda <- coda.samples(jm, variable.names = params,
                        n.iter = 15000, thin = 15)

# Visualize chains
mcmcplot(jm_coda, parms = c("deviance", "Dsum", "R2", "sig",
                            "mu.log.a", "mu.a", 
                            "mu.log.b", "mu.b", 
                            "mu.log.c", "mu.c", 
                            "mu.cp",
                            "sig.log.a", "sig.log.b", "sig.log.c", "sig.cp",
                            "mu.tlp", "tlp", "mean.tlp"))

mcmcplot(jm_coda, parms = c("mu.a", "mu.b", "mu.c", "mu.cp", "mu.tlp", "mean.tlp"))

# mcmcplot(jm_coda, parms = c("a", "b", "c", "cp"))

caterplot(jm_coda, parms = c("a", "mu.a"), reorder = FALSE)
caterplot(jm_coda, parms = c("b", "mu.b"), reorder = FALSE)
caterplot(jm_coda, parms = c("c", "mu.c"), reorder = FALSE)
caterplot(jm_coda, parms = c("cp", "mu.cp"), reorder = FALSE)
caterplot(jm_coda, parms = c("tlp", "mean.tlp"), reorder = FALSE)

# Restart values
newinits <- initfind(jm_coda, OpenBUGS = FALSE)
newinits[[1]]
saved_state <- removevars(newinits, variables = c(1:10, 15:16, 22))
saved_state[[1]]
save(saved_state, file = "scripts/mod-2/inits/inits_mod2.Rdata")

int <- which(colnames(jm_coda[[1]]) == "Dsum")
mean(jm_coda[[1]][,int])
mean(jm_coda[[2]][,int])
mean(jm_coda[[3]][,int])


# Check estimates
a = 0.8
b = -4
c = 0.25
cp = 0.2

x <- seq(0, 1.5, by = 0.01)
y1 <- a*exp(b*x)
y2 <- -1*c*x +(a*exp(b*cp) - -1*c*cp)

plot(pv$mass_lost, 1/pv$P.MPa)
points(x, y1, col = "coral")
points(x, y2, col = "skyblue")
