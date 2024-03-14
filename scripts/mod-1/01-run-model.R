### Control script to run model
library(tidyverse)
library(rjags)
library(mcmcplots)
library(postjags)
library(broom.mixed)

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

load("scripts/mod-1/inits/saved_state.Rdata")
# Run model
jm <- jags.model("scripts/mod-1/simple_v1.jags",
                 data = dat_list,
                 inits = saved_state[[2]],
                 n.chains = 3)

# update(jm, n.iter = 100000)
dic.samples(jm, n.iter = 3000, thin = 50)

params <- c("deviance", "Dsum", "R2_resid",
            "a", "b", "c", "d",
            "mu.a", "mu.b", "mu.c", "mu.d",
            "tau", "sig", 
            "sig.a", "sig.b", "sig.c", "sig.d")
jm_coda <- coda.samples(jm,
                        variable.names = params,
                        n.iter = 50000,
                        thin = 5)

# Save out
save(jm_coda, file = "scripts/mod-1/coda/jm_coda.Rdata")

mcmcplot(jm_coda, parms = c("deviance", "Dsum", "R2_resid",
                            "a", "b", "c", "d",
                            "mu.a", "mu.b", "mu.c", "mu.d",
                            "tau", "sig", 
                            "sig.a", "sig.b", "sig.c", "sig.d"))

caterplot(jm_coda, parms = "a", reorder = FALSE)
caterplot(jm_coda, parms = "b", reorder = FALSE)
caterplot(jm_coda, parms = "c", reorder = FALSE)
caterplot(jm_coda, parms = "d", reorder = FALSE)
caterplot(jm_coda, regex = "mu", reorder = FALSE)

ss_temp <- initfind(jm_coda, OpenBUGS = FALSE)
ss_temp[[1]]
saved_state <- removevars(ss_temp, variables = c(1:4, 9))
saved_state[[1]]

save(saved_state, file = "scripts/mod-1/inits/saved_state.Rdata")

# check convergence
gel <- gelman.diag(jm_coda, multivariate = FALSE)


gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) |grepl("deviance", rowname) | grepl("R2", rowname))
gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("mu", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^a", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^b", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^c", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^d", rowname))

# Posterior mean of Dsum and R2
coda_sum <- tidyMCMC(jm_coda,
                   conf.int = TRUE,
                   conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)
coda_sum %>%
  filter(grepl("Dsum", term) |
           grepl("R2_resid", term))


# Run replicated data
jm_rep <- coda.samples(jm,
                       variable.names = "invP.rep",
                       n.iter = 50000,
                       n.thin = 5)

# Save out
save(jm_rep, file = "scripts/mod-1/coda/jm_rep.Rdata")

rep_sum <- tidyMCMC(jm_rep,
                    conf.int = TRUE,
                    conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

pred <- cbind.data.frame(pv, rep_sum) |> 
  mutate(invP = 1/P.MPa)

m1 <- lm(pred.mean ~ invP, data = pred)
sm1 <- summary(m1) # R2 = 0.9959

ggplot(pred) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              linewidth = 1) +
  geom_abline(intercept = coef(sm1)[1,1], 
              
              slope = coef(sm1)[2,1], 
              col = "black",
              lty = 2) +
  geom_point(aes(x = invP,
                 y = pred.mean,
                 color = sample)) 
