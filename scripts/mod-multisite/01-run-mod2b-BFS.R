### Control script to run model
# Changepoint exponential to linear, hierarchical with sample as RE
# Model does not run well along the whole real line, so keep parameter constraints
# Small group size, try folded T priors for precisions
# Set standard deviation among parameters as data

library(tidyverse)
library(rjags)
load.module('dic')
library(mcmcplots)
# devtools::install_github("fellmk/PostJAGS/postjags")
library(postjags)
library(broom.mixed)

# Load data
pv <- read_csv("data_clean/BFS_pv_curve2.csv") |>
  clean_names() |> 
  select(-x1) |> 
  rename(mass_lost = mass_lost_g, ID = id, P.MPa = water_pot_mpa) |> 
  mutate(P.MPa = abs(P.MPa)) |> 
  filter(!is.na(P.MPa))

# Quick plot
# Log y
pv |> 
  ggplot(aes(x = mass_lost,
             y = log(1/P.MPa))) +
  geom_vline(xintercept = c(0.1, 0.35),
             lty = "dotted") +
  geom_point(aes(color = factor(ID))) +
  geom_line(aes(color = factor(ID)))
# Untransformed y
pv |> 
  ggplot(aes(x = mass_lost,
             y = 1/P.MPa)) +
  geom_vline(xintercept = c(0.1, 0.35),
             lty = "dotted") +
  geom_point(aes(color = factor(ID))) +
  geom_line(aes(color = factor(ID)))

# Create list of data
dat_list <- list(y = 1/pv$P.MPa,
                 x = pv$mass_lost,
                 N = nrow(pv),
                 id = factor(pv$ID),
                 Nid = length(unique(pv$ID)),
                 Sa = 2,
                 Sb = 2,
                 Sc = 2,
                 Scp = 1)

# Create list of initials
inits <- function() {
  list(mu.log.a = rnorm(1, 0, 10),
       mu.log.b = rnorm(1, 0, 5),
       mu.log.c = runif(1, 0, 5),
       mu.cp = runif(1, 0.2, 1),
       tau = runif(1, 0, 1),
       tau.eps.a = runif(1, 0, 1),
       tau.eps.b = runif(1, 0, 1),
       tau.eps.c = runif(1, 0, 1),
       tau.eps.cp = runif(1, 0, 1))
}
inits_list <- list(inits(), inits(), inits())

# Or, load saved state
# load(file = "scripts/mod-multisite/inits/inits_mod2b-BFS2.Rdata") # BFS2 with chains 1 and 2 is best?
load(file = "scripts/mod-2/inits/inits_mod2b.Rdata")

# Compile model
jm <- jags.model("scripts/mod-multisite/mod2b-BFS.JAGS",
                 data = dat_list,
                 inits = saved_state[[2]],
                 # inits = inits_list,
                 # inits = list(saved_state[[2]][[1]],
                 #              saved_state[[2]][[2]],
                 #              saved_state[[2]][[1]]),
                 n.chains = 3)
update(jm, 10000)
# dic.samples(jm, 10000)

# Monitor posterior samples
params <- c("deviance", "Dsum", "R2", # model fit parameters
            "a", "mu.a", "mu.log.a", "sig.log.a", "tau.eps.a",
            "b", "mu.b", "mu.log.b", "sig.log.b", "tau.eps.b", 
            "c", "mu.c", "mu.log.c", "sig.log.c", "tau.eps.c",
            "cp", "mu.cp", "sig.cp", "tau.eps.cp",
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

mcmcplot(jm_coda, parms = c("mu.a", "mu.b", "mu.c", 
                            "mu.cp", "mu.tlp", "mean.tlp"))

# mcmcplot(jm_coda, parms = c("a", "b", "c", "cp", "tlp))

caterplot(jm_coda, parms = c("a", "mu.a"), reorder = FALSE)
caterplot(jm_coda, parms = c("b", "mu.b"), reorder = FALSE)
caterplot(jm_coda, parms = c("c", "mu.c"), reorder = FALSE)
caterplot(jm_coda, parms = c("cp", "mu.cp"), reorder = FALSE)
caterplot(jm_coda, parms = c("tlp", "mean.tlp"), reorder = FALSE)

# Restart values
# newinits <- initfind(jm_coda, OpenBUGS = FALSE)
# newinits[[1]]
# saved_state <- removevars(newinits, variables = c(1:10, 15:20, 26))
# saved_state[[1]]
# save(saved_state, file = "scripts/mod-multisite/inits/inits_mod2b-BFS4.Rdata")

# ind <- which(colnames(jm_coda[[2]]) == "Dsum")
# mean(jm_coda[[1]][,ind]) # pink
# mean(jm_coda[[2]][,ind]) # green
# mean(jm_coda[[3]][,ind]) # blue

# save(jm_coda, file = "scripts/mod-multisite/coda/coda_mod2b-BFS.Rdata")
# load(file = "scripts/mod-multisite/coda/coda_mod2b-BFS.Rdata")

# Check convergence
gel <- gelman.diag(jm_coda, multivariate = FALSE)
gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) | grepl("deviance", rowname) | grepl("R2", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("mu", rowname)) # all converged except mu.tlp

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
  filter(grepl("^c\\[", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^cp", rowname))

# Run model for replicated data, time series of sigma and lambda
coda_rep <- coda.samples(jm, 
                         variable.names = c("y.rep"),
                         n.iter = 15000,
                         n.thin = 15)

# save(coda_rep, file = "scripts/mod-multisite/coda/rep_mod2b-BFS.Rdata")
# load(file = "scripts/mod-multisite/coda/rep_mod2b-BFS.Rdata")

# Summarize replicated output
coda_sum <- tidyMCMC(coda_rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

# Check model fit
pred <- cbind.data.frame(pv, coda_sum) |> 
  mutate(y = 1/P.MPa)

m1 <- lm(pred.mean ~ y, data = pred)
summary(m1) # R2 = , slope = 

# Fit plot
pred |> 
  ggplot(aes(x = y, y =pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper),
                alpha = 0.25) +
  geom_point() +
  scale_x_continuous("Observed") +
  scale_y_continuous("Predicted") +
  theme_bw(base_size = 12) +
  facet_wrap(~ID) +
  coord_equal()

# Summarize parameters
param_sum <- tidyMCMC(jm_coda,
                      robust = TRUE, # If TRUE, then median and MAD
                      conf.int = TRUE,
                      conf.method = "HPDinterval") |> 
  rename(pred.median = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

cps <- param_sum |> 
  filter(grepl("^cp", term)) |> 
  tidyr::separate(term, 
                  into = c("Parameter", "ID")) |> 
  mutate(ID_temp = as.numeric(ID),
         ID = case_when(ID_temp == 1 ~ 5,
                        ID_temp == 2 ~ 16,
                        ID_temp == 3 ~ 39,
                        ID_temp == 4 ~ 46,
                        ID_temp == 5 ~ 69,
                        ID_temp == 6 ~ 122,
                        ID_temp == 7 ~ 56,
                        ID_temp == 8 ~ 86,
                        ID_temp == 9 ~ 90,
                        ID_temp == 10 ~ 91))

tlps <- param_sum |> 
  filter(grepl("^tlp", term)) |> 
  tidyr::separate(term, 
                  into = c("Parameter", "ID")) |> 
  mutate(ID_temp = as.numeric(ID),
         ID = case_when(ID_temp == 1 ~ 5,
                        ID_temp == 2 ~ 16,
                        ID_temp == 3 ~ 39,
                        ID_temp == 4 ~ 46,
                        ID_temp == 5 ~ 69,
                        ID_temp == 6 ~ 122,
                        ID_temp == 7 ~ 56,
                        ID_temp == 8 ~ 86,
                        ID_temp == 9 ~ 90,
                        ID_temp == 10 ~ 91),
         x = 0.1,
         y = 0.95,
         lab = paste0("TLP = ", round(pred.median, 3)))

# Curves
pred |> 
  ggplot() +
  geom_text(data = tlps,
            aes(x = x, y = y, label = lab),
            hjust = 0) +
  geom_point(aes(x = mass_lost, y = y, color = "observed")) +
  geom_errorbar(aes(x = mass_lost, ymin = pred.lower, ymax = pred.upper,
                    color = "predicted"),
                alpha = 0.25) +
  geom_point(aes(x = mass_lost, y = pred.mean, color = "predicted"),
             alpha = 0.5) +
  geom_line(aes(x = mass_lost, y = pred.mean, color = "predicted")) +
  geom_rect(data = cps, aes(xmin = pred.lower, xmax = pred.upper,
                            ymin = -Inf, ymax = Inf), 
            alpha = 0.25) +
  geom_rect(data = tlps, aes(ymin = -1/pred.lower, ymax = -1/pred.upper,
                             xmin = -Inf, xmax = Inf), 
            alpha = 0.25) +
  geom_vline(data = cps, aes(xintercept = pred.median), lty = "dashed") +
  geom_hline(data = tlps, aes(yintercept = -1/pred.median), lty = "dashed") +
  facet_wrap(~ID) +
  scale_x_continuous(expression(paste(H[2], "O lost (g)"))) +
  scale_y_continuous(expression(paste("1/", Psi[leaf], " (-MPa)"))) +
  scale_color_manual(values = c("black", "coral")) +
  theme_bw(base_size = 12) +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.25))
