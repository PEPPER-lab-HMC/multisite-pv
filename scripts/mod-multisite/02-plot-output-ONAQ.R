# Plotting

library(tidyverse)
library(coda)
library(broom.mixed)

# Load data
pv <- read_csv("data_clean/03-pv-curve-data/ONAQ_pv_curve.csv")

# Load coda
load(file = "scripts/mod-multisite/coda/coda_mod2b-ONAQ.Rdata")

# Summarize coda output
coda_sum <- tidyMCMC(jm_coda,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

# Isolate mean.tlp
tlp <- coda_sum |> 
  filter(term == "mean.tlp") |>
  mutate(lab = paste0("bold(Psi[TLP] == ", round(pred.mean, 2), ")"))

# write_csv(tlp, file = "data_clean/JAGS-model-outputs/ONAQ-model-outputs/ONAQ_model_mean_tlp.csv")

# Untransformed y
ggplot(pv, ) +
  geom_point(data = pv,
             aes(x = rwc_plotting,
                 y = -1/water_pot_mpa,
                 color = factor(id))) +
  geom_line(data = pv,
            aes(x = rwc_plotting,
                y = -1/water_pot_mpa,
                color = factor(id))) +
  geom_rect(data = tlp,
            aes(xmin = -Inf, xmax = Inf,
                ymin = -1/pred.upper, ymax = -1/pred.lower),
            alpha = 0.5,
            fill = "gray70") +
  geom_hline(data = tlp, aes(yintercept = -1/pred.mean),
             lty = "dashed",
             size = 1) +
  geom_text(data = tlp,
            aes(x = 0.2, y = 0.85, label = lab),
            hjust = 0,
            size = 5,
            parse = TRUE) +
  scale_x_continuous("1 - RWC (%)") +
  scale_y_continuous(expression(paste("1/", Psi[leaf], " (-MPa)"))) +
  theme_bw(base_size = 14) +
  guides(color = "none") +
  theme(panel.grid.minor = element_blank())
  
