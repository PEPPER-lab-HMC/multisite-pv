# Plotting

library(tidyverse)
library(coda)
library(broom.mixed)

# Load data
pv <- read_csv("data_clean/pmass_comb_20240308.csv") |>
  filter(keep == TRUE)

# Load coda
load("scripts/mod-2/coda/coda_mod2b.Rdata")

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

# Untransformed y
ggplot(pv, ) +
  geom_point(data = pv,
             aes(x = mass_lost,
                 y = 1/P.MPa,
                 color = factor(ID))) +
  geom_line(data = pv,
            aes(x = mass_lost,
                y = 1/P.MPa,
                color = factor(ID))) +
  geom_rect(data = tlp,
            aes(xmin = -Inf, xmax = Inf,
                ymin = -1/pred.upper, ymax = -1/pred.lower),
            alpha = 0.5,
            fill = "gray70") +
  geom_hline(data = tlp, aes(yintercept = -1/pred.mean),
             lty = "dashed",
             size = 1) +
  geom_text(data = tlp,
            aes(x = 0.7, y = 0.85, label = lab),
            hjust = 0,
            size = 5,
            parse = TRUE) +
  scale_x_continuous(expression(paste(H[2], "O lost (g)"))) +
  scale_y_continuous(expression(paste("1/", Psi[leaf], " (-MPa)"))) +
  theme_bw(base_size = 14) +
  guides(color = "none") +
  theme(panel.grid.minor = element_blank())
  
