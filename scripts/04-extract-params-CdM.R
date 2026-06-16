### Extracting PV parameters
library(tidyverse)
library(janitor)
library(broom)
library(here)
library(stringr)
library(purrr)
theme_set(theme_bw())

#### Read in model predicted TLPs ####

pred <- read_csv("data_clean/CdM-model-outputs/CdM_model_pred.csv") |> 
  select(-term)
cps <- read_csv("data_clean/CdM-model-outputs/CdM_model_cps.csv") |>
  rename(cp_value = pred.median, cp_std.error = std.error, cp_pred.lower = pred.lower, cp_pred.upper = pred.upper)
tlps <- read_csv("data_clean/CdM-model-outputs/CdM_model_tlps.csv") |>
  rename(tlp_value = pred.median, tlp_std.error = std.error, tlp_pred.lower = pred.lower, tlp_pred.upper = pred.upper, lab_x = x, lab_y = y)
ds <- read_csv("data_clean/CdM-model-outputs/CdM_model_ds.csv") |> 
  rename(d_value = pred.median, d_std.error = std.error, d_pred.lower = pred.lower, d_pred.upper = pred.upper)
cs <- read_csv("data_clean/CdM-model-outputs/CdM_model_cs.csv") |> 
  rename(c_value = pred.median, c_std.error = std.error, c_pred.lower = pred.lower, c_pred.upper = pred.upper)

line_params <- merge(ds, cs, by = c("ID")) |> 
  select(-Parameter.x, -Parameter.y)

#### Plotting parameters ####
# TLPs and regression lines
pred |> 
  ggplot() +
  geom_point(aes(x = rwc_plotting, y = y, color = "observed")) +
  geom_errorbar(aes(x = rwc_plotting, ymin = pred.lower, ymax = pred.upper,
                    color = "predicted"),
                alpha = 0.25) +
  geom_point(aes(x = rwc_plotting, y = pred.mean, color = "predicted"),
             alpha = 0.5) +
  geom_line(aes(x = rwc_plotting, y = pred.mean, color = "predicted")) +
  geom_rect(data = cps, aes(xmin = cp_pred.lower, xmax = cp_pred.upper,
                            ymin = -Inf, ymax = Inf), 
            alpha = 0.25) +
  geom_rect(data = tlps, aes(ymin = -1/tlp_pred.lower, ymax = -1/tlp_pred.upper,
                             xmin = -Inf, xmax = Inf), 
            alpha = 0.25) +
  geom_vline(data = cps, aes(xintercept = cp_value), lty = "dashed") +
  geom_hline(data = tlps, aes(yintercept = -1/tlp_value), lty = "dashed") +
  geom_abline(data = line_params, aes(intercept = d_value, slope = -1 * c_value),
              color = "cornflowerblue", linewidth = 0.75) +
  geom_vline(aes(xintercept = 0)) +
  facet_wrap(~ID, scales = "free_y") +
  scale_x_continuous(expression("1 - RWC (%)")) +
  scale_y_continuous(expression(paste("1/", Psi[leaf], " (-MPa)"))) +
  scale_color_manual(values = c("black", "coral")) +
  theme_bw(base_size = 12) +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.25))

# TLP for each shrub
tlps |> 
  ggplot() +
  geom_errorbar(aes(x = ID, y = tlp_value, ymin = tlp_value - tlp_std.error,
                    ymax = tlp_value + tlp_std.error),
                width = 0, linewidth = 0.75, alpha = 0.7, color = "gray44") +
  geom_point(aes(x = ID, y = tlp_value), size = 5) +
  labs(x = "ID",
       y = expression(paste(Psi[TLP], "(MPa)"))) +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Osmotic potential and capacitance over time (facet ID)
line_params |> 
  ggplot() +
  geom_errorbar(aes(x = ID, y = d_value, ymin = d_value - d_std.error,
                    ymax = d_value + d_std.error),
                color = "cornflowerblue", alpha = 0.7, width = 0) +
  geom_point(aes(x = ID, y = d_value),
             color = "cornflowerblue", size = 3) +
  scale_y_continuous(expression(paste(pi[o])),
                     sec.axis = sec_axis(~ . + 0, name = "Capacitance")) +
  geom_errorbar(aes(x = ID, y = c_value, ymin = c_value - c_std.error,
                    ymax = c_value + c_std.error),
                color = "coral", alpha = 0.7, width = 0) +
  geom_point(aes(x = ID, y = c_value), color = "coral", size = 3) +
  # facet_wrap(~ shrub_id) +
  labs(x = "ID") +
  theme(axis.title.y = element_text(color = "cornflowerblue"),
        axis.title.y.right = element_text(color = "coral"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Meinzer 2014 figure 3, TLP over initial WP
pred |> 
  left_join(select(tlps, tlp_value, ID), by = c("ID")) |> 
  ggplot() +
  geom_point(aes(x = start_wp, y = tlp_value, color = ID), size = 5) +
  scale_color_gradient(low = "navyblue", high = "coral") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed") +
  xlim(-4, 0) +
  ylim(-4, 0) +
  labs(x = "Initial Water Potential (MPa)",
       y = expression(paste(Psi[TLP], "(MPa)")),
       color = "ID") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Meinzer 2014 figure 5, TLP over osmotic potential
tlps |> 
  left_join(select(ds, d_value, ID), by = c("ID")) |>
  mutate(d_value = -1 / d_value) |> 
  ggplot() +
  geom_point(aes(x = d_value, y = tlp_value, color = ID), size = 5) +
  scale_color_gradient(low = "navyblue", high = "coral") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed") +
  xlim(-3.5, -1) +
  ylim(-3.5, -1) +
  labs(x = expression(paste(pi[o]), ", MPa)"),
       y = expression(paste(Psi[TLP]), ", MPa)"),
       color = "ID") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))


