### Extracting PV parameters
library(tidyverse)
library(janitor)
library(broom)
library(here)
library(stringr)
library(purrr)
theme_set(theme_bw())

#### Read in model predicted TLPs ####
pred <- read_csv("data_clean/BFS-model-outputs/BFS_model_pred.csv") |> 
  select(-term)
cps <- read_csv("data_clean/BFS-model-outputs/BFS_model_cps.csv") |>
  rename(cp_value = pred.median, cp_std.error = std.error, cp_pred.lower = pred.lower, cp_pred.upper = pred.upper)
tlps <- read_csv("data_clean/BFS-model-outputs/BFS_model_tlps.csv") |>
  rename(tlp_value = pred.median, tlp_std.error = std.error, tlp_pred.lower = pred.lower, tlp_pred.upper = pred.upper, lab_x = x, lab_y = y)
ds <- read_csv("data_clean/BFS-model-outputs/BFS_model_ds.csv") |>
  rename(d_value = pred.median, d_std.error = std.error, d_pred.lower = pred.lower, d_pred.upper = pred.upper)
cs <- read_csv("data_clean/BFS-model-outputs/BFS_model_cs.csv") |>
  rename(c_value = pred.median, c_std.error = std.error, c_pred.lower = pred.lower, c_pred.upper = pred.upper)

line_params <- merge(ds, cs, by = c("ID", "week", "shrub_id", "superID")) |> 
  select(-Parameter.x, -Parameter.y)

#### Plotting parameters ####
# TLPs and regression lines
pred |> 
  arrange(desc(week)) |> 
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
  facet_wrap(~factor(superID, c(unique(pred$superID))), scales = "free_y") +
  scale_x_continuous(expression("1 - RWC (%)")) +
  scale_y_continuous(expression(paste("1/", Psi[leaf], " (-MPa)"))) +
  scale_color_manual(values = c("black", "coral")) +
  theme_bw(base_size = 12) +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.25))

# TLP for each shrub over time
tlps |> 
  left_join(select(pred, shrub_id, week, state), by = c("shrub_id", "week")) |> 
  ggplot() +
  geom_errorbar(aes(x = week, y = tlp_value, ymin = tlp_value - tlp_std.error,
                    ymax = tlp_value + tlp_std.error),
                width = 0, linewidth = 0.75, alpha = 0.7, color = "gray44") +
  geom_point(aes(x = week, y = tlp_value, shape = state), size = 5) +
  geom_line(aes(x = week, y = tlp_value), linewidth = 1) +
  facet_wrap(~ shrub_id) +
  labs(x = "Week",
       y = expression(paste(Psi[TLP], " (MPa)")),
       shape = "State") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# TLP over time for each state
tlps |> 
  left_join(select(pred, shrub_id, week, state), by = c("shrub_id", "week")) |> 
  ggplot() +
  geom_errorbar(aes(group = interaction(week, shrub_id), x = week, y = tlp_value,
                    ymin = tlp_value - tlp_std.error, ymax = tlp_value - tlp_std.error,
                    color = week),
                width = 0, linewidth = 2, alpha = 0.7) +
  geom_point(aes(x = week, y = tlp_value, shape = state, color = week), size = 5) +
  scale_color_gradient(low = "navyblue", high = "coral") +
  facet_wrap(~ state) +
  labs(x = "Week",
       y = expression(paste(Psi[TLP], "(MPa)")),
       shape = "State", 
       color = "Week") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Osmotic potential and capacitance over time (facet ID)
line_params |> 
  left_join(select(pred, shrub_id, week, state), by = c("shrub_id", "week")) |> 
  ggplot() +
  geom_errorbar(aes(x = week, y = d_value, ymin = d_value - d_std.error,
                    ymax = d_value + d_std.error),
                color = "cornflowerblue", alpha = 0.7, width = 0) +
  geom_point(aes(x = week, y = d_value, shape = state),
             color = "cornflowerblue", size = 3) +
  geom_line(aes(x = week, y = d_value),
            color = "cornflowerblue", linewidth = 0.75) +
  scale_y_continuous(expression(paste(pi[o])),
                     sec.axis = sec_axis(~ . + 0, name = "Capacitance")) +
  geom_errorbar(aes(x = week, y = c_value, ymin = c_value - c_std.error,
                    ymax = c_value + c_std.error),
                color = "coral", alpha = 0.7, width = 0) +
  geom_point(aes(x = week, y = c_value, shape = state), color = "coral", size = 3) +
  geom_line(aes(x = week, y = c_value), color = "coral", linewidth = 0.75) +
  facet_wrap(~ shrub_id) +
  labs(x = "Week", shape = "State") +
  theme(axis.title.y = element_text(color = "cornflowerblue"),
        axis.title.y.right = element_text(color = "coral"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Meinzer 2014 figure 3, TLP over initial WP
pred |> 
  left_join(select(tlps, tlp_value, shrub_id, week), by = c("shrub_id", "week")) |> 
  ggplot() +
  geom_point(aes(x = start_wp, y = tlp_value, color = week, shape = state), size = 5) +
  scale_color_gradient(low = "navyblue", high = "coral") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed") +
  xlim(-4, 0) +
  ylim(-4, 0) +
  labs(x = "Initial Water Potential (MPa)",
       y = expression(paste(Psi[TLP], "(MPa)")),
       shape = "State",
       color = "Week") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Adding chull for grouping visibility
state_polygon <- pred |>
  left_join(select(tlps, tlp_value, shrub_id, week), by = c("shrub_id", "week")) |> 
  group_by(state) |> 
  slice(chull(start_wp, tlp_value))

pred |> 
  left_join(select(tlps, tlp_value, shrub_id, week), by = c("shrub_id", "week")) |> 
  ggplot() +
  geom_polygon(data = state_polygon,
               aes(x = start_wp,y = tlp_value, fill = as.factor(state)),
               alpha = 0.2) +
  geom_point(aes(x = start_wp, y = tlp_value, color = week, shape = state), size = 5) +
  scale_color_gradient(low = "navyblue", high = "coral") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed") +
  xlim(-4, 0) +
  ylim(-4, 0) +
  labs(x = "Initial Water Potential (MPa)",
       y = expression(paste(Psi[TLP], "(MPa)")),
       shape = "State",
       color = "Week",
       fill = "State (Group)") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Meinzer 2014 figure 5, TLP over osmotic potential
tlps |> 
  left_join(select(ds, d_value, shrub_id, week), by = c("shrub_id", "week")) |>
  left_join(select(pred, state, shrub_id, week), by = c("shrub_id", "week")) |> 
  mutate(d_value = -1 / d_value) |> 
  ggplot() +
  geom_point(aes(x = d_value, y = tlp_value, color = week, shape = state), size = 5) +
  scale_color_gradient(low = "navyblue", high = "coral") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed") +
  xlim(-3.5, -1) +
  ylim(-3.5, -1) +
  labs(x = expression(paste(pi[o]), ", MPa)"),
       y = expression(paste(Psi[TLP]), ", MPa)"),
       color = "Week",
       shape = "State") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Adding chull for grouping visibility
state_polygon2 <- tlps |>
  left_join(select(ds, d_value, shrub_id, week), by = c("shrub_id", "week")) |>
  left_join(select(pred, state, shrub_id, week), by = c("shrub_id", "week")) |> 
  mutate(d_value = -1 / d_value) |>
  group_by(state) |> 
  slice(chull(d_value, tlp_value))

tlps |> 
  left_join(select(ds, d_value, shrub_id, week), by = c("shrub_id", "week")) |>
  left_join(select(pred, state, shrub_id, week), by = c("shrub_id", "week")) |> 
  mutate(d_value = -1 / d_value) |> 
  ggplot() +
  geom_polygon(data = state_polygon2,
               aes(x = d_value, y = tlp_value, fill = as.factor(state)),
               alpha = 0.2) +
  geom_point(aes(x = d_value, y = tlp_value, color = week, shape = state), size = 5) +
  scale_color_gradient(low = "navyblue", high = "coral") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed") +
  xlim(-3.5, -1) +
  ylim(-3.5, -1) +
  labs(x = expression(paste(pi[o]), ", MPa)"),
       y = expression(paste(Psi[TLP]), ", MPa)"),
       color = "Week",
       shape = "State",
       fill = "State (Group)") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Grouped state TLP
tlps |> 
  left_join(select(pred, shrub_id, week, state), by = c("shrub_id", "week")) |> 
  group_by(state, week) |> 
  summarize(
    tlp_m = mean(tlp_value),
    tlp_sd = sd(tlp_value),
    tlp_lower_m = mean(tlp_pred.lower),
    tlp_lower_sd = sd(tlp_pred.lower),
    tlp_upper_m = mean(tlp_pred.upper),
    tlp_upper_sd = sd(tlp_pred.upper),
    n = n()
  ) |> 
  ungroup() |> 
  ggplot() +
  geom_errorbar(aes(x = week, y = tlp_m, color = state, ymin = tlp_m - tlp_sd,
                    ymax = tlp_m + tlp_sd),
                width = 0, linewidth = 0.75, alpha = 0.7) +
  geom_point(aes(x = week, y = tlp_m, shape = state, color = state), size = 5) +
  geom_line(aes(x = week, y = tlp_m, group = state, color = state), linewidth = 1) +
  labs(x = "Week",
       y = expression(paste(Psi[TLP]), ", MPa)"),
       color = "State",
       shape = "State") +
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))


