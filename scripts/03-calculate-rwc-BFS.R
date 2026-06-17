#### Wrangling and plotting PV curve data for sagebrush ####
# Incorporate dry weights to get RWC

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

# Read in output of 3a
sageb_df3 <- read_csv("data_clean/02-cleaning-outputs/02-BFS-cleaning-output.csv") |> 
  mutate(superID = consecutive_id(shrubID))

# Concatenate dry weights
dry <- read_sheet("https://docs.google.com/spreadsheets/d/1KPhjMsLD-xkgNLjXhC9D4Js_FN3mPehTqr5W5n4inJg/edit?gid=0#gid=0") |> 
  mutate(shrubID = as.integer(id), 
         date = as.Date(date),
         week = as.integer(week)) |> 
  relocate(shrubID, .after = date) |> 
  select(-id)

# Combine with dry mass and calculate water content
comb <- sageb_df3 |> 
  left_join(dry, by = c("shrubID", "week", "date")) |> 
  mutate(wc = (mass_g - dry_mass) / dry_mass)

comb |> 
  ggplot(aes(x = P_mpa, y = mass_g,
             color = as.factor(shrubID))) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0) +
  facet_wrap(~ week)

comb |> 
  ggplot(aes(x = wc, y = -1/P_mpa,
             color = as.factor(shrubID), shape = as.factor(week))) +
  scale_shape_manual("Week", values = c(15, 16, 17, 18, 5, 13, 9, 1)) +
  geom_point() +
  geom_line()

#### Calculating saturated mass and rwc ####

# Optimize r^2 on the mass ~ P_mpa lm
sageb_r2_list <- list()
temp_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("r2", "n", "slope", "intercept"))
superIDs <- c(1:47)
for(i in 1:length(superIDs)) {
  current_comb <- comb |> filter(superID == i) |> arrange(desc(mass_g))
  for(j in 1:(length(current_comb$P_mpa) - 2)) {
    cur_model <- current_comb |> 
      mutate(fitkeep = FALSE, fitkeep = replace(fitkeep, 1:(j + 2), TRUE)) |> 
      filter(fitkeep == TRUE)
    
    cur_lm <- lm(mass_g ~ P_mpa, data = cur_model)
    
    temp_df <- temp_df |>
      add_row(r2 = summary(cur_lm)$r.squared, n = (j + 2), slope = coef(cur_lm)[2], intercept = coef(cur_lm)[1])
  }
  temp_df <- temp_df |> 
    mutate(superID = superIDs[i]) |> 
    slice_max(order_by = r2, n = 1)
  
  sageb_r2_list[[i]] <- bind_rows(temp_df)
  temp_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("r2", "n", "slope", "intercept"))
}
sageb_r2 <- do.call(bind_rows, sageb_r2_list)

comb1.2 <- comb |> 
  left_join(sageb_r2, by = "superID") |> 
  group_by(superID) |> 
  mutate(keep = FALSE, keep = replace(keep, 1:n, TRUE))

comb1.2 |> 
  ggplot(aes(x = P_mpa, y = mass_g,
             color = as.factor(keep), group = interaction(shrubID, week))) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0) +
  facet_wrap(~ week)

comb1.3 <- comb1.2 |> filter(keep == TRUE)

# Match to comb
comb2 <- comb1.2 |> 
  rename(sat_mass_est = intercept) |> # estimate saturated mass
  mutate(rwc = (mass_g - dry_mass) / (sat_mass_est - dry_mass),
         rwc_plotting = 1 - rwc)

comb2 |> 
  ggplot(aes(x = P_mpa, y = mass_g,
             color = as.factor(shrubID))) +
  geom_point() +
  geom_line() +
  geom_abline(aes(slope = slope, intercept = sat_mass_est, color = as.factor(shrubID))) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ week)

# PV curve (1/WP)
comb2 |> filter(superID != 8) |> 
  ggplot(aes(x = rwc_plotting, y = -1/P_mpa,
             color = as.factor(week),
             shape = as.factor(state))) +
  geom_point(size = 4) +
  geom_line(size = 0.75) +
  scale_x_continuous("(1-RWC)") +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) + 
  facet_wrap(~ shrubID, scales = "free") +
  labs(color = "Week") +
  guides(shape = "none") +
  facet_wrap(~ shrubID, scales = "free") +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

comb2 |> filter(shrubID != 51 & week != 3) |> 
  ggplot(aes(group = interaction(shrubID, date), x = rwc_plotting, y = -1/P_mpa,
             color = date, shape = state)) +
  scale_color_gradient(low = "coral", high = "navyblue") +
  geom_point(size = 3) +
  geom_line(size = 0.75) +
  scale_x_continuous("(1-RWC)") +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) + 
  facet_wrap(~ state) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# (WP)
comb2 |> 
  ggplot(aes(x = rwc_plotting, y = P_mpa,
             color = factor(shrubID), shape = as.factor(state), group = interaction(shrubID, week))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("1 - RWC (%)") +
  scale_y_continuous(expression(paste(Psi, " (-", MPa^-1, ")"))) + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

comb2 |> 
  ggplot(aes(x = rwc, y = P_mpa,
             color = factor(shrubID), shape = as.factor(state), group = interaction(shrubID, week))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("RWC (%)") +
  scale_y_continuous(expression(paste(Psi, " (-", MPa^-1, ")"))) + 
  facet_wrap(~ shrubID) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

#Save this out here!!!
write.csv(comb2,
          file = "data_clean/03-pv-curve-data/BFS_pv_curve.csv")
