#### Wrangling and plotting PV curve data for CdM ####
# Incorporate dry weights to get RWC

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

# Read in output of 3a
<<<<<<< Updated upstream
juniper_df3 <- read_csv("data_clean/02-cleaning-outputs/02-CdM-cleaning-output.csv") |> 
=======
juniper_df3 <- read_csv("data_clean/02-CdM-cleaning-output.csv") |> 
>>>>>>> Stashed changes
  mutate(superID = consecutive_id(id))

# Concatenate dry weights
dry <- data.frame(id = 1:7,
                  dry_mass = c(1.0693,  # total_dryweight
                               2.066,   # total_dryweight
                               4.5902,  # total_dryweight
                               1.4735,  # total_dryweight
                               1.8918,  # total_dryweight
                               4.6606,  # total_dryweight
                               0.9064)) # total_dryweight

# Combine with dry mass and calculate water content
comb <- juniper_df3 |> 
  left_join(dry, by = "id") |> 
  mutate(wc = (total_weight_g - dry_mass) / dry_mass) |> 
  select(-dry_weight_g)

comb |> 
  ggplot(aes(x = water_pot_mpa, y = total_weight_g,
             color = factor(superID))) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0)

comb |> 
  ggplot(aes(x = wc, y = -1/water_pot_mpa,
             color = as.factor(superID))) +
  geom_point() +
  geom_line()

#### Calculating saturated mass and rwc ####
# Optimize r^2 on the mass ~ water_pot_mpa lm
juniper_r2_list <- list()
temp_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("r2", "n", "slope", "intercept"))
superIDs <- c(1:6)
for(i in 1:length(superIDs)) {
  current_comb <- comb |> filter(superID == i) |> arrange(desc(total_weight_g))
  for(j in 1:(length(current_comb$water_pot_mpa) - 2)) {
    cur_model <- current_comb |> 
      mutate(fitkeep = FALSE, fitkeep = replace(fitkeep, 1:(j + 2), TRUE)) |> 
      filter(fitkeep == TRUE)
    
    cur_lm <- lm(total_weight_g ~ water_pot_mpa, data = cur_model)
    
    temp_df <- temp_df |> add_row(r2 = summary(cur_lm)$r.squared, n = (j + 2), slope = coef(cur_lm)[2], intercept = coef(cur_lm)[1])
  }
  temp_df <- temp_df |> 
    mutate(superID = superIDs[i]) |> 
    slice_max(order_by = r2, n = 1)
  
  juniper_r2_list[[i]] <- bind_rows(temp_df)
  temp_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("r2", "n", "slope", "intercept"))
}
juniper_r2 <- do.call(bind_rows, juniper_r2_list) |> 
  group_by(superID) |> 
  summarize(
    r2 = max(r2),
    n = max(n),
    slope = max(slope),
    intercept = max(intercept),
    superID = max(superID)
  )

comb1.2 <- comb |> 
  left_join(juniper_r2, by = "superID") |> 
  group_by(superID) |> 
  mutate(keep = FALSE, keep = replace(keep, 1:n, TRUE))

comb1.2 |> 
  ggplot(aes(x = water_pot_mpa, y = total_weight_g,
<<<<<<< Updated upstream
             color = keep, group = superID)) +
=======
             color = as.factor(keep), group = superID)) +
>>>>>>> Stashed changes
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0) +
  facet_wrap(~ superID)

comb1.3 <- comb1.2 |> filter(keep == TRUE)

# Match to comb
comb2 <- comb1.2 |> 
  rename(sat_mass_est = intercept) |> # estimate saturated mass
  mutate(rwc = (total_weight_g - dry_mass) / (sat_mass_est - dry_mass),
         rwc_plotting = 1 - rwc)

comb2 |> 
  ggplot(aes(x = water_pot_mpa, y = total_weight_g,
             color = as.factor(id))) +
  geom_point() +
  geom_line() +
  geom_abline(aes(slope = slope, intercept = sat_mass_est, color = as.factor(id))) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ id)

# PV curve (1/WP)
comb2 |> 
  ggplot(aes(x = 1-rwc, y = -1/water_pot_mpa,
             color = factor(id))) +
  geom_point(size = 4) +
  geom_line(size = 0.75) +
  scale_x_continuous("1 - RWC") +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# (WP)
comb2 |> 
  ggplot(aes(x = 1-rwc, y = water_pot_mpa,
             color = factor(id))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("1 - RWC") +
  scale_y_continuous(expression(paste(Psi, " (-", MPa^-1, ")"))) + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

#Save this out here!!!
write.csv(comb2,
<<<<<<< Updated upstream
          file = "data_clean/03-pv-curve-data/CdM_pv_curve.csv")
=======
          file = "data_clean/juniper_pv_curve.csv")
>>>>>>> Stashed changes


