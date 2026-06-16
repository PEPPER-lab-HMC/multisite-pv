#### Wrangling and plotting PV curve data for SRER ####
# Incorporate dry weights to get RWC

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

# Read in output of 3a
larrea_df <- read_csv("data_clean/02-cleaning-outputs/02-SRER-cleaning-output.csv")

# create dataframe of dry weights and calculate RWC
dry <- data.frame(ID = 1:11,
                  dry_mass = c(2.8602,
                               1.0059,
                               0.9365,
                               0.9159,
                               0.6535,
                               2.1702,
                               1.0052,
                               1.6909,
                               1.5351,
                               1.4964,
                               1.9769))

# combined with dry mass and calculate water content
comb <- larrea_df |> 
  left_join(dry, by = "ID") |> 
  mutate(wc = (mass.g - dry_mass) / dry_mass)

comb |> 
  ggplot(aes(x = P.MPa, y = mass.g,
             color = sample)) +
  geom_point() +
  geom_line() +
  theme_bw()

#### Calculating saturated mass and rwc ####
# Optimize r^2 on the mass ~ water_pot_mpa lm
larrea_r2_list <- list()
temp_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("r2", "n", "slope", "intercept"))
IDs <- c(1:11)
for(i in 1:length(IDs)) {
  current_comb <- comb |> filter(ID == i) |> arrange(desc(mass.g))
  for(j in 1:(length(current_comb$P.MPa) - 2)) {
    cur_model <- current_comb |> 
      mutate(fitkeep = FALSE, fitkeep = replace(fitkeep, 1:(j + 2), TRUE)) |> 
      filter(fitkeep == TRUE)
    
    cur_lm <- lm(mass.g ~ P.MPa, data = cur_model)
    
    temp_df <- temp_df |> add_row(r2 = summary(cur_lm)$r.squared, n = (j + 2), slope = coef(cur_lm)[2], intercept = coef(cur_lm)[1])
  }
  temp_df <- temp_df |> 
    mutate(ID = IDs[i]) |> 
    slice_max(order_by = r2, n = 1)
  
  larrea_r2_list[[i]] <- bind_rows(temp_df)
  temp_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("r2", "n", "slope", "intercept"))
}
larrea_r2 <- do.call(bind_rows, larrea_r2_list) |> 
  group_by(ID) |> 
  summarize(
    r2 = max(r2),
    n = max(n),
    slope = max(slope),
    intercept = max(intercept),
    ID = max(ID)
  )

comb1.2 <- comb |> 
  left_join(larrea_r2, by = "ID") |> 
  group_by(ID) |> 
  mutate(keep = FALSE, keep = replace(keep, 1:n, TRUE))

comb1.2 |> 
  ggplot(aes(x = P.MPa, y = mass.g,
             color = keep, group = ID)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0) +
  facet_wrap(~ ID)

# Match to original
comb2 <- comb1.2 |> 
  rename(sat_mass_est = intercept) |>  # estimate saturated mass
  mutate(rwc = (mass.g - dry_mass) / (sat_mass_est - dry_mass),
         rwc_plotting = 1 - rwc)

comb2 |> 
  ggplot(aes(x = 1-rwc, y = 1/P.MPa,
             color = sample)) +
  geom_point() +
  geom_line() +
  scale_x_continuous("1 - RWC") +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")")))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.70))

# Save out only TRUE points
comb2 |> 
  write_csv("data_clean/03-pv-curve-data/SRER_pv_curve.csv")
