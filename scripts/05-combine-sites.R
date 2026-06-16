#### Combining data from SRER, ONAQ, and CdM ####

library(tidyverse)
library(janitor)
library(broom)
theme_set(theme_bw())

#### Reading in data sets ####

cdm_df <- read_csv("data_clean/03-pv-curve-data/CdM_pv_curve.csv") |> 
  select(-...1) |> 
  mutate(site = "CdM", water_pot_mpa = abs(water_pot_mpa))

onaq_df <- read_csv("data_clean/03-pv-curve-data/ONAQ_pv_curve.csv") |> 
  select(-...1) |> 
  mutate(site = "ONAQ", water_pot_mpa = abs(water_pot_mpa))

srer_df_1 <- read_csv("data_clean/03-pv-curve-data/pv_comb_20240308.csv") |> 
  rename(water_pot_mpa = P.MPa, id = ID, weight_g = mass.g, mass_lost_g = mass_lost, notes = note) |> 
  select(-sample, -leaf.mass.g, -offset.mass.g, -keep) |> 
  mutate(site = "SRER", total_weight_g = weight_g)         # Keeping the weight_g column as the total_weight_g column

srer_df_2 <- read_csv("data_clean/03-pv-curve-data/pv_comb_20231030.csv") |> 
  rename(water_pot_mpa = P.MPa, id = ID, weight_g = mass.g, mass_lost_g = mass_lost, notes = note) |> 
  select(-sample, -leaf.mass.g, -offset.mass.g, -keep) |> 
  mutate(site = "SRER", total_weight_g = weight_g)         # Keeping the weight_g column as the total_weight_g column

srer_start_1 <- data.frame(id = 1:6,
                         start_weight = c(4.7327,
                                          2.2643,
                                          3.4973,
                                          3.2939,
                                          3.3954,
                                          4.7907))
srer_df_1 <- srer_df_1 |> 
  left_join(srer_start_1, by = "id")

#### Joining into one df ####

multisite_df <- do.call("rbind", list(onaq_df, cdm_df, srer_df_1)) |> 
  relocate(site, id, water_pot_mpa, start_weight, weight_g, total_weight_g, mass_lost_g, dry_mass, wc, sat_mass_est, rwc, notes)
  
write_csv(multisite_df,
          file = "data_clean/multisite_pv_data.csv")

