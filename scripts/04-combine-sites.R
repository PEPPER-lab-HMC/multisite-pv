#### Combining data from SRER, ONAQ, and CdM ####

library(tidyverse)
library(janitor)
library(broom)
theme_set(theme_bw())

#### Reading in data sets ####

cdm_df <- read_csv("data_clean/juniper_pv_curve.csv") |> 
  select(-...1, -dt, -dry_weight_g, -total_dryweight_g) |> 
  mutate(site = "CdM")

onaq_df <- read_csv("data_clean/sage_pv_curve.csv") |> 
  select(-...1, -dt, -total_dryweight_g) |> 
  mutate(site = "ONAQ")

srer_df <- read_csv("data_clean/pv_comb_20240308.csv") |> 
  rename(water_pot_mpa = P.MPa, id = ID, weight_g = mass.g, mass_lost_g = mass_lost, notes = note) |> 
  select(-sample, -leaf.mass.g, -offset.mass.g, -keep) |> 
  mutate(site = "SRER", total_weight_g = weight_g)         # Keeping the weight_g column as the total_weight_g column

srer_start <- data.frame(id = 1:6,
                         start_weight = c(4.7327,
                                          2.2643,
                                          3.4973,
                                          3.2939,
                                          3.3954,
                                          4.7907))
srer_df <- srer_df |> 
  left_join(srer_start, by = "id")

#### Joining into one df ####

multisite_df <- do.call("rbind", list(onaq_df, cdm_df, srer_df)) |> 
  relocate(site, id, water_pot_mpa, start_weight, weight_g, total_weight_g, mass_lost_g, dry_mass, wc, sat_mass_est, rwc, notes)
  
