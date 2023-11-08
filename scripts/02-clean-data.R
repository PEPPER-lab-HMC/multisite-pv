# loop through and combine raw data
library(tidyverse)
library(broom)

fn <- list.files("data_raw/")
fn_list<- list()
for(i in 1:length(fn)) {
  temp <- read_csv(paste0("data_raw/", fn[i])) |> 
    select(1:7) |> 
    mutate(sample = paste0("LATR_", i),
           ID = i,
           mass_init = mass.g[1],
           mass_lost = abs(mass.g - mass_init)) |> 
    relocate(sample, ID) |> 
    relocate(mass_lost, .after = mass.g) |> 
    select(-mass_init)
  fn_list[[i]] <- temp
}

fn_df <- do.call(rbind, fn_list)

# Select only TRUE points and plot together
fn_df |> 
  filter(keep == TRUE) |> 
  ggplot(aes(x = mass_lost, y = 1/P.MPa,
             color = sample)) +
  geom_point() +
  geom_line() +
  theme_bw()

# create dataframe of dry weights and calculate RWC
dry <- data.frame(ID = 1:5,
                  dry_mass = c(2.8602,
                               1.0059,
                               0.9365,
                               0.9159,
                               0.6535))

# combined with dry mass and calculate water content
comb <- fn_df |> 
  filter(keep == TRUE) |> 
  left_join(dry, by = "ID") |> 
  mutate(wc = (mass.g - dry_mass) / dry_mass)

comb |> 
  filter(P.MPa <= 3) |> 
  ggplot(aes(x = P.MPa, y = mass.g,
             color = sample)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Add lm for each sample to obtain intercepts
params_list <- comb |> 
  filter(P.MPa <= 3) |> 
  group_by(sample) |> 
  group_map(~ broom::tidy(lm(mass.g ~ P.MPa, data = .x)))

params <- do.call(rbind, params_list) |> 
  filter(term == "(Intercept)") |> 
  mutate(ID = 1:5)

# Match to original
comb2 <- comb |> 
  left_join(select(params, ID, estimate), by = "ID") |> 
  rename(sat_mass_est = estimate) |>  # estimate saturated mass
  mutate(rwc = (mass.g - dry_mass) / (sat_mass_est - dry_mass))

comb2 |> 
  ggplot(aes(x = 1-rwc, y = 1/P.MPa,
             color = sample)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Save out only TRUE points
fn_df |> 
  filter(keep == TRUE) |> 
  left_join(dry, by = "ID")
  write_csv("data_clean/pv_comb.csv")
