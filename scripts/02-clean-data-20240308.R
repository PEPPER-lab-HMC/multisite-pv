# loop through and combine raw data
library(tidyverse)
library(broom)

# For spring campaign 20240308, no rehydration (predawns were above -2 MPa)
fn <- list.files("data_raw/20240308")
fn_list<- list()
for(i in 1:length(fn)) {
  temp <- read_csv(paste0("data_raw/20240308/", fn[i])) |> 
    select(1:7) |> 
    mutate(sample = paste0("LATR_", i),
           ID = i,
           mass_init = mass.g[1],
           mass_lost = abs(mass.g - mass_init)) |> 
    relocate(sample, ID) |> 
    relocate(mass_lost, .after = mass.g) |> 
    select(-mass_init, -time)
  fn_list[[i]] <- temp
}

fn_df <- do.call(rbind, fn_list)

# Calculate number of true points per curve
fn_df |>
  filter(keep == TRUE) |>
  group_by(ID) |>
  summarize(n = n(),
            maxY = -1*min(P.MPa),
            minY = -1*max(P.MPa))

# Select only TRUE points and plot together
fn_df |> 
  filter(keep == TRUE) |>
  ggplot(aes(x = mass_lost, y = 1/P.MPa,
             color = factor(ID))) +
  geom_point() +
  geom_line() +
  theme_bw()

# Write out temporary for model

# write_csv(fn_df, file = "data_clean/pmass_comb_20240308.csv")

# create dataframe of dry weights and calculate RWC
dry <- data.frame(ID = 1:6,
                  dry_mass = c(2.1702,
                               1.0052,
                               1.6909,
                               1.5351,
                               1.4964,
                               1.9769))

# combined with dry mass and calculate water content
comb <- fn_df |> 
  filter(keep == TRUE) |> 
  left_join(dry, by = "ID") |> 
  mutate(wc = (mass.g - dry_mass) / dry_mass)

comb |> 
  # filter(P.MPa <= 3) |> 
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
  mutate(ID = 1:6)

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
  scale_x_continuous("1 - RWC") +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")")))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.70))

# Save out only TRUE points
comb2 |> 
  filter(keep == TRUE) |>
  write_csv("data_clean/pv_comb_20240308.csv")
