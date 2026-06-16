# loop through and combine raw data
library(tidyverse)
library(broom)

# For initial trial 20231030, using rehydrated branchlets
fn <- list.files("data_raw/20231030")
fn_list<- list()
for(i in 1:length(fn)) {
  temp <- read_csv(paste0("data_raw/20231030/", fn[i])) |> 
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

fn_df <- do.call(rbind, fn_list) |> 
  select(-time) |> 
  mutate(state = "Rehydrated", date = "20231030")

# For spring campaign 20240308, no rehydration (predawns were above -2 MPa)
fn2 <- list.files("data_raw/20240308")
fn2_list<- list()
for(i in 1:length(fn2)) {
  temp <- read_csv(paste0("data_raw/20240308/", fn2[i])) |> 
    select(1:7) |> 
    mutate(sample = paste0("LATR_", i),
           ID = i,
           mass_init = mass.g[1],
           mass_lost = abs(mass.g - mass_init)) |> 
    relocate(sample, ID) |> 
    relocate(mass_lost, .after = mass.g) |> 
    select(-mass_init, -time)
  fn2_list[[i]] <- temp
}

fn2_df <- do.call(rbind, fn2_list) |> 
  mutate(sample = case_when(sample == "LATR_1" ~ "LATR_6",
                            sample == "LATR_2" ~ "LATR_7",
                            sample == "LATR_3" ~ "LATR_8",
                            sample == "LATR_4" ~ "LATR_9",
                            sample == "LATR_5" ~ "LATR_10",
                            sample == "LATR_6" ~ "LATR_11"),
         ID = case_when(ID == 1 ~ 6,
                        ID == 2 ~ 7,
                        ID == 3 ~ 8,
                        ID == 4 ~ 9,
                        ID == 5 ~ 10,
                        ID == 6 ~ 11),
         state = "Fresh", date = "20240308")

srer_df <- rbind(fn_df, fn2_df)

# Select only TRUE points and plot together
srer_df |> 
  filter(keep == TRUE) |> 
  ggplot(aes(x = mass_lost, y = 1/P.MPa,
             color = sample)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Clean data for WP measurements taken too close together
to_remove <- srer_df |> 
  group_by("ID") |> 
  mutate(lag_wp  = lag(P.MPa),
         diff_wp = lag_wp - P.MPa) |> 
  ungroup() |> 
  filter(diff_wp >= 0.06) 

larrea_df <- srer_df |> 
  anti_join(to_remove) |> 
  filter(keep == TRUE) |> 
  select(-keep)

# Plotting without the removed points
larrea_df |> 
  ggplot(aes(x = mass_lost, y = 1/P.MPa, color = as.factor(ID))) + 
  geom_point(size = 4) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(expression(paste("Mass lost (", H[2], "O lost)"))) +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

write_csv(larrea_df,
          "data_clean/02-cleaning-outputs/02-SRER-cleaning-output.csv")
