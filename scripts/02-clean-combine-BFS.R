#### Wrangling and plotting PV curve data for sagebrush ####

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

#### Reading in sagebrush data ####
# Begin by concatenating data together and initial plotting

# Double for loop
# first through folders (dates)
# then by shrubID
pv_dates <- list.files("data_clean/BFS_pv/")
pv_all <- list()
for(i in 1:length(pv_dates)) {
  fns <- list.files(paste0("data_clean/BFS_pv/", pv_dates[i], "/"))
  pv_date <- c()
  for(j in 1:length(fns)) {
    temp <- read_csv(paste0("data_clean/BFS_pv/", pv_dates[i], "/", fns[j])) |> 
      mutate(P.MPa = -1 * P.MPa,
             date = as.Date(pv_dates[i]),
             shrubID = str_extract(fns[j], pattern = "\\d{1,3}") |> 
               as.integer()) |> 
      relocate(date, shrubID) |> 
      select(-time)
    
    pv_date <- pv_date |> bind_rows(temp)
  }
  pv_all[[i]] <- pv_date
}
# Combine output of for loop into single df
sageb_df <- do.call(bind_rows, pv_all) |> 
  rename(P_mpa = P.MPa, mass_g = mass.g, leaf_mass_g = leaf.mass.g, offset_mass_g = offset.mass.g) |> 
  filter(keep == TRUE) |> 
  # rename(water_pot_mpa = p_m_pa) |>
  dplyr::select(-keep, -eTLP, -xTLP, -a, -b, -c, -d) |> 
  group_by(date, shrubID) |> 
  mutate(start_weight = max(mass_g), start_wp = max(P_mpa)) |> 
  relocate(start_weight, .after = mass_g) |> 
  mutate(mass_lost_g = abs(mass_g - start_weight),
         date_char = as.character(date), 
         week = case_when(date_char == "2025-06-04" ~ 2,
                          date_char == "2025-06-11" ~ 3,
                          date_char == "2025-06-17" ~ 4,
                          date_char == "2025-06-24" ~ 5,
                          date_char == "2025-07-01" ~ 6,
                          date_char == "2025-07-08" ~ 7,
                          date_char == "2025-07-15" ~ 8,
                          date_char == "2025-07-22" ~ 9)) |> 
  mutate(state = case_when(shrubID %in% c(5, 16, 46) & week == 4 ~ "Rehydrated",
                           shrubID %in% c(69, 122, 39) & week == 4 ~ "Fresh - Wet paper towel",
                           shrubID %in% c(90, 51, 101) & week == 5 ~ "Rehydrated",
                           shrubID %in% c(86, 56, 91) & week == 5 ~ "Fresh",
                           shrubID %in% c(5, 16, 46) & week == 6 ~ "Fresh - Wet paper towel",
                           shrubID %in% c(69, 122, 39) & week == 6 ~ "Rehydrated",
                           shrubID %in% c(90, 51, 101) & week == 7 ~ "Fresh",
                           shrubID %in% c(86, 56, 91) & week == 7 ~ "Rehydrated",
                           shrubID %in% c(69, 122, 39) & week == 8 ~ "Fresh",
                           shrubID %in% c(5, 16, 46) & week == 8 ~ "Rehydrated",
                           shrubID %in% c(90, 51, 101) & week == 7 ~ "Rehydrated",
                           shrubID %in% c(86, 56, 91) & week == 7 ~ "Fresh",
                           week == 2 ~ "Fresh - Wet paper towel",
                           week == 3 ~ "Fresh - Wet paper towel",
                           .default = "Fresh")) |> 
  dplyr::select(-date_char) |> 
  ungroup()

#### Plotting PV curves ####
# Plotting inverse water pot over mass lost
sageb_df |> 
  ggplot(aes(group = interaction(shrubID, week), x = mass_lost_g, y = -1/P_mpa, color = as.factor(shrubID))) + 
  geom_point() +
  geom_line()

sageb_df |>
  ggplot(aes(group = interaction(shrubID, date), x = mass_lost_g, y = -1/P_mpa, color = date)) +
  scale_color_gradient(low = "coral", high = "navyblue") +
  # scale_shape_manual("Week", values = c(15, 16, 17, 18, 5, 13, 9)) +
  geom_point(size = 3) +
  geom_line()

# Plot one at a time
# sageb_df |> filter(shrubID == 39 & week == 8) |>
#   ggplot(aes(x = P_mpa, y = mass_g)) +
#   geom_line(color = "cornflowerblue", linewidth = 0.75) +
#   geom_point(color = "cornflowerblue", size = 3)

# Throw out plants here if needed
sageb_df <- sageb_df |>
  mutate(rem = FALSE) |>
  mutate(rem = if_else(shrubID == 69 & week == 8, TRUE, rem)) |>
  filter(rem != TRUE) |>
  dplyr::select(-rem)

# Clean data for WP measurements taken too close together
to_remove <- sageb_df |> 
  group_by(shrubID, week) |> 
  mutate(lag_wp  = lag(P_mpa),
         diff_wp = lag_wp - P_mpa) |> 
  ungroup() |> 
  filter(diff_wp <= 0.06) 

sageb_df2 <- sageb_df |> 
  anti_join(to_remove)

to_remove <- sageb_df2 |> 
  group_by(shrubID, week) |> 
  arrange(desc(mass_g)) |> 
  slice(1:4) |> 
  mutate(next_wp = lead(P_mpa),
         diff_wp = next_wp - P_mpa) |> 
  ungroup() |> 
  filter(diff_wp >= 0)

sageb_df3 <- sageb_df2 |> 
  anti_join(to_remove) |> 
  group_by(shrubID, week) |> 
  mutate(lag_wp  = lag(P_mpa),
         diff_wp = lag_wp - P_mpa) |> 
  ungroup()
  
# Plotting without the removed points
sageb_df3 |> 
  ggplot(aes(x = mass_lost_g, y = -1/P_mpa, color = as.factor(shrubID))) + 
  # scale_shape_manual("Week", values = c(15, 16, 17, 18, 5, 13, 9)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(expression(paste("Mass lost (", H[2], "O lost)"))) +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) +
  scale_color_discrete("ID") +
<<<<<<< Updated upstream
  facet_wrap(~ week, scales = "free") +
=======
  facet_wrap(~ week) +
>>>>>>> Stashed changes
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15))

write_csv(sageb_df3, "data_clean/02-BFS-cleaning-output.csv")



