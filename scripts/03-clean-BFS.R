#### Wrangling and plotting PV curve data for juniper ####

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

#### Reading in sagebrush data ####

ids <- paste0(c(5, 16, 69, 122, 46, 39))
sageb_list <- list()
for(i in 1:length(ids)) {
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-17/", ids[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 4) |> 
    relocate(id) |> 
    dplyr::select(-time)
  sageb_list[[i]] <- sageb_import
  
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-04/", ids[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 2) |> 
    relocate(id) |> 
    dplyr::select(-time)
  sageb_list[[i + length(ids)]] <- sageb_import
  
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-07-01/", ids[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 2) |> 
    relocate(id) |> 
    dplyr::select(-time)
  sageb_list[[i + length(ids)]] <- sageb_import
}

ids2 <- paste0(c(101, 51, 56, 86, 90, 91))
for(i in 1:length(ids2)) {
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-24/", ids2[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids2[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 5) |> 
    relocate(id) |> 
    dplyr::select(-time)
  sageb_list[[i + length(ids) + length(ids2)]] <- sageb_import
  
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-11/", ids2[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids2[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 3) |> 
    relocate(id) |> 
    dplyr::select(-time)
  sageb_list[[i + length(ids) + 2* length(ids2)]] <- sageb_import
}

sageb_df <- do.call(bind_rows, sageb_list) |> 
  filter(keep == TRUE) |> 
  rename(water_pot_mpa = p_m_pa) |> 
  dplyr::select(-keep, -e_tlp, -x_tlp, -a, -b, -c, -d)

lost_temp <- list()
for(i in 1:(2 * length(ids) + 2 * length(ids2))) {
  lost_temp[[i]] <- data.frame(id = sageb_list[[i]]$id[1], start_weight = sageb_list[[i]]$mass_g[1], week = sageb_list[[i]]$week[1])
}

lost <- do.call(bind_rows, lost_temp)

sageb_df2 <- sageb_df |> 
  left_join(lost, by = c("id", "week")) |> 
  mutate(mass_lost_g = abs(mass_g - start_weight)) |> 
  mutate(state = case_when(id %in% c(5, 16, 46) & week == 4 ~ "Rehydrated",
                           id %in% c(69, 122, 39) & week == 4 ~ "Fresh",
                           id %in% c(90, 51, 101) & week == 5 ~ "Rehydrated",
                           id %in% c(86, 56, 91) & week == 5~ "Fresh",
                           .default = "Fresh")) |> 
  mutate(date = as.Date(dt), date = case_when(week == 2 ~ "2025-06-04",
                                              week == 3 ~ "2025-06-11",
                                              week == 4 ~ "2025-06-17",
                                              week == 5 ~ "2025-06-24"))

# Plotting inverse water pot over mass lost
sageb_df2 |>
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id), shape = as.factor(week))) +
  scale_shape_manual("Week", values = c(15, 16, 17, 18)) +
  geom_point() +
  geom_line()

# Throw out plants here if needed
sageb_df2 <- sageb_df2 |>
  mutate(rem = FALSE) |>
  mutate(rem = if_else(id == 51 & week == 3, TRUE, rem)) |>
  filter(rem != TRUE) |>
  dplyr::select(-rem)

# Clean data for WP measurements taken too close together
to_remove <- sageb_df2 |> 
  group_by(id, week) |> 
  mutate(lag_wp  = lag(water_pot_mpa),
         diff_wp = lag_wp - water_pot_mpa) |> 
  ungroup() |> 
  filter(diff_wp <= 0.06) 

sageb_df3 <- sageb_df2 |> 
  anti_join(to_remove)

# Plotting without the removed points
sageb_df3 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id), shape = as.factor(week))) + 
  scale_shape_manual("Week", values = c(15, 16, 17, 18)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(expression(paste("Mass lost (", H[2], "O lost)"))) +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) +
  scale_color_discrete("ID") +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15))

# Concatenate dry weights
dry <- read_sheet("https://docs.google.com/spreadsheets/d/1KPhjMsLD-xkgNLjXhC9D4Js_FN3mPehTqr5W5n4inJg/edit?gid=0#gid=0") |> 
  mutate(id = as.character(id), date = as.character(date), week = as.integer(week))

# Combine with dry mass and calculate water content
comb <- sageb_df3 |> 
  left_join(dry, by = c("id", "week", "date")) |> 
  mutate(wc = (mass_g - dry_mass) / dry_mass)

comb |> 
  ggplot(aes(x = water_pot_mpa, y = mass_g,
             color = as.factor(id), shape = as.factor(week))) +
  scale_shape_manual("Week", values = c(15, 16, 17, 18)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0)

comb |> 
  ggplot(aes(x = wc, y = -1/water_pot_mpa,
             color = as.factor(id), shape = as.factor(week))) +
  scale_shape_manual("Week", values = c(15, 16, 17, 18)) +
  geom_point() +
  geom_line()

# Add lm for each bush to obtain intercepts
params_list <- comb |>
  group_by(id, week) |>
  group_map(~ broom::tidy(lm(mass_g ~ water_pot_mpa, data = .x)))

params <- do.call(rbind, params_list) |>
  filter(term == "(Intercept)") |>
  mutate(id = c(ids, ids, ids2, 101, 56, 86, 90, 91)) |>
  mutate(week = c(rep(c(4, 2, 5), each = 6), rep(3, times = 5)))

# Match to comb
comb2 <- comb |> 
  left_join(select(params, id, week, estimate), by = c("id", "week")) |> 
  rename(sat_mass_est = estimate) |> # estimate saturated mass
  mutate(rwc = (mass_g - dry_mass) / (sat_mass_est - dry_mass)) |> 
  mutate(holder = if_else(sat_mass_est > dry_mass, TRUE, FALSE))

# PV curve (1/WP)
comb2 |> 
  ggplot(aes(x = 1-rwc, y = -1/water_pot_mpa,
             color = as.factor(id), shape = as.factor(week))) +
  geom_point(size = 4) +
  geom_line(size = 0.75) +
  scale_x_continuous("(1-RWC)") +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# (WP)
comb2 |> 
  ggplot(aes(x = mass_lost_g, y = water_pot_mpa,
             color = factor(id), shape = as.factor(week))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("Mass lost") +
  scale_y_continuous(expression(paste(Psi, " (-", MPa^-1, ")"))) + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

#Save this out here!!!
write.csv(sageb_df3,
          file = "data_clean/BFS_pv_curve.csv")



