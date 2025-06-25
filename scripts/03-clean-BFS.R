#### Wrangling and plotting PV curve data for juniper ####

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

#### Reading in sagebrush data, weeks 2 and 3 ####

ids <- paste0(c(5, 16, 69, 122, 46, 39))
sageb_list <- list()
for(i in 1:length(ids)) {
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-04/", ids[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 1) |> 
    relocate(id) |> 
    select(-time)
  sageb_list[[i]] <- sageb_import
}

ids2 <- paste0(c(101, 51, 56, 86, 90, 91))
for(i in 1:length(ids2)) {
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-11/", ids2[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids2[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 2) |> 
    relocate(id) |> 
    select(-time)
  sageb_list[[i + 6]] <- sageb_import
}

sageb_df <- do.call(bind_rows, sageb_list) |> 
  filter(keep == TRUE) |> 
  rename(water_pot_mpa = p_m_pa) |> 
  select(-keep, -e_tlp, -x_tlp, -a, -b, -c, -d)

# Calculate number of points for each bush, min and max WP
sageb_df |>
  group_by(id) |>
  summarize(n = n(),
            maxY = -1*min(water_pot_mpa),
            minY = -1*max(water_pot_mpa))

# Create mass_lost column
lost <- data.frame(c(ids, ids2),
                   start_weight = c(sageb_list[[1]]$mass_g[1],
                                    sageb_list[[2]]$mass_g[1],
                                    sageb_list[[3]]$mass_g[1],
                                    sageb_list[[4]]$mass_g[1],
                                    sageb_list[[5]]$mass_g[1],
                                    sageb_list[[6]]$mass_g[1],
                                    sageb_list[[7]]$mass_g[1],
                                    sageb_list[[8]]$mass_g[1],
                                    sageb_list[[9]]$mass_g[1],
                                    sageb_list[[10]]$mass_g[1],
                                    sageb_list[[11]]$mass_g[1],
                                    sageb_list[[12]]$mass_g[1])) |> 
  clean_names() |> 
  rename(id = c_ids_ids2)

sageb_df2 <- sageb_df |> 
  left_join(lost, by = "id") |> 
  mutate(mass_lost_g = abs(mass_g - start_weight))

# Plotting inverse water pot over mass lost
sageb_df2 |>
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id))) + 
  geom_point() +
  geom_line()

# Throw out plants here if needed
sageb_df2 <- sageb_df2 |> 
  filter(id != "51" & id != "101")

# Clean data for WP measurements taken too close together

# Plotting without the removed points


sageb_df2 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(week), shape = as.factor(id))) + 
  scale_shape_manual("Shrub", values = c(15, 16, 17, 18, 0, 1, 2, 5, 6, 4)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(expression(paste("Mass lost (", H[2], "O lost)"))) +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Concatenate dry weights

# Combine with dry mass and calculate water content

# Add lm for each bush to obtain intercepts
params_list <- sageb_df2 |> 
  group_by(id) |> 
  group_map(~ broom::tidy(lm(mass_g ~ water_pot_mpa, data = .x)))

params <- do.call(rbind, params_list) |> 
  filter(term == "(Intercept)") |> 
  mutate(id = c("5", "16", "69", "122", "46", "39", "56", "86", "90", "91"))

# Match to comb
comb2 <- sageb_df2 |> 
  left_join(select(params, id, estimate), by = "id") |> 
  rename(sat_mass_est = estimate)  # estimate saturated mass
  # mutate(rwc = (total_weight_g - dry_mass) / (sat_mass_est - dry_mass))

# PV curve (1/WP)
comb2 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa,
             color = factor(id))) +
  geom_point(size = 4) +
  geom_line(size = 0.75) +
  scale_x_continuous("Mass lost") +
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
             color = factor(id))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("Mass lost") +
  scale_y_continuous(expression(paste(Psi, " (-", MPa^-1, ")"))) + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

#Save this out here!!!
write.csv(comb2,
          file = "data_clean/BFS_pv_curve.csv")



#### Reading in sagebrush data, week 4 ####

ids <- paste0(c(5, 16, 69, 122, 46, 39))
sageb_list <- list()
for(i in 1:length(ids)) {
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-17/", ids[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 1) |> 
    relocate(id) |> 
    select(-time)
  sageb_list[[i]] <- sageb_import
}

sageb_df <- do.call(bind_rows, sageb_list) |> 
  filter(keep == TRUE) |> 
  rename(water_pot_mpa = p_m_pa) |> 
  select(-keep, -e_tlp, -x_tlp, -a, -b, -c, -d)

# Calculate number of points for each bush, min and max WP
sageb_df |>
  group_by(id) |>
  summarize(n = n(),
            maxY = -1*min(water_pot_mpa),
            minY = -1*max(water_pot_mpa))

# Create mass_lost column
lost <- data.frame(c(ids, ids2),
                   start_weight = c(sageb_list[[1]]$mass_g[1],
                                    sageb_list[[2]]$mass_g[1],
                                    sageb_list[[3]]$mass_g[1],
                                    sageb_list[[4]]$mass_g[1],
                                    sageb_list[[5]]$mass_g[1],
                                    sageb_list[[6]]$mass_g[1])) |> 
  clean_names() |> 
  rename(id = c_ids_ids2)

sageb_df2 <- sageb_df |> 
  left_join(lost, by = "id") |> 
  mutate(mass_lost_g = abs(mass_g - start_weight))

# Plotting inverse water pot over mass lost
sageb_df2 |>
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id))) + 
  geom_point(size = 3) +
  geom_line(linewidth = 0.75)

# Throw out plants here if needed
# sageb_df2 <- sageb_df2 |>
#   filter(id != "16")

# Clean data for WP measurements taken too close together
to_remove <- sageb_df2 |> 
  group_by(id) |> 
  mutate(lag_wp  = lag(water_pot_mpa),
         diff_wp = lag_wp - water_pot_mpa) |> 
  ungroup() |> 
  filter(diff_wp <= 0.06) 

sageb_df3 <- sageb_df2 |> 
  anti_join(to_remove)

# Plotting without the removed points

sageb_df3 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id))) + 
  geom_point(size = 3) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(expression(paste("Mass lost (", H[2], "O lost)"))) +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Concatenate dry weights

# Combine with dry mass and calculate water content

# Add lm for each bush to obtain intercepts
params_list <- sageb_df2 |> 
  group_by(id) |> 
  group_map(~ broom::tidy(lm(mass_g ~ water_pot_mpa, data = .x)))

params <- do.call(rbind, params_list) |> 
  filter(term == "(Intercept)") |> 
  mutate(id = c("5", "16", "69", "122", "46", "39"))

# Match to comb
comb2 <- sageb_df2 |> 
  left_join(select(params, id, estimate), by = "id") |> 
  rename(sat_mass_est = estimate)  # estimate saturated mass
# mutate(rwc = (total_weight_g - dry_mass) / (sat_mass_est - dry_mass))

# PV curve (1/WP)
comb2 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa,
             color = factor(id))) +
  geom_point(size = 4) +
  geom_line(size = 0.75) +
  scale_x_continuous("Mass lost") +
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
             color = factor(id))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("Mass lost") +
  scale_y_continuous(expression(paste(Psi, " (-", MPa^-1, ")"))) + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

#Save this out here!!!
write.csv(comb2,
          file = "data_clean/BFS_pv_curve2.csv")






#### Comparing rehydrated / fresh ####

ids <- paste0(c(5, 16, 69, 122, 46, 39))
sageb_list <- list()
for(i in 1:length(ids)) {
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-17/", ids[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 1) |> 
    relocate(id) |> 
    select(-time)
  sageb_list[[i]] <- sageb_import
}

ids2 <- paste0(c(101, 51, 56, 86, 90, 91))
for(i in 1:length(ids2)) {
  sageb_import <- read_csv(paste0("data_clean/BFS_pv/2025-06-24/", ids2[i], " pvdata.csv")) |>
    clean_names() |> 
    mutate(id = ids2[i],
           dt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
           p_m_pa = -1 * p_m_pa,
           week = 2) |> 
    relocate(id) |> 
    select(-time)
  sageb_list[[i + 6]] <- sageb_import
}

sageb_df <- do.call(bind_rows, sageb_list) |> 
  filter(keep == TRUE) |> 
  rename(water_pot_mpa = p_m_pa) |> 
  select(-keep, -e_tlp, -x_tlp, -a, -b, -c, -d)

lost <- data.frame(c(ids, ids2),
                   start_weight = c(sageb_list[[1]]$mass_g[1],
                                    sageb_list[[2]]$mass_g[1],
                                    sageb_list[[3]]$mass_g[1],
                                    sageb_list[[4]]$mass_g[1],
                                    sageb_list[[5]]$mass_g[1],
                                    sageb_list[[6]]$mass_g[1],
                                    sageb_list[[7]]$mass_g[1],
                                    sageb_list[[8]]$mass_g[1],
                                    sageb_list[[9]]$mass_g[1],
                                    sageb_list[[10]]$mass_g[1],
                                    sageb_list[[11]]$mass_g[1],
                                    sageb_list[[12]]$mass_g[1])) |> 
  clean_names() |> 
  rename(id = c_ids_ids2)

sageb_df2 <- sageb_df |> 
  left_join(lost, by = "id") |> 
  mutate(mass_lost_g = abs(mass_g - start_weight)) |> 
  mutate(state = case_when(id %in% c(5, 16, 46) ~ "Rehydrated",
                           id %in% c(69, 122, 39) ~ "Fresh",
                           id %in% c(90, 51, 101) ~ "Rehydrated",
                           id %in% c(86, 56, 91) ~ "Fresh")) |> 
  mutate(date = as.Date(dt), date = if_else(id %in% c(86, 90, 51, 56, 101, 91), as.Date("2025-06-24"), as.Date("2024-06-17")))

# Clean data for WP measurements taken too close together
to_remove <- sageb_df2 |> 
  group_by(id) |> 
  mutate(lag_wp  = lag(water_pot_mpa),
         diff_wp = lag_wp - water_pot_mpa) |> 
  ungroup() |> 
  filter(diff_wp <= 0.06) 

sageb_df3 <- sageb_df2 |> 
  anti_join(to_remove)

# Plotting inverse water pot over mass lost
sageb_df3 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(state), group = id)) + 
  geom_point(aes(shape = as.factor(date)), size = 3) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(expression(paste("Mass lost (", H[2], "O lost)"))) +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) +
  scale_color_discrete("State") +
  scale_shape_discrete("Week of") +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  facet_wrap(~ id, scales = "free")
