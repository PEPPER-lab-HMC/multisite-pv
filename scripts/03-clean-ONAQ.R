#### Wrangling and plotting PV curve data for sagebrush ####

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

#### Reading in sagebrush ####

ids <- paste0(c("2b", "3a", "4a", "5b", "6", "8"))
sage_list <- list()
for(i in 1:length(ids)) {
  sage_import <- read_sheet("https://docs.google.com/spreadsheets/d/1ycVFmL2Df4YFGp70lT7WdsF4Xn1dGllPtAnFksf5IHU/edit?usp=sharing",
                            sheet = paste0("Sage ", ids[i])) |>
    clean_names() |> 
    mutate(id = ids[i],
           time_c = as.character(time) |> 
             str_extract(pattern = "\\d\\d:\\d\\d"),
           hour = hm(time_c),
           date_c = if_else(hour < hour[1], "03/13/2025", "03/12/2025"),
           dt = as.POSIXct(paste(date_c, time_c), format = "%m/%d/%Y %H:%M", tz = "America/Denver")) |> 
    relocate(dt, id) |> 
    select(-time, -time_c, -hour, -date_c)
  
  sage_list[[i]] <- sage_import
}
sage_df <- do.call(bind_rows, sage_list) |> 
  select(-fallen_leaf_weight_g, -fallen_bark_dryweight_g, -fallen_leaf_dryweight_g, -fallen_leaf_2_weight_g, -fallen_leaf_3_weight_g, -fallen_leaf_4_weight_g, -fallen_leaf_5_weight_g, -dry_weight_g, -fallen_leaf_2_dryweight_g, -fallen_leaf_3_dryweight_g, -fallen_leaf_4_dryweight_g, -fallen_leaf_5_dryweight_g, -fallen_branch_dryweight_g, -fallen_branch_1_dryweight_g, -fallen_branch_2_dryweight_g)

# Calculate number of points for each bush, min and max WP
sage_df |>
  group_by(id) |>
  summarize(n = n(),
            maxY = -1*min(water_pot_mpa),
            minY = -1*max(water_pot_mpa))

# Create mass_lost column
lost <- data.frame(ids,
                   start_weight = c(1.7372,
                                    2.4743,
                                    1.7507,
                                    3.5279,
                                    4.9429,
                                    2.1356)) |> 
  rename("id" = "ids")

sage_df2 <- sage_df |> 
  left_join(lost, by = "id") |> 
  mutate(mass_lost_g = abs(total_weight_g - start_weight))

# Plotting inverse water pot over mass lost
sage_df2 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id))) + 
  geom_point() +
  geom_line()

# Clean data for WP measurements taken too close together
to_remove <- sage_df2 |> 
  group_by(id) |> 
  mutate(lag_wp  = lag(water_pot_mpa),
         diff_wp = lag_wp - water_pot_mpa) |> 
  ungroup() |> 
  filter(diff_wp <= 0.06) 

sage_df3 <- sage_df2 |> 
  anti_join(to_remove)

# Plotting without the removed points
sage_df3 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id))) + 
  geom_point(size = 4) +
  geom_line(size = 0.75) +
  scale_x_continuous(expression(paste("Mass lost (", H[2], "O lost)"))) +
  scale_y_continuous(expression(paste("1/", Psi, " (-", MPa^-1, ")"))) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 19),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 15))

# Concatenate dry weights
dry <- data.frame(id = c("2b", "3a", "4a", "5b", "6", "8"),
                  dry_mass = c(0.8676,
                               1.2583,
                               1.0066,
                               2.4518,
                               2.6039,
                               1.059))

# Combine with dry mass and calculate water content
comb <- sage_df3 |> 
  left_join(dry, by = "id") |> 
  mutate(wc = (total_weight_g - dry_mass) / dry_mass)

comb |> 
  ggplot(aes(x = water_pot_mpa, y = total_weight_g,
             color = factor(id))) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0)

comb |> 
  ggplot(aes(x = wc, y = -1/water_pot_mpa,
             color = factor(id))) +
  geom_point() +
  geom_line()

# Add lm for each bush to obtain intercepts
params_list <- comb |> 
  group_by(id) |> 
  group_map(~ broom::tidy(lm(total_weight_g ~ water_pot_mpa, data = .x)))

params <- do.call(rbind, params_list) |> 
  filter(term == "(Intercept)") |> 
  mutate(id = c("2b", "3a", "4a", "5b", "6", "8"))

# Match to comb
comb2 <- comb |> 
  left_join(select(params, id, estimate), by = "id") |> 
  rename(sat_mass_est = estimate) |>  # estimate saturated mass
  mutate(rwc = (total_weight_g - dry_mass) / (sat_mass_est - dry_mass))

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
          file = "data_clean/sage_pv_curve.csv")
