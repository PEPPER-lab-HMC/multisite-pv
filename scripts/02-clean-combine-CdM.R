#### Wrangling and plotting PV curve data for juniper ####

library(googlesheets4)
library(googledrive)
library(tidyverse)
library(janitor)
library(broom)
library(here)
theme_set(theme_bw())

#### Reading in juniper data ####
ids <- paste0(1:7)
juniper_list <- list()
for(i in 1:length(ids)) {
  juniper_import <- read_sheet("https://docs.google.com/spreadsheets/d/1ycVFmL2Df4YFGp70lT7WdsF4Xn1dGllPtAnFksf5IHU/edit?usp=sharing",
                               sheet = paste0("Juniper ", ids[i])) |>
    clean_names() |> 
    mutate(id = i,
           time_c = as.character(time) |> 
             str_extract(pattern = "\\d\\d:\\d\\d"),
           hour = hm(time_c),
           date_c = if_else(hour < hour[1], "03/13/2025", "03/12/2025"),
           dt = as.POSIXct(paste(date_c, time_c), format = "%m/%d/%Y %H:%M", tz = "America/Denver")) |> 
    relocate(dt, id) |> 
    select(-time, -time_c, -hour, -date_c)
  
  juniper_list[[i]] <- juniper_import
}
juniper_df <- do.call(bind_rows, juniper_list) |> 
  select(-fallen_leaf_weight_g, -fallen_bark_dryweight_g, -fallen_leaf_dryweight_g, -fallen_leaf_number_2_weight_g, -fallen_leaf_2_dryweight_g, -dryweight_g, -fallen_leaf_2_weight_g, -maybe_fallen_leaf_dryweight_g) |> 
  group_by(id) |> 
  summarize(
    dt = dt,
    id = id,
    water_pot_mpa = water_pot_mpa,
    total_weight_g = total_weight_g,
    dry_weight_g = max(dry_weight_g, na.rm = TRUE),
    start_wp = max(water_pot_mpa, na.rm = TRUE)
  ) |> 
  ungroup()

# Create mass_lost column
lost <- data.frame(id = 1:7,
                   start_weight = c(1.9925,
                                    3.7049,
                                    8.2263,
                                    2.8165,
                                    3.3955,
                                    8.4901,
                                    1.7764))
juniper_df2 <- juniper_df |> 
  left_join(lost, by = "id") |> 
  mutate(mass_lost_g = abs(total_weight_g - start_weight))

#### Plotting PV Curves ####
# Plotting inverse water pot over mass lost
juniper_df2 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id))) + 
  geom_point() +
  geom_line()

# Throw out Juniper 2 for unreliable data
juniper_df2 <- juniper_df2 |> 
  filter(id != "2")

# Clean data for WP measurements taken too close together
to_remove <- juniper_df2 |> 
  group_by("id") |> 
  mutate(lag_wp  = lag(water_pot_mpa),
         diff_wp = lag_wp - water_pot_mpa) |> 
  ungroup() |> 
  filter(diff_wp <= 0.06) 

juniper_df3 <- juniper_df2 |> 
  anti_join(to_remove) |> 
  mutate(superID = consecutive_id(id))

# Plotting without the removed points
juniper_df3 |> 
  ggplot(aes(x = mass_lost_g, y = -1/water_pot_mpa, color = as.factor(id))) + 
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

write_csv(juniper_df3,
          "data_clean/02-CdM-cleaning-output.csv")
