# loop through and combine raw data
library(tidyverse)

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
  geom_line()

# Save out only TRUE points
fn_df |> 
  filter(keep == TRUE) |> 
  write_csv("data_clean/pv_comb.csv")
