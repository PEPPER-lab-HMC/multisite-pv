# Test run of Duncan Smith's pv function

source("source/pv - function2.R")
source("source/corner - function.R")
source("source/c.loc - function.R")

#### 20231027 run ####
# collected 20231026, bagged at predawn, wrapped cut end in damp paper towel
# left overnight in dark cooler

pv("latr_1")
pv("latr_2")
pv("latr_3")
pv("latr_4")
pv("latr_5")

pv("latr_1", plotonly = T)
pv("latr_2", plotonly = T)
pv("latr_3", plotonly = T)
pv("latr_4", plotonly = T)
pv("latr_5", plotonly = T)

# When completed, csv's moved to data_raw/

#### 20240308 run ####
# collected 20230308, bagged at predawn in ziptop bag
# was actively raining during predawn, very high values

pv("latr_1")
pv("latr_2")
pv("latr_3")
pv("latr_4")
pv("latr_5")
pv("latr_6")

pv("latr_1", plotonly = T)
pv("latr_2", plotonly = T)
pv("latr_3", plotonly = T)
pv("latr_4", plotonly = T)
pv("latr_5", plotonly = T)
pv("latr_6", plotonly = T)
