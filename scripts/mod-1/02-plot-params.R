#
library(tidyverse)
library(broom.mixed)
library(coda) # needed to recognize mcmc.list type

# Load both codas
load(file = "scripts/mod-1/coda/jm_coda.Rdata")
load(file = "scripts/mod-1/coda/jm_rep.Rdata")

# Load in original data
pv <- read_csv("data_clean/pv_comb_20240308.csv")

# Function to predict total water potential
tpf <- function(x,vec) {
  # vec = vector of a, b, c, d
  # a*exp(-x*b)+cc+d*x
  
  y = vec[1]*exp(-x*vec[2]) + vec[3] + vec[4]*x
  return(y)
}

# Function to predict mass_loss of TLP
xTLP <- function(vec) {
  # From minimizing f(x) = a*exp(-b*x) - d*x
  # by setting f'(x) = 0
  xtlp <- -log(-vec[4]/(vec[1]*vec[2]))/vec[2]
  return(xtlp)
}

# summarize parameters
coda_sum <- tidyMCMC(jm_coda,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

sample_params <- coda_sum |> 
  filter(grepl("^[abcd]\\[\\d", term)) |> 
  mutate(param = str_extract(term, "^[abcd]"),
         sample = str_extract(term, "\\d")) |> 
  arrange(sample)

pop_params <- coda_sum |> 
  filter(grepl("^mu\\.[abcd]", term)) |> 
  mutate(param = str_extract(term, "[abcd]"))

# Create dataframe for predictions and TLP for each sample
pred_list <- list()
tlp_df <- data.frame(ID = numeric(0),
                     xTLP = numeric(0),
                     inv_TLP = numeric(0))
for(i in 1:5) {
  mass_lost <- seq(0, max(pv$mass_lost[pv$ID == i]), 0.01)
  pred_invP <- tpf(mass_lost,
                   sample_params |> 
                     filter(sample == i) |> 
                     select(pred.mean) |> 
                     pull())
  
  x <- xTLP(sample_params |> 
              filter(sample == i) |> 
              select(pred.mean) |> 
              pull())
  tlp_df[i,] <- c(i, x, tpf(x, sample_params |> 
                              filter(sample == i) |> 
                              select(pred.mean) |> 
                              pull()))
  
  pred_list[[i]] <- data.frame(ID = i, pred_invP, mass_lost)
}
pred_df <- do.call(rbind, pred_list)
tlp_df$TLP <- -1/tlp_df$inv_TLP

# Calculate population-level curve and TLP
mass_lost <- seq(0, max(pv$mass_lost), 0.01)
pred_invP <- tpf(mass_lost,
                 pop_params |> 
                   select(pred.mean) |> 
                   pull())

pop_pred_df <- data.frame(mass_lost, pred_invP)

x <- xTLP(pop_params |> 
            select(pred.mean) |> 
            pull())

inv_TLP <- tpf(x,
               pop_params |> 
                 select(pred.mean) |> 
                 pull())
pop_tlp_df <- data.frame(xTLP = x, 
                         inv_TLP = inv_TLP,
                         TLP = -1/inv_TLP)


# Plot of data + posterior mean predictions
ggplot() +
  geom_line(data = pop_pred_df, 
            aes(x = mass_lost, y = pred_invP),
            color = "gray40",
            show.legend = FALSE) +
  geom_vline(data = pop_tlp_df, 
             aes(xintercept = xTLP),
             lty = 3,
             color = "gray40",
             show.legend = FALSE) +
  geom_hline(data = pop_tlp_df, 
             aes(yintercept = inv_TLP),
             lty = 3,
             color = "gray40",
             show.legend = FALSE) +
  geom_text(dat = pop_tlp_df, 
            aes(x= 1.5, y = 2,
                label = round(TLP, 3)),
            color = "gray40",
            hjust = 1,
            show.legend = FALSE) +
  geom_line(data = pred_df, 
            aes(x = mass_lost, y = pred_invP, 
                color = factor(ID)),
            show.legend = FALSE) +
  geom_point(data = pv, 
             aes(x = mass_lost, y = 1/P.MPa,
                 color = factor(ID)),
             size = 2) +
  geom_vline(data = tlp_df, aes(xintercept = xTLP,
                                color = factor(ID)),
             lty = 2,
             show.legend = FALSE) +
  geom_hline(data = tlp_df, aes(yintercept = inv_TLP,
                                color = factor(ID)),
             lty = 2,
             show.legend = FALSE) +
  geom_text(dat = tlp_df, aes(x= 1.5, y = 0.5,
                              label = round(TLP, 3),
                              color = factor(ID)),
            hjust = 1,
            show.legend = FALSE) +
  facet_wrap(~ID) +
  scale_x_continuous("Mass lost (g)") +
  scale_y_continuous(expression(paste("1 / ", Psi, " (-MPa)"))) +
  theme_bw(base_size = 12) +
  labs(color = "Sample") +
  theme(legend.position = c(0.8, 0.2),
        panel.grid = element_blank()) +
  guides(color = "none")

            