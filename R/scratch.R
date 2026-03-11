library(tidyverse)
library(patchwork)

load(file = here('data/modsreslagk.RData'))

toplo <- modsreslagk |>
  filter(bay_segment == 'OTB') |>
  filter(tn_load_lag == 'tn_load3') |>
  pull(prds) |>
  pluck(1) |>
  as_tibble() |>
  mutate(
    yr = lubridate::year(date)
  ) |>
  summarise(
    chla = median(chla, na.rm = T),
    btfit = median(btfit, na.rm = T),
    .by = yr
  )

ggplot(toplo, aes(x = yr)) +
  geom_line(aes(y = chla, color = 'Observed')) +
  geom_line(aes(y = btfit, color = 'Simulated')) +
  theme_minimal()
