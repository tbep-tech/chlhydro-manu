library(tidyverse)
library(mgcv)

load(file = 'data/mods.RData')
load(file = 'data/wqdat.RData')
load(file = 'data/lddat.RData')

source('R/funcs.R')

bayseg <- 'OTB'
prdgrd <- wqdat |>
  dplyr::select(date, sal, bay_segment) |>
  mutate(
    yr = year(date)
  ) |>
  filter(yr >= 1995) |>
  summarise(
    sal = mean(sal, na.rm = T),
    .by = c(yr, bay_segment)
  ) |>
  mutate(
    anngrp = case_when(
      yr < 2010 ~ '1975-2009',
      yr >= 2010 ~ '2010-2024'
    )
  ) |>
  mutate(
    md = quantile(sal, 0.5, na.rm = T),
    hi = quantile(sal, 0.9, na.rm = T),
    lo = quantile(sal, 0.1, na.rm = T),
    .by = c(bay_segment, anngrp)
  )

ggplot(prdgrd, aes(x = yr, y = sal, group = anngrp)) +
  geom_point() +
  geom_line(aes(y = md), color = 'blue', size = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  facet_wrap(~bay_segment, scales = 'free_y') +
  theme_bw()

# fit salinity trends out to 2050, created 1000 predictions using prediction intervals
# use three loading scenarios: constant at TMDL and 2x increase, 50% decrease
# how do I do this on a monthly basis - remember models use three month cumulative lag load
# maybe split into approximate percentages based on historcal monthly distribution? then apply lags

# then plot likelihood of exceeding thresholds by time for each scenario

trgs <- tbeptools::targets |>
  filter(bay_segment %in% mods$bay_segment) |>
  dplyr::select(bay_segment, thresh = chla_thresh) |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) |>
  tibble()


salmods <- wqdat |>
  dplyr::select(date, sal, bay_segment) |>
  mutate(
    yr = year(date),
    mo = month(date)
  ) |>
  summarise(
    sal = mean(sal, na.rm = T),
    .by = c(yr, mo, date, bay_segment)
  ) |>
  group_nest(bay_segment) |>
  mutate(
    salmod = map(
      data,
      ~ lm(sal ~ yr + factor(mo), data = .x, na.action = 'na.exclude')
    )
  )

prds <- salmods |>
  select(-data) |>
  left_join(mods, by = 'bay_segment') |>
  select(-data, -prds, -annsum) |>
  rename(
    gmmod = mod
  ) |>
  mutate(
    ldscale = map(
      bay_segment,
      ~ ldscale_fun(lddat, ldfac = c(0.5, 1, 2), bay_segment = .x)
    ),
    salfore = map(salmod, ~ salfore_fun(.x, nsim = 100, yrs = c(2025, 2050))),
    gamfore = pmap(list(salfore, ldscale, gmmod), ~ gamfore_fun(..1, ..2, ..3))
  )

toplo <- prds |>
  select(bay_segment, salfore) |>
  unnest('salfore') |>
  summarise(
    sal = mean(sal),
    .by = c(bay_segment, yr, sim)
  )

ggplot(toplo, aes(x = yr, y = sal, group = sim)) +
  geom_point(alpha = 0.3) +
  stat_smooth(
    geom = 'line',
    method = 'lm',
    color = 'blue',
    alpha = 0.3,
    formula = y ~ x
  ) +
  facet_wrap(~bay_segment, ncol = 2) +
  theme_bw()

toplo <- prds |>
  unnest('gamfore') |>
  summarise(
    chla = mean(chla),
    .by = c(bay_segment, yr, sim, ldfac)
  )

ggplot(toplo, aes(x = yr, y = chla, group = sim)) +
  geom_line(alpha = 0.3) +
  facet_grid(bay_segment ~ ldfac, scales = 'free_y') +
  theme_bw()


toplo2 <- toplo |>
  left_join(trgs, by = 'bay_segment') |>
  summarise(
    percexceeds = mean(chla > thresh) * 100,
    sdexceeds = sd(chla > thresh) * 100,
    .by = c(bay_segment, ldfac, yr)
  )

ggplot(toplo2, aes(x = yr, y = percexceeds, color = factor(ldfac))) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = percexceeds - sdexceeds,
      ymax = percexceeds + sdexceeds,
      fill = factor(ldfac)
    ),
    alpha = 0.2
  ) +
  facet_grid(bay_segment ~ ldfac) +
  theme_bw()

tmp <- prds |>
  filter(bay_segment == 'LTB') |>
  select(bay_segment, gamfore) |>
  unnest('gamfore') |>
  filter(ldfac == 'TMDL Load Factor: 1')

toplo <- tmp |>
  filter(sim == 1) |>
  filter(yr <= 2030)

ggplot(toplo, aes(x = date, y = chla)) +
  geom_line() +
  theme_bw()

tmp2 <- simulate(mods$mod[[4]], nsim = 100) |>
  bind_cols(mods$data[[4]]) |>
  pivot_longer(
    cols = starts_with('sim'),
    names_to = 'sim',
    values_to = 'chla_sim'
  ) |>
  mutate(
    sim = as.integer(str_remove(sim, 'sim_'))
  ) |>
  arrange(sim, date)

ggplot(tmp2, aes(x = date, y = chla_sim, group = sim)) +
  geom_line(alpha = 0.3) +
  geom_point(aes(y = chla), color = 'red') +
  theme_bw()

tmp3 <- tmp2 |>
  mutate(
    yr = year(date)
  ) |>
  summarise(
    chla_mean = mean(chla_sim, na.rm = T),
    chla = mean(chla, na.rm = T),
    .by = c(yr, sim)
  )

ggplot(tmp3, aes(x = yr, y = chla_mean, group = sim)) +
  geom_line(alpha = 0.3) +
  geom_point(aes(y = chla), color = 'red') +
  theme_bw()
