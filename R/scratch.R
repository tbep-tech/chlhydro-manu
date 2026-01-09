library(tidyverse)
library(mgcv)
library(gratia)

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

yrs <- 2017:2021

trgs <- tbeptools::targets |>
  filter(bay_segment %in% mods$bay_segment) |>
  dplyr::select(bay_segment, thresh = chla_thresh) |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) |>
  tibble()

salyrrates <- wqdat |>
  dplyr::select(date, sal, bay_segment) |>
  mutate(
    yr = year(date)
  ) |>
  summarise(
    sal = mean(sal, na.rm = T),
    .by = c(yr, bay_segment)
  ) |>
  group_nest(bay_segment) |>
  mutate(
    salmod = map(
      data,
      ~ lm(sal ~ yr, data = .x, na.action = 'na.exclude')
    ),
    slopeyr = map_dbl(salmod, ~ coef(.x)[2])
  ) |>
  select(-salmod, -data) |>
  unnest(c('slopeyr')) |>
  crossing(
    yrs = c(1:50)
  ) |>
  mutate(
    futureyr = 2024 + yrs,
    salforeyr = slopeyr * yrs,
    .by = bay_segment
  )

# get salinity avgs by yr, mo for test period
saltstyrmo <- wqdat |>
  dplyr::select(date, sal, bay_segment) |>
  mutate(
    yr = year(date),
    mo = month(date)
  ) |>
  filter(yr %in% yrs) |>
  summarise(
    sal = mean(sal, na.rm = T),
    .by = c(yr, mo, bay_segment)
  )

# add to forecasts
salrates <- salyrrates |>
  left_join(saltstyrmo, by = c('bay_segment'), relationship = 'many-to-many') |>
  mutate(
    salforeyrmo = sal + salforeyr,
    .by = bay_segment
  ) |>
  mutate(
    future_mo = mo,
    date = make_date(year = yr, month = mo, day = 1),
    doy = yday(date),
    dec_time = decimal_date(date)
  )

toplo <- salrates |>
  filter(yrs %in% c(1, 25, 50))
ggplot(toplo, aes(x = date, y = salforeyrmo, color = bay_segment, group = yr)) +
  geom_line() +
  facet_grid(bay_segment ~ yrs, scales = 'free_x') +
  # geom_line(aes(y = chk), linetype = 'dashed') +
  theme_bw()

ldscale <- tibble(
  bay_segment = c('OTB', 'HB', 'MTB', 'LTB')
) |>
  mutate(
    ldscale = map(
      bay_segment,
      ~ ldscale_fun(lddat, ldfac = c(0.5, 1, 2), bay_segment = .x)
    )
  ) |>
  unnest('ldscale')

newdat <- salrates |>
  arrange(bay_segment, yrs, date) |>
  left_join(
    ldscale,
    by = c('bay_segment', 'mo'),
    relationship = 'many-to-many'
  ) |>
  select(-sal) |>
  rename(sal = salforeyrmo)

toplo <- newdat |>
  filter(bay_segment == 'OTB') |>
  filter(yrs %in% c(1, 25, 50))
ggplot(toplo, aes(x = date, y = sal, color = yrs)) +
  geom_line() +
  facet_wrap(ldfac ~ yrs, ncol = 3) +
  theme_bw()
ggplot(toplo, aes(x = date, y = tn_loadlag)) +
  geom_line() +
  facet_grid(ldfac ~ yrs) +
  theme_bw()

tst <- newdat |>
  group_nest(bay_segment) |>
  rename(tst = data) |>
  left_join(trgs, by = c('bay_segment'), relationship = 'one-to-one')

p <- mods |>
  left_join(tst, by = c('bay_segment'), relationship = 'one-to-one') |>
  mutate(
    pyr = pmap(list(mod, tst, thresh), function(mod, tst, thresh) {
      p <- simulate(mod, data = tst, nsim = 500) |>
        bind_cols(tst) |>
        pivot_longer(
          cols = starts_with('sim'),
          names_to = 'sim',
          values_to = 'chla_sim'
        ) |>
        mutate(
          sim = as.integer(str_remove(sim, 'sim_'))
        ) |>
        arrange(sim, date)

      pyr <- p |>
        summarise(
          chla_sim = mean(chla_sim, na.rm = T),
          .by = c(yrs, yr, ldfac, sim)
        ) |>
        summarise(
          chla_sim = mean(chla_sim, na.rm = T),
          .by = c(yrs, ldfac, sim)
        )
      return(pyr)
    }),
    pyr2 = map(pyr, function(pyr) {
      pyr2 <- pyr |>
        summarise(
          percexceeds = mean(chla_sim > thresh) * 100,
          sdexceeds = sd(chla_sim > thresh) * 100,
          .by = c(yrs, ldfac)
        )

      return(pyr2)
    })
  )

toplo <- p |>
  select(bay_segment, pyr) |>
  unnest('pyr') |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  )

ggplot(toplo, aes(x = yrs, y = chla_sim, group = yrs)) +
  geom_boxplot() +
  facet_grid(bay_segment ~ ldfac, scale = 'free_y') +
  theme_bw()

toplo <- p |>
  select(bay_segment, pyr2) |>
  unnest('pyr2') |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  )

ggplot(toplo, aes(x = yrs, y = percexceeds, color = factor(ldfac))) +
  geom_line() +
  coord_cartesian(ylim = c(0, 100)) +
  facet_grid(bay_segment ~ ldfac) +
  geom_ribbon(
    aes(
      ymin = percexceeds - sdexceeds,
      ymax = percexceeds + sdexceeds,
      fill = factor(ldfac)
    ),
    alpha = 0.2
  ) +
  theme_bw()
