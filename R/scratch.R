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
    slopeyr = map_dbl(salmod, ~ coef(.x)[2]),
    meanporsal = map_dbl(
      data,
      ~ mean(.x |> filter(yr >= 1995) |> pull(sal), na.rm = T)
    ),
    lastyr = map(
      salmod,
      ~ {
        lastyr <- max(.x$model$yr)
        lastsal <- predict(.x, newdata = data.frame(yr = lastyr))
        data.frame(yr = lastyr, sal = lastsal)
      }
    )
  ) |>
  select(-salmod, -data) |>
  unnest(c('meanporsal', 'lastyr', 'slopeyr')) |>
  crossing(
    yrs = c(1:50)
  ) |>
  mutate(
    futureyr = 2024 + yrs,
    salforeyr = sal + (slopeyr * yrs),
    .by = bay_segment
  )

# get salinity monthly differences
modiffs <- wqdat |>
  dplyr::select(date, sal, bay_segment) |>
  mutate(
    yr = year(date),
    mo = month(date)
  ) |>
  filter(yr >= 1995) |>
  summarise(
    sal = mean(sal, na.rm = T),
    .by = c(mo, bay_segment)
  ) |>
  mutate(
    modiff = sal - mean(sal, na.rm = T),
    .by = bay_segment
  ) |>
  select(-sal)

# add to forecasts
salrates <- salyrrates |>
  left_join(modiffs, by = c('bay_segment'), relationship = 'many-to-many') |>
  mutate(
    salforeyrmo = salforeyr + modiff,
    .by = bay_segment
  ) |>
  mutate(
    date = make_date(year = futureyr, month = mo, day = 1),
    doy = yday(date),
    dec_time = decimal_date(date)
  )

ggplot(salrates, aes(x = date, y = salforeyr, color = bay_segment)) +
  geom_line() +
  geom_line(aes(y = salforeyrmo)) +
  facet_wrap(~bay_segment) +
  # geom_line(aes(y = chk), linetype = 'dashed') +
  theme_bw()

timevec <- crossing(
  yr = c(2010:2024),
  mo = 1:12
) |>
  mutate(
    date = make_date(year = yr, month = mo, day = 1),
    doy = yday(date),
    dec_time = decimal_date(date)
  )

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


modin <- mods$mod[[4]]

newdat <- salrates |>
  select(bay_segment, date, yrs, sal = salforeyrmo, doy, dec_time) |>
  group_nest(yrs, bay_segment) |>
  mutate(
    data = map(
      data,
      ~ {
        left_join(
          .x |> mutate(mo = month(date)) |> select(mo, sal),
          timevec,
          by = c('mo')
        )
      }
    )
  ) |>
  unnest('data') |>
  arrange(bay_segment, yrs, date) |>
  left_join(ldscale, by = c('bay_segment', 'mo'), relationship = 'many-to-many')

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
  filter(bay_segment == 'OTB')

p <- simulate(mods$mod[[1]], data = tst, nsim = 200) |>
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

ggplot(pyr, aes(x = yrs, y = chla_sim, group = yrs)) +
  geom_boxplot() +
  facet_wrap(~ldfac) +
  theme_bw()

pyr2 <- pyr |>
  summarise(
    percexceeds = mean(chla_sim > 9.8) * 100,
    sdexceeds = sd(chla_sim > 9.8) * 100,
    .by = c(yrs, ldfac)
  )
ggplot(pyr2, aes(x = yrs, y = percexceeds, color = factor(ldfac))) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = percexceeds - sdexceeds,
      ymax = percexceeds + sdexceeds,
      fill = factor(ldfac)
    ),
    alpha = 0.2
  ) +
  theme_bw()
