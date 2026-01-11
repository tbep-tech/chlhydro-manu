# fit salinity trends out to 2050, created 1000 predictions using prediction intervals
# use four loading scenarios: observed, constant at TMDL, TMDL 2x increase, TMDL 50% decrease
# simulations are selected for conditions for a given year period
# where salinity and loads are changed proportionally given observed monthly/annual conditions for the year period
# then plot likelihood of exceeding thresholds by time for each scenario

# must be sequential
yrs <- c(2021:2024)

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

# toplo <- salrates |>
#   filter(yrs %in% c(1, 25, 50))
# ggplot(toplo, aes(x = date, y = salforeyrmo, color = bay_segment, group = yr)) +
#   geom_line() +
#   facet_grid(bay_segment ~ yrs, scales = 'free_x') +
#   # geom_line(aes(y = chk), linetype = 'dashed') +
#   theme_bw()

ldscale <- tibble(
  bay_segment = c('OTB', 'HB', 'MTB', 'LTB')
) |>
  mutate(
    ldscale = map(
      bay_segment,
      ~ ldscale_fun(lddat, ldfac = c(0.5, 1, 2), yrs = yrs, bay_segment = .x)
    )
  ) |>
  unnest('ldscale')

# toplo <- ldscale |>
#   mutate(
#     tn_loadannchk = sum(tn_load),
#     .by = c(bay_segment, yr, ldfac)
#   )
# ggplot (toplo, aes(x = date, y = tn_load, color = factor(ldfac))) +
#   geom_line() +
#   facet_wrap(~bay_segment, scales = 'free_y') +
#   theme_minimal()

# ggplot (toplo, aes(x = date, y = tn_loadann, color = factor(ldfac))) +
#   geom_line(aes(linetype = 'Actual tn_load')) +
#   geom_line(aes(y = tn_loadannchk, linetype = 'Est sum tn_load')) +
#   facet_wrap(~bay_segment, scales = 'free_y') +
#   theme_minimal()

newdat <- salrates |>
  arrange(bay_segment, yrs, date) |>
  left_join(
    ldscale,
    by = c('bay_segment', 'date', 'yr', 'mo'),
    relationship = 'many-to-many'
  ) |>
  select(-sal) |>
  rename(sal = salforeyrmo)

# toplo <- newdat |>
#   filter(bay_segment == c('OTB')) |>
#   filter(yrs %in% c(1, 25, 50))
# ggplot(toplo, aes(x = date, y = sal, color = yrs)) +
#   geom_line() +
#   facet_wrap(ldfac ~ yrs, ncol = 3) +
#   theme_bw()
# ggplot(toplo, aes(x = date, y = tn_load)) +
#   geom_line() +
#   facet_grid(ldfac ~ yrs) +
#   theme_bw()

tst <- newdat |>
  group_nest(bay_segment) |>
  rename(tst = data) |>
  left_join(trgs, by = c('bay_segment'), relationship = 'one-to-one')

p <- mods |>
  left_join(tst, by = c('bay_segment'), relationship = 'one-to-one') |>
  mutate(
    p = pmap(list(mod, tst), function(mod, tst) {
      simulate(mod, data = tst, nsim = 500) |>
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
    }),
    pyr = map(p, function(p) {
      p |>
        summarise(
          chla_sim = mean(chla_sim, na.rm = T),
          .by = c(yrs, yr, ldfac, sim)
        )
    }),
    pyr2 = pmap(list(pyr, thresh), function(pyr, thresh) {
      pyr |>
        mutate(
          exceeds = chla_sim > thresh
        ) |>
        summarise(
          percexceeds = mean(exceeds, na.rm = T) * 100,
          .by = c(yrs, sim, ldfac)
        )
    }),
    pyr3 = map(pyr2, function(pyr2) {
      pyr2 |>
        summarise(
          avexceeds = mean(percexceeds, na.rm = T),
          sdexceeds = sd(percexceeds, na.rm = T),
          .by = c(yrs, ldfac)
        )
    })
  )

toplo <- p |>
  select(bay_segment, p) |>
  unnest('p') |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) |>
  filter(yrs %in% c(1, 25, 50)) |>
  filter(bay_segment == 'OTB')
chlact <- wqdat |>
  filter(bay_segment == 'OTB') |>
  select(date, chla) |>
  filter(date >= min(toplo$date) & date <= max(toplo$date))
ggplot(toplo, aes(x = date, y = chla_sim, group = sim)) +
  geom_line(alpha = 0.2) +
  geom_point(
    data = chlact,
    aes(x = date, y = chla),
    color = 'red',
    size = 1,
    inherit.aes = F
  ) +
  facet_grid(ldfac ~ yrs) +
  theme_bw()

toplo <- p |>
  select(bay_segment, pyr) |>
  unnest('pyr') |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) |>
  filter(yrs %in% c(1, 25, 50)) |>
  filter(bay_segment == 'OTB')
chlact <- wqdat |>
  filter(bay_segment == 'OTB') |>
  select(date, chla) |>
  mutate(
    yr = year(date)
  ) |>
  filter(yr >= min(toplo$yr) & yr <= max(toplo$yr)) |>
  summarise(
    chla = mean(chla, na.rm = T),
    .by = c(yr)
  )

ggplot(toplo, aes(x = yr, y = chla_sim, group = sim)) +
  geom_line(alpha = 0.2) +
  geom_point(
    data = chlact,
    aes(x = yr, y = chla),
    color = 'red',
    size = 1,
    inherit.aes = F
  ) +
  facet_grid(ldfac ~ yrs) +
  theme_bw()

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

ggplot(
  toplo,
  aes(x = yrs, y = percexceeds, color = factor(ldfac), group = yrs)
) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 100)) +
  facet_grid(bay_segment ~ ldfac) +
  theme_bw()

toplo <- p |>
  select(bay_segment, pyr3) |>
  unnest('pyr3') |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  )

ggplot(toplo, aes(x = yrs, y = avexceeds, color = factor(ldfac))) +
  geom_line() +
  coord_cartesian(ylim = c(0, 100)) +
  facet_grid(bay_segment ~ ldfac) +
  geom_ribbon(
    aes(
      ymin = avexceeds - sdexceeds,
      ymax = avexceeds + sdexceeds,
      fill = factor(ldfac)
    ),
    alpha = 0.2
  ) +
  theme_bw()


# check predictions for ld / 2
mod <- mods$mod[[1]]

dat <- mods$data[[1]]
toplo <- dat |>
  mutate(
    prds = predict(mod, type = 'response'),
    prds2 = predict(
      mod,
      newdata = dat |> mutate(tn_load = tn_load / 2),
      type = 'response'
    )
  )

ggplot(toplo, aes(x = date, y = chla)) +
  geom_point() +
  geom_line(aes(y = prds), color = 'blue') +
  geom_line(aes(y = prds2), color = 'red') +
  theme_minimal()
