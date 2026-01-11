library(tbeptools)
library(tidyverse)
library(mgcv)
library(here)

source(here('R/funcs.R'))

# data prep --------------------------------------------------------------

# water quality
wqdat <- epcdata |>
  filter(yr > 1975 & yr <= 2024) |>
  select(
    bay_segment,
    SampleTime,
    tn,
    chla,
    Sal_Top_ppth,
    Sal_Mid_ppth,
    Sal_Bottom_ppth
  ) |>
  mutate(
    sal = Sal_Top_ppth, #mean(c(Sal_Top_ppth, Sal_Mid_ppth, Sal_Bottom_ppth), na.rm = T),
    .by = c(bay_segment, SampleTime)
  ) |>
  mutate(
    date = floor_date(as.Date(SampleTime), unit = 'month')
  ) |>
  summarise(
    tn = mean(tn, na.rm = T),
    chla = mean(chla, na.rm = T),
    sal = mean(sal, na.rm = T),
    .by = c(bay_segment, date)
  ) |>
  mutate(
    dec_time = decimal_date(date),
    doy = yday(date),
    chla = ifelse(chla == 0, NA, chla),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) |>
  arrange(bay_segment, date) |>
  mutate(
    salorig = sal,
    sal = zoo::na.approx(salorig, x = date, na.rm = F, maxgap = 3),
    .by = c(bay_segment)
  )

# loading data
lddat <- rdataload(
  'https://github.com/tbep-tech/load-estimates/raw/refs/heads/main/data/mosdat.RData'
) |>
  filter(
    !bay_segment %in% c('All Segments (- N. BCB)', 'Remainder Lower Tampa Bay')
  ) |>
  rename(
    yr = year,
    mo = month
  ) |>
  mutate(
    bay_segment = factor(
      bay_segment,
      levels = c(
        'Old Tampa Bay',
        'Hillsborough Bay',
        'Middle Tampa Bay',
        'Lower Tampa Bay'
      ),
      labels = c('OTB', 'HB', 'MTB', 'LTB')
    ),
    date = make_date(year = yr, month = mo, day = 1)
  ) |>
  summarise(
    tn_load = sum(tn_load),
    .by = c(bay_segment, yr, mo, date, source)
  )

save(lddat, file = here('data/lddat.RData'))

load(file = here('data/lddat.RData'))

lddatsum <- lddat |>
  summarise(
    tn_load = sum(tn_load),
    .by = c(bay_segment, date)
  )

wqdat <- wqdat |>
  left_join(lddatsum, by = c('bay_segment', 'date'))

save(wqdat, file = here('data/wqdat.RData'))

# salinity NA by bay segment
wqdat |>
  summarise(
    n_sal_na = sum(is.na(salorig)),
    n_tot = n(),
    .by = c(bay_segment)
  ) |>
  mutate(
    perc_na = n_sal_na / n_tot * 100
  )

# all bay segments -------------------------------------------------------

load(file = here('data/wqdat.RData'))

mods <- wqdat |>
  filter(dec_time >= 1995) |>
  group_nest(bay_segment) |>
  mutate(
    mod = map(
      data,
      ~ gam(
        chla ~ s(dec_time, k = 40, bs = 'tp') +
          s(doy, k = 10, bs = 'cc') +
          s(sal, k = 10, bs = 'tp') +
          s(tn_load, k = 10, bs = 'tp') +
          ti(dec_time, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(dec_time, sal, k = c(5, 5), bs = c('tp', 'tp')) +
          ti(sal, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(dec_time, tn_load, k = c(5, 5), bs = c('tp', 'tp')) +
          ti(tn_load, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(tn_load, sal, k = c(5, 5), bs = c('tp', 'tp')),
        data = .x,
        knots = list(doy = c(0, 366)),
        family = Gamma(link = 'log'),
        na.action = na.exclude,
        method = 'REML'
      )
    ),
    prds = map2(data, mod, pred_fun),
    annsum = map(
      prds,
      ~ data.frame(.x) |>
        mutate(
          yr = lubridate::year(date)
        ) |>
        summarise(
          chla = mean(chla, na.rm = T),
          btfit = mean(btfit, na.rm = T),
          btnorm = mean(btnorm, na.rm = T),
          btfitmd = mean(btfitmd, na.rm = T),
          btnormmd = mean(btnormmd, na.rm = T),
          btfithi = mean(btfithi, na.rm = T),
          btnormhi = mean(btnormhi, na.rm = T),
          btfitlo = mean(btfitlo, na.rm = T),
          btnormlo = mean(btnormlo, na.rm = T),
          .by = c(yr)
        )
    )
  )

save(mods, file = here('data/mods.RData'))

load(file = here('data/mods.RData'))

# run gam.check, edf should be < k', high p-values might suggest auto-correlation
gam.check(mods$mod[[1]])
gam.check(mods$mod[[2]])
gam.check(mods$mod[[3]])
gam.check(mods$mod[[4]])

# check autocor at lag 1, should be < 0.2
acf(residuals(mods$mod[[1]], type = "deviance"), plot = F, lag.max = 1)$acf[2]
acf(residuals(mods$mod[[2]], type = "deviance"), plot = F, lag.max = 1)$acf[2]
acf(residuals(mods$mod[[3]], type = "deviance"), plot = F, lag.max = 1)$acf[2]
acf(residuals(mods$mod[[4]], type = "deviance"), plot = F, lag.max = 1)$acf[2]

summary(mods$mod[[1]])$dev.expl
summary(mods$mod[[2]])$dev.expl
summary(mods$mod[[3]])$dev.expl
summary(mods$mod[[4]])$dev.expl

toplo <- mods |>
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds)

ggplot(toplo, aes(x = log(chla), y = fit)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm', color = 'blue', se = F, formula = y ~ x) +
  facet_wrap(~bay_segment, scales = 'free') +
  theme_minimal()

# use toplo to get r2 between logchla and fit
r2vals <- toplo |>
  mutate(
    logchla = log(chla)
  ) |>
  summarise(
    r2 = cor(logchla, fit, use = 'complete.obs')^2,
    .by = bay_segment
  )

toplo <- mods |>
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |>
  mutate(
    mo = lubridate::month(date)
  ) |>
  filter(dec_time >= 2004 & dec_time <= 2024) |>
  filter(mo %in% 6:11)

ggplot(toplo, aes(x = chla, y = btfit)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  facet_grid(bay_segment ~ mo, scales = 'free') +
  theme_minimal()

# residuals and hydro load -----------------------------------------------

mohydatraw <- rdataload(
  'https://github.com/tbep-tech/load-estimates/raw/refs/heads/main/data/mohydat.RData'
)

hydat <- mohydatraw |>
  filter(
    !bay_segment %in% c('All Segments (- N. BCB)', 'Remainder Lower Tampa Bay')
  ) |>
  rename(
    yr = year,
    mo = month
  ) |>
  mutate(
    bay_segment = factor(
      bay_segment,
      levels = c(
        'Old Tampa Bay',
        'Hillsborough Bay',
        'Middle Tampa Bay',
        'Lower Tampa Bay'
      ),
      labels = c('OTB', 'HB', 'MTB', 'LTB')
    ),
    date = make_date(year = yr, month = mo, day = 1),
    qrt = factor(
      lubridate::quarter(date),
      levels = 1:4,
      labels = c('JFM', 'AMJ', 'JAS', 'OND')
    ),
  )

save(hydat, file = here('data/hydat.RData'))
