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

mosdat <- rdataload(
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
    .by = c(bay_segment, date)
  ) |>
  mutate(
    lag1 = dplyr::lag(tn_load, n = 1),
    lag2 = dplyr::lag(tn_load, n = 2),
    lag3 = dplyr::lag(tn_load, n = 3),
    .by = c(bay_segment)
  ) |>
  fill(
    lag1,
    lag2,
    lag3,
    .by = c(bay_segment),
    .direction = 'up'
  ) |>
  mutate(
    tn_load = tn_load + lag1 + lag2 + lag3
  ) |>
  select(-lag1, -lag2, -lag3)

wqdat <- wqdat |>
  left_join(mosdat, by = c('bay_segment', 'date'))

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

# OTB model --------------------------------------------------------------

tomod <- wqdat |>
  filter(bay_segment == 'OTB') |>
  filter(!is.na(tn_load))

mod <- gam(
  chla ~ s(dec_time, k = 40, bs = 'tp') +
    s(doy, k = 10, bs = 'cc') +
    s(sal, k = 10) +
    s(tn_load, k = 10) +
    ti(dec_time, doy, k = c(5, 5), bs = c('tp', 'cc')) +
    ti(dec_time, sal, k = c(5, 5)) +
    ti(sal, doy, k = c(5, 5), bs = c('tp', 'cc')) +
    ti(dec_time, tn_load, k = c(5, 5), bs = c('tp', 'tp')) +
    ti(tn_load, doy, k = c(5, 5), bs = c('tp', 'cc')) +
    ti(tn_load, sal, k = c(5, 5), bs = c('tp', 'tp')),
  data = tomod,
  family = Gamma(link = 'log'),
  knots = list(doy = c(0, 366)),
  method = 'REML'
)

toprd1 <- crossing(
  tn_load = c(max(tomod$tn_load, na.rm = T), min(tomod$tn_load, na.rm = T)),
  sal = c(max(tomod$sal, na.rm = T), min(tomod$sal, na.rm = T)),
  tomod |> select(-tn_load, -sal)
)

prd1 <- toprd1 |>
  mutate(
    prd = predict(mod, newdata = toprd1, type = 'response')
  ) |>
  mutate(
    loadcond = case_when(
      tn_load == max(tomod$tn_load, na.rm = T) ~ 'highload',
      tn_load == min(tomod$tn_load, na.rm = T) ~ 'loload'
    ),
    salcond = case_when(
      sal == max(tomod$sal, na.rm = T) ~ 'highsal',
      sal == min(tomod$sal, na.rm = T) ~ 'losal'
    ),
    yr = year(date)
  ) |>
  summarise(
    prd = mean(prd, na.rm = T),
    .by = c(yr, loadcond, salcond)
  )

ggplot(prd1, aes(x = yr, y = prd, color = loadcond)) +
  geom_line() +
  facet_wrap(~salcond) +
  theme_minimal()


summary(mod)$dev.expl
gam.check(mod)
acf(residuals(mod, type = "deviance")) # check for autocorrelation in residds
prds <- pred_fun(tomod, mod)

ggplot(prds, aes(x = dec_time, y = chla)) +
  geom_point() +
  geom_line(aes(y = btfit), color = 'red') +
  geom_line(aes(y = btfithi), color = 'blue')


grid_plo(prds, month = 'all', allsal = T)
grid_plo(prds, month = 'all', allsal = T, years = c(2000, 2024))
grid_plo(prds)
grid_plo(
  prds,
  month = c(5:10),
  years = c(2000, 2024),
  allsal = F,
  ncol = 6,
  sal_fac = 6
)

toplo <- data.frame(prds) |>
  mutate(
    yr = lubridate::year(date)
  ) |>
  select(
    date,
    yr,
    btfit,
    btfithi,
    btfitlo,
    btfitmd,
    btnorm,
    btnormhi,
    btnormlo,
    btnormmd
  ) |>
  pivot_longer(
    cols = -c(date, yr),
    names_to = 'var',
    values_to = 'value'
  ) |>
  mutate(
    nrmval = ifelse(grepl('norm', var, ignore.case = T), 'norm', 'fit'),
    est = case_when(
      grepl('hi', var, ignore.case = T) ~ 'hi',
      grepl('lo', var, ignore.case = T) ~ 'lo',
      grepl('md', var, ignore.case = T) ~ 'md',
      TRUE ~ 'obs'
    )
  ) |>
  summarise(
    val = mean(value, na.rm = T),
    .by = c(yr, nrmval, est)
  ) |>
  pivot_wider(
    names_from = c(nrmval),
    values_from = val
  )

ggplot(toplo, aes(x = yr, y = fit)) +
  geom_point() +
  geom_line(aes(y = norm)) +
  facet_wrap(~est, scales = 'free') +
  theme_minimal()

toplo <- data.frame(prds)

ggplot(toplo, aes(x = yr, y = fit, group = est, color = est)) +
  geom_point() +
  geom_line(aes(y = norm)) +
  facet_wrap(~nrmval, scales = 'free') +
  theme_minimal()

# all bay segments -------------------------------------------------------

mods <- wqdat |>
  # filter(dec_time >= 2000) |>
  group_nest(bay_segment) |>
  mutate(
    mod = map(
      data,
      ~ gam(
        chla ~ s(dec_time, k = 40, bs = 'tp') +
          s(doy, k = 10, bs = 'cc') +
          s(sal, k = 10) +
          ti(dec_time, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(dec_time, sal, k = c(5, 5)) +
          ti(sal, doy, k = c(5, 5), bs = c('tp', 'cc')),
        data = .x,
        knots = list(doy = c(0, 366)),
        family = Gamma(link = 'log'),
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
          meansalfit = mean(meansalfit, na.rm = T),
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
