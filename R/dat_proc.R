library(tbeptools)
library(tidyverse)
library(mgcv)
library(here)
library(patchwork)

source(here('R/funcs.R'))

# data prep --------------------------------------------------------------

# water quality
wqdat <- epcdata |> 
  filter(yr > 1975 & yr <= 2024) |> 
  select(bay_segment, SampleTime, tn, chla, Sal_Top_ppth, Sal_Mid_ppth, Sal_Bottom_ppth) |>
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
  filter(bay_segment == 'HB')

mod <- gam(chla ~ s(dec_time, k = 40, bs = 'tp') + 
  s(doy, k = 10, bs = 'cc') + 
  s(sal, k = 10) +
  ti(dec_time, doy, k = c(5, 5), bs = c('tp', 'cc')) +
  ti(dec_time, sal, k = c(5, 5)) +
  ti(sal, doy, k = c(5, 5), bs = c('tp', 'cc')),
  data = tomod,
  family = Gamma(link = 'log'),
  knots = list(doy=c(0,366)),
  method = 'REML'
)

summary(mod)$dev.expl
gam.check(mod)
acf(residuals(mod, type = "deviance")) # check for autocorrelation in residds
prds <- pred_fun(tomod, mod)

ggplot(prds, aes(x = dec_time, y = chla)) +
  geom_point() +
  geom_line(aes(y = btfit), color = 'red')


grid_plo(prds, month = 'all', allsal = T)
grid_plo(prds, month = 'all', allsal = T, years = c(2000, 2024))
grid_plo(prds)
grid_plo(prds, month = c(5:10), years = c(2000, 2024), allsal = F, ncol = 6, sal_fac = 6)

toplo <- data.frame(prds) |> 
  mutate(
    yr = lubridate::year(date)
  ) |> 
  summarise(
    btfit = mean(btfit, na.rm = T), 
    btnorm = mean(btnorm, na.rm = T),
    .by = c(yr)
  )

ggplot(toplo, aes(x = yr, y = btfit)) +
  geom_point() + 
  geom_line(aes(y = btnorm), color = 'red') + 
  theme_minimal()

toplo <- data.frame(prds)

ggplot(toplo, aes(x = date, y = chla)) +
  geom_point() +
  geom_line(aes(y = btfit), color = 'red') +
  geom_line(aes(y = btnorm), color = 'blue') +
  theme_minimal()

# all bay segments -------------------------------------------------------

mods <- wqdat |> 
  # filter(dec_time >= 2000) |>
  group_nest(bay_segment) |> 
  mutate(
    mod = map(data, ~ gam(chla ~ s(dec_time, k = 40, bs = 'tp') + 
          s(doy, k = 10, bs = 'cc') + 
          s(sal, k = 10) +
          ti(dec_time, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(dec_time, sal, k = c(5, 5)) +
          ti(sal, doy, k = c(5, 5), bs = c('tp', 'cc')),
        data = .x,
        knots = list(doy=c(0,366)),
        family = Gamma(link = 'log'),
        method = 'REML'
      )
    ),
    prds = map2(data, mod, pred_fun), 
    annsum = map(prds, ~ data.frame(.x) |>
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
  facet_wrap(~ bay_segment, scales = 'free') +
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

mohydatraw <- rdataload('https://github.com/tbep-tech/load-estimates/raw/refs/heads/main/data/mohydat.RData')

hydat <- mohydatraw |> 
  filter(!bay_segment %in% c('All Segments (- N. BCB)', 'Remainder Lower Tampa Bay')) |> 
  rename(
    yr = year,
    mo = month
  ) |> 
  mutate(
    bay_segment = factor(bay_segment,
      levels = c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay'), 
     labels = c('OTB', 'HB', 'MTB', 'LTB')), 
    date = make_date(year = yr, month = mo, day = 1),
    qrt = factor(lubridate::quarter(date), levels = 1:4, labels = c('JFM', 'AMJ', 'JAS', 'OND')),
  )

save(hydat, file = here('data/hydat.RData'))


load(file = here('data/mods.RData'))
load(file = here('data/hydat.RData'))

qrthydat <- hydat |> 
  mutate(
    date = floor_date(date, unit = 'quarter')
  ) |> 
  summarise(
    hyqrt = sum(hy_load_106_m3_mo, na.rm = T),
    .by = c(bay_segment, date, yr, qrt)
  )

qrtprds <- mods |> 
  select(bay_segment, prds) |> 
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |> 
  mutate(
    yr = lubridate::year(date),
    qrt = factor(lubridate::quarter(date), levels = 1:4, labels = c('JFM', 'AMJ', 'JAS', 'OND')),
    date = floor_date(date, unit = 'quarter'),
  ) |> 
  summarise(
    chla = mean(chla, na.rm = T),
    btfit = mean(btfit, na.rm = T), 
    btnorm = mean(btnorm, na.rm = T),
    .by = c(bay_segment, date, yr, qrt)
  ) |> 
  inner_join(qrthydat, by = c('bay_segment', 'yr', 'qrt', 'date')) |> 
  mutate(
    prdresid = btfit - btnorm
  )

toplo <- qrtprds |> 
  filter(yr >= 2004 & yr <= 2024) |> 
  filter(bay_segment == 'OTB') |> 
  mutate(
    hyresid = hyqrt - mean(hyqrt, na.rm = T), 
    .by = qrt
  )

# convert hyresid and prdresid to z-score
toplo <- toplo |> 
  mutate(
    # hyresid = (hyresid - mean(hyresid, na.rm = T)) / sd(hyresid, na.rm = T), 
    # prdresid = (prdresid - mean(prdresid, na.rm = T)) / sd(prdresid, na.rm = T),
    rsq = round(summary(lm(prdresid ~ hyresid))$r.squared, 2), 
    rsq = paste(qrt, ' (R² = ', rsq, ')', sep = ''), 
    .by = qrt
  )

fac <- toplo |> 
  select(qrt, rsq) |>
  distinct()

toplo <- toplo |>
  mutate(
    qrtlab = factor(qrt, levels = fac$qrt, labels = fac$rsq)
  )

p1 <- ggplot(toplo, aes(x = yr)) + 
  geom_point(aes(y = btfit)) + 
  geom_line(aes(y = btnorm), color = 'red') +
  facet_wrap(~ qrt, ncol = 4) + 
  theme_minimal() +
  labs(
    y = 'Chl-a (µg/L)',
    x = NULL,
    title = '(a) Quarterly Mean Predicted Chl-a'
  )

p2 <- ggplot(toplo, aes(x = yr)) +
  geom_col(aes(y = hyresid), fill = 'lightblue') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'darkgrey') +
  facet_wrap(~ qrt, ncol = 4) +
  theme_minimal() +
  labs(
    y = expression('Hydrologic Load ('*10^6*' m'^3*')'),
    x = NULL,
    title = '(b) Quarterly Hydrologic Load'
  )

p3 <- ggplot(toplo, aes(x = hyresid, y = prdresid)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkgrey') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'darkgrey') +
  geom_point() +
  geom_smooth(method = 'lm', color = 'red', se = F, formula = y ~ x) +
  facet_wrap(~ qrtlab, ncol = 4) +
  theme_minimal() +
  labs(
    x = expression('Hydrologic Load Residuals ('*10^6*' m'^3*')'),
    y = 'Predicted Chl-a Residuals (µg/L)',
    title = '(c) Residuals Relationship'
  )

p1 + p2 + p3 + 
  plot_layout(ncol = 1)

# need to extend this to other bay segments