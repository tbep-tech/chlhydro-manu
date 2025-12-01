library(tbeptools)
library(tidyverse)
library(mgcv)
library(here)

source(here('R/funcs.R'))

# data prep --------------------------------------------------------------

# water quality
wqdat <- epcdata |> 
  filter(yr > 1974 & yr <= 2024) |> 
  select(bay_segment, SampleTime, tn, chla, Sal_Top_ppth, Sal_Mid_ppth, Sal_Bottom_ppth) |>
  mutate(
    sal = mean(c(Sal_Top_ppth, Sal_Mid_ppth, Sal_Bottom_ppth), na.rm = T), 
    .by = c(bay_segment, SampleTime)
  ) |> 
  mutate(
    date = floor_date(as.Date(SampleTime), unit = 'month')
  ) |> 
  summarise(
    tn = median(tn, na.rm = T), 
    chla = median(chla, na.rm = T), 
    sal = median(sal, na.rm = T),
    .by = c(bay_segment, date)
  ) |> 
  arrange(bay_segment, date) |> 
  mutate(
    dec_time = decimal_date(date), 
    doy = yday(date),
    chla = ifelse(chla == 0, NA, chla)
  )

# OTB model --------------------------------------------------------------

tomod <- wqdat |> 
  filter(bay_segment == 'HB')

mod <- gam(chla ~ te(dec_time, doy, sal, bs = c('tp', 'cc', 'tp')), data = tomod, 
knots = list(doy=c(1,366)), family = Gamma(link = 'log'))

mod <- gam(chla ~ s(dec_time, bs = 'tp') + 
  s(doy, bs = 'cc') +                       
  s(sal, bs = 'tp') +
  ti(dec_time, doy, sal, bs = c('tp', 'cc', 'tp')),
  data = tomod,
  knots = list(doy=c(1,366)),
  family = Gamma(link = 'log'), 
)

summary(mod)$dev.expl

prds <- pred_fun(tomod, mod)

ggplot(prds, aes(x = dec_time, y = chla)) +
  geom_point() +
  geom_line(aes(y = btfit), color = 'red')


grid_plo(prds, month = 'all', allsal = T)
grid_plo(prds)
grid_plo(prds, month = c(5:10), years = c(2000, 2024), allsal = T, ncol = 6, sal_fac = 6)


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
  # filter(dec_time >= 2000) |
  group_nest(bay_segment) |> 
  mutate(
    mod = map(data, ~ gam(chla ~ s(dec_time, k = 100, bs = 'tp') + 
      s(doy, bs = 'cc') +                       
      s(sal, bs = 'tp') +
      ti(dec_time, doy, sal, bs = c('tp', 'cc', 'tp')),
      data = .x,
      knots = list(doy=c(1,366)),
      family = Gamma(link = 'log')
    )),
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


toplo <- mods |>
  select(bay_segment, annsum) |> 
  unnest(annsum)

ggplot(toplo, aes(x = yr, y = btfit)) +
  geom_point() + 
  geom_line(aes(y = btnorm), color = 'red') + 
  geom_line(aes(y = meansalfit), color = 'blue', linetype = 'dashed') +
  # coord_cartesian(xlim = c(2000, 2024)) +
  facet_wrap(~ bay_segment, scales = 'free_y') +
  theme_minimal()
