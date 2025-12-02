library(tbeptools)
library(tidyverse)
library(mgcv)
library(here)
library(patchwork)

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
  mutate(
    dec_time = decimal_date(date), 
    doy = yday(date),
    chla = ifelse(chla == 0, NA, chla),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) |> 
  arrange(bay_segment, date)

# OTB model --------------------------------------------------------------

tomod <- wqdat |> 
  filter(bay_segment == 'HB')

mod <- gam(chla ~ te(dec_time, doy, sal, bs = c('tp', 'cc', 'tp')), data = tomod, 
knots = list(doy=c(1,366)), family = Gamma(link = 'log'))

mod <- gam(chla ~ s(dec_time, k = 150, bs = 'tp') + 
  s(doy, bs = 'cc') +                       
  s(sal, bs = 'tp') +
  ti(dec_time, doy, bs = c('tp', 'cc')) +
  ti(dec_time, sal, bs = c('tp', 'tp')) +
  ti(sal, doy, bs = c('tp', 'cc')),
  data = tomod,
  knots = list(doy=c(1,366)),
  family = Gamma(link = 'log'), 
)

mod <- gam(chla ~ s(dec_time, k = 300, bs = 'tp') +                   
  s(sal, k = 50, bs = 'tp') +
  ti(dec_time, sal, bs = c('tp', 'tp')),
  data = tomod,
  family = Gamma(link = 'log') 
)

summary(mod)$dev.expl

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
  # filter(dec_time >= 2000) |
  group_nest(bay_segment) |> 
  mutate(
    mod = map(data, ~ gam(chla ~ s(dec_time, k = 150, bs = 'tp') + 
      s(doy, k = 20, bs = 'cc') +                       
      s(sal, bs = 'tp') +
      ti(dec_time, doy, bs = c('tp', 'cc')) +
      ti(dec_time, sal, bs = c('tp', 'tp')) +
      ti(sal, doy, bs = c('tp', 'cc')),
      data = .x,
      knots = list(doy=c(1,366)),
      family = Gamma(link = 'log'), 
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


toplo <- mods |>
  select(bay_segment, annsum) |> 
  unnest(annsum) |> 
  
ggplot(toplo, aes(x = yr, y = btfit)) +
  geom_point() + 
  geom_line(aes(y = btnorm), color = 'red') + 
  geom_line(aes(y = meansalfit), color = 'blue', linetype = 'dashed') +
  # coord_cartesian(xlim = c(2000, 2024)) +
  facet_wrap(~ bay_segment, scales = 'free_y') +
  theme_minimal()

toplo <- mods |>
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds)

ggplot(toplo, aes(x = date, y = chla)) +
  geom_point() +
  geom_line(aes(y = btfit), color = 'red') +
  # geom_line(aes(y = btnorm), color = 'blue') +
  facet_wrap(~ bay_segment, scales = 'free_y') +
  scale_x_date(limits = c(as.Date('2004-01-01'), as.Date('2024-12-31'))) +
  theme_minimal()

actplo <- wqdat |> 
  mutate(
    yr = lubridate::year(date),
    month = lubridate::month(date, label = T, abbr = F)
  ) |> 
  filter(yr >= 2004 & yr <= 2024) |>
  summarise(
    sal = mean(sal, na.rm = T), 
    res = mean(chla, na.rm = T),
    .by = c(bay_segment, yr, month)
  ) |> 
  filter(month %in% month.name[6:11])

grds <- mods |> 
  mutate(
    plo = map2(bay_segment, prds, ~ grid_plo(.y, month = c(6:11), years = c(2004, 2024), col_lim = c(1, 35), allsal = F, ncol = 6, sal_fac = 6) + labs(title = .x) + 
      geom_line(data = actplo |> filter(bay_segment == .x), aes(x = yr, y = sal), color = 'black', linetype = 'solid', linewidth = 0.5, inherit.aes = F) + 
      geom_point(data = actplo |> filter(bay_segment == .x), aes(x = yr, y = sal, fill = res), color = 'black', size = 3, shape = 21)
    )
  )

grds$plo[[1]] + grds$plo[[2]] + grds$plo[[3]] + grds$plo[[4]] + 
  plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = 'bottom')
