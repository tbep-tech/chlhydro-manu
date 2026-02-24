library(tbeptools)
library(tidyverse)
library(mgcv)
library(here)
library(gratia)

source(here('R/funcs.R'))

# data prep --------------------------------------------------------------

load(file = here('data/wqdat.RData'))

tomod <- wqdat |>
  filter(dec_time >= 1985) |>
  group_by(bay_segment) |>
  group_nest() |>
  mutate(
    data = purrr::map(
      data,
      function(x) {
        salmod <- gam(
          sal ~ s(doy, bs = 'cc'),
          data = x,
          method = 'REML',
          knots = list(doy = c(0, 366))
        )
        x$sal <- residuals(salmod)
        tnmod <- gam(
          tn_load ~ s(doy, bs = 'cc') + s(sal),
          data = x,
          method = 'REML',
          knots = list(doy = c(0, 366))
        )
        x$tn_load <- residuals(tnmod)
        x
      }
    )
  )
mods <- tomod |>
  mutate(
    mod = purrr::map(
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
    prds = purrr::map2(data, mod, pred_fun),
    annsum = purrr::map(
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
          .by = c(yr)
        )
    )
  )


toplo <- mods |>
  select(bay_segment, annsum) |>
  unnest(annsum)

p <- ggplot(toplo, aes(x = yr, y = btfit)) +
  geom_point(aes(fill = "Predicted"), color = 'black') +
  geom_line(aes(y = btnorm, color = "Normalized", linetype = "Normalized")) +
  geom_line(aes(
    y = btnormmd,
    color = "Normalized, Mean Load",
    linetype = "Normalized, Mean Load"
  )) +
  scale_linetype_manual(
    values = c("Normalized" = "solid", "Normalized, Mean Load" = "dashed")
  ) +
  scale_color_manual(
    values = c(
      "Normalized" = "tomato1",
      "Normalized, Mean Load" = "dodgerblue3"
    )
  ) +
  scale_fill_manual(values = c("Predicted" = "black")) +
  facet_wrap(~bay_segment, scales = 'free_y') +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    color = NULL,
    fill = NULL,
    linetype = NULL,
    y = 'Annual Chl-a (µg/L)',
  )
