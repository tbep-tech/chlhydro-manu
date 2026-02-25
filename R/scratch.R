library(tbeptools)
library(tidyverse)
library(mgcv)
library(here)
library(gratia)
library(ggrepel)
library(patchwork)
library(flextable)

source(here('R/funcs.R'))

# data prep --------------------------------------------------------------

load(file = here('data/wqdat.RData'))
load(file = here('data/lddat.RData'))

lddatlag <- lddat |>
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
    tn_load0 = tn_load,
    tn_load1 = tn_load + lag1,
    tn_load2 = tn_load + lag1 + lag2,
    tn_load3 = tn_load + lag1 + lag2 + lag3
  ) |>
  select(-tn_load, -lag1, -lag2, -lag3)

wqdat <- wqdat |>
  select(-tn_load) |>
  left_join(lddatlag, by = c('bay_segment', 'date'))

mods <- wqdat |>
  filter(dec_time >= 1985) |>
  pivot_longer(
    cols = starts_with('tn_load'),
    names_to = 'tn_load_lag',
    values_to = 'tn_load'
  ) |>
  group_nest(bay_segment, tn_load_lag) |>
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
    ),
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


totab <- mods |>
  select(bay_segment, tn_load_lag, mod) |>
  mutate(
    summ = purrr::map(mod, function(x) {
      summmod <- summary(x)
      GCV <- round(summmod$sp.criterion[[1]], 0)
      devexpl <- round(summmod$dev.expl, 2)

      out <- bind_cols(GCV = GCV, devexpl = devexpl) |>
        rename(
          `Dev. Expl.` = devexpl
        )

      return(out)
    })
  ) |>
  select(-mod) |>
  unnest('summ') |>
  mutate(
    bay_segment = ifelse(
      duplicated(bay_segment),
      '',
      as.character(bay_segment)
    ),
    tn_load_lag = case_when(
      tn_load_lag == 'tn_load0' ~ 'None',
      tn_load_lag == 'tn_load1' ~ '1',
      tn_load_lag == 'tn_load2' ~ '2',
      tn_load_lag == 'tn_load3' ~ '3'
    )
  )

suppgamtab <- totab |>
  flextable() |>
  set_header_labels(bay_segment = 'Bay Segment', tn_load_lag = 'TN Load Lag') |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  colformat_double(
    j = 3,
    big.mark = '',
    digits = 0
  ) |>
  autofit() |>
  hline(i = 4) |>
  hline(i = 8) |>
  hline(i = 12)

toplo <- mods |>
  filter(tn_load_lag == 'tn_load3') |>
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


toplo1 <- mods |>
  filter(tn_load_lag == 'tn_load3') |>
  select(bay_segment, annsum) |>
  unnest(annsum) |>
  # filter(bay_segment == 'OTB') |>
  select(bay_segment, yr, btfit, btnorm, btnormmd) |>
  mutate(
    saleff = btfit - btnorm,
    ldeff = btnorm - btnormmd,
    yrgrp = case_when(
      yr < 2000 ~ '1985 - 1999',
      yr >= 2000 & yr < 2015 ~ '2000 - 2014',
      yr >= 2015 ~ '2015 - 2024'
    )
  )

toplo2 <- toplo1 |>
  summarise(
    saleffmn = mean(saleff),
    saleffhi = t.test(saleff)$conf.int[2],
    salefflo = t.test(saleff)$conf.int[1],
    ldeffmn = mean(ldeff),
    ldeffhi = t.test(ldeff)$conf.int[2],
    ldefflo = t.test(ldeff)$conf.int[1],
    .by = c(yrgrp, bay_segment)
  )

thm <- theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )

# create path by year
set.seed(123)
p1 <- ggplot(toplo1, aes(x = saleff, y = ldeff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_path(aes(group = 1), color = 'gray80', linewidth = 1) +
  geom_point(aes(color = yrgrp), size = 2) +
  geom_text_repel(
    aes(label = yr, color = yrgrp),
    size = 2,
    max.overlaps = 20
  ) +
  scale_color_viridis_d(option = "D", end = 0.8) +
  # scale_fill_viridis_d(option = "D", end = 0.8) +
  # coord_equal() +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  facet_wrap(~bay_segment, scales = 'free') +
  labs(
    color = "Year",
    x = 'Salinity Effect µg/L (+ is lower salinity)',
    y = 'Load Effect µg/L (+ is higher load)',
    subtitle = '(a) By Year'
  )

# create yrgrp means
set.seed(123)
p2 <- ggplot(toplo1, aes(x = saleff, y = ldeff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(alpha = 0) +
  geom_point(
    data = toplo2,
    aes(x = saleffmn, y = ldeffmn, color = yrgrp),
    size = 4,
    show.legend = F
  ) +
  geom_errorbar(
    data = toplo2,
    aes(
      x = saleffmn,
      y = ldeffmn,
      xmin = salefflo,
      xmax = saleffhi,
      color = yrgrp
    ),
    orientation = 'y',
    show.legend = F,
    lineend = 'square',
    width = 0
  ) +
  geom_errorbar(
    data = toplo2,
    aes(
      x = saleffmn,
      y = ldeffmn,
      ymin = ldefflo,
      ymax = ldeffhi,
      color = yrgrp
    ),
    show.legend = F,
    lineend = 'square',
    width = 0
  ) +
  geom_text_repel(
    aes(label = yr),
    alpha = 0,
    max.overlaps = 20
  ) +
  scale_color_viridis_d(option = "D", end = 0.8) +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  facet_wrap(~bay_segment, scales = 'free') +
  labs(
    color = "Year",
    x = 'Salinity Effect µg/L (+ is lower salinity)',
    y = 'Load Effect µg/L (+ is higher load)',
    subtitle = '(b) Means (+/- 95% Confidence) by Year Group'
  )

p <- p1 +
  p2 +
  plot_layout(ncol = 1, guides = 'collect', axis_titles = 'collect') &
  thm
