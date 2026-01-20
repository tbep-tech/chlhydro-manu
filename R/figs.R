library(tidyverse)
library(here)
library(patchwork)
library(mgcv)
library(ggrepel)

source(here('R/funcs.R'))

load(file = here('data/mods.RData'))
load(file = here('data/wqdat.RData'))
load(file = here('data/hydat.RData'))

# observed and predicted -------------------------------------------------

toplo <- mods |>
  mutate(
    devexpl = purrr::map(mod, ~ summary(.x)$dev.expl)
  ) |>
  select(bay_segment, prds, devexpl) |>
  unnest(devexpl) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |>
  mutate(
    striplab = paste0(round(devexpl * 100, 0), '% dev. expl.')
  )

fac <- toplo |>
  select(bay_segment, striplab) |>
  distinct()

toplo <- toplo |>
  mutate(
    striplab = factor(
      bay_segment,
      levels = fac$bay_segment,
      labels = fac$striplab
    ),
    doydum = make_date(2000, month = month(date), day = day(date)),
    yr = year(date)
  )

p1 <- ggplot(toplo, aes(x = date, y = chla)) +
  geom_point(size = 1, aes(fill = 'Observed'), pch = 21, color = 'darkgrey') +
  geom_line(aes(y = btfit, color = yr)) +
  # coord_cartesian(xlim = c(2000, 2024)) +
  scale_y_log10() +
  scale_color_viridis_c() +
  scale_fill_manual(values = c('Observed' = 'darkgrey')) +
  facet_wrap(~bay_segment, scales = 'free', ncol = 1) +
  guides(
    color = guide_colorbar(
      barheight = unit(0.25, "cm"),
      barwidth = unit(6, "cm")
    )
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(
    x = NULL,
    fill = NULL,
    color = 'Year',
    y = 'µg/L',
    title = '(a) Observed and Predicted Chl-a'
  )

p2 <- ggplot(toplo, aes(x = doydum, y = btfit, group = yr, color = yr)) +
  geom_point(
    size = 1,
    aes(fill = 'Observed', y = chla),
    pch = 21,
    color = 'darkgrey'
  ) +
  geom_line() +
  # geom_smooth(se = F, method = "gam", formula = y ~ s(x, bs = "tp", k = 10)) +
  scale_x_date(date_labels = '%b', date_breaks = '1 month') +
  # coord_cartesian(xlim = c(2000, 2024)) +
  scale_y_log10() +
  scale_color_viridis_c() +
  scale_fill_manual(values = c('Observed' = 'darkgrey')) +
  facet_wrap(~striplab, scales = 'free', ncol = 1) +
  theme_minimal() +
  guides(
    color = guide_colorbar(
      barheight = unit(0.25, "cm"),
      barwidth = unit(6, "cm")
    )
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(
    x = NULL,
    color = "Year",
    fill = NULL,
    y = 'µg/L',
    title = '(b) Predicted Chl-a by Day of Year'
  )

p <- p1 +
  p2 +
  plot_layout(ncol = 2, axis_titles = 'collect_y', guides = 'collect') &
  theme(legend.position = 'bottom')

png(here('figs/obsprd.png'), width = 8, height = 9, units = 'in', res = 300)
print(p)
dev.off()

# predicted and normalized annual ----------------------------------------

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
  # coord_cartesian(xlim = c(2000, 2024)) +
  scale_linetype_manual(
    values = c("Normalized" = "solid", "Normalized, Mean Load" = "dashed")
  ) +
  scale_color_manual(
    values = c("Normalized" = "tomato1", "Normalized, Mean Load" = "blue")
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

png(here('figs/prdnrm.png'), width = 7, height = 5, units = 'in', res = 300)
print(p)
dev.off()

# salinity vs load effects ------------------------------------------------

toplo1 <- mods |>
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
p2 <- ggplot(toplo, aes(x = saleff, y = ldeff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(alpha = 0) +
  geom_point(
    data = mns,
    aes(x = saleffmn, y = ldeffmn, color = yrgrp),
    size = 4,
    show.legend = F
  ) +
  geom_errorbar(
    data = mns,
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
    data = mns,
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

png(
  here('figs/salvload.png'),
  width = 7.5,
  height = 9.5,
  units = 'in',
  res = 300
)
print(p)
dev.off()

# grid plot --------------------------------------------------------------

thresh <- tbeptools::targets |>
  select(bay_segment, thresh = chla_thresh)

minyr <- 2010

prdplo <- mods |>
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |>
  select(bay_segment, date, sal, res = btfitmd) |>
  mutate(
    month = lubridate::month(date, label = T, abbr = F),
    yr = lubridate::year(date)
  ) |>
  filter(month %in% month.name[6:11] & yr >= minyr) |>
  left_join(thresh, by = 'bay_segment') |>
  mutate(
    exceeds = ifelse(res > thresh, 'above', 'below')
  )

grds <- mods |>
  left_join(thresh, by = 'bay_segment') |>
  mutate(
    plo = pmap(
      list(bay_segment, prds, thresh),
      function(x, y, z) {
        grid_plo(
          y,
          month = c(6:11),
          col_lim = c(1, 25),
          years = c(minyr, 2024),
          allsal = F,
          ncol = 6,
          sal_fac = 6,
          thresh = z,
          ldmod = 'btfitsmd'
        ) +
          labs(title = x, shape = 'Threshold') +
          geom_line(
            data = prdplo |> filter(bay_segment == x),
            aes(x = yr, y = sal),
            color = 'black',
            linetype = 'solid',
            linewidth = 0.5,
            inherit.aes = F
          ) +
          geom_point(
            data = prdplo |> filter(bay_segment == x),
            aes(x = yr, y = sal, fill = res, shape = exceeds),
            color = 'black',
            size = 3
          ) +
          scale_shape_manual(
            values = c('above' = 22, 'below' = 21)
          )
      }
    )
  )

p <- grds$plo[[1]] +
  grds$plo[[2]] +
  grds$plo[[3]] +
  grds$plo[[4]] +
  plot_layout(ncol = 1, guides = 'collect', axis_titles = 'collect_y') &
  theme(legend.position = 'bottom', axis.text.x = element_text(size = 7))

png(here('figs/gridplo.png'), width = 10, height = 8, units = 'in', res = 300)
print(p)
dev.off()

# salinity response by year ----------------------------------------------

load(file = here('data/simprddat1721.RData'))
load(file = here('data/simprddat9294.RData'))

toploa <- simprddat9294 |>
  select(bay_segment, exceedssum) |>
  unnest('exceedssum') |>
  mutate(
    Period = '1992 - 1994'
  )
toplob <- simprddat1721 |>
  select(bay_segment, exceedssum) |>
  unnest('exceedssum') |>
  mutate(
    Period = '2017 - 2021'
  )
toplo <- bind_rows(toploa, toplob) |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  )

p <- ggplot(toplo, aes(x = yrs, y = avexceeds, color = Period, fill = Period)) +
  geom_line() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_grid(bay_segment ~ ldfac, scales = 'free_y') +
  geom_ribbon(
    aes(
      ymin = avexceeds - sdexceeds,
      ymax = avexceeds + sdexceeds,
    ),
    alpha = 0.2
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(
    x = 'Years from simulation',
    y = 'Likelihood of exceeding threshold'
  )

png(here('figs/simplo.png'), width = 8, height = 8, units = 'in', res = 300)
print(p)
dev.off()

# hydrology vs norm ------------------------------------------------------

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
    qrt = factor(
      lubridate::quarter(date),
      levels = 1:4,
      labels = c('JFM', 'AMJ', 'JAS', 'OND')
    ),
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
  facet_wrap(~qrt, ncol = 4) +
  theme_minimal() +
  labs(
    y = 'µg/L',
    x = NULL,
    title = '(a) Quarterly Mean Predicted Chl-a'
  )

p2 <- ggplot(toplo, aes(x = yr)) +
  geom_col(aes(y = hyresid), fill = 'lightblue') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'darkgrey') +
  facet_wrap(~qrt, ncol = 4) +
  theme_minimal() +
  labs(
    y = expression(10^6 ~ m^3 * ' / Quarter'),
    x = NULL,
    title = '(b) Quarterly Hydrologic Anomalies'
  )

p3 <- ggplot(toplo, aes(x = hyresid, y = prdresid)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkgrey') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'darkgrey') +
  geom_point(color = 'darkgrey') +
  geom_smooth(method = 'lm', color = 'black', se = F, formula = y ~ x) +
  facet_wrap(~qrtlab, ncol = 4) +
  theme_minimal() +
  labs(
    x = expression('Hydrologic Anomalies (' * 10^6 ~ m^3 * ' / quarter)'),
    y = 'Predicted - Normalized Chl-a (µg/L)',
    title = '(c) Chlorophyll vs Hydrology'
  )

p <- p1 +
  p2 +
  p3 +
  plot_layout(ncol = 1) &
  theme(
    panel.grid.minor = element_blank()
  )

png(here('figs/hydnrm.png'), width = 8, height = 8, units = 'in', res = 300)
print(p)
dev.off()

# sim example for obs conditions -----------------------------------------

load(file = here('data/mods.RData'))
load(file = here('data/simdat1721.RData'))

act <- mods |>
  filter(bay_segment == 'OTB') |>
  unnest('data') |>
  filter(date >= as.Date('2017-01-01') & date <= as.Date('2021-12-31')) |>
  select(date, tn_load, sal, dec_time, doy, chla) |>
  mutate(
    yr = year(date)
  ) |>
  summarise(
    chla = mean(chla, na.rm = T),
    .by = yr
  )

toplo <- simprd_fun(mods, simdat1721, nsims = 100, all = T) |>
  filter(bay_segment == 'OTB') |>
  unnest('simsyr') |>
  filter(yrs %in% c(1, 50) & ldfac == 'Actual Load') |>
  mutate(
    yrs = paste('Year', yrs)
  )

p <- ggplot(toplo, aes(x = yr, y = chla_sim)) +
  geom_line(alpha = 0.2, aes(group = sim, color = 'Simulations')) +
  geom_smooth(
    formula = y ~ x,
    aes(color = 'Simulations Mean'),
    se = F,
    method = 'loess',
    linewidth = 2
  ) +
  geom_hline(aes(yintercept = 9.3, color = 'Threshold'), inherit.aes = F) +
  geom_line(data = act, aes(y = chla, color = 'Actual'), linewidth = 2) +
  scale_color_manual(
    values = c(
      'Simulations' = 'black',
      'Simulations Mean' = 'black',
      'Actual' = 'cornflowerblue',
      'Threshold' = 'red'
    )
  ) +
  facet_wrap(~yrs, ncol = 1) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(
    x = NULL,
    y = 'Annual Chl-a (µg/L)',
    color = NULL
  )

png(here('figs/simex.png'), width = 5, height = 6, units = 'in', res = 300)
print(p)
dev.off()
