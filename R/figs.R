library(tidyverse)
library(here)
library(patchwork)
library(mgcv)
library(ggrepel)
library(gratia)
library(maps)
library(ggspatial)
library(sf)
library(purrr)

source(here('R/funcs.R'))

load(file = here('data/mods.RData'))
load(file = here('data/wqdat.RData'))

# map --------------------------------------------------------------------

data("tbseg", package = "tbeptools")
data("tbseglines", package = "tbeptools")

flpoly <- map_data('state', 'florida') %>%
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

chnls <- st_read(here('data/data-raw/Dredge_channels.shp')) %>%
  st_transform(4326)

bbox <- st_bbox(tbseg)

minset <- ggplot() +
  geom_sf(data = flpoly, fill = 'grey', color = NA) +
  geom_sf(data = st_as_sfc(bbox), fill = NA, color = 'black', linewidth = 0.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = '#FFFFFF', colour = 'white'),
    panel.border = element_rect(colour = 'black', fill = 'transparent')
  )

segcent <- tbseg %>%
  st_centroid()

epcpts <- tbeptools::stations %>%
  select(long = Longitude, lat = Latitude) %>%
  unique() %>%
  st_as_sf(coords = c('long', 'lat'), crs = 4326)

thm <- theme(
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text.y = element_text(size = 6),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
  axis.ticks = element_blank()
)

m <- ggplot() +
  ggspatial::annotation_map_tile(
    zoom = 11,
    type = 'cartolight',
    cachedir = system.file("rosm.cache", package = "ggspatial")
  ) +
  geom_sf(
    data = chnls,
    color = 'tomato1',
    fill = 'tomato1',
    alpha = 1,
    inherit.aes = F
  ) +
  geom_sf(data = epcpts, color = 'black', inherit.aes = F) +
  annotation_north_arrow(
    location = 'tl',
    style = north_arrow_orienteering(fill = c('black', 'black'), text_col = NA),
    height = unit(0.5, "cm"),
    width = unit(0.5, "cm")
  ) +
  annotation_scale(location = 'br', text_cex = 1) +
  geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
  geom_sf(data = tbseglines, color = 'black', inherit.aes = F) +
  geom_sf_text(
    data = segcent,
    aes(label = bay_segment),
    size = 5,
    color = 'white',
    inherit.aes = F
  ) +
  # annotation_custom(ggplotGrob(minset), xmin = -9.185e6, xmax = -9.17e6, ymin = 3.22e6, ymax = 3.28e6) +
  annotation_custom(
    ggplotGrob(minset),
    xmin = bbox[3] - 0.1,
    xmax = bbox[3] + 0.015,
    ymin = bbox[4] - 0.1,
    ymax = bbox[4] + 0.06
  ) +
  coord_sf(
    xlim = bbox[c('xmin', 'xmax')],
    ylim = bbox[c('ymin', 'ymax')],
    crs = 4326
  ) +
  thm

png(here('figs/map.png'), width = 4, height = 5, units = 'in', res = 300)
print(m)
dev.off()

# chl, sal, loading over time --------------------------------------------

chlthresh <- tbeptools::targets |>
  select(bay_segment, thresh = chla_thresh) |>
  filter(bay_segment %in% c('OTB', 'HB', 'MTB', 'LTB')) |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  )

ldthresh <- tibble(
  bay_segment = c('OTB', 'HB', 'MTB', 'LTB'),
  thresh = c(486, 1451, 799, 349)
) |>
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  )

toplo <- wqdat |>
  filter(dec_time >= 1985) |>
  select(bay_segment, date, chla, sal, tn_load) |>
  mutate(
    yr = year(date)
  )

toploa <- toplo |>
  select(bay_segment, yr, chla) |>
  summarise(
    val = mean(chla, na.rm = T),
    hiv = t.test(chla)$conf.int[2],
    lov = t.test(chla)$conf.int[1],
    .by = c(bay_segment, yr)
  )

pa <- ggplot(toploa, aes(x = yr, y = val)) +
  geom_point(color = 'darkgrey') +
  geom_errorbar(aes(ymin = lov, ymax = hiv), width = 0, color = 'darkgrey') +
  geom_hline(
    data = chlthresh,
    aes(yintercept = thresh),
    linetype = 'dashed',
    color = 'black'
  ) +
  facet_wrap(~bay_segment, ncol = 4) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()
  ) +
  geom_smooth(method = 'lm', se = F, color = 'black', formula = y ~ x) +
  labs(
    x = NULL,
    y = 'µg/L',
    subtitle = '(a) Chlorophyll-a'
  )

toplob <- toplo |>
  select(bay_segment, yr, sal) |>
  summarise(
    val = mean(sal, na.rm = T),
    hiv = t.test(sal)$conf.int[2],
    lov = t.test(sal)$conf.int[1],
    .by = c(bay_segment, yr)
  )

pb <- ggplot(toplob, aes(x = yr, y = val)) +
  geom_point(color = 'dodgerblue3') +
  geom_errorbar(aes(ymin = lov, ymax = hiv), width = 0, color = 'dodgerblue3') +
  facet_wrap(~bay_segment, ncol = 4) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_blank()
  ) +
  geom_smooth(method = 'lm', se = F, color = 'black', formula = y ~ x) +
  labs(
    x = NULL,
    y = 'ppt',
    subtitle = '(b) Salinity'
  )

toploc <- toplo |>
  select(bay_segment, yr, tn_load) |>
  summarise(
    val = sum(tn_load, na.rm = T),
    .by = c(bay_segment, yr)
  )

pc <- ggplot(toploc, aes(x = yr, y = val / 1000)) +
  geom_point(color = 'tomato1') +
  facet_wrap(~bay_segment, ncol = 4) +
  geom_hline(
    data = ldthresh,
    aes(yintercept = thresh / 1000),
    linetype = 'dashed',
    color = 'black'
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_blank()
  ) +
  geom_smooth(method = 'lm', se = F, color = 'black', formula = y ~ x) +
  labs(
    x = NULL,
    y = 'tons (x1000)',
    subtitle = '(c) Total Nitrogen Load'
  )

p <- pa + pb + pc + plot_layout(ncol = 1)

png(here('figs/obsdat.png'), width = 7, height = 7, units = 'in', res = 300)
print(p)
dev.off()

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
    title = '(b) Predicted Chl-a by Month'
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

png(here('figs/prdnrm.png'), width = 7, height = 5, units = 'in', res = 300)
print(p)
dev.off()

# salinity v load effects ------------------------------------------------

dat <- mods |>
  select(bay_segment, data) |>
  unnest(data) |>
  mutate(yr = year(date)) |>
  summarise(
    tn_load = sum(tn_load, na.rm = TRUE),
    sal = mean(sal, na.rm = TRUE),
    .by = c(bay_segment, yr)
  )

toplo <- mods |>
  select(bay_segment, annsum) |>
  unnest(annsum) |>
  mutate(
    saleff = btfit - btnorm,
    ldeff = btnorm - btnormmd
  ) |>
  select(bay_segment, yr, saleff, ldeff) |>
  left_join(dat, by = c('bay_segment', 'yr'))

seg_order <- c('OTB', 'HB', 'MTB', 'LTB')

# Primary axis: global range (total span across all segments)
rng_tn <- diff(range(toplo$tn_load / 1000, na.rm = TRUE))
rng_sal <- diff(range(toplo$sal, na.rm = TRUE))

# Secondary axis: fixed ranges
rng_ldeff <- diff(range(toplo$ldeff, na.rm = TRUE)) * 1.1 # actual range + 10% buffer
rng_saleff <- 13 # spans -6.5 to 6.5

# Single global slope per pair (same for all segments)
slope_ld <- rng_tn / rng_ldeff
slope_sal <- -rng_sal / rng_saleff # negative: inverse relationship

# Fixed secondary axis breaks — identical across all panels
ld_breaks <- pretty(c(-rng_ldeff / 2, rng_ldeff / 2))
sal_breaks <- pretty(c(-rng_saleff / 2, rng_saleff / 2))

dual_themesal <- theme(
  legend.position = "none",
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y.left = element_text(color = "black"),
  # axis.text.y.left = element_text(color = "black"),
  # axis.ticks.y.left = element_line(color = "black"),
  axis.title.y.right = element_text(color = "cornflowerblue"),
  axis.text.y.right = element_text(color = "cornflowerblue"),
  axis.ticks.y.right = element_line(color = "cornflowerblue")
)

dual_themeld <- theme(
  legend.position = "none",
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y.left = element_text(color = "black"),
  # axis.text.y.left = element_text(color = "black"),
  # axis.ticks.y.left = element_line(color = "black"),
  axis.title.y.right = element_text(color = "tomato1"),
  axis.text.y.right = element_text(color = "tomato1"),
  axis.ticks.y.right = element_line(color = "tomato1")
)

hide_left <- theme(axis.title.y.left = element_blank())
hide_right <- theme(axis.title.y.right = element_blank())

make_ld_plot <- function(seg) {
  d <- filter(toplo, bay_segment == seg) |> mutate(tn_load = tn_load / 1000)
  mid_tn <- mean(range(d$tn_load, na.rm = TRUE))
  intcp <- mid_tn
  d <- mutate(d, ldeff_scaled = slope_ld * ldeff + intcp)
  ylim <- c(mid_tn - rng_tn / 2, mid_tn + rng_tn / 2)
  p <- ggplot(d, aes(x = yr)) +
    geom_line(aes(y = tn_load), color = "grey30") +
    geom_line(aes(y = ldeff_scaled), color = "tomato1") +
    scale_y_continuous(
      name = "tons (x1000)",
      limits = ylim,
      sec.axis = sec_axis(
        ~ (. - intcp) / slope_ld,
        name = "Effect (\u03bcg/L)",
        breaks = ld_breaks
      )
    ) +
    dual_themeld
  if (seg != 'OTB') {
    p <- p + hide_left
  }
  if (seg != 'LTB') {
    p <- p + hide_right
  }
  p
}

make_sal_plot <- function(seg) {
  d <- filter(toplo, bay_segment == seg)
  mid_sal <- mean(range(d$sal, na.rm = TRUE))
  intcp <- mid_sal
  d <- mutate(d, saleff_scaled = slope_sal * saleff + intcp)
  ylim <- c(mid_sal - rng_sal / 2, mid_sal + rng_sal / 2)
  p <- ggplot(d, aes(x = yr)) +
    geom_line(aes(y = sal), color = "grey30") +
    geom_line(aes(y = saleff_scaled), color = "dodgerblue3") +
    scale_y_continuous(
      name = "Salinity (ppt)",
      limits = ylim,
      trans = "reverse",
      sec.axis = sec_axis(
        ~ (. - intcp) / slope_sal,
        name = "Effect (\u03bcg/L)",
        breaks = sal_breaks
      )
    ) +
    dual_themesal
  if (seg != 'OTB') {
    p <- p + hide_left
  }
  if (seg != 'LTB') {
    p <- p + hide_right
  }
  p
}

toplo_eff <- toplo |>
  mutate(bay_segment = factor(bay_segment, levels = seg_order)) |>
  pivot_longer(cols = c(saleff, ldeff), names_to = 'var', values_to = 'val')

pa <- ggplot(toplo_eff |> filter(var == 'saleff'), aes(x = yr, y = val)) +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') +
  geom_point(color = 'dodgerblue3', size = 1) +
  facet_wrap(~bay_segment, ncol = 4) +
  geom_smooth(method = 'lm', se = TRUE, color = 'dodgerblue3') +
  labs(subtitle = '(a) Salinity Effect', x = NULL, y = 'Effect (\u03bcg/L)') +
  theme_minimal() +
  theme(
    axis.ticks.y = element_line(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1.5, 'cm')
  )

pb <- ggplot(toplo_eff |> filter(var == 'ldeff'), aes(x = yr, y = val)) +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') +
  geom_point(color = 'tomato1', size = 1) +
  facet_wrap(~bay_segment, ncol = 4) +
  geom_smooth(method = 'lm', se = TRUE, color = 'tomato1') +
  labs(subtitle = '(b) Load Effect', x = NULL, y = 'Effect (\u03bcg/L)') +
  theme_minimal() +
  theme(
    axis.ticks.y = element_line(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1.5, 'cm')
  )

sal_plots <- purrr::map(seg_order, make_sal_plot)
ld_plots <- purrr::map(seg_order, make_ld_plot)

p <- pa /
  wrap_plots(sal_plots, nrow = 1) /
  pb /
  wrap_plots(ld_plots, nrow = 1)

png(here('figs/salvload.png'), width = 9, height = 8, units = 'in', res = 300)
print(p)
dev.off()

# state space ------------------------------------------------------------

toplo1 <- mods |>
  select(bay_segment, annsum) |>
  unnest(annsum) |>
  # filter(bay_segment == 'OTB') |>
  select(bay_segment, yr, btfit, btnorm, btnormmd) |>
  mutate(
    saleff = btfit - btnorm,
    ldeff = btnorm - btnormmd,
    yrgrp = case_when(
      yr < 1999 ~ '1985 - 1998',
      yr >= 1999 & yr < 2012 ~ '1999 - 2011',
      yr >= 2012 ~ '2012 - 2024'
    ),
    yr = substr(yr, 3, 4)
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
    max.overlaps = 30
  ) +
  scale_color_viridis_d(option = "D", end = 0.8) +
  # coord_equal() +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  facet_wrap(~bay_segment, scales = 'free') +
  labs(
    color = "Year",
    x = 'Salinity Effect (µg/L)',
    y = 'Load Effect (µg/L)',
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
  scale_color_viridis_d(option = "D", end = 0.8) +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  facet_wrap(~bay_segment, scales = 'free') +
  labs(
    color = "Year",
    x = 'Salinity Effect (µg/L)',
    y = 'Load Effect (µg/L)',
    subtitle = '(b) Means (+/- 95% Confidence) by Year Group'
  )

p <- p1 +
  p2 +
  plot_layout(ncol = 1, guides = 'collect', axis_titles = 'collect') &
  thm

png(
  here('figs/statespace.png'),
  width = 7.5,
  height = 9.5,
  units = 'in',
  res = 300
)
print(p)
dev.off()

# sim example for obs conditions -----------------------------------------

load(file = here('data/mods.RData'))
load(file = here('data/simdat.RData'))

act <- mods |>
  filter(bay_segment == 'OTB') |>
  unnest('data') |>
  filter(date >= min(simdat$date) & date <= max(simdat$date)) |>
  select(date, tn_load, sal, dec_time, doy, chla) |>
  mutate(
    yr = year(date)
  ) |>
  summarise(
    chla = mean(chla, na.rm = T),
    .by = yr
  )

toplo1 <- mods |>
  simprd_fun(simdat, nsims = 100, chunk = F, all = T) |>
  filter(bay_segment == 'OTB') |>
  select(-data, -mod, -prds, -annsum, -tst, -exceedsyr, -exceedssum) |>
  unnest('simsyr') |>
  filter(yrs %in% c(1, 50)) |>
  mutate(
    yrs = paste('Year', yrs)
  )

simave <- toplo1 |>
  summarise(
    chla_sim = mean(chla_sim, na.rm = T),
    .by = c('bay_segment', 'ldfac', 'yrs', 'yr')
  )

p1 <- ggplot(toplo1, aes(x = yr, y = chla_sim)) +
  geom_line(alpha = 0.1, aes(group = sim, color = 'Simulations')) +
  geom_line(data = act, aes(y = chla, color = 'Actual'), linewidth = 2) +
  geom_line(
    data = simave,
    aes(y = chla_sim, color = 'Mean', group = NULL),
    linewidth = 2
  ) +
  geom_hline(aes(yintercept = 9.3, color = 'Threshold'), inherit.aes = F) +
  scale_color_manual(
    values = c(
      'Simulations' = 'black',
      'Mean' = 'black',
      'Actual' = 'cornflowerblue',
      'Threshold' = 'red'
    )
  ) +
  facet_grid(ldfac ~ yrs) +
  scale_x_continuous(breaks = seq(2012, 2024, by = 2)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom',
    strip.text.y = element_blank()
  ) +
  labs(
    x = NULL,
    y = 'Annual Chl-a (µg/L)',
    color = NULL,
    subtitle = '(a) Simulations Subset Years 1 and 50\n'
  )

load(file = here('data/simprddat.RData'))

toplo2 <- simprddat |>
  filter(bay_segment == 'OTB') |>
  select(bay_segment, exceedssum) |>
  unnest('exceedssum')

p2 <- ggplot(toplo2, aes(x = yrs, y = avexceeds)) +
  geom_line(aes(color = 'Mean')) +
  coord_cartesian(ylim = c(0, NA)) +
  facet_grid(ldfac ~ .) +
  geom_ribbon(
    aes(
      ymin = avexceeds - sdexceeds,
      ymax = avexceeds + sdexceeds,
      fill = 'Standard Deviation'
    ),
    alpha = 0.5
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c('Mean' = 'black')
  ) +
  scale_fill_manual(
    values = c('Standard Deviation' = 'darkgrey')
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(
    x = 'Years from Simulation Period',
    y = 'Likelihood of Exceeding Threshold',
    subtitle = '(b) Likelihood of Chl-a Exceedence\nOver 50 Years',
    color = NULL,
    fill = NULL
  )

p <- p1 +
  p2 +
  plot_layout(ncol = 2, widths = c(1, 0.5)) &
  theme(legend.position = 'bottom')

png(here('figs/simplo.png'), width = 8, height = 8, units = 'in', res = 300)
print(p)
dev.off()

# simloads ---------------------------------------------------------------

load(file = 'data/simprddat.RData')

toplo <- simprddat |>
  filter(bay_segment == 'OTB') |>
  select(bay_segment, exceedssum) |>
  unnest('exceedssum') |>
  filter(ldfac %in% c('TMDL Load Factor: 0.5', 'TMDL Load Factor: 1')) |>
  mutate(
    ldfac = factor(
      ldfac,
      levels = c('TMDL Load Factor: 0.5', 'TMDL Load Factor: 1')
    )
  )

cols <- c(
  'TMDL Load Factor: 0.5' = '#FC4E2A',
  'TMDL Load Factor: 1' = '#800026'
)

# Year-50 endpoints for horizontal arrow/label annotations
endpoints <- toplo |>
  filter(yrs == 50) |>
  select(ldfac, avexceeds)

# Vertical segment showing reduction between the two scenarios at the y-axis
vert_seg <- tibble(
  y = endpoints$avexceeds[endpoints$ldfac == 'TMDL Load Factor: 1'],
  yend = endpoints$avexceeds[endpoints$ldfac == 'TMDL Load Factor: 0.5'],
) |>
  mutate(
    ymid = (y + yend) / 2,
    label = paste0('-', c(round(y, 1) - round(yend, 1)), '%')
  )

p <- ggplot(
  toplo,
  aes(x = yrs, y = avexceeds, color = ldfac, fill = ldfac)
) +
  geom_ribbon(
    aes(ymin = avexceeds - sdexceeds, ymax = avexceeds + sdexceeds),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_segment(
    data = endpoints,
    aes(x = 50, xend = 0, y = avexceeds, yend = avexceeds),
    linetype = 'dashed',
    linewidth = 0.8,
    arrow = arrow(length = unit(0.25, 'cm'), type = 'closed'),
    show.legend = FALSE
  ) +
  geom_text(
    data = endpoints,
    aes(x = 0, y = avexceeds, label = paste0(round(avexceeds, 1), '%')),
    hjust = -0.3,
    vjust = -0.5,
    size = 3.5,
    show.legend = FALSE
  ) +
  geom_segment(
    data = vert_seg,
    aes(x = 0, xend = 0, y = y, yend = yend),
    color = 'black',
    linetype = 'dashed',
    linewidth = 0.8,
    arrow = arrow(length = unit(0.25, 'cm'), type = 'closed'),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = vert_seg,
    aes(x = 0, y = ymid, label = label),
    hjust = -0.2,
    size = 3.5,
    color = 'black',
    inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  # coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(
    x = 'Years from Simulation Period',
    y = 'Likelihood of Exceeding Threshold',
    color = NULL,
    fill = NULL
  )

png(
  here('figs/simload.png'),
  width = 4.5,
  height = 4.5,
  units = 'in',
  res = 300
)
print(p)
dev.off()

# # grid plot --------------------------------------------------------------

# thresh <- tbeptools::targets |>
#   select(bay_segment, thresh = chla_thresh)

# minyr <- 2010

# prdplo <- mods |>
#   select(bay_segment, prds) |>
#   mutate(prds = purrr::map(prds, as_tibble)) |>
#   unnest(prds) |>
#   select(bay_segment, date, sal, res = btfitmd) |>
#   mutate(
#     month = lubridate::month(date, label = T, abbr = F),
#     yr = lubridate::year(date)
#   ) |>
#   filter(month %in% month.name[6:11] & yr >= minyr) |>
#   left_join(thresh, by = 'bay_segment') |>
#   mutate(
#     exceeds = ifelse(res > thresh, 'above', 'below')
#   )

# grds <- mods |>
#   left_join(thresh, by = 'bay_segment') |>
#   mutate(
#     plo = pmap(
#       list(bay_segment, prds, thresh),
#       function(x, y, z) {
#         grid_plo(
#           y,
#           month = c(6:11),
#           col_lim = c(1, 25),
#           years = c(minyr, 2024),
#           allsal = F,
#           ncol = 6,
#           sal_fac = 6,
#           thresh = z
#         ) +
#           labs(title = x, shape = 'Threshold') +
#           geom_line(
#             data = prdplo |> filter(bay_segment == x),
#             aes(x = yr, y = sal),
#             color = 'black',
#             linetype = 'solid',
#             linewidth = 0.5,
#             inherit.aes = F
#           ) +
#           geom_point(
#             data = prdplo |> filter(bay_segment == x),
#             aes(x = yr, y = sal, fill = res, shape = exceeds),
#             color = 'black',
#             size = 3
#           ) +
#           scale_shape_manual(
#             values = c('above' = 22, 'below' = 21)
#           )
#       }
#     )
#   )

# p <- grds$plo[[1]] +
#   grds$plo[[2]] +
#   grds$plo[[3]] +
#   grds$plo[[4]] +
#   plot_layout(ncol = 1, guides = 'collect', axis_titles = 'collect_y') &
#   theme(legend.position = 'bottom', axis.text.x = element_text(size = 7))

# png(here('figs/gridplo.png'), width = 10, height = 8, units = 'in', res = 300)
# print(p)
# dev.off()

#  supp normalization ----------------------------------------------------

seg <- 'HB'

prds <- mods |>
  filter(bay_segment == seg) |>
  pull(prds) |>
  pluck(1)

salgrd <- attr(prds, 'salgrd')

# reshape a fits attribute to long, aggregate to annual means per salinity level
fits_long <- function(fits_attr, salgrd) {
  fits_attr |>
    pivot_longer(
      cols = starts_with('X'),
      names_to = 'sal_lev',
      values_to = 'chla_pred'
    ) |>
    mutate(sal = salgrd[as.integer(str_remove(sal_lev, 'X'))]) |>
    summarise(chla_pred = mean(chla_pred, na.rm = TRUE), .by = c(year, sal))
}

toplo_act <- fits_long(attr(prds, 'btfits'), salgrd)
toplo_md <- fits_long(attr(prds, 'btfitsmd'), salgrd)

# annual summaries of predicted and normalized values
normsann <- data.frame(prds) |>
  mutate(yr = lubridate::year(date)) |>
  summarise(
    btfit = mean(btfit, na.rm = TRUE),
    btnorm = mean(btnorm, na.rm = TRUE),
    btnormmd = mean(btnormmd, na.rm = TRUE),
    .by = yr
  )

# shared y range across panels A and B
ylim_sal <- range(c(toplo_act$chla_pred, toplo_md$chla_pred), na.rm = TRUE)

thm_top <- theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel A: spaghetti at actual loading, mean = sal-normalized value
p1 <- ggplot(
  toplo_act,
  aes(x = year, y = chla_pred, group = sal, color = sal)
) +
  geom_line(linewidth = 0.7, alpha = 0.7) +
  geom_line(
    data = normsann,
    aes(x = yr, y = btnorm, linetype = 'Mean Predictions'),
    inherit.aes = FALSE,
    color = 'black',
    linewidth = 1.5
  ) +
  scale_color_viridis_c(
    name = 'Salinity (ppt)',
    option = 'mako',
    direction = -1
  ) +
  scale_linetype_manual(
    values = c('Mean Predictions' = 'solid'),
    name = NULL
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(color = 'black', linewidth = 1.5)
    )
  ) +
  scale_y_continuous(limits = ylim_sal) +
  labs(
    y = 'Chl-a (µg/L)',
    subtitle = 'a) Normalized Predictions'
  ) +
  thm_top

# Panel B: spaghetti at mean loading, mean = sal+load-normalized value
p2 <- ggplot(toplo_md, aes(x = year, y = chla_pred, group = sal, color = sal)) +
  geom_line(linewidth = 0.7, alpha = 0.7) +
  geom_line(
    data = normsann,
    aes(x = yr, y = btnormmd, linetype = 'Mean Predictions'),
    inherit.aes = FALSE,
    color = 'black',
    linewidth = 1.5
  ) +
  scale_color_viridis_c(
    name = 'Salinity (ppt)',
    option = 'mako',
    direction = -1
  ) +
  scale_linetype_manual(
    values = c('Mean Predictions' = 'solid'),
    name = NULL
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(color = 'black', linewidth = 1.5)
    )
  ) +
  scale_y_continuous(limits = ylim_sal) +
  labs(
    y = 'Chl-a (µg/L)',
    subtitle = 'b) Normalized Predictions, Mean Load'
  ) +
  thm_top

# Panel C: resulting normalized time series vs predicted
p3 <- ggplot(normsann, aes(x = yr)) +
  geom_point(aes(fill = 'Predicted'), color = 'black', shape = 21) +
  aes(y = btfit) +
  geom_line(aes(y = btnorm, color = 'Normalized', linetype = 'Normalized')) +
  geom_line(aes(
    y = btnormmd,
    color = 'Normalized, Mean Load',
    linetype = 'Normalized, Mean Load'
  )) +
  scale_linetype_manual(
    values = c('Normalized' = 'solid', 'Normalized, Mean Load' = 'dashed')
  ) +
  scale_color_manual(
    values = c(
      'Normalized' = 'tomato1',
      'Normalized, Mean Load' = 'dodgerblue3'
    )
  ) +
  scale_fill_manual(values = c('Predicted' = 'black')) +
  labs(
    x = NULL,
    y = 'Chl-a (µg/L)',
    subtitle = 'c) Final Normalized Time Series',
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank()
  )

# combine: 1 column, 3 rows; p1 legend hidden (p2 repeats it), p3 legend separate
p <- (p1 + theme(legend.position = 'none')) /
  (p2 + theme(legend.position = 'bottom')) /
  (p3 + theme(legend.position = 'bottom')) +
  plot_layout(axis_titles = 'collect_y')

png(here('figs/suppnorm.png'), width = 6, height = 8, units = 'in', res = 300)
print(p)
dev.off()
