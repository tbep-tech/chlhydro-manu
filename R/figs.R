library(tidyverse)
library(here)
library(patchwork)
library(mgcv)
library(ggrepel)
library(gratia)
library(maps)
library(ggspatial)
library(sf)

source(here('R/funcs.R'))

load(file = here('data/mods.RData'))
load(file = here('data/wqdat.RData'))

# map --------------------------------------------------------------------

flpoly <- map_data('state', 'florida') %>%
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

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

toplo <- mods |>
  select(bay_segment, annsum) |>
  unnest(annsum) |>
  mutate(
    saleff = btfit - btnorm,
    ldeff = btnorm - btnormmd
  ) |>
  select(bay_segment, yr, saleff, ldeff) |>
  pivot_longer(
    cols = c(saleff, ldeff),
    names_to = 'var',
    values_to = 'val'
  )

pa <- ggplot(toplo |> filter(var == 'saleff'), aes(x = yr, y = val)) +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') +
  geom_point(color = 'dodgerblue3', size = 1) +
  facet_wrap(~bay_segment, ncol = 4) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'dodgerblue3'
  ) +
  labs(
    subtitle = '(a) Salinity Effect',
    x = NULL,
    y = 'µg/L (+ is lower salinity)',
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank()
  )
pb <- ggplot(toplo |> filter(var == 'ldeff'), aes(x = yr, y = val)) +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') +
  geom_point(color = 'tomato1', size = 1) +
  facet_wrap(~bay_segment, ncol = 4) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'tomato1'
  ) +
  labs(
    subtitle = '(b) Load Effect',
    x = NULL,
    y = 'µg/L (+ is higher load)',
  ) +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    panel.grid.minor = element_blank()
  )

p <- pa + pb + plot_layout(ncol = 1, guides = 'collect')

png(here('figs/salvload.png'), width = 7, height = 6, units = 'in', res = 300)
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
load(file = here('data/simdat1524.RData'))

act <- mods |>
  filter(bay_segment == 'OTB') |>
  unnest('data') |>
  filter(date >= as.Date('2015-01-01') & date <= as.Date('2024-12-31')) |>
  select(date, tn_load, sal, dec_time, doy, chla) |>
  mutate(
    yr = year(date)
  ) |>
  summarise(
    chla = mean(chla, na.rm = T),
    .by = yr
  )

toplo1 <- simprd_fun(mods, simdat1524, nsims = 100, chunk = F, all = T) |>
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
  scale_x_continuous(breaks = seq(2015, 2024, by = 2)) +
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

load(file = here('data/simprddat1524.RData'))

toplo2 <- simprddat1524 |>
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
