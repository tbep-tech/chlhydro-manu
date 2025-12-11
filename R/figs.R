library(tidyverse)
library(here)
library(patchwork)

load(file = here('data/mods.RData'))

# observed and predicted -------------------------------------------------


toplo <- mods |> 
  mutate(
    devexpl = purrr::map(mod, ~summary(.x)$dev.expl)
  ) |> 
  select(bay_segment, prds, devexpl) |> 
  unnest(devexpl) |> 
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |> 
  mutate(
    striplab = paste0(round(devexpl*100, 0), '% dev. expl.')
  )

fac <- toplo |> 
  select(bay_segment, striplab) |>
  distinct()

toplo <- toplo |>
  mutate(
    striplab = factor(bay_segment, levels = fac$bay_segment, labels = fac$striplab),
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
  facet_wrap(~ bay_segment, scales = 'free', ncol = 1) +
  guides(color = guide_colorbar(barheight = unit(0.25, "cm"),
                                barwidth = unit(4, "cm"))) +
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
  geom_point(size = 1, aes(fill = 'Observed', y = chla), pch = 21, color = 'darkgrey') + 
  geom_line() + 
  # geom_smooth(se = F, method = "gam", formula = y ~ s(x, bs = "tp", k = 10)) +
  scale_x_date(date_labels = '%b', date_breaks = '1 month') +
  # coord_cartesian(xlim = c(2000, 2024)) +
  scale_y_log10() + 
  scale_color_viridis_c() +
  scale_fill_manual(values = c('Observed' = 'darkgrey')) +
  facet_wrap(~ striplab, scales = 'free', ncol = 1) +
  theme_minimal() +
  guides(color = guide_colorbar(barheight = unit(0.25, "cm"),
                                barwidth = unit(4, "cm"))) +
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

p <- p1 + p2 + plot_layout(ncol = 2, axis_titles = 'collect_y', guides = 'collect') & 
  theme(legend.position = 'bottom')

png(here('figs/obsprd.png'), width = 8, height = 10, units = 'in', res = 300)
print(p)
dev.off()