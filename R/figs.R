library(tidyverse)
library(here)
library(patchwork)
library(mgcv)

source(here('R/funcs.R'))

load(file = here('data/mods.RData'))
load(file = here('data/wqdat.RData'))

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
  geom_line(aes(y = meansalfit, color = "Mean Salinity", linetype = "Mean Salinity")) +
  # coord_cartesian(xlim = c(2000, 2024)) +
  scale_linetype_manual(values = c("Normalized" = "solid", "Mean Salinity" = "dashed")) +
  scale_color_manual(values = c("Normalized" = "tomato1", "Mean Salinity" = "blue")) +
  scale_fill_manual(values = c("Predicted" = "black")) +
  facet_wrap(~ bay_segment, scales = 'free_y') +
  theme_minimal() + 
  theme(legend.position = 'bottom') +
  labs(
    x = NULL, 
    color = NULL, 
    fill = NULL,
    linetype = NULL,
    y  = 'Annual Chl-a (µg/L)',
  )

png(here('figs/prdnrm.png'), width = 7, height = 5, units = 'in', res = 300)
print(p)
dev.off()


# grid plot --------------------------------------------------------------

prdplo <- mods |> 
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |> 
  select(bay_segment, date, sal, res = btfit) |> 
  mutate(
    month = lubridate::month(date, label = T, abbr = F),
    yr = lubridate::year(date)
  ) |> 
  filter(yr >= 2004 & yr <= 2024 & month %in% month.name[6:11])

grds <- mods |> 
  mutate(
    plo = map2(bay_segment, prds, ~ grid_plo(.y, month = c(6:11), years = c(2004, 2024), col_lim = c(1, 30), allsal = F, ncol = 6, sal_fac = 6) + labs(title = .x) + 
      geom_line(data = prdplo |> filter(bay_segment == .x), aes(x = yr, y = sal), color = 'black', linetype = 'solid', linewidth = 0.5, inherit.aes = F) + 
      geom_point(data = prdplo |> filter(bay_segment == .x), aes(x = yr, y = sal, fill = res), color = 'black', size = 3, shape = 21)
    )
  )

p <- grds$plo[[1]] + grds$plo[[2]] + grds$plo[[3]] + grds$plo[[4]] + 
  plot_layout(ncol = 1, guides = 'collect', axis_titles = 'collect_y') & theme(legend.position = 'bottom', axis.text.x = element_text(size = 7))

png(here('figs/gridplo.png'), width = 11, height = 9, units = 'in', res = 300)
print(p)
dev.off()
