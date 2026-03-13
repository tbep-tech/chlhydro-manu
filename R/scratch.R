library(tidyverse)
library(patchwork)
library(here)

load(file = here('data/mods.RData'))

seg <- 'OTB'

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

# combine: A and B share the viridis legend; C keeps its own color legend
(p1 + p2 + plot_layout(guides = 'collect')) /
  p3 &
  theme(legend.position = 'bottom')
