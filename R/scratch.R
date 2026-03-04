library(tidyverse)
library(patchwork)

# simplo (b) overlay — TMDL and half-TMDL superimposed -------------------

toplo2_sub <- simprddat1524 |>
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

load_colors2 <- c(
  'TMDL Load Factor: 0.5' = '#FC4E2A',
  'TMDL Load Factor: 1' = '#800026'
)

# Year-50 endpoints for horizontal arrow/label annotations
endpoints <- toplo2_sub |>
  filter(yrs == 50) |>
  select(ldfac, avexceeds)

# Vertical segment showing reduction between the two scenarios at the y-axis
vert_seg <- tibble(
  y = endpoints$avexceeds[endpoints$ldfac == 'TMDL Load Factor: 1'],
  yend = endpoints$avexceeds[endpoints$ldfac == 'TMDL Load Factor: 0.5'],
) |>
  mutate(
    ymid = (y + yend) / 2,
    label = paste0('-', round((y - yend) * 100), '%')
  )

p_overlay2 <- ggplot(
  toplo2_sub,
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
    linewidth = 0.4,
    arrow = arrow(length = unit(0.2, 'cm'), type = 'closed'),
    show.legend = FALSE
  ) +
  geom_text(
    data = endpoints,
    aes(x = 0, y = avexceeds, label = paste0(round(avexceeds * 100), '%')),
    hjust = 0,
    vjust = -0.5,
    size = 4,
    show.legend = FALSE
  ) +
  geom_segment(
    data = vert_seg,
    aes(x = 0, xend = 0, y = y, yend = yend),
    color = 'black',
    linetype = 'dashed',
    linewidth = 0.4,
    arrow = arrow(length = unit(0.2, 'cm'), type = 'closed'),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = vert_seg,
    aes(x = 0, y = ymid, label = label),
    hjust = -0.2,
    size = 4,
    color = 'black',
    inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = load_colors2) +
  scale_fill_manual(values = load_colors2) +
  theme_minimal() +
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

p_overlay2
