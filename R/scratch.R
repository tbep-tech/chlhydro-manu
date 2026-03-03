library(tidyverse)
library(patchwork)

load(file = 'data/mods.RData')

dat <- mods |>
  select(bay_segment, data) |>
  unnest(data) |>
  mutate(
    yr = year(date)
  ) |>
  summarise(
    tn_load = sum(tn_load, na.rm = T),
    sal = mean(sal, na.rm = T),
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

dual_theme <- theme(
  legend.position = "none",
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y.left = element_text(color = "tomato"),
  axis.text.y.left = element_text(color = "tomato"),
  axis.ticks.y.left = element_line(color = "tomato"),
  axis.title.y.right = element_text(color = "cornflowerblue"),
  axis.text.y.right = element_text(color = "cornflowerblue"),
  axis.ticks.y.right = element_line(color = "cornflowerblue")
)

hide_left <- theme(axis.title.y.left = element_blank())
hide_right <- theme(axis.title.y.right = element_blank())

make_ld_plot <- function(seg) {
  d <- filter(toplo, bay_segment == seg) |> mutate(tn_load = tn_load / 1000)
  mid_tn <- mean(range(d$tn_load, na.rm = TRUE))
  intcp <- mid_tn # effect = 0 maps to primary center
  d <- mutate(d, ldeff_scaled = slope_ld * ldeff + intcp)
  ylim <- c(mid_tn - rng_tn / 2, mid_tn + rng_tn / 2)
  p <- ggplot(d, aes(x = yr)) +
    geom_line(aes(y = tn_load), color = "tomato") +
    geom_line(aes(y = ldeff_scaled), color = "cornflowerblue") +
    scale_y_continuous(
      name = "tons (x1000)",
      limits = ylim,
      sec.axis = sec_axis(
        ~ (. - intcp) / slope_ld,
        name = "Effect (\u03bcg/L)",
        breaks = ld_breaks
      )
    ) +
    dual_theme
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
  intcp <- mid_sal # effect = 0 maps to primary center
  d <- mutate(d, saleff_scaled = slope_sal * saleff + intcp)
  ylim <- c(mid_sal - rng_sal / 2, mid_sal + rng_sal / 2)
  p <- ggplot(d, aes(x = yr)) +
    geom_line(aes(y = sal), color = "tomato") +
    geom_line(aes(y = saleff_scaled), color = "cornflowerblue") +
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
    dual_theme
  if (seg != 'OTB') {
    p <- p + hide_left
  }
  if (seg != 'LTB') {
    p <- p + hide_right
  }
  p
}

pA <- map(seg_order, make_ld_plot) |> wrap_plots()
pA_sal <- map(seg_order, make_sal_plot) |> wrap_plots()

# Combined 4x4 plot -----------------------------------------------------------

# Rebuild pa and pb from figs.R with correct segment ordering
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
    panel.spacing.x = unit(1.2, "cm")
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
    panel.spacing.x = unit(1.2, "cm")
  )

# Individual plots in segment order
sal_plots <- purrr::map(seg_order, make_sal_plot)
ld_plots <- purrr::map(seg_order, make_ld_plot)

pcombined <- pa /
  wrap_plots(sal_plots, nrow = 1) /
  pb /
  wrap_plots(ld_plots, nrow = 1)
