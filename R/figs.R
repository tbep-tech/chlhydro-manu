library(tidyverse)
library(here)
library(patchwork)
library(WRTDStidal)

source(here('R/funcs.R'))

load(file = here('data/mods.RData'))
load(file = here('data/wqdat.RData'))
load(file = here('data/hydat.RData'))

# observed and predicted -------------------------------------------------

toplo <- mods |>
  mutate(
    devexpl = purrr::map(mod, ~ rsq_fun(.x))
  ) |>
  select(bay_segment, mod, devexpl) |>
  unnest(devexpl) |>
  unnest(mod) |>
  mutate(
    chla = exp(res),
    btfit = exp(fit0.5),
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
      barwidth = unit(4, "cm")
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
      barwidth = unit(4, "cm")
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
  mutate(
    annsum = map(mod, function(x) {
      fit <- prdnrmplot(x, annuals = T, plot = F, logspace = F)$fits |>
        mutate(var = 'btfit') |>
        rename(value = fits_value)
      nrm <- prdnrmplot(x, annuals = T, plot = F, logspace = F)$nrms |>
        mutate(var = 'btnorm') |>
        rename(value = nrms_value)
      bind_rows(fit, nrm) |>
        pivot_wider(names_from = var, values_from = value)
    })
  ) |>
  select(-data, -mod) |>
  unnest(annsum) |>
  mutate(
    yr = year(date)
  )

p <- ggplot(toplo, aes(x = yr, y = btfit)) +
  geom_point(aes(fill = "Predicted"), color = 'black') +
  geom_line(aes(y = btnorm, color = "Normalized")) +
  # coord_cartesian(xlim = c(2000, 2024)) +
  scale_color_manual(
    values = c("Normalized" = "tomato1")
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
    y = 'Annual Chl-a (µg/L)',
  )

png(here('figs/prdnrm.png'), width = 7, height = 5, units = 'in', res = 300)
print(p)
dev.off()

# grid plot --------------------------------------------------------------

prdplo <- mods |>
  select(bay_segment, mod) |>
  mutate(
    mod = map(mod, function(x) {
      floobs_rng <- attr(x, 'floobs_rng')
      floscl_rng <- range(x$flo, na.rm = TRUE)
      x$sal <- (x$flo - floscl_rng[1]) /
        diff(floscl_rng) *
        diff(floobs_rng) +
        floobs_rng[1]
      x
    })
  ) |>
  unnest(mod) |>
  select(bay_segment, date, sal, res = fit0.5) |>
  mutate(
    res = exp(res),
    month = lubridate::month(date, label = T, abbr = F),
    yr = lubridate::year(date)
  ) |>
  filter(yr >= 2004 & yr <= 2024 & month %in% month.name[6:11])

grds <- mods |>
  mutate(
    plo = map2(
      bay_segment,
      mod,
      ~ grid_plo(
        .y,
        month = c(6:11),
        years = c(2004, 2024),
        allsal = F,
        ncol = 6,
        sal_fac = 6,
        logspace = F,
        salscl = T,
        col_lim = c(2, 32)
      ) +
        labs(title = .x) +
        geom_line(
          data = prdplo |> filter(bay_segment == .x),
          aes(x = yr, y = sal),
          color = 'black',
          linetype = 'solid',
          linewidth = 0.5,
          inherit.aes = F
        ) +
        geom_point(
          data = prdplo |> filter(bay_segment == .x),
          aes(x = yr, y = sal, fill = res),
          color = 'black',
          size = 3,
          shape = 21
        )
    )
  )

p <- grds$plo[[1]] +
  grds$plo[[2]] +
  grds$plo[[3]] +
  grds$plo[[4]] +
  plot_layout(ncol = 1, guides = 'collect', axis_titles = 'collect') &
  theme(legend.position = 'bottom', axis.text.x = element_text(size = 7))

png(here('figs/gridplo.png'), width = 11, height = 9, units = 'in', res = 300)
print(p)
dev.off()

# salinity response by year ----------------------------------------------

dec_time <- c(1980, 2020)
dys <- c('-02-15', '-05-15', '-08-15', '-11-15')
slc <- expand.grid(dec_time, dys) |>
  unite('doy', Var1, Var2, sep = '') |>
  pull(doy) |>
  as.Date() |>
  sort() |>
  decimal_date()

plos <- mods |>
  mutate(
    salplo = pmap(
      list(mod, data, bay_segment),
      function(mod, data, bay_segment) {
        # get salinity ranges
        prdgrd <- wqdat |>
          mutate(
            yr = year(date),
            anngrp = case_when(
              yr %in% 1975:2004 ~ '1975 - 2004',
              # yr %in% 1986:1995 ~ '1986 - 1995',
              # yr %in% 1991:2009 ~ '1991 - 2009',
              # yr %in% 2006:2015 ~ '2006 - 2015',
              yr %in% 2005:2024 ~ '2005 - 2024',
              T ~ NA_character_
            ),
            mo = month(date),
            qrt = quarter(date)
          ) |>
          filter(mo %in% c(2, 5, 8, 11)) |>
          filter(!is.na(anngrp)) |>
          summarise(
            salmin = quantile(sal, 0.05, na.rm = T),
            salmax = quantile(sal, 0.95, na.rm = T),
            .by = c(mo, anngrp)
          )

        # reformat grid predictions in salinity space
        salgrd <- attr(mod, 'flo_grd')
        prds <- attr(mod, 'fits')$`fit0.5`
        names(prds)[grep('^X', names(prds))] <- paste('sal', salgrd)
        prds <- tidyr::gather(prds, 'sal', 'res', 5:ncol(prds)) |>
          mutate(sal = as.numeric(gsub('^sal ', '', sal))) |>
          select(-month, -day)

        salobs_rng <- attr(mod, 'floobs_rng')
        salscl_rng <- range(prds$sal, na.rm = TRUE)
        prds$sal <- (prds$sal - salscl_rng[1]) /
          diff(salscl_rng) *
          diff(salobs_rng) +
          salobs_rng[1]

        toplo <- prds |>
          mutate(
            anngrp = case_when(
              year %in% 1975:2004 ~ '1975 - 2004',
              # year %in% 1986:1995 ~ '1986 - 1995',
              # year %in% 1991:2009 ~ '1991 - 2009',
              # year %in% 2006:2015 ~ '2006 - 2015',
              year %in% 2005:2024 ~ '2005 - 2024',
              T ~ NA_character_
            ),
            qrt = quarter(date),
            mo = month(date)
          ) |>
          inner_join(prdgrd, by = c('anngrp', 'mo')) |>
          filter(sal >= salmin & sal <= salmax) |>
          mutate(
            res = exp(res),
            qrt = factor(
              qrt,
              levels = 1:4,
              labels = c('JFM', 'AMJ', 'JAS', 'OND')
            )
          ) |>
          summarise(
            btfit = mean(res, na.rm = T),
            btlwr = t.test(res)$conf.int[1], #quantile(res, 0.05, na.rm = T),
            btupr = t.test(res)$conf.int[2], #quantile(res, 0.95, na.rm = T),
            .by = c(sal, anngrp, mo)
          )

        p <- ggplot(
          toplo,
          aes(
            x = sal,
            y = btfit,
            group = anngrp,
            color = anngrp,
            fill = anngrp
          )
        ) +
          geom_ribbon(
            aes(ymin = btlwr, ymax = btupr),
            alpha = 0.3,
            color = NA
          ) +
          geom_line() +
          facet_wrap(~mo, ncol = 4, scales = 'free_y') +
          labs(
            x = 'Salinity (ppth)',
            y = 'Chl-a (µg/L)',
            fill = 'Annual Group',
            color = 'Annual Group',
            subtitle = bay_segment
          ) +
          theme_minimal()

        return(p)
      }
    )
  )

p <- plos$salplo[[1]] +
  plos$salplo[[2]] +
  plos$salplo[[3]] +
  plos$salplo[[4]] +
  plot_layout(ncol = 1, guides = 'collect', axis_titles = 'collect') &
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank()
  )

png(here('figs/salresp.png'), width = 8, height = 8, units = 'in', res = 300)
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
