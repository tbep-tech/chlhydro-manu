library(tidyverse)
library(mgcv)
library(here)
library(flextable)
library(broom)

source(here('R/funcs.R'))

load(file = here('data/mods.RData'))
load(file = here('data/hydat.RData'))

# GAM fits ---------------------------------------------------------------

totab <- mods |>
  select(bay_segment, mod) |>
  mutate(
    summ = purrr::map(mod, function(x) {
      summmod <- summary(x)
      n <- summmod$n
      GCV <- round(summmod$sp.criterion[[1]], 0)
      devexpl <- round(summmod$dev.expl, 2)

      smths <- x |>
        tidy() |>
        rowwise() |>
        mutate(
          p.value = p_txt(p.value, addp = F)
        ) |>
        ungroup() |>
        mutate_if(is.numeric, round, 2) |>
        mutate_if(is.numeric, as.character) |>
        rename(
          'Smoother' = term,
          'Ref.df' = ref.df,
          'F' = statistic,
          'p' = p.value
        )

      out <- bind_cols(n = n, GCV = GCV, devexpl = devexpl, smths) |>
        mutate(
          n = ifelse(duplicated(n), '', n),
          GCV = ifelse(duplicated(GCV), '', GCV),
          devexpl = ifelse(duplicated(devexpl), '', devexpl)
        ) |>
        rename(
          `Num. Obs.` = n,
          `Dev. Expl.` = devexpl
        )

      return(out)
    })
  ) |>
  select(-mod) |>
  unnest('summ')

gamtab <- totab |>
  as_grouped_data(groups = 'bay_segment') |>
  flextable() |>
  set_header_labels(bay_segment = 'Bay Segment') |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  autofit()

save(gamtab, file = here('tabs/gamtab.RData'))


# GAM fits by time slides ------------------------------------------------

datprp <- mods |>
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |>
  mutate(
    yr = lubridate::year(date),
    yrcat = case_when(
      yr < 2004 ~ '1985-2003',
      yr >= 2004 ~ '2004-2024'
    ),
    qrt = factor(
      lubridate::quarter(date),
      levels = 1:4,
      labels = c('JFM', 'AMJ', 'JAS', 'OND')
    )
  )

totabqrt <- datprp |>
  summarise(
    cr = cor(chla, btfit, use = 'complete.obs'), #summary(lm(btfit ~ chla))$r.squared,
    .by = c(bay_segment, yrcat, qrt)
  )
totabann <- datprp |>
  summarise(
    cr = cor(chla, btfit, use = 'complete.obs'), #summary(lm(btfit ~ chla))$r.squared,
    .by = c(bay_segment, yrcat)
  ) |>
  mutate(
    qrt = 'Annual'
  )

totab <- bind_rows(totabqrt, totabann) |>
  mutate(
    qrt = factor(
      qrt,
      levels = c('Annual', 'JFM', 'AMJ', 'JAS', 'OND')
    ),
    cr = round(cr, 2)
  ) |>
  pivot_wider(
    names_from = qrt,
    values_from = cr
  ) |>
  select(bay_segment, yrcat, Annual, JFM, AMJ, JAS, OND)

# color function
colfun <- function(x) {
  x[is.na(x)] <- 0
  pal <- RColorBrewer::brewer.pal(9, 'Greys')
  colorRampPalette(pal)(100)[as.numeric(cut(x, breaks = 100))]
}

# font size function
fntfun <- function(x) {
  x[is.na(x)] <- min(x, na.rm = T)
  sizes <- seq(7, 11, length.out = 50)
  sizes[as.numeric(cut(x, breaks = 50))]
}

# Pre-calculate font sizes for each column
cols <- c('Annual', 'JFM', 'AMJ', 'JAS', 'OND')
ft_data <- totab |> as_grouped_data(groups = 'bay_segment')

gamcrtab <- ft_data |>
  flextable() |>
  set_header_labels(
    bay_segment = 'Bay Segment',
    yrcat = 'Year group'
  ) |>
  color(j = ~ Annual + JFM + AMJ + JAS + OND, color = colfun)

# Apply font sizes for each cell
font_sizes <- fntfun(unlist(ft_data[, cols])) |>
  matrix(nrow = nrow(ft_data), ncol = length(cols)) |>
  as.data.frame()
names(font_sizes) <- cols
for (col in cols) {
  font_col <- font_sizes[[col]]
  for (i in seq_along(font_col)) {
    gamcrtab <- fontsize(
      gamcrtab,
      i = i,
      j = col,
      size = font_col[i],
      part = 'body'
    )
  }
}

gamcrtab <- gamcrtab |>
  padding(padding = 0, part = 'all') |>
  bold(j = cols, part = 'body') |>
  align(align = 'center', part = 'all', j = cols) |>
  font(part = 'all', fontname = 'Times New Roman') |>
  autofit()

save(gamcrtab, file = here('tabs/gamcrtab.RData'))

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

modsumm <- qrtprds |>
  mutate(
    anngrp = case_when(
      yr < 2004 ~ '1985-2003',
      yr >= 2004 ~ '2004-2024'
    )
  ) |>
  group_by(bay_segment, qrt, anngrp) |>
  nest() |>
  mutate(
    mod = purrr::map(
      data,
      ~ lm(prdresid ~ I(hyqrt - mean(hyqrt, na.rm = T)), data = .x)
    ),
    summ = purrr::map(mod, function(x) {
      out <- tibble(
        n = NA,
        slo = NA,
        slose = NA,
        pval = NA,
        rsq = NA
      )

      n <- nrow(x$model)
      coef <- summary(x)$coefficients

      # output
      out$n <- n
      out$slo <- coef[2, 1]
      out$slose <- coef[2, 2]
      out$pval <- coef[2, 4]
      out$rsq <- summary(x)$r.squared

      return(out)
    })
  )

totab <- modsumm |>
  select(bay_segment, qrt, anngrp, summ) |>
  unnest(summ) |>
  ungroup() |>
  mutate(
    n = as.character(n),
    slo = as.character(round(slo, 3)),
    slose = paste0('(', round(slose, 3), ')'),
    pval = p_ast(pval),
    rsq = as.character(round(rsq, 2))
  ) |>
  unite('slo', slo, slose, sep = ' ') |>
  unite('slo', slo, pval, sep = '') |>
  pivot_longer(
    cols = c(n, slo, rsq),
    names_to = 'Metric',
    values_to = 'Value'
  ) |>
  unite('Dataset', anngrp, Metric, sep = ' ') |>
  pivot_wider(names_from = Dataset, values_from = Value) |>
  select(-`1985-2003 n`, -`2004-2024 n`)

hydnrmtab <- totab |>
  as_grouped_data(groups = 'bay_segment') |>
  flextable() |>
  add_header_row(
    values = c('', '', rep('1985 - 2003', 2), rep('2004 - 2024', 2))
  ) |>
  merge_at(i = 1, j = c(3:4), part = 'header') |>
  merge_at(i = 1, j = c(5:6), part = 'header') |>
  set_header_labels(
    i = 2,
    values = c(
      'Bay Segment',
      'Quarter',
      'Slope (SE)',
      'R-Squared',
      'Slope (SE)',
      'R-Squared'
    )
  ) |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  autofit()

save(hydnrmtab, file = here('tabs/hydnrmtab.RData'))
