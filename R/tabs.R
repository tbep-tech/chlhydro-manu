library(tidyverse)
library(mgcv)
library(here)
library(flextable)
library(broom)

source(here('R/funcs.R'))

load(file = here('data/mods.RData'))

# GAM fits ---------------------------------------------------------------

totab <- mods |>
  select(bay_segment, mod) |>
  mutate(
    summ = purrr::map(mod, function(x) {
      summmod <- summary(x)
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

      out <- bind_cols(GCV = GCV, devexpl = devexpl, smths) |>
        mutate(
          GCV = ifelse(duplicated(GCV), '', GCV),
          devexpl = ifelse(duplicated(devexpl), '', devexpl)
        ) |>
        rename(
          `Dev. Expl.` = devexpl
        )

      return(out)
    })
  ) |>
  select(-mod) |>
  unnest('summ') |>
  mutate(
    bay_segment = ifelse(duplicated(bay_segment), '', as.character(bay_segment))
  )

gamtab <- totab |>
  flextable() |>
  set_header_labels(bay_segment = 'Bay Segment') |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  autofit() |>
  bold(~ as.numeric(gsub('<', '', p)) < 0.05, j = 8) |>
  fontsize(size = 9, part = "body") |>
  hline(i = 10) |>
  hline(i = 20) |>
  hline(i = 30)

save(gamtab, file = here('tabs/gamtab.RData'))

# GAM fits by time slices ------------------------------------------------

datprp <- mods |>
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |>
  mutate(
    yr = lubridate::year(date),
    yrcat = case_when(
      yr < 2000 ~ '1985 - 1999',
      yr >= 2000 & yr < 2015 ~ '2000 - 2014',
      yr >= 2015 ~ '2015 - 2024'
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
  select(bay_segment, yrcat, Annual, JFM, AMJ, JAS, OND) |>
  mutate(
    bay_segment = ifelse(duplicated(bay_segment), '', as.character(bay_segment))
  )

# color function
colfun <- function(x) {
  x[is.na(x)] <- min(x, na.rm = T)
  pal <- RColorBrewer::brewer.pal(9, 'Greys')[5:9]
  colorRampPalette(pal)(100)[as.numeric(cut(x, breaks = 100))]
}

# font size function
fntfun <- function(x) {
  x[is.na(x)] <- min(x, na.rm = T)
  sizes <- seq(7, 11, length.out = 50)
  sizes[as.numeric(cut(x, breaks = 50))]
}

# Pre-calculate colors across all columns
cols <- c('Annual', 'JFM', 'AMJ', 'JAS', 'OND')

# Calculate colors for all values at once
cell_colors <- colfun(unlist(totab[, cols])) |>
  matrix(nrow = nrow(totab), ncol = length(cols)) |>
  as.data.frame()
names(cell_colors) <- cols

gamcrtab <- totab |>
  flextable() |>
  set_header_labels(
    bay_segment = 'Bay Segment',
    yrcat = 'Year Group'
  )

# Apply colors for each cell
for (col in cols) {
  color_col <- cell_colors[[col]]
  for (i in seq_along(color_col)) {
    gamcrtab <- color(
      gamcrtab,
      i = i,
      j = col,
      color = color_col[i],
      part = 'body'
    )
  }
}

# Apply font sizes for each cell
font_sizes <- fntfun(unlist(totab[, cols])) |>
  matrix(nrow = nrow(totab), ncol = length(cols)) |>
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
  hline(i = 3) |>
  hline(i = 6) |>
  hline(i = 9) |>
  autofit()

save(gamcrtab, file = here('tabs/gamcrtab.RData'))

# pred - norm by time slices ---------------------------------------------

datprp <- mods |>
  select(bay_segment, prds) |>
  mutate(prds = purrr::map(prds, as_tibble)) |>
  unnest(prds) |>
  mutate(
    yr = lubridate::year(date),
    yrcat = case_when(
      yr < 2000 ~ '1985 - 1999',
      yr >= 2000 & yr < 2015 ~ '2000 - 2014',
      yr >= 2015 ~ '2015 - 2024'
    ),
    qrt = factor(
      lubridate::quarter(date),
      levels = 1:4,
      labels = c('JFM', 'AMJ', 'JAS', 'OND')
    )
  )

# mean absolute difference
totabqrt <- datprp |>
  summarise(
    salmad = mean(btfit - btnorm, na.rm = T),
    ldmad = mean(btnorm - btnormmd, na.rm = T),
    .by = c(bay_segment, yrcat, qrt)
  )
totabann <- datprp |>
  summarise(
    salmad = mean(btfit - btnorm, na.rm = T),
    ldmad = mean(btnorm - btnormmd, na.rm = T),
    .by = c(bay_segment, yrcat)
  ) |>
  mutate(
    qrt = 'Annual'
  )

levs <- c('Annual', 'JFM', 'AMJ', 'JAS', 'OND')
totab <- bind_rows(totabqrt, totabann) |>
  mutate(
    qrt = factor(
      qrt,
      levels = levs
    ),
    salmad = round(salmad, 2),
    ldmad = round(ldmad, 2)
  ) |>
  pivot_longer(
    names_to = 'var',
    values_to = 'mad',
    cols = c('salmad', 'ldmad')
  ) |>
  mutate(
    var = gsub('mad$', '', var),
    qrt = as.character(qrt)
  ) |>
  unite('var', var, qrt) |>
  mutate(
    var = factor(var, levels = c(paste0('sal_', levs), paste0('ld_', levs)))
  ) |>
  pivot_wider(
    names_from = var,
    values_from = mad,
    names_sort = TRUE
  ) |>
  mutate(
    bay_segment = ifelse(duplicated(bay_segment), '', as.character(bay_segment))
  )

# Pre-calculate font sizes for each column
cols <- c(paste0('sal_', levs), paste0('ld_', levs))

# color function
colfun <- leaflet::colorNumeric(
  palette = c('blue', 'blue', 'dodgerblue3', 'grey', 'tomato1', 'red', 'red'),
  domain = c(
    -max(abs(totab[, cols]), na.rm = T),
    max(abs(totab[, cols]), na.rm = T)
  )
)

# font size function
fntfun <- function(x) {
  x[is.na(x)] <- min(x, na.rm = T)
  sizes <- seq(7, 11, length.out = 50)
  sizes[as.numeric(cut(x, breaks = 50))]
}

prdnrmtab <- totab |>
  flextable()

# Calculate colors for all values at once
cell_colors <- colfun(unlist(totab[, cols])) |>
  matrix(nrow = nrow(totab), ncol = length(cols)) |>
  as.data.frame()
names(cell_colors) <- cols

# Apply colors for each cell
for (col in cols) {
  color_col <- cell_colors[[col]]
  for (i in seq_along(color_col)) {
    prdnrmtab <- color(
      prdnrmtab,
      i = i,
      j = col,
      color = color_col[i],
      part = 'body'
    )
  }
}

# Apply font sizes for each cell
font_sizes <- fntfun(abs(unlist(totab[, cols]))) |>
  matrix(nrow = nrow(totab), ncol = length(cols)) |>
  as.data.frame()
names(font_sizes) <- cols
for (col in cols) {
  font_col <- font_sizes[[col]]
  for (i in seq_along(font_col)) {
    prdnrmtab <- fontsize(
      prdnrmtab,
      i = i,
      j = col,
      size = font_col[i],
      part = 'body'
    )
  }
}

prdnrmtab <- prdnrmtab |>
  padding(padding = 0, part = 'all') |>
  bold(j = cols, part = 'body') |>
  align(align = 'center', part = 'all', j = cols) |>
  align(align = 'left', part = 'header', j = cols, i = 1) |>
  font(part = 'all', fontname = 'Times New Roman') |>
  add_header_row(
    values = c('', '', rep('Salinity', 5), rep('Load', 5))
  ) |>
  merge_at(i = 1, j = c(3:7), part = 'header') |>
  merge_at(i = 1, j = c(8:12), part = 'header') |>
  set_header_labels(
    i = 2,
    values = c(
      'Bay Segment',
      'Year Group',
      rep(c('Annual', 'JFM', 'AMJ', 'JJA', 'OND'), 2)
    )
  ) |>
  hline(i = 3) |>
  hline(i = 6) |>
  hline(i = 9) |>
  autofit()

save(prdnrmtab, file = here('tabs/prdnrmtab.RData'))

# simulation summaries ---------------------------------------------------

load(file = here('data/simprddat1524.RData'))

salests <- simprddat1524 |>
  select(bay_segment, tst) |>
  unnest(tst) |>
  filter(yrs %in% c(25, 50)) |>
  select(bay_segment, Time = yrs, salforeyr) |>
  distinct() |>
  mutate(
    Time = paste('Year', Time)
  )

rho <- 0.7
totab <- simprddat1524 |>
  select(bay_segment, exceedssum) |>
  unnest(exceedssum) |>
  filter(yrs %in% c(1, 25, 50)) |>
  mutate(yrs = paste('yr', yrs)) |>
  pivot_longer(names_to = 'var', values_to = 'val', avexceeds:sdexceeds) |>
  pivot_wider(names_from = 'yrs', values_from = 'val') |>
  mutate(
    `yr 25diff` = case_when(
      var == 'avexceeds' ~ `yr 25` - `yr 1`,
      var == 'sdexceeds' ~ sqrt(
        `yr 25`^2 + `yr 1`^2 - 2 * rho * `yr 25` * `yr 1`
      )
    ),
    `yr 50diff` = case_when(
      var == 'avexceeds' ~ `yr 50` - `yr 1`,
      var == 'sdexceeds' ~ sqrt(
        `yr 50`^2 + `yr 1`^2 - 2 * rho * `yr 50` * `yr 1`
      )
    )
  ) |>
  pivot_wider(
    names_from = var,
    values_from = c(`yr 1`, `yr 25`, `yr 50`, `yr 25diff`, `yr 50diff`)
  ) |>
  select(-contains('diff_sd')) |>
  mutate_if(is.numeric, round, 1) |>
  mutate(across(contains('sdexceeds'), ~ paste0('\n(', .x))) |>
  mutate(across(
    contains('diff'),
    ~ ifelse(
      sign(.x) == 1,
      paste0('\u2191', .x),
      ifelse(
        sign(.x) == 0,
        paste0('\u2248', .x),
        paste0('\u2193', .x)
      )
    )
  )) |>
  mutate(across(contains('diff'), ~ paste0(', ', .x, ')'))) |>
  unite(
    'yr 25_sdexceeds',
    `yr 25_sdexceeds`,
    `yr 25diff_avexceeds`,
    sep = ''
  ) |>
  unite(
    'yr 50_sdexceeds',
    `yr 50_sdexceeds`,
    `yr 50diff_avexceeds`,
    sep = ''
  ) |>
  mutate(
    `yr 1_sdexceeds` = paste0(`yr 1_sdexceeds`, ')')
  ) |>
  unite('Present', `yr 1_avexceeds`, `yr 1_sdexceeds`, sep = ' ') |>
  unite('Year 25', `yr 25_avexceeds`, `yr 25_sdexceeds`, sep = ' ') |>
  unite('Year 50', `yr 50_avexceeds`, `yr 50_sdexceeds`, sep = ' ') |>
  mutate(
    ldfac = factor(
      ldfac,
      levels = c(
        'Actual Load',
        'TMDL Load Factor: 0.5',
        'TMDL Load Factor: 1',
        'TMDL Load Factor: 2'
      )
    )
  ) |>
  pivot_longer(names_to = 'Time', values_to = 'val', Present:`Year 50`) |>
  pivot_wider(names_from = 'ldfac', values_from = 'val', names_sort = TRUE) |>
  left_join(salests, by = c('bay_segment', 'Time')) |>
  mutate(
    salforeyr = ifelse(is.na(salforeyr), '\u2014', round(salforeyr, 2))
  ) |>
  select(
    bay_segment,
    Time,
    `Salinity Change (ppth)` = salforeyr,
    everything()
  ) |>
  mutate(bay_segment = ifelse(duplicated(bay_segment), '', bay_segment))

liktab <- totab |>
  flextable() |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  add_header_row(
    values = c(rep('', 4), rep('TMDL Load Factor', 3))
  ) |>
  merge_at(i = 1, j = c(5:7), part = 'header') |>
  set_header_labels(
    i = 2,
    values = c(
      'Bay\nSegment',
      'Time',
      'Salinity\nChange (ppth)',
      'Actual Load',
      '0.5',
      '1',
      '2'
    )
  ) |>
  valign(valign = "top", part = "all") |>
  hline(i = 3) |>
  hline(i = 6) |>
  hline(i = 9) |>
  autofit()

save(liktab, file = here('tabs/liktab.RData'))

# supp tab lag comp ------------------------------------------------------

# water quality
load(file = here('data/wqdat.RData'))
load(file = here('data/lddat.RData'))

lddatlag <- lddat |>
  summarise(
    tn_load = sum(tn_load),
    .by = c(bay_segment, date)
  ) |>
  mutate(
    lag1 = dplyr::lag(tn_load, n = 1),
    lag2 = dplyr::lag(tn_load, n = 2),
    lag3 = dplyr::lag(tn_load, n = 3),
    .by = c(bay_segment)
  ) |>
  fill(
    lag1,
    lag2,
    lag3,
    .by = c(bay_segment),
    .direction = 'up'
  ) |>
  mutate(
    tn_load0 = tn_load,
    tn_load1 = tn_load + lag1,
    tn_load2 = tn_load + lag1 + lag2,
    tn_load3 = tn_load + lag1 + lag2 + lag3
  ) |>
  select(-tn_load, -lag1, -lag2, -lag3)

wqdat <- wqdat |>
  select(-tn_load) |>
  left_join(lddatlag, by = c('bay_segment', 'date'))

mods <- wqdat |>
  filter(dec_time >= 1985) |>
  pivot_longer(
    cols = starts_with('tn_load'),
    names_to = 'tn_load_lag',
    values_to = 'tn_load'
  ) |>
  group_nest(bay_segment, tn_load_lag) |>
  mutate(
    mod = purrr::map(
      data,
      ~ gam(
        chla ~ s(dec_time, k = 40, bs = 'tp') +
          s(doy, k = 10, bs = 'cc') +
          s(sal, k = 10, bs = 'tp') +
          s(tn_load, k = 10, bs = 'tp') +
          ti(dec_time, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(dec_time, sal, k = c(5, 5), bs = c('tp', 'tp')) +
          ti(sal, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(dec_time, tn_load, k = c(5, 5), bs = c('tp', 'tp')) +
          ti(tn_load, doy, k = c(5, 5), bs = c('tp', 'cc')) +
          ti(tn_load, sal, k = c(5, 5), bs = c('tp', 'tp')),
        data = .x,
        knots = list(doy = c(0, 366)),
        family = Gamma(link = 'log'),
        na.action = na.exclude,
        method = 'REML'
      )
    ),
    prds = purrr::map2(data, mod, pred_fun),
    annsum = purrr::map(
      prds,
      ~ data.frame(.x) |>
        mutate(
          yr = lubridate::year(date)
        ) |>
        summarise(
          chla = mean(chla, na.rm = T),
          btfit = mean(btfit, na.rm = T),
          btnorm = mean(btnorm, na.rm = T),
          btfitmd = mean(btfitmd, na.rm = T),
          btnormmd = mean(btnormmd, na.rm = T),
          .by = c(yr)
        )
    )
  )

totab <- mods |>
  select(bay_segment, tn_load_lag, mod) |>
  mutate(
    summ = purrr::map(mod, function(x) {
      summmod <- summary(x)
      GCV <- round(summmod$sp.criterion[[1]], 0)
      devexpl <- round(summmod$dev.expl, 2)

      out <- bind_cols(GCV = GCV, devexpl = devexpl) |>
        rename(
          `Dev. Expl.` = devexpl
        )

      return(out)
    })
  ) |>
  select(-mod) |>
  unnest('summ') |>
  mutate(
    bay_segment = ifelse(
      duplicated(bay_segment),
      '',
      as.character(bay_segment)
    ),
    tn_load_lag = case_when(
      tn_load_lag == 'tn_load0' ~ 'None',
      tn_load_lag == 'tn_load1' ~ '1',
      tn_load_lag == 'tn_load2' ~ '2',
      tn_load_lag == 'tn_load3' ~ '3'
    )
  )

suppgamtab <- totab |>
  flextable() |>
  set_header_labels(bay_segment = 'Bay Segment', tn_load_lag = 'TN Load Lag') |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  colformat_double(
    j = 3,
    big.mark = '',
    digits = 0
  ) |>
  autofit() |>
  hline(i = 4) |>
  hline(i = 8) |>
  hline(i = 12)

save(suppgamtab, file = here('tabs/suppgamtab.RData'))
