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
      yr < 1999 ~ '1985 - 1998',
      yr >= 1999 & yr < 2012 ~ '1999 - 2011',
      yr >= 2012 ~ '2012 - 2024'
    ),
    seas = factor(
      ifelse(lubridate::month(date) %in% 6:9, 'Wet', 'Dry'),
      levels = c('Wet', 'Dry')
    )
  )

totabseas <- datprp |>
  summarise(
    cr = cor(chla, btfit, use = 'complete.obs'), #summary(lm(btfit ~ chla))$r.squared,
    .by = c(bay_segment, yrcat, seas)
  )
totabann <- datprp |>
  summarise(
    cr = cor(chla, btfit, use = 'complete.obs'), #summary(lm(btfit ~ chla))$r.squared,
    .by = c(bay_segment, yrcat)
  ) |>
  mutate(
    seas = 'Annual'
  )

totab <- bind_rows(totabseas, totabann) |>
  mutate(
    seas = factor(
      seas,
      levels = c('Annual', 'Wet', 'Dry')
    ),
    cr = round(cr, 2)
  ) |>
  pivot_wider(
    names_from = seas,
    values_from = cr
  ) |>
  select(bay_segment, yrcat, Annual, Wet, Dry) |>
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
cols <- c('Annual', 'Wet', 'Dry')

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
      yr < 1999 ~ '1985 - 1998',
      yr >= 1999 & yr < 2012 ~ '1999 - 2011',
      yr >= 2012 ~ '2012 - 2024'
    ),
    seas = factor(
      ifelse(lubridate::month(date) %in% 6:9, 'Wet', 'Dry'),
      levels = c('Wet', 'Dry')
    )
  )

# mean absolute difference
totabseas <- datprp |>
  summarise(
    salmad = mean(btfit - btnorm, na.rm = T),
    ldmad = mean(btnorm - btnormmd, na.rm = T),
    .by = c(bay_segment, yrcat, seas)
  )
totabann <- datprp |>
  summarise(
    salmad = mean(btfit - btnorm, na.rm = T),
    ldmad = mean(btnorm - btnormmd, na.rm = T),
    .by = c(bay_segment, yrcat)
  ) |>
  mutate(
    seas = 'Annual'
  )

levs <- c('Annual', 'Wet', 'Dry')
totab <- bind_rows(totabseas, totabann) |>
  mutate(
    seas = factor(
      seas,
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
    seas = as.character(seas)
  ) |>
  unite('var', var, seas) |>
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
    values = c('', '', rep('Salinity', 3), rep('Load', 3))
  ) |>
  merge_at(i = 1, j = c(3:5), part = 'header') |>
  merge_at(i = 1, j = c(6:8), part = 'header') |>
  set_header_labels(
    i = 2,
    values = c(
      'Bay Segment',
      'Year Group',
      rep(c('Annual', 'Wet', 'Dry'), 2)
    )
  ) |>
  hline(i = 3) |>
  hline(i = 6) |>
  hline(i = 9) |>
  autofit()

save(prdnrmtab, file = here('tabs/prdnrmtab.RData'))

# simulation summaries ---------------------------------------------------

load(file = here('data/simprddat.RData'))

salests <- simprddat |>
  select(bay_segment, tst) |>
  unnest(tst) |>
  filter(yrs %in% c(25, 50)) |>
  select(bay_segment, Time = yrs, salforeyr) |>
  distinct() |>
  mutate(
    Time = paste('Year', Time)
  )

rho <- 0.7
totab <- simprddat |>
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

# supp gamma vs gauss log models -----------------------------------------

load(file = here('data/wqdat.RData'))
load(file = here('data/mods.RData'))

modsgamma <- mods |>
  select(-prds, -annsum)

modsgauss <- wqdat |>
  filter(dec_time >= 1985) |>
  group_nest(bay_segment) |>
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
        family = gaussian(link = 'log'),
        na.action = na.exclude,
        method = 'REML'
      )
    )
  )

AIC(modsgamma$mod[[1]], modsgauss$mod[[1]])
# AIC(modsgamma$mod[[2]], modsgauss$mod[[2]])
# AIC(modsgamma$mod[[3]], modsgauss$mod[[3]])
# AIC(modsgamma$mod[[4]], modsgauss$mod[[4]])

aiccomp <- modsgamma |>
  select(bay_segment, mod_gamma = mod) |>
  left_join(
    modsgauss |> select(bay_segment, mod_gauss = mod),
    by = 'bay_segment'
  ) |>
  mutate(
    aic = purrr::map2(mod_gamma, mod_gauss, ~ AIC(.x, .y)),
    gamma_df = purrr::map_dbl(aic, ~ round(.x[1, 'df'], 1)),
    gamma_aic = purrr::map_dbl(aic, ~ round(.x[1, 'AIC'], 1)),
    gauss_df = purrr::map_dbl(aic, ~ round(.x[2, 'df'], 1)),
    gauss_aic = purrr::map_dbl(aic, ~ round(.x[2, 'AIC'], 1)),
    delta_aic = round(gauss_aic - gamma_aic, 1)
  ) |>
  select(bay_segment, gamma_df, gamma_aic, gauss_df, gauss_aic, delta_aic)

suppaictab <- aiccomp |>
  flextable() |>
  add_header_row(
    values = c('', 'Gamma', 'Gaussian (log)', ''),
    colwidths = c(1, 2, 2, 1)
  ) |>
  set_header_labels(
    i = 2,
    values = c(
      'Bay Segment',
      'df',
      'AIC',
      'df',
      'AIC',
      'ΔAIC'
    )
  ) |>
  merge_at(i = 1, j = 2:3, part = 'header') |>
  merge_at(i = 1, j = 4:5, part = 'header') |>
  align(align = 'center', part = 'header') |>
  align(align = 'left', part = 'header', i = 1) |>
  align(align = 'left', part = 'all', j = 2:6) |>
  bold(i = which(aiccomp$gamma_aic < aiccomp$gauss_aic), j = 'gamma_aic') |>
  bold(i = which(aiccomp$gauss_aic < aiccomp$gamma_aic), j = 'gauss_aic') |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  fontsize(size = 9, part = 'body') |>
  add_footer_lines(
    'ΔAIC = Gaussian (log) AIC - Gamma AIC'
  ) |>
  font(part = 'footer', fontname = 'Times New Roman') |>
  fontsize(size = 8, part = 'footer') |>
  autofit()

save(suppaictab, file = here('tabs/suppaictab.RData'))

# supp tab lag comp ------------------------------------------------------

# load model data
load(file = here('data/modslag.RData'))
load(file = here('data/modslagk.RData'))
load(file = here('data/modsreslag.RData'))
load(file = here('data/modsreslagk.RData'))

totab <- list(
  `No seasonal correction, flexible K` = modslag,
  `No seasonal correction, fixed K` = modslagk,
  `Seasonal correction, flexible K` = modsreslag,
  `Seasonal correction, fixed K` = modsreslagk
) |>
  tibble::enframe() |>
  unnest(value) |>
  select(name, bay_segment, tn_load_lag, mod) |>
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
  pivot_longer(
    names_to = 'metric',
    values_to = 'value',
    c(GCV, `Dev. Expl.`)
  ) |>
  pivot_wider(names_from = name, values_from = value) |>
  select(metric, everything()) |>
  arrange(desc(metric)) |>
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
    ),
    .by = c(metric)
  )

suppgamtab <- totab |>
  as_grouped_data('metric') |>
  flextable() |>
  set_header_labels(
    metric = 'Summary',
    bay_segment = 'Bay Segment',
    tn_load_lag = 'TN Load Lag',
    `No seasonal correction, flexible K` = 'Flexible K',
    `No seasonal correction, fixed K` = 'Fixed K',
    `Seasonal correction, flexible K` = 'Flexible K',
    `Seasonal correction, fixed K` = 'Fixed K'
  ) |>
  add_header_row(
    values = c('', '', '', 'No Seasonal Correction', 'Seasonal Correction'),
    colwidths = c(1, 1, 1, 2, 2)
  ) |>
  padding(padding = 0, part = 'all') |>
  font(part = 'all', fontname = 'Times New Roman') |>
  align(align = 'center', part = 'all', j = c(4:7)) |>
  colformat_double(
    j = 4:7,
    i = 2:17,
    big.mark = '',
    digits = 0
  ) |>
  autofit() |>
  hline(i = c(5, 9, 13, 17, 22, 26, 30))

save(suppgamtab, file = here('tabs/suppgamtab.RData'))

# supp concurvity --------------------------------------------------------

# load model data
load(file = here('data/modslag.RData'))
load(file = here('data/modslagk.RData'))
load(file = here('data/modsreslag.RData'))
load(file = here('data/modsreslagk.RData'))

totab <- list(
  `No seasonal correction, flexible K` = modslag,
  `No seasonal correction, fixed K` = modslagk,
  `Seasonal correction, flexible K` = modsreslag,
  `Seasonal correction, fixed K` = modsreslagk
) |>
  tibble::enframe() |>
  unnest(value) |>
  select(name, bay_segment, tn_load_lag, mod) |>
  mutate(
    summ = purrr::map(mod, function(x) {
      nms <- cols <- c(
        's(dec_time)',
        's(doy)',
        's(sal)',
        's(tn_load)',
        'ti(dec_time,doy)',
        'ti(dec_time,sal)',
        'ti(sal,doy)',
        'ti(dec_time,tn_load)',
        'ti(tn_load,doy)',
        'ti(tn_load,sal)'
      )
      mat <- concurvity(x)[3, nms] |>
        data.frame() |>
        t() |>
        data.frame()
      names(mat) <- nms
      row.names(mat) <- 1
      return(mat)
    })
  ) |>
  select(-mod) |>
  unnest('summ') |>
  mutate(
    tn_load_lag = case_when(
      tn_load_lag == 'tn_load0' ~ 'None',
      tn_load_lag == 'tn_load1' ~ '1',
      tn_load_lag == 'tn_load2' ~ '2',
      tn_load_lag == 'tn_load3' ~ '3'
    )
  )

tabfun <- function(x, bayseg) {
  totab <- x |>
    filter(bay_segment == bayseg) |>
    select(-bay_segment) |>
    separate(name, into = c('season', 'k'), sep = ', ') |>
    mutate(
      k = case_when(
        k == 'flexible K' ~ 'Flexible',
        k == 'fixed K' ~ 'Fixed',
      ),
      season = case_when(
        season == 'No seasonal correction' ~ 'No',
        season == 'Seasonal correction' ~ 'Yes'
      )
    ) |>
    mutate(
      k = ifelse(duplicated(k), '', k),
      .by = season
    ) |>
    mutate(
      season = ifelse(duplicated(season), '', season)
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
    sizes <- seq(7, 9, length.out = 50)
    sizes[as.numeric(cut(x, breaks = 50))]
  }

  # Pre-calculate colors across all columns
  cols <- c(
    's(dec_time)',
    's(doy)',
    's(sal)',
    's(tn_load)',
    'ti(dec_time,doy)',
    'ti(dec_time,sal)',
    'ti(sal,doy)',
    'ti(dec_time,tn_load)',
    'ti(tn_load,doy)',
    'ti(tn_load,sal)'
  )

  # Calculate colors for all values at once
  cell_colors <- colfun(unlist(totab[, cols])) |>
    matrix(nrow = nrow(totab), ncol = length(cols)) |>
    as.data.frame()
  names(cell_colors) <- cols

  tab <- totab |>
    flextable() |>
    set_header_labels(
      k = 'K',
      season = 'Seasonal\nCorrection',
      tn_load_lag = 'TN Load Lag'
    )

  # Apply colors for each cell
  for (col in cols) {
    color_col <- cell_colors[[col]]
    for (i in seq_along(color_col)) {
      tab <- color(
        tab,
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
      tab <- fontsize(
        tab,
        i = i,
        j = col,
        size = font_col[i],
        part = 'body'
      )
    }
  }

  out <- tab |>
    padding(padding = 0, part = 'all') |>
    bold(j = cols, part = 'body') |>
    align(align = 'center', part = 'all', j = cols) |>
    font(part = 'all', fontname = 'Times New Roman') |>
    fontsize(size = 9, part = 'header') |>
    fontsize(size = 9, part = 'body', j = 1:3) |>
    colformat_double(
      digits = 2
    ) |>
    hline(i = 4, j = 2:13) |>
    hline(i = 8) |>
    hline(i = 12, j = 2:13) |>
    width(j = 4:13, width = 0.5) |>
    width(j = 1:3, width = 0.5)

  return(out)
}

suppconcurvotbtab <- tabfun(totab, 'OTB')
suppconcurvhbtab <- tabfun(totab, 'HB')
suppconcurvmtbtab <- tabfun(totab, 'MTB')
suppconcurvltbtab <- tabfun(totab, 'LTB')

save(suppconcurvotbtab, file = here('tabs/suppconcurvotbtab.RData'))
save(suppconcurvhbtab, file = here('tabs/suppconcurvhbtab.RData'))
save(suppconcurvmtbtab, file = here('tabs/suppconcurvmtbtab.RData'))
save(suppconcurvltbtab, file = here('tabs/suppconcurvltbtab.RData'))
