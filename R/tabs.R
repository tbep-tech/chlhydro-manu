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
  bold(~ as.numeric(gsub('<', '', p)) < 0.05, j = 9) |>
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
  select(bay_segment, yrcat, Annual, JFM, AMJ, JAS, OND)

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

# Pre-calculate font sizes for each column
cols <- c('Annual', 'JFM', 'AMJ', 'JAS', 'OND')
ft_data <- totab

gamcrtab <- ft_data |>
  flextable() |>
  set_header_labels(
    bay_segment = 'Bay Segment',
    yrcat = 'Year Group'
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
    salmad = mean(abs(btfit - btnorm), na.rm = T),
    ldmad = mean(abs(btnorm - btnormmd), na.rm = T),
    .by = c(bay_segment, yrcat, qrt)
  )
totabann <- datprp |>
  summarise(
    salmad = mean(abs(btfit - btnorm), na.rm = T),
    ldmad = mean(abs(btnorm - btnormmd), na.rm = T),
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

# Pre-calculate font sizes for each column
cols <- c(paste0('sal_', levs), paste0('ld_', levs))
ft_data <- totab |>
  mutate(
    bay_segment = ifelse(duplicated(bay_segment), '', as.character(bay_segment))
  )

prdnrmtab <- ft_data |>
  flextable() |>
  color(j = cols, color = colfun)

# Apply font sizes for each cell
font_sizes <- fntfun(unlist(ft_data[, cols])) |>
  matrix(nrow = nrow(ft_data), ncol = length(cols)) |>
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

load(file = here('data/simprddat1721.RData'))
