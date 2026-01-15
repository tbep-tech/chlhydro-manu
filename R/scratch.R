library(tidyverse)


tst <- newdat |>
  group_nest(bay_segment) |>
  rename(tst = data) |>
  left_join(trgs, by = c('bay_segment'), relationship = 'one-to-one')

p <- mods |>
  left_join(tst, by = c('bay_segment'), relationship = 'one-to-one') |>
  mutate(
    p = pmap(list(mod, tst), function(mod, tst) {
      simulate(mod, data = tst, nsim = 100) |>
        bind_cols(tst) |>
        pivot_longer(
          cols = starts_with('sim'),
          names_to = 'sim',
          values_to = 'chla_sim'
        ) |>
        mutate(
          sim = as.integer(str_remove(sim, 'sim_'))
        ) |>
        arrange(sim, date)
    }),
    pyr = map(p, function(p) {
      p |>
        summarise(
          chla_sim = mean(chla_sim, na.rm = T),
          .by = c(yrs, yr, ldfac, sim)
        )
    }),
    pyr2 = pmap(list(pyr, thresh), function(pyr, thresh) {
      pyr |>
        mutate(
          exceeds = chla_sim > thresh
        ) |>
        summarise(
          percexceeds = mean(exceeds, na.rm = T) * 100,
          .by = c(yrs, sim, ldfac)
        )
    }),
    pyr3 = map(pyr2, function(pyr2) {
      pyr2 |>
        summarise(
          avexceeds = mean(percexceeds, na.rm = T),
          sdexceeds = sd(percexceeds, na.rm = T),
          .by = c(yrs, ldfac)
        )
    })
  )
