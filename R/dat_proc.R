library(tbeptools)
library(tidyverse)
library(mgcv)
library(here)

source(here('R/funcs.R'))

# data prep --------------------------------------------------------------

# water quality
wqdat <- epcdata |>
  filter(yr > 1975 & yr <= 2024) |>
  select(
    bay_segment,
    SampleTime,
    tn,
    chla,
    Sal_Top_ppth,
    Sal_Mid_ppth,
    Sal_Bottom_ppth
  ) |>
  mutate(
    sal = mean(c(Sal_Top_ppth, Sal_Mid_ppth, Sal_Bottom_ppth), na.rm = T),
    .by = c(bay_segment, SampleTime)
  ) |>
  mutate(
    date = floor_date(as.Date(SampleTime), unit = 'month')
  ) |>
  summarise(
    tn = mean(tn, na.rm = T),
    chla = mean(chla, na.rm = T),
    sal = mean(sal, na.rm = T),
    .by = c(bay_segment, date)
  ) |>
  mutate(
    dec_time = decimal_date(date),
    doy = yday(date),
    chla = ifelse(chla == 0, NA, chla),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) |>
  arrange(bay_segment, date) |>
  mutate(
    salorig = sal,
    sal = zoo::na.approx(salorig, x = date, na.rm = F, maxgap = 3),
    .by = c(bay_segment)
  )

save(wqdat, file = here('data/wqdat.RData'))

# salinity NA by bay segment
wqdat |>
  summarise(
    n_sal_na = sum(is.na(salorig)),
    n_tot = n(),
    .by = c(bay_segment)
  ) |>
  mutate(
    perc_na = n_sal_na / n_tot * 100
  )

# all bay segments -------------------------------------------------------

mods <- wqdat |>
  # filter(dec_time >= 2000) |>
  group_nest(bay_segment) |>
  mutate(
    mod = map(
      data,
      function(dat) {
        tomod <- dat |>
          mutate(
            res = log(chla),
            res = ifelse(is.infinite(res), NA, res),
            lim = 0
          ) |>
          select(date, res, flo = sal, lim) |>
          na.omit()

        out <- as.data.frame(tomod) |>
          tidal() |>
          modfit(flo_div = 40, fill_empty = T, tau = 0.5)

        return(out)
      }
    )
  )

save(mods, file = here('data/mods.RData'))

# residuals and hydro load -----------------------------------------------

mohydatraw <- rdataload(
  'https://github.com/tbep-tech/load-estimates/raw/refs/heads/main/data/mohydat.RData'
)

hydat <- mohydatraw |>
  filter(
    !bay_segment %in% c('All Segments (- N. BCB)', 'Remainder Lower Tampa Bay')
  ) |>
  rename(
    yr = year,
    mo = month
  ) |>
  mutate(
    bay_segment = factor(
      bay_segment,
      levels = c(
        'Old Tampa Bay',
        'Hillsborough Bay',
        'Middle Tampa Bay',
        'Lower Tampa Bay'
      ),
      labels = c('OTB', 'HB', 'MTB', 'LTB')
    ),
    date = make_date(year = yr, month = mo, day = 1),
    qrt = factor(
      lubridate::quarter(date),
      levels = 1:4,
      labels = c('JFM', 'AMJ', 'JAS', 'OND')
    ),
  )

save(hydat, file = here('data/hydat.RData'))
