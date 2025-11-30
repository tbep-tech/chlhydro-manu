library(tbeptools)
library(tidyverse)
library(mgcv)

# playing with models ----------------------------------------------------

# water quality
wqdat <- epcdata |> 
  filter(yr > 1974 & yr <= 2024) |> 
  select(bay_segment, SampleTime, tn, chla, Sal_Top_ppth, Sal_Mid_ppth, Sal_Bottom_ppth) |>
  mutate(
    sal = mean(c(Sal_Top_ppth, Sal_Mid_ppth, Sal_Bottom_ppth), na.rm = T), 
    .by = c(bay_segment, SampleTime)
  ) |> 
  mutate(
    date = floor_date(as.Date(SampleTime), unit = 'month')
  ) |> 
  summarise(
    tn = median(tn, na.rm = T), 
    chla = median(chla, na.rm = T), 
    sal = median(sal, na.rm = T),
    .by = c(bay_segment, date)
  ) |> 
  arrange(bay_segment, date)


tomod <- wqdat |> 
  filter(bay_segment == 'OTB') |> 
  mutate(
    dec_time = decimal_date(date), 
    doy = yday(date),
    chla = ifelse(chla == 0, NA, chla)
  )


mod <- gam(log10(chla) ~ te(dec_time, doy, sal, bs = c('tp', 'cc', 'tp')), data = tomod, knots = list(doy=c(1,366)), na.action = na.exclude)

summary(mod)$dev.expl

ggplot(tomod, aes(x = dec_time, y = chla)) +
  geom_point() +
  geom_line(aes(y = predict(mod)), color = 'red')

toprd <- crossing(
    date = seq.Date(as.Date('2000-01-01'), max(tomod$date), by = '1 month'),
    sal = seq(min(tomod$sal, na.rm = T), max(tomod$sal, na.rm = T), length.out = 10)
  ) |> 
  mutate(
    doy = yday(date),
    dec_time = lubridate::decimal_date(date)
  )

prds <- tibble(
    pred_chla = 10^predict(mod, newdata = toprd)
  ) |> 
  bind_cols(toprd)

gridgam_plo(prds, tomod, month = 'all', allsal = T)
gridgam_plo(prds, tomod, month = c(5:10), allsal = T, ncol = 6)
