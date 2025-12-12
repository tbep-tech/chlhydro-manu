# get datasets from repo
# try simple load, download if fail
rdataload <- function(dataurl = NULL){
  
  x <- gsub('\\.RData', '', basename(dataurl))
  
  # try simple load
  ld <- try(load(url(dataurl)), silent = T)
  
  # return x if load worked
  if(!inherits(ld, 'try-error')){
    out <- get(x)
  }

  # download x if load failed
  if(inherits(ld, 'try-error')){
    
    fl <- paste(tempdir(), basename(dataurl), sep = '/')
    download.file(flurl, destfile = fl, quiet = T)
    load(file = fl)
    out <- get(x)
    suppressMessages(file.remove(fl))
    
  }
  
  return(out)
    
}

# model rsq, same as summary(mod)$dev.expl
rsq_fun <- function(mod){
  
  # get complete cases
  toeval <- data.frame(resid = mod$residuals, obs = mod$y)
  toeval <- na.omit(toeval)
  
  ssr <- sum(toeval$resid^2)
  sst <- sum((toeval$obs - mean(toeval$obs))^2)
  
  out <- 1 - (ssr/sst)

  return(out)

}

# get gam predictions
pred_fun <- function(datin, modin){
  
  moddat <- datin
  moddat$fit <- predict(modin, type = 'link', newdata = moddat)
  moddat$btfit <- predict(modin, type = 'response', newdata = moddat)

  salgrd <- moddat |> 
    pull(sal) |>
    range(na.rm = TRUE)
  salgrd <- seq(salgrd[1], salgrd[2], length.out = 10)

  # make prediction grid
  toprd <- moddat |> 
    select(date, dec_time, doy) |> 
    crossing(
      sal = salgrd
    )

  # get link predictions, wide format
  fits <- tibble(
      res = as.numeric(predict(modin, newdata = toprd, type = 'link'))
    ) |> 
    bind_cols(toprd) |> 
    mutate(
      year = lubridate::year(date), 
      month = lubridate::month(date), 
      day = lubridate::day(date), 
      sal = factor(sal, levels = salgrd, labels = paste0('X', seq(1:length(salgrd))))
    ) |> 
    select(-dec_time, -doy) |> 
    pivot_wider(names_from = sal, values_from = res)

  # get response predictions, wide format
  btfits <- tibble(
      res = as.numeric(predict(modin, newdata = toprd, type = 'response'))
    ) |> 
    bind_cols(toprd) |> 
    mutate(
      year = lubridate::year(date), 
      month = lubridate::month(date), 
      day = lubridate::day(date), 
      sal = factor(sal, levels = salgrd, labels = paste0('X', seq(1:length(salgrd))))
    ) |> 
    select(-dec_time, -doy) |> 
    pivot_wider(names_from = sal, values_from = res)

  # normalized results
  moddat <- norm_fun(moddat, fits, btfits, salgrd)

  # predictions at mean salinity
  mod <- lm(sal ~ dec_time, data = moddat)
  toprd <- data.frame(
    dec_time = moddat$dec_time,
    doy = moddat$doy,
    sal = predict(mod, newdata = data.frame(dec_time = moddat$dec_time))
  )
  moddat$meansalfit <- predict(modin, newdata = toprd, type = 'response')

  out <- structure(
    .Data = moddat, 
    class = c('data.frame', 'tibble'),
    fits = fits,
    btfits = btfits,
    salgrd = salgrd
  )

  return(out)
  
}

# gridplot for gam
# ... additional arguments to pass to gridpred_fun
grid_plo <- function(prds, month = c(1:12), years = NULL, col_vec = NULL, col_lim = NULL, salscl = TRUE, allsal = FALSE, sal_fac = 3, yr_fac = 3, ncol = NULL, grids = FALSE, pretty = TRUE, ...){

  # convert month vector to those present in data
  allmo <- FALSE
  if('all' %in% month){ 
    allmo <- TRUE
    month <- c(1:12)
  }

  # format predictions as wide
  to_plo <- attr(prds, 'btfits')
   
  # model data with date column
  moddat <- data.frame(prds) |> 
    filter(lubridate::year(date) >= min(to_plo$year) & lubridate::year(date) <= max(to_plo$year))

  # salinity grid values
  salgrd <- attr(prds, 'salgrd')
  
  # convert month vector to those present in data
  month <- month[month %in% to_plo$month]
  if(length(month) == 0) stop('No observable data for the chosen month')
  
  # axis labels
  ylabel <- 'Chlorophyll-a (µg/L)'
  xlabel <- 'Salinity (ppth)'

  # subset years to plot
  if(!is.null(years)){
   
    if(length(years) != 2)
      stop('years argument must have two values for first and last')
  
  } else{

    years <- range(to_plo$year)

  }
  
  years <- seq(years[1], years[2])
  to_plo <- to_plo[to_plo$year %in% years, ]
    
  if(nrow(to_plo) == 0) stop('No data to plot for the date range')

  # reshape data frame
  to_plo <- to_plo[to_plo$month %in% month, , drop = FALSE]
  names(to_plo)[grep('^X', names(to_plo))] <- paste('sal', salgrd)
  to_plo <- tidyr::gather(to_plo, 'sal', 'res', 5:ncol(to_plo)) |> 
    mutate(sal = as.numeric(gsub('^sal ', '', sal))) |> 
    select(-date, -day) |>  
    summarize(
      res = mean(res, na.rm = TRUE), 
      .by = c(year, month, sal)
    )
  
  # change sal to original scale
  if(!salscl){
   
    # grid data
    salobs_rng <- range(moddat$sal, na.rm = TRUE)
    salscl_rng <- range(to_plo$sal, na.rm = TRUE)
    to_plo$sal <- (to_plo$sal - salscl_rng[1]) / diff(salscl_rng) * diff(salobs_rng) + salobs_rng[1]
    
    #input data
    salscl_rng <- range(prddat$sal, na.rm = TRUE)
    prddat$sal <- (prddat$sal - salscl_rng[1]) / diff(salscl_rng) * diff(salobs_rng) + salobs_rng[1]
    
    # salgrd to raw scale
    salgrd <- seq(salobs_rng[1], salobs_rng[2], length = length(salgrd))
    
  }
    
  ## use linear interpolation to make a smoother plot
  if(!allmo){
    
    # these are factors by which salinity and years are multiplied for interpolation
    sal_fac <- length(salgrd) * sal_fac
    sal_fac <- seq(min(salgrd), max(salgrd), length.out = sal_fac)
    yr_fac <- length(unique(to_plo$year)) * yr_fac
    yr_fac <- seq(min(to_plo$year), max(to_plo$year), length.out = yr_fac)
    
    # separately by month
    to_plo <- split(to_plo, to_plo$month)

    to_plo <- lapply(to_plo, function(x){

      # interp across salinity first
      interped <- lapply(
        split(x, x$year), 
        function(y){
          out <- approx(y$sal, y$res, xout = sal_fac, rule = 2)
          out <- data.frame(year = unique(y$year), month = unique(y$month), out)
          return(out)
        })
      interped <- do.call('rbind', interped)
      names(interped) <- c('year', 'month', 'sal', 'res')
      
      # interp across years
      interped <- lapply(
        split(interped, interped$sal), 
        function(y){
          out <- approx(y$year, y$res, xout = yr_fac, rule = 2)
          out <- data.frame(year = out$x, month = unique(y$month), sal = unique(y$sal), res = out$y)
          return(out)
        })
      interped <- do.call('rbind', interped)
      names(interped) <- c('year', 'month', 'sal', 'res')
    
      return(interped)
    
    })
    
    to_plo <- do.call('rbind', to_plo)
    row.names(to_plo) <- 1:nrow(to_plo)
      
  }
  
  ## use linear interpolation to make a smoother plot, handle differently if allmo
  if(allmo){

    # format to_plo for interp (wide, as matrix)
    to_interp <- to_plo
    to_interp$date <- with(to_interp, as.Date(paste(year, month, '1', sep = '-')))
    to_interp <- ungroup(to_interp) |> 
      select(date, sal, res) |> 
      tidyr::spread(sal, res)
    
    # values to pass to interp
    dts <- lubridate::decimal_date(to_interp$date)
    fit_grd <- select(to_interp, -date)
    sal_fac <- length(salgrd) * sal_fac
    sal_fac <- seq(min(salgrd), max(salgrd), length.out = sal_fac)
    yr_fac <- seq(min(dts), max(dts), length.out = length(dts) *  yr_fac)
    to_norm <- expand.grid(yr_fac, sal_fac)
          
    # bilinear interpolation of fit grid with data to average for norms
    norms <- fields::interp.surface(
      obj = list(
        y = salgrd,
        x = dts,
        z = data.frame(fit_grd)
      ), 
      loc = to_norm
    )
    
    to_plo <- data.frame(to_norm, norms)
    names(to_plo) <- c('year', 'sal', 'res')
    
  }
  
  # constrain plots to salinity limits for the selected month
  if(!allsal & !allmo){
    
    #min, max salinity values to plot
    lim_vals<- moddat |> 
      mutate(
        month = lubridate::month(date),
        yr = lubridate::year(date)
      ) |> 
      filter(yr %in% years) |> 
      summarize(
        Low = quantile(sal, 0, na.rm = TRUE),
        High = quantile(sal, 1, na.rm = TRUE), 
        .by = month
      )
  
    # month sal ranges for plot
    lim_vals <- lim_vals[lim_vals$month %in% month, ]
    
    # merge limts with months
    to_plo <- left_join(to_plo, lim_vals, by = 'month')
    
    # reduce data
    sel_vec <- with(to_plo, 
      sal >= Low &
      sal <= High
      )
    to_plo <- to_plo[sel_vec, !names(to_plo) %in% c('Low', 'High')]
    to_plo <- arrange(to_plo, year, month)
    
  }
  
  # contstrain all data by quantiles if not separated by month    
  if(!allsal & allmo){
   
    quants <- moddat |> 
      pull(sal) |> 
      quantile(c(0.05, 0.95), na.rm = TRUE)
    to_plo <- to_plo[with(to_plo, sal >= quants[1] & sal <= quants[2]), ]
     
  }
  
  # change month vector of not plotting all months in same plot
  if(!allmo){
    # months labels as text
    mo_lab <- data.frame(
      num = seq(1:12), 
      txt = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')
    )
    mo_lab <- mo_lab[mo_lab$num %in% month, ]
    to_plo$month <- factor(to_plo$month, levels =  mo_lab$num, labels = mo_lab$txt)
  } 
  
  # make plot
  p <- ggplot(to_plo, aes(x = year, y = sal, fill = res)) + 
    geom_tile(data = subset(to_plo, !is.na(to_plo$res)), aes(fill = res)) +
    geom_tile(data = subset(to_plo,  is.na(to_plo$res)), fill = 'black', alpha = 0)
  
  if(!allmo) p <- p + facet_wrap(~month, ncol = ncol)
  
  # return bare bones if FALSE
  if(!pretty) return(p)
  
  # get colors
  cols <- cols <- RColorBrewer::brewer.pal(11, 'Spectral')
  
  p <- p +
    theme_minimal() +
    theme(
      legend.position = 'top',
      axis.title.x = element_blank(), 
      plot.title = element_text(hjust = 0.5)
      )  +
    # scale_x_continuous(expand = c(0.05, 0.05)) +
    scale_y_continuous(xlabel) + #, expand = c(0.05,0.05)) +
    scale_fill_gradientn(ylabel, colours = rev(cols), limits = col_lim) +
    guides(fill = guide_colourbar(barwidth = 10)) 
    
  # add grid lines
  if(!grids) 
    p <- p + 
      theme(      
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  
  return(p)
    
}

# normalize predictions
norm_fun <- function(dat_in, fits, btfits, salgrd){

  num_obs <- nrow(dat_in)

  # prep interp grids by adding month, year columns
  dts <- fits$date
  fits <- select(fits, -year, -month, -day, -date)
  btfits <- btfits |> 
    select(-year, -date, -month, -day)

  # sal values occuring by month, used for interpolation
  sal_mon <- data.frame(dat_in) |> 
    select(date, sal) |> 
    mutate(
      yr = lubridate::year(date), 
      mo = lubridate::month(date)
    ) |> 
    select(-date) |> 
    summarize(
      sal = mean(sal, na.rm = TRUE), 
      .by = c(yr, mo)
    ) |> 
    tidyr::spread(yr, sal)
  
  # values to interpolate for normalization
  to_norm <- data.frame(dat_in) |> 
    select(date) |> 
    mutate(mo = lubridate::month(date)) |> 
    left_join(sal_mon, by = 'mo') |> 
    tidyr::gather('yr', 'sal', -date, -mo) |> 
    arrange(date) |> 
    select(date, sal) |> 
    na.omit()

  # bilinear interpolatoin of fit grid with data to average for norms
  norms <- fields::interp.surface(
    obj = list(
      y = salgrd,
      x = dts,
      z = data.frame(fits)
    ), 
    loc = to_norm
  )
  
  # append to normalization data, then average for unique dates
  # averaging happens only if more than 80% of the predictions for a date are filled
  norms <- data.frame(date = to_norm$date, norms) |> 
    summarise(
      norm = ifelse(sum(is.na(norms))/length(norms) > 0.2, NA, mean(norms, na.rm = TRUE)), 
      .by = date
    )
   
  # bilinear interpolatoin of fit grid with data to average for norms
  btnorms <- fields::interp.surface(
    obj = list(
      y = salgrd,
      x = dts,
      z = data.frame(btfits)
    ), 
    loc = to_norm
  )
    
  # append to normalization data, then average for unique dates
  # averaging happens only if more than 80% of the predictions for a date are filled
  btnorms <- data.frame(date = to_norm$date, btnorms) |> 
    summarise(
      btnorm = ifelse(sum(is.na(btnorms))/length(btnorms) > 0.2, NA, mean(btnorms, na.rm = TRUE)),
      .by = date
    )
   
  # append to dat_in object
  dat_in <- dplyr::left_join(dat_in, norms, by = 'date')
  dat_in <- dplyr::left_join(dat_in, btnorms, by = 'date')

  # exit function
  return(dat_in)
    
}
