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

# gridplot for gam
gridgam_plo <- function(prddat, wqdat, month = c(1:12), years = NULL, col_vec = NULL, col_lim = NULL, salscl = TRUE, allsal = FALSE, sal_fac = 3, yr_fac = 3, ncol = NULL, grids = FALSE, pretty = TRUE, ...){
 
  # convert month vector to those present in data
  allmo <- FALSE
  if('all' %in% month){ 
    allmo <- TRUE
    month <- c(1:12)
  }

  # salinity grid values
  sal_grd <- unique(prddat$sal)

  # get the selected months
  to_plo <- prddat |> 
    mutate(
      date = as.Date(date), 
      year = lubridate::year(date), 
      month = lubridate::month(date), 
      day = lubridate::day(date), 
      sal = factor(sal, levels = sal_grd, labels = paste0('X', seq(1:length(sal_grd))))
    ) |> 
    select(-doy, -dec_time) |> 
    pivot_wider(names_from = sal, values_from = pred_chla)

  # convert month vector to those present in data
  month <- month[month %in% to_plo$month]
  if(length(month) == 0) stop('No observable data for the chosen month')
  
  # axis labels
  ylabel <- 'Chlorophyll-a (µg/L)'
  xlabel <- NULL

  # subset years to plot
  if(!is.null(years)){
   
    if(length(years) != 2)
      stop('years argument must have two values for first and last')
  
    years <- seq(years[1], years[2])
    to_plo <- to_plo[to_plo$year %in% years, ]
     
    if(nrow(to_plo) == 0) stop('No data to plot for the date range')
  
  }
  
  # reshape data frame
  to_plo <- to_plo[to_plo$month %in% month, , drop = FALSE]
  names(to_plo)[grep('^X', names(to_plo))] <- paste('sal', sal_grd)
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
    salobs_rng <- range(wqdat$sal, na.rm = TRUE)
    salscl_rng <- range(to_plo$sal, na.rm = TRUE)
    to_plo$sal <- (to_plo$sal - salscl_rng[1]) / diff(salscl_rng) * diff(salobs_rng) + salobs_rng[1]
    
    #input data
    salscl_rng <- range(prddat$sal, na.rm = TRUE)
    prddat$sal <- (prddat$sal - salscl_rng[1]) / diff(salscl_rng) * diff(salobs_rng) + salobs_rng[1]
    
    # sal_grd to raw scale
    sal_grd <- seq(salobs_rng[1], salobs_rng[2], length = length(sal_grd))
    
  }
    
  ## use linear interpolation to make a smoother plot
  if(!allmo){
    
    # these are factors by which salinity and years are multiplied for interpolation
    sal_fac <- length(sal_grd) * sal_fac
    sal_fac <- seq(min(sal_grd), max(sal_grd), length.out = sal_fac)
    yr_fac <- length(unique(to_plo$year)) * yr_fac
    yr_fac <- seq(min(to_plo$year), max(to_plo$year), length.out = yr_fac)
    
    # separately by month
    to_plo <- split(to_plo, to_plo$month)

    to_plo <- lapply(to_plo, function(x){
      
      # interp across salinity first
      interped <- lapply(
        split(x, x$year), 
        function(y){
          out <- approx(y$sal, y$res, xout = sal_fac)
          out <- data.frame(year = unique(y$year), month = unique(y$month), out)
          return(out)
        })
      interped <- do.call('rbind', interped)
      names(interped) <- c('year', 'month', 'sal', 'res')
      
      # interp across years
      interped <- lapply(
        split(interped, interped$sal), 
        function(y){
          out <- approx(y$year, y$res, xout = yr_fac)
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
    dts <- dec_time(to_interp$date)$dec_time
    fit_grd <- select(to_interp, -date)
    sal_fac <- length(sal_grd) * sal_fac
    sal_fac <- seq(min(sal_grd), max(sal_grd), length.out = sal_fac)
    yr_fac <- seq(min(dts), max(dts), length.out = length(dts) *  yr_fac)
    to_norm <- expand.grid(yr_fac, sal_fac)
          
    # bilinear interpolation of fit grid with data to average for norms
    norms <- fields::interp.surface(
      obj = list(
        y = sal_grd,
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
    lim_vals<- tomod |> 
      mutate(
        month = lubridate::month(date)
      ) |> 
      filter(lubridate::year(date) >= min(to_plo$year) & lubridate::year(date) <= max(to_plo$year)) |>
      summarize(
        Low = quantile(sal, 0.05, na.rm = TRUE),
        High = quantile(sal, 0.95, na.rm = TRUE), 
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
   
    quants <- wqdat |> 
      filter(lubridate::year(date) >= min(to_plo$year) & lubridate::year(date) <= max(to_plo$year)) |>
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
      axis.title.x = element_blank()
      )  +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(xlabel, expand = c(0,0)) +
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
