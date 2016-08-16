#' Generate summary warm/cold rain mask files from 3D ECHAM output
#'
#' Prerequisite is ECHAM model output that contains rain and snow
#' vertical profiles as variables \code{aprlv_na} and \code{aprsv_na}.
#' Files should be created at monthly time frequency (but data time
#' frequency is arbitrary).  As a side effect of calling this
#' function, one RDS file per ECHAM file is created containing the 2D
#' precip mask at each time in the input file.  The precipitation mask
#' can take the values \code{dry}, \code{warm}, \code{cold}, and
#' \code{snow}.
#' 
#' At the end of processing, a summary file containing the number of
#' occurrences of dry/warm/cold/snow in each grid box is produced.
#' This file is named \code{experiment}.rds and is written in the
#' working directory.
#' 
#' @param datadir Path to data; this directory should containt \code{experiment} as a subdirectory.
#' @param experiment Subdirectory of \code{datadir} that contains the experiment to process.
#' @param years Numeric.  Vector of years to process.
#' @param ncores Single number.  Number of children to use in processing the data.
#' @return \code{NULL}
#' @export
process.precip.profile.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             experiment = "amip-rain-15", years = 1979:1983,
             ncores = 36) {
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            ## expand.grid(year = 2000, month = 1) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname <- sprintf("%s/%s/%s_%d%02d.01_rain3d.nc",
                                     datadir, experiment, experiment, year, month)
                    out.name <- gsub(".nc", ".rds", fname)
                    df <- try(readRDS(out.name))
                    if (class(df) == "try-error") {
                        nc <- ncdf4::nc_open(fname)
                        t <- ncdf4::ncvar_get(nc, "time")
                        lon <- ncdf4::ncvar_get(nc, "lon")
                        lat <- ncdf4::ncvar_get(nc, "lat")
                        df <- plyr::ldply(1:length(t), function(i) {
                            qr <- ncdf4::ncvar_get(nc, "aprlv_na", start = c(1,1,1,i), count = c(-1,-1,-1,1)) ## vertically resolved liquid precip mixing ratio
                            qs <- ncdf4::ncvar_get(nc, "aprsv_na", start = c(1,1,1,i), count = c(-1,-1,-1,1)) ## vertically resolved solid precip mixing ratio
                            rain.mask <- apply(qr, c(1,2), function(x) any(x > 1e-7)) ## 1e-7 is the cutoff in CESM... ECHAM appears not to have a cutoff
                            snow.mask <- apply(qs, c(1,2), function(x) any(x > 1e-7))
                            mask <- factor(ifelse(rain.mask, ifelse(snow.mask, "cold", "warm"),
                                           ifelse(snow.mask, "snow", "dry")),
                                           levels = c("dry", "warm", "cold", "snow"))
                            df <- expand.grid(lon = as.vector(lon),
                                              lat = as.vector(lat)) %>%
                                cbind(time = t[i], mask = as.vector(mask))
                        }, .parallel = FALSE, .progress = "text")
                        ncdf4::nc_close(nc)
                        saveRDS(df, out.name)
                    }
                    df %>%
                        plyr::ddply(~ lon + lat, function(x) table(x$mask))
                })
            }, .parallel = TRUE) -> df
        
        saveRDS(df, sprintf("%s.rds", experiment))
    }

#' @export
process.pr.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             experiment = "amip-rain-15", years = 1979:1983,
             ncores = 12) {
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname <- sprintf("%s/%s/%s_%d%02d.01_rain2d.nc",
                                     datadir, experiment, experiment, year, month)
                    nc <- ncdf4::nc_open(fname)
                    t <- ncdf4::ncvar_get(nc, "time")
                    lon <- ncdf4::ncvar_get(nc, "lon")
                    lat <- ncdf4::ncvar_get(nc, "lat")
                    pr <- ncdf4::ncvar_get(nc, "aprl_na")
                    prc <- ncdf4::ncvar_get(nc, "aprc_na")
                    ncdf4::nc_close(nc)
                    df <- expand.grid(lon = as.vector(lon),
                                      lat = as.vector(lat),
                                      t = as.vector(t)) %>%
                        cbind(pr = as.vector(pr),
                              prc = as.vector(prc)) %>%
                        plyr::ddply(~ lon + lat, function(x) {
                            x %>%
                                dplyr::mutate(pr.class = cut(3600 * (pr - prc),
                                                             c(0, 0.01, 0.1, 1, 1e8), right = FALSE)) %>%
                                dplyr::group_by(pr.class) %>%
                                dplyr::summarize(n = n(),
                                                 sum.pr = sum(pr - prc))
                        })
                    ##         filter(as.integer(pr.class) > 1 & !is.na(pr.class)) %>%
                    
                    ##         ## calculate totals
                    ##         mutate(sum.occ = n(),
                    ##                sum.pr = sum(pr - prc)) %>%
                    ##         ## calculate by class
                    ##         group_by(pr.class) %>%
                    ##         summarize(N = sum.occ[1],
                    ##                   Pr = sum.pr[1],
                    ##                   frac.occ = n() / N,
                    ##                   frac.pr = sum(pr - prc) / Pr)
                    ## })
                })
            }, .progress = "text", .parallel = TRUE) -> df
        saveRDS(df, sprintf("pr-hist-%s.rds", experiment))
    }

#' @export
process.rad.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             experiment = "amip-rain-15", years = 1979:1983,
             ncores = 12) {
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname <- sprintf("%s/%s/%s_%d%02d.01_echamm.nc",
                                     datadir, experiment, experiment, year, month)
                    nc <- ncdf4::nc_open(fname)
                    t <- ncdf4::ncvar_get(nc, "time")
                    lon <- ncdf4::ncvar_get(nc, "lon")
                    lat <- ncdf4::ncvar_get(nc, "lat")
                    srad0d <- ncdf4::ncvar_get(nc, "srad0d")
                    srad0  <- ncdf4::ncvar_get(nc, "srad0")
                    trad0  <- ncdf4::ncvar_get(nc, "trad0")
                    sraf0  <- ncdf4::ncvar_get(nc, "sraf0") ## clear-sky 
                    traf0  <- ncdf4::ncvar_get(nc, "traf0")
                    aprl  <- ncdf4::ncvar_get(nc, "aprl") ## precip
                    aprc  <- ncdf4::ncvar_get(nc, "aprc")
                    aprs  <- ncdf4::ncvar_get(nc, "aprs")
                    ncdf4::nc_close(nc)
                    df <- expand.grid(lon = as.vector(lon),
                                      lat = as.vector(lat),
                                      t = as.vector(t)) %>%
                        cbind(srad0d = as.vector(srad0d ),
                              srad0  = as.vector(srad0  ),
                              trad0  = as.vector(trad0  ),
                              sraf0  = as.vector(sraf0  ),
                              traf0  = as.vector(traf0  ),
                              aprl  = as.vector(aprl  ),
                              aprc  = as.vector(aprc  ),
                              aprs  = as.vector(aprs  ))

                })
            }, .progress = "text", .parallel = FALSE) -> df
        saveRDS(df, sprintf("rad-%s.rds", experiment))
    }
