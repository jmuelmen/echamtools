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
             experiment = "amip-rain-15", out.prefix = "",
             years = 1979:1983,
             ncores = 12,
             flux = TRUE,
             subsample = NULL) {
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            ## expand.grid(year = 2000, month = 1) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname <- sprintf("%s/%s/%s_%d%02d.01_rain3d.nc",
                                     datadir, experiment, experiment, year, month)
                    out.name <- gsub(".nc", ifelse(flux, ".rds", "-mr.rds"), fname)
                    df <- try(readRDS(out.name), silent = TRUE)
                    if (class(df) == "try-error") {
                        nc <- try(ncdf4::nc_open(fname), silent = TRUE)
                        if (class(nc) == "try-error")
                            return(NULL)
                        t <- ncdf4::ncvar_get(nc, "time")
                        lon <- ncdf4::ncvar_get(nc, "lon")
                        lat <- ncdf4::ncvar_get(nc, "lat")
                        df <- plyr::ldply(c(1e-7, 0.01 / 24 / 3600, 0.01 / 3600, 0.1 / 3600), function(flux.thresh) {
                            plyr::ldply(1:length(t), function(i) {
                                if (flux) { ## mask based on fluxes
                                    qr <- ncdf4::ncvar_get(nc, "aprlv_na", start = c(1,1,1,i), count = c(-1,-1,-1,1)) ## vertically resolved liquid precip flux
                                    qs <- ncdf4::ncvar_get(nc, "aprsv_na", start = c(1,1,1,i), count = c(-1,-1,-1,1)) ## vertically resolved solid precip flux
                                    rain.mask <- apply(qr, c(1,2), function(x) any(x > flux.thresh)) ## 1e-7 is the cutoff in CESM... ECHAM appears not to have a cutoff
                                    snow.mask <- apply(qs, c(1,2), function(x) any(x > flux.thresh))
                                } else { ## mask based on mixing ratios
                                    qr <- ncdf4::ncvar_get(nc, "xrl_na", start = c(1,1,1,i), count = c(-1,-1,-1,1)) ## vertically resolved liquid precip mixing ratio
                                    qs <- ncdf4::ncvar_get(nc, "xsl_na", start = c(1,1,1,i), count = c(-1,-1,-1,1)) ## vertically resolved solid precip mixing ratio
                                    rain.mask <- apply(qr, c(1,2), function(x) any(x > 1e-12))
                                    snow.mask <- apply(qs, c(1,2), function(x) any(x > 1e-12))
                                }
                                mask <- factor(ifelse(rain.mask, ifelse(snow.mask, "cold", "warm"),
                                               ifelse(snow.mask, "snow", "dry")),
                                               levels = c("dry", "warm", "cold", "snow"))
                                df <- expand.grid(lon = as.vector(lon),
                                                  lat = as.vector(lat)) %>%
                                    dplyr::mutate(time = t[i],
                                                  mask = mask,
                                                  flux.thresh = flux.thresh)
                            }, .parallel = FALSE, .progress = "none")
                        }, .parallel = FALSE, .progress = "none")
                        ncdf4::nc_close(nc)
                        saveRDS(df, out.name)
                    }
                    print(out.name)
                    str(df)
                    ## if subsampling is requested, do it here
                    (if (is.null(subsample))
                         df
                     else {
                         t <- unique(df$time)
                         dplyr::filter(df, time %in% t[seq(1, length(t), by = subsample)])
                     }) %>%
                        plyr::ddply(~ lon + lat + flux.thresh, function(x) table(x$mask))
                })
            }, .parallel = TRUE) -> df

        str(df)
        saveRDS(df, sprintf("%s%s%s.rds", out.prefix, experiment, ifelse(flux, "", "-mr")))
        df
    }

#' @export
process.precip.cosp.profile.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             experiment = "amip-rain-15", out.prefix = "",
             years = 1979:1983,
             ncores = 36,
             flux = TRUE,
             subsample = NULL) {
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            ## expand.grid(year = 2000, month = 1) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname.lsrain   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_lsrain"      )
                    fname.lssnow   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_lssnow"      )
                    fname.tau      <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_cisccp_tau3d")
                    fname.reffl    <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_reffl"       )
                    fname.reffi    <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_reffi"       )
                    fname.aclc     <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_aclc"       )
                    fname.tm1      <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_tm1"       )
                    fname.tm1_cosp <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_tm1_cosp"       )
                    fname.xl       <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_xl"       )
                    fname.xi       <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_xi"       )
                    fname.rain3d   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "rain3d"       )
                    fname.rain2d   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "rain2d"       )
                    fname.dbze     <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_001"       )

                    out.name <- gsub(".nc", "-cosp.rds", fname)

                    ## if (any(class(nc.lssnow   <- try(ncdf4::nc_open(fname.lssnow  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.lssnow   ), add = TRUE)
                    ## if (any(class(nc.lsrain   <- try(ncdf4::nc_open(fname.lsrain  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.lsrain   ), add = TRUE)
                    ## if (any(class(nc.tau      <- try(ncdf4::nc_open(fname.tau     ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.tau      ), add = TRUE)
                    ## if (any(class(nc.reffl    <- try(ncdf4::nc_open(fname.reffl   ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.reffl    ), add = TRUE)
                    ## if (any(class(nc.reffi    <- try(ncdf4::nc_open(fname.reffi   ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.reffi    ), add = TRUE)
                    if (any(class(nc.aclc     <- try(ncdf4::nc_open(fname.aclc    ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.aclc     ), add = TRUE)
                    ## if (any(class(nc.tm1      <- try(ncdf4::nc_open(fname.tm1     ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.tm1      ), add = TRUE)
                    if (any(class(nc.tm1_cosp <- try(ncdf4::nc_open(fname.tm1_cosp), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.tm1_cosp ), add = TRUE)
                    if (any(class(nc.xl       <- try(ncdf4::nc_open(fname.xl      ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.xl       ), add = TRUE)
                    if (any(class(nc.xi       <- try(ncdf4::nc_open(fname.xi      ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.xi       ), add = TRUE)
                    if (any(class(nc.rain3d   <- try(ncdf4::nc_open(fname.rain3d  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.rain3d   ), add = TRUE)
                    if (any(class(nc.rain2d   <- try(ncdf4::nc_open(fname.rain2d  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.rain2d   ), add = TRUE)
                    if (any(class(nc.dbze     <- try(ncdf4::nc_open(fname.dbze    ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.dbze     ), add = TRUE)

                    lon  <- ncdf4::ncvar_get(nc.rain3d, "lon")
                    lat  <- ncdf4::ncvar_get(nc.rain3d, "lat")
                    lev  <- ncdf4::ncvar_get(nc.rain3d, "mlev")
                    time <- ncdf4::ncvar_get(nc.rain3d, "time")

                    df <- plyr::ldply(1:length(time), function(i) {
                        ## get aclc first for further processing
                        aclc  <- ncdf4::ncvar_get(nc.aclc  , "aclc", start = c(1,1,1,i), count = c(-1,-1,-1,1))
                        layer <- apply(aclc > 0, 1:2, label.vertical.features) %>% aperm(c(2,3,1))
                        aprl <- rep(ncdf4::ncvar_get(nc.rain2d  , "aprl_na", start = c(1,1,i), count = c(-1,-1,1)), length(lev))
                        aprs <- rep(ncdf4::ncvar_get(nc.rain2d  , "aprs_na", start = c(1,1,i), count = c(-1,-1,1)), length(lev))

                        df <- expand.grid(lon = as.vector(lon),
                                          lat = as.vector(lat),
                                          lev = as.vector(lev)) %>%
                            dplyr::mutate(aprlv    = as.vector(ncdf4::ncvar_get(nc.rain3d  , "aprlv_na", start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          aprsv    = as.vector(ncdf4::ncvar_get(nc.rain3d  , "aprsv_na", start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          dbze     = as.vector(ncdf4::ncvar_get(nc.dbze    , "dbze94_001"    , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          aclc     = as.vector(aclc),
                                          tm1_cosp = as.vector(ncdf4::ncvar_get(nc.tm1_cosp, "tm1_cosp", start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          xl       = as.vector(ncdf4::ncvar_get(nc.xl      , "xl"      , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          xi       = as.vector(ncdf4::ncvar_get(nc.xi      , "xi"      , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          layer    = as.vector(layer),
                                          fracout  = as.vector(ncdf4::ncvar_get(nc.dbze    , "frac_out_001" , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          aprl = aprl,
                                          aprs = aprs) %>%
                            dplyr::group_by(lon, lat) %>%
                            dplyr::summarize(aprl = aprl[1],
                                             aprs = aprs[1],
                                             dbze.max = max(dbze),
                                             aprlv.max = max(aprlv),
                                             aprsv.max = max(aprsv),
                                             aprlv.max.dbze = aprlv[which.max(aprlv)],
                                             highest.drizzle.layer = ifelse(any(layer > 0 & dbze > -15),
                                                                            min(layer[layer > 0 & dbze > -15]),
                                                                            NA),
                                             highest.rain.layer = ifelse(any(layer > 0 & dbze > 0),
                                                                         min(layer[layer > 0 & dbze > 0]),
                                                                         NA),
                                             cold.drizzle = ifelse(!is.na(highest.drizzle.layer),
                                                                   any(xi[layer == highest.drizzle.layer] > 1e-7),
                                                                   NA),
                                             cold.rain = ifelse(!is.na(highest.rain.layer),
                                                                any(xi[layer == highest.rain.layer] > 1e-7),
                                                                NA)) %>%
                            dplyr::ungroup() %>%
                            dplyr::mutate(time = time[i])
                        ## print(dim(ncdf4::ncvar_get(nc.rain2d  , "aprl_na", start = c(1,1,i), count = c(-1,-1,1))))
                        df
                    }, .parallel = FALSE, .progress = "none")
                    saveRDS(df, out.name)
                    df
                })
            }, .parallel = TRUE) -> df
        
        saveRDS(df, sprintf("%scosp-%s.rds", out.prefix, experiment))
    }

#' @export
process.cfodd.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             experiment = "amip-rain-15", out.prefix = "",
             years = 1979:1983,
             ncores = 12,
             flux = TRUE,
             subsample = NULL) {
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            ## expand.grid(year = 2000, month = 1) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname.lsrain   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_lsrain"      )
                    fname.lssnow   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_lssnow"      )
                    fname.tau      <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_cisccp_tau3d")
                    fname.reffl    <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_reffl"       )
                    fname.reffi    <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_reffi"       )
                    fname.aclc     <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_aclc"       )
                    fname.tm1      <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_tm1"       )
                    fname.tm1_cosp <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_tm1_cosp"       )
                    fname.xl       <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_xl"       )
                    fname.xi       <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_xi"       )
                    fname.rain3d   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "rain3d"       )
                    fname.rain2d   <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "rain2d"       )
                    fname.dbze     <- sprintf("%s/%s/%s_%d%02d.01_%s.nc", datadir, experiment, experiment, year, month, "cosp_001"       )
                    
                    if (any(class(nc.lssnow   <- try(ncdf4::nc_open(fname.lssnow  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.lssnow   ), add = TRUE)
                    if (any(class(nc.lsrain   <- try(ncdf4::nc_open(fname.lsrain  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.lsrain   ), add = TRUE)
                    if (any(class(nc.tau      <- try(ncdf4::nc_open(fname.tau     ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.tau      ), add = TRUE)
                    if (any(class(nc.reffl    <- try(ncdf4::nc_open(fname.reffl   ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.reffl    ), add = TRUE)
                    if (any(class(nc.reffi    <- try(ncdf4::nc_open(fname.reffi   ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.reffi    ), add = TRUE)
                    if (any(class(nc.aclc     <- try(ncdf4::nc_open(fname.aclc    ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.aclc     ), add = TRUE)
                    ## if (any(class(nc.tm1      <- try(ncdf4::nc_open(fname.tm1     ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.tm1      ), add = TRUE)
                    if (any(class(nc.tm1_cosp <- try(ncdf4::nc_open(fname.tm1_cosp), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.tm1_cosp ), add = TRUE)
                    ## if (any(class(nc.xl       <- try(ncdf4::nc_open(fname.xl      ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.xl       ), add = TRUE)
                    ## if (any(class(nc.xi       <- try(ncdf4::nc_open(fname.xi      ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.xi       ), add = TRUE)
                    if (any(class(nc.rain3d   <- try(ncdf4::nc_open(fname.rain3d  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.rain3d   ), add = TRUE)
                    ## if (any(class(nc.rain2d   <- try(ncdf4::nc_open(fname.rain2d  ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.rain2d   ), add = TRUE)
                    if (any(class(nc.dbze     <- try(ncdf4::nc_open(fname.dbze    ), silent = TRUE)) == "try-error")) return(NULL); on.exit(ncdf4::nc_close(nc.dbze     ), add = TRUE)

                    lon  <- ncdf4::ncvar_get(nc.reffi, "lon")
                    lat  <- ncdf4::ncvar_get(nc.reffi, "lat")
                    lev  <- ncdf4::ncvar_get(nc.reffi, "mlev")
                    time <- ncdf4::ncvar_get(nc.reffi, "time")
                    
                    df <- plyr::ldply(1:length(time), function(i) {
                        ## get fracout and reffi first for masking
                        reffi   <- ncdf4::ncvar_get(nc.reffi  , "reffi", start = c(1,1,1,i), count = c(-1,-1,-1,1))
                        fracout <- ncdf4::ncvar_get(nc.dbze   , "frac_out_001" , start = c(1,1,1,i), count = c(-1,-1,-1,1))
                        mask <- (apply(reffi, c(1,2), function(x) rep(all(x == 0 | x == 4), 31)) &
                                 apply(fracout, c(1,2), function(x) rep(any(x == 1), 31)))  %>%
                            aperm(c(2,3,1))
                        ## compute layers based on fracout or on aclc > 0?
                        layer <- apply(fracout, 1:2, label.vertical.features) %>% aperm(c(2,3,1))

                        ## try to get pre-precip effective radius (which we did not archive for all runs)
                        reffl.pre <- try(ncdf4::ncvar_get(nc.rain3d, "reffl_pre_na", start = c(1,1,1,i), count = c(-1,-1,-1,1)),
                                         silent = TRUE) 
                        
                        df <- expand.grid(lon = as.vector(lon),
                                          lat = as.vector(lat),
                                          lev = as.vector(lev)) %>%
                            dplyr::mutate(aclc     = as.vector(ncdf4::ncvar_get(nc.aclc, "aclc"     , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          layer    = as.vector(layer),
                                          tm1_cosp = as.vector(ncdf4::ncvar_get(nc.tm1_cosp, "tm1_cosp"     , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          lssnow   = as.vector(ncdf4::ncvar_get(nc.lssnow  , "lssnow"     , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          lsrain   = as.vector(ncdf4::ncvar_get(nc.lsrain  , "lsrain"     , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          tau      = as.vector(ncdf4::ncvar_get(nc.tau     , "cisccp_tau3d"     , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          dbze     = as.vector(ncdf4::ncvar_get(nc.dbze    , "dbze94_001"     , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          reffl    = as.vector(ncdf4::ncvar_get(nc.reffl   , "reffl"     , start = c(1,1,1,i), count = c(-1,-1,-1,1))),
                                          reffi    = as.vector(reffi),
                                          fracout  = as.vector(fracout))

                        if (!any(class(reffl.pre) == "try-error")) {
                            df %<>%
                                dplyr::mutate(reffl.pre     = as.vector(reffl.pre))
                        }

                        df %>%
                            dplyr::filter(as.vector(mask)) %>%
                            dplyr::group_by(lon, lat) %>%
                            dplyr::filter(layer == max(layer)) %>%
                            dplyr::mutate(tautot = cumsum(tau)) %>%
                            dplyr::mutate(refftop = reffl[1],
                                          refftop.pre = reffl.pre[1]) %>%  ## check this with Kenta
                            dplyr::ungroup() %>%
                            dplyr::mutate(time = time[i])
                    }, .parallel = FALSE, .progress = "none")
                })
            }, .parallel = TRUE) -> df
        saveRDS(df, sprintf("%scfodd-%s.rds", out.prefix, experiment))
    }

#' @export
process.pr.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             experiment = "amip-rain-15", out.prefix = "",
             years = 1979:1983,
             ncores = 12,
             flux = TRUE,
             subsample = NULL) {
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname <- sprintf("%s/%s/%s_%d%02d.01_rain2d.nc",
                                     datadir, experiment, experiment, year, month)
                    nc <- try(ncdf4::nc_open(fname), silent = TRUE)
                    if (class(nc) == "try-error")
                        return(NULL)
                    t <- ncdf4::ncvar_get(nc, "time")
                    lon <- ncdf4::ncvar_get(nc, "lon")
                    lat <- ncdf4::ncvar_get(nc, "lat")
                    pr <- ncdf4::ncvar_get(nc, "aprl_na")
                    prc <- ncdf4::ncvar_get(nc, "aprc_na")
                    ncdf4::nc_close(nc)
                    df <- expand.grid(lon = as.vector(lon),
                                      lat = as.vector(lat),
                                      time = as.vector(t)) %>%
                        cbind(pr = as.vector(pr),
                              prc = as.vector(prc))
                    ## if subsampling is requested, do it here
                    (if (is.null(subsample))
                         df
                     else dplyr::filter(df, time %in% t[seq(1, length(t), by = subsample)])) %>%
                        plyr::ddply(~ lon + lat, function(x) {
                            rbind(x %>%
                                  dplyr::mutate(pr.class = cut(3600 * pr,
                                                               c(0, 0.01, 0.1, 1, 1e8), right = FALSE)) %>%
                                  dplyr::group_by(pr.class) %>%
                                  dplyr::summarize(n = n(),
                                                   sum.pr = sum(pr),
                                                   pr.type = "prl"),
                                  x %>%
                                  dplyr::mutate(pr.class = cut(3600 * prc,
                                                               c(0, 0.01, 0.1, 1, 1e8), right = FALSE)) %>%
                                  dplyr::group_by(pr.class) %>%
                                  dplyr::summarize(n = n(),
                                                   sum.pr = sum(prc),
                                                   pr.type = "prc"),
                                  x %>%
                                  dplyr::mutate(pr.class = cut(3600 * (pr + prc),
                                                               c(0, 0.01, 0.1, 1, 1e8), right = FALSE)) %>%
                                  dplyr::group_by(pr.class) %>%
                                  dplyr::summarize(n = n(),
                                                   sum.pr = sum(pr + prc),
                                                   pr.type = "prl + prc"))
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
            }, .progress = "none", .parallel = TRUE) -> df
        saveRDS(df, sprintf("%spr-hist-%s.rds", out.prefix, experiment))
    }

#' @export
process.rad.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             experiment = "amip-rain-15", out.prefix = "",
             years = 1979:1983,
             ncores = 12,
             flux = TRUE,
             subsample = NULL) { ## monthly data, so subsampling is ignored
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                with(x, {
                    fname <- sprintf("%s/%s/%s_%d%02d.01_echamm.nc",
                                     datadir, experiment, experiment, year, month)
                    nc <- try(ncdf4::nc_open(fname), silent = TRUE)
                    if (class(nc) == "try-error") { ## try echam instead of echamm
                        fname <- sprintf("%s/%s/%s_%d%02d.01_echam.nc",
                                         datadir, experiment, experiment, year, month)
                        nc <- try(ncdf4::nc_open(fname), silent = TRUE)
                        if (class(nc) == "try-error")
                            return(NULL)
                    }
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
                    aclcov <- ncdf4::ncvar_get(nc, "aclcov") ## cloud cover
                    xivi   <- ncdf4::ncvar_get(nc, "xivi") ## IWP
                    xlvi   <- ncdf4::ncvar_get(nc, "xlvi") ## LWP
                    qvi    <- ncdf4::ncvar_get(nc, "qvi") ## water vapor path
                    ncdf4::nc_close(nc)
                    df <- expand.grid(lon = as.vector(lon),
                                      lat = as.vector(lat),
                                      time = as.vector(t)) %>%
                        cbind(srad0d = as.vector(srad0d ),
                              srad0  = as.vector(srad0  ),
                              trad0  = as.vector(trad0  ),
                              sraf0  = as.vector(sraf0  ),
                              traf0  = as.vector(traf0  ),
                              aprl  = as.vector(aprl  ),
                              aprc  = as.vector(aprc  ),
                              aprs  = as.vector(aprs  ),
                              aclcov = as.vector(aclcov),
                              xivi   = as.vector(xivi  ),
                              xlvi   = as.vector(xlvi  ),
                              qvi    = as.vector(qvi   ))
                    ## ## if subsampling is requested, do it here
                    ## ifelse(is.null(subsample),
                    ##        df,
                    ##        dplyr::slice(df, seq(1, nrow(df), subsample)))
                })
            }, .progress = "none", .parallel = FALSE) -> df
        saveRDS(df, sprintf("%srad-%s.rds", out.prefix, experiment))
    }

#' @export
process.forcing.echam <-
    function(datadir = "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments",
             exp_pd = "rain_4_1_-1_no-cosp_nudged", exp_pi = "rain_4_1_-1_pi_no-cosp_nudged",
             out.prefix = "",
             years = 1979:1983,
             ncores = 24,
             flux = TRUE,
             subsample = NULL) { ## monthly data, so subsampling is ignored
        doParallel::registerDoParallel(cores = ncores)
        expand.grid(year = years, month = 1:12) %>%
            plyr::ddply(~ year + month, function(x) {
                gc()
                ## first try to find cached data frame
                fname.cache <- sprintf("%sforcing-%s-%d-%d.rds", out.prefix, exp_pd, x$year, x$month)
                df2 <- try(readRDS(fname.cache))
                ## if it's not cached, generate it
                if (class(df2) == "try-error") {
                    plyr::ddply(cbind(x, pi_pd = c("PI", "PD"), exp = c(exp_pi, exp_pd)),
                                ~ pi_pd,
                                function(x) {
                                    with(x, {
                                        fname <- sprintf("%s/%s/%s_%d%02d.01_forcing.nc",
                                                         datadir, exp, exp, year, month)
                                        nc <- try(ncdf4::nc_open(fname), silent = TRUE)
                                        if (class(nc) == "try-error")
                                            return(NULL)
                                        t <- ncdf4::ncvar_get(nc, "time")
                                        lon <- ncdf4::ncvar_get(nc, "lon")
                                        lat <- ncdf4::ncvar_get(nc, "lat")
                                        xlvi <- ncdf4::ncvar_get(nc, "XLVI")
                                        fsw_diff <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_LWP")
                                        fsw_total_top_unpert <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_UNPERT")
                                        flw_total_top_unpert <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_UNPERT")
                                        fsw_total_top_lwp    <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_LWP")
                                        flw_total_top_lwp    <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_LWP")
                                        fsw_total_top_cdnc   <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_CDNC")
                                        flw_total_top_cdnc   <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_CDNC")
                                        ## fsw_total_top_cldfra   <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_CLDFRA")
                                        ## flw_total_top_cldfra   <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_CLDFRA")
                                        xlvi       <- ncdf4::ncvar_get(nc, "XLVI")
                                        cdnc       <- ncdf4::ncvar_get(nc, "CDNC")
                                        cldfra     <- ncdf4::ncvar_get(nc, "CLDFRA")
                                        cldfra_liq <- ncdf4::ncvar_get(nc, "CLDFRA_LIQ")
                                        ncdf4::nc_close(nc)
                                        ## finite difference approximation to logarithmic derivatives
                                        dlog.fsw.dlog.lwp <-
                                            2 * (fsw_total_top_lwp - fsw_total_top_unpert) /
                                            (fsw_total_top_lwp + fsw_total_top_unpert) /
                                            (0.1 / 1.05)
                                        dlog.flw.dlog.lwp <-
                                            2 * (flw_total_top_lwp - flw_total_top_unpert) /
                                            (flw_total_top_lwp + flw_total_top_unpert) /
                                            (0.1 / 1.05)
                                        dlog.fsw.dlog.cdnc <-
                                            2 * (fsw_total_top_cdnc - fsw_total_top_unpert) /
                                            (fsw_total_top_cdnc + fsw_total_top_unpert) /
                                            (0.1 / 1.05)
                                        dlog.flw.dlog.cdnc <-
                                            2 * (flw_total_top_cdnc - flw_total_top_unpert) /
                                            (flw_total_top_cdnc + flw_total_top_unpert) /
                                            (0.1 / 1.05)
                                        expand.grid(lon = as.vector(lon),
                                                    lat = as.vector(lat),
                                                    time = as.vector(t)) %>%
                                            dplyr::mutate(dlog.fsw.dlog.lwp  = as.vector(dlog.fsw.dlog.lwp ),
                                                          dlog.flw.dlog.lwp  = as.vector(dlog.flw.dlog.lwp ),
                                                          dlog.fsw.dlog.cdnc = as.vector(dlog.fsw.dlog.cdnc),
                                                          dlog.flw.dlog.cdnc = as.vector(dlog.flw.dlog.cdnc),
                                                          fsw_total_top_unpert = as.vector(fsw_total_top_unpert),
                                                          flw_total_top_unpert = as.vector(flw_total_top_unpert),
                                                          xlvi               = as.vector(xlvi              ),
                                                          cdnc               = as.vector(cdnc              ),
                                                          cldfra             = as.vector(cldfra            ),
                                                          cldfra_liq         = as.vector(cldfra_liq        ))
                                    })
                                }) -> df
                    ## if (0) {
                    df.sens <- df %>%
                        dplyr::select(-(xlvi : cldfra_liq)) %>%
                        tidyr::gather(var, val, dlog.fsw.dlog.lwp : flw_total_top_unpert) %>%
                        dplyr::group_by(lon, lat, time, var) %>%
                        ## calculate mean between PD and PI for each combination of grouping variables
                        dplyr::summarize(val = mean(val)) %>%
                        tidyr::spread(var, val) %>%
                        dplyr::ungroup() 
                    
                    df.pert <- df %>%
                        dplyr::select(-(dlog.fsw.dlog.lwp : flw_total_top_unpert)) %>%
                        tidyr::gather(var, val, xlvi : cldfra_liq) %>%
                        dplyr::group_by(lon, lat, time, var) %>%
                        tidyr::spread(pi_pd, val) %>%
                        dplyr::summarize(dlog = 2 * (PD - PI) / (PD + PI)) %>%
                        tidyr::spread(var, dlog) %>%
                        dplyr::ungroup()

                    df2 <- dplyr::full_join(df.sens, df.pert,
                                            by = c("lon", "lat", "time"))

                    try(saveRDS(df2, fname.cache))
                }
                
                df2 %<>%
                    dplyr::mutate(lon = ifelse(lon <= 180, lon, lon - 360)) %>%
                    dplyr::group_by(lon, lat) %>%
                    dplyr::summarize(
                               ## note that any NAs in the following calculation
                               ## indicate an error in the method
                               sw_lwp =  mean(ifelse(dlog.fsw.dlog.lwp == 0, 0,  dlog.fsw.dlog.lwp * xlvi * fsw_total_top_unpert), na.rm = FALSE),
                               sw_cdnc = mean(ifelse(dlog.fsw.dlog.cdnc == 0, 0, dlog.fsw.dlog.cdnc * cdnc * fsw_total_top_unpert), na.rm = FALSE),
                               lw_lwp =  mean(ifelse(dlog.flw.dlog.lwp == 0, 0,  dlog.flw.dlog.lwp * xlvi * flw_total_top_unpert), na.rm = FALSE),
                               lw_cdnc = mean(ifelse(dlog.flw.dlog.cdnc == 0, 0, dlog.flw.dlog.cdnc * cdnc * flw_total_top_unpert), na.rm = FALSE),
                               sw_lwp_mean =  mean(dlog.fsw.dlog.lwp  * fsw_total_top_unpert) * mean(ifelse(dlog.fsw.dlog.lwp == 0, 0,  replace(xlvi, !is.finite(xlvi), 0))),
                               sw_cdnc_mean = mean(dlog.fsw.dlog.cdnc * fsw_total_top_unpert) * mean(ifelse(dlog.fsw.dlog.cdnc == 0, 0, replace(cdnc, !is.finite(cdnc), 0))),
                               lw_lwp_mean =  mean(dlog.flw.dlog.lwp  * flw_total_top_unpert) * mean(ifelse(dlog.flw.dlog.lwp == 0, 0,  replace(xlvi, !is.finite(xlvi), 0))),
                               lw_cdnc_mean = mean(dlog.flw.dlog.cdnc * flw_total_top_unpert) * mean(ifelse(dlog.flw.dlog.cdnc == 0, 0, replace(cdnc, !is.finite(cdnc), 0))))

                if (0) {
                    df2 %>%
                        tidyr::gather(pert, val, sw_lwp : lw_cdnc) %>%
                        ggplot( aes(lon, lat, fill = val)) +
                        geom_raster() +
                        scale_x_geo(facet = TRUE) + scale_y_geo() +
                        coord_fixed(xlim = c(-180, 180), ylim = c(-80, 80), expand = FALSE) +
                        ## scale_x_continuous("", labels = NULL, breaks = NULL) +
                        ## scale_y_continuous("", labels = NULL, breaks = NULL) +
                        ## scale_fill_manual(values = col.frac, name = expression(f[liq])) +
                        ## scale_fill_warmfrac() +
                        ## scale_fill_brewer("PD$-$PI", palette = "RdBu", drop = FALSE, direction = -1) +
                        scale_fill_distiller(#labels = tikz_sanitize,
                            "PD$-$PI", palette = "RdYlBu"## , 
                            ## limits = c(-1, 1) * 10
                        ) +
                        geom_world_polygon(highres = FALSE) +
                        theme_bw()  +
                        theme(legend.position = "bottom", legend.box = "horizontal") +
                        guides(fill = guide_colorbar(direction = "horizontal", title.vjust = 1)) +
                        facet_grid(pert ~ month)

                    df2 %>%
                        tidyr::gather(pert, val, sw_lwp : lw_cdnc) %>%
                        group_by(lon, lat, pert) %>%
                        summarize(val = mean(val)) %>%
                        ggplot( aes(lon, lat, fill = val)) +
                        geom_raster() +
                        scale_x_geo(facet = TRUE) + scale_y_geo() +
                        coord_fixed(xlim = c(-180, 180), ylim = c(-80, 80), expand = FALSE) +
                        ## scale_x_continuous("", labels = NULL, breaks = NULL) +
                        ## scale_y_continuous("", labels = NULL, breaks = NULL) +
                        ## scale_fill_manual(values = col.frac, name = expression(f[liq])) +
                        ## scale_fill_warmfrac() +
                        ## scale_fill_brewer("PD$-$PI", palette = "RdBu", drop = FALSE, direction = -1) +
                        scale_fill_distiller(#labels = tikz_sanitize,
                            "PD$-$PI", palette = "RdYlBu", 
                            limits = c(-1, 1) * 10
                        ) +
                        geom_world_polygon(highres = FALSE) +
                        theme_bw()  +
                        theme(legend.position = "bottom", legend.box = "horizontal") +
                        guides(fill = guide_colorbar(direction = "horizontal", title.vjust = 1)) +
                        facet_wrap(~ pert)

                    df2 %>%
                        tidyr::gather(pert, val, sw_lwp : lw_cdnc) %>%
                        mutate(cos.lat = cos(lat * pi / 180)) %>%
                        group_by(pert) %>%
                        summarize(flux = sum(val * cos.lat, na.rm = FALSE) / sum(cos.lat))

                }

                df2
            }, .progress = "none", .parallel = TRUE) -> df
        saveRDS(df, sprintf("%sforcing-%s.rds", out.prefix, exp_pd))
    }
