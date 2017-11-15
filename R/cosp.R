
#' @export 
cosp.process <- function(ccraut,
                         ccauloc,
                         creth,
                         pi = FALSE,
                         path = "/work/bb0839/b380126") {
    experiment <- sprintf("%s%g%s%s%s",
                          "rain_",
                          ccraut,
                          ifelse(is.na(ccauloc), "", sprintf("_%g", ccauloc)),
                          ifelse(is.na(creth), "", sprintf("_%g", creth)),
                          ifelse(pi, "_pi", ""))
    print(sprintf("Experiment: %s", experiment))

    path <- sprintf("%s/%s/", path, experiment)
    fname.lsrain   <- sprintf("%s_%s", experiment, "200405.01_cosp_lsrain.nc"      )
    fname.lssnow   <- sprintf("%s_%s", experiment, "200405.01_cosp_lssnow.nc"      )
    fname.tau      <- sprintf("%s_%s", experiment, "200405.01_cosp_cisccp_tau3d.nc")
    fname.reffl    <- sprintf("%s_%s", experiment, "200405.01_cosp_reffl.nc"       )
    fname.reffi    <- sprintf("%s_%s", experiment, "200405.01_cosp_reffi.nc"       )
    fname.aclc     <- sprintf("%s_%s", experiment, "200405.01_cosp_aclc.nc"       )
    fname.tm1      <- sprintf("%s_%s", experiment, "200405.01_cosp_tm1.nc"       )
    fname.tm1_cosp <- sprintf("%s_%s", experiment, "200405.01_cosp_tm1_cosp.nc"       )
    fname.xl       <- sprintf("%s_%s", experiment, "200405.01_cosp_xl.nc"       )
    fname.xi       <- sprintf("%s_%s", experiment, "200405.01_cosp_xi.nc"       )
    fname.rain3d   <- sprintf("%s_%s", experiment, "200405.01_rain3d.nc"       )
    fname.rain2d   <- sprintf("%s_%s", experiment, "200405.01_rain2d.nc"       )

    nc.lssnow <- ncdf4::nc_open(paste(path, fname.lssnow, sep = ""))
    nc.lsrain <- ncdf4::nc_open(paste(path, fname.lsrain, sep = ""))
    nc.tau <-    ncdf4::nc_open(paste(path, fname.tau   , sep = ""))
    nc.reffl <-  ncdf4::nc_open(paste(path, fname.reffl , sep = ""))
    nc.reffi <-  ncdf4::nc_open(paste(path, fname.reffi , sep = ""))
    nc.aclc  <-  ncdf4::nc_open(paste(path, fname.aclc  , sep = ""))
    nc.tm1      <-  ncdf4::nc_open(paste(path, fname.tm1     , sep = ""))
    nc.tm1_cosp <-  ncdf4::nc_open(paste(path, fname.tm1_cosp, sep = ""))
    nc.xl       <-  ncdf4::nc_open(paste(path, fname.xl      , sep = ""))
    nc.xi       <-  ncdf4::nc_open(paste(path, fname.xi      , sep = ""))
    nc.rain3d <- ncdf4::nc_open(paste(path, fname.rain3d, sep = ""))
    nc.rain2d <- ncdf4::nc_open(paste(path, fname.rain2d, sep = ""))

    lon  <- ncdf4::ncvar_get(nc.lssnow, "lon")
    lat  <- ncdf4::ncvar_get(nc.lssnow, "lat")
    lev  <- ncdf4::ncvar_get(nc.lssnow, "mlev")
    time  <- ncdf4::ncvar_get(nc.lssnow, "time")

    lssnow  <- ncdf4::ncvar_get(nc.lssnow, "lssnow")
    lsrain  <- ncdf4::ncvar_get(nc.lsrain, "lsrain")
    tau     <- ncdf4::ncvar_get(nc.tau   , "cisccp_tau3d"   )
    reffl   <- ncdf4::ncvar_get(nc.reffl , "reffl"  )
    reffi   <- ncdf4::ncvar_get(nc.reffi , "reffi"  )
    aclc    <- ncdf4::ncvar_get(nc.aclc  , "aclc"  )
    tm1     <- ncdf4::ncvar_get(nc.tm1     , "tm1" )
    tm1_cosp<- ncdf4::ncvar_get(nc.tm1_cosp, "tm1_cosp"  )
    xl      <- ncdf4::ncvar_get(nc.xl      , "xl" )
    xi      <- ncdf4::ncvar_get(nc.xi      , "xi"  )
    xrl     <- ncdf4::ncvar_get(nc.rain3d, "xrl_na"  )
    xsl     <- ncdf4::ncvar_get(nc.rain3d, "xsl_na"  )
    qaut    <- ncdf4::ncvar_get(nc.rain3d, "autoconv_na"  )
    aprlv   <- ncdf4::ncvar_get(nc.rain3d, "aprlv_na")
    aprsv   <- ncdf4::ncvar_get(nc.rain3d, "aprsv_na")
    aprl    <- ncdf4::ncvar_get(nc.rain2d, "aprl_na") %>%
        apply(1:3, rep, 31) %>%
        aperm(c(2,3,1,4)) 
    aprs    <- ncdf4::ncvar_get(nc.rain2d, "aprs_na") %>%
        apply(1:3, rep, 31) %>%
        aperm(c(2,3,1,4)) 
    
    label.vertical.features <- function(vfm) {
        x <- vfm
        if (length(x) == 0)
            return(x)
        diff.x <- diff(c(-1, x)) ## guarantee that the first group of 1's is preceded by a transition
        labels <- cumsum(diff.x != 0 & x != 0) * (x != 0) ## count up the edges
        labels
    }

    layer <- apply(aclc > 0, c(1,2,4), label.vertical.features) %>%
        aperm(c(2,3,1,4))
        
    ## somewhat laborious cloud fraction calculation
    cldfrac <- array(0, dim(lssnow))
    plyr::a_ply(1:100, 1, function(subcol) {
        gc()
        fname.dbze <- sprintf("%s_200405.01_cosp_%03d.nc", experiment, subcol)
        nc.dbze <-    ncdf4::nc_open(paste(path, fname.dbze  , sep = ""))
        fracout <- ncdf4::ncvar_get(nc.dbze  , sprintf("frac_out_%03d", subcol ))
        ncdf4::nc_close(nc.dbze)
        cldfrac <<- cldfrac + fracout
    }, .progress = "text")
    cldfrac <- 1e-2 * cldfrac

    expand.grid(lon = as.vector(lon),
                lat = as.vector(lat),
                lev = as.vector(lev),
                time = as.vector(time)) %>%
        dplyr::mutate(cosp.cldfrac = as.vector(cldfrac),
                      aclc = as.vector(aclc)) %>%
        saveRDS(sprintf("cosp-cldfrac-%s.rds", experiment))

    ## hist(cldfrac[cldfrac > 0])
    
    df <- plyr::adply(1:100, 1, .id = "subcol", function(subcol) {
        gc()
        fname.dbze <- sprintf("%s_200405.01_cosp_%03d.nc", experiment, subcol)
        nc.dbze <-    ncdf4::nc_open(paste(path, fname.dbze  , sep = ""))
        dbze    <- ncdf4::ncvar_get(nc.dbze  , sprintf("dbze94_%03d" , subcol ))
        fracout <- ncdf4::ncvar_get(nc.dbze  , sprintf("frac_out_%03d", subcol ))
        ncdf4::nc_close(nc.dbze)
        
        df <- expand.grid(lon = as.vector(lon),
                          lat = as.vector(lat),
                          lev = as.vector(lev),
                          time = as.vector(time)) %>%
            dplyr::mutate(lssnow  = as.vector(lssnow),
                          lsrain  = as.vector(lsrain),
                          xrl     = as.vector(xrl  ),
                          xsl     = as.vector(xsl  ),
                          qaut    = as.vector(qaut ),
                          aprlv   = as.vector(aprlv),
                          aprsv   = as.vector(aprsv),
                          dbze    = as.vector(dbze  ),
                          aclc    = as.vector(aclc  ),
                          tm1      = as.vector(tm1     ),
                          tm1_cosp = as.vector(tm1_cosp),
                          xl       = as.vector(xl      ),
                          xi       = as.vector(xi      ),
                          ## cldfrac = as.vector(cldfrac ),
                          aprl    = as.vector(aprl),
                          aprs    = as.vector(aprs),
                          layer = as.vector(layer),
                          fracout = as.vector(fracout))

        rm(lssnow  )
        rm(lsrain  )
        rm(xrl     )
        rm(xsl     )
        rm(qaut    )
        rm(aprlv   )
        rm(aprsv   )
        rm(dbze    )
        rm(aclc    )
        rm(tm1     )
        rm(tm1_cosp)
        rm(xl      )
        rm(xi      )
        rm(aprl    )
        rm(aprs    )
        rm(layer   )
        rm(fracout )
        gc()
        
        if (subcol == 1) {
            df %>%
                dplyr::filter(fracout == 1)  %>%
                ## snow or rain, not both
                dplyr::filter(lssnow < 1e-8 | lsrain < 1e-8) %>%
                dplyr::filter(dbze > -50) %>%
                saveRDS(sprintf("cosp-teaser-%s.rds", experiment))
        }

        if (subcol == 1) {
            df %>%
                dplyr::select(aprlv, aprsv) %>%
                tidyr::gather(type, flux, aprlv : aprsv) %>%
                dplyr::mutate(type = factor(type)) %>%
                dplyr::group_by(type) %>%
                dplyr::mutate(frac.0 = mean(flux == 0)) %>%
                dplyr::ungroup() %>%
                dplyr::filter(flux != 0) %>%
                saveRDS(sprintf("cosp-aprlv-aprsv-%s.rds", experiment))
        }

        if (subcol == 1) {
            df %>%
                dplyr::group_by(lon, lat, time) %>%
                dplyr::summarize(aprlv.max = max(aprlv),
                                 aprsv.max = max(aprsv)) %>%
                dplyr::ungroup() %>%
                tidyr::gather(type, flux, aprlv.max : aprsv.max) %>%
                dplyr::mutate(type = factor(type)) %>%
                dplyr::group_by(type) %>%
                dplyr::mutate(frac.0 = mean(flux == 0)) %>%
                dplyr::ungroup() %>%
                saveRDS(sprintf("cosp-aprlv.max-aprsv.max-%s.rds", experiment))
        }
        
        df <- readRDS(sprintf("~/cosp-teaser-%s.rds", experiment))

        df %>%
            group_by(lon, lat, time) %>%
            summarize(dbze.max = max(dbze),
                      aprlv.max = max(aprlv),
                      aprsv.max = max(aprsv),
                      aprlv.max.dbze = aprlv[which.max(aprlv)],
                      aprl = aprl[1],
                      aprs = aprs[1],
                      temp.highest.drizzle = ifelse(any(dbze > -15),
                                                    tm1_cosp[min(lev[dbze > -15])],
                                                    NA),
                      temp.highest.rain = ifelse(any(dbze > 0),
                                                 tm1_cosp[min(lev[dbze > 0])],
                                                 NA),
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
            saveRDS(sprintf("cosp2d-%s.rds", experiment))


        df %<>% dplyr::mutate(dbze = replace(dbze, dbze < -1e29, NA))
        df %<>% tidyr::gather(type, q_precip, lssnow : aprsv) %>%
            dplyr::mutate(type = factor(type))
        ## df %<>% discretize(q_precip, c(0, 10 ^ seq(-7, -2, 0.2)))
        ## df %<>% discretize(dbze, c(-1000, seq(-50, -20, 10), seq(-15, 10, 5), 100))
        ## df %<>% discretize(q_precip, 30, as_factor = TRUE, equal_contents = TRUE)
        ## df %<>% discretize(dbze, 30, as_factor = TRUE, equal_contents = TRUE)
        ## df %<>% plotutils::discretize(q_precip, 30, as_factor = FALSE, equal_contents = TRUE)
        df %<>% dplyr::mutate(q_precip = log10(q_precip)) %>%
            plotutils::discretize(q_precip, seq(-6.75,-2.25,0.5), as_factor = FALSE) %>%
                dplyr::filter(!is.na(q_precip))

        df %>%
            dplyr::group_by(type, q_precip) %>%
            dplyr::filter(dbze > -50) %>%
            ## dplyr::summarize(sd = sd(dbze), dbze = mean(dbze)) %>%
            ggplot(aes(x = factor(q_precip), y = dbze, fill = type)) +
            #geom_line() +
            #geom_pointrange(aes(ymin = dbze - 1.96 * sd, ymax = dbze + 1.96 * sd )) +
            geom_violin() +
            ## scale_x_log10() +
            scale_fill_brewer(palette = "Set1") +
            ## coord_cartesian(ylim = c(-50, 10)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

        df %>%
            dplyr::group_by(type, q_precip) %>%
            dplyr::summarize(sd = sd(dbze), dbze = mean(dbze)) %>%
            ggplot(aes(x = q_precip, y = dbze, col = type)) +
            geom_line() +
            geom_pointrange(aes(ymin = dbze - 1.96 * sd, ymax = dbze + 1.96 * sd )) +
            ## geom_violin() +
            scale_x_log10() +
            scale_color_brewer(palette = "Set1") +
            ## coord_cartesian(ylim = c(-50, 10)) +
            theme_bw()





        df %<>% discretize(dbze, 30, as_factor = FALSE, equal_contents = TRUE)
        ## df %<>%
        ##     dplyr::group_by(lon, lat, lev, type, q_precip, dbze) %>%
        ##     dplyr::summarize(n = as.numeric(n()))

        df %<>%
            dplyr::group_by(type, q_precip, dbze) %>%
            dplyr::summarize(n = as.numeric(n()))

        df %>%
            dplyr::group_by(type, q_precip) %>%
            dplyr::summarize(dbze = weighted.mean(dbze, n)) %>%
            ggplot(aes(x = q_precip, y = dbze, col = type)) +
            geom_line() +
            geom_point() +
            scale_x_log10() +
            scale_color_brewer(palette = "Set1") +
            coord_cartesian(ylim = c(-50, 10)) +
            theme_bw()
        
        df %>%
            group_by(type, q_precip) %>%
            mutate(n = n / sum(n)) %>%
            ggplot(aes(x = q_precip, y = dbze, fill = n)) +
            scale_fill_distiller(palette = "Spectral") +
            geom_tile() +
            facet_wrap(~ type, nrow = 2) +
            theme_bw()

        df %>%
            group_by(type, q_precip) %>%
            mutate(n = n / sum(n)) %>%
            ggplot(aes(x = q_precip, y = dbze, fill = n)) +
            scale_fill_distiller(palette = "Spectral") +
            geom_tile() +
            facet_wrap(~ type, nrow = 2) +
            theme_bw()

        df
    }, .progress = "text")
}

    
