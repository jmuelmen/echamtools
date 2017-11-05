
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
    fname.lsrain <- sprintf("%s_%s", experiment, "200405.01_cosp_lsrain.nc"      )
    fname.lssnow <- sprintf("%s_%s", experiment, "200405.01_cosp_lssnow.nc"      )
    fname.tau    <- sprintf("%s_%s", experiment, "200405.01_cosp_cisccp_tau3d.nc")
    fname.reffl  <- sprintf("%s_%s", experiment, "200405.01_cosp_reffl.nc"       )
    fname.reffi  <- sprintf("%s_%s", experiment, "200405.01_cosp_reffi.nc"       )

    nc.lssnow <- ncdf4::nc_open(paste(path, fname.lssnow, sep = ""))
    nc.lsrain <- ncdf4::nc_open(paste(path, fname.lsrain, sep = ""))
    nc.tau <-    ncdf4::nc_open(paste(path, fname.tau   , sep = ""))
    nc.reffl <-  ncdf4::nc_open(paste(path, fname.reffl , sep = ""))
    nc.reffi <-  ncdf4::nc_open(paste(path, fname.reffi , sep = ""))

    lon  <- ncdf4::ncvar_get(nc.lssnow, "lon")
    lat  <- ncdf4::ncvar_get(nc.lssnow, "lat")
    lev  <- ncdf4::ncvar_get(nc.lssnow, "mlev")
    time  <- ncdf4::ncvar_get(nc.lssnow, "time")

    lssnow  <- ncdf4::ncvar_get(nc.lssnow, "lssnow")
    lsrain  <- ncdf4::ncvar_get(nc.lsrain, "lsrain")
    tau     <- ncdf4::ncvar_get(nc.tau   , "cisccp_tau3d"   )
    reffl   <- ncdf4::ncvar_get(nc.reffl , "reffl"  )
    reffi   <- ncdf4::ncvar_get(nc.reffi , "reffi"  )

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

    hist(cldfrac[cldfrac > 0])
    
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
                          dbze    = as.vector(dbze  ),
                          fracout = as.vector(fracout)) %>%
            dplyr::filter(fracout == 1)  %>%
            ## snow or rain, not both
            dplyr::filter(lssnow < 1e-8 | lsrain < 1e-8)
        

        df %<>% dplyr::mutate(dbze = replace(dbze, dbze < -1e29, NA))
        df %<>% tidyr::gather(type, q_precip, lssnow.ic : lsrain) %>%
            dplyr::mutate(type = factor(type))
        ## df %<>% discretize(q_precip, c(0, 10 ^ seq(-7, -2, 0.2)))
        ## df %<>% discretize(dbze, c(-1000, seq(-50, -20, 10), seq(-15, 10, 5), 100))
        ## df %<>% discretize(q_precip, 30, as_factor = TRUE, equal_contents = TRUE)
        ## df %<>% discretize(dbze, 30, as_factor = TRUE, equal_contents = TRUE)
        ## df %<>% plotutils::discretize(q_precip, 30, as_factor = FALSE, equal_contents = TRUE)
        df %<>% dplyr::mutate(q_precip = log10(q_precip)) %>%
            plotutils::discretize(q_precip, seq(-8.75,-2.25,0.5), as_factor = FALSE) %>%
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
            scale_color_brewer(palette = "Set1") +
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

    
