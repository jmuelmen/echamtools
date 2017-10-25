
#' @export 
cfodd.process <- function(ccraut, ccauloc, creth, pi = FALSE) {
    experiment <- sprintf("%s%g%s%s%s",
                          "rain_",
                          ccraut,
                          ifelse(is.na(ccauloc), "", sprintf("_%g", ccauloc)),
                          ifelse(is.na(creth), "", sprintf("_%g", creth)),
                          ifelse(pi, "_pi", ""))
    print(sprintf("Experiment: %s", experiment))

    path <- sprintf("/work/bb0839/b380126/%s/", experiment)
    fname.lsrain <- sprintf("%s_%s", experiment, "200405.01_cosp_lsrain.nc"      )
    fname.lssnow <- sprintf("%s_%s", experiment, "200405.01_cosp_lssnow.nc"      )
    fname.tau    <- sprintf("%s_%s", experiment, "200405.01_cosp_cisccp_tau3d.nc")
    fname.dbze   <- sprintf("%s_%s", experiment, "200405.01_cosp_037.nc"         )
    fname.reffl  <- sprintf("%s_%s", experiment, "200405.01_cosp_reffl.nc"       )
    fname.reffi  <- sprintf("%s_%s", experiment, "200405.01_cosp_reffi.nc"       )

    nc.lssnow <- ncdf4::nc_open(paste(path, fname.lssnow, sep = ""))
    nc.lsrain <- ncdf4::nc_open(paste(path, fname.lsrain, sep = ""))
    nc.tau <-    ncdf4::nc_open(paste(path, fname.tau   , sep = ""))
    nc.dbze <-   ncdf4::nc_open(paste(path, fname.dbze  , sep = ""))
    nc.reffl <-  ncdf4::nc_open(paste(path, fname.reffl , sep = ""))
    nc.reffi <-  ncdf4::nc_open(paste(path, fname.reffi , sep = ""))

    lon  <- ncdf4::ncvar_get(nc.lssnow, "lon")
    lat  <- ncdf4::ncvar_get(nc.lssnow, "lat")
    lev  <- ncdf4::ncvar_get(nc.lssnow, "mlev")
    time  <- ncdf4::ncvar_get(nc.lssnow, "time")

    lssnow  <- ncdf4::ncvar_get(nc.lssnow, "lssnow")
    lsrain  <- ncdf4::ncvar_get(nc.lsrain, "lsrain")
    tau     <- ncdf4::ncvar_get(nc.tau   , "cisccp_tau3d"   )
    dbze    <- ncdf4::ncvar_get(nc.dbze  , "dbze94_037"  )
    fracout <- ncdf4::ncvar_get(nc.dbze  , "frac_out_037"  )
    reffl   <- ncdf4::ncvar_get(nc.reffl , "reffl"  )
    reffi   <- ncdf4::ncvar_get(nc.reffi , "reffi"  )

    mask <- (apply(reffi, c(1,2,4), function(x) rep(all(x == 0 | x == 4), 31)) &
             apply(fracout, c(1,2,4), function(x) rep(any(x == 1), 31)))  %>%
        aperm(c(2,3,1,4))

    label.vertical.features <- function(vfm) {
        x <- vfm
        if (length(x) == 0)
            return(x)
        diff.x <- diff(c(-1, x)) ## guarantee that the first group of 1's is preceded by a transition
        labels <- cumsum(diff.x != 0 & x != 0) * (x != 0) ## count up the edges
        labels
    }

    expand.grid(lon = as.vector(lon),
                lat = as.vector(lat),
                lev = as.vector(lev),
                time = as.vector(time))  %>%
        dplyr::filter(as.vector(mask)) %>%
        dplyr::mutate(lssnow  = as.vector(lssnow  [mask]),
                      lsrain  = as.vector(lsrain  [mask]),
                      tau     = as.vector(tau     [mask]),
                      dbze    = as.vector(dbze    [mask]),
                      fracout = as.vector(fracout [mask]),
                      reffl   = as.vector(reffl   [mask]),
                      reffi   = as.vector(reffi   [mask])) -> df

    df %<>%
        dplyr::group_by(lon, lat, time) %>%
        dplyr::mutate(layer = label.vertical.features(fracout)) %>%
        dplyr::filter(layer == max(layer)) %>%
        dplyr::mutate(tautot = cumsum(tau)) %>%
        dplyr::mutate(refftop = reffl[1]) %>%  ## check this
        dplyr::ungroup() %>%
        dplyr::filter(dbze > -100) %>%
        dplyr::mutate(dbze = round(dbze),
                      tautot = round(tautot)) 

    saveRDS(df, sprintf("%s.rds", experiment))
}

#' @export 
cfodd.plot <- function(df) {
    df %>%
        dplyr::filter(dbze > -30, tautot < 60) %>%
        plotutils::discretize(refftop, seq(5, 20, 5), as_factor = TRUE) %>%
        dplyr::filter(!is.na(refftop)) %>%
        plotutils::discretize(dbze, seq(-30, 20, by = 2)) %>%
        plotutils::discretize(tautot, seq(0, 60, by = 2)) %>%
        dplyr::group_by(dbze, tautot, refftop, ccraut, creth, ccauloc) %>%
        dplyr::summarize(count = n()) %>%
        dplyr::group_by(tautot, refftop, ccraut, creth, ccauloc) %>%
        dplyr::mutate(rel.count = count / sum(as.numeric(count))) %>%
        dplyr::ungroup()
}

bind_rows(readRDS("200405.rds") %>%
          mutate(ccraut = 4, creth = -1, ccauloc = 1),
          ## readRDS("rain_0.0001.rds") %>%
          ## mutate(ccraut = 1e-4, creth = -1, ccauloc = 1),
          readRDS("rain_4_1_15.rds") %>%
          mutate(ccraut = 4, creth = 15, ccauloc = 1),
          readRDS("rain_4_1_17.rds") %>%
          mutate(ccraut = 4, creth = 17, ccauloc = 1)) %>%
    cfodd.plot() %>%
    ggplot(aes(x = dbze, y = tautot)) +
    geom_raster(aes(fill = (rel.count))) +
    facet_grid(creth ~ refftop) +
    scale_y_reverse() +
    scale_fill_distiller(palette = "Blues", direction = -1) +
    labs(x = "Reflectivity (dB$Z_e$)", y = "Optical depth") +
    guides(fill = "none") +
    theme_bw(24)
    
