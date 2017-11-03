
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

    df <- plyr::adply(1:100, 1, .id = "subcol", function(subcol) {
        gc()
        fname.dbze <- sprintf("%s_200405.01_cosp_%03d.nc", experiment, subcol)
        nc.dbze <-    ncdf4::nc_open(paste(path, fname.dbze  , sep = ""))
        dbze    <- ncdf4::ncvar_get(nc.dbze  , sprintf("dbze94_%03d" , subcol ))
        fracout <- ncdf4::ncvar_get(nc.dbze  , sprintf("frac_out_%03d", subcol ))
        ncdf4::nc_close(nc.dbze)
        
        expand.grid(lon = as.vector(lon),
                    lat = as.vector(lat),
                    lev = as.vector(lev),
                    time = as.vector(time)) %>%
            dplyr::mutate(lssnow  = as.vector(lssnow),
                          lsrain  = as.vector(lsrain),
                          tau     = as.vector(tau   ),
                          dbze    = as.vector(dbze  ),
                          fracout = as.vector(fracout),
                          reffl   = as.vector(reffl ),
                          reffi   = as.vector(reffi )) %>%
            dplyr::filter(fracout == 1)
    }, .progress = "text")
}

    
