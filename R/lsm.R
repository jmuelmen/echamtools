#' Generate T63 land/sea mask
#' @export
generate_lsm <- function(fname = "inst/T63GR15_jan_surf.nc") {
    nc <- ncdf4::nc_open(fname)
    lon  <- ncdf4::ncvar_get(nc, "lon")
    lat  <- ncdf4::ncvar_get(nc, "lat")
    df <- expand.grid(lon = as.vector(lon),
                      lat = as.vector(lat)) %>%
        dplyr::mutate(lsm = as.vector(ncdf4::ncvar_get(nc, "SLM"))) %>%
        dplyr::mutate(lsm = factor(lsm, levels = 0:1, labels = c("Sea", "Land")))
    ncdf4::nc_close(nc)
    df
}
