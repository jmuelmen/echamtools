## get.phase.at.locations.echam <- function(locations) {
##     df <- readRDS("cfsites/echam-warm-frac-tbl.rds") 
##     df %>% group_by(site) %>% summarize(warm.frac = (function(mask) {
##         tbl <- table(mask)
##         tbl["warm"] / (tbl["warm"] + tbl["cold"])
##     })(mask)) %>%
##         rename(id = site) %>% right_join(locations) %>%
##         mutate(center = "MPI-M", model = "MPI-ESM-LR", ens = "r1i1p1") %>%
##         rename(station = name, liquid = warm.frac)
## }

#' @export
get.pr.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60),
                         ccaulocs = NA,
                         creth = NA,
                         amip = TRUE,
                         pi = TRUE,
                         path = "/home/jmuelmen/wcrain/echam-ham") {
    plyr::ddply(expand.grid(ccraut = ccrauts,
                            ccauloc = ccaulocs,
                            creth = creth),
                ~ ccraut + ccauloc + creth,
                function(df)
                    with(df, {
                        experiment <- sprintf("%s%g%s%s%s",
                                              ifelse(amip, "amip-rain-", "rain_"),
                                              ccraut,
                                              ifelse(is.na(ccauloc), "", sprintf("_%g", ccauloc)),
                                              ifelse(is.na(creth), "", sprintf("_%g", creth)),
                                              ifelse(pi, "_pi", ""))
                        readRDS(sprintf("%s/pr-hist-%s.rds", path, experiment)) %>%
                            filter_whole_years() %>%
                            cbind(ccraut = ccraut,
                                  ccauloc = ccauloc,
                                  creth = creth)
                    }))
}

#' @export
get.rad.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60),
                          ccaulocs = NA,
                          creth = NA,
                          amip = TRUE,
                          pi = FALSE,
                          path = "%s") {
    plyr::ddply(expand.grid(ccraut = ccrauts,
                            ccauloc = ccaulocs,
                            creth = creth),
                ~ ccraut + ccauloc + creth,
                function(df) 
                    with(df, {
                        experiment <- sprintf("%s%g%s%s%s",
                                              ifelse(amip, "amip-rain-", "rain_"),
                                              ccraut,
                                              ifelse(is.na(ccauloc), "", sprintf("_%g", ccauloc)),
                                              ifelse(is.na(creth), "", sprintf("_%g", creth)),
                                              ifelse(pi, "_pi", ""))
                        readRDS(sprintf("%s/rad-%s.rds", path, experiment)) %>%
                            filter_whole_years() %>%
                            cbind(ccraut = ccraut,
                                  ccauloc = ccauloc,
                                  creth = creth)
                    }))
}

#' @export
get.mask.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60),
                           ccaulocs = NA,
                           flux = TRUE,
                           creth = NA,
                           amip = TRUE,
                           pi = FALSE,
                           path = "%s") {
    plyr::ddply(expand.grid(ccraut = ccrauts,
                            ccauloc = ccaulocs,
                            creth = creth),
                ~ ccraut + ccauloc + creth,
                function(df)
                    with(df, {
                        experiment <- sprintf("%s%g%s%s",
                                              ifelse(amip, "amip-rain-", "rain_"),
                                              ccraut,
                                              ifelse(is.na(ccauloc), "", sprintf("_%g", ccauloc)),
                                              ifelse(is.na(creth), "", sprintf("_%g", creth)),
                                              ifelse(pi, "_pi", ""),
                                              ifelse(flux, "", "-mr"))
                        readRDS(sprintf("%s/%s.rds", path, experiment)) %>%
                            filter_whole_years() %>%
                            cbind(ccraut = ccraut,
                                  ccauloc = ccauloc,
                                  creth = creth)
                }))
}

filter_whole_years <- function(df) {
    plyr::ddply(df, ~ year, function(df) {
        tbl <- table(df$month)
        if (length(tbl) == 12 && all(tbl == tbl[1]))
            df
        else
            NULL
    })
}
