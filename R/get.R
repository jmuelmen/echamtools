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
                         cautalpha = NA,
                         cautbeta = NA,
                         amip = FALSE,
                         pi = FALSE,
                         path = "/home/jmuelmen/wcrain/echam-ham") {
    plyr::ddply(expand.grid(ccraut = ccrauts,
                            ccauloc = ccaulocs,
                            creth = creth,
                            cautalpha = cautalpha,
                            cautbeta = cautbeta),
                ~ ccraut + ccauloc + creth + cautalpha + cautbeta,
                function(df)
                    with(df, {
                        experiment <- expname(ccraut, ccauloc, creth,
                                              cautalpha, cautbeta, amip, pi)
                        readRDS(sprintf("%s/pr-hist-%s.rds", path, experiment)) %>%
                            filter_whole_years(exact = FALSE) %>%
                            mutate(pi_pd = ifelse(pi, "PI", "PD"))
                    }))
}

#' @export
get.cosp.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60),
                         ccaulocs = NA,
                         creth = NA,
                         cautalpha = NA,
                         cautbeta = NA,
                         amip = FALSE,
                         pi = FALSE,
                         path = "/home/jmuelmen/wcrain/echam-ham") {
    plyr::ddply(expand.grid(ccraut = ccrauts,
                            ccauloc = ccaulocs,
                            creth = creth,
                            cautalpha = cautalpha,
                            cautbeta = cautbeta),
                ~ ccraut + ccauloc + creth + cautalpha + cautbeta,
                function(df)
                    with(df, {
                        experiment <- expname(ccraut, ccauloc, creth,
                                              cautalpha, cautbeta, amip, pi)
                        readRDS(sprintf("%s/cosp-%s.rds", path, experiment)) %>%
                            filter_whole_years(exact = TRUE) %>%
                            dplyr::mutate(ccraut = ccraut,
                                          ccauloc = ccauloc,
                                          creth = creth,
                                          cautalpha = cautalpha,
                                          cautbeta = cautbeta,
                                          pi_pd = ifelse(pi, "PI", "PD"))
                    }))
}

#' @export
get.rad.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60),
                          ccaulocs = NA,
                          creth = NA,
                          cautalpha = NA,
                          cautbeta = NA,
                          amip = FALSE,
                          pi = FALSE,
                          path = "/home/jmuelmen/wcrain/echam-ham") {
    plyr::ddply(expand.grid(ccraut = ccrauts,
                            ccauloc = ccaulocs,
                            creth = creth,
                            cautalpha = cautalpha,
                            cautbeta = cautbeta),
                ~ ccraut + ccauloc + creth + cautalpha + cautbeta,
                function(df)
                    with(df, {
                        experiment <- expname(ccraut, ccauloc, creth,
                                              cautalpha, cautbeta, amip, pi)
                        df <- try(readRDS(sprintf("%s/rad-%s.rds", path, experiment)))
                        if (any(class(df) == "try-error"))
                            df <- try(readRDS(sprintf("%s/rad-%s_no-cosp.rds", path, experiment)))
                        df %>%
                            filter_whole_years() %>%
                            dplyr::mutate(ccraut = ccraut,
                                  ccauloc = ccauloc,
                                  creth = creth,
                                  cautalpha = cautalpha,
                                  cautbeta = cautbeta,
                                  pi_pd = ifelse(pi, "PI", "PD"))
                    }))
}

#' @export
get.mask.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60),
                           ccaulocs = NA,
                           creth = NA,
                           cautalpha = NA,
                           cautbeta = NA,
                           amip = FALSE,
                           pi = FALSE,
                           flux = TRUE,
                           path = "/home/jmuelmen/wcrain/echam-ham") {
    plyr::ddply(expand.grid(ccraut = ccrauts,
                            ccauloc = ccaulocs,
                            creth = creth,
                            cautalpha = cautalpha,
                            cautbeta = cautbeta),
                ~ ccraut + ccauloc + creth + cautalpha + cautbeta,
                function(df)
                    with(df, {
                        experiment <- expname(ccraut, ccauloc, creth,
                                              cautalpha, cautbeta, amip, pi)
                        ## print(experiment)
                        readRDS(sprintf("%s/%s.rds", path, experiment)) %>%
                            filter_whole_years() %>%
                            dplyr::mutate(ccraut = ccraut,
                                  ccauloc = ccauloc,
                                  creth = creth,
                                  cautalpha = cautalpha,
                                  cautbeta = cautbeta,
                                  pi_pd = ifelse(pi, "PI", "PD"))
                }))
}

filter_whole_years <- function(df, exact = TRUE) {
    plyr::ddply(df, ~ year, function(df) {
        tbl <- table(df$month)
        if (exact) {
            if (length(tbl) == 12 && all(tbl == tbl[1])) 
                df
            else
                NULL
        } else {
            if (length(tbl) != 12) {
                warning(sprintf("filter_whole_years: only %d month(s) in year %d",
                              length(tbl), df$year[1]))
                NULL
            } else if (any(abs(tbl - mean(tbl)) / mean(tbl) > 0.1)) {
                warning(sprintf("filter_whole_years: max deviation %.2f%% in %2d/%d",
                              100 * max(abs(tbl - mean(tbl)) / mean(tbl)),
                              which.max(abs(tbl - mean(tbl)) / mean(tbl)), df$year[1]))
                NULL
            }
            else {
                ## warning(sprintf("Max deviation %.2f%% in %2d/%d",
                ##               100 * max(abs(tbl - mean(tbl)) / mean(tbl)),
                ##               which.max(abs(tbl - mean(tbl)) / mean(tbl)), df$year[1]))
                df
            }
        }
    })
}

expname <- function(ccraut, ccauloc, creth,
                    cautalpha, cautbeta, amip, pi) {
    experiment <- sprintf("%s%g%s%s%s%s%s",
                          ifelse(amip, "amip-rain-", "rain_"),
                          ccraut,
                          ifelse(is.na(ccauloc), "", sprintf("_%g", ccauloc)),
                          ifelse(is.na(creth), "", sprintf("_%g", creth)),
                          ifelse(is.na(cautalpha), "", sprintf("_cautalpha_%g", cautalpha)),
                          ifelse(is.na(cautbeta), "", sprintf("_cautbeta_%g", cautbeta)),
                          ifelse(pi, "_pi", ""))
}
