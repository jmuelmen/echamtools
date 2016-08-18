get.phase.at.locations.echam <- function(locations) {
    df <- readRDS("cfsites/echam-warm-frac-tbl.rds") 
    df %>% group_by(site) %>% summarize(warm.frac = (function(mask) {
        tbl <- table(mask)
        tbl["warm"] / (tbl["warm"] + tbl["cold"])
    })(mask)) %>%
        rename(id = site) %>% right_join(locations) %>%
        mutate(center = "MPI-M", model = "MPI-ESM-LR", ens = "r1i1p1") %>%
        rename(station = name, liquid = warm.frac)
}

#' @export
get.pr.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60)) {
    ldply(sprintf("amip-rain-%g", ccrauts),
          function(experiment) {
              readRDS(sprintf("pr-hist-%s.rds", experiment)) %>%
                  cbind(experiment = experiment)
          })
}

#' @export
get.rad.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60)) {
    ldply(sprintf("amip-rain-%g", ccrauts),
          function(experiment) {
              readRDS(sprintf("rad-%s.rds", experiment)) %>%
                  cbind(experiment = experiment)
          })
}

#' @export
get.mask.echam <- function(ccrauts = c(0, 0.01, 0.1, 1, 2, 3.75, 7.5, 15, 30, 60)) {
    ldply(sprintf("amip-rain-%g", ccrauts),
          function(experiment) {
              readRDS(sprintf("%s.rds", experiment)) %>%
                  cbind(experiment = experiment)
          })
}
