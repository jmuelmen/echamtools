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

    
