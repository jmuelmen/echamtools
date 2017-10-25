## ---- cfodd ---------------------
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
