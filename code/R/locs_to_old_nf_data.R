library(dplyr)
library(here)
library(sf)

LCMAP_PROJ <- readRDS(here("data/AEA_WGS84.rds"))

crosswalk <- read.csv("/Volumes/big_bag/data/CUI/crosswalk_11122020_CUI.csv")
true_locations <- read.csv(
    "/Volumes/big_bag/data/CUI/Beier_ESF_NYS_ActualCoordinates_AVERAGED_201021.csv"
) |> select(PLOT, COUNTYCD, UNITCD, STATECD, LAT_AVG_ACTUAL, LON_AVG_ACTUAL) |>
    na.omit() |>
    group_by(PLOT, COUNTYCD, UNITCD, STATECD) |>
    summarize(
        LAT_AVG_ACTUAL = mean(LAT_AVG_ACTUAL, na.rm = T),
        LON_AVG_ACTUAL = mean(LON_AVG_ACTUAL, na.rm = T)
    ) |>
    ungroup() |>
    right_join(crosswalk, by = c("PLOT", "COUNTYCD", "UNITCD", "STATECD"))

add_locs_to_plots <- function(plots, true_locations) {
    plots <- left_join(
        plots, 
        select(true_locations, PLOT_PUBLIC, LAT_AVG_ACTUAL, LON_AVG_ACTUAL),
        by = c("PLOT" = "PLOT_PUBLIC")
    ) |>
        st_as_sf(
            coords = c("LON_AVG_ACTUAL", "LAT_AVG_ACTUAL"),
            crs = 4326,
            agr = "constant"
        ) |>
        st_transform(LCMAP_PROJ) %>%
        mutate(x = sf::st_coordinates(.)[, 1], y = sf::st_coordinates(.)[, 2]) |>
        as.data.frame() |>
        select(-geometry)
}
    
    
training_plots <- read.csv(here("data/training.csv")) |>
    add_locs_to_plots(true_locations) |>
    select(PLOT, year, fia, ensemble, x, y) |>
    write.csv(
        "/Volumes/big_bag/data/CUI/ERE_621/lucas_stuff/training.csv",
        row.names = F
    )

testing_plots <- read.csv(here("data/testing.csv")) |>
    add_locs_to_plots(true_locations) |>
    select(PLOT, year, fia, ensemble, x, y) |>
    write.csv(
        "/Volumes/big_bag/data/CUI/ERE_621/lucas_stuff/testing.csv", 
        row.names = F
    )