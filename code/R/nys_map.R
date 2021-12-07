library(labrador.client)
library(ggplot2)
library(ggspatial)
library(sf)
library(here)
library(dplyr)

field_plots <- rbind(
        st_read("/Volumes/big_bag/data/CUI/ERE_621/lucas_stuff/training.csv"),
        st_read("/Volumes/big_bag/data/CUI/ERE_621/lucas_stuff/testing.csv")
    ) |>
    st_as_sf(coords = c("x", "y"), crs = readRDS("data/AEA_WGS84.rds")) |>
    filter(year == 2019)
    
counties_shoreline <- st_read(here("data/nys_shape/Counties_Shoreline.shp"))

nys_plot <- counties_shoreline |>
    ggplot() +
    geom_sf(fill = "gray", color = "black") +
    theme_void() +
    annotation_north_arrow(
        location = "tr",
        style = north_arrow_orienteering(
            fill = c("black", "black")
        ),
        height = unit(0.3, "in"),
        width = unit(0.3, "in")
    ) + 
    annotation_scale(location = "bl", text_cex = 1)

ggsave(
    here("figures/nys_county_map.png"),
    nys_plot,
    width = 5,
    height = 5,
    units = "in"
)

field_plots <- st_transform(field_plots, st_crs(counties_shoreline))

fia_dist <- nys_plot + geom_sf(data = field_plots, color = "black")

ggsave(
    here("figures/fia_spat_dist.png"),
    fia_dist, 
    width = 5, 
    height = 5, 
    units = "in"
)
