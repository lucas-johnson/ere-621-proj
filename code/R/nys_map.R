library(labrador.client)
library(ggplot2)
library(ggspatial)
library(sf)
library(here)
library(dplyr)

nys_plot <- st_read(here("data/nys_shape/Counties_Shoreline.shp")) |>
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