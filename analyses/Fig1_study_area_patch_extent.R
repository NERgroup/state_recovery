
#jogsmith@ucsc.edu

#data required 
#1. Benthic survey locations
#2. Site patch polygons
#3. LDA-predicted patch types
#4. California county boundaries

rm(list=ls())

################################################################################
#Step 0: set paths and load data

require(librarian)
librarian::shelf(tidyverse, lubridate, sf,
                 rnaturalearth, rnaturalearthdata, ggspatial)

datadir <- "/Volumes/enhydra/data/kelp_recovery"

#load benthic survey data
load(file.path(datadir, "MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda"))

#load site patches
site_patches <- st_read(here::here("output","gis_data","processed","site_patch_polygons.shp"))

#load LDA-predicted patch types
load(here::here("output","lda_patch_transitions.rda")) 


################################################################################
#Step 1: Retain survey locations shared by kelp and quadrat data
quad_zone <- kelp_data %>%
  dplyr::distinct(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::inner_join(
    quad_data %>%
      dplyr::distinct(site, site_type, latitude, longitude, zone, survey_date),
    by = c("site","site_type","latitude","longitude","zone","survey_date")
  ) %>%
  dplyr::arrange(latitude, longitude, site, site_type, survey_date, zone)


################################################################################
#Step 2: assign model-predicted patch types

quad_zone_with_pred <- quad_zone %>%
  left_join(
    transitions_tbl_constrained %>%
      dplyr::select(site, site_type, zone, patch_2024, patch_2025),
    by = c("site","site_type", "zone")
  ) %>%
  mutate(
    pred_patch = case_when(
      format(survey_date, "%Y") == "2024" ~ as.character(patch_2024),
      format(survey_date, "%Y") == "2025" ~ as.character(patch_2025),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-patch_2024, -patch_2025)


################################################################################
#Step 3: join with patch geometry

# Convert to sf using lat/lon
quad_zone_sf <- st_as_sf(
  quad_zone_with_pred,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)


site_patches_single <- site_patches %>%
  st_cast("POLYGON") %>%              # split MULTIPOLYGON into indiv POLYGON
  mutate(patch_id = row_number())     # assign unique polygon ID

#join points to polygons
site_patches_with_points <- site_patches_single %>%
  st_join(quad_zone_sf, join = st_intersects, left = TRUE)


#tidy up and add classifier
quad_build3 <- site_patches_with_points %>%
                dplyr::select(-site_type.x) %>%
                dplyr::select(patch_id, latitude, longitude, 
                       survey_date,site, site_type = site_type.y, pred_patch, 
                       everything()) %>%
                mutate(pred_patch = ifelse(is.na(pred_patch),site_type,pred_patch)) %>%
                filter(!(is.na(pred_patch)))

################################################################################
#Build Figure1

my_theme <- theme(axis.text=element_text(size=6),
                  axis.text.y = element_text(angle = 90, hjust = 0.5),
                  axis.title=element_text(size=8),
                  plot.tag=element_text(size=8),
                  plot.title=element_text(size=7, face="bold"),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key = element_blank(),
                  legend.background = element_rect(fill=alpha('blue', 0)),
                  legend.key.height = unit(1, "lines"), 
                  legend.text = element_text(size = 6),
                  legend.title = element_text(size = 7),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 6 ,face="bold"),
                  panel.background = element_rect(fill="white"))


# Inset map (California with MB bounding box)
usa     <- ne_states(country="United States of America", returnclass="sf")
foreign <- ne_countries(country=c("Canada","Mexico"), returnclass="sf")
#read state
ca_counties <- st_read(file.path(datadir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 


# Build inset
g1_inset <-  ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data=foreign, fill="grey80", color="white", lwd=0.3) +
    geom_sf(data=usa, fill="grey80", color="white", lwd=0.3) +
    # Plot box
    annotate("rect", xmin=-122.6, xmax=-121, ymin=36.2, ymax=37.1, color="black", fill=NA, lwd=0.8) +
    # Labels
    labs(x="", y="") +
    # Crop
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Theme
    theme_bw() + my_theme +
    theme( plot.margin = unit(rep(0, 4), "null"),
           panel.margin = unit(rep(0, 4), "null"),
           panel.background = element_rect(fill='transparent'), #transparent panel bg
           axis.ticks = element_blank(),
           axis.ticks.length = unit(0, "null"),
           axis.ticks.margin = unit(0, "null"),
           axis.text = element_blank(),
           axis.title=element_blank(),
           axis.text.y = element_blank(),
           plot.title =  element_text(hjust = 0, size=14, face = "bold"))
)



g1 <-ggplot(quad_build3 %>% filter(year(survey_date) == 2024)) +
  # plot patches
  geom_sf(aes(fill = pred_patch), color = "black") +
  # add survey sites (subtidal sites)
  geom_sf(
    data = st_centroid(quad_build3 %>% filter(year(survey_date) == 2024)),
    aes(shape = "Subtidal survey \nsite"),
    fill = "yellow",
    color = "black",
    size = 1.5
  ) +
  # plot land
  geom_sf(data = ca_counties, color = "grey70", fill = "grey85") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.53, 36.64)) +
  theme_bw() +
  # Add plot inset
  annotation_custom(
    grob = g1_inset,
    xmin = -122.001,
    xmax = -121.97,
    ymin = 36.62
  ) +
  labs(
    title = "2024",
    fill = "Patch type",
    shape = "",  # legend title for shapes,
    tag = "A"
  ) +
  scale_fill_manual(
    name   = "Patch type",
    values = c(
      "BAR"   = "#7570B3",
      "FOR"   = "#1B9E77",
      "INCIP" = "#D95F02"
    ),
    breaks = c("BAR","FOR","INCIP"),
    labels = c("Barren","Forest","Incipient")
  )+
  # custom shape legend without otters
  scale_shape_manual(
    values = c(
      "Subtidal survey \nsite"        = 24
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 1, override.aes = list(
      fill = c("yellow"),
      color = c( "black"),
      alpha = c(0.6),
      size = c(2)
    ))
  ) +
  # north arrow (upper right)
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(0.5, "cm"),
    width  = unit(0.5, "cm"),
    style  = north_arrow_orienteering(text_col = NA)
  )+
  # scale bar (lower right)
  annotation_scale(
    location = "br",
    width_hint = 0.25,
    text_cex = 0.8
  ) +
  my_theme 


g2 <-ggplot(quad_build3 %>% filter(year(survey_date) == 2025)) +
  # plot patches
  geom_sf(aes(fill = pred_patch), color = "black") +
  # add survey sites (subtidal sites)
  geom_sf(
    data = st_centroid(quad_build3 %>% filter(year(survey_date) == 2025)),
    aes(shape = "Subtidal survey \nsite"),
    fill = "yellow",
    color = "black",
    size = 1.5
  ) +
  # plot land
  geom_sf(data = ca_counties, color = "grey70", fill = "grey85") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.53, 36.64)) +
  theme_bw() +
  labs(
    title = "2025",
    fill = "Patch type",
    shape = "",  # legend title for shapes
    tag = "B"
  ) +
  scale_fill_manual(
    name   = "Patch type",
    values = c(
      "BAR"   = "#7570B3",
      "FOR"   = "#1B9E77",
      "INCIP" = "#D95F02"
    ),
    breaks = c("BAR","FOR","INCIP"),
    labels = c("Barren","Forest","Incipient")
  )+
  # custom shape legend without otters
  scale_shape_manual(
    values = c(
      "Subtidal survey \nsite"        = 24
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 1, override.aes = list(
      fill = c("yellow"),
      color = c( "black"),
      alpha = c(0.6),
      size = c(2)
    ))
  ) +
  my_theme 

g <- ggpubr::ggarrange(g1, g2, common.legend = TRUE, legend = "right")
g

ggsave(
   here::here("figures", "Fig1_map_figure.png"),
   g,
   width = 7, height = 4, dpi = 600, bg = "white"
 )
