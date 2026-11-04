# ============================================================
# Tetrapodomorpha Publication History Analysis (Species & Genus)
# ============================================================

setwd("C:/Users/mlove/Desktop/MSc-Paleobio/semester-3/Research/Data sets")

## Packages
library(stringr)
library(rnaturalearth)
library(sf)
library(countrycode)
library(dplyr)
library(ggplot2)
library(scales)
library(rgplates)
library(chronosphere)
library(iNEXT)

## Load data
stages <- read.csv("2026_stages.csv")  # redated to analyze 2026 stages vs 2026 data

species_data <- read.csv("pbdb_data_species_overlap_v3.csv", header = TRUE, sep = ",", skip = 22)
genus_data <- read.csv("pbdb_data_genus_overlap_v3.csv", header = TRUE, sep = ",", skip = 22)

species_data$ref_year <- as.integer(str_extract(species_data$primary_reference, "\\d{4}"))
genus_data$ref_year <- as.integer(str_extract(genus_data$primary_reference, "\\d{4}"))

## ============================================================
## Accumulation curves through publication history
## ============================================================
build_accumulation_curve <- function(subset_data) {
  subset_data <- subset_data[!is.na(subset_data$max_ma) & !is.na(subset_data$min_ma), ]   # Remove rows that lack dating

  subset_data[subset_data$max_ma == subset_data$min_ma, "min_ma"] <- subset_data[subset_data$max_ma == subset_data$min_ma, "min_ma"] - 1e-6  # avoid zero-width ranges
    min_end     <- outer(subset_data$max_ma, stages$bottom, pmin)  # youngest max age per occurrence subset_data interval
    max_start   <- outer(subset_data$min_ma, stages$top, pmax)     # oldest min age per occurrence subset_data interval
    overlap_mat <- min_end - max_start                    # overlap = the difference
    best_bin    <- max.col(overlap_mat, ties.method = "first")
    subset_data$stage <- stages$stage[best_bin]
    subset_data$bin   <- stages$stg[best_bin]

  subset_data <- subset_data %>% arrange(ref_year)
  subset_data <- subset_data %>%
    mutate(
      cum_occurrences = row_number(),
      new_taxon = !duplicated(accepted_name),   # TRUE the first time a taxon appears
      cum_taxa = cumsum(new_taxon)
    )
  subset_data %>%
    group_by(ref_year) %>%
    summarise(
      cum_occurrences = max(cum_occurrences),
      cum_taxa = max(cum_taxa),
      .groups = "drop"
    )

  subset_data
}

publication_analysis <- function(data, resolution_label) {
  only_US   <- data[data$cc == "US", ]
  only_ZA   <- data[data$cc == "ZA", ]
  no_US     <- data[data$cc != "US" | is.na(data$cc), ]   # global without US (ZA kept)
  no_ZA     <- data[data$cc != "ZA" | is.na(data$cc), ]   # global without ZA (US kept)

  full_curve   <- build_accumulation_curve(data)    %>% mutate(method = "Full Global")
  onlyUS_curve <- build_accumulation_curve(only_US) %>% mutate(method = "USA Exclusive")
  onlyZA_curve <- build_accumulation_curve(only_ZA) %>% mutate(method = "South Africa Exclusive")
  noUS_curve   <- build_accumulation_curve(no_US)   %>% mutate(method = "Global without USA")
  noZA_curve   <- build_accumulation_curve(no_ZA)   %>% mutate(method = "Global without South Africa")

  line_data <- bind_rows(full_curve, onlyUS_curve, onlyZA_curve, noUS_curve, noZA_curve)
  line_data$method <- factor(line_data$method,
    levels = c("USA Exclusive", "South Africa Exclusive", "Full Global",
               "Global without USA", "Global without South Africa"))

  p <- ggplot(line_data, aes(x = ref_year, y = cum_taxa, color = method, linetype = method)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(
      values = c("Full Global" = "black",
                 "South Africa Exclusive" = "#D55E00",
                 "USA Exclusive" = "#0072B2",
                 "Global without USA" = "gray50",
                 "Global without South Africa" = "gray70"),
      name = NULL) +
    scale_linetype_manual(values = c(1, 1, 1, 4, 2), name = NULL) +
    labs(y = paste0("Cumulative unique ", resolution_label, " recorded"), x = "Year",
         title = paste0("Cumulative Unique ", resolution_label, " Through Publication History")) +
    theme_classic(base_size = 12) +
    theme(legend.position = c(0.2, 0.85), legend.background = element_blank(),
          plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

  list(plot = p, data = line_data,
       only_US = onlyUS_curve, only_ZA = onlyZA_curve,
       no_US = noUS_curve, no_ZA = noZA_curve)
}


species_result <- publication_analysis(species_data, "Species")
genus_result <- publication_analysis(genus_data, "Genus")
genus_result$plot  


## ============================================================
## Publications-by-country summary table
## ============================================================

build_publication_summary_table <- function(data) {
  data %>%
    group_by(cc) %>%
    summarise(
      first_publication = min(ref_year, na.rm = TRUE),
      last_publication  = max(ref_year, na.rm = TRUE),
      total_publications = n_distinct(primary_reference),
      total_occurrences  = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(total_publications)) %>%
    mutate(iso_a3 = countrycode(cc, origin = "iso2c", destination = "iso3c",
                                 custom_match = c("UK" = "GBR")))
}

## ============================================================
## Modern-day choropleth maps
## ============================================================

world <- ne_countries(scale = "medium", returnclass = "sf")

build_publication_summary <- function(data) {
  publication_summary <- build_publication_summary_table(data)

  map_data <- world %>%
    left_join(publication_summary, by = "iso_a3")

  number_map <- ggplot(map_data) +
    geom_sf(aes(fill = total_publications), color = "gray99", linewidth = 0.1) +
    scale_fill_viridis_c(name = "Total\npublications", na.value = "gray85", option = "cividis") +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90),expand = FALSE) +
    labs(x = "lng", y = "lat") +
    theme(panel.background = element_rect(fill = "lightblue2"),panel.grid.major = element_line(
      color = "gray90",linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "right")
  
  first_pub_map <- ggplot(map_data) +
    geom_sf(aes(fill = first_publication), color = "gray99", linewidth = 0.1) +
    scale_fill_viridis_c(name = "First\npublication\nyear", na.value = "gray85", option = "cividis") +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90),expand = FALSE) +
    labs(x = "lng", y = "lat") +
    theme(panel.background = element_rect(fill = "lightblue2"),panel.grid.major = element_line(
      color = "gray90",linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "right")
    
  list(summary = publication_summary, number_map = number_map, first_pub_map = first_pub_map)
}

species_maps <- build_publication_summary(species_data)
genus_maps   <- build_publication_summary(genus_data)

species_maps$number_map
species_maps$first_pub_map

genus_maps$number_map
genus_maps$first_pub_map


## ============================================================
## Carboniferous occurrence point map
## ============================================================

build_kasimovian_map <- function(data, lon_col = "lng", lat_col = "lat") {

  carb_data <- data %>%
    filter(max_ma < 307 & min_ma > 303.07) %>%          
    filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))
 
   carb_localities <- carb_data %>%
    group_by(collection_no, .data[[lon_col]], .data[[lat_col]]) %>%
    summarise(n_occurrences = n(), .groups = "drop")

  # Convert to sf points
  carb_sf <- st_as_sf(
    carb_localities,
    coords = c(lon_col, lat_col),
    crs = 4326,
    remove = FALSE
  )

  ggplot() +
    geom_sf(data = world, fill = "gray95", color = "gray80", linewidth = 0.1) +
    geom_sf(data = carb_sf, aes(fill = n_occurrences),
            shape = 21, color = "gray40", stroke = 0.4, size = 2.2, alpha = 0.9) +
    scale_fill_distiller(
      name = "Occurrences\nper locality",
      palette = "Reds",
      direction = 1,
      trans = "log1p"
    ) +
    labs(title = "Carboniferous-Early Permian (358.9-272.9my) sample localities") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}

kas_map <- build_kasimovian_map(species_data)
kas_map

## ============================================================
## Paleogeographic reconstruction map
## ============================================================

# Fetch a full-Phanerozoic plate rotation model -- used only for the
# background coastline outline, since paleolng/paleolat are already
# present on each occurrence (per-occurrence reconstruction from PBDB)
paleo_model <- fetch(src = "GPlates", ser = "MERDITH2021")

build_paleogeographic_map <- function(data, coastline_age = NULL, model = paleo_model, fill_limits = NULL) {
  if (!is.null(coastline_age)) {
    data <- data %>% filter((min_ma + max_ma)/2 >= coastline_age - 2.5, (min_ma + max_ma)/2 <= coastline_age + 2.5)
  }
  publication_summary <- data %>%
    group_by(collection_no) %>%
    summarise(
      total_publications = n_distinct(primary_reference),
      total_occurrences  = n_distinct(occurrence_no),
      mean_paleolng = mean(paleolng, na.rm = TRUE),
      mean_paleolat = mean(paleolat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_paleolng), !is.na(mean_paleolat))

  publication_summary <- sort_by(publication_summary, publication_summary$total_occurrences)
  p <- ggplot()
  if (!is.null(coastline_age)) {
    paleocoast <- reconstruct("coastlines", age = coastline_age, model = model)
    p <- p + geom_sf(data = paleocoast, fill = "gray85", color = "gray50", linewidth = 0.2)
  }
  p <- p +
    geom_point(data = publication_summary, shape = 21, size = 4,
           aes(x = mean_paleolng, y = mean_paleolat, fill = total_occurrences),
           color = "black", alpha = 0.85) +
    scale_fill_viridis_b(breaks = c(1, 5, 10, 20, 30, 40),
                      limits = fill_limits,
                      oob = scales::squish,
                      name = "Number of Occurrences") +
    coord_sf() +
    labs(title = paste0("Tetrapodomorpha occurrences by collection, ", coastline_age, "mya")) +
    theme_minimal() +
theme(axis.title = element_blank(),      
      axis.text = element_text(size = 8), 
      axis.ticks = element_line(color = "gray40", linewidth = 0.2),
      panel.grid = element_blank()) 
  p
}

global_summary <- species_data %>%
  group_by(collection_no) %>%
  summarise(total_occurrences = n_distinct(occurrence_no), .groups = "drop")

fill_limits <- range(global_summary$total_occurrences, na.rm = TRUE)

for (i in seq(400, 250, by = -5)) {
  png(paste0("Paleomap_5myr/Occurrences_", i, ".png"), width = 900, height = 600)
  
  print(build_paleogeographic_map(species_data, coastline_age = i, fill_limits = fill_limits))
  
  dev.off()
}


## =======================================
## Sample Completeness through Time
## =======================================

initial_completeness <- function(data){
  data <- data[!is.na(data$ref_year),]
  data <- data[!is.na(data$accepted_name),]
  data <- data[!is.na(data$collection_no),]

  completeness_through_time <- c()
    
    for (i in min(data$ref_year):max(data$ref_year)) {
      ## call each age interval
      target_age <- data[data$ref_year <= i,]
      
      ## Transform into incidence data
      incidence_data <- target_age[, c("collection_no", "accepted_name")] %>% table() %>% t()
      incidence_raw <- as.matrix(incidence_data)
      incidence_raw <- (incidence_raw > 0) * 1L 
      
      if (ncol(incidence_raw) < 2 || nrow(incidence_raw) < 1) {
        completeness_through_time <- append(completeness_through_time, NA)
        next
      }
      
      info <- DataInfo(incidence_raw, "incidence_raw")
      
      completeness_through_time <- append(completeness_through_time, info$SC)
    }
    names(completeness_through_time) <- min(data$ref_year):max(data$ref_year)
    completeness_through_time
}

completeness_analysis <- function(data, resolution_label){
  
  no_US <- data[!data$cc == "US", ]
  no_ZA <- data[!data$cc == "ZA", ]
  
  only_US <-  data[data$cc == "US", ]
  only_ZA <- data[data$cc == "ZA", ]
  
  rest_global <- data[!data$cc %in% c("ZA", "US"),]
  
  ## calculate completeness for all country subsets
  no_US_completeness <- initial_completeness(no_US)
  no_ZA_completeness <- initial_completeness(no_ZA)
  
  only_US_completeness <- initial_completeness(only_US)
  only_ZA_completeness <- initial_completeness(only_ZA)
  
  global_completeness <- initial_completeness(data)
  no_USZA_completeness <- initial_completeness(rest_global)
  
  ##return
  list(
    no_US_completeness = no_US_completeness,
    no_ZA_completeness = no_ZA_completeness,
    only_US_completeness = only_US_completeness,
    only_ZA_completeness = only_ZA_completeness,
    global_completeness = global_completeness,
    no_USZA_completeness = no_USZA_completeness
      )
}

species_completeness <- completeness_analysis(species_data, "Species")

## plotting
plot(as.numeric(names(species_completeness$global_completeness)), 
  species_completeness$global_completeness, ylab = "Estimated Sample Coverage",
  xlab = "Year", ylim=c(0,1), type="l")
lines(as.numeric(names(species_completeness$only_US_completeness)), 
  species_completeness$only_US_completeness, col="red")
lines(as.numeric(names(species_completeness$only_ZA_completeness)), 
  species_completeness$only_ZA_completeness, col="blue")


## Background graphs
ZA_pub_per_year <- tapply(only_ZA$primary_reference, only_ZA$ref_year, function(x) length(unique(x)))
  ZA_occ_per_year <- tapply(only_ZA$occurrence_no, only_ZA$ref_year, function(x) length(unique(x)))
  
   # Plot the first line
   plot(as.numeric(names(ZA_occ_per_year)), ZA_occ_per_year,
        type = "l", col = "blue", lwd = 2,
        xlab = "Reference year", ylab = "Count",
        ylim = range(c(ZA_occ_per_year, ZA_pub_per_year), na.rm = TRUE))
   
   # Add the second line to the same plot
   lines(as.numeric(names(ZA_pub_per_year)), ZA_pub_per_year,
         type = "l", col = "orange", lwd = 2)
   
   # Add a legend
   legend("topleft", legend = c("Unique occurrences", "Unique publications"),
          col = c("blue", "orange"), lwd = 2)
   # Set up margin space on the right for the second axis
   par(mar = c(5, 4, 4, 4),  0.1)
   
   # Plot the first line (left axis)
   plot(as.numeric(names(ZA_occ_per_year)), ZA_occ_per_year,
        type = "l", col = "blue", lwd = 2,
        xlab = "Reference year", ylab = "Unique occurrences")
   
   # Overlay the second plot without redrawing axes/labels
   par(new = TRUE)
   plot(as.numeric(names(ZA_pub_per_year)), ZA_pub_per_year,
        type = "l", col = "orange", lwd = 2,
        axes = FALSE, xlab = "", ylab = "")
   
   # Add the right-side axis and its label
   axis(side = 4)
   mtext("Unique publications", side = 4, line = 3)
   
   # Legend
   legend("topleft", legend = c("Unique occurrences", "Unique publications"),
          col = c("blue", "orange"), lwd = 2)
   
   
   
US_pub_per_year <- tapply(only_US$primary_reference, only_US$ref_year, function(x) length(unique(x)))
  US_occ_per_year <- tapply(only_US$occurrence_no, only_US$ref_year, function(x) length(unique(x)))
  
   # Plot the first line
   plot(as.numeric(names(US_occ_per_year)), US_occ_per_year,
        type = "l", col = "blue", lwd = 2,
        xlab = "Reference year", ylab = "Count",
        ylim = range(c(US_occ_per_year, US_pub_per_year), na.rm = TRUE))
   
   # Add the second line to the same plot
   lines(as.numeric(names(US_pub_per_year)), ZA_pub_per_year,
         type = "l", col = "orange", lwd = 2)
   
   # Add a legend
   legend("topleft", legend = c("Unique occurrences", "Unique publications"),
          col = c("blue", "orange"), lwd = 2)
   # Set up margin space on the right for the second axis
   par(mar = c(5, 4, 4, 4),  0.1)
   
   # Plot the first line (left axis)
   plot(as.numeric(names(US_occ_per_year)), US_occ_per_year,
        type = "l", col = "blue", lwd = 2,
        xlab = "Reference year", ylab = "Unique occurrences")
   
   # Overlay the second plot without redrawing axes/labels
   par(new = TRUE)
   plot(as.numeric(names(US_pub_per_year)), US_pub_per_year,
        type = "l", col = "orange", lwd = 2,
        axes = FALSE, xlab = "", ylab = "")
   
   # Add the right-side axis and its label
   axis(side = 4)
   mtext("Unique publications", side = 4, line = 3)
   
   # Legend
   legend("topleft", legend = c("Unique occurrences", "Unique publications"),
          col = c("blue", "orange"), lwd = 2)
   

## ==================================================
##  Publication history information
## ==================================================
   species_data$ref_year <- as.integer(str_extract(species_data$primary_reference, "\\d{4}"))
print(sort(sakart[sakart$ref_author == "Romer", "primary_reference"]))

artkun<- species_result$only_US[species_result$only_US$stage == "Artinskian" | species_result$only_US$stage == "Kungurian",]
artkun_ref <- as.data.frame(table(artkun$primary_reference))

kas <- species_result$only_US[species_result$only_US$stage == "Kasimovian",]
kas_ref <- as.data.frame(table(kas$county))


sakart <- species_result$only_US[species_result$only_US$stage == "Sakmarian" | species_result$only_US$stage == "Artinskian" | species_result$only_US$stage == "Kungurian"
  ,]
sakart_ref <- as.data.frame(table(sakart$collection_name))
sakart_rest <- species_result$no_US[species_result$no_US$stage == "Sakmarian" | species_result$no_US$stage == "Artinskian" | species_result$only_US$stage == "Kungurian"
  ,]
sakart_rest_ref <- as.data.frame(table(sakart_rest$state))


tx_ok_by_year <- rowo %>%
  filter(!is.na(ref_year), state %in% c("Texas", "Oklahoma")) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "Texas + Oklahoma Red Beds")

# Group 2: Rest of USA (everything in sakart NOT Texas/Oklahoma)
rest_usa_by_year <- rowo %>%
  filter(!is.na(ref_year), !state %in% c("Texas", "Oklahoma")) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "Rest of USA")

# Group 3: Rest of world (sakart_rest = species_result$no_US)
rest_world_by_year <- rowo_rest %>%
  filter(!is.na(ref_year)) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "Rest of world")

three_groups <- bind_rows(tx_ok_by_year, rest_usa_by_year, rest_world_by_year)

three_groups_cumulative <- three_groups %>%
  group_by(group) %>%
  arrange(ref_year) %>%
  mutate(cum_n = cumsum(n)) %>%
  filter(ref_year >= min(ref_year[n > 0])) %>%   # drop any pre-first-paper padding
  ungroup()
three_points <- three_groups_cumulative %>%
  group_by(group) %>%
  filter(ref_year == min(ref_year)) %>%
  ungroup()

ggplot(three_groups_cumulative, aes(x = ref_year, y = cum_n, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(data = three_points, size = 2.5) +
  labs(title = "Cumulative Kungurian-Roadian occurrences",
       x = "Publication year", y = "Cumulative occurrences", color = NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_color_manual(values = c("Texas + Oklahoma Red Beds" = "darkorange",
                                 "Rest of USA" = "#1f77b4",
                                 "Rest of world" = "gray40"))


indole <- species_result$only_ZA[species_result$only_ZA$stage == "Induan" | species_result$only_ZA$stage == "Olenekian",]
indole_ref <- as.data.frame(table(indole$ref_author))
indole_rest <- species_result$no_ZA[species_result$no_ZA$stage == "Induan" | species_result$no_ZA$stage == "Olenekian",]
indole_rest_ref <- as.data.frame(table(indole_rest$cc))

print(sort(indole_rest[indole_rest$cc == "RU", "ref_pubyr"]))


# Count occurrences per year for each group
sa_by_year <- indole %>%
  filter(!is.na(ref_year)) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "South Africa\n(Beaufort Group)")

rus_by_year <- indole_rest %>%
  filter(!is.na(ref_year), cc %in% c("RU")) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "Russia")

# Group 3: Rest of world (sakart_rest = species_result$no_US)
rest_world_by_year <- indole_rest %>%
  filter(!is.na(ref_year), !cc %in% c("RU")) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "Rest of world")

three_groups <- bind_rows(sa_by_year, rus_by_year, rest_world_by_year)

three_groups_cumulative <- three_groups %>%
  group_by(group) %>%
  arrange(ref_year) %>%
  mutate(cum_n = cumsum(n)) %>%
  filter(ref_year >= min(ref_year[n > 0])) %>%   # drop any pre-first-paper padding
  ungroup()
three_points <- three_groups_cumulative %>%
  group_by(group) %>%
  filter(ref_year == min(ref_year)) %>%
  ungroup()

ggplot(three_groups_cumulative, aes(x = ref_year, y = cum_n, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(data = three_points, size = 2.5) +
  labs(title = "Cumulative Induan-Olenekian occurrences",
       x = "Publication year", y = "Cumulative occurrences", color = NULL) +
theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_color_manual(values = c("South Africa\n(Beaufort Group)" = "darkgreen",
                                "Russia" = "darkblue",
                                 "Rest of world" = "gray40"))

rowo <- species_result$only_US[species_result$only_US$stage == "Kungurian" | species_result$only_US$stage == "Roadian",]
rowo_ref <- as.data.frame(table(rowo$stage))
rowo_rest <- species_result$no_US[species_result$no_US$stage == "Kungurian" | species_result$no_US$stage == "Roadian",]
rowo_rest_ref <- as.data.frame(table(rowo_rest$stage))



sa_by_year <- rowo %>%
  filter(!is.na(ref_year)) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "South Africa\n(Beaufort Group)")

row_by_year <- rowo_rest %>%
  filter(!is.na(ref_year)) %>%
  count(ref_year, name = "n") %>%
  mutate(group = "Rest of world")

combined <- bind_rows(sa_by_year, row_by_year)
combined_cumulative <- combined %>%
  group_by(group) %>%
  arrange(ref_year) %>%
  mutate(cum_n = cumsum(n)) %>%
  filter(ref_year >= min(ref_year[n > 0])) %>%   # drop any pre-first-paper padding
  ungroup()
first_points <- combined_cumulative %>%
  group_by(group) %>%
  filter(ref_year == min(ref_year)) %>%
  ungroup()

ggplot(combined_cumulative, aes(x = ref_year, y = cum_n, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(data = first_points, size = 2.5) +
  labs(title = "Cumulative Kungurian-Roadian occurrences",
       x = "Publication year", y = "Cumulative occurrences", color = NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_color_manual(values = c("South Africa\n(Beaufort Group)" = "darkgreen",
                                 "Rest of world" = "gray40"))

print(sort(indole_rest[indole_rest$cc == "RU", "ref_pubyr"]))
roa <- species_result$only_ZA[species_result$only_ZA$stage == "Roadian",]
roa_ref <- as.data.frame(table(roa$state))
roa_rest <- species_result$no_ZA[species_result$no_ZA$stage == "Roadian",]
roa_rest_ref <- as.data.frame(table(roa_rest$state))