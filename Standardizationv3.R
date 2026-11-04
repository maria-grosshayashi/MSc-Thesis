# ============================================================
# Tetrapodomorpha Standardized Diversity Analysis (Species & Genus)
# ============================================================

setwd("C:/Users/mlove/Desktop/MSc-Paleobio/semester-3/Research/Data sets")

##  Packages 
library(divDyn)
library(tidyverse)
  library(ggtext)
  library(scales)
  library(scatterpie)
  library(deeptime)
library(iNEXT)
library(devtools)


##  Load data 
stages <- read.csv("divDyn_2026.csv")[4:27, ]  # redated to analyze 2026 stages vs 2026 data

species_data <- read.csv("pbdb_data_species_overlap_v2.csv", header = TRUE, sep = ",")
genus_data <- read.csv("pbdb_data_genus_overlap_v2.csv", header = TRUE, sep = ",")


## ============================================================
## METRICS FUNCTIONS (binning, standardization, permutation)
## ============================================================

## bin each occurrence to the stage with maximum temporal overlap
binning <- function(data) {
  data <- data[!is.na(data$max_ma) & !is.na(data$min_ma), ]   # Remove rows that lack dating
  
  data[data$max_ma == data$min_ma, "min_ma"] <- data[data$max_ma == data$min_ma, "min_ma"] - 1e-6  # avoid zero-width ranges
    min_end     <- outer(data$max_ma, stages$bottom, pmin)  # youngest max age per occurrence data interval
    max_start   <- outer(data$min_ma, stages$top, pmax)     # oldest min age per occurrence data interval
    overlap_mat <- min_end - max_start                    # overlap = the difference
    best_bin    <- max.col(overlap_mat, ties.method = "first")
    data$stage <- stages$stage[best_bin]
    data$bin   <- stages$stg[best_bin]
  
    data
}

## Create standardization function
standardize <- function(data, quorum_level){
  
  ## Transform into incidence data
  incidence_data <- lapply(1:nrow(stages), function(i) { #intervals is a dataframe containing interval names, midpoints, and max and min Ma
    tmp_data <- data %>% filter(bin == stages$stg[i])
    incidence_raw <- tmp_data %>% .[,c("collection_no","accepted_name")] %>% table %>% t
    incidence_raw <- as.data.frame.matrix(incidence_raw)
    incidence_raw[incidence_raw > 1] <- 1
    if(nrow(incidence_raw) > 1) {incidence_raw}
    })

  names(incidence_data) <- stages$stg

  incidence_data <- incidence_data[!sapply(incidence_data,is.null)] #remove stages without incidence data

  ## Standardize
  estD <- estimateD(incidence_data, q=0, datatype="incidence_raw", base="coverage", level=quorum_level)
  estD$reference_t <- sapply(incidence_data, sum) 
  estD[which(estD$t >= 2 * estD$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) #no more than twice reference sample size

  estD$Assemblage <- as.integer(estD$Assemblage)
  estD <- estD %>% rename(stg = Assemblage) %>% full_join(., stages, by = "stg") %>% rename(bin = stg) #intervals is a dataframe containing interval names, midpoints, and max and min Ma

  estD  <- estD[!is.na(estD$Method),] #remove NAs

  ## Returns
  estD[estD$bin %in% stages$stg, ]
}


## Country-exclusion permutation test: repeatedly removes a random sample of occurrences (matched per-bin to the target country's actual contribution) and recomputes log_relative_change, building a null distribution per bin.
permutation <- function(data, country_code, ref_bins, n_perm = 1000, quorum_level) {
  set.seed(2)
  counts_per_bin <- data %>% filter(cc == country_code) %>% count(bin, name = "n_country")
  counts_left    <- data %>% filter(cc != country_code) %>% count(bin, name = "n_bin")

  null_matrix <- matrix(nrow = length(ref_bins), ncol = n_perm)

  for (i in seq_len(n_perm)) {
   perm_slice <- data %>%
      group_by(bin) %>%
      group_modify(~ {
      n_remove <- counts_per_bin$n_country[counts_per_bin$bin == .y$bin]
      if (length(n_remove) == 0) n_remove <- 0
      n_remove <- min(n_remove, nrow(.x))
      slice_sample(.x, n = nrow(.x) - n_remove, replace = FALSE)
      }) %>%
    ungroup()

    perm_div <- perm_slice %>%
      binning() %>%
      standardize(quorum_level = quorum_level) %>%
      arrange(bin) %>%
      mutate(log_relative_change = log10(qD / lag(qD)))

    null_matrix[, i] <- perm_div$log_relative_change[match(ref_bins, perm_div$bin)]
  }

  null_matrix[!is.finite(null_matrix)] <- NA
  list(null_matrix = null_matrix, counts_per_bin = counts_per_bin, counts_left = counts_left)
}

## ============================================================
## ANALYSIS FUNCTION
## ============================================================

run_standardized_analysis <- function(data, resolution_label, n_perm = 1000) {

  ##  bin 
  data <- binning(data)
  
  ##  clean, subset
  data <- data[(data$max_ma - data$min_ma) < 25 &
                           !is.na(data$max_ma) & !is.na(data$min_ma), ]
 
  no_US <- data[data$cc != "US", ]
  no_ZA <- data[data$cc != "ZA", ]

  ##  Standardize 
  estD_full <- standardize(data, quorum_level = 0.4)
  estD_noUS <- standardize(no_US, quorum_level = 0.4)
  estD_noZA <- standardize(no_ZA, quorum_level = 0.4)

  
  ## Stage-to-stage log10 relative change 
  estD_full <- estD_full %>% arrange(bin) %>%
    mutate(log_relative_change = log10(qD / lag(qD)), n = qD + lag(qD))
  estD_noUS <- estD_noUS %>% arrange(bin) %>%
    mutate(log_relative_change = log10(qD / lag(qD)), n = qD + lag(qD))
  estD_noZA <- estD_noZA %>% arrange(bin) %>%
    mutate(log_relative_change = log10(qD / lag(qD)), n = qD + lag(qD))

  first_valid_index <- min(which(is.finite(estD_full$log_relative_change) &
                                    is.finite(estD_noUS$log_relative_change) &
                                    is.finite(estD_noZA$log_relative_change)))
  valid_bins <- stages$stg[first_valid_index:length(stages$stg)]
  
  ## Permutation tests (country-exclusion null distributions) 
  us_result <- permutation(data, "US", ref_bins = valid_bins, n_perm = n_perm, quorum_level = 0.4)
  za_result <- permutation(data, "ZA", ref_bins = valid_bins, n_perm = n_perm, quorum_level = 0.4)

  ##  Stage lookup: pairs each bin with the proper comparison 
  stage_lookup <- data.frame(bin = valid_bins) %>%
    rowwise() %>%
    mutate(stage_comp = {
      idx <- which(stages$stg == bin)
      if (length(idx) == 0 || idx == 1) NA_character_
      else paste(stages$short[idx - 1], "-", stages$short[idx])
    }) %>%
    ungroup()

  ## Null distributions for both countries 
  us_null_df <- as.data.frame(us_result$null_matrix)
  colnames(us_null_df) <- paste0("perm_", seq_len(ncol(us_null_df)))
  us_null_df$bin <- valid_bins
  us_null_long <- us_null_df %>%
    filter(bin %in% valid_bins) %>%
    pivot_longer(cols = starts_with("perm_"), names_to = "perm_id", values_to = "null_value") %>%
    mutate(country = "USA") %>%
    left_join(stage_lookup, by = "bin")

  za_null_df <- as.data.frame(za_result$null_matrix)
  colnames(za_null_df) <- paste0("perm_", seq_len(ncol(za_null_df)))
  za_null_df$bin <- valid_bins
  za_null_long <- za_null_df %>%
    filter(bin %in% valid_bins) %>%
    pivot_longer(cols = starts_with("perm_"), names_to = "perm_id", values_to = "null_value") %>%
    mutate(country = "South Africa") %>%
    left_join(stage_lookup, by = "bin")

  null_long <- bind_rows(us_null_long, za_null_long)

  ## Observed relative-change values for both countries 
  observed_df <- bind_rows(
    estD_noUS %>% filter(bin %in% valid_bins) %>% select(bin, log_relative_change) %>%
      rename(observed = log_relative_change) %>% mutate(country = "USA") %>%
      left_join(stage_lookup, by = "bin"),
    estD_noZA %>% filter(bin %in% valid_bins) %>% select(bin, log_relative_change) %>%
      rename(observed = log_relative_change) %>% mutate(country = "South Africa") %>%
      left_join(stage_lookup, by = "bin")
  )
  observed_df$observed[!is.finite(observed_df$observed)] <- NA

  ## P-values: compare observed change to the null distribution, per bin 
  us_p <- sapply(valid_bins, function(b) {
    obs <- estD_noUS$log_relative_change[estD_noUS$bin == b]
    if (length(obs) == 0) return(NA)
    null_vals <- us_result$null_matrix[which(valid_bins == b), ]
    med <- median(null_vals, na.rm = TRUE)
    mean(c(abs(null_vals - med) >= abs(obs - med), TRUE), na.rm = TRUE)
  })

  za_p <- sapply(valid_bins, function(b) {
    obs <- estD_noZA$log_relative_change[estD_noZA$bin == b]
    if (length(obs) == 0) return(NA)
    null_vals <- za_result$null_matrix[which(valid_bins == b), ]
    med <- median(null_vals, na.rm = TRUE)
    mean(c(abs(null_vals - med) >= abs(obs - med), TRUE), na.rm = TRUE)
  })

  p_df <- bind_rows(
    data.frame(bin = valid_bins, p_value = us_p,
               n_country = us_result$counts_per_bin$n_country[match(valid_bins, us_result$counts_per_bin$bin)],
               n_bin = us_result$counts_left$n_bin[match(valid_bins, us_result$counts_left$bin)],
               country = "USA"),
    data.frame(bin = valid_bins, p_value = za_p,
               n_country = za_result$counts_per_bin$n_country[match(valid_bins, za_result$counts_per_bin$bin)],
               n_bin = za_result$counts_left$n_bin[match(valid_bins, za_result$counts_left$bin)],
               country = "South Africa")
  ) %>%
    filter(n_country > 0) %>%                          # bins with 0 occurrences of that country: no test possible
    mutate(p_adj = p.adjust(p_value, method = "BH"))    # BH FDR correction across all bins x both countries

  observed_df <- observed_df %>%
    left_join(p_df %>% select(bin, country, p_adj), by = c("bin", "country")) %>%
    mutate(significance_cat = case_when(
      is.na(p_adj) ~ "No occurrences in bin",
      p_adj < 0.05 ~ "p < 0.05",
      TRUE ~ "not significant"
    ))

  p_df$n_bin[is.na(p_df$n_bin)] <- 0
  p_df$n <- p_df$n_country + p_df$n_bin

  ## Return all 
  list(
    resolution_label = resolution_label,
    valid_bins = valid_bins,
    data = data, no_US = no_US, no_ZA = no_ZA,
    estD_full = estD_full, estD_noUS = estD_noUS, estD_noZA = estD_noZA,
    us_result = us_result, za_result = za_result,
    stage_lookup = stage_lookup,
    p_df = p_df, observed_df = observed_df, null_long = null_long
  )
}


##  Run Analyses  
species_result  <- run_standardized_analysis(species_data, "species")
genus_result  <- run_standardized_analysis(genus_data, "genus")



## ============================================================
## VISUALIZATION
## ============================================================

visual_stages <- read.csv("Visual_stages.csv",)
visual_periods <- read.csv("Visual_periods.csv",)

make_standardized_plots <- function(results, top_country = "South Africa") {

  resolution_label <- results$resolution_label
  stage_lookup <- results$stage_lookup

  ##  Plot 1: diversity-through-time line plot 
  line_data <- bind_rows(
    results$estD_full %>% select(bin, qD) %>% mutate(method = "Full Global"),
    results$estD_noUS %>% select(bin, qD) %>% mutate(method = "USA Excluded"),
    results$estD_noZA %>% select(bin, qD) %>% mutate(method = "ZA Excluded")
  ) %>%
    left_join(stages %>% select(stg, mid), by = c("bin" = "stg"))
  line_data$method <- factor(line_data$method, levels = c("USA Excluded", "ZA Excluded", "Full Global"))

  line_plot <- ggplot(line_data, aes(x = mid, y = qD, color = method)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(
      values = c("Full Global" = "black", "ZA Excluded" = "#D55E00", "USA Excluded" = "#0072B2"),
      name = NULL
    ) +
    scale_x_reverse(limits = c(375, 240)) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
    labs(y = paste0("Standardized diversity (", resolution_label, ")"), x = "Time (Ma)",
         title = paste0("Standardized *Tetrapodomorpha* ",
                         str_to_title(resolution_label), " Gamma Diversity")) +
    coord_geo(
      dat = list(visual_stages, visual_periods), xlim = c(375, 240), pos = list("b", "b"),
      height = unit(1.2, "line"), size = list(2.5, 3), abbrv = list(TRUE, FALSE)
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = c(0.2, 0.85), legend.background = element_blank(),
          plot.title = element_markdown(face = "bold", size = 11, hjust = 0.5))

  ##  Plot 2: violin plot of null distributions vs. observed change 
  violin_plot <- ggplot(results$null_long, aes(x = factor(bin), y = null_value)) +
    geom_violin(fill = "gray60", colour = NA, width = 1, scale = "width") +
    geom_point(data = results$observed_df, aes(x = factor(bin), y = observed, colour = significance_cat),
               size = 2, stroke = 1.2) +
    scale_x_discrete(labels = setNames(stage_lookup$stage_comp, stage_lookup$bin)) +
    scale_colour_manual(values = c("p < 0.05" = "red", "not significant" = "gray20",
                                   "No occurrences in bin" = "gray80")) +
    facet_wrap(~ country, ncol = 1, scales = "free") +
    labs(x = "Stages Compared", y = "Log10 Relative Change", colour = "Significance",
         title = paste0("Relative change of ", resolution_label, " richness through time")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 10, face = "bold"),
          plot.caption = element_text(hjust = 0, size = 10))

  ##  Plot 3: scatterpie summary of adjusted p-values by bin 
  p_df <- results$p_df
  spacing_factor <- 2
  scale_factor   <- 8

  bin_levels <- levels(factor(p_df$bin))
  p_df$bin <- as.numeric(factor(p_df$bin))
  stage_lookup$bin_numeric <- match(as.character(stage_lookup$bin), bin_levels)

  p_df$bin_spaced <- p_df$bin * spacing_factor
  stage_lookup$bin_numeric_spaced <- stage_lookup$bin_numeric * spacing_factor
  p_df$p_adj_scaled <- sqrt(p_df$p_adj) * scale_factor

  ## drop bins with no computable p-value (NA radius/y crashes geom_scatterpie's stat)
  p_df <- p_df %>% filter(is.finite(p_adj_scaled))

  y_range <- range(p_df$p_adj_scaled, na.rm = TRUE)
  y_pad   <- diff(y_range) * 0.05

  legend_radii  <- log(c(50, 900)) / 8
  legend_layers <- geom_scatterpie_legend(
    r = legend_radii, breaks = legend_radii, n = length(legend_radii),
    x = max(p_df$bin_spaced), y = max(p_df$p_adj_scaled) - 0.5,
    labeller = function(x) round(exp(x * 8), 0), label_position = "left"
  )
  legend_layers <- lapply(legend_layers, function(l) { l$data$country <- top_country; l })

  summary_plot <- ggplot(p_df, aes(x = bin_spaced, y = p_adj_scaled)) +
    geom_scatterpie(aes(x = bin_spaced, y = p_adj_scaled, r = log(n) / 8, colour = country),
                     cols = c("n_bin", "n_country")) +
    legend_layers +
    coord_fixed(ratio = 1, xlim = c(min(p_df$bin_spaced), max(p_df$bin_spaced)),
                ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)) +
    geom_hline(yintercept = sqrt(0.05) * scale_factor, linetype = "dashed", colour = "red") +
    facet_wrap(~ country, ncol = 1) +
    scale_y_continuous(breaks = sqrt(c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00)) * scale_factor,
                       labels = c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00)) +
    scale_x_continuous(breaks = stage_lookup$bin_numeric_spaced, labels = stage_lookup$stage_comp,
                        guide = guide_axis(angle = 45)) +
    scale_fill_manual(name = "Occurrence type",
                       values = c("n_country" = "#66CCEE", "n_bin" = "orange"),
                       labels = c("n_country" = "Occurrences in target country", "n_bin" = "Other global occurrences")) +
    scale_colour_manual(values = setNames(rep("black", length(unique(p_df$country))), unique(p_df$country)),
                         guide = "none") +
    labs(x = "Bin", y = "Adjusted p-value (BH)",
         title = paste0("Significance of country-specific ", resolution_label,
                         " removal effect by bin (FDR-corrected)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.caption = element_text(hjust = 0, size = 10),
          legend.position = "right", legend.box = "vertical",
          panel.spacing.x = unit(1, "lines"))

  list(line_plot = line_plot, violin_plot = violin_plot, summary_plot = summary_plot)
}

species_plots <- make_standardized_plots(species_result)
genus_plots <- make_standardized_plots(genus_result)

##  View plots 
species_plots$line_plot
species_plots$violin_plot
species_plots$summary_plot

genus_plots$line_plot
genus_plots$violin_plot
genus_plots$summary_plot