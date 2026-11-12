
#jogsmith@ucsc.edu

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, readxl, here,
  ggridges, patchwork, sf, lubridate, emmeans, sandwich
)

################################################################################
# 1) Load all data --------------------------------------------------------------
################################################################################

# SOFA results
sofa_path <- here::here("output", "sofa", "SOFA_Results_2025-10-31.xlsx")
eta_df    <- read_xlsx(sofa_path, sheet = "propforagingeffort")

# Zone-level data
load(here::here("output", "survey_data", "processed", "zone_level_data.rda"))


# Size-frequency data
sizefq_with_pred <- readRDS(
  here::here("output", "survey_data", "processed", "sizefq_with_pred.rds")
)


################################################################################
#Prep SOFA data 

patch_levels <- c("FOR", "INCIP", "BAR")

recode_patches <- function(df) {
  df %>%
    filter(prey_name %in% c("urchin_b", "urchin_f", "urchin_i")) %>%
    mutate(
      patch_type = case_when(
        prey_name == "urchin_b" ~ "BAR",
        prey_name == "urchin_f" ~ "FOR",
        prey_name == "urchin_i" ~ "INCIP"
      ),
      patch_type = factor(patch_type, levels = patch_levels)
    )
}

eta_focus <- recode_patches(eta_df)

expand_posterior <- function(tbl, n = 5000) {
  tbl %>%
    rowwise() %>%
    mutate(draws = list(rnorm(n, mean, sd))) %>%
    unnest(draws) %>%
    ungroup()
}

eta_draws <- expand_posterior(eta_focus)

posterior_summary <- function(df) {
  df %>%
    group_by(patch_type) %>%
    summarise(
      med = median(draws),
      l90 = quantile(draws, 0.05),
      u90 = quantile(draws, 0.95),
      .groups = "drop"
    )
}

eta_sum <- posterior_summary(eta_draws)


################################################################################
#Prep quad data 

quad_plot_dat <- quad_build3 %>%
  st_drop_geometry() %>%
  mutate(
    patch_type = factor(pred_patch, levels = patch_levels),
    biomass_density_g_m2 = total_biomass_g / 80,
    exposed_prop = if_else(
      purple_urchin_densitym2 > 0,
      1 - purple_urchin_conceiledm2 / purple_urchin_densitym2,
      NA_real_
    )
  ) %>%
  filter(!is.na(patch_type))

################################################################################
#Prep size-frequency data 

sizefq_with_pred <- sizefq_with_pred %>%
  mutate(
    site_type_predicted = factor(site_type_predicted, levels = patch_levels)
  )

size_curve <- sizefq_with_pred %>%
  filter(
    species == "Purple",
    !is.na(site_type_predicted),
    !is.na(size_cm),
    !is.na(count)
  ) %>%
  group_by(site_type_predicted, size_cm) %>%
  summarise(n = sum(count), .groups = "drop")

################################################################################
#Plot functions 

patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

theme_panel <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      axis.text        = element_text(size = 12, color = "black"),
      axis.title       = element_text(size = 13, color = "black"),
      legend.title     = element_text(size = 13),
      legend.text      = element_text(size = 12),
      plot.tag         = element_text(size = 16, color = "black", face = "plain"),
      strip.background = element_blank(),
      strip.text       = element_text(size = 12, color = "black")
    )
}

plot_posterior <- function(draw_df, sum_df, xlab, tag_title, show_legend = FALSE) {
  
  ggplot() +
    geom_density_ridges(
      data = draw_df,
      aes(x = draws, y = patch_type, fill = patch_type),
      alpha = 0.85,
      scale = 0.9,
      color = "black",
      linewidth = 0.3
    ) +
    geom_errorbar(
      data = sum_df,
      aes(y = patch_type, xmin = l90, xmax = u90),
      orientation = "y",
      height = 0.25,
      linewidth = 1.0,
      color = "black"
    ) +
    geom_point(
      data = sum_df,
      aes(y = patch_type, x = med),
      color = "white",
      fill = "black",
      shape = 21,
      size = 2.8,
      stroke = 0.9
    ) +
    scale_fill_manual(values = patch_colors, guide = "none") +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_panel() +
    theme(
      axis.text.y     = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )
}

plot_box <- function(df, xvar, xlab, tag_title, xlimits = NULL) {
  
  p <- ggplot(
    df,
    aes(x = .data[[xvar]], y = patch_type, fill = patch_type)
  ) +
    geom_boxplot(
      width = 0.62,
      outlier.shape = NA,
      alpha = 0.82,
      linewidth = 0.35,
      color = "black"
    ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.25,
      size = 1.8,
      alpha = 0.52,
      position = position_jitter(height = 0.12, width = 0, seed = 6831),
      show.legend = FALSE
    ) +
    scale_fill_manual(values = patch_colors, guide = "none") +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_panel() +
    theme(
      axis.text.y = element_blank(),
      legend.position = "none"
    )
  
  if (!is.null(xlimits)) {
    p <- p + coord_cartesian(xlim = xlimits)
  }
  
  p
}

fit_pairwise <- function(data, response, log_response = FALSE) {
  model_data <- data %>%
    filter(
      is.finite(.data[[response]]),
      if (log_response) .data[[response]] > 0 else TRUE
    ) %>%
    mutate(
      year = factor(year),
      test_response = if (log_response) log(.data[[response]]) else .data[[response]]
    )

  model <- lm(test_response ~ patch_type + year, data = model_data)
  emm <- emmeans::emmeans(
    model,
    ~ patch_type,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )

  list(
    model = model,
    contrasts = as.data.frame(
      summary(pairs(emm, adjust = "holm"))
    )
  )
}

add_pairwise_brackets <- function(plot, contrasts, x_max, x_limits,
                                  multipliers, label_add = 0.045) {
  pair_groups <- strsplit(as.character(contrasts$contrast), " - ", fixed = TRUE)
  bracket_dat <- tibble(
    group1 = vapply(pair_groups, `[`, character(1), 1),
    group2 = vapply(pair_groups, `[`, character(1), 2),
    p_value = contrasts$p.value,
    x = x_max * multipliers
  ) %>%
    mutate(
      y1 = match(group1, levels(quad_plot_dat$patch_type)),
      y2 = match(group2, levels(quad_plot_dat$patch_type)),
      cap_x = x - x_max * 0.025,
      label_x = x + x_max * label_add,
      label_y = (y1 + y2) / 2 + c(-0.14, 0, 0.14),
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )

  plot +
    geom_segment(
      data = bracket_dat,
      aes(x = x, y = y1, xend = x, yend = y2),
      inherit.aes = FALSE,
      linewidth = 0.45
    ) +
    geom_segment(
      data = bracket_dat,
      aes(x = cap_x, y = y1, xend = x, yend = y1),
      inherit.aes = FALSE,
      linewidth = 0.45
    ) +
    geom_segment(
      data = bracket_dat,
      aes(x = cap_x, y = y2, xend = x, yend = y2),
      inherit.aes = FALSE,
      linewidth = 0.45
    ) +
    geom_text(
      data = bracket_dat,
      aes(x = label_x, y = label_y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4.4,
      fontface = "bold"
    ) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    theme(plot.margin = margin(5.5, 24, 5.5, 5.5))
}

################################################################################
#build figure


# Panel A
pA <- plot_posterior(
  draw_df     = eta_draws,
  sum_df      = eta_sum,
  xlab        = "Sea otter foraging effort (proportion of \ntime on sea urchins in patch vs. all other prey)",
  tag_title   = "A",
  show_legend = FALSE
) +
  scale_y_discrete(
    limits = c("FOR", "INCIP", "BAR"),
    labels = c(FOR = "Forest", INCIP = "Incipient", BAR = "Barren"),
    expand = expansion(add = c(0.1, 0.6))
  )

gi_pairwise <- fit_pairwise(quad_plot_dat, "mean_gi", log_response = TRUE)
biomass_pairwise <- fit_pairwise(
  quad_plot_dat, "biomass_density_g_m2", log_response = TRUE
)
exposure_pairwise <- fit_pairwise(quad_plot_dat, "exposed_prop")

# Panel B
pB <- quad_plot_dat %>%
  filter(is.finite(mean_gi), mean_gi > 0) %>%
  plot_box(
    xvar      = "mean_gi",
    xlab      = "Mean purple urchin gonad index \n (ratio of gonad mass to animal mass)",
    tag_title = "B"
  ) %>%
  add_pairwise_brackets(
    gi_pairwise$contrasts,
    x_max = max(quad_plot_dat$mean_gi, na.rm = TRUE),
    x_limits = c(0, 28),
    multipliers = c(1.05, 1.30, 1.55)
  )

# Panel C
pC <- quad_plot_dat %>%
  filter(is.finite(biomass_density_g_m2), biomass_density_g_m2 > 0) %>%
  plot_box(
    xvar      = "biomass_density_g_m2",
    xlab      = expression("Purple urchin biomass density (g " * m^{-2} * ")"),
    tag_title = "C"
  ) +
  scale_x_log10(breaks = c(10, 100, 1000))

pC <- add_pairwise_brackets(
  pC,
  biomass_pairwise$contrasts,
  x_max = max(quad_plot_dat$biomass_density_g_m2, na.rm = TRUE),
  x_limits = c(3, 12000),
  multipliers = c(1.10, 1.75, 2.80),
  label_add = 0.06
)

# Panel D
pD <- quad_plot_dat %>%
  filter(is.finite(exposed_prop)) %>%
  plot_box(
    xvar = "exposed_prop",
    xlab = "Proportion of purple urchins exposed",
    tag_title = "D"
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = scales::label_percent()
  )

pD <- add_pairwise_brackets(
  pD,
  exposure_pairwise$contrasts,
  x_max = 1,
  x_limits = c(0, 1.7),
  multipliers = c(1.05, 1.25, 1.45),
  label_add = 0.06
)

pE <- ggplot(
  size_curve,
  aes(x = size_cm, y = n, color = site_type_predicted)
) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(
    values = patch_colors,
    breaks = c("BAR", "INCIP", "FOR"),
    labels = c("Barren", "Incipient", "Forest"),
    name = "Patch state",
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  scale_x_continuous(limits = c(1, 8), breaks = 1:8) +
  labs(x = "Purple urchin size frequency (cm)", y = NULL, tag = "E") +
  theme_panel()

################################################################################
# merge

energetics_combined <- wrap_plots(
  pA, pB, pC, pD, pE,
  design = "ABC\nADE",
  widths = c(1.2, 0.9, 1.2),
  guides = "collect"
) &
  theme(
    plot.margin = unit(c(2, 2, 2, 2), "pt"),
    legend.position = "bottom"
  )

energetics_combined

################################################################################
#export

ggsave(
  here::here("figures", "Fig4_effort_allocationv2.png"),
  energetics_combined,
  width  = 16,
  height = 9.5,
  dpi    = 600,
  bg     = "white"
)


################################################################################
# Summary stats

eta_sum %>%
  mutate(
    med_pct = med * 100,
    l90_pct = l90 * 100,
    u90_pct = u90 * 100
  ) %>%
  select(patch_type, med_pct, l90_pct, u90_pct) %>%
  arrange(patch_type) %>%
  print()

gi_summary <- quad_plot_dat %>%
  dplyr::filter(!is.na(mean_gi)) %>%
  dplyr::group_by(patch_type) %>%
  dplyr::summarise(
    n      = sum(is.finite(mean_gi)),
    mean   = mean(mean_gi, na.rm = TRUE),
    sd     = sd(mean_gi, na.rm = TRUE),
    se     = sd / sqrt(n),
    median = median(mean_gi, na.rm = TRUE),
    q25    = quantile(mean_gi, 0.25, na.rm = TRUE),
    q75    = quantile(mean_gi, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(patch_type, levels = c("BAR", "INCIP", "FOR")))

gi_summary

biomass_summary <- quad_plot_dat %>%
  dplyr::filter(!is.na(biomass_density_g_m2)) %>%
  dplyr::group_by(patch_type) %>%
  dplyr::summarise(
    n      = sum(is.finite(biomass_density_g_m2)),
    mean   = mean(biomass_density_g_m2, na.rm = TRUE),
    sd     = sd(biomass_density_g_m2, na.rm = TRUE),
    se     = sd / sqrt(n),
    median = median(biomass_density_g_m2, na.rm = TRUE),
    q25    = quantile(biomass_density_g_m2, 0.25, na.rm = TRUE),
    q75    = quantile(biomass_density_g_m2, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(patch_type, levels = c("BAR", "INCIP", "FOR")))

biomass_summary
