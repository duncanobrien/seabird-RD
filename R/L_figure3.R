##########################################################################################
#### Preamble ####
##########################################################################################

library(ggplot2)
library(ggdist)
library(patchwork)
library(dplyr)

##########################################################################################
#### Create figure ####
##########################################################################################

behav_imbalance_a <- readRDS("outputs/stability/behav_imbalance_atpu.RDS")
behav_imbalance_b <- readRDS("outputs/stability/behav_imbalance_blki.RDS")
behav_imbalance_cor <- lapply(behav_imbalance_a, `[[`, "rho") |>
  do.call("rbind", arg = _) |>
  dplyr::mutate(species = "Puffins") |>
  rbind(
    lapply(behav_imbalance_b, `[[`, "rho") |>
      do.call("rbind", arg = _) |>
      dplyr::mutate(species = "Kittiwakes")
  ) |>
  dplyr::select(rho_dev_detrend,rho_ts_detrend,species,mc_sim) |>
  dplyr::slice_sample(n = 1000, by = c(species)) |>
  tidyr::pivot_longer(
    dplyr::starts_with("rho"),
    names_to = "response",
    values_to = ".rho"
  ) |>
  dplyr::mutate(response = gsub("rho_|_detrend", "", response)) |>
  dplyr::mutate(
    response = dplyr::if_else(response == "dev", "Dominant eigenvalue", "Temporal stability")
  )

behav_dissim_a <- readRDS("outputs/stability/behav_dissim_atpu.RDS")
behav_dissim_b <- readRDS("outputs/stability/behav_dissim_blki.RDS")
behav_dissim_cor <- lapply(behav_dissim_a, `[[`, "rho") |>
  do.call("rbind", arg = _) |>
  dplyr::mutate(species = "Puffins") |>
  rbind(
    lapply(behav_dissim_b, `[[`, "rho") |>
      do.call("rbind", arg = _) |>
      dplyr::mutate(species = "Kittiwakes")
  ) |>
  dplyr::select(rho_dev_detrend,rho_ts_detrend,species,mc_sim) |>
  dplyr::slice_sample(n = 1000, by = c(species)) |>
  tidyr::pivot_longer(
    dplyr::starts_with("rho"),
    names_to = "response",
    values_to = ".rho"
  ) |>
  dplyr::mutate(response = gsub("rho_|_detrend", "", response)) |>
  dplyr::mutate(
    response = dplyr::if_else(response == "dev", "Dominant eigenvalue", "Temporal stability")
  )

dee_imbalance_a <- readRDS("outputs/stability/dee_imbalance_atpu.RDS")
dee_imbalance_b <- readRDS("outputs/stability/dee_imbalance_blki.RDS")
dee_imbalance_cor <- lapply(dee_imbalance_a, `[[`, "rho") |>
  do.call("rbind", arg = _) |>
  dplyr::mutate(species = "Puffins") |>
  rbind(
    lapply(dee_imbalance_b, `[[`, "rho") |>
      do.call("rbind", arg = _) |>
      dplyr::mutate(species = "Kittiwakes")
  ) |>
  dplyr::select(rho_dev_detrend,rho_ts_detrend,species,mc_sim) |>
  dplyr::slice_sample(n = 1000, by = c(species)) |>
  tidyr::pivot_longer(
    dplyr::starts_with("rho"),
    names_to = "response",
    values_to = ".rho"
  ) |>
  dplyr::mutate(response = gsub("rho_|_detrend", "", response)) |>
  dplyr::mutate(
    response = dplyr::if_else(response == "dev", "Dominant eigenvalue", "Temporal stability")
  )

dee_dissim_a <- readRDS("outputs/stability/dee_dissim_atpu.RDS")
dee_dissim_b <- readRDS("outputs/stability/dee_dissim_blki.RDS")
dee_dissim_cor <- lapply(dee_dissim_a, `[[`, "rho") |>
  do.call("rbind", arg = _) |>
  dplyr::mutate(species = "Puffins") |>
  rbind(
    lapply(dee_dissim_b, `[[`, "rho") |>
      do.call("rbind", arg = _) |>
      dplyr::mutate(species = "Kittiwakes")
  ) |>
  dplyr::select(rho_dev_detrend,rho_ts_detrend,species,mc_sim) |>
  dplyr::slice_sample(n = 1000, by = c(species)) |>
  tidyr::pivot_longer(
    dplyr::starts_with("rho"),
    names_to = "response",
    values_to = ".rho"
  ) |>
  dplyr::mutate(response = gsub("rho_|_detrend", "", response)) |>
  dplyr::mutate(
    response = dplyr::if_else(response == "dev", "Dominant eigenvalue", "Temporal stability")
  )

p1 <- ggplot(
  behav_imbalance_cor |>
    dplyr::mutate(term = "Imbalance") |>
    rbind(behav_dissim_cor  |>
            dplyr::mutate(term = "Dissimilarity")),
  aes(x = term, y = .rho)
) +
  stat_halfeye(aes(
    term = term,
    response = response,
    fill = after_stat(ifelse((term == "Dissimilarity" &
                                response == "Temporal stability") |
                               (term == "Imbalance" &
                                  response == "Dominant eigenvalue"),
                             y > 0,
                             y < 0
    ))
  ),
  .width = 0.95,
  scale = 0.5,
  fatten_point = 1,
  linewidth = 1,
  normalize = "groups",
  check.aes = FALSE) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0)) +
  coord_cartesian(ylim = c(-1.0, 1.0)) +
  ylab("Pearson correlation") +
  xlab("Foraging effort response diversity") +
  scale_fill_manual(values = c("gray85", "#CCCCFF"), name = "Aligned with\nhypothesised\ndirection") +
  facet_grid(species ~ response, scales = "free_y") +
  ggplot(
    dee_imbalance_cor |>
      dplyr::mutate(term = "Imbalance") |>
      rbind(dee_dissim_cor  |>
              dplyr::mutate(term = "Dissimilarity")) ,
    aes(x = term, y = .rho)
  ) +
  stat_halfeye(aes(
    term = term,
    response = response,
    fill = after_stat(ifelse((term == "Dissimilarity" &
                                response == "Temporal stability") |
                               (term == "Imbalance" &
                                  response == "Dominant eigenvalue"),
                             y > 0,
                             y < 0
    ))
  ),
  .width = 0.95,
  scale = 0.5,
  fatten_point = 1,
  linewidth = 1,
  normalize = "groups",
  check.aes = FALSE) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0)) +
  coord_cartesian(ylim = c(-1.0, 1.0)) +
  ylab("Pearson correlation") +
  xlab("Energy expenditure response diversity") +
  scale_fill_manual(values = c("gray85", "#CCCCFF"), name = "Aligned with\nhypothesised\ndirection") +
  facet_grid(species ~ response, scales = "free_y") +
  plot_layout(nrow = 2, guides = "collect") &
  guides(fill = guide_legend(override.aes = list(linetype = 0))) &
  theme_classic() &
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0,
                                hjust = 0),
    strip.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    legend.margin=margin(l = 0, unit='cm'),
    legend.box.margin=margin(l = -1, unit='cm'),
    legend.title = element_text(face = "bold.italic"),
  )  &
  plot_annotation(tag_levels = "A")

ggsave("outputs/figure3_spp_cor.png", p1, dpi = 300, width = 6, height = 6)
