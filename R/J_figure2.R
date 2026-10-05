##########################################################################################
#### Preamble ####
##########################################################################################

library(ggplot2)
library(patchwork)
library(dplyr)

model_files <- fs::dir_ls(path = "outputs/models", glob = "*.RDS")
behav_models <- lapply(grep("behav",model_files, value = TRUE),readRDS)
behav_info <- lapply(behav_models,`[[`,"sst_info")
behav_derivs <- lapply(behav_models,`[[`,"derivs") |>
  lapply(function(X) dplyr::mutate(X, median = median(estimate),.by = c(species,colony,bird1,sst)) |>
           subset(drawid %in% sample(unique(drawid), 50, replace = FALSE))) |>
  do.call("rbind", args = _) |>
  tibble::remove_rownames() |>
  dplyr::left_join(do.call("rbind",behav_info),by = c("species","colony"), relationship = "many-to-one") |>
  dplyr::mutate(sst_raw = (sst*sd_sst) + mean_sst) |>
  dplyr::mutate(
    colony = dplyr::case_when(
      colony %in% "Bjornoya" ~ "Bjørnøya",
      colony %in% "Hjelmsoya" ~ "Hjelmsøya",
      colony %in% "Hornoya" ~ "Hornøya",
      colony %in% "Rost" ~ "Røst",
      TRUE ~ colony
    )
  ) |>
  dplyr::mutate(colony = dplyr::case_when(
    colony == "Runde" & grepl("BlKi",species) ~  "Ålesund",
    colony == "Sklinna" ~ "Sør-Gjæslingan",
    TRUE ~ colony
  )) |>
  dplyr::mutate(species = ifelse(species == "AtPu", "Fratercula arctica", "Rissa tridactyla"))

dee_models <- lapply(grep("energy",model_files, value = TRUE),readRDS)
dee_info <- lapply(dee_models,`[[`,"sst_info")
dee_derivs <- lapply(dee_models,`[[`,"derivs") |>
  lapply(function(X) dplyr::mutate(X, median = median(estimate),.by = c(species,colony,bird1,sst)) |>
           subset(drawid %in% sample(unique(drawid), 50, replace = FALSE))) |>
  do.call("rbind", args = _) |>
  tibble::remove_rownames() |>
  dplyr::left_join(do.call("rbind",dee_info),by = c("species","colony"), relationship = "many-to-one") |>
  dplyr::mutate(sst_raw = (sst*sd_sst) + mean_sst) |>
  dplyr::mutate(
    colony = dplyr::case_when(
      colony %in% "Bjornoya" ~ "Bjørnøya",
      colony %in% "Hjelmsoya" ~ "Hjelmsøya",
      colony %in% "Hornoya" ~ "Hornøya",
      colony %in% "Rost" ~ "Røst",
      TRUE ~ colony
    )
  ) |>
  dplyr::mutate(colony = dplyr::case_when(
    colony == "Runde" & grepl("BlKi",species) ~  "Ålesund",
    colony == "Sklinna" ~ "Sør-Gjæslingan",
    TRUE ~ colony
  )) |>
  dplyr::mutate(species = ifelse(species == "AtPu", "Fratercula arctica", "Rissa tridactyla"))

##########################################################################################
#### Plot ####
##########################################################################################

p1 <- ggplot(behav_derivs |>
         subset(species == "Fratercula arctica"),
       aes(x = sst_raw, y = estimate)) +
  geom_line(aes(col = colony,group = interaction(bird1,drawid)),alpha = 0.2) +
  geom_line(aes(y = median, group = bird1), colour = "black") +
  facet_wrap(~colony) +
  ylab("Partial foraging effort\nfirst derivative") + xlab("Sea surface temperature (°C)") +
  guides(color = "none", fill = "none") +
  scale_color_manual(values = c("#5F4E6B","#B7EB86","#87CEEB","#EBE586","#EB8699","#5C6B4E","#969586","#4E636B"), breaks = c("Ålesund","Anda", "Hjelmsøya", "Hornøya", "Røst", "Runde","Bjørnøya", "Sør-Gjæslingan")) +  theme_classic() +
  theme_classic(base_size = 10) +
  theme(
    axis.line = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ))  +
  ggplot(behav_derivs |>
           subset(species == "Rissa tridactyla"),
         aes(x = sst_raw, y = estimate)) +
  geom_line(aes(col = colony,group = interaction(bird1,drawid)),alpha = 0.1, linewidth = 0.5) +
  geom_line(aes(y = median, group = bird1), colour = "black", linewidth = 0.5) +
  facet_wrap(~colony) +
  ylab("Partial foraging effort\nfirst derivative") + xlab("Sea surface temperature (°C)") +
  guides(color = "none", fill = "none") +
  scale_color_manual(values = c("#5F4E6B","#B7EB86","#87CEEB","#EBE586","#EB8699","#5C6B4E","#969586","#4E636B"), breaks = c("Ålesund","Anda", "Hjelmsøya", "Hornøya", "Røst", "Runde","Bjørnøya", "Sør-Gjæslingan")) +  theme_classic() +
  theme_classic(base_size = 10) +
  theme(
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ))  

p2 <- ggplot(dee_derivs |>
               subset(species == "Fratercula arctica"),
             aes(x = sst_raw, y = estimate)) +
  geom_line(aes(col = colony,group = interaction(bird1,drawid)),alpha = 0.2) +
  geom_line(aes(y = median, group = bird1), colour = "black") +
  facet_wrap(~colony) +
  ylab("Partial energy expenditure\nfirst derivative") + xlab("Sea surface temperature (°C)") +
  coord_cartesian(ylim = c(-1,1)) +
  guides(color = "none", fill = "none") +
  scale_color_manual(values = c("#5F4E6B","#B7EB86","#87CEEB","#EBE586","#EB8699","#5C6B4E","#969586","#4E636B"), breaks = c("Ålesund","Anda", "Hjelmsøya", "Hornøya", "Røst", "Runde","Bjørnøya", "Sør-Gjæslingan")) +  theme_classic() +
  theme_classic(base_size = 10) +
  theme(
    axis.line = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ))  +
  ggplot(dee_derivs |>
           subset(species == "Rissa tridactyla"),
         aes(x = sst_raw, y = estimate)) +
  geom_line(aes(col = colony,group = interaction(bird1,drawid)),alpha = 0.1, linewidth = 0.5) +
  geom_line(aes(y = median, group = bird1), colour = "black", linewidth = 0.5) +
  facet_wrap(~colony) +
  ylab("Partial energy expenditure\nfirst derivative") + xlab("Sea surface temperature (°C)") +
  coord_cartesian(ylim = c(-1,1)) +
  guides(color = "none", fill = "none") +
  scale_color_manual(values = c("#5F4E6B","#B7EB86","#87CEEB","#EBE586","#EB8699","#5C6B4E","#969586","#4E636B"), breaks = c("Ålesund","Anda", "Hjelmsøya", "Hornøya", "Røst", "Runde","Bjørnøya", "Sør-Gjæslingan")) +  theme_classic() +
  theme_classic(base_size = 10) +
  theme(
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ))  

ggsave("outputs/figure2.png",p1 / p2 & patchwork::plot_annotation(tag_levels = "A"),dpi = 300, width = 8, height = 6)
