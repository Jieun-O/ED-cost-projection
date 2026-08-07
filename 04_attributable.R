###############################################################
# 04. Attributable costs (AC) and attributable fractions (AF)
###############################################################

library(dplyr); library(ggplot2); library(ggh4x); library(gridExtra)

full_state_name  <- "heat_age_pop_cpc"
baseline_periods <- c("2010", "2015", "2020")
get_yrs <- function(x) ifelse(x == "2020", 4, ifelse(x == "2024", 6, 5))
krw_usd <- 1305.6625

###############################################################
# 1. Baseline as a pseudo-scenario ##############
###############################################################

# Baseline results are reported as the mean of the three SSP scenarios
full_costsim_mod_list[["baseline"]] <- setNames(lapply(mod_names, function(m) {
  arrs <- lapply(scnlist, function(s) full_costsim_mod_list[[s]][[m]][[full_state_name]])
  setNames(list(Reduce(`+`, arrs) / length(arrs)), full_state_name)
}), mod_names)

full_doy_cost_mod_list[["baseline"]] <- setNames(lapply(mod_names, function(m) {
  arrs <- lapply(scnlist, function(s) full_doy_cost_mod_list[[s]][[m]][[full_state_name]])
  setNames(list(Reduce(`+`, arrs) / length(arrs)), full_state_name)
}), mod_names)

# Scenario and period combinations
task_df <- rbind(
  expand.grid(scenario = scnlist, period = period_names, stringsAsFactors = FALSE),
  expand.grid(scenario = "baseline", period = c(baseline_periods, "baseline"),
              stringsAsFactors = FALSE))

cause_ext <- c(cause_names, "all")   # "all" pools the seven causes

###############################################################
# 2. Attributable cost and attributable fraction ##############
###############################################################

# Costs are summed over age groups
af_period_df <- data.frame()

for (t in 1:nrow(task_df)) {

  scn <- task_df$scenario[t]
  p   <- task_df$period[t]

  p_sel <- if (p == "baseline") baseline_periods else p
  p_idx <- which(period_names %in% p_sel)
  if (length(p_idx) == 0) next
  yrs <- sum(get_yrs(period_names[p_idx]))

  py_val <- py_df %>%
    filter(period %in% p_sel, age == "total") %>%
    summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>% pull(pop)
  if (length(py_val) == 0 || is.na(py_val) || py_val == 0) next

  for (cc in cause_ext) {

    c_idx <- if (cc == "all") seq_along(cause_names) else which(cause_names == cc)

    sim_total <- numeric(0)
    sim_af    <- numeric(0)
    doy_total <- numeric(0)

    for (mod in mod_names) {

      arr0 <- full_costsim_mod_list[[scn]][[mod]][[full_state_name]]
      darr <- full_doy_cost_mod_list[[scn]][[mod]][[full_state_name]]

      # Attributable cost, 1,000 simulations for this GCM
      sim_mod <- apply(arr0[p_idx, c_idx, , sim_idx, drop = FALSE], 4, sum, na.rm = TRUE)

      # Total cost for this GCM (no simulation dimension)
      doy_mod <- sum(darr[p_idx, c_idx, , drop = FALSE], na.rm = TRUE)

      sim_total <- c(sim_total, sim_mod)
      sim_af    <- c(sim_af, sim_mod / doy_mod * 100)
      doy_total <- c(doy_total, doy_mod)
    }

    sim_annual <- sim_total / yrs
    sim_rate   <- sim_total / py_val

    tmp <- data.frame(
      scenario = scn, period = p, cause = cc,

      est_total = mean(sim_total, na.rm = TRUE) / 1000,
      lci_total = unname(quantile(sim_total, 0.025, na.rm = TRUE)) / 1000,
      uci_total = unname(quantile(sim_total, 0.975, na.rm = TRUE)) / 1000,

      est_annual = mean(sim_annual, na.rm = TRUE) / 1000,
      lci_annual = unname(quantile(sim_annual, 0.025, na.rm = TRUE)) / 1000,
      uci_annual = unname(quantile(sim_annual, 0.975, na.rm = TRUE)) / 1000,

      est_py = mean(sim_rate, na.rm = TRUE) * 1000000 / 1000,
      lci_py = unname(quantile(sim_rate, 0.025, na.rm = TRUE)) * 1000000 / 1000,
      uci_py = unname(quantile(sim_rate, 0.975, na.rm = TRUE)) * 1000000 / 1000,

      # Total cost: GCM mean
      doy_total  = mean(doy_total, na.rm = TRUE) / 1000,
      doy_annual = mean(doy_total, na.rm = TRUE) / yrs / 1000,

      est_af = mean(sim_af, na.rm = TRUE),
      lci_af = unname(quantile(sim_af, 0.025, na.rm = TRUE)),
      uci_af = unname(quantile(sim_af, 0.975, na.rm = TRUE)))

    af_period_df <- rbind(af_period_df, tmp)
  }
}

rownames(af_period_df) <- NULL

# AF is a ratio and is already scale free
af_period_df <- af_period_df %>%
  mutate(across(c(matches("est_total|lci_total|uci_total|annual|_py"), starts_with("doy")),
                ~ . / krw_usd))

write.csv(af_period_df, "output/attributable_cost_fraction.csv", row.names = FALSE)

###############################################################
# 3. Figure: attributable cost and attributable fraction  #####
###############################################################

cause_labels <- c(all = "Total", i = "Circulatory", j = "Respiratory",
                  n = "Genitourinary", f = "Mental", AB = "Infectious",
                  E = "Endocrine & metabolic", G = "Nervous")
cause_order <- c("all", "i", "j", "n", "f", "AB", "E", "G")

plot_bar <- af_period_df %>%
  filter((scenario == "baseline" & period == "baseline") |
           (scenario %in% scnlist & period %in% c("2024", "2030", "2035"))) %>%
  mutate(grp   = factor(scenario, levels = c("baseline", scnlist)),
         cause = factor(cause, levels = cause_order),
         x = case_when(period == "baseline" ~ 1.3,
                       period == "2024" ~ 2.2,
                       period == "2030" ~ 3.2,
                       period == "2035" ~ 4.2))

fill_vals <- c(baseline = "grey50", ssp245 = "#f4a3a3", ssp585 = "#d73027")
fill_labs <- c(baseline = "Baseline", ssp245 = "SSP2-4.5", ssp585 = "SSP5-8.5")
x_breaks  <- c(1.3, 2.2, 3.2, 4.2)
x_labs    <- c("2010-23", "2024-29", "2030-34", "2035-39")

# Panel A: annual cost, panel B: attributable fraction,
# supplementary panel: rate per 1,000 person-years
metric_list <- c("annual", "af", "py")
ylab_list <- c(
  annual = "Annual heat-attributable ED admission costs\n(million USD)",
  af     = "Heat-attributable fraction of ED admission costs\n(%)",
  py     = "Heat-attributable ED admission costs\n(USD per 1,000 person-years)")
title_list <- c(annual = "A", af = "B", py = "")

p_list <- list()

for (mt in metric_list) {

  yv <- paste0("est_", mt); lv <- paste0("lci_", mt); uv <- paste0("uci_", mt)

  p_list[[mt]] <- ggplot(plot_bar, aes(x = x, y = .data[[yv]], fill = grp)) +

    geom_col(data = plot_bar %>% filter(grp == "baseline"),
             width = 0.3, color = "grey20", linewidth = 0.25) +
    geom_errorbar(data = plot_bar %>% filter(grp == "baseline"),
                  aes(ymin = .data[[lv]], ymax = .data[[uv]]),
                  position = position_dodge(width = 0.7),
                  width = 0.15, linewidth = 0.3, color = "grey20") +

    geom_col(data = plot_bar %>% filter(grp != "baseline"),
             position = position_dodge(width = 0.7),
             width = 0.55, color = "grey20", linewidth = 0.25) +
    geom_errorbar(data = plot_bar %>% filter(grp != "baseline"),
                  aes(ymin = .data[[lv]], ymax = .data[[uv]]),
                  position = position_dodge(width = 0.7),
                  width = 0.25, linewidth = 0.3, color = "grey20") +

    geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = fill_vals, labels = fill_labs) +
    ggh4x::facet_wrap2(~ cause, scales = "free_y",
                       labeller = labeller(cause = cause_labels),
                       axes = "x", remove_labels = "none", ncol = 4) +
    scale_x_continuous(breaks = x_breaks, labels = x_labs) +
    labs(title = title_list[[mt]], x = NULL, y = ylab_list[[mt]], fill = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_cartesian(ylim = c(0, NA))
}

pall <- grid.arrange(p_list[["annual"]] + theme(legend.position = "none"),
                     p_list[["af"]], ncol = 1, heights = c(1, 1.1))

ggsave("output/fig_ac_af_proj.png", pall, width = 12.5, height = 10)
ggsave("output/suppl_fig_ac_py.png", p_list[["py"]], width = 12.5, height = 6)
