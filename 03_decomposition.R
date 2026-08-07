###############################################################
# 03. Shapley decomposition of the projected change in cost
# Factors: temperature change (heat), population ageing (age), population size (pop), cost per admission (cpc)
###############################################################

library(dplyr); library(ggplot2); library(ggh4x); library(gridExtra)

rm(list = ls())

scnlist <- c("ssp245", "ssp585")

fac_df <- expand.grid(heat = c(FALSE, TRUE), age = c(FALSE, TRUE),
                      pop = c(FALSE, TRUE), cpc = c(FALSE, TRUE))
state_lab <- apply(fac_df, 1, function(x) {
  nm <- c("heat", "age", "pop", "cpc")[as.logical(x)]
  if (length(nm) == 0) "baseline" else paste(nm, collapse = "_")
})

# KRW to USD (2023 average)
krw_usd <- 1305.6625

get_yrs <- function(x) ifelse(x == "2020", 4, ifelse(x == "2024", 6, 5))
baseline_periods <- c("2010", "2015", "2020")
future_periods   <- c("2024", "2030", "2035")
cost_factors     <- c("heat", "pop", "age", "cpc")
full_state_name  <- "heat_age_pop_cpc"

###############################################################
# 1. Population denominator ##############
###############################################################

load("data/KOSIS_pop_proj.RData")

# Observed
popdf <- poplist[["obs"]]
popdf <- popdf %>% group_by(year) %>%
  summarise(pop_tot = sum(pop_tot), pop_u65 = sum(pop_u65),
            pop_o65 = sum(pop_o65), pop_6574 = sum(pop_6574),
            pop_o75 = sum(pop_o75))

# Projected
popdf2 <- poplist[["med"]]
popdf2 <- popdf2 %>% group_by(year) %>%
  summarise(pop_tot = sum(pop_tot), pop_u65 = sum(pop_u65),
            pop_o65 = sum(pop_o65), pop_6574 = sum(pop_6574),
            pop_o75 = sum(pop_o75))
popdf2 <- popdf2 %>% filter(!(year %in% c(2023)))

popdf <- rbind(popdf, popdf2) %>% filter(year %in% 2010:2039)

popdf <- popdf %>% mutate(
  period = case_when(
    year %in% 2010:2014 ~ "2010",
    year %in% 2015:2019 ~ "2015",
    year %in% 2020:2023 ~ "2020",
    year %in% 2024:2029 ~ "2024",
    year %in% 2030:2034 ~ "2030",
    year %in% 2035:2039 ~ "2035"))

popdf <- popdf %>% group_by(period) %>%
  summarise(pop_tot = sum(pop_tot), pop_u65 = sum(pop_u65),
            pop_o65 = sum(pop_o65), pop_6574 = sum(pop_6574),
            pop_o75 = sum(pop_o75)) %>% ungroup()

# Person-years by period and age group
py_df <- popdf %>%
  tidyr::pivot_longer(cols = c(pop_tot, pop_u65, pop_o65, pop_6574, pop_o75),
                      names_to = "age", values_to = "pop") %>%
  mutate(age = recode(age,
                      "pop_tot" = "total", "pop_u65" = "u65",
                      "pop_o65" = "o65", "pop_6574" = "6574",
                      "pop_o75" = "o75"))

###############################################################
# 2. Merge projection output ##############
###############################################################

# Containers: [scenario][GCM][state]
# GCM-level arrays are kept throughout; simulations of all GCMs are pooled
# at the summary stage rather than averaged
full_ansim_mod_list    <- setNames(lapply(scnlist, function(s) list()), scnlist)
full_costsim_mod_list  <- setNames(lapply(scnlist, function(s) list()), scnlist)
full_doy_n_mod_list    <- setNames(lapply(scnlist, function(s) list()), scnlist)
full_doy_cost_mod_list <- setNames(lapply(scnlist, function(s) list()), scnlist)

file_list <- list.files("output", pattern = "^proj_.*\\.RData$", full.names = TRUE)

for (f in file_list) {

  rm(list = intersect(c("ansim_list", "costsim_list", "doy_n_list", "doy_cost_list"), ls()))

  load(f)
  cat("Loaded:", f, "\n")

  for (scn in names(costsim_list)) {
    for (mod in names(costsim_list[[scn]])) {

      keep_state <- names(costsim_list[[scn]][[mod]])[
        !vapply(costsim_list[[scn]][[mod]], is.null, logical(1))]
      if (length(keep_state) == 0) next

      if (is.null(full_costsim_mod_list[[scn]][[mod]])) {
        full_ansim_mod_list[[scn]][[mod]]    <- setNames(vector("list", length(state_lab)), state_lab)
        full_costsim_mod_list[[scn]][[mod]]  <- setNames(vector("list", length(state_lab)), state_lab)
        full_doy_n_mod_list[[scn]][[mod]]    <- setNames(vector("list", length(state_lab)), state_lab)
        full_doy_cost_mod_list[[scn]][[mod]] <- setNames(vector("list", length(state_lab)), state_lab)
      }

      full_ansim_mod_list[[scn]][[mod]][keep_state]    <- ansim_list[[scn]][[mod]][keep_state]
      full_costsim_mod_list[[scn]][[mod]][keep_state]  <- costsim_list[[scn]][[mod]][keep_state]
      full_doy_n_mod_list[[scn]][[mod]][keep_state]    <- doy_n_list[[scn]][[mod]][keep_state]
      full_doy_cost_mod_list[[scn]][[mod]][keep_state] <- doy_cost_list[[scn]][[mod]][keep_state]
    }
  }
  gc()
}

rm(list = intersect(c("ansim_list", "costsim_list", "doy_n_list", "doy_cost_list"), ls()))
gc()

# Number of GCMs and of states retrieved per GCM
sapply(full_costsim_mod_list, length)
sapply(full_costsim_mod_list, function(x) sum(!vapply(x[[1]], is.null, logical(1))))

# The GCM-level lists are large (3 scenarios x 19 GCMs x 16 states);
# drop full_ansim_mod_list and full_doy_n_mod_list if memory is limited

# GCM order, kept identical across states so that the pooled simulation
# vectors remain aligned
mod_names <- names(full_costsim_mod_list[[scnlist[1]]])

# Dimension names (identical across scenarios, GCMs and states)
arr_tmp      <- full_costsim_mod_list[[scnlist[1]]][[1]][[full_state_name]]
period_names <- dimnames(arr_tmp)[[1]]
cause_names  <- dimnames(arr_tmp)[[2]]
sim_names    <- dimnames(arr_tmp)[[4]]
sim_idx      <- setdiff(seq_along(sim_names), which(sim_names == "est"))
rm(arr_tmp)

###############################################################
# 3. Shapley decomposition ##############
###############################################################

# State name in canonical factor order
normalize_state <- function(x, factors = c("heat", "age", "pop", "cpc")) {
  if (is.null(x) || x == "" || x == "baseline") return("baseline")
  sp <- unlist(strsplit(x, "_", fixed = TRUE))
  sp <- intersect(factors, sp)
  sp <- factors[factors %in% sp]
  if (length(sp) == 0) "baseline" else paste(sp, collapse = "_")
}

subset_to_state <- function(subset_vec, factors = c("heat", "age", "pop", "cpc")) {
  subset_vec <- intersect(factors, subset_vec)
  subset_vec <- factors[factors %in% subset_vec]
  if (length(subset_vec) == 0) "baseline" else paste(subset_vec, collapse = "_")
}

# Shapley value of each factor: weighted mean marginal contribution
# over all subsets of the remaining factors
decomp_one_scn <- function(state_list, factors = c("heat", "age", "pop", "cpc")) {

  state_list <- state_list[!vapply(state_list, is.null, logical(1))]

  names(state_list) <- vapply(names(state_list), normalize_state,
                              character(1), factors = factors)
  state_list <- state_list[!duplicated(names(state_list), fromLast = TRUE)]

  all_states <- unlist(lapply(0:length(factors), function(k) {
    combn(factors, k, FUN = function(z) subset_to_state(z, factors), simplify = TRUE)
  }))
  all_states[all_states == ""] <- "baseline"
  all_states <- unique(all_states)

  miss_states <- setdiff(all_states, names(state_list))
  if (length(miss_states) > 0) {
    stop(paste("Missing states:", paste(miss_states, collapse = ", ")))
  }

  template <- state_list[[1]] * 0
  out <- setNames(vector("list", length(factors)), factors)
  for (ff in factors) out[[ff]] <- template

  p <- length(factors)

  for (ff in factors) {

    others <- setdiff(factors, ff)

    for (k in 0:length(others)) {

      subs_k <- combn(others, k, simplify = FALSE)
      if (k == 0) subs_k <- list(character(0))

      for (S in subs_k) {
        state_S  <- subset_to_state(S, factors)
        state_Si <- subset_to_state(c(S, ff), factors)

        w <- factorial(length(S)) * factorial(p - length(S) - 1) / factorial(p)
        out[[ff]] <- out[[ff]] + w * (state_list[[state_Si]] - state_list[[state_S]])
      }
    }
  }

  out
}

# Decomposition is run separately for each GCM
decomp_costsim_list <- setNames(vector("list", length(scnlist)), scnlist)
for (scn in scnlist) {
  decomp_costsim_list[[scn]] <- setNames(vector("list", length(mod_names)), mod_names)
  for (mod in mod_names) {
    decomp_costsim_list[[scn]][[mod]] <- decomp_one_scn(full_costsim_mod_list[[scn]][[mod]])
    # Net change, stored alongside the four factors
    decomp_costsim_list[[scn]][[mod]][["net"]] <-
      Reduce(`+`, decomp_costsim_list[[scn]][[mod]][cost_factors])
  }
}

###############################################################
# 4. Change relative to the baseline period ##############
###############################################################

# Cost unit: million USD (2023 prices); rate: per 1,000 person-years

# The Shapley values give the difference between the full state and the baseline state within the same period. 
# To express the change relative to the baseline period, 
# the same-period baseline state is added back and the baseline-period absolute value is subtracted:
#   full_state(p) = baseline_state(p) + sum of Shapley values
#   change = full_state(p) - baseline reference (14 years)
# When splitting by factor, baseline_state(p) and the baseline reference are multiplied by w = 1 / n_fac 
# so that each is counted exactly once in the sum across factors.

cause_ext  <- c(cause_names, "all")   # "all" pools the seven causes
factor_ext <- c(cost_factors, "net")  # "net" is the sum of the four factors

p_base   <- which(period_names %in% baseline_periods)
yrs_base <- sum(get_yrs(period_names[p_base]))
py_base  <- py_df %>%
  filter(period %in% baseline_periods, age == "total") %>%
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>% pull(pop)

# Baseline reference: full state averaged over the SSP scenarios, by GCM
arr_base_mean_mod <- setNames(vector("list", length(mod_names)), mod_names)
for (mod in mod_names) {
  arrs <- lapply(scnlist, function(s) full_costsim_mod_list[[s]][[mod]][[full_state_name]])
  arr_base_mean_mod[[mod]] <- Reduce(`+`, arrs) / length(arrs)
}

# Pooled baseline simulations by cause (19 GCMs x 1,000 samples)
base_sim_list <- setNames(vector("list", length(cause_ext)), cause_ext)
for (cc in cause_ext) {
  c_idx <- if (cc == "all") seq_along(cause_names) else which(cause_names == cc)
  v <- numeric(0)
  for (mod in mod_names) {
    v <- c(v, apply(arr_base_mean_mod[[mod]][p_base, c_idx, , sim_idx, drop = FALSE],
                    4, sum, na.rm = TRUE))
  }
  base_sim_list[[cc]] <- v
}

# Single loop over scenario, cause, period and factor
decomp_change_df <- data.frame()
n_fac <- length(cost_factors)

for (scn in scnlist) {
  for (cc in cause_ext) {

    c_idx <- if (cc == "all") seq_along(cause_names) else which(cause_names == cc)
    base_sim_total <- base_sim_list[[cc]]

    for (p in future_periods) {

      p_idx <- which(period_names == p)
      if (length(p_idx) == 0) next
      yrs <- get_yrs(p)

      fut_py <- py_df %>%
        filter(period == p, age == "total") %>%
        summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>% pull(pop)
      if (length(fut_py) == 0 || is.na(fut_py) || fut_py == 0) next

      # Same-period baseline state, pooled over GCMs (shared by all factors)
      fut_base_state_sim <- numeric(0)
      for (mod in mod_names) {
        fut_base_state_sim <- c(fut_base_state_sim,
          apply(full_costsim_mod_list[[scn]][[mod]][["baseline"]][p_idx, c_idx, , sim_idx, drop = FALSE],
                4, sum, na.rm = TRUE))
      }

      for (fac in factor_ext) {

        w <- if (fac == "net") 1 else 1 / n_fac

        fut_change_sim <- numeric(0)
        for (mod in mod_names) {
          fut_change_sim <- c(fut_change_sim,
            apply(decomp_costsim_list[[scn]][[mod]][[fac]][p_idx, c_idx, , sim_idx, drop = FALSE],
                  4, sum, na.rm = TRUE))
        }

        n_sim <- min(length(fut_base_state_sim), length(fut_change_sim), length(base_sim_total))
        if (n_sim == 0) next

        fut_sim_total <- fut_base_state_sim[1:n_sim] * w + fut_change_sim[1:n_sim]

        sim_diff_total  <- fut_sim_total - base_sim_total[1:n_sim] * w
        sim_diff_annual <- (fut_sim_total / yrs) - (base_sim_total[1:n_sim] * w / yrs_base)
        sim_diff_rate   <- (fut_sim_total / fut_py) - (base_sim_total[1:n_sim] * w / py_base)

        tmp <- data.frame(
          scenario = scn, factor = fac, period = p, cause = cc, age = "total",

          est_total = mean(sim_diff_total, na.rm = TRUE) / 1000,
          lci_total = unname(quantile(sim_diff_total, 0.025, na.rm = TRUE)) / 1000,
          uci_total = unname(quantile(sim_diff_total, 0.975, na.rm = TRUE)) / 1000,

          est_annual = mean(sim_diff_annual, na.rm = TRUE) / 1000,
          lci_annual = unname(quantile(sim_diff_annual, 0.025, na.rm = TRUE)) / 1000,
          uci_annual = unname(quantile(sim_diff_annual, 0.975, na.rm = TRUE)) / 1000,

          est_py = mean(sim_diff_rate, na.rm = TRUE) * 1000000 / 1000,
          lci_py = unname(quantile(sim_diff_rate, 0.025, na.rm = TRUE)) * 1000000 / 1000,
          uci_py = unname(quantile(sim_diff_rate, 0.975, na.rm = TRUE)) * 1000000 / 1000)

        decomp_change_df <- rbind(decomp_change_df, tmp)
      }
    }
  }
}

rownames(decomp_change_df) <- NULL

decomp_change_df <- decomp_change_df %>%
  mutate(across(matches("est|lci|uci"), ~ . / krw_usd))

write.csv(decomp_change_df, "output/decomp_change.csv", row.names = FALSE)

###############################################################
# 5. Figure: stacked factor contributions ##############
###############################################################

cause_labels <- c(all = "Total", i = "Circulatory", j = "Respiratory",
                  n = "Genitourinary", f = "Mental", AB = "Infectious",
                  E = "Endocrine & metabolic", G = "Nervous")
factor_cols <- c(heat = "#d73027", age = "#91bfdb", pop = "#4575b4",
                 cpc = "#fc8d59", net = "grey50")
period_labels <- c("2024" = "2024-29", "2030" = "2030-34", "2035" = "2035-39")
cause_order <- c("all", "i", "j", "n", "f", "AB", "E", "G")

# Stacked bars: the four factors; error bars: eCI of the net change
bar_all <- decomp_change_df %>% filter(factor != "net") %>%
  mutate(cause  = factor(cause, levels = cause_order),
         factor = factor(factor, levels = c("heat", "age", "pop", "cpc")))
ci_all  <- decomp_change_df %>% filter(factor == "net") %>%
  mutate(cause = factor(cause, levels = cause_order))

# The supplementary figure of the rate per 1,000 person-years uses the same
# code with est_py, lci_py and uci_py in place of the annual columns
p1 <- ggplot() +
  geom_col(data = bar_all,
           aes(x = period, y = est_annual, fill = factor),
           width = 0.7, color = "grey20", linewidth = 0.2) +
  geom_errorbar(data = ci_all,
                aes(x = period, ymin = lci_annual, ymax = uci_annual),
                width = 0.2, linewidth = 0.3, color = "grey20") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = factor_cols,
                    breaks = c("heat", "age", "pop", "cpc"),
                    labels = c("Temperature change", "Population ageing",
                               "Population size", "Cost per admission")) +
  scale_x_discrete(labels = period_labels) +
  ggh4x::facet_grid2(cause ~ scenario, scales = "free_y", independent = "y",
                     labeller = labeller(cause = cause_labels)) +
  labs(x = NULL,
       y = "Change in annual heat-attributable ED admission costs\n(million USD)",
       fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("output/fig_decomp_annual.png", p1, width = 10, height = 14)
