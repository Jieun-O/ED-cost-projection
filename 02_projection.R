###############################################################
# 02. Projection of attributable admissions and costs
# Loop: SSP scenario > GCM > decomposition state > cause > district
###############################################################

library(dplyr); library(readxl); library(dlnm); library(mvmeta)
library(MASS); library(lubridate); library(data.table)
library(splines); library(stringr)

rm(list = ls())

###############################################################
# 1. Settings ##############
###############################################################

# Age groups
age_grps_sel <- c("u65", "6574", "o75")

# Ranges for running the script in parallel sessions
scn_strt_num   <- 1; scn_end_num   <- 2   # SSP scenario
mod_strt_num   <- 1; mod_end_num   <- 19  # GCM
state_strt_num <- 1; state_end_num <- 16  # decomposition state

# Causes
sublist <- c("i", "j", "n", "f", "AB", "E", "G")
subname <- c("Circulatory", "Respiratory", "Genitourinary", "Mental",
             "Infectious", "Endocrine & metabolic", "Nervous")
subdf <- data.frame(sublist, subname)
subdf$var <- subdf$sublist
subdf$lab <- subdf$subname

# SSP scenarios
scnlist <- c("ssp245", "ssp585")

# Periods (label = first year of the period; 2024 period covers 2024-2029)
periodinfo <- c(2010, 2015, 2020, 2024, 2030, 2035)

# Decomposition states: 2^4 combinations of the four factors
fac_df <- expand.grid(heat = c(FALSE, TRUE), age = c(FALSE, TRUE),
                      pop = c(FALSE, TRUE), cpc = c(FALSE, TRUE))
state_lab <- apply(fac_df, 1, function(x) {
  nm <- c("heat", "age", "pop", "cpc")[as.logical(x)]
  if (length(nm) == 0) "baseline" else paste(nm, collapse = "_")
})

# Cost per admission trend used in the main analysis
cpc_sel <- "lin1023"   # alternatives: "linex21", "damp"

# Exposure-response settings (identical to 01_twostage.R)
varper   <- c(50, 90)
varfun   <- "ns"
lag      <- 5
lagknots <- logknots(lag, 1)

# GCM list, taken from the projection file names
modlist <- list.files("data", pattern = "^tmean_scn_245_cal_.*\\.RData$")
modlist <- sub("^tmean_scn_245_cal_", "", sub("\\.RData$", "", modlist))

###############################################################
# 2. Observed temperature ##############
###############################################################

load("data/era5_noaa_kor_1023.RData")
de_ts <- do.call(rbind, weather_list_sgg); rm(list = c("weather_list_sgg"))
de_ts <- de_ts %>% filter(substr(date, 1, 4) %in% c(2010:2023))
de_ts$date <- as.Date(de_ts$date)
de_ts$mday <- format(de_ts$date, "%m-%d")
de_ts <- de_ts %>% filter(month(date) %in% 5:10)
de_ts$sgg_h_v2019 <- de_ts$region

sgginfo <- sort(unique(de_ts$sgg_h_v2019))

# District list for faster subsetting
de_ts_split <- split(de_ts, de_ts$sgg_h_v2019)

# District code key
sggkey <- read_xlsx("data/sgg_info_bh_v2019.xlsx")
sggkey <- sggkey %>% filter(year == 2019)

###############################################################
# 3. Admission, cost and population inputs ##############
###############################################################

load("data/doy_list.RData")

# (1) District-level admissions by cause and age group (for admission rates)
temp_list <- list()
for (i in 1:length(doy_list)) {
  temp <- doy_list[[i]] %>%
    group_by(sgg_h_v2019, icd, sub) %>%
    summarise(n = sum(n, na.rm = TRUE),
              cost = sum(cost, na.rm = TRUE),
              .groups = "drop")
  temp_list[[i]] <- temp
}
temp_df <- do.call(rbind, temp_list)
temp_df$sub <- paste0("icd_", temp_df$icd, "_", temp_df$sub)
prop <- temp_df

# (2) Day-of-year admission and cost profile by district
temp_list <- list()
for (i in 1:length(doy_list)) {
  temp <- doy_list[[i]] %>%
    group_by(month, sgg_h_v2019, icd, sub) %>%
    summarise(n = mean(n, na.rm = TRUE),
              cost = mean(cost, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(cost_pc = cost / n) %>%
    tidyr::pivot_wider(names_from = sub,
                       values_from = c(n, cost, cost_pc),
                       names_glue = "{.value}_{sub}") %>%
    arrange(month)
  temp_list[[i]] <- temp
}
temp_df <- do.call(rbind, temp_list)

month_days <- data.table(month = 1:12,
                         days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
temp_df <- merge(temp_df, month_days, by = "month", all.x = TRUE)

# Monthly totals to daily values, then expanded to day of year
temp_df <- temp_df %>%
  mutate(n_u65        = n_u65 / days,
         n_6574       = n_6574 / days,
         n_o75        = n_o75 / days,
         cost_u65     = cost_u65 / days,
         cost_6574    = cost_6574 / days,
         cost_o75     = cost_o75 / days,
         cost_pc_u65  = cost_u65 / n_u65,
         cost_pc_6574 = cost_6574 / n_6574,
         cost_pc_o75  = cost_o75 / n_o75) %>%
  tidyr::uncount(weights = days, .id = "day_in_month")
temp_df <- temp_df %>%
  mutate(date = as.Date(sprintf("2020-%02d-%02d", month, day_in_month)),
         doy  = yday(date)) %>%
  arrange(sgg_h_v2019, icd) %>%
  dplyr::select(-date, -day_in_month)
doy_list <- split(temp_df, temp_df$icd)

# (3) Nationwide monthly cost per admission
cost_list <- read_xlsx("data/doy_mth_kor.xlsx", guess_max = 50000)
cost_list <- cost_list %>% filter(sub %in% c(age_grps_sel))
cost_list <- cost_list %>%
  group_by(month, icd, sub) %>%
  summarise(n = sum(n, na.rm = TRUE),
            cost = sum(cost, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(cost_pc = cost / n) %>%
  tidyr::pivot_wider(names_from = sub,
                     values_from = c(n, cost, cost_pc),
                     names_glue = "{.value}_{sub}") %>%
  arrange(month)
cost_list <- split(cost_list, cost_list$icd)

rm(list = c("temp", "temp_df", "temp_list"))

# (4) Projected cost per admission
load("data/cost_est.RData")

# (5) Projected population (medium variant)
load("data/pop_est.RData")
popdf <- popestlist[["med"]]
popdf <- popdf %>% filter(year %in% 2010:2039)

# (6) Pooled coefficients and covariance matrices from 01_twostage.R
for (ag in age_grps_sel) {
  assign(paste0("meta_coef_", ag, "_list"),
         setNames(lapply(subdf$var, function(v) {
           read_xlsx("output/meta.xlsx", sheet = paste0("n_icd_", v, "_", ag))[, 1] %>%
             t() %>% as.vector()
         }), subdf$var))

  assign(paste0("meta_vcov_", ag, "_list"),
         setNames(lapply(subdf$var, function(v) {
           as.matrix(read_xlsx("output/meta.xlsx",
                               sheet = paste0("n_icd_", v, "_", ag))[, 2:4]) %>% t()
         }), subdf$var))
}

###############################################################
# 4. Output containers ##############
###############################################################

# [district, period, cause, age, simulation]
ansim <- array(NA,
               dim = c(length(sgginfo), length(periodinfo), length(subdf$var),
                       length(age_grps_sel), 1000 + 1),
               dimnames = list(sgginfo, periodinfo, subdf$var, age_grps_sel,
                               c("est", paste0("sim", 1:1000))))
costsim <- ansim

# Total (not attributable) admissions and costs, used as the AF denominator
doy_n <- array(NA,
               dim = c(length(sgginfo), length(periodinfo), length(subdf$var),
                       length(age_grps_sel)),
               dimnames = list(sgginfo, periodinfo, subdf$var, age_grps_sel))
doy_cost <- doy_n

# Templates for re-initialising within the loop
ansim_df    <- ansim
costsim_df  <- costsim
doy_n_df    <- doy_n
doy_cost_df <- doy_cost

# Results by SSP scenario > GCM > decomposition state
empty_state <- function() setNames(vector("list", length(unique(state_lab))), unique(state_lab))
empty_mod   <- function() setNames(lapply(modlist, function(m) empty_state()), modlist)
ansim_list    <- setNames(lapply(scnlist, function(s) empty_mod()), scnlist)
costsim_list  <- setNames(lapply(scnlist, function(s) empty_mod()), scnlist)
doy_n_list    <- setNames(lapply(scnlist, function(s) empty_mod()), scnlist)
doy_cost_list <- setNames(lapply(scnlist, function(s) empty_mod()), scnlist)

###############################################################
# 5. Projection ##############
###############################################################

# 5-1. SSP scenario
for (scn in scn_strt_num:scn_end_num) {

  scn_temp <- scnlist[[scn]]
  print(scn_temp)

  # 5-2. GCM
  for (m in mod_strt_num:mod_end_num) {

    mod_temp <- modlist[[m]]
    print(mod_temp)

    # Bias-corrected daily temperature for this GCM
    load(paste0("data/tmean_scn_", substr(scn_temp, 4, 7), "_cal_", mod_temp, ".RData"))

    idx_ws <- (year(tmean_scn_mod$date) %in% 2010:2039) &
      (month(tmean_scn_mod$date) %in% 5:10)
    tmean_scn_mod_split <- split(tmean_scn_mod[idx_ws, ], tmean_scn_mod$sgg_h[idx_ws])
    gc()

    # 5-3. Decomposition state
    for (es in state_strt_num:state_end_num) {

      state_temp <- state_lab[[es]]
      print(state_temp)

      # Re-initialise arrays for each state
      ansim    <- ansim_df
      costsim  <- costsim_df
      doy_n    <- doy_n_df
      doy_cost <- doy_cost_df

      # National totals (district dimension dropped)
      ansim_kor    <- array(0, dim = dim(ansim)[-1],    dimnames = dimnames(ansim)[-1])
      costsim_kor  <- array(0, dim = dim(costsim)[-1],  dimnames = dimnames(costsim)[-1])
      doy_n_kor    <- array(0, dim = dim(doy_n)[-1],    dimnames = dimnames(doy_n)[-1])
      doy_cost_kor <- array(0, dim = dim(doy_cost)[-1], dimnames = dimnames(doy_cost)[-1])

      # 5-4. Cause
      for (i in 1:nrow(subdf)) {

        sub_temp <- subdf[i, "var"]

        dt_doy  <- doy_list[[sub_temp]]
        dt_cost <- cost_list[[sub_temp]]

        # Replace district cost per admission with the national monthly average
        dt_doy[c("cost_u65", "cost_6574", "cost_o75",
                 "cost_pc_u65", "cost_pc_6574", "cost_pc_o75")] <- NULL
        dt_cost[c("n_u65", "n_6574", "n_o75", "icd")] <- NULL
        dt_doy <- dt_doy %>% left_join(dt_cost, by = "month")

        dt_doy_split <- split(dt_doy[dt_doy$doy %in% 122:305, ],
                              dt_doy$sgg_h_v2019[dt_doy$doy %in% 122:305])

        # Coefficient draws by age group
        set.seed(1234)
        for (ag in age_grps_sel) {
          assign(paste0("mv_", ag), list(
            coef = get(paste0("meta_coef_", ag, "_list"))[[sub_temp]],
            vcov = get(paste0("meta_vcov_", ag, "_list"))[[sub_temp]]
          ))
          mv <- get(paste0("mv_", ag))
          assign(paste0("coefsim_", ag), mvrnorm(1000, mv$coef, mv$vcov))
        }

        # 5-5. District
        for (j in 1:length(sgginfo)) {

          sgg_temp <- sgginfo[j]

          d_ts_temp         <- de_ts_split[[as.character(sgg_temp)]]          # observed
          tmean_scn_mod_sgg <- tmean_scn_mod_split[[as.character(sgg_temp)]]  # projected
          tmean_scn_mod_sgg$mday <- format(tmean_scn_mod_sgg$date, "%m-%d")
          dt_doy_temp       <- dt_doy_split[[as.character(sgg_temp)]]

          # District-specific percentiles from the observed series
          pred <- seq(1, 99, by = 0.5)
          perc <- quantile(d_ts_temp$tmean, probs = pred / 100, na.rm = TRUE)

          # Factor "heat" off: future temperature replaced by the GCM baseline
          # day-of-year mean
          tmean_temp <- tmean_scn_mod_sgg[year(tmean_scn_mod_sgg$date) %in% 2010:2023, ]
          tmean_scn_mod_sgg$tmean_scn <- if (str_detect(state_temp, "heat") == FALSE) {
            rep(tapply(tmean_temp$tmean, tmean_temp$mday, function(x) mean(x, na.rm = TRUE)),
                length = 184 * length(unique(year(tmean_scn_mod_sgg$date))))
          } else tmean_scn_mod_sgg$tmean

          # Basis centred at the 75th percentile
          b1 <- onebasis(tmean_scn_mod_sgg$tmean_scn, fun = varfun,
                         knots = perc[paste0(sprintf("%.1f", varper), "%")],
                         Bound = c(min(perc), max(perc)))
          cen <- perc["75.0%"]
          ind <- match(c("fun", names(formals(attr(b1, "fun")))),
                       names(attributes(b1)), nomatch = 0)
          cenvec <- do.call("onebasis", c(list(x = cen), attributes(b1)[ind]))
          b1cen  <- scale(b1, center = cenvec, scale = FALSE)

          # Period label (184 warm-season days per year)
          years     <- c(5, 5, 4, 6, 5, 5)
          periodlab <- rep(periodinfo, times = years * 184)

          # Only days above the centring value contribute
          keep <- tmean_scn_mod_sgg$tmean_scn >= perc["75.0%"]
          g    <- as.integer(factor(periodlab[keep], levels = periodinfo))

          # District-level inputs (shared across age groups)
          cpcdf_temp      <- cpcdf %>% filter(var == sub_temp)
          cpc_growth_temp <- cpc_growth_sum_df %>% filter(var == sub_temp)
          popdf_temp      <- popdf %>% filter(sgg_h == sgg_temp)

          # Population scenario implied by the decomposition state
          scn_pick <- NA
          if (!(state_temp %in% c("baseline", "heat"))) {
            scn_pick <- if (grepl("pop", state_temp) & grepl("age", state_temp)) "pop+age"
            else if (grepl("pop", state_temp)) "pop"
            else if (grepl("age", state_temp)) "age"
            else NA
          }

          # 5-6. Age group
          for (ag in age_grps_sel) {

            doy_n   <- rep(dt_doy_temp[[paste0("n_", ag)]],       length = nrow(tmean_scn_mod_sgg))
            doy_cpc <- rep(dt_doy_temp[[paste0("cost_pc_", ag)]], length = nrow(tmean_scn_mod_sgg))

            # Baseline admission rate (2010-2023 total / 14)
            n_ag <- prop %>%
              filter(sub == paste0("icd_", sub_temp, "_", ag), sgg_h_v2019 == sgg_temp) %>%
              dplyr::select(n) %>% pull() / 14
            stopifnot(length(n_ag) == 1)
            ir_ag <- n_ag / popdf_temp[[paste0("pop_", ag)]][1]

            # Factors "pop", "age", "cpc" on
            if (!(state_temp %in% c("baseline", "heat"))) {
              if (!is.na(scn_pick)) {
                popdf_temp_scn <- popdf_temp %>% filter(scn_dec == scn_pick)
                rate_n <- (ir_ag * rep(popdf_temp_scn[[paste0("pop_", ag)]], each = 184)) / n_ag
                doy_n  <- doy_n * rate_n
              }
              if (grepl("cpc", state_temp)) {
                ratio   <- cpcdf_temp[[paste0("cpc_ratio_", ag, "_", cpc_sel)]]
                doy_cpc <- doy_cpc * rep(as.numeric(ratio), each = 184)
              }
            }

            mv      <- get(paste0("mv_", ag))
            coefsim <- get(paste0("coefsim_", ag))

            # Point estimate, negative values truncated to zero
            an      <- pmax((1 - exp(-b1cen %*% mv$coef)) * doy_n, 0)
            an_cost <- pmax(an * doy_cpc, 0)

            ansim_kor[,   sub_temp, ag, 1] <- ansim_kor[,   sub_temp, ag, 1] +
              tapply(an[keep], periodlab[keep], sum)
            costsim_kor[, sub_temp, ag, 1] <- costsim_kor[, sub_temp, ag, 1] +
              tapply(an_cost[keep], periodlab[keep], sum)

            doy_n_kor[,    sub_temp, ag] <- doy_n_kor[,    sub_temp, ag] +
              tapply(doy_n[keep], periodlab[keep], sum)
            doy_cost_kor[, sub_temp, ag] <- doy_cost_kor[, sub_temp, ag] +
              tapply((doy_n * doy_cpc)[keep], periodlab[keep], sum)

            # Monte Carlo uncertainty (1,000 draws at once)
            eta_mat  <- b1cen %*% t(coefsim)
            an_mat   <- pmax((1 - exp(-eta_mat)) * doy_n, 0)
            cost_mat <- pmax(an_mat * doy_cpc, 0)

            an_sum   <- rowsum(an_mat[keep, , drop = FALSE],   group = g, reorder = FALSE)
            cost_sum <- rowsum(cost_mat[keep, , drop = FALSE], group = g, reorder = FALSE)

            ansim_kor[,   sub_temp, ag, 2:1001] <- ansim_kor[,   sub_temp, ag, 2:1001] + an_sum
            costsim_kor[, sub_temp, ag, 2:1001] <- costsim_kor[, sub_temp, ag, 2:1001] + cost_sum
          }
        }
      }

      ansim_list[[scn_temp]][[mod_temp]][[state_temp]]    <- ansim_kor
      costsim_list[[scn_temp]][[mod_temp]][[state_temp]]  <- costsim_kor
      doy_n_list[[scn_temp]][[mod_temp]][[state_temp]]    <- doy_n_kor
      doy_cost_list[[scn_temp]][[mod_temp]][[state_temp]] <- doy_cost_kor

      rm(list = c("ansim_kor", "costsim_kor", "doy_n_kor", "doy_cost_kor"))
      gc()
    }

    rm(list = c("tmean_scn_mod", "tmean_scn_mod_split"))
    gc()
  }

  save(ansim_list, costsim_list, doy_n_list, doy_cost_list,
       file = paste0("output/proj_", scn_temp,
                     "_mod", mod_strt_num, "_", mod_end_num,
                     "_dec", state_strt_num, "_", state_end_num, ".RData"))
}
