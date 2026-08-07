###############################################################
# 01. Two-stage analysis
# Stage 1: district-specific DLNM, time-stratified case-crossover, conditional quasi-Poisson regression
# Stage 2: multivariate meta-regression (mixmeta)
###############################################################

library(dplyr); library(tidyr); library(data.table); library(lubridate)
library(readxl); library(writexl); library(dlnm); library(mixmeta); library(gnm)
library(splines); library(scales); library(AER)

rm(list = ls())

###############################################################
# 1. Data and settings ##############
###############################################################

# Daily district-level time series, 2010-2023 (see README)
load("data/d_ts_cause_all_2010_2023.RData")

# Calendar variables
d_ts$year  <- year(d_ts$date)
d_ts$month <- month(d_ts$date)
d_ts$dow   <- wday(d_ts$date)
d_ts$doy   <- yday(d_ts$date)

# Warm season only
d_ts <- d_ts %>% filter(month(date) %in% 5:10)

# Model settings
varper   <- c(50, 90)          # internal knots for the exposure dimension
varfun   <- "ns"
lag      <- 5
lagknots <- logknots(lag, 1)   # internal knot for the lag dimension

# District list
sgginfo <- sort(unique(d_ts$sgg_h_v2019))

# Outcomes: cause x age group (primary diagnosis only)
causelist <- c("i", "j", "n", "f", "AB", "E", "G")
agelist   <- c("u65", "6574", "o75")
sublist   <- paste0("n_icd_",
                    rep(causelist, each = length(agelist)), "_", agelist)

# Temperature and humidity correlation (reported in the appendix)
cor.test(d_ts$tmean, d_ts$shum)

# Containers
coef_list <- list()
meta_list <- list()

# Loop for each subgroup and district
for (k in 1:length(sublist)) {

  varsel <- sublist[k]
  print(varsel)

  # Container for coefficients and (vectorised) covariances
  coef_tab <- data.frame(sgg_h_v2019 = sgginfo)
  ncoef <- 3
  coef_tab[, paste0("coef", seq(ncoef))] <- NA
  coef_tab[, paste0("vcov", seq(ncoef * (ncoef + 1) / 2))] <- NA

  ###############################
  # 1-1. First stage
  ###############################

  for (i in 1:nrow(coef_tab)) {

    sgg_temp <- coef_tab[i, "sgg_h_v2019"]

    d_ts_temp <- d_ts %>% filter(sgg_h_v2019 == sgg_temp) %>% as.data.frame()
    print(paste0("Sgg: ", sgg_temp))

    d_ts_temp$n_sel <- d_ts_temp[, varsel]

    # Case-crossover strata: year x month x day of week
    d_ts_temp$stratum <- factor(paste(d_ts_temp$year, d_ts_temp$month,
                                      d_ts_temp$dow, sep = ":"))
    d_ts_temp <- d_ts_temp %>% mutate(keep = sum(n_sel) > 0, .by = stratum)

    # Cross-basis (year groups prevent lags crossing seasons)
    group <- factor(paste(d_ts_temp$sgg_h_v2019, d_ts_temp$year, sep = "-"))
    cb <- crossbasis(d_ts_temp$tmean,
                     lag = lag,
                     argvar = list(fun = varfun,
                                   knots = quantile(d_ts_temp$tmean, varper / 100, na.rm = TRUE)),
                     arglag = list(knots = lagknots),
                     group = group)

    # Model and reduction to the overall cumulative association
    model <- gnm(n_sel ~ cb + ns(shum, df = 3) + factor(holiday),
                 data = d_ts_temp, family = quasipoisson,
                 subset = keep, eliminate = stratum)
    red <- crossreduce(cb, model)
    coef_tab[i, -c(1)] <- t(c(coef(red), vechMat(vcov(red))))

  }

  coef_list[[varsel]] <- coef_tab

  ###############################
  # 1-2. Second stage
  ###############################

  coef <- coef_tab[, grep("coef", colnames(coef_tab))] %>% as.matrix()
  vcov <- coef_tab[, grep("vcov", colnames(coef_tab))] %>% as.matrix()

  # Pooled association (used for projection)
  meta1 <- mixmeta(coef, vcov, random = ~ 1 | sgg_h_v2019, data = coef_tab)

  meta_list[[varsel]] <- meta1

  ###############################
  # 1-3. Pooled exposure-response curve
  ###############################

  # Nationwide average of the district-specific percentiles
  pred <- seq(1, 99, by = 0.5)
  perc <- do.call(rbind, tapply(d_ts$tmean, d_ts$sgg_h_v2019,
                                function(x) quantile(x, pred / 100)))
  perc <- apply(perc, 2, mean)

  b1 <- onebasis(perc, fun = varfun,
                 knots = perc[paste0(sprintf("%.1f", varper), "%")],
                 Bound = c(min(perc), max(perc)))
  cen <- perc["75.0%"]
  cp <- crosspred(b1, coef = coef(meta1), vcov = vcov(meta1),
                  model.link = "log", cen = cen)

  plot(cp, main = varsel, ylim = c(0.9, 1.2),
       xlab = expression(paste("Temperature ("*degree, "C)")),
       ylab = "RR (95% CI)",
       col = "black", lwd = 1.6,
       ci.arg = list(col = alpha("black", 0.05)))
  abline(v = perc["75.0%"], lty = 3, col = "grey50")
  abline(v = perc["99.0%"], lty = 3, col = "red")
}

dev.off()

###############################################################
# 2. Export ##############
###############################################################

save(coef_list, meta_list, dp_list,
     file = "output/result_gnm.RData")

# Pooled coefficients and covariance matrix, one sheet per outcome
# (read by 02_projection.R)
meta_export <- lapply(meta_list, function(m) {
  data.frame(coef = coef(m), vcov(m), check.names = FALSE)
})
write_xlsx(meta_export, "output/meta.xlsx")

