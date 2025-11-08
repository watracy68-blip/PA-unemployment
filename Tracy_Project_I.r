#Code for downloading state level data from FRED 
library(quantmod)
library(lubridate)
series <- "PAUR"
end_date <- lubridate::ymd("2025-06-30")
start_date <- end_date %m-% years(30)

suppressWarnings(
  getSymbols(Symbols = series, src = "FRED",
             from = start_date, to = end_date, auto.assign = TRUE)
)

x <- get(series); colnames(x) <- "value"

t0 <- as.Date(head(index(x),1))
ts_start <- c(year(t0), month(t0))
y_ts <- ts(as.numeric(x$value), start = ts_start, frequency = 12)

plot(y_ts)

# ================================================================
# BOX–JENKINS: STEP 1 — IDENTIFICATION
# Purpose: Determine (p, d, q) and, if seasonal, (P, D, Q)[m]
# Input:  y_ts  (monthly ts built in your preceding code)
# ================================================================

# ---- A) Packages for Identification (tests, plots) ----
needed <- c("forecast", "tseries", "urca")
to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
library(forecast)  # ACF/PACF, seasonplot, ndiffs, nsdiffs
library(tseries)   # kpss.test(), adf.test()

# Ensure numeric, finite series
y_ts <- stats::ts(as.numeric(y_ts), start = start(y_ts), frequency = frequency(y_ts))
m   <- frequency(y_ts)  # 12 for monthly

# ================================================================
# 1) STATIONARITY CHECK
#    Visual inspection + statistical tests (ADF, KPSS)
# ================================================================
par(mfrow = c(2,2))
plot(y_ts, main = "Level series y_t (visual stationarity check)", ylab = "Percent", xlab = "")
Acf(y_ts, lag.max = 72, main = "ACF (level)")
Pacf(y_ts, lag.max = 72, main = "PACF (level)")

# Quick preview of first difference and seasonal difference (if m>1)
if (m > 1) {
  plot(diff(y_ts, lag = m), main = expression("Preview: seasonal difference " * nabla[m] * " y_t"), ylab = "", xlab = "")
} else {
  plot(diff(y_ts), main = expression("Preview: first difference " * nabla * " y_t"), ylab = "", xlab = "")
}
par(mfrow = c(1,1))

###########Stationarity testting

cat("\n--- Stationarity tests on LEVEL ---\n")
# KPSS (H0: stationarity). Small p-value -> evidence AGAINST stationarity (suggest differencing)
print(kpss.test(y_ts, null = "Level"))
# ADF test using R, H0: unit root (nonstationary)
print(adf.test(y_ts))  

# NOTE:
# Combine visual + tests:
#  • If KPSS rejects stationarity and ADF fails to reject unit root → include difference.
#  • Prefer the SMALLEST differencing that yields stationarity to avoid overdifferencing.

# ================================================================
# 2) AUTOCORRELATION & PARTIAL AUTOCORRELATION ANALYSIS (LEVEL)
#    Use ACF/PACF to sense AR vs MA structure, but final read is on the
#    stationarized (differenced) series in Section 4.
# ================================================================
par(mfrow = c(1,2))
Acf(y_ts, lag.max = 72, main = "ACF (level): AR/MA cues") #repeat PACF and ACF
Pacf(y_ts, lag.max = 72, main = "PACF (level): AR/MA cues")
par(mfrow = c(1,1))

# NOTE:
#  • AR(p): PACF tends to “cut off” near p; ACF tails off.
#  • MA(q): ACF tends to “cut off” near q; PACF tails off.
#  • Mixed ARMA: both ACF & PACF often tail off.

# ================================================================
# 3) SEASONALITY CHECK
#    Visual inspection + seasonal decomposition (STL)
# ================================================================


# STL decomposition (This R function is more robust than the decompose() function, 
#but it only allows for additive decomposition. Additive is reasonable for rates)
stl_fit <- try(stl(y_ts, s.window = "periodic"), silent = TRUE)
if (!inherits(stl_fit, "try-error")) plot(stl_fit, main = "STL decomposition (Trend/Seasonal/Remainder)")

if (m > 1) {
  # Month/season means pattern also visible via monthplot (base)
  # If this plot looks like a unit root for each month, include a seasonal difference D=1
  suppressWarnings(monthplot(y_ts, main = "Monthplot: within-year pattern", ylab = "Percent"))
}

# NOTE:
# Clear repeating within-year pattern ⇒ include seasonal terms in ARIMA:
# ARIMA(p,d,q)(P,D,Q)[m]. 
# For NSA unemployment rate, D often = 1; for SA rate, D often = 0.

# ================================================================
# 4) DIFFERENCING ORDER SELECTION
#    Choose seasonal D first, then nonseasonal d. Build working series w_t.
# ================================================================
D_hat <- 0 #Try your option here based on seasonal plots above
y_seas <- if (D_hat > 0) diff(y_ts, lag = m, differences = D_hat) else y_ts

d_hat <- 1 #Try your option here based on stationarity tests above
w_ts  <- if (d_hat > 0) diff(y_seas, differences = d_hat) else y_seas
cat(sprintf("\nTrial differencing orders: D = %d (seasonal), d = %d (nonseasonal)\n", D_hat, d_hat))

# Re-check stationarity on working series
cat("\n--- Stationarity tests on WORKING SERIES w_t ---\n")
print(kpss.test(w_ts, null = "Level"))
print(adf.test(w_ts))


# ACF/PACF of w_t to guide p, q (and seasonal P, Q via spikes at multiples of m)
par(mfrow = c(1,2))
Acf(w_ts, lag.max = 72, main = expression(paste("ACF of working series  ", w[t])))
Pacf(w_ts, lag.max = 72, main = expression(paste("PACF of working series  ", w[t])))
par(mfrow = c(1,1))

# Simple candidate ranges to carry into estimation step. Choose based on above ACF/PACF
p_candidates <- 0:3; q_candidates <- 0 # Set max value based on ACF/PACF
P_candidates <- 0  # Max is one if Seasonality present in ACF/PACF, zero otherwise
Q_candidates <- 0  # Max is one if Seasonality present in ACF/PACF, zero otherwise  

cat("\n--- Identification summary ---\n")
cat(sprintf("Period m = %d | D = %d | d = %d\n", m, D_hat, d_hat))
cat("Use ACF/PACF of w_t to bound:\n")
cat("  Nonseasonal: p in {0,1,2,3}, q in {0,1,2,3}\n")
if (m > 1) cat("  Seasonal:    P in {0,1,2},   Q in {0,1,2}\n")
cat("Proceed to ESTIMATION with a small grid using these bounds.\n")

# Object to pass forward
identification <- list(
  frequency = m,
  d = d_hat, D = D_hat,
  w_ts = w_ts, y_seas = y_seas,
  p_candidates = p_candidates,
  q_candidates = q_candidates,
  P_candidates = P_candidates,
  Q_candidates = Q_candidates
)
# Inspect structure
str(identification)


# ================================================================
# BOX–JENKINS: STEP 2 — ESTIMATION
# Purpose: Fit candidate ARIMA models using orders suggested by Identification,
#          compare by information criteria (AICc/BIC), and select a final model.
# Inputs expected:
#   - y_ts : original monthly series
#   - identification list from Step 1 with:
#       d, D, frequency=m, w_ts, y_seas,
#       p_candidates, q_candidates, P_candidates, Q_candidates
# ================================================================
###########Initial Code to Run

# ---- A) Setup Packages & helpers ----
needed <- c("forecast")
to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
library(forecast)

# Safe AICc getter (avoids "could not find function 'AICc'" errors)
get_AICc <- function(fit) {
  if ("AICc" %in% ls("package:forecast")) return(forecast::AICc(fit))
  if (exists("AICc")) return(AICc(fit))
  # Fallback formula: AICc = AIC + 2k(k+1)/(n - k - 1), k = df of logLik (incl. sigma^2)
  k <- attr(logLik(fit), "df")
  n <- fit$nobs
  AIC(fit) + (2 * k * (k + 1)) / (n - k - 1)
}

#Define coefficient table function

coef_table <- function(fit) {
  se <- sqrt(diag(fit$var.coef))
  est <- fit$coef
  tval <- est / se
  pval <- 2 * pt(-abs(tval), df = length(fit$residuals) - length(est))
  data.frame(Estimate = est, Std.Error = se, t = tval, `Pr(>|t|)` = pval, check.names = FALSE)
}

# ---- B) Pull Identification outputs (with safe fallbacks) ----
m <- frequency(y_ts)
if (!exists("identification")) {
  # Fallback: recompute minimal pieces if Step 1 object is missing
  D_hat <- if (m > 1) forecast::nsdiffs(y_ts, m = m) else 0
  y_seas <- if (D_hat > 0) diff(y_ts, lag = m, differences = D_hat) else y_ts
  d_hat <- forecast::ndiffs(y_seas)
  identification <- list(
    frequency = m, d = d_hat, D = D_hat,
    w_ts = if (d_hat > 0) diff(y_seas, differences = d_hat) else y_seas,
    y_seas = y_seas,
    p_candidates = 0:3, q_candidates = 0:3,
    P_candidates = if (m > 1) 0:2 else 0:0,
    Q_candidates = if (m > 1) 0:2 else 0:0
  )
}

d_hat <- identification$d
D_hat <- identification$D
p_grid <- identification$p_candidates
q_grid <- identification$q_candidates
P_grid <- identification$P_candidates
Q_grid <- identification$Q_candidates

# Choose information criterion for selection:
criterion <- "AICc"   # options: "AICc" or "BIC"

# Include mean / drift rules:
include_mean  <- (d_hat + D_hat) == 0
include_drift <- (d_hat + D_hat) > 0   # allow deterministic drift when differenced

# ================================================================
# (1) MODEL SELECTION
#     Build a small, interpretable grid over (p,q,P,Q), fit each candidate,
#     and rank by chosen information criterion (AICc default).
# ================================================================
cand_fits <- list()
scoreboard <- data.frame()

### This function estimates the values of the coefficients for all values in the grid
for (p in p_grid) for (q in q_grid) {
  for (P in P_grid) for (Q in Q_grid) {
    # Skip the totally empty model unless a mean/drift is present
    if ((p + q + P + Q) == 0 && !(include_mean || include_drift)) next
    fit <- try(
      Arima(y_ts,
            order = c(p, d_hat, q),
            seasonal = if (m > 1) list(order = c(P, D_hat, Q), period = m) else NULL,
            include.mean  = include_mean,
            include.drift = include_drift,
            method = "ML"),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    
    # Scores
    aicc <- get_AICc(fit)
    bic  <- BIC(fit)
    
    cand_fits[[length(cand_fits) + 1]] <- fit
    scoreboard <- rbind(scoreboard, data.frame(
      idx = length(cand_fits),
      p = p, d = d_hat, q = q, P = P, D = D_hat, Q = Q, m = m,
      AICc = aicc, BIC = bic,
      stringsAsFactors = FALSE
    ))
  }
}

stopifnot(nrow(scoreboard) > 0)

# Rank by the default criterion (AICc unless you changed `criterion`) and set label
ranking_table <- scoreboard[order(scoreboard[[criterion]]), ]
best_row      <- ranking_table[1, , drop = FALSE]
final_fit     <- cand_fits[[best_row$idx]]
selection_label <- criterion  # will be updated if an override is used

cat("\n=== MODEL SELECTION (top 10 by ", selection_label, ") ===\n", sep = "")
print(head(ranking_table, 10))


cat("\nChosen model (", selection_label, "):\n", sep = "")
print(best_row)

# NOTE:
# We deliberately use a SMALL grid bound by Identification to avoid overfitting.
# Selection principle: minimize information criterion (AICc default) to balance fit vs complexity.

# ================================================================
# (2) PARAMETER ESTIMATION
#     Once the model is specified, estimate AR/MA (and seasonal) parameters
#     via Maximum Likelihood (using the Arima() function). Report estimates & standard errors.
# ================================================================
cat("\n=== PARAMETER ESTIMATION: Coefficient table (t-stats & p-values) ===\n")
print(round(coef_table(final_fit), 4))

cat("\nModel summary (including logLik, AIC/AICc/BIC):\n")
print(final_fit)
cat("\nInformation criteria for the chosen model:\n")
cat(sprintf("AIC  = %.3f | AICc = %.3f | BIC = %.3f\n", AIC(final_fit), get_AICc(final_fit), BIC(final_fit)))

# NOTE:
# ML finds θ = (ϕ, θ, Φ, Θ, drift/mean, σ^2) maximizing log-likelihood.
# Standard errors come from the observed information matrix (inverse Hessian).

# ================================================================
# (3) MODEL FITTING
#     Use the estimated model to obtain fitted values; assess in-sample fit
#     (plots + simple accuracy metrics).
# ================================================================
# One-step-ahead fitted values (in-sample)
fitted_vals <- fitted(final_fit)
resids      <- residuals(final_fit)

# Accuracy metrics (in-sample)
rmse <- sqrt(mean(resids^2, na.rm = TRUE))
mae  <- mean(abs(resids), na.rm = TRUE)
mape <- mean(abs(resids / y_ts), na.rm = TRUE) * 100

cat("\n=== MODEL FITTING: In-sample accuracy ===\n")
cat(sprintf("RMSE = %.4f | MAE = %.4f | MAPE = %.2f%%\n", rmse, mae, mape))

# Plot actual vs fitted (one-step ahead estimate)
op <- par(no.readonly = TRUE)
par(mfrow = c(2,1))
plot(y_ts, main = "Actual vs Fitted (one-step ahead)", ylab = "Percent", xlab = "")
lines(fitted_vals, col = "blue")
legend("topleft", lty = c(1,1), col = c("black", "blue"), bty = "n", c("Actual", "Fitted"))

plot(resids, main = "Residuals from chosen model", ylab = "", xlab = "")
abline(h = 0, lty = 2)
par(op)

# NOTE:
# “Model fitting” here refers to applying the estimated parameters to produce
# fitted values and inspecting basic fit metrics. Formal residual diagnostics
# (residual analysis, Ljung-Box test) belong to Step 3 (Diagnostics).


#================================================================
# BOX–JENKINS: STEP 3 — DIAGNOSTICS
# Purpose: Validate that residuals from the chosen ARIMA are white noise:
#          no autocorrelation, mean ~ 0, roughly constant variance.
# Inputs expected:
#   - final_fit : the chosen Arima() model from Step 2
#   - y_ts      : original series (for plots/scale)
# ================================================================

# ---- A) Packages & helpers ----
needed <- c("forecast", "tseries")
to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
library(forecast)
library(tseries)

# Safe access to residuals/fitted
stopifnot(exists("final_fit"))
resids <- residuals(final_fit)
fits   <- fitted(final_fit)

# helper: count ARMA-type parameters for Ljung–Box df (exclude mean/drift/sigma)
arma_df <- {
  nm <- names(coef(final_fit))
  if (is.null(nm)) 0L else sum(grepl("ar|ma", nm, ignore.case = TRUE) &
                                 !grepl("mean|drift|intercept|sigma", nm, ignore.case = TRUE))
}

# ================================================================
# (1) RESIDUAL ANALYSIS
#     Visual inspection: no pattern over time, roughly constant variance,
#     residual ACF/PACF ~ 0
# ================================================================
op <- par(no.readonly = TRUE)
par(mfrow = c(2,2))
plot(resids, main = "Residuals over time", ylab = "", xlab = "")
abline(h = 0, lty = 2)

hist(resids, breaks = "FD", main = "Histogram of residuals", xlab = "Residual")
lines(density(na.omit(resids)), lwd = 2)

Acf(resids, lag.max = 72, main = "Residual ACF")
Pacf(resids, lag.max = 72, main = "Residual PACF")
par(op)

# NOTE:
# • White-noise residuals show no visible trend/seasonal pattern; ACF/PACF spikes lie within bands.

# ================================================================
# (2) LJUNG–BOX TEST
#     H0: no autocorrelation up to lag L (residuals are white noise).
#     We test multiple horizons (e.g., 12, 24, 36) and adjust df by ARMA params.
# ================================================================
LB_lags <- c(12, 24, 36)
cat("\n=== Ljung–Box tests (H0: no residual autocorrelation) ===\n")
for (L in LB_lags) {
  lb <- Box.test(resids, lag = L, type = "Ljung-Box", fitdf = arma_df)
  cat(sprintf("Lag %2d: X^2=%.3f, df=%d, p=%.4f\n", L, lb$statistic, lb$parameter, lb$p.value))
}

cat("\n--- Combined check (forecast::checkresiduals) ---\n")
# Produces residual plot + ACF + LB test in one view
suppressWarnings(checkresiduals(final_fit))

# NOTE:
# • Prefer p-values > 0.05 (no evidence of residual autocorrelation).
# • If p < 0.05 at seasonal multiples (e.g., 12, 24), consider seasonal AR/MA adjustments.


###############################################
# ARIMA Forecasting in R (Post Box–Jenkins)
# Series: PAUR (SA monthly PA unemployment rate)
# Window: 1995-07 to 2025-06
###############################################

# ---- 0) Package setup ---------------------------------------------------------
needed <- c("quantmod", "forecast", "lubridate", "zoo")
to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(quantmod)   # getSymbols() from FRED
library(forecast)   # Arima(), forecast(), autoplot(), accuracy()
library(lubridate)  # date handling
library(zoo)        # as.yearmon(), plotting helpers

# ---- 1) Download data from FRED -----------------------------------------------
getSymbols("PAUR", src = "FRED", auto.assign = TRUE, warnings = FALSE)

# ---- 2) Restrict to the study window and convert to ts ------------------------
start_date <- ymd("1995-07-01")
end_date   <- ymd("2025-06-30")

# xts subset by date
un_sa_xts <- window(PAUR, start = start_date, end = end_date)

# Convert to a monthly 'ts' object with frequency m = 12
# start=c(1995,7) means: year 1995, month 7 in ts indexing
y_ts <- ts(as.numeric(un_sa_xts),
           start = c(1995, 7),
           frequency = 12)

# Sanity check: basic plot of the level series
plot(y_ts, main = "Pennsylvania Unemployment Rate (SA)\nPAUR: Jul 1995 – Jun 2025",
     ylab = "Percent", xlab = "")

# ---- 3) (Assumed) Box–Jenkins has been done; plug in YOUR chosen orders -------
# Replace the next six integers with the orders you selected during identification
# and estimation. The seasonal period for monthly data is m = 12.

m <- frequency(y_ts)  # should be 12 for monthly

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# EDIT THESE to your chosen orders from Box–Jenkins:
p <- 2; d <- 1; q <- 0     # nonseasonal ARIMA(p,d,q)
P <- 0; D <- 0; Q <- 0     # seasonal ARIMA(P,D,Q)[m]
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Notes on means/drifts:
# - If d + D == 0 (no differencing at all), you typically include a mean (intercept).
# - If d + D == 1, you often include a drift (linear trend in levels).
# - If d + D >= 2, drift is usually not identifiable/stable in practice.

include_mean  <- (d + D == 0)
include_drift <- (d + D == 1)

# ---- 5) Forecasting -----------------------------------------------------------
# H = forecast horizon in months; choose as needed (e.g., 12, 24)
H <- 24

# level = confidence levels for prediction intervals
fc <- forecast(fit, h = H, level = c(80, 95))

# Console summary of the forecast
print(fc)

# ---- 6) Neat forecast table with calendar months ------------------------------
# Convert the ts index of fc$mean into "yearmon" and then a YYYY-MM Date
f_dates <- as.yearmon(time(fc$mean))            # "2025 Jul", "2025 Aug", ...
f_dates <- as.Date(f_dates)                     # becomes first day of each month
forecast_tbl <- data.frame(
  date  = f_dates,
  mean  = as.numeric(fc$mean),
  lo80  = as.numeric(fc$lower[,"80%"]),
  hi80  = as.numeric(fc$upper[,"80%"]),
  lo95  = as.numeric(fc$lower[,"95%"]),
  hi95  = as.numeric(fc$upper[,"95%"])
)
head(forecast_tbl, 6)
# You can write to CSV if desired:
# write.csv(forecast_tbl, "forecast_PAUR_ARIMA.csv", row.names = FALSE)

# ---- 7) Plot: full sample + forecast and PIs ---------------------------------
autoplot(fc) +
  ggplot2::labs(title = "ARIMA Forecast for Pennsylvania Unemployment (SA)",
                subtitle = sprintf("Model: (%d,%d,%d)(%d,%d,%d)[%d],  h = %d months",
                                   p,d,q,P,D,Q,m,H),
                x = "", y = "Percent")

# ---- 8) Zoomed plot: last N years + forecast ---------------------------------
N_years <- 10  # show last 10 years of history for context
start_zoom <- c(end(time(y_ts)) - (N_years - 1), 1) # rough; or use window()

autoplot(window(y_ts, start = time(y_ts)[length(y_ts) - N_years*m + 1])) +
  autolayer(fc$mean, series = "Forecast") +
  autolayer(fc$lower[,"95%"], series = "95% PI (lower)") +
  autolayer(fc$upper[,"95%"], series = "95% PI (upper)") +
  ggplot2::labs(title = "Zoomed: Recent History and ARIMA Forecast (SA)",
                x = "", y = "Percent")

# ---- 9) (Optional) Backtesting snippet (rolling-origin or holdout) -----------
# If you want to *verify* forecast skill out-of-sample, reserve the last k months
# as a test set (here, k = 18 for example) and compare forecasts to actuals.

k <- 18
y_train <- head(y_ts, length(y_ts) - k)
y_test  <- tail(y_ts, k)

fit_bt <- Arima(
  y = y_train,
  order    = c(p, d, q),
  seasonal = list(order = c(P, D, Q), period = m),
  include.mean  = include_mean,
  include.drift = include_drift,
  method  = "ML"
)

fc_bt <- forecast(fit_bt, h = k)
print(accuracy(fc_bt, y_test))  # MAE, RMSE, MAPE, etc.

autoplot(fc_bt) +
  autolayer(y_test, series = "Actual (Holdout)") +
  ggplot2::labs(title = "Holdout Backtest: ARIMA Forecast vs Actual",
                x = "", y = "Percent")

###############################################
# End of script
###############################################

