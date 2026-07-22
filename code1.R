# =====================================================================
# Realized Volatility Forecasting -- Walk-Forward Study (R), v2
# Organised into functions: acquisition / forecasting / evaluation /
# visualisation / discussion. Fixes applied per code review.
#   install.packages(c("quantmod","forecast","tseries","ggplot2"))
#
# METHOD NOTE (intended equivalence -- NOT yet verified against this R run):
#   fitted(Arima(rv, model = fit0)) is INTENDED to reproduce Python's
#   append(refit=FALSE).predict() (fixed-parameter one-step-ahead forecasts).
#   The equivalence was checked in Python on a synthetic case, but the R and
#   Python outputs have not been numerically compared. Treat parity as an
#   intended design, and confirm it by running validate_parity.py on the two
#   wf_results.csv files (see README). Do not cite parity as established
#   until that comparison passes.
# =====================================================================
suppressPackageStartupMessages({
  library(quantmod); library(forecast); library(tseries); library(ggplot2)
})
set.seed(42)

CFG <- list(
  indices     = c(NIFTY50 = "^NSEI", SENSEX = "^BSESN", BANKNIFTY = "^NSEBANK"),
  windows     = c(5, 10, 20, 30),
  start_date  = "2015-01-01",
  main_frac   = 0.60,
  train_fracs = c(0.60, 0.70, 0.80),
  IC          = "bic",
  n_boot      = 1000,
  run_sens    = TRUE,
  boot_all    = FALSE   # TRUE = bootstrap CIs for every cell (slow); rep cell always done
)

#acquisition
load_prices <- function(cfg) {
  out <- list()
  for (nm in names(cfg$indices)) {
    px <- tryCatch({
      s <- getSymbols(cfg$indices[[nm]], src = "yahoo",
                      from = cfg$start_date, auto.assign = FALSE)
      na.omit(Cl(s))
    }, error = function(e) { message("download failed: ", nm); NULL })
    if (!is.null(px) && nrow(px) > 0) out[[nm]] <- px
  }
  out
}

realized_vol <- function(price, w) {
  ret <- na.omit(diff(log(price)))
  as.numeric(na.omit(runSD(ret, n = w)))
}

# memoized RV: compute each (index, window) once and reuse (#5, #6)
.RV_CACHE <- new.env(parent = emptyenv())
rv_of <- function(price, nm, w) {
  key <- paste(nm, w)
  if (is.null(.RV_CACHE[[key]])) .RV_CACHE[[key]] <- realized_vol(price, w)
  .RV_CACHE[[key]]
}

#stationarity
stationarity_decision <- function(x, alpha = 0.05) {
  adf_p  <- suppressWarnings(adf.test(x)$p.value)
  kpss_p <- suppressWarnings(kpss.test(x)$p.value)
  a <- adf_p < alpha; k <- kpss_p > alpha
  if (a && k)        list(adf_p, kpss_p, "stationary", 0)
  else if (!a && !k) list(adf_p, kpss_p, "unit-root", 1)
  else               list(adf_p, kpss_p, "inconclusive", NA)
}

build_labels <- function(price_list, cfg) {
  labs <- list()
  for (nm in names(price_list)) for (w in cfg$windows)
    labs[[paste(nm, w)]] <- stationarity_decision(rv_of(price_list[[nm]], nm, w))
  labs
}

#forecasting
smape <- function(a, f) mean(2 * abs(a - f) / (abs(a) + abs(f)), na.rm = TRUE) * 100

fit_arima_ic <- function(x, d, ic) {
  args <- list(y = x, ic = ic, seasonal = FALSE)
  if (!is.na(d)) args$d <- d
  do.call(auto.arima, args)
}

apply_onestep <- function(model, rv, ntr) {
  # fixed-parameter one-step-ahead forecasts over the test region
  as.numeric(fitted(Arima(rv, model = model)))[(ntr + 1):length(rv)]
}

wf_forecasts <- function(rv, ntr, d, ic) {
  n <- length(rv); train <- rv[1:ntr]
  tryCatch({
    fit0    <- fit_arima_ic(train, d, ic)
    arima_o <- apply_onestep(fit0, rv, ntr)
    ses_o   <- as.numeric(fitted(ets(rv, model = ses(train,  h = 1)$model)))[(ntr + 1):n]
    holt_o  <- as.numeric(fitted(ets(rv, model = holt(train, h = 1)$model)))[(ntr + 1):n]
    naive_o <- rv[ntr:(n - 1)]
    # internal checks: right length, no missing forecasts (#10)
    nt <- n - ntr
    stopifnot(length(arima_o) == nt, length(ses_o) == nt,
              length(holt_o) == nt, length(naive_o) == nt,
              !anyNA(arima_o), !anyNA(ses_o), !anyNA(holt_o), !anyNA(naive_o))
    list(ok = TRUE, fit0 = fit0, order = arimaorder(fit0),
         Naive = naive_o, SES = ses_o, Holt = holt_o, ARIMA = arima_o)
  }, error = function(e) { message("wf fit failed: ", conditionMessage(e)); list(ok = FALSE) })
}

#evaluation
metrics <- function(a, f, scale) c(
  RMSE = sqrt(mean((a - f)^2)), MAE = mean(abs(a - f)),
  MASE = mean(abs(a - f)) / scale, sMAPE = smape(a, f))

dm_beats_naive <- function(e_m, e_n)
  tryCatch(dm.test(e_m, e_n, alternative = "less", h = 1, power = 2)$p.value,
           error = function(e) NA)

# bootstrap CI for ALL four metrics, resampling (actual, forecast) blocks
block_bootstrap_ci <- function(a, f, scale, n_boot, alpha = 0.05) {
  ok <- !is.na(a) & !is.na(f); a <- a[ok]; f <- f[ok]
  n <- length(a); L <- max(1, round(n^(1/3))); nb <- ceiling(n / L)
  M <- matrix(NA, n_boot, 4)
  for (b in 1:n_boot) {
    st <- sample(1:(n - L + 1), nb, replace = TRUE)
    idx <- unlist(lapply(st, function(s) s:(s + L - 1)))[1:n]
    M[b, ] <- metrics(a[idx], f[idx], scale)
  }
  ci <- apply(M, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
  colnames(ci) <- c("RMSE", "MAE", "MASE", "sMAPE"); ci
}

run_grid <- function(price_list, labels, frac, cfg, keep = FALSE, compare_ic = FALSE) {
  res <- list(); errs <- list(); acts <- list(); orders <- list(); sel_t <- 0
  for (nm in names(price_list)) for (w in cfg$windows) {
    rv <- rv_of(price_list[[nm]], nm, w)
    ntr <- floor(length(rv) * frac); n <- length(rv)
    train <- rv[1:ntr]; test <- rv[(ntr + 1):n]; scale <- mean(abs(diff(train)))
    lab <- labels[[paste(nm, w)]]
    t0 <- Sys.time(); wf <- wf_forecasts(rv, ntr, lab[[4]], cfg$IC)
    sel_t <- sel_t + as.numeric(Sys.time() - t0, units = "secs")
    if (!isTRUE(wf$ok)) next
    orders[[paste(nm, w)]] <- wf
    # optional AIC-vs-BIC forecast comparison for THIS cell (#12, global)
    aic_mase <- NA; aic_ord <- NA
    if (compare_ic) {
      fa <- tryCatch(fit_arima_ic(train, lab[[4]], "aicc"), error = function(e) NULL)
      if (!is.null(fa)) {
        aic_ord <- paste(arimaorder(fa), collapse = ",")
        aic_mase <- mean(abs(test - apply_onestep(fa, rv, ntr))) / scale
      }
    }
    for (m in c("Naive", "SES", "Holt", "ARIMA")) {
      f <- wf[[m]]
      if (keep) { errs[[paste(nm, w, m)]] <- test - f; acts[[paste(nm, w, m)]] <- test }
      mv <- metrics(test, f, scale)
      res[[length(res) + 1]] <- data.frame(
        index = nm, window = w, method = m,
        order = ifelse(m == "ARIMA", paste(wf$order, collapse = ","), ""),
        RMSE = mv["RMSE"], MAE = mv["MAE"], MASE = mv["MASE"], sMAPE = mv["sMAPE"],
        adf_p = lab[[1]], kpss_p = lab[[2]], stationarity = lab[[3]],
        aic_order = ifelse(m == "ARIMA", aic_ord, NA),
        aic_MASE  = ifelse(m == "ARIMA", aic_mase, NA),
        row.names = NULL, stringsAsFactors = FALSE)
    }
  }
  list(results = do.call(rbind, res), errors = errs, actuals = acts,
       orders = orders, sel_time = sel_t)
}

#representative
analyze_representative <- function(price_list, labels, results, errors, actuals,
                                   orders_by_cell, dm_sig, arima_ok, cfg) {
  arima_res <- results[results$method == "ARIMA", ]
  if (length(arima_ok) > 0) {
    cand <- arima_res[paste(arima_res$index, arima_res$window) %in% arima_ok, ]
    rr <- cand[which.min(cand$MASE), ]; reason <- "lowest MASE among DM-significant cells"
  } else { rr <- arima_res[which.min(arima_res$MASE), ]; reason <- "lowest MASE (none beat naive)" }
  ri <- rr$index; rw <- rr$window
  cat(sprintf("\nRepresentative: %s %dd, ARIMA(%s), MASE=%.3f -- %s\n",
              ri, rw, rr$order, rr$MASE, reason))
  
  rv <- rv_of(price_list[[ri]], ri, rw); n <- length(rv)
  ntr <- floor(n * cfg$main_frac); train <- rv[1:ntr]; test <- rv[(ntr + 1):n]
  scale <- mean(abs(diff(train)))
  rep_fit <- orders_by_cell[[paste(ri, rw)]]$fit0
  rep_d <- labels[[paste(ri, rw)]][[4]]
  Lag <- as.integer(max(10, min(40, 10 * log10(n))))
  
  png("acf_pacf.png", 1100, 450); par(mfrow = c(1, 2))
  Acf(rv, lag.max = Lag, main = sprintf("%s RV ACF (%dd)", ri, rw))
  Pacf(rv, lag.max = Lag, main = "RV PACF"); dev.off()
  
  # stationarity plot: show first difference ONLY when d = 1 (#11)
  show_diff <- identical(rep_d, 1)
  png("stationarity_series.png", 1100, ifelse(show_diff, 600, 350))
  par(mfrow = c(ifelse(show_diff, 2, 1), 1))
  plot(rv, type = "l", col = "darkred", main = paste(ri, "RV level"))
  if (show_diff) plot(diff(rv), type = "l", col = "steelblue", main = "First difference (d=1)")
  dev.off()
  
  # model-based one-step PI via rolling forecast() (#2,#3,#6)
  m <- l <- u <- numeric(n - ntr)
  for (i in 1:(n - ntr)) {
    fc <- forecast(Arima(rv[1:(ntr + i - 1)], model = rep_fit), h = 1, level = 95)
    m[i] <- fc$mean; l[i] <- fc$lower[1]; u[i] <- fc$upper[1]
  }
  coverage <- mean(test >= l & test <= u)
  png("representative_forecast_interval.png", 1100, 600)
  plot(test, type = "l", lwd = 2,
       main = sprintf("%s %dd one-step forecast + 95%% PI (coverage=%.1f%%)",
                      ri, rw, 100 * coverage))
  lines(m, col = "blue"); lines(l, col = "blue", lty = 2); lines(u, col = "blue", lty = 2)
  dev.off()
  
  # residual diagnostics: checkresiduals() ONLY for the plot; compute the
  # Ljung-Box p-value EXPLICITLY (its return value is version-dependent) (#2,#5)
  png("residual_diagnostics.png", 1100, 500); checkresiduals(rep_fit, plot = TRUE); dev.off()
  res_r <- na.omit(residuals(rep_fit))
  pq <- sum(arimaorder(rep_fit)[c(1, 3)])          # p + q (fitted df)
  lb_p <- Box.test(res_r, lag = max(10, pq + 1), type = "Ljung-Box", fitdf = pq)$p.value
  jb_p <- suppressWarnings(jarque.bera.test(res_r)$p.value)
  
  # error viz for ALL methods, histograms not density (#9,#10)
  png("error_histograms.png", 1100, 800); par(mfrow = c(2, 2))
  for (mm in c("Naive", "SES", "Holt", "ARIMA"))
    hist(errors[[paste(ri, rw, mm)]], breaks = 40, col = "grey",
         main = paste(mm, "errors"), xlab = "error")
  dev.off()
  roll <- function(e, k = max(10, length(e) %/% 15))
    sqrt(as.numeric(stats::filter(e^2, rep(1/k, k), sides = 1)))
  methods4 <- c("Naive", "SES", "Holt", "ARIMA")
  rr_list <- lapply(methods4, function(mm) roll(errors[[paste(ri, rw, mm)]]))
  ymax <- max(unlist(rr_list), na.rm = TRUE)        # scale to rolling RMSE (#7)
  png("rolling_rmse.png", 1000, 500)
  plot(NULL, xlim = c(1, length(test)), ylim = c(0, ymax * 1.05),
       xlab = "test step", ylab = "rolling RMSE",
       main = "Rolling RMSE by method (leading gap = warm-up window)")   # (#8)
  for (j in seq_along(methods4)) lines(rr_list[[j]], col = j)
  legend("topright", methods4, col = seq_along(methods4), lty = 1); dev.off()
  
  # bootstrap CIs for ALL methods and ALL four metrics (#7,#8)
  boot <- lapply(c("Naive","SES","Holt","ARIMA"), function(mm)
    block_bootstrap_ci(actuals[[paste(ri, rw, mm)]],
                       actuals[[paste(ri, rw, mm)]] - errors[[paste(ri, rw, mm)]],
                       scale, cfg$n_boot))
  names(boot) <- c("Naive","SES","Holt","ARIMA")
  
  list(ri = ri, rw = rw, order = rr$order, lb_p = lb_p, jb_p = jb_p,
       coverage = coverage, boot = boot)
}

#main
main <- function(cfg = CFG) {
  price_list <- load_prices(cfg)
  labels <- build_labels(price_list, cfg)
  
  G <- run_grid(price_list, labels, cfg$main_frac, cfg, keep = TRUE, compare_ic = TRUE)
  results <- G$results; errors <- G$errors; actuals <- G$actuals
  orders_by_cell <- G$orders
  write.csv(results, "wf_results.csv", row.names = FALSE)
  cat(sprintf("ARIMA selection time: %.1fs over %d cells\n",
              G$sel_time, length(price_list) * length(cfg$windows)))
  print(results[, c("index","window","method","order","RMSE","MASE","sMAPE","stationarity")],
        row.names = FALSE)
  
  n_combos <- length(price_list) * length(cfg$windows)
  dm_sig <- c(SES = 0, Holt = 0, ARIMA = 0); arima_ok <- c()
  for (nm in names(price_list)) for (w in cfg$windows) for (m in c("SES","Holt","ARIMA")) {
    p <- dm_beats_naive(errors[[paste(nm,w,m)]], errors[[paste(nm,w,"Naive")]])
    if (!is.na(p) && p < 0.05) { dm_sig[m] <- dm_sig[m] + 1
    if (m == "ARIMA") arima_ok <- c(arima_ok, paste(nm, w)) }
  }
  
  R <- analyze_representative(price_list, labels, results, errors, actuals,
                              orders_by_cell, dm_sig, arima_ok, cfg)
  
  # full-grid bootstrap CIs -> CSV (optional; slow) (#4,#7,#8)
  # The representative cell's CIs (all methods, all 4 metrics) are always
  # produced in analyze_representative(); this block adds every other cell.
  if (cfg$boot_all) {
    boot_rows <- list()
    for (nm in names(price_list)) for (w in cfg$windows) {
      rvw <- rv_of(price_list[[nm]], nm, w)                 # computed once (#5)
      sc <- mean(abs(diff(rvw[1:floor(length(rvw) * cfg$main_frac)])))
      for (m in c("Naive","SES","Holt","ARIMA")) {
        key <- paste(nm, w, m); a <- actuals[[key]]; e <- errors[[key]]
        ci <- block_bootstrap_ci(a, a - e, sc, cfg$n_boot)
        boot_rows[[key]] <- data.frame(index = nm, window = w, method = m,
                                       RMSE_lo = ci[1,"RMSE"], RMSE_hi = ci[2,"RMSE"], MAE_lo = ci[1,"MAE"], MAE_hi = ci[2,"MAE"],
                                       MASE_lo = ci[1,"MASE"], MASE_hi = ci[2,"MASE"], sMAPE_lo = ci[1,"sMAPE"], sMAPE_hi = ci[2,"sMAPE"],
                                       row.names = NULL)
      }
    }
    write.csv(do.call(rbind, boot_rows), "bootstrap_cis.csv", row.names = FALSE)
  }
  
  # visuals: heatmaps, ranking, order frequency
  ggsave("mase_heatmaps.png",
         ggplot(results, aes(factor(window), index, fill = MASE)) + geom_tile() +
           geom_text(aes(label = round(MASE, 2)), colour = "white") + facet_wrap(~method) +
           labs(x = "window", y = NULL) + scale_fill_viridis_c(), width = 10, height = 7)
  rank <- aggregate(MASE ~ method, results, mean)
  ggsave("method_ranking.png",
         ggplot(rank, aes(reorder(method, MASE), MASE)) + geom_col(fill = "slateblue") +
           labs(x = NULL, title = "Mean MASE by method"), width = 6, height = 4)
  ord_tab <- as.data.frame(table(sapply(orders_by_cell, function(o) paste(o$order, collapse = ","))))
  ggsave("order_frequency.png",
         ggplot(ord_tab, aes(Var1, Freq)) + geom_col(fill = "darkcyan") +
           labs(x = "ARIMA order", title = "Selected ARIMA order frequency") +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)), width = 8, height = 4)
  
  # window/accuracy trend (Spearman) + global AIC vs BIC (#12)
  ar <- results[results$method == "ARIMA", ]
  sp <- suppressWarnings(cor.test(ar$window, ar$MASE, method = "spearman"))
  trend_txt <- sprintf("Spearman rho=%.2f (p=%.3f): MASE %s with window length",
                       sp$estimate, sp$p.value,
                       if (sp$p.value < 0.05 && sp$estimate > 0) "increased significantly"
                       else if (sp$estimate > 0) "tended to increase (not significant)" else "did not increase")
  ic_dis <- sum(mapply(function(o, ao) !identical(paste(o$order, collapse=","), ao),
                       orders_by_cell, ar$aic_order[match(names(orders_by_cell),
                                                          paste(ar$index, ar$window))]), na.rm = TRUE)
  aicbic_txt <- sprintf("AIC vs BIC: mean MASE %.3f (BIC) vs %.3f (AIC) across cells; orders differed in %d/%d.",
                        mean(ar$MASE), mean(ar$aic_MASE, na.rm = TRUE), ic_dis, n_combos)
  
  # sensitivity
  sens_txt <- "Sensitivity skipped."
  if (cfg$run_sens) {
    sens <- sapply(cfg$train_fracs, function(fr) {
      r <- run_grid(price_list, labels, fr, cfg)$results
      tapply(r$MASE, r$method, mean) })
    colnames(sens) <- paste0(cfg$train_fracs * 100, "%")
    png("sensitivity_trainfrac.png", 800, 500)
    matplot(cfg$train_fracs, t(sens), type = "b", pch = 1, lty = 1,
            xlab = "training fraction", ylab = "mean MASE", main = "Sensitivity")
    legend("topright", rownames(sens), col = 1:nrow(sens), lty = 1); dev.off()
    sens_txt <- paste(capture.output(print(round(sens, 3))), collapse = "\n")
  }
  
  write_discussion(results, dm_sig, n_combos, trend_txt, aicbic_txt, sens_txt, R)
  cat("\nSaved: wf_results.csv, bootstrap_cis.csv, discussion.txt, and PNG figures.\n")
}

#discussion
write_discussion <- function(results, dm_sig, n_combos, trend_txt, aicbic_txt, sens_txt, R) {
  mm <- sort(tapply(results$MASE, results$method, mean))
  disc <- c(
    "INTERPRETATION AND DISCUSSION", strrep("=", 60),
    sprintf("* Lowest mean MASE: %s (%.3f). Ranking: %s.", names(mm)[1], mm[1],
            paste(sprintf("%s %.3f", names(mm), mm), collapse = ", ")),
    sprintf("* DM vs naive: ARIMA beat naive in %d/%d cells (SES %d, Holt %d).",
            dm_sig["ARIMA"], n_combos, dm_sig["SES"], dm_sig["Holt"]),
    sprintf("* %s", aicbic_txt),
    sprintf("* %s.", trend_txt),
    sprintf("* Representative (%s %dd): Ljung-Box p=%.3f, Jarque-Bera p=%.3f, 95%% PI coverage %.1f%%.",
            R$ri, R$rw, R$lb_p, R$jb_p, 100 * R$coverage),
    sprintf("* Representative bootstrap 95%% CI for MASE: ARIMA [%.3f, %.3f] vs Naive [%.3f, %.3f] (see bootstrap_cis.csv for all cells/metrics).",
            R$boot$ARIMA["2.5%","MASE"], R$boot$ARIMA["97.5%","MASE"],
            R$boot$Naive["2.5%","MASE"], R$boot$Naive["97.5%","MASE"]),
    "",
    "PRACTICAL INTERPRETATION", strrep("=", 60),
    "* Volatility forecasts drive Value-at-Risk / expected shortfall, option pricing (realized vs implied vol), volatility-targeting position sizing, and margin. Matching only a naive persistence forecast still supports risk MEASUREMENT but offers no directional edge.",
    "",
    "REFERENCES", strrep("=", 60),
    "* Mandelbrot, B. (1963). The Variation of Certain Speculative Prices. Journal of Business, 36(4), 394-419.",
    "* Engle, R. F. (1982). Autoregressive Conditional Heteroscedasticity with Estimates of the Variance of United Kingdom Inflation. Econometrica, 50(4), 987-1007.",
    "* Bollerslev, T. (1986). Generalized Autoregressive Conditional Heteroskedasticity. Journal of Econometrics, 31(3), 307-327.",
    "* Diebold, F. X., & Mariano, R. S. (1995). Comparing Predictive Accuracy. Journal of Business & Economic Statistics, 13(3), 253-263.",
    "* Andersen, T. G., Bollerslev, T., Diebold, F. X., & Labys, P. (2003). Modeling and Forecasting Realized Volatility. Econometrica, 71(2), 579-625.",
    "* Hansen, P. R., & Lunde, A. (2005). A Forecast Comparison of Volatility Models: Does Anything Beat a GARCH(1,1)? Journal of Applied Econometrics, 20(7), 873-889.",
    "* Corsi, F. (2009). A Simple Approximate Long-Memory Model of Realized Volatility. Journal of Financial Econometrics, 7(2), 174-196.",
    "* Hyndman, R. J., & Koehler, A. B. (2006). Another Look at Measures of Forecast Accuracy. International Journal of Forecasting, 22(4), 679-688.",
    "",
    "LIMITATIONS AND FUTURE WORK", strrep("=", 60),
    "* ARIMA models the conditional MEAN; volatility clustering is a conditional-VARIANCE effect, captured only indirectly via an RV proxy.",
    "* Residuals are heteroskedastic and heavy-tailed, violating ARIMA's assumptions.",
    "* Overlapping windows inject artificial autocorrelation, flattering the naive benchmark.",
    "* ARIMA is linear/symmetric and cannot capture the leverage effect.",
    "* Future work: GARCH/EGARCH/GJR-GARCH ('rugarch'), HAR-RV ('HARModel'), stochastic volatility, ML models.")
  writeLines(disc, "discussion.txt"); cat("\n", paste(disc, collapse = "\n"), "\n")
}

# ------------------------------- run --------------------------------- #
main()
