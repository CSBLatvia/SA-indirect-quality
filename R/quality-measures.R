# Step 1: Function for quality measures

# Options
options(encoding = "UTF-8")
options(stringsAsFactors = F)
options(max.print = 10e3)

# Packages
require(RJDemetra)
require(rJava)
require(ggplot2)
require(data.table)


# Reset
rm(list = ls())
gc()


# Extra functions
source("R/jd_seasonality.R")

# Demo data
load("data/data.Rdata")

tsda
tsia

# # Save full copy
# tsda_full <- tsda
# tsia_full <- tsia


# Subset by variables

ts.names <- c("sa", "t", "i")

tsda <- tsda[, ts.names]
tsia <- tsia[, ts.names]

# colnames(tsda) <- paste0("da_", ts.names)
# colnames(tsia) <- paste0("ia_", ts.names)

tscomb <- cbind(tsda, tsia)

start(tscomb)
end(tscomb)
frequency(tscomb)
time(tscomb)

tsdt <- data.table(time = as.numeric(time(tscomb)), tscomb)
tsdt <- melt(tsdt, id.vars = "time")
tsdt[, c("adj", "comp") := tstrsplit(variable, split = ".", fixed = T)]
tsdt[, comp := factor(comp, ts.names, ts.names)]


sapply(tsdt, class)

plot.ts(tscomb)

ggplot(tsdt, aes(x = time, y = value, colour = adj)) +
  geom_line() +
  facet_grid(rows = "comp", scales = "free_y") +
  theme_bw()


# Function to calculate quality measures for indirect and direct seasonal adjusted data
quality.measures <- function(tsda,
                             tsia,
                             years = 3L,
                             digits = 3L) {

  # Testing of input
  if (any(dim(tsda) != dim(tsia))) {
    message("The dimension of the provided data differs.")
    message(paste0("Dimension of direct: ",
                   paste(dim(tsda), collapse = " ")))
    message(paste0("Dimension of indirect: ",
                   paste(dim(tsia), collapse = " ")))
    stop("Cannot proceed.", call. = FALSE)
  }

  # Subselection by time - last n years
  if (is.na(years)) n <- nrow(tsda) else n <- years * frequency(tsda)

  tsda_n <- forecast:::subset.ts(tsda, start = nrow(tsda) - n + 1)
  tsia_n <- forecast:::subset.ts(tsia, start = nrow(tsia) - n + 1)

  # Mean Absolute Percentage Deviation (MAPE)
  MAPE <- function(tsda, tsia, digits) {
    x <- round(100 * colMeans(abs((tsda - tsia) / tsda), na.rm = TRUE),
               digits = digits)
    names(x) <- paste("MAPE", colnames(tsda))
    return(x)
  }

  value.MAPE <- MAPE(tsda[, c("sa", "t")],
                     tsia[, c("sa", "t")], digits = digits)
  value.MAPE.n <- MAPE(tsda_n[, c("sa", "t")],
                       tsia_n[, c("sa", "t")], digits = digits)


  # Max Absolute Percentage Deviation (MaxAPE)
  MaxAPE <- function(tsda, tsia, digits) {
    x <- round(100 * apply(abs((tsda - tsia) / tsda), 2, max,
                           na.rm = TRUE), digits = digits)
    names(x) <- paste("MaxAPE", colnames(tsda))
    return(x)
  }

  value.MaxAPE <- MaxAPE(tsda[, c("sa", "t")],
                         tsia[, c("sa", "t")], digits = digits)
  value.MaxAPE.n <- MaxAPE(tsda_n[, c("sa", "t")],
                           tsia_n[, c("sa", "t")], digits = digits)


  # Inconsistencies
  Incons <- function(x, y, digits = digits) {
    x <- as.vector(x)
    y <- as.vector(y)
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    z <- sign((tail(x, -1) / head(x, -1) - 1) *
              (tail(y, -1) / head(y, -1) - 1))
    return(round(100 * length(z[z == -1]) / length(z[!is.na(z)]),
                 digits = digits))
  }

  value.Incons <- Incons(tsda[, "sa"], tsia[, "sa"], digits = digits)
  value.Incons.n <- Incons(tsda_n[, "sa"], tsia_n[, "sa"], digits = digits)

  # MQ stat

  # SA X-13
  tsda_x13 <- jx13(tsda[, "sa"], spec = "RSA3")
  tsia_x13 <- jx13(tsia[, "sa"], spec = "RSA3")

  MQstat_da <- get_indicators(tsda_x13,
                             c(paste0("mstats.M(", 1:11, ")"),
                               "mstats.Q", "mstats.Q-M2"))
  MQstat_ia <- get_indicators(tsia_x13,
                              c(paste0("mstats.M(", 1:11, ")"),
                                "mstats.Q", "mstats.Q-M2"))

  # Roughness of the components
  R1_da <- round(sum(diff(tsda) ^ 2), digits = digits)
  R1_ia <- round(sum(diff(tsia) ^ 2), digits = digits)

  if (frequency(tsia) == 12) {
    tsia_trend <- get_indicators(jx13(tsia[, "sa"], spec = x13_spec(spec=c("RSA3"), x11.trendma = 13)), "t")[[1]]
    tsda_trend <- get_indicators(jx13(tsda[, "sa"], spec = x13_spec(spec=c("RSA3"), x11.trendma = 13)), "t")[[1]]
  } else if (frequency(tsia) == 4 | frequency(tsia) == 2) {
    tsia_trend <- get_indicators(jx13(tsia[, "sa"], spec = x13_spec(spec=c("RSA3"), x11.trendma =5)), "t")[[1]]
    tsda_trend <- get_indicators(jx13(tsda[, "sa"], spec = x13_spec(spec=c("RSA3"), x11.trendma =5)), "t")[[1]]
  } else {
    message("Trend cannot be computed for the indirectly adjusted series")
    tsia_trend <- NA
    tsda_trend <- NA
  }

  R2_ia <- round(sum((tsia[, "sa"] - tsia_trend) ^ 2), digits = digits)
  R2_da <- round(sum((tsda[, "sa"] - tsda_trend) ^ 2), digits = digits)

  R3_ia <- round(sum((tsia[, "sa"] - tsia[, "t"]) ^ 2), digits = digits)
  R3_da <- round(sum((tsda[, "sa"] - tsda[, "t"]) ^ 2), digits = digits)

  if (frequency(tsia) == 12) {
    MarS_ia <- round(sum((tsia[, "sa"] +
                            lag(tsia[, "sa"], 1) +
                            lag(tsia[, "sa"], 2) +
                            lag(tsia[, "sa"], 3) +
                            lag(tsia[, "sa"], 4) +
                            lag(tsia[, "sa"], 5) +
                            lag(tsia[, "sa"], 6) +
                            lag(tsia[, "sa"], 7) +
                            lag(tsia[, "sa"], 8) +
                            lag(tsia[, "sa"], 9) +
                            lag(tsia[, "sa"], 10) +
                            lag(tsia[, "sa"], 11)) ^ 2), digits = digits)
    MarS_da <- round(sum((tsda[, "sa"] +
                            lag(tsda[, "sa"], 1) +
                            lag(tsda[, "sa"], 2) +
                            lag(tsda[, "sa"], 3) +
                            lag(tsda[, "sa"], 4) +
                            lag(tsda[, "sa"], 5) +
                            lag(tsda[, "sa"], 6) +
                            lag(tsda[, "sa"], 7) +
                            lag(tsda[, "sa"], 8) +
                            lag(tsda[, "sa"], 9) +
                            lag(tsda[, "sa"], 10) +
                            lag(tsda[, "sa"], 11)) ^ 2), digits = digits)
  } else if (frequency(tsia) == 4) {
    MarS_ia <- round(sum((tsia[, "sa"] +
                            lag(tsia[, "sa"], 1) +
                            lag(tsia[, "sa"], 2) +
                            lag(tsia[, "sa"], 3)) ^ 2), digits = digits)
    MarS_da <- round(sum((tsda[, "sa"] +
                            lag(tsda[, "sa"], 1) +
                            lag(tsda[, "sa"], 2) +
                            lag(tsda[, "sa"], 3)) ^ 2), digits = digits)
  } else {
    MarS_ia <- NA
    MarS_da <- NA
  }

  Mar1_ia <- round(sum(diff(tsia[, "sa"]) ^ 2), digits = digits)
  Mar1_da <- round(sum(diff(tsda[, "sa"]) ^ 2), digits = digits)

  Mar2_ia <- round(sum(diff(tsia[, "sa"], differences = 2) ^ 2),
                   digits = digits)
  Mar2_da <- round(sum(diff(tsda[, "sa"], differences = 2) ^ 2),
                   digits = digits)

  # Idempotency
  sa_test_da <- sapply(c("QS", "FRIEDMAN", "KRUSKALWALLIS", "PERIODOGRAM"),
         jd_seasonality, s = tsda[, "sa"], simplify = F)
  sa_test_ia <- sapply(c("QS", "FRIEDMAN", "KRUSKALWALLIS", "PERIODOGRAM"),
         jd_seasonality, s = tsia[, "sa"], simplify = F)


  # Impact on the business cycle estimation
  bkf_da <- mFilter::bkfilter(tsda[, "sa"])
  bkf_ia <- mFilter::bkfilter(tsia[, "sa"])

  bkf <- cbind(da_cycle = bkf_da$cycle, ia_cycle = bkf_ia$cycle)
  plot.ts(x = bkf, main = "Impact on the business cycle estimation")

  # Incons(bkf_da$cycle, bkf_ia$cycle, digits = digits)


  # Characteristics of the irregular component

  # SA TS and X13 on i

  # tsda_i_ts <- tramoseats(tsda[, "i"], spec = "RSAfull")
  # tsia_i_ts <- tramoseats(tsia[, "i"], spec = "RSAfull")
  # tsda_i_x13 <- x13(tsda[, "i"], spec = "RSA5c")
  # tsia_i_x13 <- x13(tsia[, "i"], spec = "RSA5c")

  # tsda_y_ts <- tramoseats(tsda_full[, "y"], spec = "RSAfull")
  # tsia_y_ts <- tramoseats(tsia_full[, "y"], spec = "RSAfull")
  # tsda_y_x13 <- x13(tsda_full[, "y"], spec = "RSA5c")
  # tsia_y_x13 <- x13(tsia_full[, "y"], spec = "RSA5c")

  jtsda_i_ts <- jtramoseats(tsda[, "i"], spec = "RSAfull")
  jtsia_i_ts <- jtramoseats(tsia[, "i"], spec = "RSAfull")
  jtsda_i_x13 <- jx13(tsda[, "i"], spec = "RSA5c")
  jtsia_i_x13 <- jx13(tsia[, "i"], spec = "RSA5c")

  # jtsda_y_ts <- jtramoseats(tsda_full[, "y"], spec = "RSAfull")
  # jtsia_y_ts <- jtramoseats(tsia_full[, "y"], spec = "RSAfull")
  # jtsda_y_x13 <- jx13(tsda_full[, "y"], spec = "RSA5c")
  # jtsia_y_x13 <- jx13(tsia_full[, "y"], spec = "RSA5c")

  # get_dictionary(jtsda_i_ts)
  # grep("preprocessing", get_dictionary(jtsda_i_ts), value = T)
  # grep("preprocessing.model", get_dictionary(jtsda_i_ts), value = T)
  # grep("preprocessing.model.nout", get_dictionary(jtsda_i_ts), value = T)
  # grep("qs", get_dictionary(jtsda_i_ts), value = T)
  # grep("combined", get_dictionary(jtsda_i_ts), value = T)
  #
  # lapply(grep("combined", get_dictionary(jtsda_i_ts), value = T),
  #        get_indicators, x = jtsda_i_ts)
  # lapply(grep("summary", get_dictionary(jtsda_i_ts), value = T),
  #        get_indicators, x = jtsda_i_ts)
  #
  # lapply(grep("\\.ar", get_dictionary(tsda_i_ts), value = T),
  #        get_indicators, x = tsda_i_ts)
  #
  # el(get_indicators(jtsda_i_ts, "preprocessing.arima.p"))
  # el(get_indicators(jtsda_i_ts, "preprocessing.arima.d"))
  # el(get_indicators(jtsda_i_ts, "preprocessing.arima.q"))
  # el(get_indicators(jtsda_i_ts, "preprocessing.arima.bp"))
  # el(get_indicators(jtsda_i_ts, "preprocessing.arima.bd"))
  # el(get_indicators(jtsda_i_ts, "preprocessing.arima.bq"))
  #
  # get_indicators(jtsda_i_ts, "preprocessing.model.lp")
  # get_indicators(jtsda_i_ts, "preprocessing.model.easter")
  # get_indicators(jtsda_y_ts, "preprocessing.model.td(1)")
  # get_indicators(jtsda_y_ts, "preprocessing.model.td(2)")
  # get_indicators(jtsda_i_ts, "preprocessing.model.nout")
  # get_indicators(jtsda_i_ts, "preprocessing.model.noutao")
  # get_indicators(jtsda_i_ts, "preprocessing.model.noutls")
  # get_indicators(jtsda_i_ts, "preprocessing.model.nouttc")
  # get_indicators(jtsda_y_ts, "diagnostics.qs")[[1]][2]
  #
  # get_indicators(jtsda_y_ts, paste0("preprocessing.model.td(",1:7,")"))
  #
  # length(unname(unlist(sapply(paste0("preprocessing.model.td(", 1:7, ")"), get_indicators, x = jtsda_y_ts))))
  # length(unlist(sapply(paste0("preprocessing.model.td(", 1:7, ")"), get_indicators, x = jtsda_y_ts)))
  #
  # tsda_i_ts$regarima$arma
  # tsda_i_ts$regarima$model$spec_rslt$`Leap year`
  # tsda_i_ts$regarima$model$spec_rslt$Easter
  # tsda_i_ts$regarima$model$spec_rslt$`Trading days`
  # tsda_i_ts$regarima$model$spec_rslt$Outliers
  #
  # tsda_i_ts$diagnostics$combined_test$combined_seasonality_test
  #
  # tsda_y_ts$regarima$model$spec_rslt$`Trading days`
  # tsia_y_ts$regarima$model$spec_rslt$`Trading days`
  # tsda_y_ts$regarima$model$spec_rslt$`Trading days`
  # tsia_y_ts$regarima$model$spec_rslt$`Trading days`

  # Nice model print out
  print_model <- function(x) {
    foo <- function(x) {
      paste0("(", paste(x, collapse = ","), ")")
    }
    paste0(foo(x[1:3]), foo(x[4:6]))
  }

  # get_indicators(jtsda_i_ts, "preprocessing.model.lp")
  # get_indicators(jtsia_i_ts, "preprocessing.model.lp")
  # get_indicators(jtsda_i_x13, "preprocessing.model.lp")
  # get_indicators(jtsia_i_x13, "preprocessing.model.lp")

  # get_spec <- function(x) {
  #   data.frame(model = print_model(x$regarima$arma),
  #              leapyear = x$regarima$model$spec_rslt$`Leap year`,
  #              easter = x$regarima$model$spec_rslt$Easter,
  #              tradingdays = x$regarima$model$spec_rslt$`Trading days` -
  #                x$regarima$model$spec_rslt$`Leap year`,
  #              outliers = x$regarima$model$spec_rslt$Outliers,
  #              comb_sa_test = x$diagnostics$combined_test$combined_seasonality_test
  #   )
  # }
  #
  # get_spec(tsda_i_ts)
  #
  # x <- do.call(rbind, lapply(list(tsia_i_ts, tsia_i_x13,
  #                                 tsda_i_ts, tsda_i_x13), get_spec))
  # data.frame(series = c("Seats Indirect", "X-13 Indirect",
  #                       "Seats Direct", "X-13 Direct"), x)

  # Function for information extraction from jsa
  get_jspec <- function(x) {
    leapyear <- el(get_indicators(x, "preprocessing.model.lp"))
    if (is.null(leapyear)) leapyear <- FALSE
      else leapyear <- grepl("Leap year", leapyear)

    easter <- el(get_indicators(x, "preprocessing.model.easter"))
    if (is.null(easter)) easter <- FALSE
      else easter <- grepl("Easter", easter)

    data.frame(model = print_model(c(el(get_indicators(x, "preprocessing.arima.p")),
                                     el(get_indicators(x, "preprocessing.arima.d")),
                                     el(get_indicators(x, "preprocessing.arima.q")),
                                     el(get_indicators(x, "preprocessing.arima.bp")),
                                     el(get_indicators(x, "preprocessing.arima.bd")),
                                     el(get_indicators(x, "preprocessing.arima.bq")))),
               leapyear = leapyear,
               easter = easter,
               tradingdays = length(unlist(sapply(paste0("preprocessing.model.td(", 1:7, ")"),
                                                  get_indicators, x = x))),
               ao = el(get_indicators(x, "preprocessing.model.noutao")),
               ls = el(get_indicators(x, "preprocessing.model.noutls")),
               tc = el(get_indicators(x, "preprocessing.model.nouttc")),
               qs.pvalue = el(get_indicators(x, "diagnostics.seas-lin-qs"))[2],
               comb_sa_test = el(get_indicators(x, "diagnostics.combined.all.summary")))
  }

  # get_jspec(jtsda_i_ts)

  y <- do.call(rbind, lapply(list(jtsia_i_ts, jtsia_i_x13,
                                  jtsda_i_ts, jtsda_i_x13), get_jspec))
  i_stat <- data.frame(series = c("Seats Indirect", "X-13 Indirect",
                                  "Seats Direct", "X-13 Direct"), y)

  # x <- do.call(rbind, lapply(list(tsia_y_ts, tsia_y_x13,
  #                                 tsda_y_ts, tsda_y_x13), get_spec))
  # data.frame(series = c("Seats Indirect", "X-13 Indirect",
  #                       "Seats Direct", "X-13 Direct"), x)
  #
  # y <- do.call(rbind, lapply(list(jtsia_y_ts, jtsia_y_x13,
  #                                 jtsda_y_ts, jtsda_y_x13), get_jspec))
  # data.frame(series = c("Seats Indirect", "X-13 Indirect",
  #                       "Seats Direct", "X-13 Direct"), y)

  return(list(MAPE = list(full = value.MAPE, last_n = value.MAPE.n),
              MaxAPE = list(full = value.MaxAPE, last_n = value.MaxAPE.n),
              Incons = list(full = value.Incons, last_n = value.Incons.n),
              MQstat = list(da = unlist(MQstat_da), ia = unlist(MQstat_ia)),
              R1 = list(da = R1_da, ia = R1_ia),
              R2 = list(da = R2_da, ia = R2_ia),
              R3 = list(da = R3_da, ia = R3_ia),
              MarS = list(da = MarS_da, ia = MarS_ia),
              Mar1 = list(da = Mar1_da, ia = Mar1_ia),
              Mar2 = list(da = Mar2_da, ia = Mar2_ia),
              SA_test = list(da = sa_test_da, ia = sa_test_ia),
              Ireg_Stat = i_stat))
}

# quality.measures(tsda, tsda)
quality.measures(tsda, tsia, digits = 6)
# quality.measures(tsda, tsia, years = 1.5, digits = 10)
