# Demo data for development and testing

# Options
options(encoding = "UTF-8")
options(stringsAsFactors = F)
options(max.print = 10e3)


# Package
require(RJDemetra)


# Reset
rm(list = ls())
gc()


# Data
data("ipi_c_eu")

s1 <- ipi_c_eu[, "FR"]
s2 <- ipi_c_eu[, "IT"]
s3 <- (s1 + s2) / 2

plot(s1)
plot(s2)
plot(s3)

ts_model1 <- tramoseats(s1)
ts_model2 <- tramoseats(s2)
ts_model3 <- tramoseats(s3)

plot(ts_model1, type_chart = "sa-trend")
plot(ts_model2, type_chart = "sa-trend")
plot(ts_model3, type_chart = "sa-trend")

sa1 <- ts_model1$final$series
sa2 <- ts_model2$final$series

tsda <- ts_model3$final$series
tsia <- (sa1 + sa2) / 2

attr(tsia, "dimnames")[[2]] <- attr(tsda, "dimnames")[[2]]

str(tsda)
str(tsia)

plot(tsda)
plot(tsia)

save(tsda, tsia, file = "data/data.Rdata")
