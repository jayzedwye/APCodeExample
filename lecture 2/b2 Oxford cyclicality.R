rm(list = ls())
library(data.table)
library(plyr)
library(stargazer)
library(lfe)
library(xtable)
library(lubridate)
library(gmm)
library(MASS)
library(latex2exp)
library(zoo)

##################################################################
# Main
##################################################################
only.developed <- 1

load(file = "main.RData")
start.year <- 1980
m <- m[year(date) >= start.year]
m <- m[developed >= only.developed]


setorder(m, country, date)
m[, dq := c(diff(q), NA), by = .(country)]
m[, r := dq + fp]

countries <- sort(unique(m$country))
n.country <- length(countries)
m[, year := year(date)]


setorder(m, country, date)
lag.d <- 4

m[, surplus := surplus.oxford / debt.oxford]
m[, dsurplus := surplus - c(rep(NA, lag.d), head(surplus, -lag.d)), by = .(country)]

m <- m[yrqtr >= 198001 & yrqtr <= 201803]

count <- m[, .(a = sum(!is.na(dsurplus))), by = .(country)]
countries <- unique(sort(count[a > 10]$country))
m <- m[country %in% countries]
n.country <- length(countries)

##################################################################
# Construct portfolio
##################################################################

dates <- sort(unique(m$yrqtr))
n.date <- length(dates)

hml <- NULL
for (i in 1:n.date)
{
  tmp <- m[yrqtr == dates[i] & !is.na(fp)]
  hml <- rbind(hml, data.table(
    yrqtr = dates[i],
    hml = mean(tmp[fp >= median(fp) & !is.na(r)]$r) - mean(tmp[fp < median(fp) & !is.na(r)]$r)
  ))
}
m <- merge(m, hml, by = c("yrqtr"))

fac <- NULL
for (i in 1:n.country)
{
  tmp <- m[country != countries[i], .(fsurplus = mean(dsurplus, na.rm = T)), by = .(yrqtr)]
  tmp[, country := countries[i]]
  fac <- rbind(fac, tmp)
}
m <- merge(m, fac, by = c("country", "yrqtr"))


##################################################################
# Carry Beta
##################################################################

report <- NULL

# GMM
m[, yrqtr := as.yearqtr(date)]
sq <- m[!is.na(dsurplus) & !is.na(r)]
sq <- sq[, count := length(country), by = .(yrqtr)]
sq <- sq[count == max(sq$count)]


##################################################################
# Standard AP GMM
##################################################################
n.GMM.date <- length(unique(sq$date))
GMM.countries <- sort(unique(sq$country))
n.GMM.country <- length(GMM.countries) - 1

X <- matrix(0, nrow = n.GMM.date, ncol = 1 + n.GMM.country)
X[, 1] <- sq[, .(hml = mean(hml)), by = .(date)]$hml
for (i in 1:n.GMM.country)
{
  X[, 1 + i] <- sq[country == GMM.countries[i]]$r
}
X <- X * 100
colMeans(X)


g1 <- function(tet, x) {
  iii <- 1:n.GMM.country
  beta <- tet[iii]
  lambda <- tet[n.GMM.country + 1]
  lambda0 <- tet[n.GMM.country + 2]

  f <- matrix(0, nrow = n.GMM.date, ncol = 2 * n.GMM.country)
  for (t in 1:n.GMM.date)
  {
    f[t, iii] <- ((x[t, 1 + iii] - mean(x[, 1 + iii])) - beta * (x[t, 1] - mean(x[, 1]))) *
       (x[t, 1] - mean(x[, 1])) # return ~ carry
    f[t, n.GMM.country + iii] <- x[t, 1 + iii] - (lambda * beta + lambda0) # return - lambda*b
  }
  return(f)
}

init.tet <- (1:(1 * n.GMM.country)) * 0
for (i in 1:n.GMM.country)
{
  init.tet[i] <- summary(lm(X[, 1 + i] ~ X[, 1]))$coefficients[[2]]
}
init.tet[1 * n.GMM.country + 1] <- 0
init.tet[1 * n.GMM.country + 2] <- 0

res1 <- gmm(g = g1, x = X, t0 = init.tet, wmatrix = "ident", vcov = "iid", optfct = "nlminb")
summary(res1)


##################################################################
# Jiang (2022) GMM
##################################################################
n.GMM.date <- length(unique(sq$date))
GMM.countries <- sort(unique(sq$country))
n.GMM.country <- length(GMM.countries)

X <- matrix(0, nrow = n.GMM.date, ncol = 1 + n.GMM.country * 3)
X[, 1] <- sq[, .(hml = mean(hml)), by = .(date)]$hml
for (i in 1:n.GMM.country)
{
  X[, 1 + i] <- sq[country == GMM.countries[i]]$fsurplus
  X[, 1 + n.GMM.country + i] <- sq[country == GMM.countries[i]]$dsurplus
  X[, 1 + 2 * n.GMM.country + i] <- sq[country == GMM.countries[i]]$r
}
for (i in 1:(1 + n.GMM.country * 3))
{
  X[, i] <- X[, i] - mean(X[, i])
}


g1 <- function(tet, x) {
  iii <- 1:n.GMM.country
  b <- tet[iii]
  zeta <- tet[n.GMM.country + 1]
  zeta0 <- tet[n.GMM.country + 2]

  f <- matrix(0, nrow = n.GMM.date, ncol = 2 * n.GMM.country)
  for (t in 1:n.GMM.date)
  {
    f[t, iii] <- (x[t, 1 + n.GMM.country + iii] - b * x[t, 1 + iii]) * x[t, 1 + iii] # dsurplus ~ fsurplus
    f[t, n.GMM.country + iii] <- (x[t, 1 + 2 * n.GMM.country + iii] - (zeta * b + zeta0) * x[t, 1]) * x[t, 1] # return ~ carry
  }
  return(f)
}

init.tet <- (1:(1 * n.GMM.country + 2)) * 0
for (i in 1:n.GMM.country)
{
  init.tet[i] <- summary(lm(X[, 1 + n.GMM.country + i] ~ X[, 1 + i]))$coefficients[[2]]
}
init.tet[1 * n.GMM.country + 1] <- 0
init.tet[1 * n.GMM.country + 2] <- 0

res1 <- gmm(g = g1, x = X, t0 = init.tet, wmatrix = "ident", vcov = "iid", optfct = "nlminb")
summary(res1)


output <- rbind(report, data.table(
  sample = "Oxford Economics",
  n.t = n.GMM.date,
  slope = summary(res1)$coef[n.GMM.country + 1, 1],
  se = summary(res1)$coef[n.GMM.country + 1, 2]
))

output[, se := sprintf("(%.2f)", se)]
print(xtable(output), include.rownames = F, math.style.negative = T)
