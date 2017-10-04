library("arfima")


# ======== Controls ======== #
for1 <- arfima.sim(1000)       # 1st forecast series, here generated
for2 <- arfima.sim(1000)       # 2nd forecast series, here generated
baseline <- arfima.sim(1000)   # baseline series


# ======== Initialization ======== #
suma <- NULL
sumam <- NULL
i <- sqrt(as.complex(-1))
len = length(for1)   # length of the series
m <- round(len/(log(len)^2),0)


# ======== Calculate MAC and DM ======== #
lossfunction <- (baseline-for1)^2-(baseline-for2)^2   # can be changed to a different one

arfimafit <- arfima(lossfunction, order = c(0, 0, 0), quiet = TRUE)
d <-  get("dfrac",(get("modes",arfimafit))[[1]])

pd <- (2*gamma(1-2*d)*sin(pi*d))/(d*(1+2*d))

for (j in 1:(m)){
lambda <- 2*pi*j/len

for (t in 1:(len)){
  suma[t] <- exp(i*t*lambda)*lossfunction[t]}

pgram <- (abs(sum(suma))^2)/(2*pi*len)
sumam[j] <- lambda^(2*d)*pgram
}
bm <- sum(sumam)/m
DM <- (sqrt(len)*sum((lossfunction))/len)/sd(lossfunction)
MAC <- (len^(0.5-d))*(sum(lossfunction)/len)/sqrt(bm*pd)

print(paste0('Durbin-Mariano statistic: ',round(DM,5)))
print(paste0('MAC estimator           : ',round(MAC,5)))
