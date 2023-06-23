library(ggplot2)
library(astsa)
library(TSA)
library(fGarch)
library(forecast)

#Reads the data file in variable
calls = read.csv("out.csv", header = T)
#Type cast to Time Series
myts <- ts(log(calls[[1]]), frequency=365)
#Auto Arima to find best model
auto.arima(myts, seasonal = T)

#Fitting the best model
m4 <- sarima(log(calls[[1]]),5,1,1,0,1,0,365)
par(mfrow=c(1,1))




#Fitting a TBATS Model
myts %>%
  tbats(seasonal.periods = c(365.25)) -> fit3
#Forecasting the Same
fc3 <- forecast(fit3)
#Plotting the same
autoplot(fc3, fcol = "red")+ theme_bw()
#Print AIC
fit3$AIC
#Plotting the Residuals
plot((as.numeric(fit3$errors)), type = "l", col = "#CC0033", ylab = "Residuals", main = "Residual Plot TBATS")
#ACF of Residuals
acf(as.numeric(fit3$errors), col = "#CC0033", ylab = "ACF", main = "ACF Plot TBATS")
#ACF for 350-400 lags
acf(as.numeric(fit3$errors), lag.max = 400, xlim = c(350,400))
#PACF of the same
pacf(as.numeric(fit3$errors), col = "#CC0033", ylab = "PACF", main = "PACF Plot TBATS")
#PACF of Residuals from 350 to 400
pacf(as.numeric(fit3$errors), lag.max = 400, xlim = c(350,400))
#ACF of Residual Squared
acf((as.numeric(fit3$errors))^2, lag.max = 100)
#PACF of Residual Squared
pacf((as.numeric(fit3$errors))^2, lag.max = 100)
#Ljung Box Test Residuals
Box.test(as.numeric(fit3$errors), type = "Ljung-Box")
#Ljung Box Test Squared Residuals
Box.test((as.numeric(fit3$errors))^2, type = "Ljung-Box")
#ARCH Test
FinTS::ArchTest(fit3$errors, lag = 12)
#Fitting GARCH
a2 = garchFit(~garch(1, 0), data = fit3$errors, trace = F, include.mean = F)
summary(a2)
plot(a2@sigma.t, type = "l")
sresi2 = a2@residuals/a2@sigma.t
Box.test(sresi2^2, lag = 12, type = "Ljung", fitdf = length(coef(a2)))

#x for Seasonal Pattern of TBATS
x = seq(as.Date("2015-01-01"),as.Date("2015-12-31"),by='day')
#x for Trend Pattern of TBATS
x2 = seq(as.Date("2015-01-01"),as.Date("2020-07-01"),by='day')


library(RColorBrewer)  # This is a nice package to give  you great plotting colors. You will need to install it before using for the first time.
temp.color = c(rep(rev(brewer.pal(11, "RdBu")), each = 11), rep(brewer.pal(11, "RdYlBu"),each = 23))  # creating a color vector for the 12 months
plot(x,tbats.components(fit3)[,3][1:365], type = "o", col = temp.color, xlab = "Time", ylab = "Trend", main = "Seasonal Pattern")
plot(x2,tbats.components(fit3)[,2], type = "l", xlab = "Time", ylab = "Trend", main = "Trend Observed", col = "#CC0033")



plots <- list()
for (i in seq(6)) {
  fit <- auto.arima(myts, xreg = fourier(myts, K = i),
                    seasonal = FALSE, lambda = 0)
  plots[[i]] <- autoplot(forecast(fit,
                                  xreg=fourier(myts, K=i, h=730)),
    xlab= (paste("K=",i,"   AICC=",round(fit[["aicc"]],2))),
    ylab =(""))
}
gridExtra::grid.arrange(
  plots[[1]],plots[[2]],plots[[3]],
  plots[[4]],plots[[5]],plots[[6]], nrow=3)

# Extract the seasonal component
fourier_terms = fourier(myts, K = 4)
model <- auto.arima(myts, xreg=fourier_terms, seasonal=FALSE, lambda=0)
autoplot(forecast(model, xreg=fourier(myts, K=4, h=730)), fcol = "red", main = "Forecast with DHR", ylab = "Log PM2.5") + theme_bw()
#seasonal_component <- fourier_terms %*%tail(coef(model), 12)
plot(as.numeric(model$residuals), type = "l", col = "#CC0033", ylab = "Residuals", main = "Residual Plot DHR")
acf(as.numeric(model$residuals), col = "#CC0033", ylab = "ACF", main = "ACF Plot DHR")
acf(as.numeric(model$residuals), lag.max = 400, xlim = c(350,400))
pacf(as.numeric(model$residuals), col = "#CC0033", ylab = "PACF", main = "PACF Plot DHR")
pacf(as.numeric(model$residuals), lag.max = 400, xlim = c(350,400))
acf((as.numeric(model$residuals))^2, lag.max = 100)
pacf((as.numeric(model$residuals))^2, lag.max = 100)
Box.test(as.numeric(model$residuals), type = "Ljung-Box")
Box.test((as.numeric(model$residuals))^2, type = "Ljung-Box")
FinTS::ArchTest(model$residuals, lag = 12)
a1 = garchFit(~garch(1, 0), data = model$residuals, trace = F, include.mean = F)
summary(a1)
plot(a1@sigma.t, type = "l")
sresi = a1@residuals/a1@sigma.t
Box.test(sresi^2, lag = 12, type = "Ljung", fitdf = length(coef(a1)))




