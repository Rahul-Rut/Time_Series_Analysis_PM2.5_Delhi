library(ggplot2)
library(astsa)
library(TSA)
library(fGarch)
library(forecast)
library(feather)
library(data.table)
library(mgcv)
library(car)
library(ggplot2)
library(grid)
library(zoo)
library(lubridate)
library(mboost)

#Reads the data file in variable
calls = read.csv("out.csv", header = T)
#Type cast to Time Series
myts <- ts(log(calls[[1]]), frequency=365)
out = tsoutliers(exp(myts)) 
myts = replace(exp(myts), out$index, out$replacements)


new_data = read.csv("out2.csv", header = T)
new_data$X0 = ts(log(new_data$X0), frequency=365)
new_data$DT = as.Date(new_data$DT, format = "%Y-%m-%d")
new_data[,'quarter'] = quarter(new_data$DT)
new_data[,'week'] = week(new_data$DT)
new_data[,'year'] = year(new_data$DT)
new_data[,'month'] = month(new_data$DT)
new_data[,'cont'] = seq(length(myts))/length(myts) + 1
new_data[,'day_num'] = yday(new_data$DT)
pred_set = new_data[(length(myts)-328):(length(myts)),c('year','day_num')]
ctrl<-list(niterEM = 0, optimMethod = "L-BFGS-B", maxIter = 100, msMaxIter = 100, msVerbose=TRUE)

gam_model_2 = gam(X0 ~ s(year, k=5, bs="cr")+ s(day_num,k=16, bs="cp") ,
                 data = new_data[1:(length(myts)-329),],
                 family = gaussian,
                 method = 'REML')

preds_gam_2 = predict(gam_model_2, newdata = pred_set, type = "response")
plot(preds_gam_2, type="l")
ar_errors = auto.arima(gam_model_2$residuals)
resids_gam = new_data[1:(length(myts)-329), 'X0'] - ar_errors$fitted - gam_model_2$fitted.values
fitted_gam = ar_errors$fitted + gam_model_2$fitted.values
ar_for = forecast(ar_errors, h = length(pred_set[[1]]))
forecast_gam = (as.numeric(preds_gam_2) + ar_for$mean)

gamboost(X0 ~ bbs(day_num, knots = 20), data = new_data)

## Outlier Testing!!!!!!


calls = read.csv("out.csv", header = T)
#Type cast to Time Series
myts <- ts( (calls[[1]]), frequency=365)
out = tsoutliers((myts)) 
out$index
myts2 = replace((myts), out$index, out$replacements)
plot(myts2)
plot(myts, col = "#CC0033", ylab="AQI", main = "Outliers")
points((out$index/365)+1,myts[out$index],col="blue")

fourier_terms = fourier(myts, K = 4)
model <- auto.arima(myts, xreg=fourier_terms, seasonal=FALSE, lambda=0)
seas = ((fourier_terms[,1]*(model$coef)[4]) + (fourier_terms[,2]*(model$coef)[5]) + (fourier_terms[,3]*(model$coef)[6]) + (fourier_terms[,4]*(model$coef)[7]) + (fourier_terms[,5]*(model$coef)[8]) + (fourier_terms[,6]*(model$coef)[9]) + (fourier_terms[,7]*(model$coef)[10]) + (fourier_terms[,8]*(model$coef)[11]))
plot(seas)


fourier_terms2 = fourier(myts2, K = 4)
model2 <- auto.arima(myts2, xreg=fourier_terms2, seasonal=FALSE, lambda=0)
seas2 = ((fourier_terms2[,1]*(model2$coef)[4]) + (fourier_terms2[,2]*(model2$coef)[5]) + (fourier_terms2[,3]*(model2$coef)[6]) + (fourier_terms2[,4]*(model2$coef)[7]) + (fourier_terms2[,5]*(model2$coef)[8]) + (fourier_terms2[,6]*(model2$coef)[9]) + (fourier_terms2[,7]*(model2$coef)[10]) + (fourier_terms2[,8]*(model2$coef)[11]) )
plot(seas2)
plot(log(myts2)-seas2)


#Outlier

out = tsoutliers(log(myts))
plot(log(myts), col = "#CC0033", ylab="AQI", main = "Outliers")
points((out$index/365)+1,log(myts)[out$index],col="blue")
legend(1,6.5,legend=c("Outliers"), col=c('blue'), pch=c(1), cex=0.75)
te = log(myts)
te[out$index] = out$replacements
plot(te, col = "#CC0033", ylab="AQI", main = "Outliers Replaced")
points((out$index/365)+1,te[out$index],col="blue")
legend(1,6.5,legend=c("Replacements"), col=c('blue'), pch=c(1), cex=0.75)


fourier7 = fourier(log(myts), K = 5)
model7 <- auto.arima(log(myts), xreg=fourier7, seasonal=FALSE, lambda=0)
plot(model7$fitted)
seas_comp = ((fourier7[,1]*(model7$coef)[4]) + (fourier7[,2]*(model7$coef)[5]) + (fourier7[,3]*(model7$coef)[6]) + (fourier7[,4]*(model7$coef)[7]) + (fourier7[,5]*(model7$coef)[8]) + (fourier7[,6]*(model7$coef)[9]) + (fourier7[,7]*(model7$coef)[10]) + (fourier7[,8]*(model7$coef)[11])+ (fourier7[,9]*(model7$coef)[12]) + (fourier7[,10]*(model7$coef)[13]))
plot(seas_comp, type="l")
autoplot((forecast(model7, xreg=fourier(log(myts), K=5, h=365))))
trend_comp = log(myts) - (4*seas_comp)
plot(trend_comp)
detrend = supsmu(tt,trend_comp)
plot(detrend)
detrend = trend_comp - detrend$y
plot(detrend)
qs2 = quantile(detrend)
iqr_resids2 = qs2[4] - qs2[2]
outliers_ind2 = which((detrend < (qs2[2] - (2*iqr_resids2))) | (detrend > (qs2[4]+(2*iqr_resids2)) ) ,arr.ind = TRUE)
plot(log(myts))
points((outliers_ind2)/365+1, log(myts)[outliers_ind2], col = 'red')
plot(resids)



plot(exp(new_ts))
plot(resids)

log(myts) %>%
  tbats(seasonal.periods = c(365.25)) -> temp_fit
resids = log(myts) - temp_fit$fitted.values
qs = quantile(resids)
iqr_resids = qs[4] - qs[2]
outliers_ind = which((resids < (qs[2] - (3*iqr_resids))) | (resids > (qs[4]+(3*iqr_resids)) ) ,arr.ind = TRUE)
plot(log(myts))
points((outliers_ind)/365+1, log(myts)[outliers_ind], col = 'red')
plot(resids)
points((outliers_ind)/365+1, resids[outliers_ind], col = 'red')
resids[outliers_ind] <- NA
missng <- is.na(resids)
tt = 1:length(resids)
idx <- tt[!missng]
new_resids = (ts(approx(idx, resids[idx], tt, rule = 2)$y, frequency=365))
plot(new_resids)
new_ts = (new_resids+temp_fit$fitted.values)
new_ts %>%
  tbats(seasonal.periods = c(365.25)) -> fit4
plot(fit4)
plot(tbats.components(fit4)[,3])
pacf(new_ts - fit4$fitted.values, lag.max = 400)


fouriers = fourier(new_ts, K = 5)
model4 <- auto.arima(new_ts, xreg=fouriers, seasonal=FALSE, lambda=0)
plot(model4$fitted)
seas_comp = ((fouriers[,1]*(model4$coef)[4]) + (fouriers[,2]*(model4$coef)[5]) + (fouriers[,3]*(model4$coef)[6]) + (fouriers[,4]*(model4$coef)[7]) + (fouriers[,5]*(model4$coef)[8]) + (fouriers[,6]*(model4$coef)[9]) + (fouriers[,7]*(model4$coef)[10]) + (fouriers[,8]*(model4$coef)[11]))
plot(seas_comp, type="l")
autoplot(forecast(model4, xreg=fourier(new_ts, K=5, h=365)))
plot(exp(new_ts))
plot(resids)



#IDK
# stl_ts = (mstl(log(myts)))
# plot(seasadj(stl_ts))
# xx<-(seasadj(stl_ts))
# plot(supsmu(tt,xx))
# trend = (supsmu(tt,xx))
# plot(xx-trend$y)
# acf(xx-trend$y)
# plot(xx-trend$y)
# stl_ts = (mstl(log(myts)))
# plot(stl_ts)
# remainder(stl_ts)



