rm(list=ls())
library(astsa)
library(ggplot2)
library(forecast)
library(animation)
library(itsmr) #przeciąża forecast!
library(tsoutliers)
library(lmtest)


#wczytujemy dane
szereg <- ts(gnp, start=start(gnp), end=end(gnp), frequency=4)
length(szereg)
ggtsdisplay(szereg)
anyNA(szereg)
#(a)
help('window')
train <- window(szereg, end=c(1996, 4))
test <- window(szereg, start=c(1997, 1))
length(train)
length(test)
#(b)
#obserwujemy trend
autoplot(cbind(test, train), lwd=1, main='Podział szeregu gnp') +
  labs(x='Czas', y='Szereg gnp')

#przeksztalcamy do postaci stacjonarnej
#czy potrzebna jest transformacja Boxa-Coxa?

#przeanalizujmy wybór lambdy używając animacji
lambda <- seq(from=0, to=2, by=0.2)
lambda <- lambda[lambda != 1] #usuwamy lambda=1, bo to wyjściowy szereg
oopt = ani.options(interval=0.1)

for (p in 1:length(lambda)) {
  transformacja <- BoxCox(train, lambda=lambda[p])
  p <- forecast::autoplot(transformacja, lwd=1) + 
    ggtitle(paste0("lambda = ", lambda[p])) + 
    theme(plot.title=element_text(size=22), 
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +  
    labs(x="Czas", y="Szereg gnp") 
  plot(p)
  
  ani.pause()
}
ani.options(oopt)

#wybieramy lambda = 0
train_log <- BoxCox(train, lambda=0)

#różnicowanie jednokrotne z jest już satysfakcjonujące
ggtsdisplay(diff(train_log, differences=1), main='Różnicowanie jednokrotne')
ggtsdisplay(diff(train_log, differences = 2),
            main='Różnicowanie 2-krotne')
train_log_diff <- diff(train_log, differences = 1)

#sprawdźmy jeszcze wybór parametru lambda oraz krotności różnicowania przy użyciu wbudowanych
#funkcji
BoxCox.lambda(train)
ndiffs(train_log)

#(c)
#na bazie estymatorów funkcji ACF i PACF możemy zidentyfikować modele:
# - MA(2) lub MA(5)
# - AR(1) lub AR(12)

help('auto.arima')
#szukamy modelu ARMA - kryterium AIC - ARIMA(2,2) z dryfem, czyli ARMA(2,2)
#po 1-krotnym zróżnicowaniu zlogarytmowanego szeregu
fit_arma_aic <- auto.arima(
    y=train,
    trace=TRUE,
    stepwise=FALSE,
    lambda=0,
    ic='aic',
    method='ML'
    )

#szukamy modelu ARMA - kryterium BIC - ARIMA(1,1,0) z dryfem, czyli ARMA(1,0)
#po 1-krotnym zróżnicowaniu zlogarytmowanego szeregu
fit_arma_bic <- auto.arima(
  y=train,
  trace=TRUE,
  stepwise=FALSE,
  lambda=0,
  ic='bic',
  method='ML'
)

#szukamy modelu ARMA - kryterium AICc - ARIMA(2,1,2) z dryfem, czyli ARMA(2,2) po
#1-krotnym zróżnicowaniu i użyciu transformacji logarytmicznej
fit_arma_aicc <- auto.arima(
  y=train,
  trace=TRUE,
  stepwise=FALSE,
  lambda=0,
  ic='aicc',
  method='ML'
)

#z ciekawości sprawdzamy jeszcze dopasowanie modelu AR i MA do train_log_diff, używając do tego
#kryteria AIC, BIC, AICc
#najpierw model AR - AIC wybralo AR(3) a BIC i AICc AR(1)
#po 1-krotnym zróżnicowaniu i użyciu transformacji logarytmicznej
fit_ar_aic <- auto.arima(train, trace=T,lambda=0, stepwise=F, max.q=0, ic='aic')
fit_ar_aicc <- auto.arima(train, trace=T, lambda=0, stepwise=F, max.q=0, ic='aicc')
fit_ar_bic <- auto.arima(train, trace=T, lambda=0, stepwise=F, max.q=0, ic='bic')

#teraz model MA - kryterium AIC wybralo model ARIMA(0,1,2)(0,0,1)[4] z dryfem, AICc MA(2),
#BIC też MA(2)
fit_ma_aic <- auto.arima(train, trace=T,
                         lambda=0, stepwise=F, max.p=0, ic='aic')
fit_ma_aicc <- auto.arima(train, trace=T,
                           lambda=0, stepwise=F, max.p=0, ic='aicc')
fit_ma_bic <- auto.arima(train, trace=T,
                          lambda=0, stepwise=F, max.p=0, ic='bic')

#w dalszej analizie uwzglednimy modele:
# -Arma(2,2) 
# -AR(1) 
# -MA(2) 
fit_arma <- Arima(train,
                  order=c(2,1,2),
                  seasonal=c(0,0,0),
                  lambda=0,
                  include.drift=TRUE
                  )

fit_ar <- Arima(train,
                order=c(1,1,0),
                seasonal=c(0,0,0),
                lambda=0,
                include.drift = TRUE
                )
fit_ma <- Arima(train,
                order=c(0,1,2),
                seasonal=c(0,0,0),
                lambda=0,
                include.drift = TRUE
)

#(d)
#wyznaczamy reszty
res_ma <- fit_ma$residuals
res_ar <- fit_ar$residuals
res_arma <- fit_arma$residuals

#najlepiej prezentuje sie model ARMA(2,2)
ggtsdisplay(res_ma, main="Reszty dla modelu Ma(2)")
ggtsdisplay(res_ar, main="Reszty dla modelu AR(1)")
ggtsdisplay(res_arma, main='Reszty dla modelu ARMA(2,2)')
forecast::checkresiduals(fit_ar_aicc)
forecast::checkresiduals(fit_ma_aicc)
forecast::checkresiduals(fit_arma_aicc)


#graficzny test bialoszumowości
graph_test <- function(ts, conf.int=1.96/sqrt(NROW(ts)),lag_max=floor(NROW(ts)/4),plot=TRUE){
  if(plot==TRUE){
    ts_name <- deparse(substitute(ts))
    plot<- ggAcf(ts,ci=NA,col='white',lag.max = lag_max, main=paste('Test graficzny dla szeregu: ', ts_name))+
      geom_hline(yintercept = conf.int, linetype='dashed',color='lightblue')+
      geom_hline(yintercept = -conf.int,linetype='dashed',color='lightblue')+
      theme_dark()
    cond_1 <- sum(abs(plot$data$Freq)<=conf.int)/lag_max #warunek pierwszy testu graficznego białoszumowości (95% znajdujących się wewnątrz przedziału)
    cond_2 <- any(abs(plot$data$Freq)>=conf.int+0.05*conf.int) #warunek drugi testu graficznego białoszumowości (żaden *ISTOTNIE* nie przekracza przedziału )
    lista <- list(cond_1,cond_2,plot)
    return(lista)}
  else{ #dodatkowe, w celu przyspieszenia obliczeń dla symulacji - wykres niepotrzebny
    acf <-  ggAcf(ts,plot=FALSE,lag.max=lag_max)$acf[-1]
    cond_1 <- sum(abs(acf)<=conf.int)/lag_max
    cond_2 <- any(abs(acf)>=conf.int+0.05*conf.int)
    lista <- list(cond_1, cond_2,acf)
    return(lista)}}

#sprawdzamy wyniki testu graficznego - najlepiej wypada model ARMA(2,2) wraz z MA(2)
graph_test(res_ar)[1]
graph_test(res_ma)[1]
graph_test(res_arma)[1]

#rysujemy p-wartości testu Ljungi-Boxa dla modelów - najlepiej wypada model ARMA(2,2)
#ggfortify::ggtsdiag(fit_arma)
#ggfortify::ggtsdiag(fit_ar)
#ggfortify::ggtsdiag(fit_ma)
tsdiag(fit_arma)
tsdiag(fit_ar)
tsdiag(fit_ma)

#testy losowości - wszystkie testy twierdzą przyjmują H0: residuals are iid white noise
itsmr::test(res_ar)
itsmr::test(res_ma)
itsmr::test(res_arma)

#testy normalności
JarqueBera.test(res_ar)
JarqueBera.test(res_ma)
JarqueBera.test(res_arma)

shapiro.test(res_ar)
shapiro.test(res_ma)  
shapiro.test(res_arma)

#obserwujemy odstepstwa od normalności
df <- data.frame(res_ar, res_ma, res_arma)
ggplot(df, aes(sample = res_ar)) +
  stat_qq() +
  stat_qq_line(col='red') +
  ggtitle('Q-Q plot dla modelu AR(1)')

ggplot(df, aes(sample = res_ma)) +
  stat_qq() +
  stat_qq_line(col='red') +
  ggtitle('Q-Q plot dla modelu MA(2)')

ggplot(df, aes(sample = res_arma)) +
  stat_qq() +
  stat_qq_line(col='red') +
  ggtitle('Q-Q plot dla modelu ARMA(2,2)')


#(e)
#możemy porównać modele miedzy sobą na bazie krteriów informacyjnych
#AIC
#najlepiej wypada model MA(2), potem AR(1), potem ARMA(2,2)
fit_ar$aic
fit_ma$aic
fit_arma$aic

#AICc
#najlepiej wypada model MA(2), potem AR(1), potem ARMA(2,2)
fit_ar$aicc
fit_ma$aicc
fit_arma$aicc

#BIC
#najlepiej wypada model ARMA(2,2), potem MA(2), potem AR(1)
fit_ar$bic
fit_ma$bic
fit_arma$bic

#(f)
#Badamy istotność współczynników
#wszystkie istotne na poziomie ufności 0.05
lmtest::coeftest(fit_ar)

#wszystkie istotne na poziomie ufności 0.05
lmtest::coeftest(fit_ma)

#wszystkie istotne na poziomie ufności 0.05
lmtest::coeftest(fit_arma_aic)

#(g)
#za najlepiej dopasowany mimo wszystko uznajemy model ARMA(2,2)?

#(h)
h <- length(test)

prediction_arma <- forecast::forecast(fit_arma, h=h, fan=TRUE)
prediction_ar <- forecast::forecast(fit_ar, h=h, fan=TRUE)
prediction_ma <- forecast::forecast(fit_ma, h=h, fan=TRUE)
help('forecast')
autoplot(prediction_ar) +
  labs(x='Czas', y='Szereg gnp') +
  ggtitle('Prognoza 23-krokowa, ARIMA(1,1,0) z dryfem')
autoplot(prediction_ma) +
  labs(x='Czas', y='Szereg gnp') +
  ggtitle('Prognoza 23-krokowa, ARIMA(0,1,2) z dryfem')
autoplot(prediction_arma) +
  labs(x='Czas', y='Szereg gnp') +
  ggtitle('Prognoza 23-krokowa, ARIMA(2,1,2) z dryfem')

autoplot(ts.union(
  prognoza=prediction_arma$mean,
  `zbiór treningowy`=train,
  `zbiór testowy`=test
  ), lwd=1) +
  labs(x='Czas', y='Szereg gnp') +
  ggtitle('Prognoza 23-krokowa, ARIMA(2,1,2) z dryfem')

autoplot(ts.union(
  prognoza=prediction_ar$mean,
  `zbiór treningowy`=train,
  `zbiór testowy`=test
), lwd=1) +
  labs(x='Czas', y='Szereg gnp') +
  ggtitle('Prognoza 23-krokowa, ARIMA(1,1,0) z dryfem')

autoplot(ts.union(
  prognoza=prediction_ma$mean,
  `zbiór treningowy`=train,
  `zbiór testowy`=test
), lwd=1) +
  labs(x='Czas', y='Szereg gnp') +
  ggtitle('Prognoza 23-krokowa, ARIMA(0,1,2) z dryfem')

#spójrzmy jeszcze na metody referencyjne
#prosta średnia ruchoma
simple_ma_fit <- meanf(train, h = h, level = c(80,95), fan = TRUE)
autoplot(simple_ma_fit) +
  labs(x='Czas', y='gnp') +
  ggtitle('Prognoza 23-krokowa metodą prostej średniej ruchomej')

#metoda naiwna
naive_fit <- naive(train, h=h, level=c(80,95), fan=TRUE)
autoplot(naive_fit) +
  labs(x='Czas', y='gnp') +
  ggtitle('Prognoza 23-krokowa; metoda naiwna')

#metoda naiwna uwzgledniająca dryf
naive_drift_fit <- rwf(train, h = h, drift = TRUE, level = c(80,95), fan = T)
autoplot(naive_drift_fit) +
  labs(x='Czas', y='gnp') +
  ggtitle('Prognoza 23-krokowa; metoda naiwna z dryfem')

autoplot(ts.union(
  prognoza=naive_drift_fit$mean,
  `zbiór treningowy`=train,
  `zbiór testowy`=test
), lwd=1) +
  labs(x='Czas', y='Szereg gnp') +
  ggtitle('Prognoza 23-krokowa; metoda naiwna z dryfem')

#(i)
#najlepiej wypada AR(1)
accuracy(prediction_arma, test)
accuracy(prediction_ar, test)
accuracy(prediction_ma, test)

#robimy porównanie dokladności dla prostych metod prognozowania
accuracy(simple_ma_fit, test)
accuracy(naive_fit, test)
accuracy(naive_drift_fit, test)
