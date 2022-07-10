rm(list=ls())
library(astsa)
library(ggplot2)
library(xts)
library(forecast)
library(animation)
library(signal)
library(pracma)
library(gridExtra)
library(ggpubr)
library(timsac)
library(ts.extend)
library(xtable)
library(randtests)
library(tsoutliers)
library(quantmod)
szereg <- ts(globtemp, frequency=1, start= start(globtemp), end=end(globtemp))
szereg
#a
globtemp
length(szereg)
setwd('C:/Users/Lenovo/Desktop/Szeregi czasowe')
#Szereg zawiera wyraźny trend bez sezonowości, ponadto przyjmujemy model X_t = f(t) + Z_t,
#gdzie f(t) to deterministyczna skladowa trendu, a Z_t to losowe fluktuacje
autoplot(szereg, main='Szereg dla globtemp', lwd=1) +
  labs(x='Czas')
ggtsdisplay(szereg, main='Szereg globtemp')

#czy użyć transformacji Boxa-Coxa? Dodajmy 1 do szeregu, by nie mieć problemów z ujemnymi
#obserwacjami i prześledźmy animacje

#przeanalizujmy wybór lambdy używając animacji
lambda <- seq(from=0, to=2, by=0.2)
lambda <- lambda[lambda != 1] #usuwamy lambda=1, bo to wyjściowy szereg
oopt = ani.options(interval=0.1)

#użycie transformacji potegowej nie poprawia zachowania szeregu
for (p in 1:length(lambda)) {
  transformacja <- BoxCox(szereg+1, lambda=lambda[p])
  p <- forecast::autoplot(transformacja, lwd=1) + 
    ggtitle(paste0("lambda = ", lambda[p])) + 
    theme(plot.title=element_text(size=22), 
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +  
    labs(x="Czas", y="szereg globtemp") 
  plot(p)
  
  ani.pause()
}
ani.options(oopt) 
#przeanalizujemy teraz, jaka metoda estymacji/eliminacji daje "najlepszą stacjonarność"
#dla różnicowania oraz wielomianów (w tym trygonometrycznych)
#różnicujemy szereg
#po zróżnicowaniu jest o wiele lepiej, wariancja zostala uregulowana, a ACF zanika szybko
ggtsdisplay(diff(szereg,differences = 1), main='Różnicowanie z krotnością 1')
ggtsdisplay(diff(szereg,differences = 2), main='Różnicowanie z krotnością 2')
ggtsdisplay(diff(szereg,differences = 3), main='Różnicowanie z krotnością 3')
ggtsdisplay(diff(szereg,differences = 4), main='Różnicowanie z krotnością 4')




#teraz metody parametryczne
wielomian_tryg <- function(x,n,s) {
  x <- ts(x)
  y <- 0
  for (k in 1:n) {
    y <- y + sin(k*x/s) + cos(k*x/s)
  }
  return(y)
}

trend_st1 <- tslm(szereg~poly(trend, degree=1, raw=TRUE))$fitted
trend_st2 <- tslm(szereg~poly(trend, degree=2, raw=TRUE))$fitted
trend_st3 <- tslm(szereg~poly(trend, degree=3, raw=TRUE))$fitted
trend_tryg <- tslm(szereg~I(wielomian_tryg(trend,3, 100)))$fitted

ggtsdisplay(diff(szereg,differences = 1), main='Różnicowanie z krotnością 1')

ggtsdisplay(szereg-trend_st1, main='Eliminacja trendu wielomianowego')
ggtsdisplay(szereg-trend_st2, main='Eliminacja trendu wielomianowego')
ggtsdisplay(szereg-trend_st3, main='Eliminacja trendu wielomianowego')
ggtsdisplay(szereg-trend_tryg,
                   main='Eliminacja trendu wielomianowego trygonometrycznego')

autoplot(ts.union(szereg, trend_tryg))

#za szereg stacjonarny wybieramy 1-krotnie zróżnicowany szereg

#(b)
#wybór rzedu modelu autoregresji na podstawie estymatora funkcji PACF
#przyjmujemy rząd = 3
ggtsdisplay(diff(szereg,differences = 1), main='Różnicowanie z krotnością 1')
szereg_diff <- diff(szereg, differences = 1)

#dopasowujemy model autoregresji, aic=TRUE, ponieważ wybieramy model na podstawie tego
#kryterium
model_1 <- ar(szereg_diff, aic=TRUE, method='yule-walker', se.fit=TRUE)
df_1<- data.frame(c(seq(1, model_1$order.max, by=1)), model_1$aic[2:22])
colnames(df_1) <- c('x', 'AIC')

#Uwaga: na wykresie są przedstawione różnice pomiedzy AIC, a nie ich wartości!
#AIC także wskazuje na rząd równy 2
ggplot(df_1, aes(x = x, y = AIC)) +
  geom_point() +
  ggtitle('Wybór rzędu modelu autoregresji poprzez kryterium AIC') +
  scale_x_continuous(breaks=1:21)
model_1$order
#kryterium FPE - trzeba zrobić
fpe <- function(time_series, p.max, method) {
  var_hat <- rep(0, p.max)
  fpe_values <- rep(0, p.max)
  for(p in 1:p.max) {
    var_hat <- ar(time_series, FALSE, p, method=method)$var.pred
    fpe_values[p] <- var_hat*(length(time_series)+p)/(length(time_series)-p)
  }
  return(fpe_values)
}
#rząd modelu wybrany przez kryterium FPE to 2, choć niewiele brakuje, by wybralo rząd równy 1
which.min(fpe(szereg_diff, p.max=model_1$order.max, method='yule-walker'))
fpe(szereg_diff, p.max=model_1$order.max, method='yule-walker')[3]
df_2 <- data.frame(c(seq(1,model_1$order.max, by=1)),
                   fpe(szereg_diff, model_1$order.max, 'yule-walker'))
colnames(df_2) <- c('x', 'FPE')
ggplot(df_2, aes(x = x, y = FPE)) +
  geom_point() +
  ggtitle('Wybór rzędu modelu autoregresji poprzez kryterium FPE') +
  scale_x_continuous(breaks=1:21)

#(c)
#wyznaczamy estymatory parametrów modelu autoregresji uzyskanego metodą Yule'a-Walkera
yw_par <- model_1$ar
model_2 <- ar(szereg_diff, aic=TRUE, method='burg', se.fit=TRUE)
burg_par <- model_2$ar
model_3 <- ar(szereg_diff, aic=TRUE, method='mle', se.fit=TRUE)
mle_par <- model_3$ar
xtable(data.frame(yw_par, burg_par, mle_par))

#rząd modelu to także 3 dla model_2 i model_3
cov_yw <- model_1$asy.var.coef
cov_burg <- model_2$asy.var.coef
cov_mle <- model_3$asy.var.coef

#porównanie wariancji uzysaknych estymatorów
df <- data.frame(c(cov_yw[1,1], cov_yw[2,2], cov_yw[3,3]),
           c(cov_burg[1,1], cov_burg[2,2], cov_burg[3,3]),
           c(cov_mle[1,1], cov_mle[2,2], cov_mle[3,3]))

colnames(df)[1] <- 'war_yw'
colnames(df)[2] <- 'war_burg'
colnames(df)[3] <- 'war_mle'
xtable(df, digits=6)

#(d)
#X_t - przyczynowy szereg o średniej 0, średnia jest odjeta automatyczne podczas dopasowywania
#modelu autoregresji, przyczynowość sprawdzamy poniżej


#liczymy statystyki lewo- i prawostronne dla przedzialów ufności dla trzech wspólczynników
#modelu autoregresji
alfa <- 0.05
L_1 <- yw_par[1] - qnorm(1-alfa/2, 0,1)*sqrt(model_1$asy.var.coef[1,1])
R_1 <- yw_par[1] + qnorm(1-alfa/2, 0,1)*sqrt(model_1$asy.var.coef[1,1])
L_2 <- yw_par[2] - qnorm(1-alfa/2, 0,1)*sqrt(model_1$asy.var.coef[2,2])
R_2 <- yw_par[2] + qnorm(1-alfa/2, 0,1)*sqrt(model_1$asy.var.coef[2,2])
L_3 <- yw_par[3] - qnorm(1-alfa/2, 0,1)*sqrt(model_1$asy.var.coef[3,3])
R_3 <- yw_par[3] + qnorm(1-alfa/2, 0,1)*sqrt(model_1$asy.var.coef[3,3])

#wpisujemy przedzialy ufności w data frame
#estymowane wspolczynniki należą do przedzialów ufności -> wspólczynniki są istotne!
CI <- data.frame(c(L_1, L_2, L_3), c(R_1, R_2,R_3), rep(0,3))
colnames(CI) <- c('L', "R", 'wsp.')

#zera nie należą do przedzialów, wiec wszystkie wspólczynniki są istotne
CI
xtable(CI, digits=4)

#(e)
reszty <- model_1$resid

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

#Uwaga:zmieniamy funkcje theor_test, by uwzglednić liczbe stopni swobody
theor_test <- function(ts, p){
  if(NROW(ts)<=500){
    test <-Box.test(ts, lag=log(NROW(ts)), type='Ljung-Box', fitdf = p)
  }
  else{
    if(NROW(ts)>750){
      test <-Box.test(ts,lag=20,type='Ljung-Box', fitdf=p)
    }
    else{
      test <-Box.test(ts,lag=10,type='Ljung-Box', fitdf=p)
    }
  }
  return(test)
}

#możemy uznać, że reszty są bialym szumem, wiec model dobrze dopasowany
graph_test(reszty)
theor_test(reszty, model_1$order)

#test oparty na dopasowaniu do reszt modelu autoregresji - jeżeli rząd = 0, to zakladamy, że
#mamy bialy szum, uwaga: 3 pierwsze wartości reszt to NA
ar(reszty[4:length(reszty)], aic=TRUE, method='yule-walker')$order
ar(reszty[4:length(reszty)], aic=TRUE, method='mle')$order

#test znaków
difference.sign.test(reszty)
help('difference.sign.test')

#sprawdźmy jeszcze normalność reszt
#test Jarque'a-Bera
JarqueBera.test(reszty)

#test Shapiro-Wilka
shapiro.test(reszty)

#wykres kwantylowy
d <- data.frame(Group=rep(1:2, each=length(reszty[4:length(reszty)])),
            Sample=c(rnorm(length(reszty[4:length(reszty)]), 0, sqrt(model_1$var.pred)),
                     reszty[4:length(reszty)]))

qplot(sample=Sample, data=d, color=as.factor(Group)) +
  ggtitle('Porównanie wykresów kwantylowych')

#(f)

#prognoza 5-krokowa
prognoza_diff_5 <- predict(model_1, n.ahead = 5)$pred

#odwracamy prognoze
prognoza_5 <- diffinv(prognoza_diff_5, xi = globtemp[length(globtemp)])
autoplot(cbind(globtemp, prognoza_5), lwd=1) +
  labs(x='Czas', y='Szereg globtemp + prognoza') +
  ggtitle('Prognoza 5-krokowa')
as.numeric(df_1$L)
#prognoza 10-krokowa
prognoza_diff_10 <- predict(model_1, n.ahead = 10)$pred

#odwracamy prognoze
prognoza_10 <- diffinv(prognoza_diff_10, xi = globtemp[length(globtemp)])
autoplot(cbind(globtemp, prognoza_10), lwd=1) +
  labs(x='Czas', y='Szereg globtemp + prognoza') +
  ggtitle('Prognoza 10-krokowa')

#prognoza 15-krokowa
prognoza_diff_15 <- predict(model_1, n.ahead = 15)$pred

#odwracamy prognoze
prognoza_15 <- diffinv(prognoza_diff_15, xi = globtemp[length(globtemp)])
autoplot(cbind(globtemp, prognoza_15), lwd=1) +
  labs(x='Czas', y='Szereg globtemp + prognoza') +
  ggtitle('Prognoza 15-krokowa')


#prognoza 20-krokowa
prognoza_diff_20 <- predict(model_1, n.ahead = 20)$pred

#odwracamy prognoze
prognoza_20 <- diffinv(prognoza_diff_20, xi = globtemp[length(globtemp)])
autoplot(cbind(globtemp, prognoza_20), lwd=1) +
  labs(x='Czas', y='Szereg globtemp + prognoza') +
  ggtitle('Prognoza 20-krokowa')

#porównujemy wyniki otrzymane powyżej z prostymi metodami prognozowania

#prosta średnia ruchoma
autoplot(meanf(szereg, h = 10, level = c(80,95), fan = TRUE, lambda = NULL)) +
  labs(x='Czas', y='Globtemp') +
  ggtitle('Prognoza 10-krokowa metodą prostej średniej ruchomej')

#metoda naiwna
help('naive')
autoplot(naive(szereg, h=10, level=c(80,95), fan=TRUE, lambda=NULL)) +
  labs(x='Czas', y='Globtemp') +
  ggtitle('Prognoza 10-krokowa; metoda naiwna')

#metoda naiwna uwzgledniająca dryf
autoplot(rwf(szereg, h = 10, drift = TRUE, level = c(80,95), fan = FALSE,
              lambda = NULL, tidy = FALSE)) +
  labs(x='Czas', y='Globtemp') +
  ggtitle('Prognoza 10-krokowa metodą naiwną z dryfem')
