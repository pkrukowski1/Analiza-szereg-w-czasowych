rm(list=ls())
#potrzebne biblioteki
library(xts)
library(forecast)
library(ggplot2)
library(animation)
library(expsmooth)
library(signal)
library(smoots)
library(pracma)
library(mFilter)
library(splines)
library(gridExtra)
library(aTSA)

#Wczytujemy dane dotyczące naszej grupy wraz z odpowiednim losowaniem dekady
losuj.dekade <- function(album)
{
  set.seed(album)
  dekady <- seq(from=30, to=70, by=10)
  sample(dekady, 1)
}

losuj.dekade(249823 + 249824) #nasza dekada to 60
setwd('C:/Users/Lenovo/Desktop/Szeregi czasowe') #Ustawiamy ścieżkę do wczytania danych
dane <- read.table("Dow_60-69.html", blank.lines.skip = T, fill = T, sep="\t", header=T, skip=2)
dane <- dane[-nrow(dane),] # usuwamy ostatni wiersz (znacznik </PRE>)
kurs <- dane$closing.values
daty <- as.POSIXlt(paste0("19",as.character(dane$Date)), format = "%Y%m%d")

# obiekt klasy 'xts' (dokładne daty notowań)
szereg.xts <- xts(kurs, order.by=daty)
anyNA(szereg.xts) #Brak NA
plot(szereg.xts, main='Notowania indeksu D&J w latach 1960-1969')

# obiekt klasy 'ts' (standardowy format w R dla szeregów czasowych)
szereg.ts <- ts(kurs, start=1960, end=1969,  frequency=292)
head(szereg.ts)
autoplot(szereg.ts, main='Notowania indeksu D&J w latach 1960-1969')

#(a) dla pelnej dekady
#przeksztalcamy dane, używając transformacji logarytmicznej
length(szereg.ts)
szereg.ts.BoxCox <- BoxCox(szereg.ts, lambda=0)
autoplot(szereg.ts.BoxCox, main='Szereg notowań po transformacji BoxaCoxa, lambda=0')
BoxCox.lambda(szereg.ts) #lambdą wybraną przez R jest ~ 1.89
#przeanalizujmy wybór lambdy używając animacji
lambda <- seq(from=0, to=2, by=0.2)
lambda <- lambda[lambda != 1] #usuwamy lambda=1, bo to wyjściowy szereg

for (p in 1:length(lambda)) {
  transformacja <- BoxCox(szereg.ts, lambda=lambda[p])
  p <- forecast::autoplot(transformacja, lwd=1) + 
    ggtitle(paste0("lambda = ", lambda[p])) + 
    theme(plot.title=element_text(size=22), 
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +  
    labs(x="Czas", y="szereg DJ") 
  plot(p)
  
  ani.pause()
}
#Zatem lambda=0 to najlepszy wybór, ponieważ wariancja jest najbardziej jednorodna

#(b) dla pelnej dekady
#zakladamy, że X_t = f(t) + Z_t
#Estymujemy trend

#####
#różnicowanie
#####
#Różnicowanie z krotnością 1 jest już satysfakcjonujące
ggtsdisplay(diff(szereg.ts.BoxCox,differences = 1), main='Różnicowanie z krotnością 1')
ggtsdisplay(diff(szereg.ts.BoxCox,differences = 2), main='Różnicowanie z krotnością 2')
ggtsdisplay(diff(szereg.ts.BoxCox,differences = 3), main='Różnicowanie z krotnością 3')
ggtsdisplay(diff(szereg.ts.BoxCox,differences = 4), main='Różnicowanie z krotnością 4')
help('diff')
#szereg.ts.BoxCox.diff to szereg czasowy z wyeliminowanym trendem
szereg.ts.BoxCox.diff <- diff(szereg.ts.BoxCox,differences = 1) 

#####
#Wygładzanie symetryczną ruchomą średnią
#####

order <- seq(from=1, to=51, by=2)

for (p in 1:length(order)) {
  ruchoma.srednia <- forecast::ma(szereg.ts.BoxCox,order=order[p])
  szeregi <- ts.union(szereg.ts.BoxCox, ruchoma.srednia)
  p <- forecast::autoplot(szeregi, lwd=1) + 
    ggtitle(paste0("Ruchoma średnia: order = ", order[p])) + 
    theme(plot.title=element_text(size=22), 
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +  
    ylab("szereg DJ") 
  plot(p)
  
  ani.pause()
}
#Rząd = 31 dość dobrze wygładził szereg, zachowując informacje o lokalnych wahaniach, zatem za
#estymator trendu przyjmujemy
trend.ma.31 <- forecast::ma(szereg.ts.BoxCox, order=31)
autoplot(ts.union(trend.ma.31, szereg.ts.BoxCox),
         main='Wyestymowany trend metodą symetrycznej średniej ruchomej z q=31')

#Eliminujemy trend, nie trzeba tego robić w zadaniach
resid.ma.31 <- szereg.ts.BoxCox-trend.ma.31 #Uwaga na NA!
ggtsdisplay(resid.ma.31, main='Szereg z wyeliminowanym trendem przy użyciu sym. średniej ruchomej')
####
#Filtr Spencera
####
spencer.DJ <- spencer(szereg.ts.BoxCox) #estymowany trend
# Uwaga: konieczna jest konwersja na obiekt 'ts'
spencer.DJ <- ts(spencer.DJ, frequency=292, start=start(szereg.ts.BoxCox), 
                 end=end(szereg.ts.BoxCox))

autoplot(ts.union(szereg.ts.BoxCox, spencer.DJ), lwd=1, main='Wygładzanie szeregu DJ przy użyciu filtru Spencera') +
  ylab('Szereg DJ')
autoplot(spencer.DJ)

#Eliminujemy trend
resid.spencer <- szereg.ts.BoxCox - spencer.DJ
#Możemy przyjąć, że jest to szereg stacjonarny
ggtsdisplay(resid.spencer, main='Szereg z wyeliminowanym trendem przy użyciu filtru Spencera')
Box.test(resid.spencer) #na poziomie istotności alpha=0.05 przyjmujemy, że nie mamy do czynienia z białym
                        #szumem

####
# Wygładzanie wykładnicze (exponential smoothing)
####
ses.alfa01 <- ses(szereg.ts.BoxCox, alpha=0.1)$fitted
ses.alfa04 <- ses(szereg.ts.BoxCox, alpha=0.4)$fitted
ses.alfa.optim <- ses(szereg.ts.BoxCox)$fitted
dane.ses <- ts.union(szereg.ts.BoxCox, ses.alfa.0.1=ses.alfa01, ses.alfa.0.4=ses.alfa04, ses.alfa.optymalne=ses.alfa.optim)
autoplot(dane.ses, main="Wygładzanie wykładnicze", lwd=.5)
autoplot(ts.union(szereg.ts.BoxCox, ses.alfa01), main='Wygładzanie wykładnicze')
autoplot(ses.alfa01) #Dla alfa = 0.1 jest najlepsze wygładzanie

###
#Jądrowy estymator typu Nadaraya-Watsona
###
help('ksmooth')
help('knsmooth')
#Jakie jądro bedzie najlepsze? Wybieramy poprzez zmiane parametru mu
#Oznaczenia 0 - jądro, 0.15 - bandwith
trend.nw.0_0.15 <- knsmooth(szereg.ts.BoxCox, mu=0, b=0.15)
trend.nw.0_0.15 <- trend.nw.0_0.15$ye
trend.nw.0_0.15 <- ts(trend.nw.0_0.15,start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)


trend.nw.1_0.15 <- knsmooth(szereg.ts.BoxCox, mu=1, b=0.15)
trend.nw.1_0.15 <- trend.nw.1_0.15$ye
trend.nw.1_0.15 <- ts(trend.nw.1_0.15, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)

trend.nw.2_0.15 <- knsmooth(szereg.ts.BoxCox, mu=2, b=0.15)
trend.nw.2_0.15 <- trend.nw.2_0.15$ye
trend.nw.2_0.15 <- ts(trend.nw.2_0.15, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)

trend.nw.3_0.15 <- knsmooth(szereg.ts.BoxCox, mu=3, b=0.15)
trend.nw.3_0.15 <- trend.nw.3_0.15$ye
trend.nw.3_0.15 <- ts(trend.nw.3_0.15, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)

szeregi_trend_kernel_0.15 <- ts.union(szereg.ts.BoxCox,
                                      Uniform=trend.nw.0_0.15,
                                      Epanechnikov=trend.nw.1_0.15,
                                      Bisquare=trend.nw.2_0.15,
                                      Triweight=trend.nw.3_0.15)
                                      

autoplot(szeregi_trend_kernel_0.15, lwd=1, main='Estymacja jądrowa') +
  labs(x='Czas', y='Szereg D-J')
autoplot(trend.nw.0_0.15) #musimy lepiej kontrolować parametr bandwith


trend.nw.0_0.05 <- knsmooth(szereg.ts.BoxCox, mu=0, b=0.05)
trend.nw.0_0.05 <- trend.nw.0_0.05$ye
trend.nw.0_0.05 <- ts(trend.nw.0_0.05,start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)


trend.nw.1_0.05 <- knsmooth(szereg.ts.BoxCox, mu=1, b=0.05)
trend.nw.1_0.05 <- trend.nw.1_0.05$ye
trend.nw.1_0.05 <- ts(trend.nw.1_0.05, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)

trend.nw.2_0.05 <- knsmooth(szereg.ts.BoxCox, mu=2, b=0.05)
trend.nw.2_0.05 <- trend.nw.2_0.05$ye
trend.nw.2_0.05 <- ts(trend.nw.2_0.05, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)

trend.nw.3_0.05 <- knsmooth(szereg.ts.BoxCox, mu=3, b=0.05)
trend.nw.3_0.05 <- trend.nw.3_0.05$ye
trend.nw.3_0.05 <- ts(trend.nw.3_0.05, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)


szeregi_trend_kernel_0.05 <- ts.union(szereg.ts.BoxCox, 
                                      Uniform=trend.nw.0_0.05, 
                                      Epanechnikov=trend.nw.1_0.05,
                                      Bisquare=trend.nw.2_0.05,
                                      Triweight=trend.nw.3_0.05)

autoplot(szeregi_trend_kernel_0.05, lwd=1, main='Estymacja jądrowa') +
  labs(x='Czas', y='Szereg D-J')

trend.nw.0_0.01 <- knsmooth(szereg.ts.BoxCox, mu=0, b=0.01)
trend.nw.0_0.01 <- trend.nw.0_0.01$ye
trend.nw.0_0.01 <- ts(trend.nw.0_0.01,start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)


trend.nw.1_0.01 <- knsmooth(szereg.ts.BoxCox, mu=1, b=0.01)
trend.nw.1_0.01 <- trend.nw.1_0.01$ye
trend.nw.1_0.01 <- ts(trend.nw.1_0.01, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)

trend.nw.2_0.01 <- knsmooth(szereg.ts.BoxCox, mu=2, b=0.01)
trend.nw.2_0.01 <- trend.nw.2_0.01$ye
trend.nw.2_0.01 <- ts(trend.nw.2_0.01, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)

trend.nw.3_0.01 <- knsmooth(szereg.ts.BoxCox, mu=3, b=0.01)
trend.nw.3_0.01 <- trend.nw.3_0.01$ye
trend.nw.3_0.01 <- ts(trend.nw.3_0.01, start=start(szereg.ts.BoxCox), 
                      end=end(szereg.ts.BoxCox), frequency = 292)



szeregi_trend_kernel_0.01 <- ts.union(szereg.ts.BoxCox, 
                                      Uniform=trend.nw.0_0.01, 
                                      Epanechnikov=trend.nw.1_0.01,
                                      Bisquare=trend.nw.2_0.01,
                                      Triweight=trend.nw.3_0.01)
                                     
autoplot(szeregi_trend_kernel_0.01, lwd=1, main='Estymacja jądrowa') +
  labs(x='Czas', y='Szereg D-J')
#za estymator trendu przyjmujemy trend.nw.3_0.01
autoplot(trend.nw.3_0.01, main = 'Estymator trendu metodą jądrową typu Nadaraya-Watsona',lwd=1) +
  labs(x='Czas', y='Szereg D-J')

###
#suma funkcji bazowych B-sklejanych 
###
help('ns')
for (p in seq(10,70,by=10)) {
  trend.funkcja.sklejana <- tslm(szereg.ts.BoxCox~ns(trend, p))$fitted
  modele <- ts.union(trend.funkcja.sklejana, szereg.ts.BoxCox)
  g <- autoplot(modele) +
    ggtitle(paste0("Bazowe funkcje B-sklejane: order = ", p))
  plot(g)
  ani.pause()
}
#przyjmujemy order=30
trend.funkcja.sklejana <- tslm(szereg.ts.BoxCox~ns(trend, 30))$fitted

######
#c dla pelnej dekady
######


model.liniowy <- tslm(szereg.ts.BoxCox~trend)
summary(model.liniowy) 
trend.liniowy <- model.liniowy$fitted
trend.kwadrat <- tslm(szereg.ts.BoxCox~I(trend^2)+trend)$fitted
trend.st5 <- tslm(szereg.ts.BoxCox~poly(trend, degree=5, raw=TRUE))$fitted
trend.st6 <- tslm(szereg.ts.BoxCox~poly(trend, degree=6, raw=TRUE))$fitted

wielomian.tryg <- function(x,n,s) {
  x <- ts(x)
  y <- 0
  for (k in 1:n) {
    y <- y + sin(k*x/s) + cos(k*x/s)
  }
  return(y)
}

trend.funkcja.wiel.tryg <- tslm(szereg.ts.BoxCox~I(wielomian.tryg(trend,4, 800)))$fitted
dane.trend <- ts.union(szereg.ts.BoxCox, trend.liniowy, trend.kwadrat, trend.st5, trend.st6,
                       trend.funkcja.wiel.tryg)
autoplot(dane.trend, lwd=.75, main="Trend wielomianowy")

#wielomian trygonometryczny nie dziala najlepiej
#sprawdźmy, czy nie mamy do czynienia z przeparametryzowaniem
tslm.1 <- tslm(szereg.ts.BoxCox~poly(trend, degree=5, raw=TRUE))
tslm.2 <- tslm(szereg.ts.BoxCox~I(wielomian.tryg(trend,4, 800)))

# wyznaczamy prognozy na podstawie dopasowanego trendu wielomianowego
prognozy.1 <- forecast::forecast(tslm.1, h=5, level = c(80, 95))
prognozy.2 <- forecast::forecast(tslm.2, h=5, level = c(80, 95))

p1 <- autoplot(prognozy.1, main="Prognozy na podstawie trendu wielomianu stopnia 5.", flwd = 1)
p2 <- autoplot(prognozy.2, main="Prognozy na podstawie trendu wielomianu trygonometrycznego", flwd = 1)

#Estymowany trend wielomianem stopnia 5. jest lepszy
gridExtra::grid.arrange(p1,p2, nrow=2)

############################################
# Przeprowadźmy jeszcze raz powyższą analize dla roku 1965
############################################
daty.1965 <- daty[as.Date('1965-01-01 CET') <= as.Date(daty) & as.Date(daty)
                  <= as.Date("1965-12-31 CET")]
kurs.1965 <- kurs[as.Date('1965-01-01 CET') <= as.Date(daty) & as.Date(daty)
                  <= as.Date("1965-12-31 CET")]
szereg.ts.1965 <- ts(kurs.1965, frequency = 292, start=c(1965,1), end=c(1966,1))

szereg.xts.1965 <- xts(kurs.1965, order.by=daty.1965, frequency = 292)
autoplot(szereg.ts.1965, main="Notowania indeksu D&J w roku 1965", lwd=1) +
  labs(x="Czas", y="Szereg D-J w 1965 roku")

#(a) dla roku 1965
#przeksztalcamy dane, używając transformacji logarytmicznej

szereg.ts.1965.BoxCox <- BoxCox(szereg.ts.1965, lambda=0)
autoplot(szereg.ts.1965.BoxCox, 
         main='Szereg notowań po transformacji BoxaCoxa, lambda=0') +
  labs(x="Czas", y="Szereg D-J w 1965 roku")
BoxCox.lambda(szereg.ts.1965) #lambdą wybraną przez R jest ~ 1, czyli wyjściowy szereg
#przeanalizujmy wybór lambdy używając animacji


for (p in 1:length(lambda)) {
  transformacja.1965 <- BoxCox(szereg.ts.1965, lambda=lambda[p])
  g <- forecast::autoplot(transformacja.1965, lwd=1) + 
    ggtitle(paste0("lambda = ", lambda[p])) + 
    theme(plot.title=element_text(size=22), 
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +  
    labs(x = "Czas", y="Szereg D-J w roku 1965") 
  plot(g)
  
  ani.pause()
}
r#Zatem lambda=0 to najlepszy wybór, ponieważ wariancja jest najbardziej jednorodna

#(b) dla roku 1965
#zakladamy, że X_t = f(t) + Z_t
#Estymujemy trend

#####
#różnicowanie
#####
#Różnicowanie z krotnością 1 jest już satysfakcjonujące
ggtsdisplay(diff(szereg.ts.1965.BoxCox,differences = 1), main='Różnicowanie z krotnością 1')
ggtsdisplay(diff(szereg.ts.1965.BoxCox,differences = 2), main='Różnicowanie z krotnością 2')
ggtsdisplay(diff(szereg.ts.1965.BoxCox,differences = 3), main='Różnicowanie z krotnością 3')
ggtsdisplay(diff(szereg.ts.1965.BoxCox,differences = 4), main='Różnicowanie z krotnością 4')

#szereg.ts.BoxCox.diff to szereg czasowy z wyeliminowanym trendem
szereg.ts.1965.BoxCox.diff <- diff(szereg.ts.1965.BoxCox,differences = 1) 

#test Boxa mówi nam o tym, że mamy do czynienia z bialym szumem
Box.test(szereg.ts.1965.BoxCox.diff) 


#####
#Wygładzanie symetryczną ruchomą średnią
#####

order <- seq(from=1, to=51, by=2)

for (p in 1:length(order)) {
  ruchoma.srednia.1965 <- forecast::ma(szereg.ts.1965.BoxCox,order=order[p])
  szeregi.1965 <- ts.union(szereg.ts.1965.BoxCox, ruchoma.srednia.1965)
  p <- forecast::autoplot(szeregi.1965, lwd=1) + 
    ggtitle(paste0("Ruchoma średnia: order = ", order[p])) + 
    theme(plot.title=element_text(size=18), 
          legend.text=element_text(size=10),
          legend.title=element_text(size=10)) +  
    labs(x="CZas", y="Szereg D-J w 1965 roku")
  plot(p)
  
  ani.pause()
}
#Rząd = 25 dość dobrze wygładził szereg, choć metoda ta nie dziala tak dobrze, jak przy dluższym okresie
#czasowym - tracimy sporo informacji o wahaniach lokalnych
trend.ma.1965.25 <- forecast::ma(szereg.ts.1965.BoxCox, order=25)
autoplot(trend.ma.1965.25, 
         main='Wyestymowany trend metodą symetrycznej średniej ruchomej z q=25',lwd=1) +
  labs(x="Czas",y='Szereg DJ w 1965 roku')

#Eliminujemy trend
resid.ma.1965.25 <- szereg.ts.1965.BoxCox-trend.ma.1965.25 #Uwaga na NA!

ggtsdisplay(resid.ma.1965.25, main='Szereg z wyeliminowanym trendem przy użyciu sym. średniej ruchomej')

####
#Filtr Spencera
####

spencer.DJ.1965 <- spencer(szereg.ts.1965.BoxCox) #estymowany trend
# Uwaga: konieczna jest konwersja na obiekt 'ts'
spencer.DJ.1965 <- ts(spencer.DJ.1965, frequency=292, start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965))

autoplot(ts.union(szereg.ts.1965.BoxCox, spencer.DJ.1965), lwd=1,
         main='Wygładzanie szeregu DJ przy użyciu filtru Spencera') +
  labs(y='Szereg D-J w 1965 roku', x='Czas')

#Eliminujemy trend
resid.spencer.1965 <- szereg.ts.1965.BoxCox - spencer.DJ.1965


####
# Wygładzanie wykładnicze (exponential smoothing)
####
ses.alfa01.1965 <- ses(szereg.ts.1965.BoxCox, alpha=0.1)$fitted
ses.alfa02.1965 <- ses(szereg.ts.1965.BoxCox, alpha=0.2)$fitted
ses.alfa04.1965 <- ses(szereg.ts.1965.BoxCox, alpha=0.4)$fitted
ses.alfa.optim.1965 <- ses(szereg.ts.1965.BoxCox)$fitted
dane.ses.1965 <- ts.union(szereg.ts.1965.BoxCox, ses.alfa.0.1=ses.alfa01.1965,
                          ses.alfa.0.4=ses.alfa04.1965,
                          ses.alfa.optymalne=ses.alfa.optim.1965,
                          ses.alfa.0.2=ses.alfa02.1965)

autoplot(dane.ses.1965, main="Wygładzanie wykładnicze", lwd=.8)
autoplot(ts.union(szereg.ts.1965.BoxCox, ses.alfa01.1965),
         main='Wygładzanie wykładnicze',
         lwd=1) +
  labs(x='Czas', y='Szereg D-J w 1965 roku')
###
#Jądrowy estymator typu Nadaraya-Watsona
###

#Jakie jądro bedzie najlepsze? Wybieramy poprzez zmiane parametru mu
#Oznaczenia 0 - jądro, 0.15 - bandwith
trend.1965.nw.0_0.15 <- knsmooth(szereg.ts.1965.BoxCox, mu=0, b=0.15)
trend.1965.nw.0_0.15 <- trend.1965.nw.0_0.15$ye
trend.1965.nw.0_0.15 <- ts(trend.1965.nw.0_0.15,start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)


trend.1965.nw.1_0.15 <- knsmooth(szereg.ts.1965.BoxCox, mu=1, b=0.15)
trend.1965.nw.1_0.15 <- trend.1965.nw.1_0.15$ye
trend.1965.nw.1_0.15 <- ts(trend.1965.nw.1_0.15, start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)

trend.1965.nw.2_0.15 <- knsmooth(szereg.ts.1965.BoxCox, mu=2, b=0.15)
trend.1965.nw.2_0.15 <- trend.1965.nw.2_0.15$ye
trend.1965.nw.2_0.15 <- ts(trend.1965.nw.2_0.15, start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)

trend.1965.nw.3_0.15 <- knsmooth(szereg.ts.1965.BoxCox, mu=3, b=0.15)
trend.1965.nw.3_0.15 <- trend.1965.nw.3_0.15$ye
trend.1965.nw.3_0.15 <- ts(trend.1965.nw.3_0.15, start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)


szeregi_trend_kernel_0.15_1965 <- ts.union(szereg.ts.1965.BoxCox,
                                      Uniform=trend.1965.nw.0_0.15,
                                      Epanechnikov=trend.1965.nw.1_0.15,
                                      Bisquare=trend.1965.nw.2_0.15,
                                      Triweight=trend.1965.nw.3_0.15)
                                      

#musimy lepiej kontrolować parametr bandwith
autoplot(szeregi_trend_kernel_0.15_1965, lwd=1) +
  labs(x='Czas', y='Szereg DJ w 1965 roku')


trend.1965.nw.0_0.05 <- knsmooth(szereg.ts.1965.BoxCox, mu=0, b=0.05)
trend.1965.nw.0_0.05 <- trend.1965.nw.0_0.05$ye
trend.1965.nw.0_0.05 <- ts(trend.1965.nw.0_0.05,start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)


trend.1965.nw.1_0.05 <- knsmooth(szereg.ts.1965.BoxCox, mu=1, b=0.05)
trend.1965.nw.1_0.05 <- trend.1965.nw.1_0.05$ye
trend.1965.nw.1_0.05 <- ts(trend.1965.nw.1_0.05, start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)

trend.1965.nw.2_0.05 <- knsmooth(szereg.ts.1965.BoxCox, mu=2, b=0.05)
trend.1965.nw.2_0.05 <- trend.1965.nw.2_0.05$ye
trend.1965.nw.2_0.05 <- ts(trend.1965.nw.2_0.05, start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)

trend.1965.nw.3_0.05 <- knsmooth(szereg.ts.1965.BoxCox, mu=3, b=0.05)
trend.1965.nw.3_0.05 <- trend.1965.nw.3_0.05$ye
trend.1965.nw.3_0.05 <- ts(trend.1965.nw.3_0.05, start=start(szereg.ts.1965), 
                      end=end(szereg.ts.1965), frequency = 292)


szeregi_trend_kernel_0.05_1965 <- ts.union(szereg.ts.1965.BoxCox,
                                           Uniform=trend.1965.nw.0_0.05,
                                           Epanechnikov=trend.1965.nw.1_0.05,
                                           Bisquare=trend.1965.nw.2_0.05,
                                           Triweight=trend.1965.nw.3_0.05)
                                          

autoplot(szeregi_trend_kernel_0.05_1965, lwd=1,
         main='Estymacja trendu estymatorem jądrowym typu N-W') +
  labs(x='Czas', y='Szereg D-J w 1965 roku')
  
#za estymator trendu w tym przypadku przyjmujemy trend.1965.nw.3_0.05
autoplot(ts.union(Triweight=trend.1965.nw.3_0.05, szereg.ts.1965.BoxCox),
         main = 'Estymator jądrowy typu Nadaraya-Watsona', lwd=1) +
  labs(x='Czas', y='Szereg D-J w 1965 roku')

###
#suma funkcji bazowych B-sklejanych 
###

for (p in seq(10,70,by=10)) {
    trend.funkcja.sklejana.1965 <- tslm(szereg.ts.1965.BoxCox~ns(trend, p))$fitted
    modele <- ts.union(trend.funkcja.sklejana.1965, szereg.ts.1965.BoxCox)
    g <- autoplot(modele, lwd=1) +
      ggtitle(paste0("Bazowe funkcje B-sklejane: order = ", p)) +
      labs(x="Czas", y="Szereg D-J w 1965 roku")
    plot(g)
    
    ani.pause()
    }
#przyjmujemy order=30
trend.funkcja.sklejana.1965 <- tslm(szereg.ts.1965.BoxCox~ns(trend, 30))$fitted

######
#c dla roku 1965
######

trend.cubic.1965 <- tslm(szereg.ts.1965.BoxCox~I(trend^3))$fitted
trend.st5.1965 <- tslm(szereg.ts.1965.BoxCox~poly(trend, degree=5,
                                                      raw=TRUE))$fitted

trend.st6.1965 <- tslm(szereg.ts.1965.BoxCox~poly(trend, degree=6,
                                                      raw=TRUE))$fitted

trend.st7.1965 <- tslm(szereg.ts.1965.BoxCox~poly(trend, degree=7,
                                                      raw=TRUE))$fitted
trend.funkcja.wiel.tryg.1965 <- tslm(szereg.ts.1965.BoxCox~I(wielomian.tryg(trend,2, 40)))$fitted
dane.trend.1965 <- ts.union(szereg.ts.1965.BoxCox, trend.cubic.1965, trend.st5.1965, 
                       trend.funkcja.wiel.tryg.1965, trend.st6.1965, trend.st7.1965)
autoplot(dane.trend.1965, lwd=.75, main="Trend wielomianowy") +
  labs(x='Czas', y='Szereg DJ w roku 1965')

#wielomian trygonometryczny nie dziala najlepiej
#sprawdźmy, czy nie mamy do czynienia z przeparametryzowaniem
tslm.1.1965 <- tslm(szereg.ts.1965.BoxCox~poly(trend, degree=5, raw=TRUE))
tslm.2.1965 <- tslm(szereg.ts.1965.BoxCox~poly(trend, degree=6, raw=TRUE))
tslm.3.1965 <- tslm(szereg.ts.1965.BoxCox~poly(trend, degree=7, raw=TRUE))

# wyznaczamy prognozy na podstawie dopasowanego trendu wielomianowego
prognozy.1.1965 <- forecast::forecast(tslm.1.1965, h=5)
prognozy.2.1965 <- forecast::forecast(tslm.2.1965, h=5)
prognozy.3.1965 <- forecast::forecast(tslm.3.1965, h=5)

p1.1965 <- autoplot(prognozy.1.1965,
                    main="Prognozy na podstawie trendu wielomianu stopnia 5.", flwd = 1) +
  labs(x='Czas', y='Szereg D-J')
p2.1965 <- autoplot(prognozy.2.1965,
                    main="Prognozy na podstawie trendu wielomianu stopnia 6.", flwd = 1) +
  labs(x='Czas', y='Szereg D-J')
p3.1965 <- autoplot(prognozy.3.1965,
                    main="Prognozy na podstawie trendu wielomianu stopnia 7.", flwd = 1) +
  labs(x='Czas', y='Szereg DJ')
#Estymowany trend wielomianem stopnia 6. jest najlepszy
p3.1965
# (d)

Zt_differentation_1965 <- szereg.ts.1965.BoxCox.diff
Box.test(Zt_differentation_1965,lag=floor(NROW(Zt_differentation_1965)/4),
         type='Ljung-Box')$p.value>=0.05



