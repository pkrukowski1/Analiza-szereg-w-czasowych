rm(list=ls())

library(ggplot2)
library(forecast)
library(fitdistrplus)
library(xtable)

####################################
#Zadanie 2

###################################

gen_arima_sim <- function(n, nsim, model){
  result <- matrix(nrow=n, ncol=nsim)
  for(i in 1:nsim){
    simul <- arima.sim(n=n, model=model)
    result[,i] <- simul}
  return(result)
}
gen_arima_sim(50, 500, list(ar=-0.2))

##########
#Druga kropka
##########
k <- 80 #liczba realizacji

#średnie z definicje muszą być równe 0, dla uproszczenia przeskalujmy dane (tam, gdzie to potrzebne)

# 1) realizacje bialego szumu ze zmiennymi losowymi z rozkladu N(0,1)
n_1 <- 100 #dlugość szeregu
realizacje_1 <- matrix(rnorm(n_1*k),n_1,k)

n_2 <- 200 #dlugość szeregu
realizacje_2 <- matrix(rnorm(n_2*k),n_2,k)

#2) realizacje bialego szumu ze zmiennymi losowymi z rozkladu Exp(5)
n_3 <- 50
realizacje_3 <-matrix(rexp(n_3*k,rate=5), nrow=n_3, ncol=k)
realizacje_3 <- realizacje_3-colMeans(realizacje_3)

n_4 <- 150
realizacje_4 <- matrix(rexp(n_4*k,rate=5), n_4,k)
realizacje_4 <- realizacje_4-colMeans(realizacje_4)


n_5 <- 300
realizacje_5 <- matrix(rexp(n_5*k, rate=5),n_5,k)
realizacje_5 <- realizacje_5-colMeans(realizacje_5)


########
#Trzecia kropka
########

######Wyznaczenie wartości estymatora wartości oczekiwanej - średnia próbkowa
srednie_1 <- apply(realizacje_1, MARGIN=2, FUN=mean)
srednie_2 <- apply(realizacje_2, MARGIN=2, FUN=mean)
srednie_3 <- apply(realizacje_3, MARGIN=2, FUN=mean)
srednie_4 <- apply(realizacje_4, MARGIN=2, FUN=mean)
srednie_5 <- apply(realizacje_5, MARGIN=2, FUN=mean)
#Wlasności rozkladów asymptotycznych uzyskanych estymatorów 

########
#Analiza srednie_1

#Estymator jądrowy
dens_1 <- density(srednie_1)
plot(dens_1, lwd=3, col='green', main='Porównanie estymatorów jądrowych dla srednie_1',xlab='x',ylab='y')
curve(dnorm(x,mean=0,sd=1/sqrt(n_1)), col="red", add=T, lwd=3)
grid(col='lightblue')
legend('topleft', legend=c("Estymator", "Gęstość rozkładu asymp."),
       col=c("green", "red"), lty=c(1,1), cex=0.8)

#dystrybuanta empiryczna

ecdf_1 <- ecdf(srednie_1)
plot(ecdf_1, main='Porównanie dystrybuant empirycznych dla srednie_1')
curve(pnorm(x,mean=0,sd=1/sqrt(n_1)), col="red", add=T, lwd=2)
grid(col='lightblue')
legend('bottomright', legend=c("Dystrybuanta teoretyczna", "Dystrybuanta empiryczna"),
       col=c("red", "black"), lty=1:2, cex=1)

########
#Analiza srednie_2


#Estymator jądrowy
dens_2 <- density(srednie_2)
plot(dens_2, lwd=3, col='green', main='Porównanie estymatorów jądrowych dla srednie_2',xlab='x',ylab='y')
curve(dnorm(x,mean=0,sd=1/sqrt(n_2)), col="red", add=T, lwd=3)
grid(col='lightblue')
legend('topleft', legend=c("Estymator", "Gęstość rozkładu asymp."),
       col=c("green", "red"), lty=c(1,1), cex=0.8)

#dystrybuanta empiryczna

ecdf_2 <- ecdf(srednie_2)
plot(ecdf_2, main='Porównanie dystrybuant empirycznych dla srednie_2')
curve(pnorm(x,mean=0,sd=1/sqrt(n_2)), col="red", add=T, lwd=2)
grid(col='lightblue')
legend('topleft', legend=c("Dystrybuanta teoretyczna", "Dystrybuanta empiryczna"),
       col=c("red", "black"), lty=1:2, cex=1)



########
#Analiza srednie_3

#Estymator jądrowy
dens_3 <- density(srednie_3)
plot(dens_3, lwd=3, col='green', main='Porównanie estymatorów jądrowych dla srednie_3',xlab='x',ylab='y',
    )
curve(dnorm(x,mean=0,sd=1/(5*sqrt(n_3))), col="red", add=T, lwd=3)
grid(col='lightblue')
legend('topright', legend=c("Estymator", "Gęstość rozkładu asymp."),
       col=c("green", "red"), lty=1:1, cex=0.8)

#dystrybuanta empiryczna
ecdf_3 <- ecdf(srednie_3)
plot(ecdf_3, main='Porównanie dystrybuant empirycznych dla srednie_3')
curve(pnorm(x,mean=0,sd=1/(5*sqrt(n_3))), col="red", add=T, lwd=2)
grid(col='lightblue')
legend('topleft', legend=c("Dystrybuanta teoretyczna", "Dystrybuanta empiryczna"),
       col=c("red", "black"), lty=1:2, cex=1)

#########
#Analiza srednie_4
#Estymator jądrowy
dens_4 <- density(srednie_4)
plot(dens_4, lwd=3, col='green', main='Porównanie estymatorów jądrowych dla srednie_4',xlab='x',ylab='y',
ylim=c(1,30))
curve(dnorm(x,mean=0,sd=1/(5*sqrt(n_4))), col="red", add=T, lwd=3)
grid(col='lightblue')
legend('topleft', legend=c("Estymator", "Gęstość rozkładu asymp."),
       col=c("green", "red"), lty=1:1, cex=0.8)

#dystrybuanta empiryczna
ecdf_4 <- ecdf(srednie_4)
plot(ecdf_4, main='Porównanie dystrybuant empirycznych dla srednie_4')
curve(pnorm(x,mean=0,sd=1/(5*sqrt(n_4))), col="red", add=T, lwd=2)
grid(col='lightblue')
legend('topleft', legend=c("Dystrybuanta teoretyczna", "Dystrybuanta empiryczna"),
       col=c("red", "black"), lty=1:2, cex=1)


##########
#Analiza srednie_5


#Estymator jądrowy
dens_5 <- density(srednie_5)
plot(dens_5, lwd=3, col='green', main='Porównanie estymatorów jądrowych dla srednie_5',xlab='x',ylab='y',
     ylim=c(0,40))
curve(dnorm(x,mean=0,sd=1/(5*sqrt(n_5))), col="red", add=T, lwd=3)
grid(col='lightblue')
legend('topleft', legend=c("Estymator", "Gęstość rozkładu asymp."),
       col=c("green", "red"), lty=1:1, cex=0.8)

#dystrybuanta empiryczna
ecdf_5 <- ecdf(srednie_5)
plot(ecdf_5, main='Porównanie dystrybuant empirycznych dla srednie_5')
curve(pnorm(x,mean=0,sd=1/(5*sqrt(n_5))), col="red", add=T, lwd=2)
grid(col='lightblue')
legend('topleft', legend=c("Dystrybuanta teoretyczna", "Dystrybuanta empiryczna"),
       col=c("red", "black"), lty=1:2, cex=1)

#######
#rozkład asymptotyczny estymatora autokorelacji
#dla realizacji_1

# wyznaczamy macierz ACF (w wierszach mamy realizacje acf(h) dla h=1,2,3,...,h.max)
h.max_1 <- floor(n_1/4) # maksymalne opóźnienie (można uzależnić od n)
acf.matrix_1 <- apply(realizacje_1, 2, function(x) acf(x, lag.max=h.max_1, type="correlation", plot=FALSE)$acf)

# usuwamy ACF(0)
acf.matrix_1 <- acf.matrix_1[-1,]
h.wybrane_1 <- c(1, 5, 10, 25)  # można uzależnić od n

# histogramy dla ACF(h)
par(mfrow=c(2,2))

for (h in h.wybrane_1) {
  print(paste0("h=",h))
  tytul_1 <- paste0("Histogram dla acf(",h,")")
  acf.h_1 <- acf.matrix_1[h,]
  hist(acf.h_1, freq=FALSE, col="lightblue", main=tytul_1, xlab="", ylim=c(0,6))
  curve(dnorm(x,mean=mean(acf.h_1), sd=sd(acf.h_1)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/sqrt(n_1)), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_1))
}

h.max_1
############
#dla realizacji_2

h.max_2 <- floor(n_2/4) # maksymalne opóźnienie (można uzależnić od n)
h.max_2
acf.matrix_2 <- apply(realizacje_2, 2, function(x) acf(x, lag.max=h.max_2, type="correlation", plot=FALSE)$acf)

# usuwamy ACF(0)
acf.matrix_2 <- acf.matrix_2[-1,]
h.wybrane_2 <- c(1, 10,30 ,50 )  # można uzależnić od n

# histogramy dla ACF(h)
par(mfrow=c(2,2))

for (h in h.wybrane_2) {
  print(paste0("h=",h))
  tytul_2 <- paste0("Histogram dla acf(",h,")")
  acf.h_2 <- acf.matrix_2[h,]
  hist(acf.h_2, freq=FALSE, col="lightblue", main=tytul_2, xlab="", ylim=c(0,8))
  curve(dnorm(x,mean=mean(acf.h_2), sd=sd(acf.h_2)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/sqrt(n_2)), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_2))
}

############
#dla realizacji_3

h.max_3 <- floor(n_3/4) # maksymalne opóźnienie (można uzależnić od n)
acf.matrix_3 <- apply(realizacje_3, 2, function(x) acf(x, lag.max=h.max_3, type="correlation", plot=FALSE)$acf)

# usuwamy ACF(0)
acf.matrix_3 <- acf.matrix_3[-1,]

h.wybrane_3 <- c(1, 5,8 ,12)  # można uzależnić od n

# histogramy dla ACF(h)
par(mfrow=c(2,2))

for (h in h.wybrane_3) {
  print(paste0("h=",h))
  tytul_3 <- paste0("Histogram dla acf(",h,")")
  acf.h_3 <- acf.matrix_3[h,]
  hist(acf.h_3, freq=FALSE, col="lightblue", main=tytul_3, xlab="", ylim=c(0,8))
  curve(dnorm(x,mean=mean(acf.h_3), sd=sd(acf.h_3)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/sqrt(n_3)), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_3))
}

##########
#dla realizacji_4

h.max_4 <- floor(n_4/4) # maksymalne opóźnienie (można uzależnić od n)
acf.matrix_4 <- apply(realizacje_4, 2, function(x) acf(x, lag.max=h.max_4, type="correlation", plot=FALSE)$acf)
h.max_4
# usuwamy ACF(0)
acf.matrix_4 <- acf.matrix_4[-1,]

h.wybrane_4 <- c(1, 13,28 ,34)  # można uzależnić od n

# histogramy dla ACF(h)
par(mfrow=c(2,2))

for (h in h.wybrane_4) {
  print(paste0("h=",h))
  tytul_4 <- paste0("Histogram dla acf(",h,")")
  acf.h_4 <- acf.matrix_4[h,]
  hist(acf.h_4, freq=FALSE, col="lightblue", main=tytul_4, xlab="", ylim=c(0,8))
  curve(dnorm(x,mean=mean(acf.h_4), sd=sd(acf.h_4)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/sqrt(n_4)), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_4))
}

##########
#dla realizacji_5

h.max_5 <- floor(n_5/4) # maksymalne opóźnienie (można uzależnić od n)
acf.matrix_5 <- apply(realizacje_5, 2, function(x) acf(x, lag.max=h.max_5, type="correlation", plot=FALSE)$acf)
h.max_5
# usuwamy ACF(0)
acf.matrix_5 <- acf.matrix_5[-1,]

h.wybrane_5 <- c(4, 32,55 ,70)  # można uzależnić od n

# histogramy dla ACF(h)
par(mfrow=c(2,2))

for (h in h.wybrane_5) {
  print(paste0("h=",h))
  tytul_5 <- paste0("Histogram dla acf(",h,")")
  acf.h_5 <- acf.matrix_5[h,]
  hist(acf.h_5, freq=FALSE, col="lightblue", main=tytul_5, xlab="", ylim=c(0,10))
  curve(dnorm(x,mean=mean(acf.h_5), sd=sd(acf.h_5)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/sqrt(n_5)), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_5))
}


#####
#Co sie dzieje z rozkladem dla h bliskich n? Rozważmy BUO jedynie estymator autokorelacji

acf.matrix_4_prim <- apply(realizacje_4, 2, function(x) acf(x, lag.max=n_4-1, type="correlation", 
                                                       plot=FALSE)$acf)


h.wybrane_6 <- c(n_4-20, n_4-10, n_4-5, n_4-2)

for (h in h.wybrane_6) {
  print(paste0("h=",h))
  tytul_4 <- paste0("Histogram dla acf(",h,")")
  acf.h_4_prim <- acf.matrix_4_prim[h,]
  hist(acf.h_4_prim, freq=FALSE, col="lightblue", main=tytul_4, xlab="", ylim=c(0,30))
  curve(dnorm(x,mean=mean(acf.h_4_prim), sd=sd(acf.h_4_prim)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/sqrt(n_4)), add=T, col="red", lwd=2)
  legend('topright', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.65)
}



h.wybrane_7 <- c(n_5-15, n_5-8, n_5-4, n_5-1)

acf.matrix_5_prim <- apply(realizacje_5, 2, function(x) acf(x, lag.max=n_5-1, type="correlation", 
                                                           plot=FALSE)$acf)


for (h in h.wybrane_7) {
  print(paste0("h=",h))
  tytul_5 <- paste0("Histogram dla acf(",h,")")
  acf.h_5_prim <- acf.matrix_5_prim[h,]
  hist(acf.h_5_prim, freq=FALSE, col="lightblue", main=tytul_5, xlab="", ylim=c(0,50))
  curve(dnorm(x,mean=mean(acf.h_5_prim), sd=sd(acf.h_5_prim)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/sqrt(n_5)), add=T, col="red", lwd=2)
  legend('topright', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
}


############
#rozkład asymptotyczny estymatora autokowariancji
############
#dla realizacje_1
matrix_1_autokow <- apply(realizacje_1, 2, function(x) acf(x, lag.max=h.max_1,
                                                              type="covariance", plot=FALSE)$acf)

matrix_1_autokow <- matrix_1_autokow[-1,]

par(mfrow=c(2,2))

for (h in h.wybrane_1) {
  print(paste0("h=",h))
  tytul_1 <- paste0("Histogram estymatora autokowariancji
              dla h = ",h)
  acf.h_1_autocor <- matrix_1_autokow[h,]
  hist(acf.h_1_autocor, freq=FALSE, col="lightblue", main=tytul_1, xlab="")
  curve(dnorm(x,mean=mean(acf.h_1_autocor), sd=sd(acf.h_1_autocor)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/(sqrt(n_1))), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_1_autocor))
}


#Dla realizacje_2
matrix_2_autokow <- apply(realizacje_2, 2, function(x) acf(x, lag.max=h.max_2,
                                                              type="covariance", plot=FALSE)$acf)
matrix_2_autokow <- matrix_2_autokow[-1,]
par(mfrow=c(2,2))

for (h in h.wybrane_2) {
  print(paste0("h=",h))
  tytul_2 <- paste0("Histogram estymatora autokowariancji
              dla h = ",h)
  acf.h_2_autocor <- matrix_2_autokow[h,]
  hist(acf.h_2_autocor, freq=FALSE, col="lightblue", main=tytul_2, xlab="",ylim=c(0,9))
  curve(dnorm(x,mean=mean(acf.h_2_autocor), sd=sd(acf.h_2_autocor)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/(sqrt(n_2))), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_2_autocor))
}

#Dla realizacje_3

matrix_3_autokow <- apply(realizacje_3, 2, function(x) acf(x, lag.max=h.max_3,
                                                        type="covariance", plot=FALSE)$acf)

matrix_3_autokow <- matrix_3_autokow[-1,]

par(mfrow=c(2,2))
for (h in h.wybrane_3) {
  print(paste0("h=",h))
  tytul_3 <- paste0("Histogram estymatora autokowariancji
              dla h = ",h)
  acf.h_3_autocor <- matrix_3_autokow[h,]
  hist(acf.h_3_autocor, freq=FALSE, col="lightblue", main=tytul_3, xlab="")
  curve(dnorm(x,mean=mean(acf.h_3_autocor), sd=sd(acf.h_3_autocor)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/(5^2*sqrt(n_3))), add=T, col="red", lwd=2)
  legend('topleft', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.6)
  print(summary(acf.h_3_autocor))
}


#########
#Dla realizacje_4
matrix_4_autokow <- apply(realizacje_4, 2, function(x) acf(x, lag.max=h.max_4,
                                                              type="covariance", plot=FALSE)$acf)
par(mfrow=c(2,2))
matrix_4_autokow <- matrix_4_autokow[-1,]
for (h in h.wybrane_4) {
  print(paste0("h=",h))
  tytul_4 <- paste0("Histogram estymatora autokowariancji
              dla h = ",h)
  acf.h_4_autocor <- matrix_4_autokow[h,]
  hist(acf.h_4_autocor, freq=FALSE, col="lightblue", main=tytul_4, xlab="")
  curve(dnorm(x,mean=mean(acf.h_4_autocor), sd=sd(acf.h_4_autocor)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/(5^2*sqrt(n_4))), add=T, col="red", lwd=2)
  legend('topright', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.5)
  print(summary(acf.h_4_autocor))
}

########
#Dla realizacje_5
matrix_5_autokow <- apply(realizacje_5, 2, function(x) acf(x, lag.max=h.max_5,
                                                              type="covariance", plot=FALSE)$acf)
par(mfrow=c(2,2))
matrix_5_autokow <- matrix_5_autokow[-1,]
for (h in h.wybrane_5) {
  print(paste0("h=",h))
  tytul_5 <- paste0("Histogram estymatora autokowariancji
              dla h = ",h)
  acf.h_5_autocor <- matrix_5_autokow[h,]
  hist(acf.h_5_autocor, freq=FALSE, col="lightblue", main=tytul_5, xlab="")
  curve(dnorm(x,mean=mean(acf.h_5_autocor), sd=sd(acf.h_5_autocor)), add=T, col="blue", lwd=2)
  curve(dnorm(x,mean=0, sd=1/(5^2*sqrt(n_5))), add=T, col="red", lwd=2)
  legend('topright', legend=c("Gęstość emp.", "Gęstość teoret."),
         col=c("blue", "red"), lty=c(1,1), cex=0.5)
  print(summary(acf.h_5_autocor))
}

########
# Czwarta kropka + piąta kropka naraz, testujemy wraz z uwzglednieniem losowości
########
###średnie
#dla realizaje_1
ile.powtorz <- 1000

# deklaracja wektorow wynikowych
wynik.sw_srednie_1 <- numeric(ile.powtorz)
wynik.ks_srednie_1 <- numeric(ile.powtorz)

for (i in 1:ile.powtorz) {
  # generujemy k realizacji białego szumu i wyznaczamy średnie
  realizacje_losowe_1 <- matrix(rnorm(n_1*k, mean=0, sd=1),n_1,k)
  srednie_losowe_1    <- apply(realizacje_losowe_1, 2, mean)
  
  # testujemy zgodność z rozkładem normalnym
  wynik.sw_srednie_1[i] <- shapiro.test(srednie_losowe_1)$p.value>0.05
  wynik.ks_srednie_1[i] <- ks.test(srednie_losowe_1,"pnorm", mean=0, sd=1/sqrt(n_1))$p.value>0.05
}


# częstości przyjęcia H0
sum(wynik.ks_srednie_1)/ile.powtorz
sum(wynik.sw_srednie_1)/ile.powtorz

#dla realizacje_2
wynik.sw_srednie_2 <- numeric(ile.powtorz)
wynik.ks_srednie_2 <- numeric(ile.powtorz)

for (i in 1:ile.powtorz) {
  # generujemy k realizacji białego szumu i wyznaczamy średnie
  realizacje_losowe_2 <- matrix(rnorm(n_2*k, mean=0, sd=1),n_2,k)
  srednie_losowe_2    <- apply(realizacje_losowe_2, 2, mean)
  
  # testujemy zgodność z rozkładem normalnym
  wynik.sw_srednie_2[i] <- shapiro.test(srednie_losowe_2)$p.value>0.05
  wynik.ks_srednie_2[i] <- ks.test(srednie_losowe_2,"pnorm", mean=0, sd=1/sqrt(n_2))$p.value>0.05
}


# częstości przyjęcia H0
sum(wynik.ks_srednie_2)/ile.powtorz
sum(wynik.sw_srednie_2)/ile.powtorz

#dla realizacje_3
wynik.sw_srednie_3 <- numeric(ile.powtorz)
wynik.ks_srednie_3 <- numeric(ile.powtorz)

for (i in 1:ile.powtorz) {
  # generujemy k realizacji białego szumu i wyznaczamy średnie
  realizacje_losowe_3 <- matrix(rexp(n_3*k, rate=5),n_3,k)
  realizacje_losowe_3 <- realizacje_losowe_3 - colMeans(realizacje_losowe_3)
  srednie_losowe_3    <- apply(realizacje_losowe_3, 2, mean)
  
  # testujemy zgodność z rozkładem normalnym
  wynik.sw_srednie_3[i] <- shapiro.test(srednie_losowe_3)$p.value>0.05
}


# częstości przyjęcia H0
sum(wynik.ks_srednie_3)/ile.powtorz
sum(wynik.sw_srednie_3)/ile.powtorz

#dla realizacje_4
wynik.sw_srednie_4 <- numeric(ile.powtorz)
wynik.ks_srednie_4 <- numeric(ile.powtorz)

for (i in 1:ile.powtorz) {
  # generujemy k realizacji białego szumu i wyznaczamy średnie
  realizacje_losowe_4 <- matrix(rexp(n_4*k, rate=5),n_4,k)
  realizacje_losowe_4 <- realizacje_losowe_4 - colMeans(realizacje_losowe_4)
  srednie_losowe_4    <- apply(realizacje_losowe_4, 2, mean)
  
  # testujemy zgodność z rozkładem normalnym
  wynik.sw_srednie_4[i] <- shapiro.test(srednie_losowe_4)$p.value>0.05
}


# częstości przyjęcia H0
sum(wynik.ks_srednie_4)/ile.powtorz
sum(wynik.sw_srednie_4)/ile.powtorz


#dla realizacje_5
wynik.sw_srednie_5 <- numeric(ile.powtorz)
wynik.ks_srednie_5 <- numeric(ile.powtorz)

for (i in 1:ile.powtorz) {
  # generujemy k realizacji białego szumu i wyznaczamy średnie
  realizacje_losowe_5 <- matrix(rexp(n_5*k, rate=5),n_5,k)
  realizacje_losowe_5 <- realizacje_losowe_5-colMeans(realizacje_losowe_5)
  srednie_losowe_5    <- apply(realizacje_losowe_5, 2, mean)
  
  # testujemy zgodność z rozkładem normalnym
  wynik.sw_srednie_5[i] <- shapiro.test(srednie_losowe_5)$p.value>0.05
  wynik.ks_srednie_5[i] <- ks.test(srednie_losowe_5,"pnorm", mean=0, sd=1/(5*sqrt(n_5)))$p.value>0.05
}


# częstości przyjęcia H0
sum(wynik.ks_srednie_5)/ile.powtorz
sum(wynik.sw_srednie_5)/ile.powtorz

######
#Teraz testujemy autokorelacje
######

ile.powtorz_1 <- 100

#dla realizacje_1
wynik.sw_autokor_1 <- matrix(rep(0, ile.powtorz_1*h.max_1), nrow=ile.powtorz_1,
                             ncol=h.max_1)
wynik.ks_autokor_1 <- matrix(rep(0, ile.powtorz_1*h.max_1), nrow=ile.powtorz_1,
                             ncol=h.max_1)

for (h in h.wybrane_1) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_1 <- matrix(rnorm(n_1*k),n_1,k)
    matrix_losowe_1 <- apply(realizacje_losowe_1, 2, function(x) acf(x, lag.max=h.max_1,
                                                            type="correlation", plot=FALSE)$acf)
    matrix_losowe_1 <- matrix_losowe_1[-1,]
    h_losowe_1 <- matrix_losowe_1[h,]
    wynik.sw_autokor_1[j,h] <- shapiro.test(h_losowe_1)$p.value>0.05
    wynik.ks_autokor_1[j,h] <- ks.test(h_losowe_1, 'pnorm', mean=0, sd=1/sqrt(n_1))$p.value>0.05
  }
}

#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokor_1[,1])/ile.powtorz_1
sum(wynik.sw_autokor_1[,5])/ile.powtorz_1
sum(wynik.sw_autokor_1[,10])/ile.powtorz_1
sum(wynik.sw_autokor_1[,25])/ile.powtorz_1
sum(wynik.ks_autokor_1[,1])/ile.powtorz_1
sum(wynik.ks_autokor_1[,5])/ile.powtorz_1
sum(wynik.ks_autokor_1[,10])/ile.powtorz_1
sum(wynik.ks_autokor_1[,25])/ile.powtorz_1

#dla realizacje_2


wynik.sw_autokor_2 <- matrix(rep(0, ile.powtorz_1*h.max_2), nrow=ile.powtorz_1,
                             ncol=h.max_2)
wynik.ks_autokor_2 <- matrix(rep(0, ile.powtorz_1*h.max_2), nrow=ile.powtorz_1,
                             ncol=h.max_2)

for (h in h.wybrane_2) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_2 <- matrix(rnorm(n_2*k),n_2,k)
    matrix_losowe_2 <- apply(realizacje_losowe_2, 2, function(x) acf(x, lag.max=h.max_2,
                                                                   type="correlation", plot=FALSE)$acf)
    matrix_losowe_2 <- matrix_losowe_2[-1,]
    h_losowe_2 <- matrix_losowe_2[h,]
    wynik.sw_autokor_2[j,h] <- shapiro.test(h_losowe_2)$p.value>0.05
    wynik.ks_autokor_2[j,h] <- ks.test(h_losowe_2, 'pnorm', mean=0, sd=1/sqrt(n_2))$p.value>0.05
    
  }
}

#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokor_2[,1])/ile.powtorz_1
sum(wynik.sw_autokor_2[,10])/ile.powtorz_1
sum(wynik.sw_autokor_2[,30])/ile.powtorz_1
sum(wynik.sw_autokor_2[,50])/ile.powtorz_1
sum(wynik.ks_autokor_2[,1])/ile.powtorz_1
sum(wynik.ks_autokor_2[,10])/ile.powtorz_1
sum(wynik.ks_autokor_2[,30])/ile.powtorz_1
sum(wynik.ks_autokor_2[,50])/ile.powtorz_1

#dla realizacje_3

wynik.sw_autokor_3 <- matrix(rep(0, ile.powtorz_1*h.max_3), nrow=ile.powtorz_1,
                             ncol=h.max_3)
wynik.ks_autokor_3 <- matrix(rep(0, ile.powtorz_1*h.max_3), nrow=ile.powtorz_1,
                             ncol=h.max_3)

for (h in h.wybrane_3) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_3 <- matrix(rexp(n_3*k, 5),n_3,k)
    realizacje_losowe_3 <- realizacje_losowe_3-colMeans(realizacje_losowe_3)
    acf.matrix_losowe_3 <- apply(realizacje_losowe_3, 2, function(x) acf(x, lag.max=h.max_3,
                                                            type="correlation", plot=FALSE)$acf)
    acf.matrix_losowe_3 <- acf.matrix_losowe_3[-1,]
    acf.h_losowe_3 <- acf.matrix_losowe_3[h,]
    wynik.sw_autokor_3[j,h] <- shapiro.test(acf.h_losowe_3)$p.value>0.05
    wynik.ks_autokor_3[j,h] <- ks.test(acf.h_losowe_3, 'pnorm', mean=0, sd=1/sqrt(n_3))$p.value>0.05
  }
}

#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokor_3[,1])/ile.powtorz_1
sum(wynik.sw_autokor_3[,5])/ile.powtorz_1
sum(wynik.sw_autokor_3[,8])/ile.powtorz_1
sum(wynik.sw_autokor_3[,12])/ile.powtorz_1
sum(wynik.ks_autokor_3[,1])/ile.powtorz_1
sum(wynik.ks_autokor_3[,5])/ile.powtorz_1
sum(wynik.ks_autokor_3[,8])/ile.powtorz_1
sum(wynik.ks_autokor_3[,12])/ile.powtorz_1

#dla realizacje_4
wynik.sw_autokor_4 <- matrix(rep(0, ile.powtorz_1*h.max_4), nrow=ile.powtorz_1,
                             ncol=h.max_4)
wynik.ks_autokor_4 <- matrix(rep(0, ile.powtorz_1*h.max_4), nrow=ile.powtorz_1,
                             ncol=h.max_4)

for (h in h.wybrane_4) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_4 <- matrix(rexp(n_4*k, 5),n_4,k)
    realizacje_losowe_4 <- realizacje_losowe_4-colMeans(realizacje_losowe_4)
    acf.matrix_losowe_4 <- apply(realizacje_losowe_4, 2, function(x) acf(x, lag.max=h.max_4,
                                                                         type="correlation", plot=FALSE)$acf)
    acf.matrix_losowe_4 <- acf.matrix_losowe_4[-1,]
    acf.h_losowe_4 <- acf.matrix_losowe_4[h,]
    wynik.sw_autokor_4[j,h] <- shapiro.test(acf.h_losowe_4)$p.value>0.05
    wynik.ks_autokor_4[j,h] <- ks.test(acf.h_losowe_4, 'pnorm', mean=0, sd=1/sqrt(n_4))$p.value>0.05
    
  }
}

#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokor_4[,1])/ile.powtorz_1
sum(wynik.sw_autokor_4[,13])/ile.powtorz_1
sum(wynik.sw_autokor_4[,28])/ile.powtorz_1
sum(wynik.sw_autokor_4[,34])/ile.powtorz_1
sum(wynik.ks_autokor_4[,1])/ile.powtorz_1
sum(wynik.ks_autokor_4[,13])/ile.powtorz_1
sum(wynik.ks_autokor_4[,28])/ile.powtorz_1
sum(wynik.ks_autokor_4[,34])/ile.powtorz_1

#dla realizacje_5

wynik.sw_autokor_5 <- matrix(rep(0, ile.powtorz_1*h.max_5), nrow=ile.powtorz_1,
                             ncol=h.max_5)
wynik.ks_autokor_5 <- matrix(rep(0, ile.powtorz_1*h.max_5), nrow=ile.powtorz_1,
                             ncol=h.max_5)


for (h in h.wybrane_5) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_5 <- matrix(rexp(n_5*k, 5),n_5,k)
    acf.matrix_losowe_5 <- apply(realizacje_losowe_5, 2, function(x) acf(x, lag.max=h.max_5,
                                                                         type="correlation", plot=FALSE)$acf)
    acf.matrix_losowe_5 <- acf.matrix_losowe_5[-1,]
    acf.h_losowe_5 <- acf.matrix_losowe_5[h,]
    wynik.sw_autokor_5[j,h] <- shapiro.test(acf.h_losowe_5)$p.value>0.05
    wynik.ks_autokor_5[j,h] <- ks.test(acf.h_losowe_5, 'pnorm', mean=0, sd=1/sqrt(n_5))$p.value>0.05
    
  }
}

#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokor_5[,4])/ile.powtorz_1
sum(wynik.sw_autokor_5[,32])/ile.powtorz_1
sum(wynik.sw_autokor_5[,55])/ile.powtorz_1
sum(wynik.sw_autokor_5[,70])/ile.powtorz_1
sum(wynik.ks_autokor_5[,4])/ile.powtorz_1
sum(wynik.ks_autokor_5[,32])/ile.powtorz_1
sum(wynik.ks_autokor_5[,55])/ile.powtorz_1
sum(wynik.ks_autokor_5[,70])/ile.powtorz_1

######
#Testy zgodności dla autokowariancji, nie trzeba ich robić, ponieważ bedzie to wyglądać, jak wyżej
######

#realizacje_3

wynik.sw_autokow_3 <- matrix(rep(0, ile.powtorz_1*h.max_3), nrow=ile.powtorz_1,
                             ncol=h.max_3)

for (h in h.wybrane_3) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_3 <- matrix(rexp(n_3*k, 5),n_3,k)
    realizacje_losowe_3 <- realizacje_losowe_3-colMeans(realizacje_losowe_3)
    acf.matrix_3_autokow <- apply(realizacje_losowe_3, 2, function(x) acf(x, lag.max=h.max_3,
                                              type="covariance", plot=FALSE)$acf)
    acf.matrix_3_autokow <- acf.matrix_3_autokow[-1,]
    acf.h_losowe_3 <- acf.matrix_3_autokow[h,]
    wynik.sw_autokow_3[j,h] <- shapiro.test(acf.h_losowe_3)$p.value>0.05
  }
}

#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokow_3[,1])/ile.powtorz_1
sum(wynik.sw_autokow_3[,5])/ile.powtorz_1
sum(wynik.sw_autokow_3[,8])/ile.powtorz_1
sum(wynik.sw_autokow_3[,12])/ile.powtorz_1

#realizacje_4


wynik.sw_autokow_4 <- matrix(rep(0, ile.powtorz_1*h.max_4), nrow=ile.powtorz_1,
                             ncol=h.max_4)

for (h in h.wybrane_4) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_4 <- matrix(rexp(n_4*k, 5),n_4,k)
    realizacje_losowe_4 <- realizacje_losowe_4-colMeans(realizacje_losowe_4)
    acf.matrix_4_autokow <- apply(realizacje_losowe_4, 2, function(x) acf(x, lag.max=h.max_4,
                                                              type="covariance", plot=FALSE)$acf)
    acf.matrix_4_autokow <- acf.matrix_4_autokow[-1,]
    acf.h_losowe_4 <- acf.matrix_4_autokow[h,]
    wynik.sw_autokow_4[j,h] <- shapiro.test(acf.h_losowe_4)$p.value>0.05
  }
}

#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokow_4[,1])/ile.powtorz_1
sum(wynik.sw_autokow_4[,13])/ile.powtorz_1
sum(wynik.sw_autokow_4[,28])/ile.powtorz_1
sum(wynik.sw_autokow_4[,34])/ile.powtorz_1

#realizacje_5
wynik.sw_autokow_5 <- matrix(rep(0, ile.powtorz_1*h.max_5), nrow=ile.powtorz_1,
                             ncol=h.max_5)

for (h in h.wybrane_5) {
  for (j in 1:ile.powtorz_1) {
    realizacje_losowe_5 <- matrix(rexp(n_5*k, 5),n_5,k)
    realizacje_losowe_5 <- realizacje_losowe_5-colMeans(realizacje_losowe_5)
    acf.matrix_5_autokow <- apply(realizacje_losowe_5, 2, function(x) acf(x, lag.max=h.max_5,
                                                  type="covariance", plot=FALSE)$acf)
    acf.matrix_5_autokow <- acf.matrix_5_autokow[-1,]
    acf.h_losowe_5 <- acf.matrix_5_autokow[h,]
    wynik.sw_autokow_5[j,h] <- shapiro.test(acf.h_losowe_5)$p.value>0.05
  }
}

h.wybrane_5
#Ile razy przyjeliśmy H0 dla poszczególnych opóźnień h?
sum(wynik.sw_autokow_5[,4])/ile.powtorz_1
sum(wynik.sw_autokow_5[,32])/ile.powtorz_1
sum(wynik.sw_autokow_5[,55])/ile.powtorz_1
sum(wynik.sw_autokow_5[,70])/ile.powtorz_1

##################################
#Zadanie 2
##################################


wn_1 <- arima.sim(model=list(order=c(0,0,0)),n=100) #przykładowy biały szum

realizacje_sd1_50 <- matrix(rnorm(50*500, mean=0, sd=1),50, 500)
realizacje_sd2_50 <- matrix(rnorm(50*500, mean=0, sd=2),50, 500)
realizacje_sd05_50 <- matrix(rnorm(50*500, mean=0, sd=0.5),50, 500)

realizacje_sd1_100 <- matrix(rnorm(100*500, mean=0, sd=1),100, 500)
realizacje_sd2_100 <- matrix(rnorm(100*500, mean=0, sd=2),100, 500)
realizacje_sd05_100 <- matrix(rnorm(100*500, mean=0, sd=0.5),100, 500)

realizacje_sd1_1000 <- matrix(rnorm(1000*500, mean=0, sd=1),1000, 500)
realizacje_sd2_1000 <- matrix(rnorm(1000*500, mean=0, sd=2),1000, 500)
realizacje_sd05_1000 <- matrix(rnorm(1000*500, mean=0, sd=0.5),1000, 500)

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

theor_test <- function(ts){
  if(NROW(ts)<=500){
    test <-Box.test(ts, lag=log(NROW(ts)), type='Ljung-Box')
  }
  else{
    if(NROW(ts)>750){
      test <-Box.test(ts,lag=20,type='Ljung-Box')
    }
    else{
      test <-Box.test(ts,lag=10,type='Ljung-Box')
    }
  }
  return(test)
}


graph_results_sd1_50 <- apply(realizacje_sd1_50, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd1_50 <- apply(realizacje_sd1_50,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd1_50 <- apply(realizacje_sd1_50,2,function(x) theor_test(x)$p.value>=0.05)

graph_results_sd1_100 <- apply(realizacje_sd1_100, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd1_100 <- apply(realizacje_sd1_100,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd1_100 <- apply(realizacje_sd1_100,2,function(x) theor_test(x)$p.value>=0.05)

graph_results_sd1_1000 <- apply(realizacje_sd1_1000, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd1_1000 <- apply(realizacje_sd1_1000,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd1_1000 <- apply(realizacje_sd1_1000,2,function(x) theor_test(x)$p.value>=0.05)


graph_results_sd2_50 <- apply(realizacje_sd2_50, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd2_50 <- apply(realizacje_sd2_50,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd2_50 <- apply(realizacje_sd2_50,2,function(x) theor_test(x)$p.value>=0.05)

graph_results_sd2_100 <- apply(realizacje_sd2_100, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd2_100 <- apply(realizacje_sd2_100,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd2_100 <- apply(realizacje_sd2_100,2,function(x) theor_test(x)$p.value>=0.05)

graph_results_sd2_1000 <- apply(realizacje_sd2_1000, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd2_1000 <- apply(realizacje_sd2_1000,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd2_1000 <- apply(realizacje_sd2_1000,2,function(x) theor_test(x)$p.value>=0.05)

graph_results_sd05_50 <- apply(realizacje_sd05_50, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd05_50 <- apply(realizacje_sd05_50,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd05_50 <- apply(realizacje_sd05_50,2,function(x) theor_test(x)$p.value>=0.05)

graph_results_sd05_100 <- apply(realizacje_sd05_100, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd05_100 <- apply(realizacje_sd05_100,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd05_100 <- apply(realizacje_sd05_100,2,function(x) theor_test(x)$p.value>=0.05)

graph_results_sd05_1000 <- apply(realizacje_sd05_1000, 2, function(x) if(graph_test(x,plot=FALSE)[[1]]>=0.95 & graph_test(x,plot=FALSE)[[2]]==FALSE) return(TRUE) else{return(FALSE)})
LB_results_sd05_1000 <- apply(realizacje_sd05_1000,2, function(x) Box.test(x,lag=floor(NROW(x)/4),type='Ljung-Box')$p.value>=0.05)
theor_results_sd05_1000 <- apply(realizacje_sd05_1000,2,function(x) theor_test(x)$p.value>=0.05)

mean(graph_results_sd1_50)
mean(LB_results_sd1_50)
mean(theor_results_sd1_50)

mean(graph_results_sd1_100)
mean(LB_results_sd1_100)
mean(theor_results_sd1_100)

mean(graph_results_sd1_1000)
mean(LB_results_sd1_1000)
mean(theor_results_sd1_1000)

data.frame('T=50'=c(mean(graph_results_sd1_50),
                    mean(LB_results_sd1_50),
                    mean(theor_results_sd1_50)),
           'T=100'=c(mean(graph_results_sd1_100),
                     mean(LB_results_sd1_100),
                     mean(theor_results_sd1_100)),
           'T=1000'=c(mean(graph_results_sd1_1000),
                      mean(LB_results_sd1_1000),
                      mean(theor_results_sd1_1000))
) -> sd_1_wn_results
rownames(sd_1_wn_results) <- c('graphic','LB_test','LB_adjusted_H')
xtable(sd_1_wn_results)

