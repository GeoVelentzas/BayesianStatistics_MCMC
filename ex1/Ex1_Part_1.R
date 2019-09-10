setwd(dirname(rstudioapi::getSourceEditorContext()$path))  #for R-studio
# Load Libraries ----------------------------------------------------------
rm(list = ls())
graphics.off()
library(ggplot2)
library(reliaR)
library(MASS)
library(fitdistrplus)
library(goftest)


# Load Data ---------------------------------------------------------------
#read data
data(dataset2)
#convert data to data frame
df <- data.frame(dataset2)



# Fig 1.1.a ---------------------------------------------------------------
ggplot(df)+
  geom_point(aes(x = 1:length(dataset2), y = dataset2), shape = 21, color = "black", fill = "black", size = 2, alpha = 0.5)+
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.1.b ---------------------------------------------------------------

qplot(dataset2, geom="histogram", bins=36, fill= I("black"), col = I("black"), alpha=I(0.15), y = ..density..) +  #historgram
  stat_density(geom = "line", linetype = "dashed", color = "blue") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Table 1.1 ---------------------------------------------------------------
fivenum(dataset2) #five number statistics
mean(dataset2)   #data sample mean
sqrt(var(dataset2)) #data sample std



# Calculation for a first approach of Gumbel ------------------------------
z <- seq(from=0, to = 1100, by = 1)      #just index
b1 <- sqrt(6*var(dataset2))/pi           #solve Var(x)=pi^2*beta^2/6 for beta
mu1 <- mean(dataset2) - b1*0.5772        #use it to find mu from E(X)=mu+beta*gamma
mu2 <- median(dataset2) + b1*log(log(2)) #use beta to find the median also
g1 <- dgumbel(z, mu1, b1)                #gumbel distribution with known mean and variance


# Fig 1.2.a ---------------------------------------------------------------

ggplot()+
  geom_area(aes(x = z, y = dgumbel(z, 150, 50)), color = "red" , fill = "red", alpha = 0.05)+
  geom_area(aes(x = z, y = dgumbel(z, 200, 100)), color = "black", fill = "black", alpha = 0.05)+
  geom_area(aes(x = z, y = dgumbel(z, 500, 150)), color = "blue", fill = "blue", alpha = 0.05)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))



# Fig 1.2.b ---------------------------------------------------------------

qplot(dataset2, geom="histogram", bins=36, fill= I("black"), col = I("black"), alpha=I(0.15), y = ..density..) + 
  geom_area(aes(x = z, y = g1), color = "black", linetype = "dashed", fill = "blue", alpha = 0.05)+
  #geom_area(aes(x = z, y = g2), color = "red", linetype = "dashed", fill = "red", alpha = 0.05)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Gumbel functions to use -------------------------------------------------
dgumbel <- function(x, mu, b) 1/b*exp((mu-x)/b)*exp(-exp((mu-x)/b))
dgumbel_log <- function(x, mu, b) -log(b) + (mu-x)/b - exp((mu-x)/b)
pgumbel <- function(q, mu, b) exp(-exp((mu-q)/b))
qgumbel <- function(p, mu, b) mu-b*log(-log(p))


# MLE Gumbel Fit ----------------------------------------------------------

fitgumbel <- fitdist(dataset2, "gumbel", start=list(mu=mu1, b=b1))
summary(fitgumbel)
muhat <- fitgumbel$estimate["mu"]
bhat <-  fitgumbel$estimate["b"]

cdf_data <- ecdf(dataset2)
cdf_data <- cdf_data(dataset2)
cdf_g <- pgumbel(dataset2, muhat, bhat)
pdf_g <- dgumbel(z, muhat, bhat)
qis <-  seq(from = 0, to = 1, length = length(dataset2) + 2)
qis <- qis[seq(from=2, to = length(qis) -1)]
qis <- seq(from = 0.5/length(dataset2), to = (length(dataset2)-0.5)/length(dataset2), by = 1/length(dataset2))
q_gumbel <- qgumbel(qis, muhat, bhat)

pdf_bayes <- dgumbel(z, 212.8, 154.6) #this was computed after the Bayesian approach



# Fig 1.3.a (can put Bayesian on top) -------------------------------------

p5 <- qplot(dataset2, geom="histogram", bins=36, fill= I("black"), col = I("black"), alpha=I(0.15), y = ..density..) + 
  geom_area(aes(x = z, y = pdf_g), color = "NA", linetype = "solid", fill = "blue",  alpha = 0.2)+
  geom_line(aes(x = z, y= pdf_g), color = "blue", size = 1)+
  #geom_line(aes(x = z, y= pdf_bayes), color = "red", size = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p5


# Fig 1.3.b (q-q plot) ----------------------------------------------------

p6 <-  ggplot() + 
  geom_point(aes(x = cdf_g, y = seq(from = 0, to = 1, length = length(dataset2))), shape = 21, 
             color = "black", fill = "black", size = 3, alpha = 0.5)+ 
  geom_abline(intercept=0, slope = 1, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p6


# Fig 1.3.c (p-p plot) ----------------------------------------------------

p7 <- ggplot() + 
  geom_point(aes(x = q_gumbel, y = dataset2), shape = 21, color = "black", fill = "black", size = 3, alpha = 0.5)+ #qqplot
  geom_abline(intercept=0, slope = 1, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p7



# Fig 1.3.d (cdf plot) ----------------------------------------------------

p8 <- ggplot() + 
  geom_point(aes(x = dataset2, y = cdf_data), shape = 21, color = "black", fill = "black", size = 3, alpha = 0.5)+ 
  geom_line(aes(x = z, y = pgumbel(z, muhat, bhat)), color = "blue", size = 1)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p8


# Goodness of fit tests ---------------------------------------------------

gofstat(fitgumbel)
ks.test(dataset2, muhat, bhat, alternative = "two.sided")
ad.test(dataset2, pgumbel, muhat, bhat)
cvm.test(dataset2, pgumbel, muhat, bhat)



# 95% confints for gumbel parameters --------------------------------------

summary(fitgumbel)
a <- confint(fitgumbel, level = 0.95)
z1 = seq(from = a[1,1], to = a[1,2], by = 0.001)
z2 = seq(from = a[2,1], to = a[2,2], by = 0.001)
dnorm(a[1,1], 0, 10^3) - dnorm(a[1,2], 0, 10^3)
dgamma(a[2,1],1, 10^(-3)) - dgamma(a[2,2], 1, 10^(-3))









# Gap space ---------------------------------------------------------------




















