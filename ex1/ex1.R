setwd(dirname(rstudioapi::getSourceEditorContext()$path))  #for R-studio
rm(list = ls())
graphics.off()
library(ggplot2)
library(reliaR)
library(MASS)
library(fitdistrplus)
library(goftest)



#### observe data 
data(dataset2)
df <- data.frame(dataset2)


### plot datapoints
p1 = ggplot(df)+
  geom_point(aes(x = 1:length(dataset2), y = dataset2), shape = 21, color = "black", fill = "black", size = 2, alpha = 0.5)+
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p1

#plot histogram...
p2 = qplot(dataset2, geom="histogram", bins=36, fill= I("black"), col = I("black"), alpha=I(0.15), y = ..density..) + 
  stat_density(geom = "line", linetype = "dashed", color = "blue") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p2

#statistical components
fivenum(dataset2)
mean(dataset2)
sqrt(var(dataset2))


z <- seq(from=0, to = 1100, by = 1)
b1 <- sqrt(6*var(dataset2))/pi
mu1 <- mean(dataset2) - b1*0.5772
mu2 <- median(dataset2) + b1*log(log(2))
g1 <- dgumbel(z, mu1, b1)
#g2 <- dgumbel(z, mu2, b1)

p3 = ggplot()+
  geom_area(aes(x = z, y = dgumbel(z, 150, 50)), color = "red" , fill = "red", alpha = 0.05)+
  geom_area(aes(x = z, y = dgumbel(z, 200, 100)), color = "black", fill = "black", alpha = 0.05)+
  geom_area(aes(x = z, y = dgumbel(z, 500, 150)), color = "blue", fill = "blue", alpha = 0.05)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p3

p4 = qplot(dataset2, geom="histogram", bins=36, fill= I("black"), col = I("black"), alpha=I(0.15), y = ..density..) + 
  geom_area(aes(x = z, y = g1), color = "black", linetype = "dashed", fill = "blue", alpha = 0.05)+
  #geom_area(aes(x = z, y = g2), color = "red", linetype = "dashed", fill = "red", alpha = 0.05)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p4



dgumbel <- function(x, mu, b) 1/b*exp((mu-x)/b)*exp(-exp((mu-x)/b))
pgumbel <- function(q, mu, b) exp(-exp((mu-q)/b))
qgumbel <- function(p, mu, b) mu-b*log(-log(p))
fitgumbel <- fitdist(dataset2, "gumbel", start=list(mu=mu1, b=b1))
summary(fitgumbel)
#plot(fitgumbel) #not visually appealing...
muhat <- fitgumbel$estimate["mu"]
bhat <-  fitgumbel$estimate["b"]


cdf_data <- ecdf(dataset2)
cdf_data <- cdf_data(dataset2)
cdf_g <- pgumbel(dataset2, muhat, bhat)
pdf_g <- dgumbel(z, muhat, bhat)

p5 = qplot(dataset2, geom="histogram", bins=36, fill= I("black"), col = I("black"), alpha=I(0.15), y = ..density..) + 
  geom_area(aes(x = z, y = pdf_g), color = "NA", linetype = "solid", fill = "blue",  alpha = 0.2)+
  geom_line(aes(x = z, y= pdf_g), color = "blue", size = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p5 #histogram best fit


p6 = ggplot() + 
  geom_point(aes(x = cdf_g, y = seq(from = 0, to = 1, length = length(dataset2))), shape = 21, 
             color = "black", fill = "black", size = 3, alpha = 0.5)+ 
  geom_abline(intercept=0, slope = 1, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))

p6 ## pp plot


qis <-  seq(from = 0, to = 1, length = length(dataset2) + 2)
qis <- qis[seq(from=2, to = length(qis) -1)]
qis <- seq(from = 0.5/length(dataset2), to = (length(dataset2)-0.5)/length(dataset2), by = 1/length(dataset2))
q_gumbel <- qgumbel(qis, muhat, bhat)

p7 = ggplot() + 
  geom_point(aes(x = q_gumbel, y = dataset2), shape = 21, color = "black", fill = "black", size = 3, alpha = 0.5)+ #qqplot
  geom_abline(intercept=0, slope = 1, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))

p7 ## qqplot

p8 = ggplot() + 
  geom_point(aes(x = dataset2, y = cdf_data), shape = 21, color = "black", fill = "black", size = 3, alpha = 0.5)+ 
  geom_line(aes(x = z, y = pgumbel(z, muhat, bhat)), color = "blue", size = 1)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))

p8 ## qqplot


gof <- gofstat(fitgumbel)
gof$kstest
ks.test(dataset2, muhat, bhat, alternative = "two.sided")
ad.test(dataset2, pgumbel, muhat, bhat)
cvm.test(dataset2, pgumbel, muhat, bhat)

