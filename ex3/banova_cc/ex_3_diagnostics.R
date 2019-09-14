setwd(dirname(rstudioapi::getSourceEditorContext()$path))  #put data in the same folder
rm(list = ls()) #clear all variables
graphics.off() #close all graphics
library(coda) #load coda
library(lattice) #load lattice
library(psych) #load psych (for scatter matrix)
library(ggplot2) #for nice plots

codamenu() #Diagnostics with CODA

################ CORNER CONSTRAINTS ##########################

df1 <- read.table("data_cc_1.txt") #load data
ind1 <- read.table("indexes_cc.txt") #load indices
m1 <- matrix(data = df1[,2], 100000, 6) #matrix with parameter sample in columns
mf1 <- as.data.frame(m1) #data frame for plots

ggplot(mf1) + #colored histograms with nice visualizations
  geom_histogram(aes(x = V1, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="red", alpha = .3) +
  geom_histogram(aes(x = V2, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="green", alpha = .3) +
  geom_histogram(aes(x = V3, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="blue", alpha = .3) +
  geom_histogram(aes(x = V4, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="black", alpha = .3) +
  geom_histogram(aes(x = V5, y = ..density..), breaks=seq(17, 37, by=0.1),  col="NA", fill="magenta", alpha = .3) +
  geom_density(aes(x = V1), col = "black", fill = "NA", alpha = .0, size = 0.2) +
  geom_density(aes(x = V2), col = "black", fill = "NA", alpha = .0, size = 0.2) +
  geom_density(aes(x = V3), col = "black", fill = "NA", alpha = .0, size = 0.2) +
  geom_density(aes(x = V4), col = "black", fill = "NA", alpha = .0, size = 0.2) +
  geom_density(aes(x = V5), col = "black", fill = "NA", alpha = .0, size = 0.2) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggplot(mf1)+ #densities estimated from histograms
  geom_density(aes(x = V1), col = "red", fill = "red", alpha = .1) +
  geom_density(aes(x = V2), col = "green", fill = "green", alpha = .1) +
  geom_density(aes(x = V3), col = "blue", fill = "blue", alpha = .1) +
  geom_density(aes(x = V4), col = "black", fill = "black", alpha = .1) +
  geom_density(aes(x = V5), col = "magenta", fill = "magenta", alpha = .1) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggplot(mf1)+ #density for variance parameter
  geom_histogram(aes(x = V6, y = ..density..), breaks=seq(8, 30, by=0.3),  col="NA", fill="black", alpha = .3) +
  geom_density(aes(x = V6), col = "black", fill = "NA", alpha = .0, size = 0.1) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


thin <- seq(from = 90000, to = 100000, by = 100) #possible thining for scatter-matrix
pairs.panels(mf1[thin, ], #scatter matrix with correlations and histograms
             method = "pearson", # correlation method
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

mean(mf1$V1>mf1$V2) - mean(mf1$V1<mf1$V2) #P(mu1>mu2) - P(mu1<mu2)
mean(mf1$V1>mf1$V3) - mean(mf1$V1<mf1$V3) #equivalent chekcs
mean(mf1$V1>mf1$V4) - mean(mf1$V1<mf1$V4)
mean(mf1$V1>mf1$V5) - mean(mf1$V1<mf1$V5)
mean(mf1$V2>mf1$V3) - mean(mf1$V2<mf1$V3)
mean(mf1$V2>mf1$V4) - mean(mf1$V2<mf1$V4)
mean(mf1$V2>mf1$V5) - mean(mf1$V2<mf1$V5)
mean(mf1$V3>mf1$V4) - mean(mf1$V3<mf1$V4)
mean(mf1$V3>mf1$V5) - mean(mf1$V3<mf1$V5)
mean(mf1$V4>mf1$V5) - mean(mf1$V4<mf1$V5)



######################## ZERO SUM CONSTRAINTS ##############################

df2 <- read.table("data_zsc_1.txt")
ind2 <- read.table("indexes_zsc.txt")
m2 <- matrix(data = df2[,2], 100000, 6)
mf2 <- as.data.frame(m2)

ggplot(mf2) + 
  geom_histogram(aes(x = V1, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="red", alpha = .3) +
  geom_histogram(aes(x = V2, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="green", alpha = .3) +
  geom_histogram(aes(x = V3, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="blue", alpha = .3) +
  geom_histogram(aes(x = V4, y = ..density..), breaks=seq(20, 37, by=0.1),  col="NA", fill="black", alpha = .3) +
  geom_histogram(aes(x = V5, y = ..density..), breaks=seq(17, 37, by=0.1),  col="NA", fill="magenta", alpha = .3) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggplot(mf2)+
  geom_density(aes(x = V1), col = "red", fill = "red", alpha = .1) +
  geom_density(aes(x = V2), col = "green", fill = "green", alpha = .1) +
  geom_density(aes(x = V3), col = "blue", fill = "blue", alpha = .1) +
  geom_density(aes(x = V4), col = "black", fill = "black", alpha = .1) +
  geom_density(aes(x = V5), col = "magenta", fill = "magenta", alpha = .1) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggplot(mf2)+
  geom_density(aes(x = V6), col = "black", fill = "black", alpha = .1) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


thin <- seq(from = 1, to = 100000, by = 100)
pairs(mf2[thin, 1:5], pch = 19)
pairs.panels(mf2[thin, ], 
             method = "pearson", # correlation method
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


mean(mf2$V1>mf2$V2) - mean(mf2$V1<mf2$V2) 
mean(mf2$V1>mf2$V3) - mean(mf2$V1<mf2$V3) 
mean(mf2$V1>mf2$V4) - mean(mf2$V1<mf2$V4)
mean(mf2$V1>mf2$V5) - mean(mf2$V1<mf2$V5)
mean(mf2$V2>mf2$V3) - mean(mf2$V2<mf2$V3)
mean(mf2$V2>mf2$V4) - mean(mf2$V2<mf2$V4)
mean(mf2$V2>mf2$V5) - mean(mf2$V2<mf2$V5)
mean(mf2$V3>mf2$V4) - mean(mf2$V3<mf2$V4)
mean(mf2$V3>mf2$V5) - mean(mf2$V3<mf2$V5)
mean(mf2$V4>mf2$V5) - mean(mf2$V4<mf2$V5)











