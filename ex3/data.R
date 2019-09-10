setwd(dirname(rstudioapi::getSourceEditorContext()$path))  #(put data in the same folder)
rm(list = ls()) #clear all other data
graphics.off()  #close existing graphics
library(ggplot2) #load library ggplot2

df <- read.table("data.txt",header=TRUE) #read data
cols <- c("black", "blue", "red", "green", "magenta") #colors
cats <- c("A1", "A2", "A3", "A4", "A5") #names of categories

p1 <- ggplot(df, aes(group=diet, x= cats[diet], y=weight)) +
  geom_boxplot(fill = cols, alpha = 0.3) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = "none")

p1 #boxplot
res.aov <- aov(weight ~ cats[diet], data = df) #anova
summary(res.aov) #display anova outcome
TukeyHSD(res.aov) #perform Tukey test
plot(res.aov, 1) #residuals vs fitted
plot(res.aov, 2) #qq plot
