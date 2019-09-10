setwd(dirname(rstudioapi::getSourceEditorContext()$path))  #for R-studio

# Load Libraries and Data -------------------------------------------------
rm(list = ls())
graphics.off()
library(ggplot2)
library(reliaR)
library(MASS)
library(fitdistrplus)
library(goftest)
library(coda)
library(KernSmooth)
library(GGally) #for some cool stuff with ggpairs
data(dataset2)




# SCMH function -----------------------------------------------------------


scmh <-  function(x, T, thetas_0, lambda, c,  m0=0, s02=1e6, a0=1, b0=1e-3) {
  n <- length(x)
  theta <-  matrix(0, nrow = T+1, ncol = 2)
  theta[1,] <-  c(thetas_0)
  accept <-matrix(0, nrow = T, ncol = 2) 
  for (t in 2:(T+1)){
    mu_star <- rnorm(1, theta[t-1,1], sqrt(lambda*theta[t-1, 2])) 
    log_L1 <- -sum((x - mu_star)/theta[t-1, 2]) -sum(exp(-(x - mu_star)/theta[t-1, 2])) +
               sum((x - theta[t-1, 1])/theta[t-1, 2]) + sum(exp(-(x - theta[t-1, 1])/theta[t-1, 2]))
    log_F1 <- dnorm(mu_star, m0, sqrt(s02), log=T) - dnorm(theta[t-1, 1], m0, sqrt(s02), log = T)
    log_a1 <- min(0, log_L1 + log_F1)
    
    if (log(runif(1)) < log_a1) {
      theta[t, 1] <- mu_star
      accept[t-1, 1] <- 1 
    }
    else {
      theta[t, 1] <- theta[t-1, 1]
    }
    b_star <- rgamma(1, theta[t-1, 2]^2/(c*theta[t, 1]), rate = theta[t-1, 2]/(c*theta[t, 1]))
    
    log_L2 <- -n*log(b_star) - sum(x - theta[t, 1])/b_star - sum(exp(-(x - theta[t, 1])/b_star)) +
               n*log(theta[t-1, 2]) + sum(x - theta[t, 1])/theta[t-1, 2] + sum(exp(-(x - theta[t, 1])/theta[t-1, 2])) 
    log_F2 <- dgamma(b_star, a0, b0, log = T) - dgamma(theta[t-1, 2], a0, b0, log = T)
    log_G2 <- dgamma(theta[t-1, 2], b_star^2/(c*theta[t, 1]), b_star/(c*theta[t, 1]), log = T) - 
              dgamma(b_star, theta[t-1, 2]^2/(c*theta[t, 1]), theta[t-1, 2]/(c*theta[t, 1]), log = T)
    log_a2 <- min(0, log_L2 + log_F2 + log_G2)
    
    if (log(runif(1)) < log_a2) {
      theta[t, 2] <- b_star
      accept[t-1, 2] <- 1 
    }
    else {
      theta[t, 2] <- theta[t-1, 2]
    }
  }
  return(list(thetas = theta, accepts = accept) )
}




# Hyperparameter Tuning ---------------------------------------------------

T <- 5000
grid_size <- 100
lambda_grid <-  seq(from = 10, to = 50, length = grid_size)
c_grid <- seq(from = 5, to = 20, length = grid_size)
Rates <- matrix(0, nrow=grid_size, ncol=4)
colnames(Rates) <- c("M_vs_lambda", "B_vs_lambda", "M_vs_c", "B_vs_c")
for (i in 1:grid_size){
  out1 <- scmh(dataset2, T, c(212.75, 151.86), lambda_grid[i], 16)
  out2 <- scmh(dataset2, T, c(212.75, 151.86), 37, c_grid[i])
  Rates[i,"M_vs_lambda"] <- mean(out1$accepts[,1])
  Rates[i,"B_vs_lambda"] <- mean(out1$accepts[,2])
  Rates[i,"M_vs_c"] <- mean(out2$accepts[,1])
  Rates[i,"B_vs_c"] <- mean(out2$accepts[,2])
}
df <- data.frame(Rates)


# Fig 1.4.a ---------------------------------------------------------------


ggplot(df, aes(x = lambda_grid, y = df$M_vs_lambda )) + #Mean_acceptvs lambda
  geom_point(size=2, stroke = 1, color="black", alpha = 0.6)+
  geom_smooth(level = 0.9999) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.4.b ---------------------------------------------------------------



ggplot(df, aes(x = lambda_grid, y = df$B_vs_lambda )) + 
  geom_point(size=2, stroke = 1, color="black", alpha = 0.6)+
  geom_smooth(level = 0.9999) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.4.c ---------------------------------------------------------------


ggplot(df, aes(x = c_grid, y = df$M_vs_c )) + 
  geom_point(size=2, stroke = 1, color="black", alpha = 0.6)+
  geom_smooth(level = 0.9999) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.4.d ---------------------------------------------------------------


ggplot(df, aes(x = c_grid, y = df$B_vs_c )) + 
  geom_point(size=2, stroke = 1, color="black", alpha = 0.6)+
  geom_smooth(level = 0.9999) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))




# Two Chains --------------------------------------------------------------

T <- 10^4
out1 <- scmh(dataset2, T, c(300, 300), 37, 16)
out2 <- scmh(dataset2, T, c(100, 100), 37, 16)
erg.mean.mu1 <- cumsum(out1$thetas[,1])/cumsum(rep(1, T+1))
erg.mean.mu2 <- cumsum(out2$thetas[,1])/cumsum(rep(1, T+1))
erg.mean.beta1 <- cumsum(out1$thetas[,2])/cumsum(rep(1, T+1))
erg.mean.beta2 <- cumsum(out2$thetas[,2])/cumsum(rep(1, T+1))
sq.dif.mus <- (erg.mean.mu1 - erg.mean.mu2)^2
sq.dif.betas <- (erg.mean.beta1 - erg.mean.beta2)^2


# Fig 1.5.a ---------------------------------------------------------------


ggplot() + geom_line(aes(x = seq(1,T+1), y = erg.mean.mu1), color = "black", size = 0.6)+
  geom_line(aes(x = seq(1,T+1), y = erg.mean.mu2), color = "blue", size = 0.6) +
  theme_classic() + ylim(190, 230) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.5.b ---------------------------------------------------------------


ggplot() + geom_line(aes(x = seq(1,T+1), y = erg.mean.beta1), color = "black", size = 0.6)+
  geom_line(aes(x = seq(1,T+1), y = erg.mean.beta2), color = "blue", size = 0.6) +
  theme_classic() + ylim(140, 180) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.5.c ---------------------------------------------------------------


ggplot() + geom_line(aes(x = seq(1,T+1), y = sq.dif.mus), color = "black", size = 0.6)+
  geom_line(aes(x = seq(1,T+1), y = sq.dif.betas), color = "darkgray", size = 0.6) +
  theme_classic() + ylim(0, 5) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.5.d ---------------------------------------------------------------


ggplot() +geom_point(aes(x = out1$thetas[,1], y = out1$thetas[,2]), color = "black", size =2, alpha = 0.2) +
  geom_point(aes(x = out2$thetas[,1], y = out2$thetas[,2]), color = "blue", size = 2, alpha = 0.2)+ 
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.6.a ---------------------------------------------------------------


ggplot()+ geom_line(aes(x = seq(1,T+1), y = out1$thetas[, 1]), color = "black", alpha =0.8) +
 geom_line(aes(x = seq(1,T+1), y = out2$thetas[, 1]), color = "blue", alpha = 0.7)+
 theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.6.b ---------------------------------------------------------------


ggplot()+ geom_line(aes(x = seq(1,T+1), y = out1$thetas[, 2]), color = "black", alpha =0.8) +
  geom_line(aes(x = seq(1,T+1), y = out2$thetas[, 2]), color = "blue", alpha = 0.7)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))




# Convert to mcmc list ----------------------------------------------------

mc1 <- mcmc(out1$thetas)
mc2 <- mcmc(out2$thetas)
mc_all <- mcmc.list(mc1, mc2)



# First Case b=10000, m=20000 ---------------------------------------------


burnin <- 10000
T <- 20000 + burnin
out1 <- scmh(dataset2, T, c(212.75, 151.86), 37, 16)
steps <- out1$thetas[seq(burnin+2, T+1), ]
mc <- mcmc(steps)
ac <- autocorr(mc, lags = seq(0, 20))



# Fig 1.8.a ---------------------------------------------------------------


ggplot()+
  geom_point(aes(x=seq(0, 20), y = ac[, 1, 1]), size=4, color = "black", alpha = 0.4)+
  geom_segment(aes(x=seq(0,20), y=ac[, 1, 1], xend=seq(0,20),yend=0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.8.b ---------------------------------------------------------------


ggplot()+
  geom_point(aes(x=seq(0, 20), y = ac[, 2, 2]), size=4, color = "black", alpha = 0.4)+
  geom_segment(aes(x=seq(0,20), y=ac[, 2, 2], xend=seq(0,20),yend=0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))



# Final Case b=10000, m=100000, thin=20 -----------------------------------

burnin <- 10^4 #CHANGE THIS TO 10^3
thin <- 1      #CHANGE THIS TO 20    
size <- 10^4   #CHANGE THIS TO 10^5
T <- thin*size + burnin
out <- scmh(dataset2, T, c(212.75, 151.86), 37, 16) #or 37,16
ind <- seq(from = burnin+2, to = T+1, by = thin)
steps <- out$thetas[ind, ]
accepts <- out$accepts[ind-1, ]
mc <- mcmc(steps)
#my_mc <- mc
#(183.1, 242.5) (133.1, 179.8)
#names(table(x))[table(x)==max(table(x))] ###mode!



# Fig 1.9.b ---------------------------------------------------------------


ggplot() + geom_line(aes(x = seq(1,size), y = steps[, 2]), color = "black", size = 0.2)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))



# Fig 1.9.a ---------------------------------------------------------------


ggplot() + geom_line(aes(x = seq(1,size), y = steps[, 1]), color = "black", size = 0.2)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.9.c ---------------------------------------------------------------


erg.mean.mu <- cumsum(steps[,1])/cumsum(rep(1, size))
erg.mean.beta <- cumsum(steps[,2])/cumsum(rep(1, size))
ggplot() + geom_line(aes(x = seq(1,size), y = erg.mean.mu), color = "black", size = 0.6)+
  theme_classic() + ylim(211, 215) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.9.d ---------------------------------------------------------------


ggplot() + geom_line(aes(x = seq(1,size), y = erg.mean.beta), color = "black", size = 0.6)+
  theme_classic() + ylim(153, 157) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.10.b --------------------------------------------------------------


ggplot()+ 
  geom_histogram(aes(x = steps[, 2], y = ..density..), breaks=seq(115, 200, by=0.9),  col="black", fill="lightblue", alpha = 1) +
  geom_density(aes(x = steps[, 2]), col = "black", linetype = "solid", size = 0.4,  fill = "NA", alpha = .1) +
  geom_point(aes(x=c(133.1, 154.6, 179.8), y = c(0,0,0) ), color = c("blue", "black", "blue"), size=4)+
  xlim(115, 200)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



# Fig 1.10.a --------------------------------------------------------------


ggplot()+ 
  geom_histogram(aes(x = steps[, 1], y = ..density..), breaks=seq(155, 275, by=1.1),  col="black", fill="red", alpha = .1) +
  geom_density(aes(x = steps[,1]), col = "black", size=.4, fill = "NA", alpha = .1) +
  geom_point(aes(x=c(182.7, 212.8, 243.7), y = c(0,0,0) ), color = c("red", "black", "red"), size=4)+
  xlim(155, 275)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14), legend.position = "none",
  axis.text.y = element_text(size = 14),
  panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.11.a --------------------------------------------------------------


ggplot() +geom_point(aes(x = steps[,1], y = steps[,2]), color = "black", size =2, alpha = 0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Fig 1.11.b --------------------------------------------------------------


z <- bkde2D(steps[, 1:2], 2)
persp(z$fhat, main="none",
      zlab = "p(μ,β|y)",
      theta = 65, phi = 19,
      col = "lightblue", shade = 0.3, xlab = "μ", ylab="β", border = "NA")









# Others ------------------------------------------------------------------



MUS <- out1$thetas[, 1]
BS <- out1$thetas[, 2]
erg.mean.mu <- cumsum(MUS)/cumsum(rep(1, length(MUS)))
erg.mean.b <- cumsum(BS)/cumsum(rep(1, length(BS)))
p1 <- ggplot() +
  geom_line(aes(x= 1: length(BS), y = erg.mean.mu)) # ergodig mean mu
p1
p2 <- ggplot() +
  geom_line(aes(x = 1:length(MUS), y = erg.mean.b)) #ergodig mean sigma
p2
p3 <- ggplot()+
  geom_line(aes(x = 1: length(out1$thetas[, 1]), y = out1$thetas[, 1]))
p3
p4 <- ggplot()+
  geom_line(aes(x = 1: length(out1$thetas[, 2]), y = out1$thetas[, 2] ) )
p4
p5 <- ggplot() + 
  geom_point(aes(x = out1$thetas[,1], y = out1$thetas[,2]))
p5

















