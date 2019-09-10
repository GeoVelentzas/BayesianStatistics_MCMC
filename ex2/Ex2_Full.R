setwd(dirname(rstudioapi::getSourceEditorContext()$path))  #for R-studio


# Load libraries ----------------------------------------------------------

rm(list = ls())
graphics.off()
library(ggplot2)

# Computations and Distributions ------------------------------------------

eps <- 1e-3
dtheta <- 0.001
theta <- seq(from = eps, to = 1-eps, by = dtheta)
alpha <- 21
beta <- 5
likelihood <- choose(alpha+beta-2, alpha -1)*theta^(alpha-1)*(1-theta)^(beta-1)
likelihood_pdf <- dbeta(theta, alpha, beta, ncp = 0, log = FALSE)
noninform1 <- dbeta(theta, 1, 1, ncp = 0, log = FALSE)
noninform2 <- dbeta(theta, 1/2, 1, ncp = 0, log = FALSE)
noninform3 <- dbeta(theta, 3/2, 1 , ncp = 0, log = FALSE)


# Plot Likelihood ---------------------------------------------------------

p1 = ggplot()+
  geom_line(aes(x = theta, y = likelihood_pdf), color = "black")+
  geom_area(aes(x = theta, y = likelihood_pdf ),  alpha = 0.1)+
  ylim(0, max(likelihood_pdf)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p1



# Plot First Prior --------------------------------------------------------


p2 = ggplot()+
  geom_line(aes(x = theta, y = likelihood_pdf), color = "black")+
  geom_area(aes(x = theta, y = likelihood_pdf ),  alpha = 0.1)+
  geom_line(aes(x = theta, y = noninform1), color = "black") +
  geom_area(aes(x = theta, y = noninform1), fill = "lightblue",  alpha = 0.3)+
  ylim(0, max(likelihood_pdf)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p2



# Plot second Prior -------------------------------------------------------


p3 = ggplot()+
  geom_line(aes(x = theta, y = likelihood_pdf), color = "black")+
  geom_area(aes(x = theta, y = likelihood_pdf ),  alpha = 0.1)+
  geom_line(aes(x = theta, y = noninform2), color = "black") +
  geom_area(aes(x = theta, y = noninform2), fill = "lightblue",  alpha = 0.3)+
  ylim(0, max(likelihood_pdf)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p3



# Plot third prior --------------------------------------------------------


p3 = ggplot()+
  geom_line(aes(x = theta, y = likelihood_pdf), color = "black")+
  geom_area(aes(x = theta, y = likelihood_pdf ),  alpha = 0.1)+
  geom_line(aes(x = theta, y = noninform3), color = "black") +
  geom_area(aes(x = theta, y = noninform3), fill = "lightblue",  alpha = 0.3)+
  ylim(0, max(likelihood_pdf)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p3









# Plot posteriors + priors ------------------------------------------------


p4 = ggplot()+
  #geom_point(aes(x = z, y = f21), shape = 21, color = "black", fill = "black", size = 4, alpha = 0.5)+
  #geom_point(aes(x = z, y = f22), shape = 21, color = "black",  fill = "blue", size  = 4, alpha = 0.5) +
  #geom_point(aes(x = z, y = f23), shape = 21, color = "black", fill = "red", size = 4, alpha = 0.5) +
  geom_line(aes(x = theta, y = likelihood_pdf), linetype = "dashed", color = "black", size = 0.5)+
  geom_area(aes(x = theta, y = likelihood_pdf), fill = "black",  alpha = 0.05)+
  geom_line(aes(x = theta, y= noninform1), color = "black", size = 0.5, linetype = "dashed") +
  geom_line(aes(x = theta, y= noninform2), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(aes(x = theta, y= noninform3), color = "red", size = 0.5, linetype = "dashed") +
  geom_line(aes(x = theta, y= dbeta(theta, 21,5)), color = "black", size = 0.8) +
  geom_line(aes(x = theta, y= dbeta(theta, 20.5,5)), color = "blue", size = 0.8) +
  geom_line(aes(x = theta, y= dbeta(theta, 21.5,5)), color = "red", size = 0.8) +
  theme_classic() + 
  ylim(0, max(likelihood_pdf)*1.05)+
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p4



# Plot posteriors (zoom in) -----------------------------------------------


p5 = ggplot()+
  geom_area(aes(x = theta, y = likelihood_pdf), fill = "black",  alpha = 0.05)+
  geom_line(aes(x = theta, y= dbeta(theta, 21,5)), color = "black", size = 0.8) +
  geom_line(aes(x = theta, y= dbeta(theta, 20.5,5)), color = "blue", size = 0.8) +
  geom_line(aes(x = theta, y= dbeta(theta, 21.5,5)), color = "red", size = 0.8) +
  theme_classic() + 
  ylim(0, max(likelihood_pdf)*1.05)+
  xlim(0.7, 0.95) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p5




# Compute predictive - Binomial -------------------------------------------


N <- 30
x <- 20
n <- 24
r <- 4
R <- 10
a <- c(1, 0.5, 1.5)
b <- c(1, 1, 1)
z <- seq(from = 0, to = N, by = 1)
f11 <- choose(N, z)*beta(a[1] + x + z, b[1] + n - x + N - z)/ beta(a[1] + x, b[1] + n - x)
f12 <- choose(N, z)*beta(a[2] + x + z, b[2] + n - x + N - z)/ beta(a[2] + x, b[2] + n - x)
f13 <- choose(N, z)*beta(a[3] + x + z, b[3] + n - x + N - z)/ beta(a[3] + x, b[3] + n - x)
#df = data.frame(x_ax = rep(z, 3), vals = c(f11, f12, f13), cols = c(rep("black", N+1), c(rep("blue", N+1)), c(rep("red", N+1))))



# Plot predictive - Binomial ----------------------------------------------


p6 = ggplot()+
  geom_point(aes(x = z, y = f11), shape = 21, color = "black", fill = "black", size = 4, alpha = 0.5)+
  geom_point(aes(x = z, y = f12), shape = 21, color = "black",  fill = "blue", size  = 4, alpha = 0.5) +
  geom_point(aes(x = z, y = f13), shape = 21, color = "black", fill = "red", size = 4, alpha = 0.6) +
  #geom_col(data = df, aes(x = x_ax, y = vals, fill = cols), position = "dodge", color = df$cols)+
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p6



# Compute predictive - Negative Binomial ----------------------------------


z2 <- seq(from = 0, to = 200, by = 1)  
f21 <- choose(R+z2-1, R-1)*beta(a[1] + x + z2, b[1] + r + R)/ beta(a[1] + x, b[1] + r)
f22 <- choose(R+z2-1, R-1)*beta(a[2] + x + z2, b[1] + r + R)/ beta(a[2] + x, b[2] + r)
f23 <- choose(R+z2-1, R-1)*beta(a[3] + x + z2, b[1] + r + R)/ beta(a[3] + x, b[3] + r)



# Plot predictive - Negative Binomial -------------------------------------


p7 = ggplot()+
  geom_point(aes(x = z2, y = f21), shape = 21, color = "black", fill = "black", size = 4, alpha = 0.5)+
  geom_point(aes(x = z2, y = f22), shape = 21, color = "black",  fill = "blue", size  = 4, alpha = 0.5) +
  geom_point(aes(x = z2, y = f23), shape = 21, color = "black", fill = "red", size = 4, alpha = 0.5) +
  #geom_line(aes(x = z, y= f21), color = "black", size = 1) +
  #geom_line(aes(x = z, y= f22), color = "blue", size = 1) +
  #geom_line(aes(x = z, y= f23), color = "red", size = 1) +
  theme_classic() + 
  xlim(0, 60) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p7



# Plot predictive - Negative Binomial (zoom out) --------------------------


p8 = ggplot()+
  geom_line(aes(x = z2, y= f21), color = "black", size = 0.8) +
  geom_line(aes(x = z2, y= f22), color = "blue", size = 0.8) +
  geom_line(aes(x = z2, y= f23), color = "red", size = 0.8) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p8




# Plot predictive - Binomial - MLE ----------------------------------------


theta_hat_1 = (alpha-1+a[1]-1)/(alpha-1 + a[1] -1  +beta -1 + b[1]-1)
pr_1 = choose(N, z) * theta_hat_1^z * (1-theta_hat_1)^(N-z)

p9  = p6 +  geom_point(aes(x = z, y = pr_1), shape = 21, color = "black", fill = "green", size = 4, alpha = 0.7)
p9



# Plot predictive - Negative Binomial - MLE (zoom out) --------------------


pr_2 = choose(R + z2 -1, R-1) * theta_hat_1^z2* (1-theta_hat_1)^R
p10  = p8 + geom_line(aes(x = z2, y= pr_2), color = "green", size = 0.8) 
p10



# Plot predictive - Negative Binomial (zoom in) ---------------------------


p11 = p7 + geom_point(aes(x = z2, y = pr_2), shape = 21, color = "black", fill = "green", size = 4, alpha = 0.5)
p11




# Computations and Bayes Factor -------------------------------------------


a1 <- alpha -1 + a[1]
a2 <- alpha -1 + a[2]
a3 <- alpha -1 + a[3]
b1 <- beta - 1 + b[1]
b2 <- beta - 1 + b[2]
b3 <- beta - 1 + b[3]

m1 <- N*a1 /(a1 + b1)
m2 <- N*a2 /(a2 + b2)
m3 <- N*a3 /(a3 + b3)
mc <- N*theta_hat_1

v1 <- N*a1*b1*(a1 + b1 + N)/( (a1+ b1)^2*(a1 + b1 + 1) )
v2 <- N*a2*b2*(a2 + b2 + N)/( (a2+ b2)^2*(a2 + b2 + 1) )
v3 <- N*a3*b3*(a3 + b3 + N)/( (a3+ b3)^2*(a3 + b3 + 1) )
vc <- N*theta_hat_1*(1-theta_hat_1)

F11 <- cumsum(f11)
F12 <- cumsum(f12)
F13 <- cumsum(f13)
FC  <- cumsum(pr_1)

cf11 <- which(F11>0.05 & F11<0.95)
cf12 <- which(F12>0.05 & F12<0.95)
cf13 <- which(F13>0.05 & F13<0.95)
cfc <-  which(FC>0.05 & FC<0.95)

mnb1 <- R*a1/(b1 - 1)
mnb2 <- R*a2/(b1 - 1)
mnb3 <- R*a3/(b3 - 1)
mnbc <- theta_hat_1*R/(1-theta_hat_1)

vnb1 <- R*(b1+r-1)*a1*(b1+a1-1)/((b1-2)*(b1-1)^2)
vnb2 <- R*(b1+r-1)*a2*(b2+a2-1)/((b2-2)*(b2-1)^2)
vnb3 <- R*(b2+r-1)*a3*(b3+a3-1)/((b3-2)*(b3-1)^2)
vnbc <- theta_hat_1*R/(1-theta_hat_1)^2


F21 <- cumsum(f21)
F22 <- cumsum(f22)
F23 <- cumsum(f23)
FNC  <- cumsum(pr_2)

cf21 <- which(F21>0.05 & F21<0.95)
cf22 <- which(F22>0.05 & F22<0.95)
cf23 <- which(F23>0.05 & F23<0.95)
cfnc <- which(FNC>0.05 & FNC<0.95)


th0 <- 0.65
bf11 <- th0^x*(1-th0)^(n-x) * beta(a[1], b[1]) / beta(a[1] + x, b[1] + n - x)
bf12 <- th0^x*(1-th0)^(n-x) * beta(a[2], b[2]) / beta(a[2] + x, b[2] + n - x)
bf13 <- th0^x*(1-th0)^(n-x) * beta(a[3], b[3]) / beta(a[3] + x, b[3] + n - x)
bf11
bf12
bf13






# Gap ---------------------------------------------------------------------








