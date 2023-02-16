################################################################################
# Topias Tolonen 16.2.2023
# Uppsala University
# Plots for article regarding Hiring and Firing: a Signaling Game
# Maths based on our article together with Prof. Erik Ekström

################################################################################
# Packages

#install.packages("rootSolve")
#install.packages("magrittr")
#install.packages("ggplot2", dependencies=T)
#install.packages("ggpubr")
library(rootSolve) # for solving \gamma
library(magrittr) # for using lean pipes
library(ggplot2) # for plotting
library(ggpubr) # for text_grob()

################################################################################
# Parameters

# Set parameters, we want c0<mu0<c1<mu1 in our model
# For applied purposes, we wish that c0 and c1 would not bee "too far apart"
# These are chosen fairly arbitrarily for plot purposes
# c0 = low salary claim
# c1 = high salary claim
# mu0 = low type productivity
# mu1 = high type productivity
c0 <- 1.2
mu0 <- 1.4
c1 <- 1.5
mu1 <- 1.7

# Second set of parameters
# r = interest rate
# sigma = volatility of the underlying Brownian Motion
# omega = signal-to-noise-ratio
r <- 0.05
sigma <- 1
omega <- (mu1-mu0)/sigma

################################################################################
# Solving for further constants

# Solve for \gamma^1 in interval (-10,0) and \gamma^2 for similar pos. interval
gamma_eqn <- function(gamma) gamma^2 - gamma - (2*r)/omega^2
gamma1 <- gamma_eqn %>% rootSolve::uniroot.all(c(-10,0))
gamma2 <- gamma_eqn %>% rootSolve::uniroot.all(c(0,10))

# Calculate the boundary level b
b <- (-(c1-mu0)*gamma1)/(mu1-c1-(mu1-mu0)*gamma1)

# Calculate unique \hat{p} for which U(\hat{p}) = c0/r
phat <- (b*(1-(c0/c1))^(1/gamma1))/(1-b+b*(1-(c0/c1))^(1/gamma1))

# Calculate the multiplier A1 spanning from V(pi)'s ODE solution 
A1 <- -((mu1-mu0)*b)/(r*(gamma1-b)*(b/(1-b))^(gamma1))

################################################################################
# Establishing the Value functions U and V

# Function for U with boundary condition U=0 when pi<b.
# Note that U(pi) here represents the Low type employee (High type doesn't mix).
U_eqn <- function(pi) { 
  output <- (c1/r)*(1-((pi*(1-b))/((1-pi)*b))^gamma1)
  output[pi < b] <- 0
  return(output)
}

# Function for V when observing C=c1.
V_eqn <- function(pi) {
  output <- A1*(1-pi)*(pi/(1-pi))^gamma1 + (mu0-c1+(mu1-mu0)*pi)/(r)
  output[pi < b | output<0] <- 0
  return(output)
}

################################################################################
# Training plots

# Grid for the training plot, this is used as dummy in final plots as well.
pi_grid <- seq(0,1,by=0.005)

# Training plots 
# You can uncomment the following if you want to tweak the parameters for quick overview

# Set the plot window to contain 2 pictures
par(mfrow = c(1,2))

# Plot U
pi_grid %>% U_eqn() %>% plot(pi_grid,.,type="l", main="Value Function for Player 1", xlab="pi", ylab="U", col="black")
abline(h=c0/r, col="blue") # low perpetuity salary
abline(h=c1/r, col="magenta") # max salary
clip(x1=0,x2=1,y1=0,y2=U_eqn(phat)) # clip the subsequent abline
abline(v=phat,lty = "dashed", col="red") # Mark \hat{p}
do.call("clip", as.list(par("usr")))  # reset to plot region
abline(v=b, col="red") # draw the boundary condition

#plot V
pi_grid %>% V_eqn() %>% plot(pi_grid,.,type="l", main="Value Function for Player 2", xlab="pi", ylab="V", col="black")
abline(v=b, col="red") # draw the boundary condition
abline(h=(mu1-c1)/r, col="magenta") # draw max salary

#set the plot window back to default
par(mfrow = c(1,1))

################################################################################
# Final plots with ggplot2()
# Note that these are tuned for the current parameter values
# Some of the annotations and textes might not appear in nice positions if parameter values are changed.

# Plotting U
base <-  ggplot() # Base of the plot
base +
  # Etablish the value function and name it
  geom_function(fun = U_eqn, n=100000, linewidth=1) +
  annotation_custom(text_grob(expression(U(pi)), col = "black"), xmin = 0.85, xmax = 0.85, ymin = U_eqn(0.85) + (c1/r)/20, ymax = U_eqn(0.85) + (c1/r)/20) +

  # Scale limits, breaks, extensions
  scale_x_continuous(limits = c(0,1), breaks = seq(from = 0, to = 1, by = 0.2), expand = expansion(mult = c(0,0))) +
  scale_y_continuous(limits = c(0,c1/r+(c1/r)/30), breaks = NULL, expand = expansion(mult = c(0,0))) +

  # Horisontal line for c0/r
  geom_segment(aes(x = -Inf, xend = phat, y = c0/r, yend = c0/r), linewidth = 0.75, linetype = "twodash", col = "black") +
  annotation_custom(text_grob(expression(c[0]/r), col = "black"), xmin = 0.04, xmax = 0.04, ymin = c0/r-(c1/r)/20, ymax = c0/r-(c1/r)/20) +

  # Horisontal line for c1/r
  geom_segment(aes(x = -Inf, xend = 1, y = c1/r, yend = c1/r), linewidth = 0.75, linetype = "dotdash", col = "black") +
  annotation_custom(text_grob(expression(c[1]/r), col = "black"), xmin = 0.04, xmax = 0.04, ymin = c1/r-(c1/r)/20, ymax = c1/r-(c1/r)/20) +

  # Vertical line for p^{hat}
  geom_segment(aes(x=phat,xend=phat,y=0, yend=U_eqn(phat)),linewidth = 0.75, linetype = "dashed", col = "red") +
  annotation_custom(text_grob(expression(hat(p)), col = "red"), xmin = phat, xmax = phat, ymin = -(c1/r)/30, ymax = -(c1/r)/30) +

  # Vertical line for pi=b
  geom_segment(aes(x=b,xend=b,y=0, yend=+Inf),linewidth = 0.75, linetype = "dotted", col = "red") +
  annotation_custom(text_grob(expression(b), col = "red"), xmin = b, xmax = b, ymin = -(c1/r)/30, ymax = -(c1/r)/30) +

  # Add plot grid, quite rigorously but still
  geom_segment(aes(x = -Inf, xend = Inf, y = c1/r+(c1/r)/30, yend = c1/r+(c1/r)/30), linewidth = 0.5, col = "black") +
  geom_segment(aes(x=1,xend=1,y=-Inf, yend=+Inf),linewidth = 0.5, col = "black") +

  # Disable text clipping outside plot area
  coord_cartesian(clip = "off") +

  # Axis names
  labs(x = expression(pi), y="") +

  # Theme for other aesthetics
  theme_classic(base_size = 12)

# Plotting V
base <-  ggplot() # Base of the plot
base +
  # Etablish the value function and name it
  geom_function(fun = V_eqn, n=10000, linewidth=1) +
  annotation_custom(text_grob(expression(V(pi)), col = "black"), xmin = 0.85, xmax = 0.85, ymin = V_eqn(0.85)+((mu1-c1)/r)/15, ymax = V_eqn(0.85)+((mu1-c1)/r)/15) +

  # Scale limits, breaks, extensions
  scale_x_continuous(limits = c(0,1), breaks = seq(from = 0, to = 1, by = 0.2), expand = expansion(mult = c(0,0))) +
  scale_y_continuous(limits = c(0,max(V_eqn(pi_grid))+((mu1-c1)/r)/30), breaks = NULL, expand = expansion(mult = c(0,0))) +

  # Horisontal line for (mu1-c1)/r
  geom_segment(aes(x = -Inf, xend = 1, y = (mu1-c1)/r, yend = (mu1-c1)/r), linewidth = 0.75, linetype = "dotdash", col = "black") +
  annotation_custom(text_grob(expression((mu[1]-c[1])/r), col = "black"), xmin = 0.075, xmax = 0.075, ymin = (mu1-c1)/r-((mu1-c1)/r)/15, ymax = (mu1-c1)/r-((mu1-c1)/r)/15) +

  # Vertical line for pi=b
  geom_segment(aes(x=b,xend=b,y=0, yend=+Inf),linewidth = 0.75, linetype = "dotted", col = "red") +
  annotation_custom(text_grob(expression(b), col = "red"), xmin = b, xmax = b, ymin = -((mu1-c1)/r)/30, ymax = -((mu1-c1)/r)/30) +

  # Add plot grid, quite rigorously but still
  geom_segment(aes(x = -Inf, xend = Inf, y = (mu1-c1)/r+((mu1-c1)/r)/30, yend = (mu1-c1)/r + ((mu1-c1)/r)/30), linewidth = 0.5, col = "black") +
  geom_segment(aes(x=1,xend=1,y=-Inf, yend=+Inf),linewidth = 0.5, col = "black") +

  # Disable text clipping outside plot area
  coord_cartesian(clip = "off") +

  # Axis names
  labs(x = expression(pi), y="") +
  
  # Theme for other aestethics
  theme_classic(base_size = 12)

########################
###### THANK YOU #######
########################

################################################################################
################################################################################
################################################################################
# Unused parameters and other calculations close to our paper, but not used.

# Solve for \eta^1 in interval (-10,0)
# We don't use eta in our example cases
#eta_eqn <- function(eta) eta^2 + eta - (2*r)/sigma^2
#eta1 <- eta_eqn %>% rootSolve::uniroot.all(c(-10,0))

# Boundary Levels
# Could be used to calculate U and V for other combinations of mu, C, pi0, pi1
# Not necessarily correct 
#b0 <- (-(c0-mu0)*gamma1)/(mu1-c0-(mu1-mu0)*gamma1)
#b1 <- (-(c1-mu0)*gamma1)/(mu1-c1-(mu1-mu0)*gamma1)

################################################################################
# Comparison of Pi0, Pi1
# Used to justify further equations t.ex. in region Pi0<pi<Pi1...

# #pies solved in interval (-10,10)
# #pi0
# pi0_eqn <- function(pi0) (c0/r)*(1-((pi0*(1-b0))/((1-pi0)*b0)*gamma1))-(c1/r)*(1-((pi0*(1-b1))/((1-pi0)*b1)*gamma1))
# pi0 <- pi0_eqn %>% rootSolve::uniroot.all(c(-10,10))
# 
# #pi1
# pi1_eqn <- function(pi1) (c0/r)*(1-((pi1*(1-b0))/((1-pi1)*b0)*eta1))-(c1/r)*(1-((pi1*(1-b1))/((1-pi1)*b1)*eta1))
# pi1 <- pi1_eqn %>% rootSolve::uniroot.all(c(-10,10))
# 
# #print values and discussion
# print(pi1)
# print(pi0)
# print("pi1 > pi0:")
# (pi1 > pi0)

#if there are no mistakes it would seem indeed that \pi^1 < \pi^0


################################################################################
# establishing alternative constants not used 
#A1a<- ((mu1-mu0)*b0)/(r*(gamma1-b0*(b0/(1-b0))^(gamma1)))
#A1b<- ((mu1-mu0)*b1)/(r*(gamma1-b1*(b0/(1-b1))^(gamma1)))
#A2a <- (-A1a*(1-b0)*(b0/(1-b0))^gamma1-(c1-mu0-(mu1-mu0)*b0)/r)/((1-b0)*(b0/(1-b0))^gamma2)
#A2b <- (-A1b*(1-b1)*(b1/(1-b1))^gamma1-(c1-mu0-(mu1-mu0)*b1)/r)/((1-b1)*(b1/(1-b1))^gamma2)
#A11 <- (((mu1-mu0)*b1)/(r*(gamma1-b1)))*(b1/(1-b1))^(-gamma1)
#A10 <- (((mu1-mu0)*b0)/(r*(gamma1-b0)))*(b0/(1-b0))^(-gamma1)

# This could be used but converges to 0 so it vanishes 
#A2 <- (-A1*(1-b)*(b/(1-b))^gamma1-(c1-mu0-(mu1-mu0)*b)/r)/((1-b)*(b/(1-b))^gamma2)

################################################################################
##establish test U's and V's equations (a bit of repetition) for odd parameters
#Vtest1_eqn <- function(pi) A1*(1-pi)*(pi/(1-pi))^gamma1 + A2b*(1-pi)*(pi/(1-pi))^gamma2 + (c1-mu0-(mu1-mu0)*pi)/(r)
#V1_eqn <- function(pi) A1a*(1-pi)*(pi/(1-pi))^gamma1 + A2a*(1-pi)*(pi/(1-pi))^gamma2 + (c1-mu0-(mu1-mu0)*pi)/(r) 
#V2_eqn <- function(pi) A1b*(1-pi)*(pi/(1-pi))^gamma1 + A2b*(1-pi)*(pi/(1-pi))^gamma2 + (c1-mu0-(mu1-mu0)*pi)/(r) 

################################################################################
## Alternative formulations for different parameters, not used in paper
# U1_eqn <- function(pi) (c1/r)*(1-((pi*(1-b0))/((1-pi)*b0))^gamma1)
# U2_eqn <- function(pi) (c1/r)*(1-((pi*(1-b1))/((1-pi)*b1))^gamma1)
# 
# #for weird parameters
# U00_eqn <- function(pi) (c0/r)*(1-((pi*(1-b0))/((1-pi)*b0)^gamma1))
# U10_eqn <- function(pi) (c1/r)*(1-((pi*(1-b1))/((1-pi)*b1)^gamma1))
# 
# U01_eqn <- function(pi) (c0/r)*(1-((pi*(1-b0))/((1-pi)*b0)^eta1))
# U11_eqn <- function(pi) (c1/r)*(1-((pi*(1-b1))/((1-pi)*b1)^eta1))
# 
# V0_eqn <- function(pi) A10*(1-pi)*(pi/(1-pi))^gamma1 + (c0-mu0-(mu1-mu0)*pi)/(r) 
# V1_eqn <- function(pi) A11*(1-pi)*(pi/(1-pi))^gamma1 + (c0-mu0-(mu1-mu0)*pi)/(r) 

################################################################################
## Alternative plot sketches for weird parameters
#pi_grid %>% U01_eqn() %>% plot(pi_grid,.,type="line", main="U^{i,j}", xlab="pi", ylab="U", col="magenta")
#pi_grid %>% U11_eqn() %>% lines(pi_grid,.,type="line",col="black")
#pi_grid %>% U00_eqn() %>% lines(pi_grid,.,type="line",col="brown")
#pi_grid %>% U10_eqn() %>% lines(pi_grid,.,type="line",col="red")

#legend
#legend(0,1000,legend=c("U01","U11","U00","U10"), col=c("magenta","black","brown","red"),lty=rep(1,4), cex=0.5)

#pi_grid %>% V0_eqn() %>% plot(pi_grid,.,type="line", main="V^{i}", xlab="pi", ylab="V", col="magenta")
#pi_grid %>% V1_eqn() %>% lines(pi_grid,.,type="line",col="black")
#legend(0,-0.3,legend=c("V0","V1"), col=c("magenta","black"),lty=rep(1,4), cex=0.5)

#par(mfrow = c(1,1))
