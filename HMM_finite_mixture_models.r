# This code is implementing a finite mixture model as a particular case of HMMs, as explained in the paper:
# Bayesian Methods for Hidden Markov Models: Recursive Computing in the 21st Century - Author(s): Steven L. Scott

# PARAMETERS
H <- 4 # number of hidden states
alpha_0 <- runif(n = H, min = 1, max = 4)
mu_0 <- 0
tau_0 <- 0.02 # precision
a_0 <- 5
b_0 <- 5

LENGTH <- 100
N <- 100 # number of samples
niter <- 200 # number of Gibbs Sampling iterations


# Defining the unknown mixture
w_real <- gtools::rdirichlet(1, alpha_0)
mu_real <- rnorm(H, mu_0, sqrt(1 / tau_0))
tau_real <- rgamma(H, shape = a_0, rate = b_0)
cat("w_real :", w_real, "\n")
cat("mu_real :", mu_real, "\n")
cat("tau_real :", tau_real, "\n")


# Transition probability matrix 
# For finite mixture models, all rows are equal
Q <- c()
for (h in 1:H) {
  Q <- rbind(Q, w_real)
}
print("\n Transition matrix\n")
print(Q)


# Generating an Hidden Markov Chain
HMC <- function(LENGTH, Q, H) {
  h <- c(1)
  for (i in 1:LENGTH) {
    h <- c(h, sample(1:H, prob = Q[h[length(h)], ], size = 1))
  }
  return(h)
}
h_real <- HMC(LENGTH, Q, H)
print(h_real)

# Generating the gaussians of the mixture
x <- seq(-20, 20, by = 0.001)
norms <- data.frame(x)
col_names <- c("x")
norms <- c()
for (i in 1:H) {
  norms <- cbind(norms, dnorm(x, mu_real[i], sqrt(1 / tau_real[i])) * w_real[i])
}
print(dim(norms))


# Plotting
library(plotly)
z_plot <- c()
y_plot <- c()
x_plot <- c()
for (i in 1:LENGTH) {
  z_plot <- c(z_plot, rep(i, length(x)))
  x_plot <- c(x_plot, x)
  y_plot <- c(y_plot, norms[h_real[i]])

  plot3d(x, norms[h_real[i]], rep(i, length(x)), type = "l", lwd = 4, col = h_real[i], zlim = c(1, LENGTH))
}
#grid3d(c("x", "y+", "z"))