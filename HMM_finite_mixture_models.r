# This code is combining the following two papers
# Bayesian Methods for Hidden Markov Models: Recursive Computing in the 21st Century - Author(s): Steven L. Scott
# Gibbs Sampling for Bayesian Finite Mixtures of Normal Distributions - Cristina Mollica, Luca Tardella


# PARAMETERS
H <- 2 # number of hidden states
alpha_0 <- runif(n = H, min = 1, max = 4)
mu_0 <- 0
tau_0 <- 0.02 # precision
a_0 <- 5
b_0 <- 5

LENGTH <- 100 # length of the chain
niter <- 100 # number of Gibbs Sampling iterations


# Defining the unknown mixture
w_real <- gtools::rdirichlet(1, alpha_0)
mu_real <- rnorm(H, mu_0, sqrt(1 / tau_0))
tau_real <- rgamma(H, shape = a_0, rate = b_0)
cat("w_real :", w_real, "\n")
cat("mu_real :", mu_real, "\n")
cat("tau_real :", tau_real, "\n")


# Transition probability matrix 
# For finite mixture models, all rows are equal
Q_real <- c()
for (i in 1:H) {
  Q_real <- rbind(Q_real, w_real)
}
cat("\n Transition matrix\n")
print(Q_real)


# Generating an Hidden Markov Chain
HMC <- function(LENGTH, Q, H) {
  h <- c(1)
  for (i in 2:LENGTH) {
    h <- c(h, sample(1:H, prob = Q[h[length(h)], ], size = 1))
  }
  return(h)
}
h_real <- HMC(LENGTH, Q_real, H)
#print(h_real)


# Generating the gaussians of the mixture
x <- seq(-20, 20, by = 0.001)
norms <- data.frame(x)
col_names <- c("x")
norms <- c()
for (i in 1:H) {
  norms <- rbind(norms, dnorm(x, mu_real[i], sqrt(1 / tau_real[i])) * w_real[i])
}

# Sampling from the unknown distribution
rmix <- function(h_real, mu_real, tau_real) {
  d <- rnorm(length(h_real), mu_real[h_real], sqrt(1 / tau_real[h_real]))
  return(d)
}
d <- rmix(h_real, mu_real, tau_real)



# Plotting the Hidden Markov Model and samples

library(rgl)
plot_HMM_samples <- function(x, d, h_real, LENGTH) {
  for (i in 1:LENGTH) {
    plot3d(x, norms[h_real[i], ], rep(i, length(x)), type = "l", lwd = 4, col = h_real[i], zlim = c(1, LENGTH))
    plot3d(d[i], 0, i, size = 4, col = h_real[i], zlim = c(1, LENGTH))
  }
  #grid3d(c("x", "y+", "z"))
}
plot_HMM_samples(x, d, h_real, LENGTH)


# Full conditionals
sample_Q <- function(alpha_0, h, H) {
  Q <- matrix(, nrow = H, ncol = H)
  NN <- matrix(0, nrow = H, ncol = H)
  for (z in 2:H) {
    NN[h[z-1], h[z]] <- NN[h[z-1], h[z]] + 1
  }
    
  for (i in 1:H) {
    Q[i, ] <- gtools::rdirichlet(1, alpha_0 + NN[i, ])
  }

  return(Q)
}


sample_tau <- function(mu, h, x, a_0, b_0, N, H, LENGTH) {
  z <- matrix(0, nrow = LENGTH, ncol = H)
  for (i in 1:LENGTH) {
    z[i, h[i]] <- 1
  }

  tau <- c()
  for (c in 1:length(mu)) {
    summation <- 0
    for (i in 1:length(x)) {
      summation <- summation + z[i, c] * (x[i] - mu[c])^2
    }
    tau <- c(tau, rgamma(1, shape = a_0 + N[c] / 2, rate = b_0 + summation / 2))
  }
  return(tau)
}


sample_mu <- function(tau, z, x, tau_0, mu_0, N) {
  mu <- c()

  for (c in 1:length(tau)) {
    # to prevent division by 0
    if (N[c] == 0) {
      N[c] <- 1
    }

    delta <- N[c] * tau[c] / (tau_0 + N[c] * tau[c])
    x_bar <- (z[, c] %*% x) / N[c]
    mu <- c(mu, rnorm(1, delta * x_bar + (1 - delta) * mu_0, sqrt(1 / (tau_0 + N[c] * tau[c]))))
  }
  return(mu)
}


# Forward-Backward 
sample_h <- function(d, Q, mu, tau, LENGTH, H) {
  # Forward recursion
  P <- array(0, dim = c(H, H, LENGTH))
  pi <- matrix(0, nrow = LENGTH, ncol = H)
  pi[1, 1] <- 1
  for (t in 2:LENGTH) {
    for (r in 1:H) {
      for (s in 1:H) {
        P[r, s, t] <- pi[t - 1, r] * Q[r, s] * dnorm(d[t], mu[s], tau[s])        
      }
    }
    #print(sum(P[, , t]))
    summation <- sum(P[, , t])
    if (summation == 0) {
      P[, , t] <- matrix(1 / H^2, nrow = H, ncol = H)
    } else {
      P[, , t] <- P[, , t] / summation
    }

    for (s in 1:H) {
      pi[t, s] <-  sum(P[, s, t])
    }
  }
  #print(P)
  # Backward recursion
  h <- sample(1:H, prob = pi[LENGTH, ], size = 1)
  for (i in (LENGTH - 1):1) {
    if (sum(P[, h[length(h)], i + 1]) == 0) {
      prob <- rep(1 / H, H)
    } else {
      prob <- P[, h[length(h)], i + 1]
    }
    #print(P[, , i+1])
    h <- c(sample(1:H, prob = prob, size = 1), h)
  }
  return(h)
}


# Gibbs Sampler
gibbs <- function(d, niter, H, alpha_0, mu_0, tau_0, a_0, b_0, LENGTH) {
  cat("\n\nGibbs Sampler\n")
  w <- gtools::rdirichlet(1, alpha_0)
  Q <- c()
  for (h in 1:H) {
    Q <- rbind(Q, w)
  }
  tau <- rgamma(H, shape = a_0, rate = b_0)
  mu <- rnorm(H, mu_0, sqrt(1 / tau_0))
  h <- HMC(LENGTH, Q, H)
  print(Q)
  #print(h)
  print(mu)
  print(tau)

  mu_GS <- matrix(, nrow = niter, ncol = H)
  mu_GS[1, ] <- mu

  cat(0, "/", niter, "\n")

  for (i in 1:niter) {
    if (i %% 50 == 0) {
      cat(i, "/", niter, "\n")
    }
    z <- matrix(0, nrow = LENGTH, ncol = H)
    for (l in 1:LENGTH) {
      z[l, h[l]] <- 1
    }
    N <- c()
    for (c in 1:H) {
      N <- c(N, sum(z[, c]))
    }
    Q <- sample_Q(alpha_0, h, H)
    tau <- sample_tau(mu, h, d, a_0, b_0, N, H, LENGTH)
    mu <- sample_mu(tau, z, d, tau_0, mu_0, N)
    h <- sample_h(d, Q, mu, tau, LENGTH, H)
    mu_GS[i, ] <- mu
  }
  #print(Q)
  #print(h)
  return(mu_GS)
}

mu_GS <- gibbs(d, niter, H, alpha_0, mu_0, tau_0, a_0, b_0, LENGTH)
#print(mu_GS)
print(mu_real)
x11()
matplot(mu_GS, main="Markov Chain for mu", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)