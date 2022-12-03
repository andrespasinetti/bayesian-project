# Implementing the following paper:
# https://epub.jku.at/obvulihs/download/pdf/5554146?originalFilename=true

library("mvtnorm")
# PARAMETERS
C <- 4 # number of classes
alpha_0 <- runif(n = C, min = 1, max = 4)
p <- 2
mu_0 <- rep(0, p)
tau_0 <- diag(rep(0.02, p)) 
a_0 <- 5
b_0 <- diag(rep(0.05, p)) 


N_SAMPLES <- 100 # number of samples
niter <- 300 # number of Gibbs Sampling iterations


# Defining the unknown mixture

w_real <- gtools::rdirichlet(1, alpha_0)
mu_real <- matrix(, nrow = C, ncol = p)
mu_real <- rmvnorm(C, mu_0, solve(tau_0))
tau_real <- array(0, dim = c(p, p, C))
tau_real[, , ] <- rWishart(C, a_0, b_0)


cat("\nw_real :", w_real, "\n")
cat("\nmu_real\n")
print(mu_real)
#cat("mu_real :", mu_real, "\n")
cat("\ntau_real\n")
print(tau_real)
#cat("tau_real :", tau_real, "\n")


# Sampling from the unknown distribution
rmix <- function(w_real, mu_real, tau_real) {
  z_real <- sample(1:C, prob = w_real, size = N_SAMPLES, replace = TRUE)
  x <- matrix(, nrow = N_SAMPLES, ncol = p)
  for (i in 1:N_SAMPLES) {
    x[i, ] <- rmvnorm(1, mu_real[z_real[i], ], solve(tau_real[, , z_real[i]]))
  }
  
  return(list(x, z_real))
}
mix <- rmix(w_real, mu_real, tau_real)
x <- mix[1]
z_real <- mix[2]

x <- data.frame(x)
z_real <- unlist(z_real)
x11()
plot(x, type="p", col=z_real, pch=20)

# Plotting
plot_mixture <- function(x_seq, w, mu, tau, title) {
  norms <- data.frame(x_seq)
  col_names <- c("x_seq")
  for (i in 1:C) {
    norms <- cbind(norms, dnorm(x_seq, mu[i], sqrt(1 / tau[i])) * w[i])
    col_names <- c(col_names, paste("norm", i, sep = ""))
  }

  norms_matrix <- data.matrix(norms[2:length(norms)])
  mixture <- rowSums(norms_matrix)
  norms <- cbind(norms, mixture)

  col_names <- c(col_names, "mixture")
  colnames(norms) <- col_names
  library("tidyr")
  long_norms <- pivot_longer(norms, cols = all_of(col_names[2:length(col_names)]))

  library(ggplot2)
  line_sizes <- rep(0.6, C + 1)
  line_sizes[1] <- 2
  line_colors <- rep(1, C + 1)
  line_colors[1] <- 2
  ggp <- ggplot(long_norms) + 
        geom_line(aes(x_seq, value, col = name, size = name)) + 
        scale_size_manual(values = line_sizes) +
        scale_color_manual(values = line_colors) +
        ggtitle(title) +
        theme(plot.title = element_text(size = 20, face = "bold"))

  # x_seq <- data.frame(x_seq)
  # ggp <- ggp +
  #   geom_histogram(data = x_seq, aes(x_seq = x_seq, y = after_stat(density)), alpha = 0.3,
  #                 bins = 300, position = "identity", lwd = 0.2) +
  #   ggtitle("Unknown mixture of gaussians + samples")

  # ggp <- x_seq[[1]]

  x11(type = "cairo")
  plot(ggp)
}
#x_seq <- seq(min(x), max(x), by = 0.001)
#title <- "Unknown mixture of gaussians"
#plot_mixture(x_seq, w_real, mu_real, tau_real, title)


# Full conditionals
sample_w <- function(alpha_0, N) {
  alpha_new <- alpha_0 + N
  w <- gtools::rdirichlet(1, alpha_new)
  return(w)
}

sample_tau <- function(mu, z, x, c_0, C_0, N) {
  tau <- array(, dim = c(p, p, C))

  for (c in 1:C) {
    summation <- t(sweep(x[z == c, ], p, mu[c, ], "-")) %*% (sweep(x[z == c, ], p, mu[c, ], "-"))
    c_new <- c_0 + N[c] / 2
    C_new <- C_0 + summation / 2
    tau <- rWishart(C, c_new, solve(C_new))
  }
  return(tau)
}

sample_mu <- function(tau, z, x, b_0, B_0, N) {
  mu <- matrix(, nrow = C, ncol = p)
  for (c in 1:C) {
    B_new <- solve(solve(B_0) + N[c] * tau[, , c])


    # avoid divison by 0
    if (N[c] == 0) {
      N[c] <- 1
    }
    b_new <- (B_new) %*% (solve(B_0) %*% b_0 + N[c] * (tau[, , c] %*% apply(x[z == c,], p, sum)/N[c]))
    mu[c, ] <- rmvnorm(C, b_new, B_new)
  }
  return(mu)
}

sample_z <- function(mu, tau, w, x) {
  z <- matrix(, nrow = N_SAMPLES, ncol = C)
  for (i in 1:N_SAMPLES) {
    prob <- c()
    summation <- 0
    for (c in 1:C) {
      summation <- summation + w[c] * dmvnorm(x[i, ], mu[c, ], solve(tau[, , c]))
    }
    # avoid division by 0
    if (summation != 0) {
      for (c in 1:C) {
        prob <- c(prob, w[c] * dmvnorm(x[i, ], mu[c, ], solve(tau[, , c])) / summation)
      }
    } else {
      prob <- runif(n = C, min = 0, max = 1)
    }
    print("-----")
    z[i, ] <- rmultinom(1, 1, prob)
  }
  return(z)
}


# Gibbs Sampler
gibbs <- function(x, niter, C, alpha_0, mu_0, tau_0, a_0, b_0) {
  x <- data.matrix(x)

  # Hyperparameters
  c_0 <- 2.5 * 0.5 * (p - 1)
  c_0 <- 2
  phi <- 1
  C_0 <- phi * cov(x)
  b_0 <- apply(x, p, median)
  B_0 <- diag((apply(x, p, max) - apply(x, p, min))^2)
  # Initialize the Markov Chain
  w <- gtools::rdirichlet(1, alpha_0) * 0.01

  mu <- matrix(, nrow = C, ncol = p)
  tau <- array(0, dim = c(p, p, C))
  tau[, , ]  <- rWishart(C, c_0, solve(C_0))
  mu <- rmvnorm(C, b_0, B_0)
  z <- matrix(, nrow = N_SAMPLES, ncol = C)
  for (i in 1:N_SAMPLES) {
    z[i, ] <- rmultinom(1, 1, w)
  }

  # Save the Markov Chain of mu
  mu_GS <- array(dim = c(C, p, niter))
  mu_GS[, , 1] <- mu

  cat("\nGibbs Sampling\n")
  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 50 == 0) {
      cat(i, "/", niter, "\n")
    }

    N <- c()
    for (c in 1:C) {
      N <- c(N, sum(z[, c]))
    }

    k <- array(0, dim = N_SAMPLES)
    for (i in 1:N_SAMPLES) {
      k[i] <- which.max(z[i, ])
    }
    z <- k

    w <- sample_w(alpha_0, N)
    tau <- sample_tau(mu, z, x, c_0, C_0, N)
    mu <- sample_mu(tau, z, x, b_0, B_0, N)
    z <- sample_z(mu, tau, w, x)

    mu_GS[, , i] <- mu
  }

  return(mu_GS)
}

mu_GS <- gibbs(x, niter, C, alpha_0, mu_0, tau_0, a_0, b_0)
print(mu_GS)
x11()
matplot(mu_GS, main="Markov Chain for mu", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)
