\texttt{rbetabinom}
rm(list = ls())

library(tikzDevice)
library(cmdstanr) # 0.5.2, with CmdStan v 2.28.2
library(posterior) # v 1.3.1

dir.create("figures")

plot_circle <- function(x, y, r = 1, prop = 1, ...) {
  thetas <- seq(0, prop * 2*pi, length.out = 100)
  xs <- x + cos(thetas) * r
  ys <- y + sin(thetas) * r
  if (prop == 1) {
    polygon(xs, ys, ...)
  } else if (prop < 1) {
    xs <- c(x, xs, x)
    ys <- c(y, ys, y)
    polygon(xs, ys, ...)
  } else {
    stop ("prop invalid!")
  }
}
  
Fsts <- c(0.1, 0.5, 0.9)
mus <- c(0.1, 0.5)
dat <- expand.grid(Fst = Fsts, mu = mus)
dat$seed <- c(847, 56628, 54623, 68953, 54680, 10612)

tikz("figures/pieExamples.tex", height = 4, width = 6)

par(mfrow = c(2, 3))

for (j in 1:nrow(dat)) {

  set.seed(dat$seed[j])
  n_groups <- 9
  mu <- dat$mu[j]
  Fst <- dat$Fst[j]
  theta <- (1 - Fst) / Fst
  p <- rbeta(n_groups, mu * theta, (1 - mu) * theta)

  par(mar = c(0, 0, 0, 0))

  plot(NULL, axes = FALSE, frame.plot = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = c(0, 4), ylim = c(0, 4))
  if ((j - 1) %% 3 == 0) {
    mtext(paste0("$\\overline{x}$: ", mu), 2, -1.5, cex = 1.0)
    mtext("$\\overbrace{\\hspace{3.2cm}}$", 2, -2.8, cex = 1.0)
  }
  if ((j-1) < 3) {
    mtext(paste0("$F_{ST}$: ", Fst), 3, -1.5, cex = 1.0)
    mtext("$\\overbrace{\\hspace{3.2cm}}$", 3, -2.8, cex = 1.0)
  }

  light <- "#FF8A89"
  light_border <- "#FF8A89"
  dark <- "#55E1FF"
  dark_border <- "#55E1FF"

  light <- "#FFCE86"
  light_border <- "#FFCE86"
  dark <- "#00AFB9"
  dark_border <- "#00AFB9"
  
  for (row in 1:3) {
    for (col in 1:3) {
      plot_circle(row, col, r = 0.4, col = light, border = light_border)
      if (p[3*(row-1) + col] > 0.01) plot_circle(row, col, r = 0.4, prop = p[3*(row-1) + col], col = dark, border = dark_border)
    }
  }

}

dev.off()


rbetabinom <- function(m, n, mu, F) {
  theta <- (1 - F) / F
  p <- rbeta(m, mu * theta, (1 - mu) * theta)
  x <- rbinom(m, n, prob = p)
  return(x)
}



model <- cmdstan_model("betabinomial_fst.stan")

pars <- expand.grid(
  n = c(10, 50),
  m = 20,
  x = seq(0.1, 0.9, 0.1),
  F = seq(0.1, 0.9, 0.1)
)

est <- data.frame(
  x_mu = rep(NA, nrow(pars)),
  F_mu = rep(NA, nrow(pars))
)


for (i in 1:nrow(pars)) {

  stan_dat <- list(
    n = rep(pars$n[i], pars$m[i]),
    m = pars$m[i],
    x = rbetabinom(pars$m[i], pars$n[i], pars$x[i], pars$F[i])
  )

  fit <- model$sample(
    data = stan_dat, 
    seed = 123,
    iter_warmup = floor(1000/2), iter_sampling = 1000,
    chains = 4, 
    parallel_chains = 4,
    refresh = 50
  )

  samples <- as_draws_rvars(fit$draws())
  est$F_mu[i] <- mean(draws_of(samples$F))
  est$x_mu[i] <- mean(draws_of(samples$mu))

  print(i)

}

saveRDS(pars, file = "figures/pars.RDS")
saveRDS(est, file = "figures/est.RDS")

pars <- readRDS("figures/pars.RDS")
est <- readRDS("figures/est.RDS")

offset <- 0.01

tikz("figures/sim_est.tex", height = 3.5, width = 6)

par(mfrow = c(1, 2))

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "true $\\overline{x}$", ylab = "estimated $\\overline{x}$")
tar <- which(pars$n == 10)
points(pars$x[tar] - offset, est$x_mu[tar], pch = 16)
tar <- which(pars$n == 50)
points(pars$x[tar] + offset, est$x_mu[tar], col = "red", pch = 16)
abline(0, 1, lty = 2)

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "true $F_{ST}$", ylab = "estimated $F_{ST}$")
tar <- which(pars$n == 10)
points(pars$F[tar] - offset, est$F_mu[tar], pch = 16)
tar <- which(pars$n == 50)
points(pars$F[tar] + offset, est$F_mu[tar], col = "red", pch = 16)
abline(0, 1, lty = 2)

dev.off()