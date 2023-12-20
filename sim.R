
rm(list = ls())

library(tikzDevice)

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
