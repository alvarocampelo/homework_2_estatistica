# Questão 3 - Simulação da temperatura da CPU (Box-Muller vs rnorm)

# 1. Imports das demais questões e parâmetros globais
source("R/parametros.R") # carrega M, SEED, n, p

q1 <- new.env()
invisible(capture.output(source("R/questao1.R", local = q1)))

# 2. Dados herdados da Questão 1 e de parametros.R
n_q1 <- n
p_q1 <- p
m    <- q1$E_teo          # E[X] calculado na Q1 (item 3)
n_T  <- round(100 * m)

# Parâmetros da temperatura da CPU (próprios da Q3)
mu    <- 62
sigma <- 3.5

# Gerador Normal(0, 1) por Box-Muller
box_muller <- function(n) {
  k <- ceiling(n / 2)
  u1 <- runif(k)
  u2 <- runif(k)
  r <- sqrt(-2 * log(u1))
  th <- 2 * pi * u2
  c(r * cos(th), r * sin(th))[1:n]
}

set.seed(SEED)
t_bm <- mu + sigma * box_muller(n_T)
t_rn <- rnorm(n_T, mu, sigma)

# Valores teóricos (usando a CDF F_x importada da Q1)
r_teo  <- 1 - q1$F_x(floor(m))
t_crit <- qnorm(1 - r_teo, mu, sigma)

p68_teo    <- 1 - pnorm(68, mu, sigma)
p60_65_teo <- pnorm(65, mu, sigma) - pnorm(60, mu, sigma)

# Extremos típicos para n_T observações (P(Z > z) = ln(2)/n_T)
z_tip   <- qnorm(1 - log(2) / n_T)
min_tip <- mu - sigma * z_tip
max_tip <- mu + sigma * z_tip

resumo <- function(x) {
  c(
    Media = mean(x),
    DP = sd(x),
    Min = min(x),
    Max = max(x),
    P_gt68 = mean(x > 68),
    P_60_65 = mean(x > 60 & x < 65),
    P_gtTcrit = mean(x > t_crit)
  )
}

teo <- c(mu, sigma, min_tip, max_tip, p68_teo, p60_65_teo, r_teo)

tab <- data.frame(
  Teorico = round(teo, 4),
  Box_Muller = round(resumo(t_bm), 4),
  rnorm = round(resumo(t_rn), 4)
)

cat("n_T =", n_T, "| r =", round(r_teo, 6), "| T_crit =", round(t_crit, 4), "\n\n")
print(tab)

# Checagem dos extremos observados (z-score e probabilidade na amostra)
for (nm in c("Box-Muller", "rnorm")) {
  x <- if (nm == "Box-Muller") t_bm else t_rn
  z_min <- (min(x) - mu) / sigma
  z_max <- (max(x) - mu) / sigma
  p_min <- 1 - (1 - pnorm(z_min))^n_T
  p_max <- 1 - (1 - pnorm(-z_max))^n_T
  cat(sprintf("\n%s -> min: %.2f (z=%.2f, P=%.3f) | max: %.2f (z=%.2f, P=%.3f)",
              nm, min(x), z_min, p_min, max(x), z_max, p_max))
}
cat("\n")

# Histogramas + PDF teórica
x_seq <- seq(mu - 4 * sigma, mu + 4 * sigma, length.out = 500)
y_pdf <- dnorm(x_seq, mu, sigma)
t_crit_lbl <- sub("\\.", ",", sprintf("%.2f", t_crit))

par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1.0), oma = c(3.0, 0, 2.2, 0))

hist(t_bm, breaks = 30, freq = FALSE, xlim = c(48, 76), ylim = c(0, 0.135),
     col = "#B0C4DE", border = "white", las = 1,
     main = "Box-Muller", xlab = "Temperatura da CPU (°C)", ylab = "Densidade")
lines(x_seq, y_pdf, col = "#8B0000", lwd = 2)
abline(v = t_crit, col = "#2F4F4F", lty = 2, lwd = 1.5)
box()

hist(t_rn, breaks = 30, freq = FALSE, xlim = c(48, 76), ylim = c(0, 0.135),
     col = "#B0C4DE", border = "white", las = 1,
     main = "Gerador nativo (rnorm)", xlab = "Temperatura da CPU (°C)", ylab = "Densidade")
lines(x_seq, y_pdf, col = "#8B0000", lwd = 2)
abline(v = t_crit, col = "#2F4F4F", lty = 2, lwd = 1.5)
box()

mtext(sprintf("Temperaturas Simuladas da CPU (n = %d)", n_T),
      outer = TRUE, cex = 1.1, font = 2, line = 0.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom",
       legend = c("PDF Teórica N(62; 3,5²)", paste0("T_crit = ", t_crit_lbl, " °C")),
       col = c("#8B0000", "#2F4F4F"), lty = c(1, 2), lwd = c(2, 1.5),
       horiz = TRUE, bty = "n", cex = 0.9, inset = c(0, 0.02))

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(5.1, 4.1, 4.1, 2.1))

# Verificações citadas no relatório
z_crit <- qnorm(1 - r_teo)
cat("\nz_crit =", round(z_crit, 4), "| T_crit =", round(mu + sigma * z_crit, 4), "\n")
cat("P(T>68) =", round(1 - pnorm(68, mu, sigma), 4),
    "| P(60<T<65) =", round(pnorm(65, mu, sigma) - pnorm(60, mu, sigma), 4), "\n")

cat("P(55<T<69): BM =", round(mean(t_bm > 55 & t_bm < 69), 4),
    "| rnorm =", round(mean(t_rn > 55 & t_rn < 69), 4),
    "| teórico =", round(pnorm(2) - pnorm(-2), 4), "\n")

ep_media <- sigma / sqrt(n_T)
ep_s     <- sigma / sqrt(2 * (n_T - 1))
cat("EP(media) =", round(ep_media, 4), "| EP(s) =", round(ep_s, 4), "\n")
cat("Média em EP: BM =", round((mean(t_bm) - mu) / ep_media, 2),
    "| rnorm =", round((mean(t_rn) - mu) / ep_media, 2), "\n")
cat("DP em EP:    BM =", round((sd(t_bm) - sigma) / ep_s, 2),
    "| rnorm =", round((sd(t_rn) - sigma) / ep_s, 2), "\n")

ep_prop <- function(p) sqrt(p * (1 - p) / n_T)
cat("Desvios das proporções (em EP):\n")
print(round(rbind(
  `T>68`    = c(BM = mean(t_bm > 68) - p68_teo, rn = mean(t_rn > 68) - p68_teo) / ep_prop(p68_teo),
  `60<T<65` = c(BM = mean(t_bm > 60 & t_bm < 65) - p60_65_teo,
                rn = mean(t_rn > 60 & t_rn < 65) - p60_65_teo) / ep_prop(p60_65_teo),
  `T>Tcrit` = c(BM = mean(t_bm > t_crit) - r_teo, rn = mean(t_rn > t_crit) - r_teo) / ep_prop(r_teo)
), 2))

cat("P(T>Tcrit) esperado com a média do rnorm =",
    round(1 - pnorm(t_crit, mean(t_rn), sigma), 4), "\n")
