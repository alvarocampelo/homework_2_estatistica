# q1.R - questao 1: X ~ binomial(n, p)
# importando parametros
source("R/parametros.R")   # carrega M, SEED, n, p

dir.create("figuras", showWarnings = FALSE)

k   <- 0:n                     # suporte de X
pmf <- dbinom(k, n, p)         # P(X = k)
cdf <- pbinom(k, n, p)         # F(k) = P(X <= k)

cat(sprintf("p = %.3f\n", p))


#q1. item 2
plot_pmf <- function() {
  plot(k, pmf, type = "h", lwd = 2, col = "steelblue", main = "PMF de X ~ Binomial(50; 0,405)", xlab = "Numero de requisicoes com processamento adicional (k)", ylab = "P(X = k)")
  points(k, pmf, pch = 19, cex = 0.6, col = "steelblue")
  abline(v = n * p, lty = 2, col = "gray40")
  legend("topright", legend = c("P(X = k)", "E[X] = 20,25"), col = c("steelblue", "gray40"), lty = c(1, 2), lwd = c(2, 1), bty = "n")
}

plot_pmf() #mostra na tela
png("figuras/q1_pmf.png", width = 1600, height = 1000, res = 200)
plot_pmf() #salva o arquivo
dev.off()

F_x <- stepfun(k, c(0, cdf)) #funcao escada continua a direita
plot_cdf <- function() {
  plot(F_x, verticals = FALSE, do.points = FALSE, lwd = 2, col = "steelblue", xlim = c(-2, 52), ylim = c(0, 1), main = "CDF de X ~ Binomial(50; 0,405)", xlab = "x (numero de requisicoes)", ylab = "F(x) = P(X <= x)")
  points(k, cdf, pch = 19, cex = 0.6, col = "steelblue") #valor no salto
  points(k, c(0, cdf[-length(cdf)]), pch = 1, cex = 0.6, col = "steelblue") #limite pela esquerda
}

plot_cdf()
png("figuras/q1_cdf.png", width = 1600, height = 1000, res = 200)
plot_cdf()
dev.off()


#q1. item 3
#teorico
E_teo   <- n * p
Var_teo <- n * p * (1 - p)
DP_teo  <- sqrt(Var_teo)

#R: pela definicao, somando sobre a PMF
E_R   <- sum(k * pmf)
Var_R <- sum((k - E_R)^2 * pmf)
DP_R  <- sqrt(Var_R)

tab_13 <- data.frame(Medida  = c("E[X]", "Var(X)", "DP(X)"), Teorico = c(E_teo, Var_teo, DP_teo), R = c(E_R, Var_R, DP_R))
print(tab_13, digits = 7)


#q1. item 4
#exatas pela CDF
pa <- 1 - pbinom(19, n, p) # P(X >= 20) = 1 - F(19)
pb <- pbinom(42, n, p) - pbinom(30, n, p) # P(30 < X < 43) = F(42) - F(30)
pc <- dbinom(31, n, p) # P(X = 31)

#conferencia: soma direta da PMF e formula da PMF manualmente
stopifnot(all.equal(pa, sum(dbinom(20:50, n, p))))
stopifnot(all.equal(pb, sum(dbinom(31:42, n, p))))
stopifnot(all.equal(pc, choose(50, 31) * p^31 * (1 - p)^19))

#aproximacao normal com correcao de continuidade
mu    <- n * p
sigma <- sqrt(n * p * (1 - p))
pa_N <- 1 - pnorm(19.5, mu, sigma)
pb_N <- pnorm(42.5, mu, sigma) - pnorm(30.5, mu, sigma)
pc_N <- pnorm(31.5, mu, sigma) - pnorm(30.5, mu, sigma)

tab_14 <- data.frame(Evento  = c("P(X >= 20)", "P(30 < X < 43)", "P(X = 31)"), Exato   = c(pa, pb, pc), Normal  = c(pa_N, pb_N, pc_N))
print(tab_14, digits = 4)