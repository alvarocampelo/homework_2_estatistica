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


#q1. item 5
#item 5 - c)

set.seed(SEED) # semente unica, definida em parametros.R

y <- rbinom(n, size = 1, prob = p) # 50 observacoes Bernoulli(p)
k_obs <- sum(y) # numero de sucessos observados
p_hat <- k_obs / n # MLE analitico: k/n

# verificacao numerica: maximizando a logverossimilhanca
loglik <- function(pp) dbinom(k_obs, n, pp, log = TRUE)
p_hat_num <- optimize(loglik, interval = c(1e-6, 1 - 1e-6), maximum = TRUE)$maximum

cat(sprintf("k = %d | p_hat = %.4f | p_hat (numerico) = %.6f | p teorico = %.3f\n", k_obs, p_hat, p_hat_num, p))
cat(sprintf("Diferenca p_hat - p = %+.4f\n", p_hat - p))

p_grid <- seq(0.15, 0.75, by = 0.001)

plot_loglik <- function() {
  plot(p_grid, loglik(p_grid), type = "l", lwd = 2, col = "steelblue", main = sprintf("Log-verossimilhanca para k = %d, n = %d", k_obs, n), xlab = "p", ylab = "Log-verossimilhanca  l(p | k)")
  abline(v = p_hat, col = "firebrick", lwd = 2)
  abline(v = p, col = "gray40", lty = 2, lwd = 2)
  legend("bottomright",
         legend = c("l(p | k)", sprintf("MLE p_hat = %.2f", p_hat), sprintf("p teorico = %.3f", p)),
         col = c("steelblue", "firebrick", "gray40"),
         lty = c(1, 1, 2), lwd = 2, bty = "n")
}
plot_loglik()
png("figuras/q1_loglik.png", width = 1600, height = 1000, res = 200)
plot_loglik()
dev.off()

#item 5 - c)
# erro padrao teorico do MLE
ep_teo <- function(nn) sqrt(p * (1 - p) / nn)
cat(sprintf("EP teorico (n = 50): %.4f\n", ep_teo(n)))

# faixa de 95% das estimativas com n = 50
cat(sprintf("95%% das estimativas com n = 50 ficam entre %.2f e %.2f\n", qbinom(0.025, n, p) / n, qbinom(0.975, n, p) / n))

# simulacao: 1000 amostras para cada tamanho amostral.
# a soma de nn Bernoullis(p) e uma Binomial(nn, p), entao rbinom(R, nn, p) / nn
# gera diretamente R estimativas p_hat, cada uma de uma amostra de tamanho nn.
R_rep    <- 1000
tamanhos <- c(50, 500, 5000)

p_hat_sim <- lapply(tamanhos, function(nn) rbinom(R_rep, nn, p) / nn)
names(p_hat_sim) <- paste0("n = ", tamanhos)

tab_15d <- data.frame(
  n = tamanhos,
  EP_teorico = ep_teo(tamanhos),
  DP_simulado = sapply(p_hat_sim, sd),
  Media_p_hat = sapply(p_hat_sim, mean)
)
print(tab_15d, digits = 4)
# histogramas comparados
plot_sim <- function() {
  par(mfrow = c(1, 3))
  for (i in seq_along(tamanhos)) {
    hist(p_hat_sim[[i]], breaks = 30, xlim = c(0.2, 0.6), col = "lightsteelblue", border = "white", main = names(p_hat_sim)[i], xlab = "Estimativa p_hat", ylab = "Frequencia")
    abline(v = p, col = "firebrick", lty = 2, lwd = 2)
  }
  par(mfrow = c(1, 1))
}

plot_sim()
png("figuras/q1_mle_tamanho_amostral.png", width = 2400, height = 900, res = 200)
plot_sim()
dev.off()