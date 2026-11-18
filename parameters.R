# =============================================================================
# parametros.R valores compartilhados por todas as questoes do HW2
# Responsavel: Alvaro
# Uso: no inicio de cada script, rodar  source("R/parametros.R")
# =============================================================================

# Entrada: maior numero de matricula do grupo
# matriculas: 578061, 585644, 588071
M <- 588071

# Semente do projeto
SEED <- 2  # TODO: confirmar com o grupo

# Tamanho do lote de requisições
n <- 50

# Probabilidade de uma requisição exigir processamento adicional — Equação (1)
p <- 0.05 + (M %% 100) / 200