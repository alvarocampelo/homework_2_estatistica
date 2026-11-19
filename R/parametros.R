
# parametros.R valores compartilhados por todas as questoes do HW2
# no inicio de cada script, rodar  source("R/parametros.R")
# =============================================================================

# entrada: maior numero de matricula do grupo
# matriculas: 578061, 585644, 588071
M <- 588071

# semente do projeto
SEED <- 2  # TODO: confirmar com o grupo

# tamanho do lote de requisicoes
n <- 50

# probabilidade de uma requisicao exigir processamento adicional - Equacao (1)
p <- 0.05 + (M %% 100) / 200