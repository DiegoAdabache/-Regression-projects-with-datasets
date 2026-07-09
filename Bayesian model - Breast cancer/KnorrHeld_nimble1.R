# =============================================================================
# Proyecto final: Modelo dinámico–espacial de Knorr-Held
# Incidencia de cáncer de mama en México por estado (2003–2024)
# Manuel McCadden | Regresión Avanzada — ITAM
#
# Conversión de R2jags -> nimble
# =============================================================================

library(nimble)
library(nimbleCarbon) # para dcar_proper si se prefiere; aquí usamos la
# formulación manual via dmnorm con Q construida en R

wdir <- "C:/Users/diego/OneDrive/Documentos/ProyectoRegres/Proyecto_final"
setwd(wdir)
# -----------------------------------------------------------------------------
# Funciones auxiliares
# -----------------------------------------------------------------------------

prob <- function(x) min(length(x[x > 0]) / length(x),
                        length(x[x < 0]) / length(x))

# Resumen de convergencia para un parámetro escalar
plot_convergence <- function(samples, par_name) {
  # samples: objeto mcmc.list de coda (salida de runMCMC con samplesAsCodaMCMC=TRUE)
  z1 <- as.vector(samples[[1]][, par_name])
  z2 <- as.vector(samples[[2]][, par_name])
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  plot(z1, type = "l", col = "grey50",
       main = paste("Traza:", par_name), ylab = par_name, xlab = "Iteración")
  lines(z2, col = "firebrick2")
  y1 <- cumsum(z1) / seq_along(z1)
  y2 <- cumsum(z2) / seq_along(z2)
  plot(y1, type = "l", col = "grey50",
       main = "Media ergódica", ylab = par_name, xlab = "Iteración")
  lines(y2, col = "firebrick2")
  acf(z1, main = "ACF (cadena 1)")
}

# DIC manual desde muestras NIMBLE
# Requiere que el modelo haya monitoreado "logProb_y" o se recalcule.
# Aquí lo calculamos directamente desde las muestras de lambda.
calc_DIC <- function(samples_mat, model_obj) {
  # Método estándar: DIC = Dbar + pD, pD = Dbar - D(theta_bar)
  # Delegamos al summary de nimble si está disponible, o lo calculamos
  # a mano si se monitorearon los log-likelihoods.
  # Ver sección de ajuste de cada modelo para el cálculo concreto.
  message("Usar calc_DIC_from_lambda() después de extraer lambda del summary.")
}

# DIC desde la matriz de log-verosimilitudes monitoreadas
calc_DIC_from_loglik <- function(loglik_samples) {
  D_samples <- -2 * rowSums(loglik_samples)
  Dbar      <- mean(D_samples)
  # Usar la varianza de la devianza muestral para estimar pD
  pD        <- var(D_samples) / 2 
  DIC       <- Dbar + pD
  c(DIC = DIC, pD = pD, Dbar = Dbar)
}
# -----------------------------------------------------------------------------
# 1. Lectura y preparación de datos  (idéntica a la versión JAGS)
# -----------------------------------------------------------------------------

d      <- read.csv("analysis_dataset.csv", stringsAsFactors = FALSE)
estados <- sort(unique(d$state))
anios   <- sort(unique(d$year))
S <- length(estados)
T <- length(anios)
cat("S =", S, "estados | T =", T, "años (", min(anios), "-", max(anios), ")\n")

y   <- matrix(0L, S, T, dimnames = list(estados, as.character(anios)))
E   <- matrix(0,  S, T, dimnames = list(estados, as.character(anios)))
pop <- matrix(0,  S, T, dimnames = list(estados, as.character(anios)))

for (i in 1:nrow(d)) {
  si <- match(d$state[i], estados); ti <- match(d$year[i], anios)
  y[si, ti]   <- d$y[i]
  E[si, ti]   <- d$E[i]
  pop[si, ti] <- d$pop_fem[i]
}
stopifnot(all(E > 0), all(pop > 0), !any(is.na(y)))

tasa <- y / pop * 1e5
rbar <- sum(y) / sum(pop)

# Adyacencia
vec <- list(
  "Aguascalientes"    = c("Zacatecas","Jalisco"),
  "Baja California"   = c("Baja California Sur","Sonora"),
  "Baja California Sur" = c("Baja California"),
  "Campeche"          = c("Tabasco","Quintana Roo","Yucatán"),
  "Chiapas"           = c("Tabasco","Veracruz","Oaxaca"),
  "Chihuahua"         = c("Sonora","Sinaloa","Durango","Coahuila"),
  "Ciudad de México"  = c("México","Morelos"),
  "Coahuila"          = c("Chihuahua","Durango","Zacatecas","San Luis Potosí","Nuevo León"),
  "Colima"            = c("Jalisco","Michoacán"),
  "Durango"           = c("Chihuahua","Sinaloa","Nayarit","Jalisco","Zacatecas","Coahuila"),
  "Guanajuato"        = c("Jalisco","Zacatecas","San Luis Potosí","Querétaro","Michoacán"),
  "Guerrero"          = c("Michoacán","México","Morelos","Puebla","Oaxaca"),
  "Hidalgo"           = c("San Luis Potosí","Querétaro","México","Tlaxcala","Puebla","Veracruz"),
  "Jalisco"           = c("Nayarit","Durango","Zacatecas","Aguascalientes","San Luis Potosí","Guanajuato","Michoacán","Colima"),
  "Michoacán"         = c("Colima","Jalisco","Guanajuato","Querétaro","México","Guerrero"),
  "Morelos"           = c("México","Ciudad de México","Puebla","Guerrero"),
  "México"            = c("Michoacán","Querétaro","Hidalgo","Tlaxcala","Puebla","Morelos","Guerrero","Ciudad de México"),
  "Nayarit"           = c("Sinaloa","Durango","Zacatecas","Jalisco"),
  "Nuevo León"        = c("Coahuila","Zacatecas","San Luis Potosí","Tamaulipas"),
  "Oaxaca"            = c("Guerrero","Puebla","Veracruz","Chiapas"),
  "Puebla"            = c("Hidalgo","México","Morelos","Guerrero","Oaxaca","Veracruz","Tlaxcala"),
  "Querétaro"         = c("Guanajuato","San Luis Potosí","Hidalgo","México","Michoacán"),
  "Quintana Roo"      = c("Campeche","Yucatán"),
  "San Luis Potosí"   = c("Coahuila","Nuevo León","Tamaulipas","Veracruz","Hidalgo","Querétaro","Guanajuato","Jalisco","Zacatecas"),
  "Sinaloa"           = c("Sonora","Chihuahua","Durango","Nayarit"),
  "Sonora"            = c("Baja California","Chihuahua","Sinaloa"),
  "Tabasco"           = c("Campeche","Chiapas","Veracruz"),
  "Tamaulipas"        = c("Nuevo León","San Luis Potosí","Veracruz"),
  "Tlaxcala"          = c("Hidalgo","Puebla","México"),
  "Veracruz"          = c("Tamaulipas","San Luis Potosí","Hidalgo","Puebla","Oaxaca","Chiapas","Tabasco"),
  "Yucatán"           = c("Campeche","Quintana Roo"),
  "Zacatecas"         = c("Durango","Coahuila","Nuevo León","San Luis Potosí","Aguascalientes","Jalisco","Nayarit","Guanajuato")
)

W <- matrix(0, S, S, dimnames = list(estados, estados))
for (a in names(vec)) for (b in vec[[a]]) { W[a, b] <- 1; W[b, a] <- 1 }
num  <- rowSums(W)
zero <- rep(0, S)
cat("Vecinos: min", min(num), "máx", max(num),
    "| simétrica:", isSymmetric(W),
    "| sin aislados:", all(num > 0), "\n")

# -----------------------------------------------------------------------------
# Nota sobre la matriz de precisión CAR en NIMBLE
# -----------------------------------------------------------------------------
# En JAGS construíamos Q dentro del modelo con un doble for. En NIMBLE es más
# eficiente construir Q.phi en R y pasarla como constante, evitando recomputarla
# en cada iteración.  Se actualiza cuando tau.phi o rho cambian usando la
# función build_Q() definida abajo.
#
# Alternativa: usar dcar_proper() de nimble >= 0.13 directamente, lo que
# elimina la necesidad de construir Q manualmente.  Se muestran ambas formas:
# la versión dmnorm (fiel al original) y la versión dcar_proper (más elegante).
# Aquí usamos dmnorm para conservar la estructura original.

# -----------------------------------------------------------------------------
# 2. Datos e iniciales comunes
# -----------------------------------------------------------------------------

constants <- list(S = S, T = T, num = num, zero = zero, W = W)
data_kh   <- list(y = y, E = E)

inits_base <- function() list(
  alpha     = 0,
  tau.phi   = 1,
  tau.theta = 1,
  tau.gamma = 1,
  rho       = 0.5,
  phi       = rep(0, S),
  theta     = rep(0, S),
  gamma     = rep(0, T)
)

inits_delta <- function() c(
  inits_base(),
  list(tau.delta = 1,
       delta     = matrix(0, S, T))
)

# Parámetros a monitorear
pars_base  <- c("alpha", "rho", "tau.phi", "tau.theta", "tau.gamma",
                "phi.c", "gamma.c", "lambda", "loglik")
pars_delta <- c(pars_base, "tau.delta")

# Configuración MCMC (misma lógica que el original)
mcmc_base  <- list(niter = 15000, nburnin = 3000,  nthin = 3,  nchains = 2)
mcmc_tipoI <- list(niter = 50000, nburnin = 6000,  nthin = 6,  nchains = 2)
mcmc_heavy <- list(niter = 70000, nburnin = 10000, nthin = 10, nchains = 2)

# Función genérica de ajuste
fit_nimble <- function(code, constants, data, inits_fn, monitors, mcmc_cfg,
                       seed = 123) {
  set.seed(seed)
  model <- nimbleModel(code       = code,
                       constants  = constants,
                       data       = data,
                       inits      = inits_fn(),
                       calculate  = FALSE)
  cmodel  <- compileNimble(model)
  conf    <- configureMCMC(model, monitors = monitors, thin = mcmc_cfg$nthin)
  mcmc    <- buildMCMC(conf)
  cmcmc   <- compileNimble(mcmc, project = model)
  samples <- runMCMC(cmcmc,
                     niter             = mcmc_cfg$niter,
                     nburnin           = mcmc_cfg$nburnin,
                     nchains           = mcmc_cfg$nchains,
                     samplesAsCodaMCMC = TRUE,
                     setSeed           = seed + seq_len(mcmc_cfg$nchains) - 1)
  samples
}

# =============================================================================
# 3. MODELO BASE  (efectos principales, sin interacción)
# =============================================================================

code_base <- nimbleCode({
  # Verosimilitud
  for (s in 1:S) {
    for (t in 1:T) {
      y[s, t]      ~ dpois(mu[s, t])
      mu[s, t]     <- E[s, t] * lambda[s, t]
      log(lambda[s, t]) <- alpha + phi.c[s] + theta.c[s] + gamma.c[t]
      loglik[s, t] <- dpois(y[s, t], mu[s, t], log = TRUE)
    }
  }
  
  # Efecto espacial estructurado: CAR propio via dmnorm con Q construida inline
  phi[1:S] ~ dmnorm(zero[1:S], prec = Q.phi[1:S, 1:S])
  for (i in 1:S) {
    for (j in 1:S) {
      Q.phi[i, j] <- tau.phi * (equals(i, j) * (num[i] + 0.001) - rho * W[i, j])
    }
  }
  for (s in 1:S) { phi.c[s] <- phi[s] - mean(phi[1:S]) }
  
  # Heterogeneidad no estructurada
  for (s in 1:S) { theta[s] ~ dnorm(0, tau = tau.theta) }
  for (s in 1:S) { theta.c[s] <- theta[s] - mean(theta[1:S]) }
  
  # Tendencia temporal RW1
  gamma[1] ~ dnorm(0, var = 100)
  for (t in 2:T) { gamma[t] ~ dnorm(gamma[t - 1], tau = tau.gamma) }
  for (t in 1:T) { gamma.c[t] <- gamma[t] - mean(gamma[1:T]) }
  
  # Previas
  alpha     ~ dnorm(0, var = 100)
  tau.phi   ~ dgamma(1, 0.1)
  tau.theta ~ dgamma(1, 0.1)
  tau.gamma ~ dgamma(1, 0.1)
  rho       ~ dunif(0, 0.99)
})

cat("\n=== Ajustando modelo BASE ===\n")
samples_base <- fit_nimble(code_base, constants, data_kh,
                           inits_base, pars_base, mcmc_base)

# DIC base
loglik_mat_base <- as.matrix(do.call(rbind, samples_base))[,
                                                           grep("^loglik\\[", colnames(as.matrix(samples_base[[1]])))]
dic_base <- calc_DIC_from_loglik(loglik_mat_base)
cat("BASE  DIC =", round(dic_base["DIC"], 1),
    "pD =",  round(dic_base["pD"],  1), "\n")

# =============================================================================
# 4. TIPO I  (no estructurada × no estructurada)
# =============================================================================

code_tipoI <- nimbleCode({
  for (s in 1:S) {
    for (t in 1:T) {
      y[s, t]      ~ dpois(mu[s, t])
      mu[s, t]     <- E[s, t] * lambda[s, t]
      log(lambda[s, t]) <- alpha + phi.c[s] + theta.c[s] + gamma.c[t] + delta[s, t]
      loglik[s, t] <- dpois(y[s, t], mu[s, t], log = TRUE)
    }
  }
  
  phi[1:S] ~ dmnorm(zero[1:S], prec = Q.phi[1:S, 1:S])
  for (i in 1:S) {
    for (j in 1:S) {
      Q.phi[i, j] <- tau.phi * (equals(i, j) * (num[i] + 0.001) - rho * W[i, j])
    }
  }
  for (s in 1:S) { phi.c[s] <- phi[s] - mean(phi[1:S]) }
  for (s in 1:S) { theta[s] ~ dnorm(0, tau = tau.theta) }
  for (s in 1:S) { theta.c[s] <- theta[s] - mean(theta[1:S]) }
  
  gamma[1] ~ dnorm(0, var = 100)
  for (t in 2:T) { gamma[t] ~ dnorm(gamma[t - 1], tau = tau.gamma) }
  for (t in 1:T) { gamma.c[t] <- gamma[t] - mean(gamma[1:T]) }
  
  # Interacción Tipo I: iid
  for (s in 1:S) {
    for (t in 1:T) { delta[s, t] ~ dnorm(0, tau = tau.delta) }
  }
  
  
  alpha     ~ dnorm(0, var = 100)
  tau.phi   ~ dgamma(1, 0.1)
  tau.theta ~ dgamma(1, 0.1)
  tau.gamma ~ dgamma(1, 0.1)
  tau.delta ~ dgamma(1, 0.1)
  rho       ~ dunif(0, 0.99)
})

cat("\n=== Ajustando modelo TIPO I ===\n")
samples_I <- fit_nimble(code_tipoI, constants, data_kh,
                        inits_delta, pars_delta, mcmc_tipoI)

loglik_mat_I <- as.matrix(do.call(rbind, samples_I))[,
                                                     grep("^loglik\\[", colnames(as.matrix(samples_I[[1]])))]
dic_I <- calc_DIC_from_loglik(loglik_mat_I)
cat("TIPO I  DIC =", round(dic_I["DIC"], 1),
    "pD =",  round(dic_I["pD"],  1), "\n")

# =============================================================================
# 5. TIPO II  (no estructurada × RW1 temporal)
# =============================================================================

code_tipoII <- nimbleCode({
  for (s in 1:S) {
    for (t in 1:T) {
      y[s, t]      ~ dpois(mu[s, t])
      mu[s, t]     <- E[s, t] * lambda[s, t]
      log(lambda[s, t]) <- alpha + phi.c[s] + theta.c[s] + gamma.c[t] + delta.c[s, t]
      loglik[s, t] <- dpois(y[s, t], mu[s, t], log = TRUE)
    }
  }
  
  phi[1:S] ~ dmnorm(zero[1:S], prec = Q.phi[1:S, 1:S])
  for (i in 1:S) {
    for (j in 1:S) {
      Q.phi[i, j] <- tau.phi * (equals(i, j) * (num[i] + 0.001) - rho * W[i, j])
    }
  }
  for (s in 1:S) { phi.c[s] <- phi[s] - mean(phi[1:S]) }
  for (s in 1:S) { theta[s] ~ dnorm(0, tau = tau.theta) }
  for (s in 1:S) { theta.c[s] <- theta[s] - mean(theta[1:S]) }
  
  gamma[1] ~ dnorm(0, var = 100)
  for (t in 2:T) { gamma[t] ~ dnorm(gamma[t - 1], tau = tau.gamma) }
  for (t in 1:T) { gamma.c[t] <- gamma[t] - mean(gamma[1:T]) }
  
  # Interacción Tipo II: RW1 independiente por estado
  for (s in 1:S) {
    delta[s, 1] ~ dnorm(0, tau = tau.delta)
    for (t in 2:T) { delta[s, t] ~ dnorm(delta[s, t - 1], tau = tau.delta) }
  }
  for (s in 1:S) {
    for (t in 1:T) { delta.c[s, t] <- delta[s, t] - mean(delta[s, 1:T]) }
  }
  
  alpha     ~ dnorm(0, var = 100)
  tau.phi   ~ dgamma(1, 0.1)
  tau.theta ~ dgamma(1, 0.1)
  tau.gamma ~ dgamma(1, 0.1)
  tau.delta ~ dgamma(1, 0.1)
  rho       ~ dunif(0, 0.99)
})

cat("\n=== Ajustando modelo TIPO II ===\n")
samples_II <- fit_nimble(code_tipoII, constants, data_kh,
                         inits_delta, pars_delta, mcmc_tipoI)

loglik_mat_II <- as.matrix(do.call(rbind, samples_II))[,
                                                       grep("^loglik\\[", colnames(as.matrix(samples_II[[1]])))]
dic_II <- calc_DIC_from_loglik(loglik_mat_II)
cat("TIPO II  DIC =", round(dic_II["DIC"], 1),
    "pD =",  round(dic_II["pD"],  1), "\n")



# =============================================================================
# 6. TIPO III  (CAR espacial × no estructurada temporal)
# =============================================================================

code_tipoIII <- nimbleCode({
  for (s in 1:S) {
    for (t in 1:T) {
      y[s, t]      ~ dpois(mu[s, t])
      mu[s, t]     <- E[s, t] * lambda[s, t]
      log(lambda[s, t]) <- alpha + phi.c[s] + theta.c[s] + gamma.c[t] + delta.c[s, t]
      loglik[s, t] <- dpois(y[s, t], mu[s, t], log = TRUE)
    }
  }
  
  # Efecto espacial principal
  phi[1:S] ~ dmnorm(zero[1:S], prec = Q.phi[1:S, 1:S])
  for (i in 1:S) {
    for (j in 1:S) {
      Q.phi[i, j] <- tau.phi * (equals(i, j) * (num[i] + 0.001) - rho * W[i, j])
    }
  }
  for (s in 1:S) { phi.c[s] <- phi[s] - mean(phi[1:S]) }
  for (s in 1:S) { theta[s] ~ dnorm(0, tau = tau.theta) }
  for (s in 1:S) { theta.c[s] <- theta[s] - mean(theta[1:S]) }
  
  gamma[1] ~ dnorm(0, var = 100)
  for (t in 2:T) { gamma[t] ~ dnorm(gamma[t - 1], tau = tau.gamma) }
  for (t in 1:T) { gamma.c[t] <- gamma[t] - mean(gamma[1:T]) }
  
  # Interacción Tipo III: CAR independiente en cada periodo
  for (i in 1:S) {
    for (j in 1:S) {
      Q.delta[i, j] <- tau.delta * (equals(i, j) * (num[i] + 0.001) - rho * W[i, j])
    }
  }
  for (t in 1:T) { delta[1:S, t] ~ dmnorm(zero[1:S], prec = Q.delta[1:S, 1:S]) }
  for (s in 1:S) {
    for (t in 1:T) { delta.c[s, t] <- delta[s, t] - mean(delta[1:S, t]) }
  }
  
  alpha     ~ dnorm(0, var = 100)
  tau.phi   ~ dgamma(1, 0.1)
  tau.theta ~ dgamma(1, 0.1)
  tau.gamma ~ dgamma(1, 0.1)
  tau.delta ~ dgamma(1, 0.1)
  rho       ~ dunif(0, 0.99)
})

cat("\n=== Ajustando modelo TIPO III ===\n")
samples_III <- fit_nimble(code_tipoIII, constants, data_kh,
                          inits_delta, pars_delta, mcmc_heavy)

loglik_mat_III <- as.matrix(do.call(rbind, samples_III))[,
                                                         grep("^loglik\\[", colnames(as.matrix(samples_III[[1]])))]
dic_III <- calc_DIC_from_loglik(loglik_mat_III)
cat("TIPO III  DIC =", round(dic_III["DIC"], 1),
    "pD =",  round(dic_III["pD"],  1), "\n")

# =============================================================================
# 7. TIPO IV  (CAR × RW1 — totalmente acoplado)
# =============================================================================

code_tipoIV <- nimbleCode({
  for (s in 1:S) {
    for (t in 1:T) {
      y[s, t]      ~ dpois(mu[s, t])
      mu[s, t]     <- E[s, t] * lambda[s, t]
      log(lambda[s, t]) <- alpha + phi.c[s] + theta.c[s] + gamma.c[t] + delta.c[s, t]
      loglik[s, t] <- dpois(y[s, t], mu[s, t], log = TRUE)
    }
  }
  
  # Efecto espacial principal
  phi[1:S] ~ dmnorm(zero[1:S], prec = Q.phi[1:S, 1:S])
  for (i in 1:S) {
    for (j in 1:S) {
      Q.phi[i, j] <- tau.phi * (equals(i, j) * (num[i] + 0.001) - rho * W[i, j])
    }
  }
  for (s in 1:S) { phi.c[s] <- phi[s] - mean(phi[1:S]) }
  for (s in 1:S) { theta[s] ~ dnorm(0, tau = tau.theta) }
  for (s in 1:S) { theta.c[s] <- theta[s] - mean(theta[1:S]) }
  
  gamma[1] ~ dnorm(0, var = 100)
  for (t in 2:T) { gamma[t] ~ dnorm(gamma[t - 1], tau = tau.gamma) }
  for (t in 1:T) { gamma.c[t] <- gamma[t] - mean(gamma[1:T]) }
  
  # Interacción Tipo IV: CAR en espacio + RW1 en tiempo (producto Kronecker)
  for (i in 1:S) {
    for (j in 1:S) {
      Q.delta[i, j] <- tau.delta * (equals(i, j) * (num[i] + 0.001) - rho * W[i, j])
    }
  }
  delta[1:S, 1] ~ dmnorm(zero[1:S], prec = Q.delta[1:S, 1:S])
  for (t in 2:T) { delta[1:S, t] ~ dmnorm(delta[1:S, t - 1], prec = Q.delta[1:S, 1:S]) }
  for (s in 1:S) {
    for (t in 1:T) { delta.c[s, t] <- delta[s, t] - mean(delta[s, 1:T]) - mean(delta[1:S, t]) + mean(delta[1:S, 1:T]) }
  }
  
  alpha     ~ dnorm(0, var = 100)
  tau.phi   ~ dgamma(1, 0.1)
  tau.theta ~ dgamma(1, 0.1)
  tau.gamma ~ dgamma(1, 0.1)
  tau.delta ~ dgamma(1, 0.1)
  rho       ~ dunif(0, 0.99)
})

cat("\n=== Ajustando modelo TIPO IV ===\n")
samples_IV <- fit_nimble(code_tipoIV, constants, data_kh,
                         inits_delta, pars_delta, mcmc_heavy)

loglik_mat_IV <- as.matrix(do.call(rbind, samples_IV))[,
                                                       grep("^loglik\\[", colnames(as.matrix(samples_IV[[1]])))]
dic_IV <- calc_DIC_from_loglik(loglik_mat_IV)
cat("TIPO IV  DIC =", round(dic_IV["DIC"], 1),
    "pD =",  round(dic_IV["pD"],  1), "\n")

# =============================================================================
# 8. Convergencia (Tipo IV)
# =============================================================================

plot_convergence(samples_IV, "rho")
plot_convergence(samples_IV, "alpha")

# Resumen de hiperparámetros y Rhat via coda
sm_IV <- summary(samples_IV)
hp    <- c("alpha", "rho", "tau.phi", "tau.theta", "tau.gamma", "tau.delta")
print(round(sm_IV$statistics[hp, ], 4))
print(round(sm_IV$quantiles[hp, ], 4))

library(coda)
gd <- gelman.diag(samples_IV, multivariate = FALSE)
cat("\nMáx Rhat (todos los nodos):", round(max(gd$psrf[, 1]), 3),
    " (deseable < 1.1)\n")

# =============================================================================
# 9. Comparación de modelos (DIC)
# =============================================================================

tab_dic <- rbind(
  Base    = dic_base,
  "Tipo I"  = dic_I,
  "Tipo II" = dic_II,
  "Tipo III"= dic_III,
  "Tipo IV" = dic_IV
)
tab_dic <- round(tab_dic, 1)
tab_dic <- tab_dic[order(tab_dic[, "DIC"]), ]
print(tab_dic)
mejor <- rownames(tab_dic)[1]
cat("\n>> Modelo seleccionado por DIC:", mejor, "<<\n")

options(repr.plot.width = 10, repr.plot.height = 5)
par(bg = "white", mar = c(5, 4, 3, 1))
bp <- barplot(tab_dic[, "DIC"], col = "steelblue", las = 2, ylab = "DIC",
              main = "Comparación por DIC",
              ylim = c(min(tab_dic[, "DIC"]) * .999,
                       max(tab_dic[, "DIC"]) * 1.001), xpd = FALSE)
text(bp, tab_dic[, "DIC"], round(tab_dic[, "DIC"]), pos = 3, cex = .8)

# =============================================================================
# 10. Resultados e interpretación
# =============================================================================

# Seleccionar muestras del modelo ganador
samples_list <- list(Base     = samples_base,
                     "Tipo I"  = samples_I,
                     "Tipo II" = samples_II,
                     "Tipo III"= samples_III,
                     "Tipo IV" = samples_IV)
best_samples <- samples_list[[mejor]]
best_mat     <- as.matrix(do.call(rbind, best_samples))

# Media posterior de lambda[s,t]
lam_cols <- grep("^lambda\\[", colnames(best_mat))
lam_mean  <- matrix(colMeans(best_mat[, lam_cols]), S, T,
                    dimnames = list(estados, as.character(anios)))
tasa_suav <- lam_mean * rbar * 1e5

# Mapa último año (requiere mapa_mexico — definir igual que en el notebook)
# mapa_mexico() requiere terra, RColorBrewer, classInt; se omite aquí porque
# depende del shapefile local, pero el código siguiente es directo drop-in.
#
# stl1 <- c("Aguascalientes","Baja California", ...)   # orden STL-1
# library(terra); library(RColorBrewer); library(classInt)
# mexico.map <- vect("Mexico/shapes/")
# poly_cnt   <- as.numeric(table(mexico.map$`STL-1`))
# mapa_mexico <- function(vals, titulo, nclr=7, pal="Reds") { ... }
#
# ult <- T; a_ult <- anios[ult]
# par(mfrow=c(1,2), bg="white")
# mapa_mexico(setNames(tasa[, ult], estados),
#             paste0("Tasa CRUDA x100k — ", a_ult))
# mapa_mexico(setNames(tasa_suav[, ult], estados),
#             paste0("Tasa SUAVIZADA (", mejor, ") x100k — ", a_ult))

# Tendencia temporal y efecto espacial
sm_best <- summary(best_samples)
g_idx   <- grep("^gamma.c\\[", rownames(sm_best$statistics))
g_mean  <- sm_best$statistics[g_idx, "Mean"]
g_q     <- sm_best$quantiles[g_idx, c("2.5%", "97.5%")]

options(repr.plot.width = 14, repr.plot.height = 6)
par(mfrow = c(1, 2), bg = "white", mar = c(4, 4, 3, 1))
plot(anios, g_mean, type = "l", col = "firebrick2", lwd = 2,
     ylim = range(g_q),
     xlab = "Año", ylab = "gamma_t", main = "Tendencia temporal (RW1)")
lines(anios, g_q[, 1], lty = 2, col = "grey50")
lines(anios, g_q[, 2], lty = 2, col = "grey50")

pe_idx <- grep("^phi.c\\[", rownames(sm_best$statistics))
pe     <- sm_best$statistics[pe_idx, "Mean"]
names(pe) <- estados
par(mar = c(9, 4, 3, 1))
barplot(sort(pe), las = 2,
        col = ifelse(sort(pe) > 0, "firebrick", "steelblue"),
        main = "Efecto espacial estructurado (phi.c)",
        ylab = "log-RR", cex.names = .5)

rho_stats <- sm_best$statistics["rho", ]
rho_q     <- sm_best$quantiles["rho", ]
cat("rho: media", round(rho_stats["Mean"], 3),
    " IC95% [", round(rho_q["2.5%"], 3), ",", round(rho_q["97.5%"], 3), "]\n")


load("ckpt_total_Tipo_I")

# 1. Cargar la librería necesaria para diagnósticos MCMC
library(coda)

# 2. Definir la función gráfica (por si están en una sesión limpia)
plot_convergence <- function(samples, par_name) {
  z1 <- as.vector(samples[[1]][, par_name])
  z2 <- as.vector(samples[[2]][, par_name])
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), bg = "white")
  plot(z1, type = "l", col = "grey50",
       main = paste("Traza:", par_name), ylab = par_name, xlab = "Iteración")
  lines(z2, col = "firebrick2")
  y1 <- cumsum(z1) / seq_along(z1)
  y2 <- cumsum(z2) / seq_along(z2)
  plot(y1, type = "l", col = "grey50",
       main = "Media ergódica", ylab = par_name, xlab = "Iteración")
  lines(y2, col = "firebrick2")
  acf(z1, main = "ACF (cadena 1)")
}


# Al hacer load(), se cargan en la memoria los objetos:
# samples_list, m_mat, dic_calc, yf_smp y t_final

# ------------------------------------------------------------------------------
# 4. ANÁLISIS VISUAL (Trazas y Media Ergódica)
# ------------------------------------------------------------------------------
cat("\nGenerando gráficas... (Usa las flechas del panel 'Plots' para navegar)\n")

# Evaluar hiperparámetros e intercepto global
plot_convergence(samples_list, "alpha")
plot_convergence(samples_list, "tau.delta")
plot_convergence(samples_list, "tau.phi")
plot_convergence(samples_list, "tau.theta") # Varianza del ruido espacial
plot_convergence(samples_list, "tau.gamma") # Varianza de la tendencia temporal
plot_convergence(samples_list, "tau.delta")

# Opcionalmente, evaluar un nodo específico de interacción
# plot_convergence(samples_list, "delta[1, 1]")

# ------------------------------------------------------------------------------
# 5. ANÁLISIS NUMÉRICO (Gelman-Rubin / R-hat)
# ------------------------------------------------------------------------------
cat("\n=== Estadístico Gelman-Rubin (R-hat) Global ===\n")
# Calculamos el R-hat para los miles de parámetros a la vez
gd_I <- gelman.diag(samples_list, multivariate = FALSE)
rh_I <- gd_I$psrf[, 1]
rh_I <- rh_I[is.finite(rh_I)] # Evitar errores si alguna varianza es cero

cat(sprintf("Max Rhat global: %.3f   Mediana Rhat: %.3f\n", max(rh_I), median(rh_I)))
cat(sprintf("Parámetros con Rhat > 1.1: %d de %d\n", sum(rh_I > 1.1), length(rh_I)))

# 6. Desglose del R-hat agrupado por familia de parámetros
cat("\n=== R-hat por grupo de parámetros ===\n")
pnames <- names(rh_I)
groups <- gsub("\\[.*", "", pnames) # Limpiar índices (ej. "delta[1,2]" -> "delta")

for (g in unique(groups)) {
  rh_g <- rh_I[groups == g]
  if (length(rh_g) == 0) next
  cat(sprintf("  %-12s n=%4d   max %.2f   q95 %.2f   med %.2f\n",
              g, length(rh_g), max(rh_g), quantile(rh_g, 0.95), median(rh_g)))
}

# 7. Resumen de tiempo y DIC
cat("\n=== Desempeño del Modelo ===\n")
cat("Tiempo de ejecución en paralelo:", round(t_final, 1), "minutos\n")
cat("DIC del Tipo I:", round(dic_calc["DIC"], 1), "\n")
