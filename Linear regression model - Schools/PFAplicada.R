wdir <- "C:/Users/diego/OneDrive/Documentos/PROYECTOAP2"
setwd(wdir)

library(corrplot)
library(lmtest)
library(car)
library(ggplot2)
library(reshape2)

# Cargar los datos
datos_crudos <- read.csv("student-mat.csv", sep = ";", stringsAsFactors = TRUE)

datos <- subset(datos_crudos, G3 > 0) # hay que quitar a los que se dieron de baja 

cat("Observaciones originales:", nrow(datos_crudos), "\n")
cat("Observaciones tras quitar G3=0:", nrow(datos), "\n\n")

# Datos de entrenamiento (para la parte de predicción, vamos a trabajar con los datos de entrenamiento)
set.seed(123)
n   <- nrow(datos) 
idx <- sample(1:n, size = 0.8 * n)

train <- datos[idx, ]
test  <- datos[-idx, ]  # Este conjunto se guarda en un "cajón" hasta la Fase 2.

cat("Observaciones de entrenamiento:", nrow(train), "\n")
cat("Observaciones de prueba:", nrow(test), "\n\n")

# 1.1 Correlación con datos de entrenamiento
datos_numericos <- train[, sapply(train, is.numeric)]
matriz_cor <- cor(datos_numericos, use = "complete.obs")

cor_melt <- melt(matriz_cor)
ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, limit = c(-1, 1), name = "Correlación") +
  geom_text(aes(label = round(value, 2)), size = 3, color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mapa de calor (Datos de Entrenamiento sin G3=0)", x = "", y = "")


cat("\n3. Multicolinealidad (VIF):\n")
# Solo calcular VIF si hay más de 1 variable predictora en el modelo completo
if(length(datos_numericos) > 1) {
  print(vif(modelo_completo))
} else {
  cat("El modelo reducido solo tiene una variable predictora. No aplica VIF.\n")
}



# Reducir el modelo
# Excluimos G1 y G2 para encontrar relaciones estructurales/sociodemográficas
modelo_completo <- lm(G3 ~ . - G1 - G2, data = train)

modelo_reducido <- step(modelo_completo, 
                        direction = "both",
                        trace = 0) # cambiar a 1 si queremos ver las iteraciones

summary(modelo_reducido)

# 1.3 Gráfica de R2 Ajustado Iterativo
variables <- all.vars(formula(modelo_reducido))[-1]
r2_adj <- numeric(length(variables))

for (i in 1:length(variables)) {
  formula_i <- as.formula(paste("G3 ~", paste(variables[1:i], collapse = " + ")))
  modelo_i  <- lm(formula_i, data = train)
  r2_adj[i] <- summary(modelo_i)$adj.r.squared
}

df_r2 <- data.frame(
  n_variables = 1:length(variables),
  variable_agregada = variables,
  R2_ajustado = r2_adj
)

grafica_r2 <- ggplot(df_r2, aes(x = n_variables, y = R2_ajustado)) +
  geom_line(color = "#2166ac", linewidth = 1.2) +
  geom_point(color = "#d73027", size = 4) +
  geom_text(aes(label = paste0(variable_agregada, "\n", round(R2_ajustado*100, 2), "%")), 
            vjust = -0.8, size = 3) +
  # Expande el eje X un 10% a la derecha e izquierda para que el primer y último texto no se corten
  scale_x_continuous(breaks = 1:length(variables), 
                     expand = expansion(mult = c(0.1, 0.1))) +
  # Expande el eje Y un 20% hacia arriba para darle espacio al vjust = -0.8
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme_minimal() +
  # Añade margen al marco de la imagen (Arriba, Derecha, Abajo, Izquierda)
  theme(plot.margin = margin(t = 20, r = 20, b = 10, l = 10)) +
  labs(title = "R² ajustado al agregar variables al modelo (entrenamiento)",
       x = "Número de variables agregadas (Columnas del dataset)", 
       y = "R² ajustado")

print(grafica_r2)

# 1.4 Diagnóstico de Residuos del Modelo Final
cat("\n--- DIAGNÓSTICO DE SUPUESTOS MATEMÁTICOS ---\n")
cat("1. Normalidad (Shapiro-Wilk):\n")
print(shapiro.test(residuals(modelo_reducido)))

cat("\n2. Homocedasticidad (Breusch-Pagan):\n")
print(bptest(modelo_reducido))

# Gráficas de diagnóstico nativas de R
par(mfrow = c(2,2))
plot(modelo_reducido, main = "Diagnóstico con datos entrenamiento")
par(mfrow = c(1,1))


# Comparación formal
anova(modelo_reducido, modelo_completo)
AIC(modelo_completo, modelo_reducido)


# FASE 2: PREDICCIÓN (ahora con datos de prueba)
predicciones <- predict(modelo_reducido, newdata = test)
residuos_test <- test$G3 - predicciones

RMSE  <- sqrt(mean(residuos_test^2))
MAE   <- mean(abs(residuos_test))
R2_pred <- 1 - sum(residuos_test^2) / sum((test$G3 - mean(test$G3))^2)

cat("\n--- MÉTRICAS DE PREDICCIÓN (TEST - DATOS NUEVOS) ---\n")
cat("RMSE:", round(RMSE, 3), "\n")
cat("MAE: ", round(MAE, 3), "\n")
cat("R² predicción:", round(R2_pred, 3), "\n")

# Gráfica final: Valores Reales vs Predichos
plot(test$G3, predicciones,
     xlab = "Calificación Real (G3)",
     ylab = "Calificación Predicha",
     main = "Capacidad Predictiva: Reales vs Predichos",
     pch  = 16, col = "steelblue",
     xlim = c(2, 20),   # Fija el límite del eje X de 0 a 20
     ylim = c(2, 20))   # Fija el límite del eje Y de 0 a 20)
abline(0, 1, col = "red", lwd = 2)

