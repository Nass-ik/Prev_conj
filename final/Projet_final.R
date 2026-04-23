## -----------------------------------------------------------------------------------------------------------------------
library(readr)
library(tidyr)
library(dplyr)    
library(tsoutliers)
library(RJDemetra)  
library(forecast)   
library(smooth)
library(tsoutliers)
MRTSSM4453USN <- read_csv("MRTSSM4453USN.csv")
# View(MRTSSM4453USN)  # désactivé pour la compilation


## -----------------------------------------------------------------------------------------------------------------------
# Filtrage exact sur 20 ans (Janvier 2006 à Décembre 2025)
df_clean <- MRTSSM4453USN |> 
  mutate(observation_date = as.Date(observation_date)) |> 
  arrange(observation_date) |>                             
  filter(observation_date >= as.Date("2006-01-01") & 
         observation_date <= as.Date("2025-12-31"))

serie_ts <- ts(df_clean$MRTSSM4453USN, start = c(2006, 1), frequency = 12)

# Affichage du graphique 
plot(serie_ts, 
     main = "Ventes de boissons alcoolisées aux US (2006-2025)", 
     ylab = "Millions de dollars", 
     xlab = "Année",
     col = "darkblue",
     lwd = 2)


## -----------------------------------------------------------------------------------------------------------------------
# Distribution, ACF et PACF (très utile pour justifier les futurs modèles ARIMA)
ggtsdisplay(serie_ts, plot.type = "histogram", 
            main = "Distribution et Autocorrélations - Ventes d'Alcool")


## -----------------------------------------------------------------------------------------------------------------------
# Effets saisonniers (permet de bien voir le pic de décembre !)
ggseasonplot(serie_ts, col = rainbow(20), year.labels = TRUE, 
             main = "Graphique Saisonnier des Ventes d'Alcool")


## -----------------------------------------------------------------------------------------------------------------------
fit_outliers <- tso(serie_ts, types = c("AO", "LS", "TC"))


## -----------------------------------------------------------------------------------------------------------------------
plot(fit_outliers)


## -----------------------------------------------------------------------------------------------------------------------
print(fit_outliers$outliers)


## -----------------------------------------------------------------------------------------------------------------------
serie_adj <- fit_outliers$yadj


## -----------------------------------------------------------------------------------------------------------------------
#| label: tbl-outliers
#| tbl-cap: "Détection des principaux points atypiques"

library(knitr)

# Création du tableau de données synchronisé avec les résultats de fit_outliers
tableau_outliers <- data.frame(
  Date = c("Mars 2020", "Mai 2020", "Avril 2021", "Janvier 2025"),
  Type_de_point = c("LS (Level Shift)", "LS (Level Shift)", "TC (Temporary Change)", "LS (Level Shift)"),
  T_stat = c(7.917, 5.576, 4.782, -5.610)
)

# Affichage 
kable(tableau_outliers, col.names = c("Date", "Type de point", "T-stat"), digits = 3)


## -----------------------------------------------------------------------------------------------------------------------
library(fBasics)
library(tibble)
library(ggplot2)


## -----------------------------------------------------------------------------------------------------------------------
statistiques <- basicStats(serie_adj)
print(statistiques)


## -----------------------------------------------------------------------------------------------------------------------
df_box <- tibble(value = as.numeric(serie_adj))

ggplot(df_box, aes(y = value)) +
  geom_boxplot(fill = "royalblue", alpha = 0.5, outlier.color = "red") +
  labs(title = "Boxplot de la série corrigée (Ventes d'alcool)", 
       y = "Millions de dollars") +
  theme_minimal()


## -----------------------------------------------------------------------------------------------------------------------
library(TSA)
library(RJDemetra)


## -----------------------------------------------------------------------------------------------------------------------
d_serie_adj <- diff(serie_adj, differences = 1)


## -----------------------------------------------------------------------------------------------------------------------
library(seastests)

print("Test isSeasonal :")
isSeasonal(serie_adj, test = "combined", freq = 12)

print("Combined Test détaillé :")
summary(combined_test(serie_adj))


## -----------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
periodogram(serie_adj, main="Périodogramme (en niveau)")
periodogram(d_serie_adj, main="Périodogramme (différence 1ère)")
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------------------------------------------------
myregx13 <- regarima_x13(serie_adj, spec ="RG5c")
summary(myregx13)
s_transform(myregx13)
plot(myregx13)


## -----------------------------------------------------------------------------------------------------------------------
myspec <- x13_spec("RSA5c")
mysax13 <- x13(serie_adj, myspec)
summary(mysax13$regarima)
plot(mysax13$final)


## -----------------------------------------------------------------------------------------------------------------------
h_prev <- 12


## -----------------------------------------------------------------------------------------------------------------------
prev_naive <- naive(serie_adj, h = h_prev)
plot(prev_naive, main = "Prévision : Méthode Naïve", 
     ylab = "Millions de dollars", xlab = "Année")


## -----------------------------------------------------------------------------------------------------------------------
prev_x13 <- myregx13$forecast
myregx13$model$effects
print("Prévisions X13-SARIMA-SEATS :")
print(prev_x13)


## -----------------------------------------------------------------------------------------------------------------------
prev_x13_12mois <- head(myregx13$forecast, 12)
print(prev_x13_12mois)


## -----------------------------------------------------------------------------------------------------------------------
fit_stl <- stlm(serie_adj)
prev_stl <- forecast(fit_stl, h = h_prev)

plot(prev_stl, 
     main = "Prévision : Méthode STL", 
     ylab = "Millions de dollars", 
     xlab = "Année")

# Astuce anti-trou
fin_x <- time(serie_adj)[length(serie_adj)]
fin_y <- serie_adj[length(serie_adj)]
debut_x <- time(prev_stl$mean)[1]
debut_y <- prev_stl$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
fitsts <- StructTS(serie_adj)

# Graphique de la série ajustée et des résidus 
plot(cbind(fitted(fitsts), residuals(fitsts)), 
     main = "Modèle STS : Ajustement et Résidus")

# Affichage des paramètres du modèle 
print(" Variances du modèle STS ")
show(fitsts)

# Génération de la prévision
fitsts <- StructTS(serie_adj)
prev_sts <- forecast(fitsts, h=12) 

# Tracer le graphique de base (avec le trou)
plot(prev_sts, 
     main = "Prévision : Méthode STS", 
     ylab = "Millions de dollars", 
     xlab = "Année")

fin_x <- time(serie_adj)[length(serie_adj)]
fin_y <- serie_adj[length(serie_adj)]
debut_x <- time(prev_sts$mean)[1]
debut_y <- prev_sts$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue", lwd = 2)


## -----------------------------------------------------------------------------------------------------------------------
prev_hw <- hw(serie_adj, h = 12, seasonal = "multiplicative")

# 2. On trace le graphique
plot(prev_hw, 
     main = "Prévision : Holt-Winters Multiplicatif", 
     ylab = "Millions de dollars", 
     xlab = "Année")

fin_x <- time(serie_adj)[length(serie_adj)]
fin_y <- serie_adj[length(serie_adj)]
debut_x <- time(prev_hw$mean)[1]
debut_y <- prev_hw$mean[1]

lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
fit_ets <- ets(serie_adj)
summary(fit_ets) 

prev_ets <- forecast(fit_ets, h = 12)

plot(prev_ets, 
     main = "Prévision : Méthode ETS", 
     ylab = "Millions de dollars", 
     xlab = "Année")

# anti-trou
fin_x <- time(serie_adj)[length(serie_adj)]
fin_y <- serie_adj[length(serie_adj)]
debut_x <- time(prev_ets$mean)[1]
debut_y <- prev_ets$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
#| cache: true
fit_tbats <- tbats(serie_adj)
print(fit_tbats) 

prev_tbats <- forecast(fit_tbats, h = 12)

plot(prev_tbats, 
     main = "Prévision : Méthode TBATS", 
     ylab = "Millions de dollars", 
     xlab = "Année")

# Astuce anti-trou
debut_x <- time(prev_tbats$mean)[1]
debut_y <- prev_tbats$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
#| cache: true
fit_adam1 <- auto.adam(serie_adj, model="MAdA", lags=c(1,12), select=TRUE)
summary(fit_adam1)

prev_adam1 <- forecast(fit_adam1, h = 12)

plot(prev_adam1, 
     main = "Prévision : ADAM ETS", 
     ylab = "Millions de dollars", 
     xlab = "Année")

# Astuce anti-trou
debut_x <- time(prev_adam1$mean)[1]
debut_y <- prev_adam1$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
#| cache: true
fit_adam2 <- auto.adam(serie_adj, model="MAdA", lags=c(1,12), 
                       orders=list(ar=c(3,3), i=c(2), ma=c(3,3)),
                       select=TRUE)
summary(fit_adam2)

# Prévision 
prev_adam2 <- forecast(fit_adam2, h=12, level = 0.90)

plot(prev_adam2, 
     main = "Prévision : ADAM ETS SARIMA", 
     ylab = "Millions de dollars", 
     xlab = "Année")

# Astuce anti-trou
debut_x <- time(prev_adam2$mean)[1]
debut_y <- prev_adam2$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
#| cache: true
# Modèle SSARIMA (State Space ARIMA)
fit_ssarima <- auto.ssarima(serie_adj, lags=c(1,12), 
                            orders=list(ar=c(3,3), i=c(2), ma=c(3,3)), 
                            select=TRUE)
summary(fit_ssarima)

# Prévision 
prev_ssarima <- forecast(fit_ssarima, h=12, level = 0.90)

plot(prev_ssarima, 
     main = "Prévision : SSARIMA", 
     ylab = "Millions de dollars", 
     xlab = "Année")

# Astuce anti-trou
debut_x <- time(prev_ssarima$mean)[1]
debut_y <- prev_ssarima$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
# 1. Recherche automatique du meilleur modèle SARIMA
fit_sarima <- auto.arima(serie_adj)

# 2. Affichage des paramètres trouvés par l'algorithme
print("Résumé du modèle SARIMA")
summary(fit_sarima)

# 3. Prévision sur 12 mois
prev_sarima <- forecast(fit_sarima, h = 12)

# 4. Graphique
plot(prev_sarima, 
     main = "Prévision : Modèle SARIMA", 
     ylab = "Millions de dollars", 
     xlab = "Année")

# Astuce anti-trou 
fin_x <- time(serie_adj)[length(serie_adj)]
fin_y <- serie_adj[length(serie_adj)]
debut_x <- time(prev_sarima$mean)[1]
debut_y <- prev_sarima$mean[1]
lines(c(fin_x, debut_x), c(fin_y, debut_y), col = "blue")


## -----------------------------------------------------------------------------------------------------------------------
# =========================================================
# Création du tableau comparatif complet (Sécurisé)
# =========================================================

# Fonction anti-crash pour transformer les NULL en NA et garder 11 lignes
safe_val <- function(x) { if (is.null(x)) NA else as.numeric(x) }

# Calcul manuel de l'AIC et AICc pour le modèle STS 
k_sts <- length(fitsts$coef)
n_obs <- length(serie_adj)
aic_sts <- 2 * k_sts - 2 * fitsts$loglik
aicc_sts <- aic_sts + (2 * k_sts^2 + 2 * k_sts) / (n_obs - k_sts - 1)

# Construction du grand tableau de données
tableau_aic <- data.frame(
  Modèle = c("Naïve", "X13-SARIMA-SEATS", "STL", "STS", "Holt-Winters", 
             "SARIMA", "ETS", "TBATS", "SSARIMA", "ADAM ETS", "ADAM SARIMA"),
  AIC = c(
    NA,
    safe_val(myregx13$model$regarima$aic), 
    safe_val(fit_stl$model$aic),           
    aic_sts,                         
    safe_val(prev_hw$model$aic),           
    safe_val(fit_sarima$aic),
    safe_val(fit_ets$aic),
    safe_val(fit_tbats$AIC), 
    safe_val(fit_ssarima$AIC),   
    safe_val(fit_adam1$AIC),     
    safe_val(fit_adam2$AIC)      
  ),
  AICc = c(
    NA,
    NA, # X13 ne fournit pas d'AICc natif
    safe_val(fit_stl$model$aicc),
    aicc_sts,                         
    safe_val(prev_hw$model$aicc),
    safe_val(fit_sarima$aicc),
    safe_val(fit_ets$aicc),
    NA, # TBATS ne fournit pas d'AICc natif
    safe_val(fit_ssarima$AICc),  
    safe_val(fit_adam1$AICc),    
    safe_val(fit_adam2$AICc)     
  )
)

# Tri du tableau par AICc croissant 
tableau_aic <- tableau_aic |> 
  dplyr::arrange(AICc)

# Affichage final
knitr::kable(tableau_aic, digits = 2, caption = "Comparaison globale des critères d'information")


## -----------------------------------------------------------------------------------------------------------------------
# 1. On crée de vraies dates lisibles pour l'axe X (de Janvier à Décembre 2026)
dates_futures <- seq(as.Date("2026-01-01"), by = "month", length.out = 12)

# 2. Création du tableau avec les belles dates
df_forecasts <- data.frame(
  Date = dates_futures,
  `X13-ARIMA` = as.numeric(prev_x13[1:12]),
  STL = as.numeric(prev_stl$mean[1:12]),
  STS = as.numeric(prev_sts$mean[1:12]),
  `Holt-Winters` = as.numeric(prev_hw$mean[1:12]),
  ETS = as.numeric(prev_ets$mean[1:12]),
  TBATS = as.numeric(prev_tbats$mean[1:12]),
  `ADAM ETS` = as.numeric(prev_adam1$mean[1:12]),
  `ADAM SARIMA` = as.numeric(prev_adam2$mean[1:12]),
  SSARIMA = as.numeric(prev_ssarima$mean[1:12]),
  SARIMA = as.numeric(prev_sarima$mean[1:12]),
  Naïve = as.numeric(prev_naive$mean[1:12]),
  check.names = FALSE 
)

# 3. Transformation en format long
df_long <- df_forecasts |>
  pivot_longer(-Date, names_to = "Méthode", values_to = "Valeur")

# 4. Couleurs
couleurs <- c(
  "X13-ARIMA" = "darkred",
  "STL" = "blue",
  "STS" = "forestgreen",
  "Holt-Winters" = "pink",
  "ETS" = "steelblue",
  "TBATS" = "yellow3",
  "ADAM ETS" = "mediumaquamarine",
  "ADAM SARIMA" = "mediumvioletred",
  "SSARIMA" = "purple",
  "SARIMA" = "red",         
  "Naïve" = "black"         
)

# 5. Création du graphique avec un format de date propre sur l'axe X
ggplot(df_long, aes(x = Date, y = Valeur, color = Méthode)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Évolution des prévisions selon l'ensemble des modèles",
    x = "Mois (Année 2026)", 
    y = "Millions de dollars"
  ) +
  scale_color_manual(values = couleurs) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") + # L'astuce magique pour les dates !
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## -----------------------------------------------------------------------------------------------------------------------
# 1. SÉPARATION DES DONNÉES (TRAIN / TEST) 
# On cache l'année 2025 pour l'évaluation
train_test <- window(serie_adj, end = c(2024, 12))
obs_test <- window(serie_adj, start = c(2025, 1), end = c(2025, 12))

#  2. RÉ-ENTRAÎNEMENT EXPRESS DES MODÈLES
# On refait tourner les algorithmes uniquement sur le passé (2006-2024)

# Méthodes basiques et classiques
prev_naive_test <- naive(train_test, h = 12)$mean
prev_hw_test <- hw(train_test, h = 12, seasonal = "multiplicative")$mean

fit_stl_test <- stlm(train_test)
prev_stl_test <- forecast(fit_stl_test, h = 12)$mean

fit_sts_test <- StructTS(train_test)
prev_sts_test <- forecast(fit_sts_test, h = 12)$mean

# Lissage exponentiel avancé
fit_ets_test <- ets(train_test)
prev_ets_test <- forecast(fit_ets_test, h = 12)$mean

fit_tbats_test <- tbats(train_test)
prev_tbats_test <- forecast(fit_tbats_test, h = 12)$mean

# Modèles ADAM
fit_adam1_test <- auto.adam(train_test, model="MAdA", lags=c(1,12), select=TRUE)
prev_adam1_test <- forecast(fit_adam1_test, h = 12)$mean

fit_adam2_test <- auto.adam(train_test, model="MAdA", lags=c(1,12), orders=list(ar=c(3,3), i=c(2), ma=c(3,3)), select=TRUE)
prev_adam2_test <- forecast(fit_adam2_test, h = 12)$mean

# Modèles ARIMA
fit_ssarima_test <- auto.ssarima(train_test, lags=c(1,12), orders=list(ar=c(3,3), i=c(2), ma=c(3,3)), select=TRUE)
prev_ssarima_test <- forecast(fit_ssarima_test, h = 12)$mean

fit_sarima_test <- auto.arima(train_test)
prev_sarima_test <- forecast(fit_sarima_test, h = 12)$mean

# (Note : On exclut X13 ici car sa fonction regarima_x13 est complexe à automatiser dans une boucle de test sans créer d'erreurs, les 10 autres modèles suffisent largement pour le classement).

# 3. REGROUPEMENT DES PRÉVISIONS TESTS 
df_test <- data.frame(
  STL = as.numeric(prev_stl_test),
  STS = as.numeric(prev_sts_test),
  `Holt-Winters` = as.numeric(prev_hw_test),
  ETS = as.numeric(prev_ets_test),
  TBATS = as.numeric(prev_tbats_test),
  `ADAM ETS` = as.numeric(prev_adam1_test),
  `ADAM SARIMA` = as.numeric(prev_adam2_test),
  SSARIMA = as.numeric(prev_ssarima_test),
  SARIMA = as.numeric(prev_sarima_test),
  Naïve = as.numeric(prev_naive_test),
  check.names = FALSE
)


## -----------------------------------------------------------------------------------------------------------------------
# Calcul du MSE de la méthode Naïve (notre référence)
mse_naive <- mean((as.numeric(obs_test) - df_test$Naïve)^2)

results <- data.frame(Modèle = character(), MSE = numeric(), R2OOS = numeric())

# Boucle d'évaluation
for (mod in colnames(df_test)) {
  pred <- df_test[[mod]]
  err <- as.numeric(obs_test) - pred
  mse <- mean(err^2)
  
  if(mod == "Naïve") {
    r2oos <- NA 
  } else {
    r2oos <- 1 - (mse / mse_naive)
  }
  
  results <- rbind(
    results,
    data.frame(Modèle = mod, MSE = round(mse, 2), R2OOS = round(r2oos, 4))
  )
}

# Tri du meilleur (MSE le plus bas) au pire
results <- results |> arrange(MSE)

# Affichage du tableau final
kable(results, caption = "Tableau des performances Out-of-Sample (Test sur l'année 2025)")






## -----------------------------------------------------------------------------------------------------------------------
cspe_data <- df_test
for(col in names(cspe_data)) {
  cspe_data[[col]] <- cumsum((cspe_data[[col]] - as.numeric(obs_test))^2)
}

# Ajout des dates de 2025
cspe_data$Date <- seq(as.Date("2025-01-01"), by = "month", length.out = 12)

# Passage en format long ET exclusion de la méthode Naïve pour corriger l'échelle !
cspe_long <- cspe_data |>
  pivot_longer(-Date, names_to = "Modèle", values_to = "CSPE") |>
  filter(Modèle != "Naïve")

# Réutilisation de tes couleurs
couleurs_test <- c(
  "STL" = "blue", "STS" = "forestgreen", "Holt-Winters" = "pink",
  "ETS" = "steelblue", "TBATS" = "yellow3", "ADAM ETS" = "mediumaquamarine",
  "ADAM SARIMA" = "mediumvioletred", "SSARIMA" = "purple",
  "SARIMA" = "red", "Naïve" = "black"
)

ggplot(cspe_long, aes(x = Date, y = CSPE, color = Modèle)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Évolution de l'Erreur Quadratique Cumulative (CSPE) sur 2025",
    subtitle = "Le meilleur modèle est celui dont la courbe termine le plus bas en décembre",
    x = "Mois de l'année 2025",
    y = "Erreur cumulée (CSPE)"
  ) +
  scale_color_manual(values = couleurs_test) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## -----------------------------------------------------------------------------------------------------------------------
# 1. Création et prévision du modèle AR(1) sur le train_test
fit_ar1 <- Arima(train_test, order=c(1,0,0))
prev_ar1 <- forecast(fit_ar1, h = 12)$mean

# 2. Calcul de l'erreur du nouveau benchmark AR(1)
error_ar1 <- as.numeric(obs_test) - as.numeric(prev_ar1)

# 3. Création du tableau des résultats
results_dm <- data.frame(Modèle = character(), DM_p_value = numeric())

# 4. Boucle de test pour chaque modèle
for (mod in colnames(df_test)) {
  if (mod == "Naïve") next  # On peut ignorer la naïve maintenant
  
  pred <- df_test[[mod]]
  err <- as.numeric(obs_test) - pred
  
  # Calcul du test par rapport à l'AR(1)
  dm_result <- tryCatch({
    dm_test <- dm.test(err, error_ar1, alternative = "less", h = length(obs_test))
    round(dm_test$p.value, 4)
  }, error = function(e) { NA })
  
  results_dm <- rbind(results_dm, data.frame(Modèle = mod, DM_p_value = dm_result))
}

# 5. Tri et affichage
results_dm <- results_dm |> arrange(DM_p_value)
kable(results_dm, caption = "Test de Diebold-Mariano (Référence : AR(1), p-value < 0.05 = meilleur)")


## -----------------------------------------------------------------------------------------------------------------------
library(forecast)
library(ggplot2)
library(tidyr)

# --- 1. CALCUL DE LA BOUCLE ---
h <- 12        
format_rolling <- numeric(h)

for(i in 1:h) {
  end_year <- 2024 + floor((10 + i) / 12)
  end_month <- (10 + i) %% 12 + 1
  
  train_rolling <- window(serie_adj, end = c(end_year, end_month))
  fit_roll <- auto.arima(train_rolling)
  prev_roll <- forecast(fit_roll, h = 1)
  
  format_rolling[i] <- as.numeric(prev_roll$mean)
}

# --- 2. CREATION DU GRAPHIQUE DIRECTEMENT APRES ---
valeurs_reelles <- window(serie_adj, start = c(2025, 1), end = c(2025, 12))
dates_2025 <- seq(as.Date("2025-01-01"), by = "month", length.out = 12)

df_final <- data.frame(
  Date = dates_2025,
  Prévision_SARIMA = format_rolling,
  Réalité = as.numeric(valeurs_reelles)
)

df_final_long <- tidyr::pivot_longer(df_final, -Date, names_to = "Série", values_to = "Valeur")

ggplot(df_final_long, aes(x = Date, y = Valeur, color = Série)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Prévision_SARIMA" = "orchid", "Réalité" = "mediumaquamarine")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(
    title = "Validation Finale : SARIMA Glissant vs Réalité (2025)",
    subtitle = "Alignement corrigé : prévision à 1 mois réestimée chaque mois",
    x = NULL, y = "Millions de dollars"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

