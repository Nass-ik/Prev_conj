
# SCRIPT AUTONOME : PREVISION A 1 PAS


# 1. Chargement des packages
library(readr)
library(dplyr)
library(tsoutliers)
library(forecast)

# 2. Importation et Filtrage (20 ans)
MRTSSM4453USN <- read_csv("MRTSSM4453USN.csv", show_col_types = FALSE)
df_clean <- MRTSSM4453USN |> 
  mutate(observation_date = as.Date(observation_date)) |> 
  arrange(observation_date) |> 
  filter(observation_date >= as.Date("2006-01-01") & observation_date <= as.Date("2025-12-31"))

serie_ts <- ts(df_clean$MRTSSM4453USN, start = c(2006, 1), frequency = 12)

# 3. Correction des outliers
fit_outliers <- tso(serie_ts, types = c("AO", "LS", "TC"))
serie_adj <- fit_outliers$yadj

# 4. Boucle de prévision Rolling Window (SARIMA)
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

# 5. Sauvegarde du vecteur de prévision dans un fichier texte
write(t(format_rolling), file = "for_sarima.out", ncolumn = 1, append = FALSE)
print("Script terminé ! Le fichier for_sarima.out a été généré.")