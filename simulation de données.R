library(tidyr)
library(dplyr)
library(ggplot2)

set.seed(1)  # Pour la reproductibilité
temps = 196

simulate_td_curve <- function(time, max_value, growth_rate, decline_rate, replicat_id, noise_factor) {
  # Générer des variations spécifiques à chaque réplicat
  noisy_max_value <- max_value * (1 + rnorm(1, mean = 0, sd = noise_factor) * replicat_id)
  noisy_growth_rate <- growth_rate * (1 + rnorm(1, mean = 0, sd = noise_factor) * replicat_id)
  noisy_decline_rate <- decline_rate * (1 + rnorm(1, mean = 0, sd = noise_factor) * replicat_id)
  
  # Combiner croissance et décroissance
  growth <- noisy_max_value * (1 - exp(-noisy_growth_rate * time))  # Phase de croissance
  decline <- exp(-((time - 24)^2) / (2 * noisy_decline_rate^2))  # Phase de décroissance
  growth * decline
}

# Générer les données simulées
simulated_data <- expand.grid(
  Time = seq(0, temps, by = 1),  # Pas de temps
  Proportion = c(1, 10, 100, 1000),  # Modalités
  Replicat = 1:3  # Réplicats
) %>%
  mutate(
    # Paramètres dépendant de la modalité
    max_value_td = case_when(
      Proportion == 1 ~ 60000,
      Proportion == 10 ~ 120000,
      Proportion == 100 ~ 150000,
      Proportion == 1000 ~ 180000
    ),
    growth_rate_td = case_when(
      Proportion == 1 ~ 0.1,
      Proportion == 10 ~ 0.08,
      Proportion == 100 ~ 0.06,
      Proportion == 1000 ~ 0.04
    ),
    decline_rate_td = case_when(
      Proportion == 1 ~ 50,
      Proportion == 10 ~ 50,
      Proportion == 100 ~ 50,
      Proportion == 1000 ~ 48
    ),
    # Simuler Td alive avec un bruit unique par réplicat
    `Td alive` = simulate_td_curve(
      Time, max_value_td, growth_rate_td, decline_rate_td, Replicat, noise_factor = 0.1
    ),
    # Simuler Sc alive (similaire, mais sans dépendance forte à la modalité)
    `Sc alive` = simulate_td_curve(
      Time, max_value_td * 0.8, growth_rate_td * 1.2, decline_rate_td * 1.2, Replicat, noise_factor = 30
    ),
    `Yeast conc` = `Td alive` + `Sc alive`,  # Total vivant
    `% dead` = (1 - (Time / temps)) * 100  # Exemples de % morts
  )

# Simuler Sc alive
simulate_sc_curve <- function(time, max_value, growth_rate, decline_rate, replicat_id, noise_factor) {
  # Générer des variations spécifiques à chaque réplicat
  noisy_max_value <- max_value * (1 + rnorm(1, mean = 0, sd = noise_factor) * replicat_id)
  noisy_growth_rate <- growth_rate * (1 + rnorm(1, mean = 0, sd = noise_factor) * replicat_id)
  noisy_decline_rate <- decline_rate * (1 + rnorm(1, mean = 0, sd = noise_factor) * replicat_id)
  
  # Courbe combinée : croissance et décroissance légères
  growth <- noisy_max_value * (1 - exp(-noisy_growth_rate * time))  # Phase de croissance
  decline <- exp(-((time - 48)^2) / (2 * noisy_decline_rate^2))  # Phase de décroissance lente
  growth * decline
}

# Mise à jour des données simulées pour inclure `Sc alive`
simulated_data <- simulated_data %>%
  mutate(
    max_value_sc = case_when(
      Proportion == 1 ~ 200000,
      Proportion == 10 ~ 150000,
      Proportion == 100 ~ 120000,
      Proportion == 1000 ~ 90000
    ),
    growth_rate_sc = case_when(
      Proportion == 1 ~ 0.1,
      Proportion == 10 ~ 0.08,
      Proportion == 100 ~ 0.06,
      Proportion == 1000 ~ 0.04
    ),
    decline_rate_sc = case_when(
      Proportion == 1 ~ 400,
      Proportion == 10 ~ 400,
      Proportion == 100 ~ 400,
      Proportion == 1000 ~ 400
    ),
    # Simuler Sc alive avec un bruit unique par réplicat
    `Sc alive` = simulate_sc_curve(
      Time, max_value_sc, growth_rate_sc, decline_rate_sc, Replicat, noise_factor = 0.05
    )
  )

ggplot(simulated_data) +
  # Courbes pour Td alive
  geom_line(aes(x = Time, y = `Td alive`, color = as.factor(Replicat), group = interaction(Proportion, Replicat)), size = 1) +
  # Courbes pour Sc alive
  geom_line(aes(x = Time, y = `Sc alive`, color = as.factor(Replicat), group = interaction(Proportion, Replicat)), size = 1, linetype = "dashed") +
  labs(
    title = "Simulation de Td et Sc alive par modalité et réplicat",
    x = "Temps (h)",
    y = "Nb cellules",
    color = "Réplicat"
  ) +
  facet_wrap(~ Proportion, ncol = 2, scales = "free_y") +  # Une facette par modalité
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))




### INTERPOLATION SPLINE ### ******script used for the article (with some variable before)********  

# Calcul des splines pour chaque Proportion
splines <- simulated_data %>%
  group_by(Proportion) %>%
  summarize(
    spline_td = list(as.data.frame(spline(Time, `Td alive`, n = 100, method = "natural"))),
    spline_sc = list(as.data.frame(spline(Time, `Sc alive`, n = 100, method = "natural")))
  )

# Transformer les splines en format utilisable pour ggplot
splines_td <- splines %>%
  select(Proportion, spline_td) %>%
  unnest(cols = spline_td) %>%
  mutate(Type = "Td alive")

splines_sc <- splines %>%
  select(Proportion, spline_sc) %>%
  unnest(cols = spline_sc) %>%
  mutate(Type = "Sc alive")

## seuil sur spline ## (plus précis)

# Combiner les splines de Td et Sc pour calculer le pourcentage
splines_combined <- splines_td %>%
  inner_join(splines_sc, by = c("x", "Proportion"), suffix = c("_Td", "_Sc")) %>%
  mutate(`% Td` = (y_Td / (y_Td + y_Sc)) * 100)

# Identifier le premier moment où Td alive reste sous les 5 %
first_below_5_percent_spline <- splines_combined %>%
  group_by(Proportion) %>% 
  filter(`% Td` < 5) %>%
  slice_head(n = 1) %>%        # Prendre le premier moment
  ungroup() %>%                # Sortir du regroupement
  rename(Time = x)             # Renommer x en Time pour ggplot

# Afficher franchissement du seuil pour spline global
print(first_below_5_percent_spline)

## test spline sur réplicat pour IC ##

# Étape 1 : Fonction pour générer des splines
generate_splines <- function(data) {
  data %>%
    group_by(Proportion) %>%
    summarize(
      spline_td = list(as.data.frame(spline(Time, `Td alive`, n = 100, method = "natural"))),
      spline_sc = list(as.data.frame(spline(Time, `Sc alive`, n = 100, method = "natural"))),
      .groups = "drop"
    )
}

# Étape 2 : Fonction pour traiter les splines

process_splines <- function(splines) {
  # Transformer les splines
  splines_td <- splines %>%
    select(Proportion, spline_td) %>%
    unnest(cols = spline_td) %>%
    rename(x = x, y_Td = y)
  
  splines_sc <- splines %>%
    select(Proportion, spline_sc) %>%
    unnest(cols = spline_sc) %>%
    rename(x = x, y_Sc = y)
  
  # Combiner les splines et calculer % Td
  splines_combined <- splines_td %>%
    inner_join(splines_sc, by = c("x", "Proportion")) %>%
    mutate(`% Td` = pmax(0, y_Td / (y_Td + y_Sc) * 100))  # Remplacer les négatifs par 0
  
  # Identifier le premier moment où Td passe sous 5 %
  splines_combined %>%
    group_by(Proportion) %>%
    filter(`% Td` < 5) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    rename(Time = x)
}

# Étape 3 : Boucle pour traiter chaque réplicat
results <- lapply(1:3, function(rep) {
  # Filtrer les données pour le réplicat
  data_replicat <- simulated_data %>% filter(Replicat == rep)
  
  # Générer et traiter les splines
  splines <- generate_splines(data_replicat)
  process_splines(splines)
})

# Combiner les résultats des réplicats dans un seul tableau
final_results <- bind_rows(results, .id = "Replicat")

# Prendre le min et max des seuils des réplicats
final_results <- final_results %>%
  group_by(Proportion) %>%  # Grouper par Proportion
  filter(`Time` == min(`Time`) | `Time` == max(`Time`)) %>%  # Garder min et max pour chaque groupe
  ungroup()  # Sortir du regroupement

# Afficher franchissements du seuil pour splines réplicats
print(final_results[order(final_results$Proportion), ])

# Afficher graphe
ggplot() +
  # Courbes spline pour Td et Sc
  geom_line(data = splines_td, aes(x = x, y = y, color = "Spline Td")) +
  geom_line(data = splines_sc, aes(x = x, y = y, color = "Spline Sc")) +
  # Lignes verticales pour les moments identifiés (tous les réplicats)
  geom_vline(data = final_results, aes(xintercept = Time, color = as.factor(Replicat)),
             linetype = "dashed") +
  # Moyenne des franchissements (par modalité)
  geom_vline(data = first_below_5_percent_spline, aes(xintercept = Time, group = Proportion), 
             linetype = "solid", color = "darkorange") +
  # Titre, labels et ajustements de légende
  labs(title = "Interpolation spline pour Sc et Td simulés par Proportion\nAvec seuil de 5% et IC",
       x = "Temps (h)", y = "Valeur",
       color = "Type", fill = "IC spline") +
  scale_color_manual(values = c("Points Td" = "blue", "Spline Td" = "blue", 
                                "Points Sc" = "red", "Spline Sc" = "red")) +
  facet_wrap(~ Proportion, ncol = 1, scales = "free_y") +  # Facetting par Proportion
  # scale_x_continuous(limits = c(110, 196)) +  # Limites de l'axe X
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



