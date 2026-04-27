library(tidyverse)
library(labelled)

# ==============================================================================
# Color Themes ----
# ==============================================================================
custom_palettes <- list(
  maritime = c(
    "indigo" = "#102542",
    "terracotta" = "#C5533B",
    "cerulean" = "#2E86AB",
    "ochre" = "#D8973C",
    "slate" = "#6B818C",
    "kelp" = "#4A6660", # Neu: Ein gedecktes Meeresgrün
    "sand" = "#E1D5BA" # Neu: Ein heller Sandton
  ),
  botanical = c(
    "forest" = "#1B4332",
    "berry" = "#A63A50",
    "moss" = "#6A994E",
    "slate" = "#4A4E69",
    "tan" = "#D4A373",
    "sage" = "#A3B19B", # Neu: Ein helles Salbeigrün
    "rose" = "#D28F9D" # Neu: Ein weiches Altrosa
  ),
  nordic = c(
    "petrol" = "#0F4C5C",
    "crimson" = "#9A031E",
    "spruce" = "#4C7273",
    "caramel" = "#CB793A",
    "steel" = "#8DA9C4",
    "pine" = "#2C3D3B", # Neu: Ein sehr dunkles Nadelgrün
    "linen" = "#E8E4D9" # Neu: Ein kühles Leinen-Weiß
  ),
  azure_coral = c(
    "azure" = "#1D70A2",
    "coral" = "#EF6461",
    "jade" = "#3BB273",
    "mustard" = "#E1BC29",
    "steel" = "#A0AAB2",
    "navy" = "#153B50", # Neu: Ein tiefes Dunkelblau
    "peach" = "#F6C1A4" # Neu: Ein heller Pfirsichton
  ),
  eucalyptus = c(
    "eucalyptus" = "#238B86",
    "raspberry" = "#D1495B",
    "celadon" = "#8DB596",
    "apricot" = "#F3A78B",
    "slate" = "#6F7A89",
    "fir" = "#1A4A47", # Neu: Ein dunkles Tannengrün
    "cream" = "#F2E8CF" # Neu: Ein weicher Cremeton
  ),
  aegean = c(
    "ocean" = "#3D5A80",
    "terracotta" = "#EE6C4D",
    "seafoam" = "#81B29A",
    "sand" = "#F2CC8F",
    "lilac" = "#9A8C98",
    "deep_sea" = "#293241", # Neu: Ein sehr dunkles Ozeanblau
    "cloud" = "#E0E1DD" # Neu: Ein helles Wolkengrau
  ),
  oasis = c(
    "teal" = "#2A9D8F",
    "grapefruit" = "#E76F51",
    "olive" = "#8AB17D",
    "maize" = "#E9C46A",
    "glaucous" = "#8D99AE",
    "midnight" = "#264653", # Neu: Ein fast schwarzes Blaugrün
    "peach" = "#F4A261" # Neu: Ein weiches, helles Orange
  ),
  iris_copper = c(
    "iris" = "#5B6C9C",
    "copper" = "#CB6D51",
    "aqua" = "#54B0A7",
    "ochre" = "#E5B25D",
    "warm_grey" = "#9C9592",
    "aubergine" = "#4A3B4C", # Neu: Ein dunkler Aubergineton
    "sand" = "#E6D5B8" # Neu: Ein warmer, heller Sandton
  )
)

usethis::use_data(custom_palettes, overwrite = TRUE)

# ==============================================================================
# Thesis ----
# ==============================================================================
set.seed(1775575805)

n_per_group <- 10000

before_words <- rnorm(n_per_group, mean = 600, sd = 150)
individual_effect <- rnorm(n_per_group, mean = 20, sd = 50)

after_words <- before_words + individual_effect

thesis <- tibble(
  id = rep(1:n_per_group, 2),
  time = factor(
    rep(c("Before", "After"), each = n_per_group),
    levels = c("Before", "After")
  ),
  words_written = c(before_words, after_words)
) |>
  set_variable_labels(
    id = "ID",
    time = "Measurement",
    words_written = "Words Written"
  )

usethis::use_data(thesis, overwrite = TRUE)

# ==============================================================================
# Imposter ----
# ==============================================================================
set.seed(1775576764)

n <- 50
pre_mean <- 70
mean_change <- -15

imposter <- tibble(
  id = seq(1, n),
  pre = rnorm(n, mean = pre_mean, sd = 10),
  change = rnorm(n, mean = mean_change, sd = 12)
) |>
  mutate(
    post = pre + change
  ) |>
  select(-change) |>
  pivot_longer(
    cols = c(pre, post),
    values_to = "imposter_score",
    names_to = "measurement"
  ) |>
  mutate(
    imposter_score = pmin(pmax(imposter_score, 0), 100)
  ) |>
  set_variable_labels(
    id = "ID",
    measurement = "Measurement",
    imposter_score = "Imposter Score"
  )

usethis::use_data(imposter, overwrite = TRUE)

# ==============================================================================
# Shoulder ----
# ==============================================================================
set.seed(1775578562)
n <- 50

# Wir simulieren realistische "Frozen Shoulder" Werte in Grad
pre_rom <- rnorm(n, mean = 90, sd = 15)
change <- rnorm(n, mean = 20, sd = 12)
post_rom <- pre_rom + change

shoulder <- tibble(
  id = 1:n,
  Pre = pre_rom,
  Post = post_rom
) |>
  pivot_longer(
    cols = -id,
    names_to = "time",
    values_to = "rotation"
  ) |>
  mutate(
    time = as_factor(time)
  ) |>
  set_variable_labels(
    id = "ID",
    time = "Measurement",
    rotation = "Range of Motion (Degrees)"
  )

usethis::use_data(shoulder, overwrite = TRUE)


set.seed(1775648491)

# ==============================================================================
# PHQ-9 ----
# ==============================================================================
set.seed(1777294530)
N <- 150 # Anzahl Probanden
I <- 9 # Anzahl Items

# Wahre Prä-Werte: z.B. hohe Symptombelastung zu Beginn
theta_pre <- rnorm(N, mean = 2.1, sd = 1)

# Wahre Post-Werte: Intervention senkt den Wert im Schnitt um 1.0, mit individueller Varianz (Schwankung)
theta_post <- theta_pre + rnorm(N, mean = -2, sd = 1.4)

# Schwierigkeiten der Items
item_diffs <- seq(-2, 2, length.out = I)

# Hilfsfunktion zur Simulation einer 4-stufigen Likert-Skala (0 bis 3)
sim_likert <- function(theta, b) {
  n_persons <- length(theta)
  n_items <- length(b)
  dat <- matrix(0, nrow = n_persons, ncol = n_items)

  for (i in 1:n_items) {
    latent_response <- (theta - b[i]) + rlogis(n_persons)

    dat[, i] <- case_when(
      latent_response > 1.5 ~ 3, # "Fast jeden Tag"
      latent_response > 0.0 ~ 2, # "An mehr als der Hälfte der Tage"
      latent_response > -1.5 ~ 1, # "An einzelnen Tagen"
      TRUE ~ 0 # "Überhaupt nicht"
    )
  }
  dat
}

# 4-stufige Likert-Antworten für Prä und Post generieren
dat_pre <- sim_likert(theta_pre, item_diffs)
dat_post <- sim_likert(theta_post, item_diffs)

colnames(dat_pre) <- colnames(dat_post) <- paste0("item_", 1:I)

# --- 2. Concurrent Calibration (IRT) ---
dat_combined <- rbind(dat_pre, dat_post)

phq_9 <- dat_combined |>
  as.data.frame() |>
  as_tibble() |>
  mutate(
    id = rep(1:N, 2),
    time = factor(rep(c("pre", "post"), each = N), levels = c("pre", "post")),
    # Summe über alle 9 Items bilden (Maximalwert = 27)
    phq_9 = rowSums(across(starts_with("Item")))
  ) |>
  relocate(id, time, .before = 1) |>
  arrange(id) |>
  set_variable_labels(
    id = "ID",
    time = "Meassurement",
    phq_9 = "PHQ-9 Score"
  )

usethis::use_data(phq_9, overwrite = TRUE)

# ==============================================================================
# Burnout ----
# ==============================================================================
set.seed(1776172078)

n_per_group <- 30

int_data <- map_df(1:n_per_group, \(a) {
  max_weeks <- sample(2:6, 1)
  tibble(id = a, group = "Intervention", time_weeks = 0:max_weeks)
})

wl_data <- map_df((n_per_group + 1):(2 * n_per_group), \(a) {
  tibble(id = a, group = "Waitlist", time_weeks = seq(0, 12, by = 4))
})

burnout <- bind_rows(int_data, wl_data) |>
  mutate(group = factor(group, levels = c("Waitlist", "Intervention")))

# Wahre Effekte simulieren
# Startwert (Intercept) bei ca. 75 (Hohes Burnout)
# Waitlist Slope: -1.0 Punkte pro Woche
# Intervention Slope: -3.5 Punkte pro Woche
burnout <- burnout |>
  mutate(
    # Random Intercept (Jeder startet anders)
    ind_intercept = rnorm(1, mean = 75, sd = 8),
    # Random Slope (Jeder reagiert anders auf die Zeit)
    ind_slope_base = rnorm(1, mean = -1.0, sd = 0.5),
    # Zusätzlicher Interventionseffekt (nur für Int-Gruppe)
    ind_int_effect = if_else(
      group == "Intervention",
      rnorm(1, mean = -2.5, sd = 1.0),
      0
    ),

    # EEI berechnen: Startwert + (Basis-Slope + Int-Slope) * Zeit + Messfehler
    eei_true = ind_intercept + (ind_slope_base + ind_int_effect) * time_weeks,
    eei = eei_true + rnorm(n(), mean = 0, sd = 3), # Messrauschen
    .by = id
  ) |>
  mutate(eei = pmin(pmax(eei, 0), 100)) |>
  select(id, group, time_weeks, eei) |>
  set_variable_labels(
    id = "ID",
    group = "Intervention Group",
    time_weeks = "Weeks Since Baseline",
    eei = "Emotional Exhaustion Inventory Score"
  )

usethis::use_data(burnout, overwrite = TRUE)


# ==============================================================================
# Espresso ----
# ==============================================================================
set.seed(1775664741)

n_patients <- 50

# 1. Patienten-Profile: Realistische Varianzen
patients <- tibble(
  id = 1:n_patients,
  age = round(runif(n_patients, 20, 60)),
  espresso = sample(0:4, n_patients, replace = TRUE),

  # Solide, aber keine absurden Unterschiede beim Start-Level (Baseline)
  u_0i = rnorm(n_patients, mean = 0, sd = 8),

  # Ein kleines bisschen individuelle Biologie bleibt erhalten (nicht jeder
  # Null-Espresso-Trinker ist ein exakter Klon).
  u_1i = rnorm(n_patients, mean = 0, sd = 1)
)

# 2. Datengenerierung: Klarer Kausalpfad, aber mit klinischem Rauschen
espresso <- patients |>
  expand_grid(week = c(0, 1, 2, 3)) |> # Deine 4 Messzeitpunkte!
  mutate(
    is_post = ifelse(week >= 2, 1, 0),

    # Baseline wird etwas durchs Alter diktiert
    baseline_energy = 60 - (0.5 * age) - (1.5 * espresso) + u_0i,

    # MODERATE INTERAKTION:
    # 0 Espressos = +12 Punkte (Klarer Erfolg)
    # 2 Espressos = +6 Punkte  (Knapp über Signifikanz)
    # 4 Espressos =  0 Punkte  (Kein Effekt)
    true_effect = 12 - (3 * espresso) + u_1i,

    # DAS GEHEIMNIS FÜR DAS BUCH: Realistisches Messrauschen (sd = 4)!
    # Das Rauschen sorgt dafür, dass das Rohmodell verwirrt ist.
    # Das Kontext-Modell hingegen nutzt "Smart Shrinkage", um das Rauschen
    # durch das Wissen über die Subgruppe auszuhebeln.
    morning_energy = baseline_energy +
      (is_post * true_effect) +
      rnorm(n(), 0, 4)
  ) |>
  select(id, age, espresso, week, is_post, morning_energy) |>
  set_variable_labels(
    id = "ID",
    age = "Age",
    espresso = "Daily Espressos",
    week = "Measurement Week",
    is_post = "Intervention Phase",
    morning_energy = "Energy Score"
  )

usethis::use_data(espresso, overwrite = TRUE)


# ==============================================================================
# Misophonia ----
# ==============================================================================
set.seed(1775673111)

n_days <- 30
intervention_start <- 11

true_intercept <- 25
true_slope_A <- 1.5
true_slope_B <- -2.0
ar_rho <- 0.3
sigma <- 4

# AR(1) Error simulation
errors <- numeric(n_days)
errors[1] <- rnorm(1, 0, sigma)
for (i in 2:n_days) {
  errors[i] <- ar_rho * errors[i - 1] + rnorm(1, 0, sigma)
}

misophonia <- tibble::tibble(day = 1:n_days) |>
  dplyr::mutate(
    phase = ifelse(day >= intervention_start, "Therapy", "Baseline"),
    phase_1 = ifelse(day < intervention_start, day, intervention_start - 1),
    phase_2 = ifelse(
      day >= intervention_start,
      day - intervention_start + 1,
      0
    ),
    severity_score = true_intercept +
      (true_slope_A * phase_1) +
      (true_slope_B * phase_2) +
      errors,

    # Cap values between 0 and 60
    severity_score = pmin(60, pmax(0, severity_score))
  ) |>
  # Keep only the relevant columns for the final dataset
  dplyr::select(day, phase, severity_score) |>
  set_variable_labels(
    day = "Day of Study",
    phase = "Study Phase",
    severity_score = "Daily Distress Level"
  )

usethis::use_data(misophonia, overwrite = TRUE)


# ==============================================================================
# Man on Fire ----
# ==============================================================================
set.seed(1775673111)

n_days_per_phase <- 10

# The biological truth (The levels of pain in each phase)
true_A1 <- 80 # High pain
true_B1 <- 35 # Drug works!
true_A2 <- 85 # Pain returns (Withdrawal with slight rebound)
true_B2 <- 30 # Drug works again!

ar_rho <- 0.4 # Moderate carry-over from yesterday's pain
sigma <- 8 # Daily fluctuation

man_on_fire <- tibble::tibble(
  day = 1:40,
  # Create the labels A1, B1, A2, B2
  phase = factor(
    rep(c("A1", "B1", "A2", "B2"), each = n_days_per_phase),
    levels = c("A1", "B1", "A2", "B2")
  )
) |>
  dplyr::mutate(
    # Set the true mean depending on the phase
    true_mean = dplyr::case_when(
      phase == "A1" ~ true_A1,
      phase == "B1" ~ true_B1,
      phase == "A2" ~ true_A2,
      phase == "B2" ~ true_B2
    )
  )

# Simulate AR(1) errors
errors <- numeric(40)
errors[1] <- rnorm(1, 0, sigma)
for (i in 2:40) {
  errors[i] <- ar_rho * errors[i - 1] + rnorm(1, 0, sigma)
}

# Final pain score
man_on_fire <- man_on_fire |>
  dplyr::mutate(
    pain_score = true_mean + errors,
    pain_score = pmin(100, pmax(0, pain_score)) # Keep within 0-100 bounds
  ) |>
  dplyr::select(day, phase, pain_score) |> # Clean up for the final dataset
  set_variable_labels(
    day = "Day of Study",
    phase = "Study Phase",
    pain_score = "Pain Severity"
  )

usethis::use_data(man_on_fire, overwrite = TRUE)


# ==============================================================================
# Return to Work ----
# ==============================================================================
set.seed(1775731636)

n_employees <- 50
total_months <- 6

employees <- tibble(
  id = 1:n_employees,
  # Baseline variability: Different starting levels of sickness
  individual_intercept = rnorm(n_employees, 0, 0.15),
  # Slope variability: Some recover fast, some stagnate
  individual_slope = rnorm(n_employees, 0, 0.20)
)

return_to_work <- employees |>
  expand_grid(month = 1:total_months) |>
  mutate(
    # TRICK: We center time at 0 so the intercept represents Month 1
    time = month - 1,

    # THE LOG SCALE MATH:
    # Start at ~12 days sick (log(12)).
    # Drop by ~20% every month (log(0.8)).
    log_mu = log(12) +
      individual_intercept +
      (log(0.80) + individual_slope) * time,

    # Translate back to the original scale (expected days)
    mu = exp(log_mu),

    # Generate actual counts with Overdispersion (Negative Binomial)
    sick_days = rnbinom(n(), mu = mu, size = 2),
    sick_days = pmin(sick_days, 31)
  ) |>
  select(id, month, time, sick_days) |>
  set_variable_labels(
    id = "ID",
    month = "Month of Observation",
    time = "Time",
    sick_days = "Number of Sick Days per Month"
  )

usethis::use_data(return_to_work, overwrite = TRUE)


# ==============================================================================
# Athletes ----
# ==============================================================================
#%%
set.seed(1775757656)

n_athletes <- 60

individual_athletes <- faux::rnorm_multi(
  n = n_athletes,
  vars = 2,
  varnames = c("strength_gain", "anxiety_drop"),
  mu = c(30, -20),
  # REDUZIERTE VARIANZ 1 (Wahres Therapie-Ansprechen):
  # Vorher SD = 10 und 8. Jetzt halbieren wir das auf 5 und 4.
  sd = c(18, 12),
  # Die Korrelation bleibt bei starken -0.75!
  r = -0.75
) |>
  dplyr::mutate(
    id = 1:dplyr::n(),
    strength_base = rnorm(dplyr::n(), 180, 10),
    anxiety_base = rnorm(dplyr::n(), 60, 5)
  )

# Wir messen 3x vor dem Camp und 3x danach!
athletes <- individual_athletes |>
  tidyr::expand_grid(
    phase = c(0, 1),
    measurement = 1:3
  ) |>
  dplyr::mutate(
    strength = strength_base +
      (strength_gain * phase) +
      rnorm(dplyr::n(), 0, 8),
    anxiety = anxiety_base + (anxiety_drop * phase) + rnorm(dplyr::n(), 0, 4)
  ) |>
  # Nur das ausgeben, was der Forscher in der Realität auch messen kann!
  dplyr::select(id, phase, measurement, strength, anxiety) |>
  set_variable_labels(
    id = "ID",
    phase = "Study Phase",
    measurement = "Measurement",
    strength = "Physical Strength",
    anxiety = "Anxiety"
  )

usethis::use_data(athletes, overwrite = TRUE)

# ==============================================================================
# Post-Cardiac Rehabilitation ----
# ==============================================================================
set.seed(1775825526)
n_patients <- 90

individual_params <- tibble(
  id = 1:n_patients,
  efficacy = rnorm(n_patients, mean = 50, sd = 15), # Kardiale Selbstwirksamkeit
  intercept = rnorm(n_patients, mean = 45, sd = 10), # Startwert Aktivität
  slope = rnorm(n_patients, mean = -1.5, sd = 2.5) # Steigung (die Trompete)
)

cardiac_rehab <- individual_params |>
  expand_grid(week = 1:12) |>
  mutate(
    active_minutes = intercept +
      (slope * week) +
      (efficacy * 0.05 * week) +
      rnorm(n(), mean = 0, sd = 8),

    active_minutes = pmax(0, active_minutes),

    # MAR-Dropout: Ab Woche 5 brechen Patienten mit niedriger Efficacy oft ab
    is_missing = ifelse(week >= 5 & efficacy < 45, rbinom(n(), 1, 0.8), 0),
    active_minutes = ifelse(is_missing == 1, NA, active_minutes)
  ) |>
  select(id, week, efficacy, active_minutes) |>
  set_variable_labels(
    id = "ID",
    week = "Week of Rehabilitation",
    efficacy = "Cardiac Self-Efficacy",
    active_minutes = "Physical Activity (Minutes)"
  )

usethis::use_data(cardiac_rehab, overwrite = TRUE)
