# ================================================================
# 📦 Synthetic Healthcare Dataset (200 Patients)
# With ~20–30% missingness in key variables
# ================================================================

library(dplyr)

set.seed(123)
n <- 200

# --- Base demographics ---
data <- data.frame(
  ID = sprintf("P%03d", 1:n),
  SEX = sample(c("Male", "Female"), n, replace = TRUE, prob = c(0.55, 0.45)),
  AGE = round(rnorm(n, mean = 45, sd = 12)),
  WT  = round(rnorm(n, mean = 70, sd = 12), 1),
  HT  = round(rnorm(n, mean = 170, sd = 10), 1)
)

# --- Derived BMI ---
data$BMI <- round(data$WT / ((data$HT/100)^2), 1)

# --- Dosing information ---
data$DOSE <- sample(c(50, 100, 150, 200, 250, 300), n, replace = TRUE)
data$ROUTE <- sample(c("Oral", "IV"), n, replace = TRUE, prob = c(0.8, 0.2))
data$FREQ <- sample(c("Once daily", "Twice daily"), n, replace = TRUE)
data$TIME_HR <- sample(seq(0.5, 12, 0.5), n, replace = TRUE)

# --- Pharmacokinetics ---
data$CMAX <- round(rnorm(n, mean = 2.5, sd = 0.6) * (data$DOSE/100) / (data$WT/70), 2)
data$TMAX <- round(runif(n, 0.5, 4.5), 2)
data$AUC <- round(rnorm(n, mean = 15, sd = 4) * (data$DOSE/100) / (data$WT/70), 2)
data$CL <- round(runif(n, 3, 9), 2)
data$T_HALF <- round(runif(n, 2, 8), 2)

# --- Clinical measurements ---
data$BP_SYS <- round(rnorm(n, mean = 120, sd = 12))
data$BP_DIA <- round(rnorm(n, mean = 80, sd = 8))
data$HR <- round(rnorm(n, mean = 75, sd = 10))
data$ALT <- round(rnorm(n, mean = 35, sd = 10), 1)
data$CREAT <- round(rnorm(n, mean = 1.0, sd = 0.2), 2)
data$GLUCOSE <- round(rnorm(n, mean = 95, sd = 15), 1)

# --- Response variable ---
data <- data %>%
  mutate(RESPONSE = ifelse(AUC > median(AUC) & BP_SYS < 130, "Responder", "Nonresponder"))

# --- Introduce Missingness (20–30%) ---
introduce_missing <- function(x, prop) {
  idx <- sample(1:length(x), size = floor(prop * length(x)))
  x[idx] <- NA
  return(x)
}

data$CMAX <- introduce_missing(data$CMAX, runif(1, 0.2, 0.3))
data$AUC <- introduce_missing(data$AUC, runif(1, 0.2, 0.3))
data$RESPONSE <- introduce_missing(data$RESPONSE, runif(1, 0.2, 0.3))
data$ALT <- introduce_missing(data$ALT, runif(1, 0.2, 0.3))
data$CREAT <- introduce_missing(data$CREAT, runif(1, 0.2, 0.3))
data$RESPONSE <-  as.factor(data$RESPONSE)
data$SEX <-  as.factor(data$SEX)
# --- Save to CSV ---
write.csv(data, "healthcare_dataset_with_missing.csv", row.names = FALSE)

cat("✅ 'healthcare_dataset_with_missing.csv' created successfully!\n")

