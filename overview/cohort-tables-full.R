data.filename <- '../../data/cohort-sanitised.csv'
require(data.table)
COHORT <- fread(data.filename)

percentMissing <- function(x, sf = 3) {
  round(sum(is.na(x))/length(x), digits = sf)*100
}

# Remove the patients we shouldn't include
COHORT <-
  COHORT[
    # remove negative times to death
    COHORT$time_death > 0 &
      # remove patients who should be excluded
      !COHORT$exclude
    ,
    ]

# Age, 5, 50, 95, %missing
print(quantile(COHORT$age, c(0.5, 0.025, 0.975)))

# Gender
print(table(COHORT$gender))
print(table(COHORT$gender)/nrow(COHORT)*100)

# Deprivation, 5, 50, 95, %missing
print(quantile(COHORT$imd_score, c(0.5, 0.025, 0.975), na.rm = TRUE))
print(percentMissing(COHORT$imd_score))

# SCAD subtype
print(table(COHORT$diagnosis)/nrow(COHORT)*100)

# PCI
print(sum(COHORT$pci_6mo)/nrow(COHORT)*100)

# CABG
print(sum(COHORT$cabg_6mo)/nrow(COHORT)*100)

# previous/recurrent MI
print(sum(COHORT$hx_mi)/nrow(COHORT)*100)

# nitrates (listed as 1 and NA not T and F)
print(sum(COHORT$long_nitrate, na.rm = TRUE)/nrow(COHORT)*100)

# Smoking, by category, %missing
print(table(COHORT$smokstatus)/nrow(COHORT)*100)
print(percentMissing(COHORT$smokstatus))

# Hypertension
print(sum(COHORT$hypertension)/nrow(COHORT)*100)

# Diabetes, yes/no
print(
  (sum(COHORT$diabetes == 'Diabetes unspecified type') +
    sum(COHORT$diabetes == 'Type 1 diabetes') +
    sum(COHORT$diabetes == 'Type 2 diabetes')) /nrow(COHORT)*100
)

# Total cholesterol
print(quantile(COHORT$total_chol_6mo, c(0.5, 0.025, 0.975), na.rm = TRUE))
print(percentMissing(COHORT$total_chol_6mo))

# HDL
print(quantile(COHORT$hdl_6mo, c(0.5, 0.025, 0.975), na.rm = TRUE))
print(percentMissing(COHORT$hdl_6mo))

# Heart failure
print(sum(COHORT$heart_failure)/nrow(COHORT)*100)

# Peripheral arterial disease
print(sum(COHORT$pad)/nrow(COHORT)*100)

# Atrial fibrillation
print(sum(COHORT$hx_af)/nrow(COHORT)*100)

# Stroke
print(sum(COHORT$hx_stroke)/nrow(COHORT)*100)

# Chronic kidney disease
print(sum(COHORT$hx_renal)/nrow(COHORT)*100)

# COPD
print(sum(COHORT$hx_copd)/nrow(COHORT)*100)

# Cancer
print(sum(COHORT$hx_cancer)/nrow(COHORT)*100)

# Chronic liver disease
print(sum(COHORT$hx_liver)/nrow(COHORT)*100)

# Depression
print(sum(COHORT$hx_depression)/nrow(COHORT)*100)

# Anxiety
print(sum(COHORT$hx_anxiety)/nrow(COHORT)*100)

# Heart rate
print(quantile(COHORT$pulse_6mo, c(0.5, 0.025, 0.975), na.rm = TRUE))
print(percentMissing(COHORT$pulse_6mo))

# Creatinine
print(quantile(COHORT$crea_6mo, c(0.5, 0.025, 0.975), na.rm = TRUE))
print(percentMissing(COHORT$crea_6mo))

# WCC
print(quantile(COHORT$total_wbc_6mo, c(0.5, 0.025, 0.975), na.rm = TRUE))
print(percentMissing(COHORT$total_wbc_6mo))

# Haemoglobin
print(quantile(COHORT$haemoglobin_6mo, c(0.5, 0.025, 0.975), na.rm = TRUE))
print(percentMissing(COHORT$haemoglobin_6mo))

# Follow-up, 5, 50, 95
print(quantile(COHORT$endpoint_death_date, c(0.5, 0.025, 0.975)))/365.25

# Death vs censored, %
print(table(COHORT$endpoint_death)) /nrow(COHORT)*100
