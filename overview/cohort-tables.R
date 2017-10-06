data.filename <- '../../data/cohort-sanitised.csv'

source('../lib/shared.R')

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- fread(data.filename)

# Remove the patients we shouldn't include
COHORT.full <-
  COHORT.full[
    # remove negative times to death
    COHORT.full$time_death > 0 &
    # remove patients who should be excluded
    !COHORT.full$exclude
    ,
    ]

# Age, 5, 50, 95, %missing
print(quantile(COHORT.full$age, c(0.5, 0.05, 0.95)))

# Gender
print(table(COHORT.full$gender))/nrow(COHORT.full)*100

# Deprivation, 5, 50, 95, %missing
print(quantile(COHORT.full$imd_score, c(0.5, 0.05, 0.95), na.rm = TRUE))
print(percentMissing(COHORT.full$imd_score))

# Smoking, by category, %missing
print(table(COHORT.full$smokstatus))/nrow(COHORT.full)*100
print(percentMissing(COHORT.full$smokstatus))

# Diabetes, yes/no
print(
  ( sum(COHORT.full$diabetes == 'Diabetes unspecified type') +
    sum(COHORT.full$diabetes == 'Type 1 diabetes') +
    sum(COHORT.full$diabetes == 'Type 2 diabetes')) /nrow(COHORT.full)*100
)

# Follow-up, 5, 50, 95
print(quantile(COHORT.full$endpoint_death_date, c(0.5, 0.05, 0.95)))/365.25

# Death vs censored, %
print(table(COHORT.full$endpoint_death)) /nrow(COHORT.full)*100
