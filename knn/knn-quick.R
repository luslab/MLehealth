# A quick and dirty script to try the k nearest neighbours approach for survival
# as in the bnnSurvival package. C-indices were about 0.75 without any real
# optimisation, which isn't that impressive, even when I used the bagging which
# is the b in bnnSurvival... What I want to do is this, but learning the weight
# to apply to each component when calculating distances, but no packages seem to
# offer that. The other question is, when your approach becomes more
# sophisticated, is there ever a theoretical reason to expect this to do better
# than random forests, which are effectively designed to partition the
# variable-space in exactly the way this algorithm would?

n.data <- NA

data.filename <- '../../data/cohort-sanitised.csv'
output.base <- '../../output/rf-continuous-try5'

endpoint <- 'death' # Change to MI to look for MI...anything else uses death

n.trees <- 500
n.threads <- 16
# What to do with missing data
continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )
untransformed.vars <- c('anonpatid', 'surv_time', 'imd_score', 'exclude')

source('../lib/shared.R')
requirePlus('bnnSurvival')


# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# Remove the patients we shouldn't include
COHORT.full <-
  COHORT.full[
    # remove negative times to death
    COHORT.full$time_death > 0 &
      # remove patients who should be excluded
      !COHORT.full$exclude
    ,
    ]

# If n.data was specified...
if(!is.na(n.data)){
  # Take a subset n.data in size
  COHORT.use <- sample.df(COHORT.full, n.data)
  rm(COHORT.full)
} else {
  # Use all the data
  COHORT.use <- COHORT.full
  rm(COHORT.full)
}

# Redefine n.data just in case we lost any rows
n.data <- nrow(COHORT.use)
# Define indices of test set
test.set <- sample(1:n.data, (1/3)*n.data)

COHORT.scaled <-
  data.frame(
    ## Time to event
    surv_time = COHORT.use[, surv.time],
    ## Death/censorship
    surv_event = COHORT.use[, surv.event] %in% surv.event.yes,
    ## Rescaled age
    age = COHORT.use$age,
    ## Gender
    gender = COHORT.use$gender,
    ## Most deprived quintile, yes vs. no
    most_deprived =
      COHORT.use$imd_score > quantile(COHORT.use$imd_score, 0.8, na.rm = TRUE),
    ### SCAD diagnosis and severity ############################################
    ## Other CHD / unstable angina / NSTEMI / STEMI vs. stable angina
    diagnosis = factorChooseFirst(factor(COHORT.use$diagnosis), 'SA'),
    ## PCI in last 6 months, yes vs. no
    pci_6mo = COHORT.use$pci_6mo,
    ## CABG in last 6 months, yes vs. no
    cabg_6mo = COHORT.use$cabg_6mo,
    ## Previous/recurrent MI, yes vs. no
    hx_mi = COHORT.use$hx_mi,
    ## Use of nitrates, yes vs. no
    long_nitrate = COHORT.use$long_nitrate,
    ### CVD risk factors #######################################################
    ## Ex-smoker vs. never / Current smoker vs. never
    smokstatus = factorChooseFirst(factor(COHORT.use$smokstatus), 'Non'),
    ## Hypertension, present vs. absent
    hypertension = COHORT.use$hypertension,
    ## Diabetes mellitus, present vs. absent
    diabetes_logical = COHORT.use$diabetes != 'No diabetes',
    ## Total cholesterol, per 1 mmol/L increase
    total_chol_6mo = (COHORT.use$total_chol_6mo - 5),
    ## HDL, per 0.5 mmol/L increase
    hdl_6mo = (COHORT.use$hdl_6mo - 1.5) / 0.5,
    ### CVD co-morbidities #####################################################
    ## Heart failure, present vs. absent
    heart_failure = COHORT.use$heart_failure,
    ## Peripheral arterial disease, present vs. absent
    pad = COHORT.use$pad,
    ## Atrial fibrillation, present vs. absent
    hx_af = COHORT.use$hx_af,
    ## Stroke, present vs. absent
    hx_stroke = COHORT.use$hx_stroke,
    ### Non-CVD comorbidities ##################################################
    ## Chronic kidney disease, present vs. absent
    hx_renal = COHORT.use$hx_renal,
    ## Chronic obstructive pulmonary disease, present vs. absent
    hx_copd = COHORT.use$hx_copd,
    ## Cancer, present vs. absent
    hx_cancer = COHORT.use$hx_cancer,
    ## Chronic liver disease, present vs. absent
    hx_liver = COHORT.use$hx_liver,
    ### Psychosocial characteristics ###########################################
    ## Depression at diagnosis, present vs. absent
    hx_depression = COHORT.use$hx_depression,
    ## Anxiety at diagnosis, present vs. absent
    hx_anxiety = COHORT.use$hx_anxiety,
    ### Biomarkers #############################################################
    ## Heart rate, per 10 b.p.m increase
    pulse_6mo = (COHORT.use$pulse_6mo - 70) / 10,
    ## Creatinine, per 30 Î¼mol/L increase
    crea_6mo = (COHORT.use$crea_6mo - 60) / 30,
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo = (COHORT.use$total_wbc_6mo - 7.5) / 1.5,
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo = (COHORT.use$haemoglobin_6mo - 13.5) / 1.5
  )

missing.cols <-
  c(
    "diagnosis", "most_deprived", "smokstatus", "total_chol_6mo", "hdl_6mo",
    "pulse_6mo", "crea_6mo", "total_wbc_6mo", "haemoglobin_6mo"
  )
COHORT.demissed <- medianImpute(COHORT.scaled, missing.cols)

COHORT.demissed$surv_time_round <-
  round_any(COHORT.demissed$surv_time, 0.1)

surv.formula <-
  Surv(surv_time_round, surv_event) ~
  ### Sociodemographic characteristics #######################################
## Age in men, per year
## Age in women, per year
## Women vs. men
# ie include interaction between age and gender!
age+gender +
  ## Most deprived quintile, yes vs. no
  most_deprived +
  ### SCAD diagnosis and severity ############################################
## Other CHD vs. stable angina
## Unstable angina vs. stable angina
## NSTEMI vs. stable angina
## STEMI vs. stable angina
diagnosis +
  #diagnosis_missing +
  ## PCI in last 6 months, yes vs. no
  pci_6mo +
  ## CABG in last 6 months, yes vs. no
  cabg_6mo +
  ## Previous/recurrent MI, yes vs. no
  hx_mi +
  ## Use of nitrates, yes vs. no
  long_nitrate +
  ### CVD risk factors #######################################################
## Ex-smoker / current smoker / missing data vs. never
smokstatus +
  ## Hypertension, present vs. absent
  hypertension +
  ## Diabetes mellitus, present vs. absent
  diabetes_logical +
  ## Total cholesterol, per 1 mmol/L increase
  total_chol_6mo +
  ## HDL, per 0.5 mmol/L increase
  hdl_6mo +
  ### CVD co-morbidities #####################################################
## Heart failure, present vs. absent
heart_failure +
  ## Peripheral arterial disease, present vs. absent
  pad +
  ## Atrial fibrillation, present vs. absent
  hx_af +
  ## Stroke, present vs. absent
  hx_stroke +
  ### Non-CVD comorbidities ##################################################
## Chronic kidney disease, present vs. absent
hx_renal +
  ## Chronic obstructive pulmonary disease, present vs. absent
  hx_copd +
  ## Cancer, present vs. absent
  hx_cancer +
  ## Chronic liver disease, present vs. absent
  hx_liver +
  ### Psychosocial characteristics ###########################################
## Depression at diagnosis, present vs. absent
hx_depression +
  ## Anxiety at diagnosis, present vs. absent
  hx_anxiety +
  ### Biomarkers #############################################################
## Heart rate, per 10 b.p.m increase
pulse_6mo +
  ## Creatinine, per 30 Î¼mol/L increase
  crea_6mo +
  ## White cell count, per 1.5 109/L increase
  total_wbc_6mo +
  ## Haemoglobin, per 1.5 g/dL increase
  haemoglobin_6mo

surv.model.fit <-
  bnnSurvival(
    surv.formula, COHORT.demissed[-test.set, ],
    k = 20, num_base_learners = 5, replace = TRUE
  )

options(mc.cores = n.threads)

predictions <- predict(surv.model.fit, COHORT.demissed[test.set[1:1000],])

COHORT.predictions <- COHORT.demissed[test.set[1:1000], c('surv_time_round', 'surv_event')]
COHORT.predictions$risk <- 1- predictions(predictions)[,50]

survConcordance(
  as.formula(paste0('Surv(surv_time_round, surv_event) ~ risk')),
  COHORT.predictions
)
