#' # Replicating Rapsomaniki _et al._ 2014 without imputation
#' 
#' Rapsomaniki and co-workers' paper creates a Cox hazard model to predict time
#' to death for patients using electronic health record data.
#' 
#' Because such data
#' are gathered as part of routine clinical care, some values may be missing.
#' For example, covariates include blood tests such as creatine level, and use
#' the most recent value no older than six months. A patient may simply not have
#' had this test in that time period.
#' 
#' The approach adopted in this situation, in common with many other similar
#' studies, is to impute missing values. This essentially involves finding
#' patients whose other values are similar to make a best guess at the likely
#' value of a missing datapoint. However, in the case of health records data,
#' the missingness of a value may be informative: the fact that a certain test
#' wasn't ordered could result from a doctor not thinking that test necessary,
#' meaning that its absence carries information.
#' 
#' Whilst there are many approaches which could be employed here, this program
#' represents arguably the simplest: we keep the Cox model as close to the
#' published version as possible but, instead of imputing the missing data, we
#' include it explicitly.
#' 
#' ## User variables
#' 
#' First, define some variables...

#+ define_vars

data.filename <- '../../data/cohort-sanitised.csv'
n.data <- 10000 # This is of full dataset...further rows may be excluded in prep

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

#' ## Setup

#+ setup, message=FALSE

source('../random-forest/shared.R')

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

#' ## Transform variables
#' 
#' The model uses variables which have been standardised in various ways, so
#' let's go through and transform our input variables in the same way...

#' ### Age
#' 
#' Risk of death is a nonlinear function of age at the start of the study, so
#' it is modelled with a cubic spline.

ageSpline <- function(x) {
  max((x-51)/10.289,0)^3 + 
    (69-51) * (max((x-84)/10.289,0)^3) -
    ((84-51) * (max(((x-69))/10.289,0))^3)/(84-69)
}

#' Let's have a look at what that function looks like...

df <- data.frame(
  x = 50:100,
  y = sapply(50:100, ageSpline)
)

ggplot(df, aes(x,y)) +
  geom_line()

#' The spline function looks approximately exponential, which makes sense.
#' 
#' ### Other variables
#' 
#' * Most other variables in the model are simply normalised by a factor with an
#'   offset, so build a data frame of these transformed variables.
#' * Variables which are factors (such as diagnosis) need to make sure that
#'   their first level is the one which was used as the baseline in the paper.
#' * The IMD (social deprivation) score is used by flagging those patients who
#'   are in the bottom quintile.

COHORT.scaled <-
  data.frame(
    ## Time to event
    time_death = COHORT.use$time_death,
    ## Death/censorship
    surv_event = COHORT.use$endpoint_death == 'Death',
    ## Rescaled age
    age_rescale = sapply(COHORT.use$age, ageSpline),
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
    total_chol_6mo_rescale = (COHORT.use$total_chol_6mo - 5),
    ## HDL, per 0.5 mmol/L increase
    hdl_6mo_rescale = (COHORT.use$hdl_6mo - 1.5) / 0.5,
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
    pulse_6mo_rescale = (COHORT.use$pulse_6mo - 70) / 10,
    ## Creatinine, per 30 μmol/L increase
    crea_6mo_rescale = (COHORT.use$crea_6mo - 60) / 30,
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo_rescale = (COHORT.use$total_wbc_6mo - 7.5) / 1.5,
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo_rescale = (COHORT.use$haemoglobin_6mo - 13.5) / 1.5
  )

#' ## Missing values
#' 
#' To incorporate missing values, we first make columns of logicals to
#' indicate whether the value was missing or not.

# Deal with missing values...go through all the columns
missing.suffix <- '_missing'
missing.cols <-
  c(
    'diagnosis',
    "most_deprived", "smokstatus", "total_chol_6mo_rescale", "hdl_6mo_rescale",
    "pulse_6mo_rescale", "crea_6mo_rescale", "total_wbc_6mo_rescale",
    "haemoglobin_6mo_rescale"
  )
for(surv.col in names(COHORT.scaled)) {
  # ...and, if it contains any missing values
  if(sum(is.na(COHORT.scaled[, surv.col])) > 0) {
    print(surv.col)
    # Create a new column which designates the missing ones
    COHORT.scaled[, paste0(surv.col, missing.suffix)] <-
      is.na(COHORT.scaled[, surv.col])
    # Then, deal with the actual values to remove their effect, depending on
    # variable type
    if(is.factor(COHORT.scaled[, surv.col])) {
      # If it's a factor, NAs can be their own level
      COHORT.scaled[, surv.col] <-
        factorNAfix(COHORT.scaled[, surv.col], NAval = 'missing')
    } else if(is.logical(COHORT.scaled[, surv.col])) {
      # Set the NA values to the baseline so they don't contribute to the model
      COHORT.scaled[is.na(COHORT.scaled[, surv.col]), surv.col] <- FALSE
    } else {
      # Set the NA values to the baseline so they don't contribute to the model
      COHORT.scaled[is.na(COHORT.scaled[, surv.col]), surv.col] <- 0
    }
  }
}

# make a survival object
COHORT.surv.train <- Surv(
  time  = COHORT.scaled[-test.set, 'time_death'],
  event = COHORT.scaled[-test.set, 'surv_event']
)

fit.exp <- survreg(
  COHORT.surv.train ~
    ### Sociodemographic characteristics #######################################
    ## Age in men, per year
    ## Age in women, per year
    ## Women vs. men
    # ie include interaction between age and gender!
    age_rescale*gender +
    ## Most deprived quintile, yes vs. no
    most_deprived +
    most_deprived_missing +
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
    ## Ex-smoker / current smoker vs. never
    smokstatus +
    smokstatus_missing +
    ## Hypertension, present vs. absent
    hypertension +
    ## Diabetes mellitus, present vs. absent
    diabetes_logical +
    ## Total cholesterol, per 1 mmol/L increase
    total_chol_6mo_rescale +
    total_chol_6mo_rescale_missing +
    ## HDL, per 0.5 mmol/L increase
    hdl_6mo_rescale +
    hdl_6mo_rescale_missing +
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
    pulse_6mo_rescale +
    pulse_6mo_rescale_missing +
    ## Creatinine, per 30 μmol/L increase
    crea_6mo_rescale +
    crea_6mo_rescale_missing +
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo_rescale +
    total_wbc_6mo_rescale_missing +
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo_rescale +
    haemoglobin_6mo_rescale_missing,
  data = COHORT.scaled[-test.set, ],
  dist = "exponential"
)

coefficients <- -fit.exp$coeff

cat(
  'Covariate                theirs         ours
  Age if man, years	        0.063230409    ', coefficients['age_rescale'],'
  Age if woman, years	      0.078075876    ', coefficients['age_rescale:genderWomen'],'
  Being a woman             -0.54888796    ', coefficients['genderWomen'],'
  DEPRIVATION: bottom 5th   0.140887623    ', coefficients['most_deprivedTRUE'],'
  DEPRIVATION: missing                     ', coefficients['most_deprived_missingTRUE'],'
  CAD DIAGNOSIS & SEVERITY
  Stable_angina (ref)	      0.03           ', 'hmm - how does ref have a beta?','
  Unstable angina	          0.023          ', coefficients['diagnosisUA'],'
  STEMI	                    0.079863191    ', coefficients['diagnosisSTEMI'],'
  NSTEMI                    0.260971786    ', coefficients['diagnosisNSTEMI'],'
  PCI last 6 months	        -0.429038846   ', coefficients['pci_6moTRUE'],'
  CABG last 6 months	      -0.661443369   ', coefficients['cabg_6moTRUE'],'
  Previous/recurrent MI	    0.128214367    ', coefficients['hx_miTRUE'],'
  Nitrates                  0.142262066    ', coefficients['long_nitrateTRUE'],'
  CVD RISK FACTORS
  SMOKING: Ex               0.104560407    ', coefficients['smokstatusEx'],'
  SMOKING: Current          0.28           ', coefficients['smokstatusCurrent'],'
  SMOKING: missing                         ', coefficients['most_deprived_missingTRUE'],'
  Hypertension              -0.035521708   ', coefficients['hypertensionTRUE'],'
  Diabetes                  0.185590982    ', coefficients['diabetes_logicalTRUE'],'
  Total cholesterol         0.012691881    ', coefficients['total_chol_6mo_rescale'],'
  Total cholesterol missing                ', coefficients['total_chol_6mo_rescale_missingTRUE'],'
  HDL, mmol/L               0.006510087    ', coefficients['hdl_6mo_rescale'],'
  HDL, missing                             ', coefficients['hdl_6mo_rescale_missingTRUE'],'
  CVD COMORBIDITIES
  Heart failure             0.43416416     ', coefficients['heart_failureTRUE'],'
  PAD                       0.251746864    ', coefficients['padTRUE'],'
  Atrial fibrillation       0.247297351    ', coefficients['hx_afTRUE'],'
  Prior stroke              0.284611916    ', coefficients['hx_strokeTRUE'],'
  CO-EXISTING MEDICAL CONDITIONS
  Chronic renal disease     0.110543627    ', coefficients['hx_renalTRUE'],'
  COPD                      0.140243473    ', coefficients['hx_copdTRUE'],'
  Cancer                    0.32014784     ', coefficients['hx_cancerTRUE'],'
  Chronic liver disease     0.48920724     ', coefficients['hx_liverTRUE'],'
  PSYCHOSOCIAL CHARACTERISTICS
  Depression                0.16530389     ', coefficients['hx_depressionTRUE'],'
  Anxiety                   0.159257697    ', coefficients['hx_anxietyTRUE'],'
  BIOMARKERS	
  Heart rate, beats/min     0.093975676    ', coefficients['pulse_6mo_rescale'],'
  Heart rate, missing                      ', coefficients['pulse_6mo_rescale_missingTRUE'],'
  Creatinine, mmol/L        0.063882169    ', coefficients['crea_6mo_rescale'],'
  Creatinine, missing                      ', coefficients['crea_6mo_rescale_missingTRUE'],'
  White cell count, 109/L   0.113975087    ', coefficients['total_wbc_6mo_rescale'],'
  White cell count, missing                ', coefficients['total_wbc_6mo_rescale_missingTRUE'],'
  Haemoglobin, g/dL         -0.26734457    ', coefficients['haemoglobin_6mo_rescale'],'
  Haemoglobin, missing                     ', coefficients['haemoglobin_6mo_rescale_missingTRUE'],
  '\n\n'
)

# Calculate C-indices on training and test sets
# Get C-indices for training and test sets
c.index.train <-
  cIndex(fit.exp, COHORT.scaled[-test.set, ], model.type = 'cph')
c.index.test <- 
  cIndex(fit.exp, COHORT.scaled[test.set, ], model.type = 'cph')

cat('C-index on training set:', c.index.train, '\n')
cat('C-index on test set:', c.index.test, '\n')
