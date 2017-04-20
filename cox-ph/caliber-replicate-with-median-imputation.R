#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

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
n.data <- NA # This is of full dataset...further rows may be excluded in prep
endpoint <- 'death'

old.coefficients.filename <- 'rapsomaniki-cox-values-from-paper.csv'
compare.coefficients.filename <-
  '../../output/caliber-replicate-median-imp-survreg-bootstrap-1.csv'
cox.var.imp.perm.filename <-
  '../../output/caliber-replicate-median-imp-survreg-bootstrap-var-imp-perm-1.csv'
cox.var.imp.perm.missing.filename <-
  '../../output/caliber-replicate-with-missing-survreg-bootstrap-var-imp-perm-1.csv'
model.filename <-
  '../../output/caliber-replicate-median-imp-model-survreg-bootstrap-1.rds'

bootstraps <- 100
n.threads <- 8

#' ## Setup

#+ setup, message=FALSE

source('../random-forest/shared.R')
require(xtable)
require(ggrepel)

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

#' OK, we've now got **`r n.data`** patients, split into a training set of
#' `r n.data - length(test.set)` and a test set of `r length(test.set)`.
#'
#'
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
#'   offset.
#' * Variables which are factors (such as diagnosis) need to make sure that
#'   their first level is the one which was used as the baseline in the paper.
#' * The IMD (social deprivation) score is used by flagging those patients who
#'   are in the bottom quintile.

COHORT.scaled <-
  data.frame(
    ## Time to event
    surv_time = COHORT.use[, surv.time],
    ## Death/censorship
    surv_event = COHORT.use[, surv.event] %in% surv.event.yes,
    ## Rescaled age
    age = sapply(COHORT.use$age, ageSpline),
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
    ## Creatinine, per 30 μmol/L increase
    crea_6mo = (COHORT.use$crea_6mo - 60) / 30,
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo = (COHORT.use$total_wbc_6mo - 7.5) / 1.5,
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo = (COHORT.use$haemoglobin_6mo - 13.5) / 1.5
  )

#' ## Missing values
#' 
#' To deal with missing values, we will do a very simple kind of imputation,
#' replacing them with the median of the nonmissing values.

# Specify missing columns - diagnosis only has a handful of missing values so
# sometimes doesn't have missing ones in the sampled training set, meaning
# prepCoxMissing wouldn't fix it.
missing.cols <-
  c(
    "diagnosis", "most_deprived", "smokstatus", "total_chol_6mo", "hdl_6mo",
    "pulse_6mo", "crea_6mo", "total_wbc_6mo", "haemoglobin_6mo"
  )
COHORT.scaled.demissed <- medianImpute(COHORT.scaled, missing.cols)

# make a survival object
COHORT.surv.train <- Surv(
  time  = COHORT.scaled.demissed[-test.set, 'surv_time'],
  event = COHORT.scaled.demissed[-test.set, 'surv_event']
)

#' ## Survival fitting
#' 
#' Fit a Cox model to the preprocessed data. The paper uses a Cox model with an
#' exponential baseline hazard, as here. The standard errors were calculated
#' with 200 bootstrap samples, which we're also doing here.

#+ fit_cox_model, cache=cacheoption

surv.formula <-
  Surv(surv_time, surv_event) ~
    ### Sociodemographic characteristics #######################################
    ## Age in men, per year
    ## Age in women, per year
    ## Women vs. men
    # ie include interaction between age and gender!
    age*gender +
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
    ## Creatinine, per 30 μmol/L increase
    crea_6mo +
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo +
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo

# Fit the main model with all data
fit.exp <- survreg(
  formula = surv.formula,
  data = COHORT.scaled.demissed[-test.set, ],
  dist = "exponential"
)

# Run a bootstrap on the model to find uncertainties
fit.exp.boot <- 
  boot(
    formula = surv.formula,
    data = COHORT.scaled.demissed[-test.set, ],
    statistic = bootstrapFitSurvreg,
    R = bootstraps,
    parallel = 'multicore',
    ncpus = n.threads,
    test.data = COHORT.scaled.demissed[test.set, ]
  )

# Save the fit, because it might've taken a while!
saveRDS(fit.exp.boot, model.filename)

# Unpackage the uncertainties from the bootstrapped data
fit.exp.boot.ests <-  bootStats(fit.exp.boot, uncertainty = '95ci')

# Save bootstrapped performance values
varsToTable(
  data.frame(
    model = 'cox',
    imputation = FALSE,
    discretised = FALSE,
    c.index = fit.exp.boot.ests['c.train', 'val'],
    c.index.lower = fit.exp.boot.ests['c.train', 'lower'],
    c.index.upper = fit.exp.boot.ests['c.train', 'upper']
  ),
  performance.file
)

#' ## Performance
#' 
#' Having fitted the Cox model, how did we do? The c-indices were calculated as
#' part of the bootstrapping, so we just need to take a look at those...
#' 
#' C-indices are **`r round(fit.exp.boot.ests['c.train', 'val'], 3)`
#' (`r round(fit.exp.boot.ests['c.train', 'lower'], 3)` - 
#' `r round(fit.exp.boot.ests['c.train', 'upper'], 3)`)** on the training set and
#' **`r round(fit.exp.boot.ests['c.test', 'val'], 3)`
#' (`r round(fit.exp.boot.ests['c.test', 'lower'], 3)` - 
#' `r round(fit.exp.boot.ests['c.test', 'upper'], 3)`)** on the test set.
#' Not too bad!
#' 
#' 
#' ## Coefficients
#' 
#' As well as getting comparable C-indices, it's also worth checking to see how
#' the risk coefficients calculated compare to those found in the original
#' paper. Let's compare...

# Load CSV of values from paper
old.coefficients <- read.csv(old.coefficients.filename)

# Get coefficients from this fit
new.coefficients <-
  bootStats(fit.exp.boot, uncertainty = '95ci', transform = negExp)
names(new.coefficients) <- c('our_value', 'our_lower', 'our_upper')
new.coefficients$quantity.level <- rownames(new.coefficients)

# Create a data frame comparing them
compare.coefficients <- merge(old.coefficients, new.coefficients)

# Kludge because age:genderWomen is the pure interaction term, not the risk for
# a woman per unit of advancing spline-transformed age
compare.coefficients[
  compare.coefficients$quantity.level == 'age:genderWomen', 'our_value'
  ] <-
  compare.coefficients[
    compare.coefficients$quantity.level == 'age:genderWomen', 'our_value'
    ] *
  compare.coefficients[
    compare.coefficients$quantity.level == 'age', 'our_value'
    ]

# Save CSV of results
write.csv(compare.coefficients, output.filename)

# Plot a graph by which to judge success
ggplot(compare.coefficients, aes(x = their_value, y = our_value)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 1, colour = 'grey') +
  geom_vline(xintercept = 1, colour = 'grey') +
  geom_point() +
  geom_errorbar(aes(ymin = our_lower, ymax = our_upper)) +
  geom_errorbarh(aes(xmin = their_lower, xmax = their_upper)) +
  geom_text_repel(aes(label = long_name)) +
  theme_classic(base_size = 8)

#+ coefficients_table, results='asis'

print(
  xtable(
    data.frame(
      variable =
        paste(
          compare.coefficients$long_name, compare.coefficients$unit, sep=', '
        ),
      compare.coefficients[c('our_value', 'their_value')]
    ),
    digits = c(0,0,3,3)
  ),
  type = 'html',
  include.rownames = FALSE
)

#' ### Variable importance
#' 
#' Let's compare the variable importance from this method with accounting for
#' missing values explicitly...
#' 
#+ cox_variable_importance

cox.var.imp.perm <- 
  generalVarImp(
    fit.exp, COHORT.scaled.demissed[test.set, ], model.type = 'survreg'
  )

write.csv(cox.var.imp.perm, cox.var.imp.perm.filename, row.names = FALSE)

cox.var.imp.perm.missing <- read.csv(cox.var.imp.perm.missing.filename)

cox.var.imp.comparison <-
  merge(
    cox.var.imp.perm,
    cox.var.imp.perm.missing,
    by = 'var',
    suffixes = c('', '.missing')
  )

ggplot(cox.var.imp.comparison, aes(x = var.imp.missing, y = var.imp)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

#' There's a good correlation!
