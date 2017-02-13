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

output.filename <- '../../output/caliber-replicate-with-missing-try1.csv'

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
    time_death = COHORT.use$time_death,
    ## Death/censorship
    surv_event = COHORT.use$endpoint_death == 'Death',
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
#' To incorporate missing values, we need to do different things depending on
#' the variable type.
#' 
#' * If the variable is a factor, we can simply add an extra factor level
#'   'missing' to account for missing values.
#' * For logical or numerical values, we first make a new column of logicals to
#'   indicate whether the value was missing or not. Those logicals will then
#'   themselves be variables with associated beta values, giving the risk
#'   associated with having a value missing. Then, set the a logical to FALSE or
#'   a numeric to 0, the baseline value, therefore not having any effect on the
#'   final survival curve.

# Specify missing columns - diagnosis only has a handful of missing values so
# sometimes doesn't have missing ones in the sampled training set, meaning
# prepCoxMissing wouldn't fix it.
missing.cols <-
  c(
    "diagnosis", "most_deprived", "smokstatus", "total_chol_6mo", "hdl_6mo",
    "pulse_6mo", "crea_6mo", "total_wbc_6mo", "haemoglobin_6mo"
  )
COHORT.scaled <- prepCoxMissing(COHORT.scaled, missing.cols)

# make a survival object
COHORT.surv.train <- Surv(
  time  = COHORT.scaled[-test.set, 'time_death'],
  event = COHORT.scaled[-test.set, 'surv_event']
)

#' ## Survival fitting
#' 
#' Fit a Cox model to the preprocessed data. The paper uses a Cox model with an
#' exponential baseline hazard, as here.

fit.exp <- survreg(
  COHORT.surv.train ~
    ### Sociodemographic characteristics #######################################
    ## Age in men, per year
    ## Age in women, per year
    ## Women vs. men
    # ie include interaction between age and gender!
    age*gender +
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
    ## Ex-smoker / current smoker / missing data vs. never
    smokstatus +
    ## Hypertension, present vs. absent
    hypertension +
    ## Diabetes mellitus, present vs. absent
    diabetes_logical +
    ## Total cholesterol, per 1 mmol/L increase
    total_chol_6mo +
    total_chol_6mo_missing +
    ## HDL, per 0.5 mmol/L increase
    hdl_6mo +
    hdl_6mo_missing +
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
    pulse_6mo_missing +
    ## Creatinine, per 30 μmol/L increase
    crea_6mo +
    crea_6mo_missing +
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo +
    total_wbc_6mo_missing +
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo +
    haemoglobin_6mo_missing,
  data = COHORT.scaled[-test.set, ],
  dist = "exponential"
)

#' ## Performance
#' 
#' Having fitted the Cox model, how did we do?
#' 

# Calculate C-indices on training and test sets
c.index.train <-
  cIndex(fit.exp, COHORT.scaled[-test.set, ], model.type = 'survreg')
c.index.test <- 
  cIndex(fit.exp, COHORT.scaled[test.set, ], model.type = 'survreg')

#' C-indices are **`r round(c.index.train, 3)`** on the training set and
#' **`r round(c.index.test, 3)`** on the test set. Not too bad!
#' 
#' 
#' ## Coefficients
#' 
#' As well as getting comparable C-indices, it's also worth checking to see how
#' the risk coefficients calculated compare to those found in the original
#' paper. Let's compare...

# Load CSV of values from paper
old.coefficients <- read.csv('rapsomaniki-cox-values-from-paper.csv')

# Get coefficients from this fit
new.coefficients <- as.data.frame(-fit.exp$coeff)
names(new.coefficients) <- 'our_value'
new.coefficients$our_value <- exp(new.coefficients$our_value)
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

#' ## Strange things about gender
#' 
#' The only huge outlier in this comparison is gender: according to the paper,
#' being a woman is significantly safer than being a man, whereas we
#' don't find a striking difference between genders.
#' 
#' Let's try some simple fits to see if this is explicable.
#' 
#' ### Fit based only on gender

print(
  exp(
    -survreg(
      COHORT.surv.train ~ gender,
      data = COHORT.scaled[-test.set, ],
      dist = "exponential"
    )$coeff
  )
)

#' According to a model based only on gender, women are at higher risk than men.
#' This suggests that there are indeed confounding factors in this dataset.
#' 
#' ### Fit based on age and gender

print(
  exp(
    -survreg(
      COHORT.surv.train ~ age + gender,
      data = COHORT.scaled[-test.set, ],
      dist = "exponential"
    )$coeff
  )
)

#' Once age is taken into account, it is (obviously) more dangerous to be older,
#' and it is once again safer to be female.

ggplot(COHORT.use, aes(x = age, group = gender, fill = gender)) +
  geom_histogram(alpha = 0.8, position = 'dodge', binwidth = 1)

#' And, as we can see from this histogram, this is explained by the fact that
#' the women in this dataset do indeed tend to be older than the men.

print(
  exp(
    -survreg(
      COHORT.surv.train ~ age * gender,
      data = COHORT.scaled[-test.set, ],
      dist = "exponential"
    )$coeff
  )
)

#' Finally, looking at a model where age and gender are allowed to interact, the
#' results are once again a little confusing: being older is worse, which makes
#' sense, but being a woman is again worse. This is then offset by the slightly
#' positive interaction between age and being a woman, meaning that the overall
#' effect of being a slightly older woman will be positive (especially because
#' the spline function gets much larger with advancing age).
#' 
#' It's important to remember that the output of any regression model with many
#' terms gives the strength of a given relationship having already controlled
#' for variation due to those other terms. In this case, even just using two
#' terms can give counter-intuitive results if you try to interpret coefficients
#' in isolation.
#' 
#' ## Coefficients for missing values
#' 
#' Let's see how coefficients for missing values compare to the range of risks
#' implied by the range of values for data...

# Value of creatinine which would result in a 25% increased risk of death
creatinine.25pc.risk <- 60 + 30 * log(1.25)/log(compare.coefficients[
  compare.coefficients$quantity.level == 'crea_6mo', 'our_value'
  ])

# Equivalent value of creatinine for patients with missing data
creatinine.missing <- 60 + 30 * 
  log(compare.coefficients[
    compare.coefficients$quantity.level == 'crea_6mo_missingTRUE', 'our_value'
    ]) /
  log(compare.coefficients[
    compare.coefficients$quantity.level == 'crea_6mo', 'our_value'
  ])

text.y.pos <- 3500

ggplot(
  # Only plot the lowest 95 percentiles of data due to outliers
  subset(COHORT.use, crea_6mo < quantile(crea_6mo, 0.95, na.rm = TRUE)),
  aes(x = crea_6mo)
) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 60) +
  annotate("text", x = 60, y = text.y.pos, angle = 270, hjust = 0, vjust = 1,
           label = "Baseline") +
  geom_vline(xintercept = creatinine.25pc.risk) +
  annotate("text", x = creatinine.25pc.risk, y = text.y.pos, angle = 270,
           hjust = 0, vjust = 1, label = "25% more risk") +
  geom_vline(xintercept = creatinine.missing) +
  annotate("text", x = creatinine.missing, y = text.y.pos, angle = 270,
           hjust = 0, vjust = 1, label = "missing data eqv")

#' So, having a missing creatinine value is slightly safer than having a
#' baseline reading, which is on the low end of observed values in this cohort.
#' Even an extremely high creatinine reading doesn't confer more than 25%
#' additional risk.

# Value which would result in a 25% increased risk of death
hdl.10pc.risk <- 1.5 + 0.5 * log(1.1)/log(compare.coefficients[
  compare.coefficients$quantity.level == 'hdl_6mo', 'our_value'
  ])

# Equivalent value for patients with missing data
hdl.missing <- 1.5 + 0.5 *
  log(compare.coefficients[
    compare.coefficients$quantity.level == 'hdl_6mo_missingTRUE', 'our_value'
    ]) /
  log(compare.coefficients[
    compare.coefficients$quantity.level == 'hdl_6mo', 'our_value'
    ])

text.y.pos <- 4000

ggplot(
  # Only plot the lowest 95 percentiles of data due to outliers
  subset(COHORT.use, hdl_6mo < quantile(hdl_6mo, 0.95, na.rm = TRUE)),
  aes(x = hdl_6mo)
) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 1.5) +
  annotate("text", x = 1.5, y = text.y.pos, angle = 270, hjust = 0, vjust = 1,
           label = "Baseline") +
  geom_vline(xintercept = hdl.10pc.risk) +
  annotate("text", x = hdl.10pc.risk, y = text.y.pos, angle = 270, hjust = 0,
           vjust = 1, label = "10% more risk") +
  geom_vline(xintercept = hdl.missing) +
  annotate("text", x = hdl.missing, y = text.y.pos, angle = 270, hjust = 0,
           vjust = 1, label = "missing data eqv")

#' Missing HDL values seem to have the opposite effect: it's substantially more
#' risky to have a blank value here. One possible explanation is that this is
#' such a common blood test that its absence indicates that the patient is not
#' seeking or receiving adequate medical care.