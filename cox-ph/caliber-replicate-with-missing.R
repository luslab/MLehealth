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

n.data <- NA # This is of full dataset...further rows may be excluded in prep
endpoint <- 'death'

old.coefficients.filename <- 'rapsomaniki-cox-values-from-paper.csv'

output.filename.base <- '../../output/caliber-replicate-with-missing-survreg-2'

bootstraps <- 100
n.threads <- 8

# Create a few filenames from the base
new.coefficients.filename <- paste0(output.filename.base, '-coeffs-2.csv')
compare.coefficients.filename <-
  paste0(output.filename.base, '-bootstrap-2.csv')
cox.var.imp.perm.filename <- paste0(output.filename.base, '-var-imp-perm-2.csv')
model.filename <- paste0(output.filename.base, '-model-bootstrap.rds')

#' ## Setup

#+ setup, message=FALSE

source('../lib/shared.R')
require(xtable)
require(ggrepel)
source('caliber-scale.R')

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
#' it is modelled with a cubic spline. Let's have a look at what that function
#' looks like...

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
  caliberScale(COHORT.use, surv.time, surv.event)

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
COHORT.scaled.demissed <- prepCoxMissing(COHORT.scaled, missing.cols)

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
    age + gender +
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
    #hx_mi +
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
    haemoglobin_6mo_missing

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
    c.index = fit.exp.boot.ests['c.test', 'val'],
    c.index.lower = fit.exp.boot.ests['c.test', 'lower'],
    c.index.upper = fit.exp.boot.ests['c.test', 'upper'],
    calibration.score = fit.exp.boot.ests['calibration.score', 'val'],
    calibration.score.lower = fit.exp.boot.ests['calibration.score', 'lower'],
    calibration.score.upper = fit.exp.boot.ests['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)

#' ## Performance
#' 
#' ### C-index
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
#' ### Calibration
#' 
#' The bootstrapped calibration score is
#' **`r round(fit.exp.boot.ests['calibration.score', 'val'], 3)`
#' (`r round(fit.exp.boot.ests['calibration.score', 'lower'], 3)` - 
#' `r round(fit.exp.boot.ests['calibration.score', 'upper'], 3)`)**.
#' 
#' Let's draw a representative curve from the unbootstrapped fit... (It would be
#' better to draw all the curves from the bootstrap fit to get an idea of
#' variability, but I've not implemented this yet.)
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(fit.exp, COHORT.scaled.demissed[test.set, ])

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score[['area']], 3)`** +/-
#' **`r round(calibration.score[['se']], 3)`**.
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
write.csv(compare.coefficients, new.coefficients.filename)

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
#' To establish how important a given variable is in determining outcome (and to
#' compare with other measures such as variable importance with random forests),
#' it would be worthwhile to calculate an equivalent measure for a Cox model.
#' The easiest way to do this across techniques is to look at it by permutation:
#' randomly permute values in each variable, and see how much worse it makes the
#' prediction.
#' 
#+ cox_variable_importance

cox.var.imp.perm <- 
  generalVarImp(
    fit.exp, COHORT.scaled.demissed[test.set, ], model.type = 'survreg'
  )

write.csv(cox.var.imp.perm, cox.var.imp.perm.filename, row.names = FALSE)

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
      Surv(surv_time, surv_event) ~ gender,
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
      Surv(surv_time, surv_event) ~ age + gender,
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
      Surv(surv_time, surv_event) ~ age * gender,
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
#' 
#' ## Systematic analysis of missing values
#' 
#' Clearly the danger (or otherwise) of having missing values differs by
#' variable, according to this model. Let's look at all of the variables, to see
#' how all missing values compare, and whether they make sense.
#'
#+ missing_values_vs_ranges

# For continuous variables, get risk values for all patients for violin plots
risk.dist.by.var <- data.frame()
# For categorical variables, let's get all levels of the factor
risk.cats <- data.frame()
# Quantities not to plot in this graph
exclude.quantities <-
  c(
    'age', # too large a risk range, and no missing values anyway
    'gender', 'diagnosis', # no missing values
    'diabetes' # is converted to diabetes_logical so causes an error if included
    )

for(quantity in surv.predict) {
  # If we're not excluding..
  if(!(quantity %in% exclude.quantities)) {
    # If it's numeric
    if(is.numeric(COHORT.scaled[, quantity])) {
      # Get risks by taking the scaled values (NOT processed for missing ones, or
      # there will be a lot of 0s in there) and taking quantity risks to that power
      risk <-
        new.coefficients$our_value[new.coefficients$quantity.level == quantity] ^
        COHORT.scaled[, quantity]
      # Discard outliers with absurd values
      inliers <- inRange(risk, quantile(risk, c(0.01, 0.99), na.rm = TRUE))
      risk <- risk[inliers]
      risk.dist.by.var <-
        rbind(
          risk.dist.by.var,
          data.frame(quantity, risk)
        )
      risk.cats <-
        rbind(
          risk.cats,
          data.frame(
            quantity,
            new.coefficients[
              new.coefficients$quantity.level ==
                paste0(quantity, '_missingTRUE'),
              ]
          )
        )
    } else if(is.factor(COHORT.scaled[, quantity])) {
      risk.cats <-
        rbind(
          risk.cats,
          data.frame(
            quantity,
            new.coefficients[
              startsWith(as.character(new.coefficients$quantity), quantity),
            ]
          )
        )
    } else {
      # We're not going to include logicals in this plot, because they cannot
      # be missing, by definition, in the logicals used here.
      # There shouldn't be any other kinds of variables.
      # If your code requires them, put something here.
    }
  }
}

# Save the results
write.csv(risk.dist.by.var, paste0(output.filename.base, '-risk-violins.csv'))
write.csv(risk.cats, paste0(output.filename.base, '-risk-cats.csv'))

# Plot the results
ggplot() +
  # First, and therefore at the bottom, draw the reference line at risk = 1
  geom_hline(yintercept = 1) +
  # Then, on top of that, draw the violin plot of the risk from the data
  geom_violin(data = risk.dist.by.var, aes(x = quantity, y = risk)) +
  geom_pointrange(
    data = risk.cats,
    aes(x = quantity, y = our_value, ymin = our_lower,
        ymax = our_upper),
    
    position = position_jitter(width = 0.1)
  ) +
  geom_text(
    data = risk.cats,
    aes(
      x = quantity,
      y = our_value,
      label = quantity.level
    )
  )

#' ## Conclusion
#' 
#' That's all folks! Just to confirm that this is reproducible, let's check that
#' our random number generator is lined up as expected after all those
#' calculations. If your seed is the same, you should expect the two random
#' strings to be the same too.
#' 
#' Original random seed: 35498 (found in ``../random-forests/shared.R``);
#' current random seed: `r random.seed`
#' 
#' 
#' Original random string: ``HPMFNEVRGMWPBYYLNUBB``
#' 
#' Current random string: ```r randomString(20, LETTERS)```