#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- TRUE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Modelling with age
#' 
#' It seems that, no matter what I do, the C-index of a model, random forest or
#' otherwise, is about 0.78. I decided to try to some simpler models, intially
#' based purely on age which is clearly the biggest factor in this dataset.
#' 
#' The sanity check works: giving the model more data does indeed result in a
#' better fit. However, on top of that, I was surprised by just how good the
#' performance can be when age alone is considered!

#+ user_variables, message=FALSE

data.filename <- '../../data/cohort-sanitised.csv'

n.trees <- 500

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )
untransformed.vars <- c('anonpatid', 'time_death', 'imd_score', 'exclude')

source('shared.R')
require(ggrepel)

#' ## Load and prepare data
#+ load_and_prepare_data

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# Define process settings; nothing for those to not transform, and missingToBig
# for the continuous ones...
process.settings <-
  list(
    var        = c(untransformed.vars, continuous.vars),
    method     =
      c(
        rep(NA, length(untransformed.vars)),
        rep('missingToBig', length(continuous.vars))
      ),
    settings   = rep(NA, length(untransformed.vars) + length(continuous.vars))
  )

COHORT.prep <-
  prepData(
    # Data for cross-validation excludes test set
    COHORT.full,
    cols.keep,
    process.settings,
    surv.time, surv.event,
    surv.event.yes,
    extra.fun = caliberExtraPrep
  )
n.data <- nrow(COHORT.prep)

# Define indices of test set
test.set <- sample(1:n.data, (1/3)*n.data)

#' ## Models
#' 
#' ### The normal model
#' 
#' All the variables, as in the vector `surv.predict`.
#+ normal_model, cache=cacheoption

# Fit random forest
surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.prep[-test.set,],
    model.type = 'ranger',
    n.trees = n.trees,
    split.rule = 'logrank',
    n.threads = 8,
    respect.unordered.factors = 'partition'
  )

print(surv.model.fit)

# Get C-index
c.index.test <- 
  cIndex(surv.model.fit, COHORT.prep[test.set, ], model.type = 'ranger')

#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.

 
#' ### Just age
#' 
#' What if all we had to go on was age?
#+ just_age_model, cache=cacheoption

# Fit random forest
surv.model.fit <-
  survivalFit(
    c('age'),
    COHORT.prep[-test.set,],
    model.type = 'ranger',
    n.trees = n.trees,
    split.rule = 'logrank',
    n.threads = 8,
    respect.unordered.factors = 'partition'
  )

print(surv.model.fit)

# Get C-index
c.index.test <- 
  cIndex(surv.model.fit, COHORT.prep[test.set, ], model.type = 'ranger')

#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.

#' ### No model, literally just age
#' 
#' What if we constructed the C-index based purely on patients' ages?
#+ just_age_cindex

c.index.age <-
  as.numeric(
    survConcordance(
      Surv(time_death, surv_event) ~ age,
      COHORT.prep
    )$concordance
  )

#' The C-index on the whole dataset based purely on age is
#' **`r round(c.index.test, 3)`**. That's most of our predictive accuracy right
#' there! Reassuringly, it's also equal to the value predicted by the random
#' forest model based purely on age...

#' ### Age and gender
#' 
#' OK, age and gender.
#+ age_gender_model, cache=cacheoption

# Fit random forest
surv.model.fit <-
  survivalFit(
    c('age', 'gender'),
    COHORT.prep[-test.set,],
    model.type = 'ranger',
    n.trees = n.trees,
    split.rule = 'logrank',
    n.threads = 8,
    respect.unordered.factors = 'partition'
  )

print(surv.model.fit)

# Get C-index
c.index.test <- 
  cIndex(surv.model.fit, COHORT.prep[test.set, ], model.type = 'ranger')

#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.

#' ### Just gender
#' 
#' Just gender, as a sanity check.
#+ just_gender_model, cache=cacheoption

# Fit random forest
surv.model.fit <-
  survivalFit(
    c('gender'),
    COHORT.prep[-test.set,],
    model.type = 'ranger',
    n.trees = n.trees,
    split.rule = 'logrank',
    n.threads = 8,
    respect.unordered.factors = 'partition'
  )

print(surv.model.fit)

# Get C-index
c.index.test <- 
  cIndex(surv.model.fit, COHORT.prep[test.set, ], model.type = 'ranger')

#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.

#' ### Everything except age
#' 
#' How do we do if we use all the variables _except_ age?
#+ no_age_model, cache=cacheoption

# Fit random forest
surv.model.fit <-
  survivalFit(
    surv.predict[surv.predict != 'age'],
    COHORT.prep[-test.set,],
    model.type = 'ranger',
    n.trees = n.trees,
    split.rule = 'logrank',
    n.threads = 8,
    respect.unordered.factors = 'partition'
  )

print(surv.model.fit)

# Get C-index
c.index.test.all.not.age <- 
  cIndex(surv.model.fit, COHORT.prep[test.set, ], model.type = 'ranger')

#' The C-index on the held-out test set is
#' **`r round(c.index.test.all.not.age, 3)`**.


#' ### Predicting age
#' 
#' So why does the model which doesn't include age do so well? Clearly the other
#' variables allow you to predict age with reasonable accuracy... So let's try
#' just that as a final test.
#+ predict_age, cache=cacheoption

age.model <-
  ranger(
    formula(
      paste0(
        # Predicting just the age
        'age ~ ',
        # Predictor variables then make up the other side
        paste(surv.predict[surv.predict != 'age'], collapse = '+')
      )
    ),
    COHORT.prep[-test.set, ],
    num.trees = n.trees,
    splitrule = 'variance',
    num.threads = 8,
    respect.unordered.factors = 'partition'
  )

age.predictions <- predict(age.model, COHORT.prep[test.set, ])

age.cor <- cor(age.predictions$predictions, COHORT.prep[test.set, 'age'])

to.plot <-
  data.frame(
    age = COHORT.prep[test.set, 'age'],
    predicted = age.predictions$predictions
  )

ggplot(sample.df(to.plot, 10000), aes(x = age, y = predicted)) +
  geom_point(alpha = 0.2)

#' It doesn't look that great, but there is some correlation...
#' r^2 = `r age.cor^2` which is unremarkable, but OK. A more relevant measure
#' would be the pure-age C-index, ie if I gave you a pair of patients, how
#' often could you predict who was older?

c.index.age.on.age <-
  1 - as.numeric(
    survConcordance(
       age ~ predicted,
       to.plot
    )$concordance
  )

#' This comes out as **`r round(c.index.age.on.age, 3)`**, as compared to
#' **`r round(c.index.test.all.not.age, 3)`**, which was the C-index for the
#' survival model based on all other variables. This is slightly worse, but
#' that is what you'd expect: in some sense it's the amount of the C-index for
#' survival accounted for by the random survival forest inferring age from other
#' variables, ie, as expected, most of it.

#' ## Conclusion
#' 
#' So, in conclusion, the sanity check does appear to show sanity: giving the
#' random forest model more data to work with improves its performance. The
#' reasonable performance on the model of all factors except age is explained in
#' part by those factors' ability to act as a proxy for age (eg perhaps people
#' get more diseases as they get older, so if you have a, b _and_ c you're
#' likely to be older).