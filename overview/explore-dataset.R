#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Simple data investigations
#' 
#+ setup, message=FALSE

data.filename <- '../../data/cohort-sanitised.csv'

source('../lib/shared.R')
requirePlus('rms')

COHORT.full <- data.frame(fread(data.filename))

COHORT.use <- subset(COHORT.full, !exclude)

#' ## Missing data
#' 
#' How much is there, and where is it concentrated?
#' 
#+ missing_data_plot

interesting.vars <- 
  c(
    'age', 'gender', 'diagnosis', 'pci_6mo', 'cabg_6mo',
    'hx_mi', 'long_nitrate', 'smokstatus', 'hypertension', 'diabetes',
    'total_chol_6mo', 'hdl_6mo', 'heart_failure', 'pad', 'hx_af', 'hx_stroke',
    'hx_renal', 'hx_copd', 'hx_cancer', 'hx_liver', 'hx_depression',
    'hx_anxiety', 'pulse_6mo', 'crea_6mo', 'total_wbc_6mo','haemoglobin_6mo',
    'imd_score'
  )

missingness <- 
  unlist(lapply(COHORT.use[, interesting.vars], function(x){sum(is.na(x))}))

missingness <- data.frame(var = names(missingness), n.missing = missingness)

missingness$pc.missing <- missingness$n.missing / nrow(COHORT.use)

ggplot(subset(missingness, n.missing > 0), aes(x = var, y = pc.missing)) +
  geom_bar(stat = 'identity') +
  ggtitle('% missingness by variable')

#' Are any variables commonly found to be jointly missing?
#' 
#+ missing_jointly_plot

COHORT.missing <-
  data.frame(
    lapply(COHORT.use[, interesting.vars], function(x){is.na(x)})
  )

COHORT.missing.cor <- data.frame(
  t(combn(1:ncol(COHORT.missing), 2)),
  var1 = NA, var2 = NA, joint.missing = NA
)

for(i in 1:nrow(COHORT.missing.cor)) {
  var1 <- sort(names(COHORT.missing))[COHORT.missing.cor[i, 'X1']]
  var2 <- sort(names(COHORT.missing))[COHORT.missing.cor[i, 'X2']]
  COHORT.missing.cor[i, c('var1', 'var2')] <- c(var1, var2)
  if(any(COHORT.missing[, var1]) & any(COHORT.missing[, var2])) {
    COHORT.missing.cor[i, 'joint.missing'] <-
      sum(!(COHORT.missing[, var1]) & !(COHORT.missing[, var2])) / 
      sum(!(COHORT.missing[, var1]) | !(COHORT.missing[, var2]))
  }
}

ggplot(subset(COHORT.missing.cor, !is.na(joint.missing)), aes(x = var1, y = var2, fill = joint.missing)) +
  geom_tile()

#' Some tests are usually ordered together. Are they missing together?
#' 
#+ missing_venn

ggplot(
  data.frame(
    x = 1,
    category = factor(c('HDL', 'both', 'total cholesterol', 'neither'),
                      levels = c('HDL', 'both', 'total cholesterol', 'neither')),
    val = c(
      sum(!COHORT.missing$hdl_6mo) - sum(!COHORT.missing$hdl_6mo & !COHORT.missing$total_chol_6mo),
      sum(!COHORT.missing$hdl_6mo & !COHORT.missing$total_chol_6mo),
      sum(!COHORT.missing$total_chol_6mo) - sum(!COHORT.missing$hdl_6mo & !COHORT.missing$total_chol_6mo),
      sum(COHORT.missing$hdl_6mo & COHORT.missing$total_chol_6mo)
    )
  ),
  aes(x = x, y = val, fill = category)
) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c("#cc0000", "#990099", "#0000cc", '#cccccc'))

ggplot(
  data.frame(
    x = 1,
    category = factor(c('haemoglobin', 'both', 'WBC', 'neither'),
                      levels = c('haemoglobin', 'both', 'WBC', 'neither')),
    val = c(
      sum(!COHORT.missing$haemoglobin_6mo) - sum(!COHORT.missing$haemoglobin_6mo & !COHORT.missing$total_wbc_6mo),
      sum(!COHORT.missing$haemoglobin_6mo & !COHORT.missing$total_wbc_6mo),
      sum(!COHORT.missing$total_wbc_6mo) - sum(!COHORT.missing$haemoglobin_6mo & !COHORT.missing$total_wbc_6mo),
      sum(COHORT.missing$hdl_6mo & COHORT.missing$total_wbc_6mo)
    )
  ),
  aes(x = x, y = val, fill = category)
) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c("#cc0000", "#990099", "#0000cc", '#cccccc'))
  

#' ### Survival and missingness
#' 
#' Let's plot survival curves for a few types of data by missingness...
#' 
#+ missingness_survival

COHORT.use$surv <- with(COHORT.use, Surv(time_death, endpoint_death == 'Death'))

plotSurvSummary <- function(df, var) {
  df$var_summary <- factorNAfix(binByQuantile(df[, var], c(0,0.1,0.9,1)))
  surv.curve <- npsurv(surv ~ var_summary, data = df)
  print(survplot(surv.curve, ylab = var))
}

plotSurvSummary(COHORT.use, 'total_wbc_6mo')

surv.curve <- npsurv(surv ~ is.na(crea_6mo), data = COHORT.use)
survplot(surv.curve)

surv.curve <- npsurv(surv ~ is.na(haemoglobin_6mo), data = COHORT.use)
survplot(surv.curve)

surv.curve <- npsurv(surv ~ is.na(hdl_6mo), data = COHORT.use)
survplot(surv.curve)

surv.curve <- npsurv(surv ~ is.na(pulse_6mo), data = COHORT.use)
survplot(surv.curve)

surv.curve <- npsurv(surv ~ is.na(smokstatus), data = COHORT.use)
survplot(surv.curve)

surv.curve <- npsurv(surv ~ is.na(total_chol_6mo), data = COHORT.use)
survplot(surv.curve)

surv.curve <- npsurv(surv ~ is.na(total_wbc_6mo), data = COHORT.use)
survplot(surv.curve)

#' Variables where it's safer to be missing data:
#'  * Creatinine
#'  * 
#' 
#' Variables where it's safer to have data:
#'  * 

#' ## Data distributions
#' 
#' How are the data distributed?
#' 
#+ data_distributions

ggplot(COHORT.full) +
  geom_histogram(aes(x = crea_6mo, fill = crea_6mo %% 10 != 0), binwidth = 1) +
  xlim(0, 300)

#' Creatinine levels are fairly smoothly distributed. The highlighted bins
#' indicate numerical values divisible by 10, and there seems to be no
#' particular bias. The small cluster of extremely low values could be
#' misrecorded somehow.

ggplot(COHORT.full) +
  geom_histogram(aes(x = pulse_6mo, fill = pulse_6mo %% 4 != 0), binwidth = 1) +
  xlim(0, 150)

#' Heart rate data have high missingness, and those few values we do have are
#' very heavily biased towards multiples of 4. This is likely because heart rate
#' is commonly measured for 15 seconds and then multiplied up to give a result
#' in beats per minute! There is also a bias towards round numbers, with large
#' peaks at 60, 80, 100 and 120...