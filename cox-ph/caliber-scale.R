ageSpline <- function(x) { 
  max((x-51)/10.289,0)^3 +  
    (69-51) * (max((x-84)/10.289,0)^3) - 
    ((84-51) * (max(((x-69))/10.289,0))^3)/(84-69) 
}

caliberScaleUnits <- function(x, quantity) {
  if(quantity == 'age') {
    ## Spline function
    x <- sapply(x, ageSpline)
  } else if(quantity == 'total_chol_6mo') {
    ## Total cholesterol, per 1 mmol/L increase
    x <- x - 5
  } else if(quantity == 'hdl_6mo') {
    ## HDL, per 0.5 mmol/L increase
    x <- (x - 1.5) / 0.5
  } else if(quantity == 'pulse_6mo') {
    ## Heart rate, per 10 b.p.m increase
    x <- (x - 70) / 10
  } else if(quantity == 'crea_6mo') {
    ## Creatinine, per 30 μmol/L increase
    x <- (x - 60) / 30
  } else if(quantity == 'total_wbc_6mo') {
    ## White cell count, per 1.5 109/L increase
    x <- (x - 7.5) / 1.5
  } else if(quantity == 'haemoglobin_6mo') {
    ## Haemoglobin, per 1.5 g/dL increase
    x <- (x - 13.5) / 1.5
  }
  
  # Return transformed values
  x
}

caliberScale <- function(df) {
  # Return a data frame with all quantities normalised/scaled/standardised to
  # ranges etc used in Rapsomaniki et al. 2014
  
  data.frame(
    ## Time to event
    surv_time = df[, surv.time],
    ## Death/censorship
    surv_event = df[, surv.event] %in% surv.event.yes,
    ## Rescaled age
    age = sapply(df$age, ageSpline),
    ## Gender
    gender = df$gender,
    ## Most deprived quintile, yes vs. no
    most_deprived =
      df$imd_score > quantile(df$imd_score, 0.8, na.rm = TRUE),
    ### SCAD diagnosis and severity ############################################
    ## Other CHD / unstable angina / NSTEMI / STEMI vs. stable angina
    diagnosis = factorChooseFirst(factor(df$diagnosis), 'SA'),
    ## PCI in last 6 months, yes vs. no
    pci_6mo = df$pci_6mo,
    ## CABG in last 6 months, yes vs. no
    cabg_6mo = df$cabg_6mo,
    ## Previous/recurrent MI, yes vs. no
    hx_mi = df$hx_mi,
    ## Use of nitrates, yes vs. no
    long_nitrate = df$long_nitrate,
    ### CVD risk factors #######################################################
    ## Ex-smoker vs. never / Current smoker vs. never
    smokstatus = factorChooseFirst(factor(df$smokstatus), 'Non'),
    ## Hypertension, present vs. absent
    hypertension = df$hypertension,
    ## Diabetes mellitus, present vs. absent
    diabetes_logical = df$diabetes != 'No diabetes',
    ## Total cholesterol, per 1 mmol/L increase
    total_chol_6mo = df$total_chol_6mo - 5,
    ## HDL, per 0.5 mmol/L increase
    hdl_6mo = (df$hdl_6mo - 1.5) / 0.5,
    ### CVD co-morbidities #####################################################
    ## Heart failure, present vs. absent
    heart_failure = df$heart_failure,
    ## Peripheral arterial disease, present vs. absent
    pad = df$pad,
    ## Atrial fibrillation, present vs. absent
    hx_af = df$hx_af,
    ## Stroke, present vs. absent
    hx_stroke = df$hx_stroke,
    ### Non-CVD comorbidities ##################################################
    ## Chronic kidney disease, present vs. absent
    hx_renal = df$hx_renal,
    ## Chronic obstructive pulmonary disease, present vs. absent
    hx_copd = df$hx_copd,
    ## Cancer, present vs. absent
    hx_cancer = df$hx_cancer,
    ## Chronic liver disease, present vs. absent
    hx_liver = df$hx_liver,
    ### Psychosocial characteristics ###########################################
    ## Depression at diagnosis, present vs. absent
    hx_depression = df$hx_depression,
    ## Anxiety at diagnosis, present vs. absent
    hx_anxiety = df$hx_anxiety,
    ### Biomarkers #############################################################
    ## Heart rate, per 10 b.p.m increase
    pulse_6mo = (df$pulse_6mo - 70) / 10,
    ## Creatinine, per 30 μmol/L increase
    crea_6mo = (df$crea_6mo - 60) / 30,
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo = (df$total_wbc_6mo - 7.5) / 1.5,
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo = (df$haemoglobin_6mo - 13.5) / 1.5
  )
}