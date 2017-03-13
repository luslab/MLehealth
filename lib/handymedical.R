source('handy.R')

requirePlus(
  'foreach',
  #'CALIBERdatamanage',
  #'CALIBERcodelists',
  #'CALIBERlookups',
  'plyr',
  'ggplot2',
  'utils',
  'reshape2',
  'GGally',
  'psych',
  'bnlearn',
  'rms',
  'survival',
  'ranger',
  'randomForestSRC',
  'e1071',
  'Rgraphviz',
  'data.table',
  'boot',
  install = FALSE
)

readMedicalData <- function(filenames, col.keep, col.class) {
  # read the file(s) into a data table
  df <- foreach(filename = filenames, .combine = 'rbind') %do% {
    fread(
      filename,
      sep = ',',
      select = col.keep,
      #colClasses = col.class,
      data.table = FALSE
    )
  }
  # go through altering the classes of the columns where specified
  for(i in 1:ncol(df)) {
    if(col.class[i] == 'factor') {
      df[,i] <- factor(df[,i])
    } else if(col.class[i] == 'date') {
      df[,i] <- as.Date(df[,i])
    }
  }
  
  # return the data
  df
}

getQuantiles <- function(x, probs, duplicate.discard = TRUE) {
  breaks <- quantile(x, probs, na.rm = TRUE)
  if(duplicate.discard) {
    breaks <- unique(breaks)
  } else if (sum(duplicated(breaks))) {
    stop(
      'Non-unique breaks and discarding of duplicates has been disabled. ',
      'Please choose different quantiles to split at.'
    )
  }
  breaks
}

binByQuantile <- function(x, probs, duplicate.discard = TRUE) {
  # discretises data by binning a vector of values x into quantile-based bins
  # defined by probs
  breaks <- getQuantiles(x, probs, duplicate.discard = duplicate.discard)
  factorNAfix(
    cut(
      x,
      breaks,
      include.lowest = TRUE
    ),
    NAval = 'missing'
  )
}

binByAbs <- function(x, breaks) {
  # discretises data by binning given absolute values of breaks, and includes
  # the minimum and maximum values so all data are included
  factorNAfix(
    cut(
      x,
      c(min(x, na.rm = TRUE), breaks, max(x, na.rm = TRUE)),
      include.lowest = TRUE
    ),
    NAval = 'missing'
  )
}

missingToBig <- function(x) {
  # Removes missing values and gives them an extreme (high) value
  
  # Get a value which is definitely far higher than the maximum value, and is
  # easy for a human to spot
  max.x <- max(x, na.rm = TRUE)
  # If the max is less than zero, zero will do
  if(max.x < 0) {
    really.big.value <- 0
  # If the max is zero, then 100 is easy to spot
  } else if(max.x == 0) {
    really.big.value <- 100
  # Finally, if the max value is positive, choose one at least 10x bigger
  } else {
    really.big.value <- 10*10^ceiling(log10(max.x))
  }
  
  # Set the NA values to that number and return
  NA2val(x, really.big.value)
}

prepData <- function(
  # surv.event cannot be 'surv_event' or will break later!
  # The fraction of the data to use as the test set (1 - this will be used as
  # the training set)
  # Default quantile boundaries for discretisation
  df, predictors, process.settings, col.time, col.event, event.yes = NA,
  default.quantiles = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
  extra.fun = NULL, random.seed = NA, NAval = 'missing', n.keep = NA
) {
  # If a random seed was provided, set it
  if(!is.na(random.seed)) set.seed(random.seed)
  
  # If we only want n.keep of the data, might as well throw it out now to make
  # all the steps from here on faster...
  if(!is.na(n.keep)) {
    # Keep rows at random to avoid bias
    df <- sample.df(df, n.keep)
  } else {
    # If there was no n.keep, we should still randomise the rows for consistency
    df <- sample.df(df, nrow(df))
  }
  
  # Add event column to predictors to create full column list
  columns <- c(col.time, col.event, predictors)
  
  # Only include the columns we actually need, and don't include any which
  # aren't in the data frame because it's possible that some predictors may be
  # calculated later, eg during extra.fun
  df <- df[, columns[columns %in% names(df)]]
  
  # Go through per predictor and process them
  for(col.name in predictors[predictors %in% names(df)]){
    # If we have a specific way to process this column, let's do it!
    if(col.name %in% process.settings$var) {
      j <- match(col.name, process.settings$var)

      # Processing method being NA means leave it alone...
      if(!is.na(process.settings$method[j])) {
        # ...so, if not NA, use the function provided
        process.fun <- match.fun(process.settings$method[j])
        
        # If there are no process settings for this, just call the function
        if(is.na(process.settings$settings[[j]])) {
          df[, col.name] <- process.fun(df[, col.name])
        # Otherwise, call the function with settings
        } else {
          df[, col.name] <-
            process.fun(
              df[, col.name],
              process.settings$settings[[j]]
            )
        }
      }
    # Otherwise, no specific processing specified, so perform defaults
    } else {
      # If it's a character column, make it a factor
      if(is.character(df[, col.name])) {
        df[, col.name] <- factor(df[, col.name])
      }
      # Then, if there are any NAs, go through and make them a level of their own
      if(is.factor(df[, col.name]) & anyNA(df[, col.name])){
        df[, col.name] <-
          factorNAfix(df[, col.name], NAval = NAval)
      }
      # If it's numerical, then it needs discretising
      if(class(df[,col.name]) %in% c('numeric', 'integer')) {
          df[,col.name] <-
              binByQuantile(df[,col.name], default.quantiles)
      # Finally, if it's logical, turn it into a two-level factor
      } else if(class(df[,col.name]) == 'logical') {
        df[,col.name] <- factor(df[,col.name])
        # If there are missing values, fix them
        if(anyNA(df[, col.name])) {
          factorNAfix(df[, col.name], NAval = NAval)
        }
      }
    }
  }
  
  # Rename the survival time column
  names(df)[names(df) == col.time] <- 'surv_time'
  # Create a column denoting censorship or otherwise of events
  df$surv_event <- df[, col.event] %in% event.yes
  # Remove the event column so we don't use it as a covariate later
  df[, col.event] <- NULL
  
  # If there's any more preprocessing to do, do it now!
  if(!is.null(extra.fun)) {
    df <- extra.fun(df)
  }
  
  # Return prepped data
  df
}

prepCoxMissing <- function(
  df, missing.cols = NA, missing.suffix = '_missing', NAval = 'missing'
){
  # If a list of columns which may contain missing data wasn't provided, then
  # find those columns which do, in fact, contain missing data.
  # (Check length == 1 or gives a warning if testing a vector.)
  if(length(missing.cols) == 1) {
    if(is.na(missing.cols)) {
      missing.cols <- c()
      for(col.name in names(df)) {
        if(sum(is.na(df[, col.name])) > 0) {
          missing.cols <- c(missing.cols, col.name)
        }
      }
    }
  }
  
  # Go through missing.cols, processing appropriately for data type
  for(col.name in missing.cols) {
    
    # If it's a factor, simply create a new level for missing values
    if(is.factor(df[, col.name])) {
      # If it's a factor, NAs can be their own level
      df[, col.name] <-
        factorNAfix(df[, col.name], NAval = NAval)
      
    } else {
      # If it isn't a factor, first create a column designating missing values
      df[, paste0(col.name, missing.suffix)] <-
        is.na(df[, col.name])
      
      # Then, deal with the actual values to remove their effect, depending on
      # variable type
      if(is.logical(df[, col.name])) {
        # Set the NA values to baseline so they don't contribute to the model
        df[is.na(COHORT.scaled[, col.name]), col.name] <- FALSE
      } else {
        # Set the NA values to baseline so they don't contribute to the model
        df[is.na(df[, col.name]), col.name] <- 0
      }
      
    }
  }
  
  df
}

plotConfusionMatrix <- function(truth, prediction, title = NA) {
  confusion.matrix <- table(truth, prediction)
  
  # normalise by columns, ie predictions sum to probability 1
  confusion.matrix.n <- sweep(confusion.matrix, 1, rowSums(confusion.matrix),
                              FUN="/")
  
  confusion.matrix.n <- melt(confusion.matrix.n)
  
  confusion.matrix.plot <-
    ggplot(confusion.matrix.n,
           aes(x=truth,
               y=prediction,
               fill=value)) +
    geom_tile()
  
  if(!is.na(title)) {
    confusion.matrix.plot <-
      confusion.matrix.plot + ggtitle(title)
  }
  
  print(confusion.matrix.plot)
  
  # return the raw confusion matrix
  confusion.matrix
}

convertFactorsToBinaryColumns<- function(data_frame){
  covariates=colnames(data_frame)
  return(model.matrix( 
    formula(paste0('~',paste0(covariates,collapse='+'))), 
    data = data_frame
  )[,-1]
  )
}

getTopStates <- function(df, n = 10) {
  # Given a data frame, find the top unique 'states', ie collections of common
  # values, and return a vector of which rows belong to each state, and NA for
  # those which aren't in the top n states.
  # df = a data frame
  # n = the number of top states
  all.states <- do.call(paste, df)
  top.states <-
    head(
      sort(
        table(all.states),
        decreasing = TRUE
      ),
      n
    )
  factor(all.states, levels=names(top.states))
}

cvFolds <- function(n.data, n.folds = 3) {
  # Return a list of n.folds vectors containing the numbers 1:n.data, scrambled
  # randomly.
  split(
    sample(1:n.data),
    ceiling((1:n.data)/(n.data/n.folds))
  )
}

calcRFCIndex <- function(model.fit, df, risk.time) {
  # Calculate the C-index for model.fit based on data in df, using the value of
  # the predicted survival curve at risk.time as a predictor
  predictions <- predict(model.fit, df)
  risk.bin <- which.min(abs(predictions$unique.death.times - risk.time))
  # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
  rf.risk.proxy <- 1 - predictions$survival[, risk.bin]
  # Return C-index
  as.numeric(
    survConcordance(
      Surv(surv_time, surv_event) ~ rf.risk.proxy,
      df
    )$concordance
  )
}

calcCoxCIndex <- function(model.fit, df) {
  # Calculate the C-index for a Cox proportional hazards model on data in df
  predictions <- predict(model.fit, df)
  as.numeric(
    survConcordance(
      Surv(surv_time, surv_event) ~ predictions,
      df
    )$concordance
  )
}

cIndex <- function(model.fit, df, model.type = 'cph', risk.time = 5) {
  # Calculate the C-index for a Cox proportional hazards model on data in df
  predictions <- predict(model.fit, df)
  # If we're dealing with a ranger model, then we need to get a proxy for risk
  if(model.type == 'ranger') {
    risk.bin <- which.min(abs(predictions$unique.death.times - risk.time))
    # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
    predictions <- 1 - predictions$survival[, risk.bin]
  } else if(model.type == 'rfsrc') {
    # If we're dealing with a randomForestSRC model, extract the 'predicted' var
    predictions <- predictions$predicted
  } else if(model.type == 'survreg') {
    # survreg type models give larger numbers for longer survival...this is a
    # hack to make this return C-indices which make sense!
    predictions <- max(predictions) - predictions
  }
  as.numeric(
    survConcordance(
      Surv(surv_time, surv_event) ~ predictions,
      df
    )$concordance
  )
}

cphCoeffs <- function(cph.model, df) {
  # Get the survival variables from the coefficients in the model, which are
  # recorded as variable=level.
  # The data frame is needed to provide the baseline levels, which don't seem to
  # be stored in the Cox model object...
  surv.vars <-
    unique(sapply(strsplit(names(cph.model$coefficients), '='), firstElement))
  # get all the possible levels by looking through the data for the levels of each
  # variable identified above
  surv.vars.levels <-
    sapply(surv.vars, function(x){levels(df[,x])})
  # create a data frame to populate with the beta values
  surv.vars.df <- data.frame(
    var   = rep(surv.vars, sapply(surv.vars.levels, length)),
    level = unlist(surv.vars.levels),
    beta  = 0 # betas are zero for all baselines, so make that the default value
  )
  # go through each coefficient in the survival fit...
  for(i in 1:length(cph.model$coefficients)) {
    # ...extract the variable name and level...
    needle <- strsplit(names(cph.model$coefficients[i]), '=')[[1]]
    surv.vars.df[
      # choose the beta value of the variable at that level
      surv.vars.df$var == needle[1] & surv.vars.df$level == needle[2], 'beta'
      ] <- cph.model$coefficients[i] # and set it to the relevant coefficient
  }
  surv.vars.df
}

# Create per-patient survival curves from a data frame and a Cox model
cphSurvivalCurves <-
  function(
    df,
    surv.model,
    surv.times = max(df$surv_time)*seq(0, 1, length.out = 100)
  ) {
    # return a large, melted data frame of the relevant curves
    data.frame(
      #anonpatid = rep(df$anonpatid, each = length(surv.times)),
      id = rep(1:nrow(df), each = length(surv.times)),
      surv_time = rep(df$surv_time, each = length(surv.times)),
      surv_event = rep(df$surv_event, each = length(surv.times)),
      t = rep(surv.times, times = nrow(df)),
      s = 
        c(
          t(
            survest(surv.model,
                    newdata=df,
                    times=surv.times,
                    conf.int = FALSE # we don't want confidence intervals
            )$surv
          )
        )
    )
  }

# Create per-patient survival curves from a data frame and a random forest
rfSurvivalCurves <-
  function(
    df,
    predict.rf
  ) {
    surv.times <- predict.rf$unique.death.times
    # return a large, melted data frame of the relevant curves
    data.frame(
      #anonpatid = rep(df$anonpatid, each = length(surv.times)),
      id = rep(1:nrow(df), each = length(surv.times)),
      surv_time = rep(df$surv_time, each = length(surv.times)),
      surv_event = rep(df$surv_event, each = length(surv.times)),
      t = rep(surv.times, times = nrow(df)),
      s = c(t(predict.rf$survival))
    )
  }

getSurvCurves <- function(
  df,
  predictions,
  model.type = 'cph',
  surv.times = max(df$surv_time)*seq(0, 1, length.out = 100)
) {
  if(model.type == 'cph') {
    # return a large, melted data frame of the relevant curves
    data.frame(
      #anonpatid = rep(df$anonpatid, each = length(surv.times)),
      id = rep(1:nrow(df), each = length(surv.times)),
      surv_time = rep(df$surv_time, each = length(surv.times)),
      surv_event = rep(df$surv_event, each = length(surv.times)),
      t = rep(surv.times, times = nrow(df)),
      s = 
        c(
          t(
            survest(surv.model,
                    newdata=df,
                    times=surv.times,
                    conf.int = FALSE # we don't want confidence intervals
            )$surv
          )
        )
    )
  } else if(model.type == 'ranger') {
    surv.times <- predictions$unique.death.times
    # return a large, melted data frame of the relevant curves
    data.frame(
      #anonpatid = rep(df$anonpatid, each = length(surv.times)),
      id = rep(1:nrow(df), each = length(surv.times)),
      surv_time = rep(df$surv_time, each = length(surv.times)),
      surv_event = rep(df$surv_event, each = length(surv.times)),
      t = rep(surv.times, times = nrow(df)),
      s = c(t(predictions$survival))
    )
  }  else if(model.type == 'rfsrc') {
    surv.times <- predictions$time.interest
    # return a large, melted data frame of the relevant curves
    data.frame(
      #anonpatid = rep(df$anonpatid, each = length(surv.times)),
      id = rep(1:nrow(df), each = length(surv.times)),
      surv_time = rep(df$surv_time, each = length(surv.times)),
      surv_event = rep(df$surv_event, each = length(surv.times)),
      t = rep(surv.times, times = nrow(df)),
      s = c(t(predictions$survival))
    )
  }
}

survivalFit <- function(
  predict.vars, df, model.type = 'cph',
  n.trees = 500, split.rule = 'logrank', n.threads = 1, tod.round = 0.1,
  bootstraps = 200, ...
) {
  
  # Depending on model.type, change the name of the variable for survival time
  if(model.type %in% c('cph', 'survreg', 'survreg.boot')) {
    # Cox models can use straight death time
    surv.time = 'surv_time'
  } else {
    # Random forests need to use rounded death time
    surv.time = 'surv_time_round'
    
    df$surv_time_round <-
      round_any(df$surv_time, tod.round)
  }
  
  # Create a survival formula with the provided variable names...
  surv.formula <-
    formula(
      paste0(
        # Survival object made in-formula
        'Surv(', surv.time,', surv_event) ~ ',
        # Predictor variables then make up the other side
        paste(predict.vars, collapse = '+')
      )
    )
  
  # Then, perform the relevant type of fit depending on the model type requested 
  if(model.type == 'cph') {
    return(
      cph(surv.formula, df, surv = TRUE)
    )
  } else if(model.type == 'survreg') {
    return(
      survreg(surv.formula, df, dist = 'exponential')
    )
  } else if(model.type == 'survreg.boot') {
    return(
      boot(
        formula = surv.formula,
        data = df,
        statistic = bootstrapFit,
        fit.fun = survreg,
        R = bootstraps,
        dist = 'exponential',
        ...
      )
    )
  } else if(model.type == 'ranger') {
    return(
      ranger(
        surv.formula,
        df,
        num.trees = n.trees,
        splitrule = split.rule,
        num.threads = n.threads,
        ...
      )
    )
  } else if(model.type == 'rfsrc') {
    return(
      rfsrc(
        surv.formula,
        df,
        ntree = n.trees
      )
    )
  }
}

bootstrapFit <- function(formula, data, indices, fit.fun) {
  # Wrapper function to pass generic fitting functions to boot for
  # bootstrapping. This is actually called by boot, so much of this isn't
  # specified manually.
  #
  # Args:
  #   formula: The formula to fit with, given by the formula argument in boot.
  #      data: The data to fit, given by the data argument in boot.
  #   indices: Used internally by boot to select each bootstrap sample.
  #   fit.fun: The function you'd like to use to fit with, eg lm, cph, survreg,
  #            etc. You pass this to boot as part of its ... arguments, so
  #            provide it as fit.fun. It must return something sensible when
  #            acted on by the coef function.
  #       ...: Other arguments to your fitting function. This is now a nested
  #            ..., since you'll put these hypothetical arguments in boot's ...
  #            to pass here, to pass to your fitting function.
  #
  # Returns:
  #   The coefficients of the fit, which are then aggregated over multiple
  #   passes by boot to construct estimates of variation in parameters.
  
  d <- data[indices,]
  fit <- fit.fun(formula, data = d, ...)
  return(coef(fit))
}

bootstrapFitSurvreg <- function(formula, data, indices, test.data) {
  # Wrapper function to pass a survreg fit with c-index calculations to boot.
  
  d <- data[indices,]
  fit <- survreg(formula, data = d, dist = 'exponential')
  # Return fit coefficients, c-index on training data, c-index on test data
  return(
    c(
      coef(fit),
      c.train = cIndex(fit, data, model.type = 'survreg'),
      c.test = cIndex(fit, test.data, model.type = 'survreg')
    )
  )
}

bootStats <- function(bootfit, transform = identity) {
  # Return a data frame with the statistics from a bootstrapped fit
  #
  # Args:
  #    bootfit: A boot object.
  #  transform: Optional transform for the statistics, defaults to identity, ie
  #             leave the values as they are. Useful if you want the value and
  #             variance of the exp(statistic), etc.
  return(
    data.frame(
      val  = transform(bootfit$t0),
      err  = sqrt(apply(transform(bootfit$t), 2, var))
    )
  )
}

negExp <- function(x) {
  exp(-x)
}