source('handy.R')

requirePlus(
  'foreach',
  #'CALIBERdatamanage',
  #'CALIBERcodelists',
  'CALIBERlookups',
  'plyr',
  'dplyr',
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

missingToAverage <- function(x) {
  if(is.factor(x)) {
    # If it's a factor, replace with the most common level
    return(NA2val(x, val = modalLevel(x)))
    
  } else {
    # If it isn't a factor, replace with the median value
    return(NA2val(x, val = median(x, na.rm = TRUE)))
  }
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

missingToZero <- function(x) {
  NA2val(x, val = 0)
}

missingToSample <- function(x) {
  NA2val(x, val = samplePlus(x, replace = TRUE))
}

prepSurvCol <- function(df, col.time, col.event, event.yes) {
  # Rename the survival time column
  names(df)[names(df) == col.time] <- 'surv_time'
  # Create a column denoting censorship or otherwise of events
  df$surv_event <- df[, col.event] %in% event.yes
  # Remove the event column so we don't use it as a covariate later
  df[, col.event] <- NULL
  
  df
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
        if(isExactlyNA(process.settings$settings[[j]])) {
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
  
  # Sort out the time to event and event class columns
  df <- prepSurvCol(df, col.time, col.event, event.yes)
  
  # If there's any more preprocessing to do, do it now!
  if(!is.null(extra.fun)) {
    df <- extra.fun(df)
  }
  
  # Return prepped data
  df
}

prepCoxMissing <- function(
  df, missing.cols = NA, missingReplace = missingToZero,
  missing.suffix = '_missing', NAval = 'missing'
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
      df[, paste0(col.name, missing.suffix)] <- is.na(df[, col.name])
      
      # If we want to replace the missing values...
      if(!isExactlyNA(missingReplace)) {
      
        # Then, deal with the actual values, depending on variable type
        if(is.logical(df[, col.name])) {
          # Set the NA values to baseline so they don't contribute to the model
          df[is.na(COHORT.scaled[, col.name]), col.name] <- FALSE
        } else {
          # Set the NA values to the desired value, eg 0 (ie baseline in a Cox
          # model with missingToZero), missingToMedian, missingToBig, etc...
          df[, col.name] <- missingReplace(df[, col.name])
        }
      }
      
    }
  }
  
  df
}

medianImpute <- function(df, missing.cols = NA) {
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
    df[, col.name] <- missingToAverage(df[, col.name])
  }
  
  df
}

modalLevel <- function(x) {
  # Return the name of the most common level in a factor x
  tt <- table(x)
  names(tt[which.max(tt)])
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

modelType <- function(model.fit) {
  # Take a model fit and return a string representing its type so as to deal
  # with it correctly
  
  # rfsrc for some reason has multiple classes associated with its fit objects
  if('rfsrc' %in% class(model.fit)) {
    return('rfsrc')
  # Other models are more sensible, and simply returning the class will do
  } else {
    return(class(model.fit))
  }
}

cIndex <- function(model.fit, df, risk.time = 5, tod.round = 0.1, ...) {
  if(modelType(model.fit) == 'rfsrc') {
    # rfsrc throws an error unless the y-values in the provided data are
    # identical to those used to train the model, so recreate the rounded ones..
    df$surv_time_round <-
      round_any(df$surv_time, tod.round)
    # This means we need to use surv_time_round in the formula
    surv.time <- 'surv_time_round'
  } else {
    # Otherwise, our survival time variable is just surv_time
    surv.time <- 'surv_time'
  }
  
  # Calculate the C-index for a Cox proportional hazards model on data in df
  
  # First, get some risks, or values proportional to them
  risk <- getRisk(model.fit, df, ...)
  
  # Then, get the C-index and, since we don't probably want to do any further
  # work with it, simply return the numerical value of the index itself.
  as.numeric(
    survConcordance(
      as.formula(paste0('Surv(', surv.time, ', surv_event) ~ risk')),
      df
    )$concordance
  )
}

generalVarImp <-
  function(
    model.fit, df, vars = NA, risk.time = 5, tod.round = 0.1, ...
    ) {
  baseline.c.index <- cIndex(model.fit, df, risk.time, tod.round, ...)
  
  # If no variables were passed, let's do it on all of the variables
  if(isExactlyNA(vars)) {
    vars <- names(model.fit$xvar)
  }
  
  var.imp <- data.frame(
     var = vars,
     var.imp = NA,
     stringsAsFactors = FALSE
  )
  for(i in 1:nrow(var.imp)) { 
    # Make a new, temporary data frame
    df2 <- df
    # Permute values of the sample in question
    df2[, var.imp[i, 'var']] <- sample(df[, var.imp[i, 'var']], replace = TRUE)
    # Calculate the C-index based on the permuted data
    var.c.index <- cIndex(model.fit, df2, risk.time, tod.round,
                          ...)
    var.imp[i, 'var.imp'] <- baseline.c.index - var.c.index
  }
  
  # Return the data frame of variable importances
  var.imp
}

modelFactorLevelName <- function(factor.name, level.name, model.type) {
  if(model.type == 'cph') {
    # factor=Level
    return(paste0(factor.name, '=', level.name))
  } else if(model.type == 'survreg') {
    # factorLevel
    return(paste0(factor.name, level.name))
  } else if(model.type == 'boot.survreg') {
    # factorLevel
    return(paste0(factor.name, level.name))
  }
}

cphCoeffs <- function(cph.model, df, surv.predict, model.type = 'cph') {
  # Depending on the model type, get a vector of the Cox coefficient names...
  if(model.type == 'cph') {
    coeff.names <- names(cph.model$coefficients)
    coeff.vals  <- cph.model$coefficients
  } else {
    # Otherwise, it will come as a data frame of some kind
    coeff.names <- rownames(cph.model)
    coeff.vals  <- cph.model$val
    coeff.lower  <- cph.model$lower
    coeff.upper  <- cph.model$upper
  }
  
  # Get the names and levels from each of the factors used to create the
  # survival model. Models by cph are good enough to separate with = (ie 
  # factor=level), but this is not universal so it's a more general solution to
  # create these coefficient names from the data in a per-model-type way.
  surv.vars.levels <- sapply(surv.predict, function(x){levels(df[,x])})
  surv.vars.df <- 
    data.frame(
      var   = rep(surv.predict, unlist(sapply(surv.vars.levels, length))),
      level = unlist(surv.vars.levels),
      val   = 0, # betas are zero for all baselines so make that the default val
      err   = 0, # uncertainty is zero for a baseline too!
      stringsAsFactors = FALSE
    )
  # go through each coefficient in the survival fit...
  for(i in 1:nrow(surv.vars.df)) {
    # ...create the factor/level coefficient name...
    needle <-
      modelFactorLevelName(
        surv.vars.df[i, 'var'], surv.vars.df[i, 'level'],
        model.type
      )
    # ...find where in the coefficients that name occurs...
    if(sum(coeff.names == needle) > 0) {
      needle.i <- which(coeff.names == needle)
      # ...and set the relevant value and error
      surv.vars.df[i, 'val'] <- coeff.vals[needle.i]
      surv.vars.df[i, 'lower'] <- coeff.lower[needle.i]
      surv.vars.df[i, 'upper'] <- coeff.upper[needle.i]
    }
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
        dist = 'exponential'
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
    # rfsrc, if you installed it correctly, controls threading by changing an
    # environment variable
    options(rf.cores = n.threads)
    
    # Fit and return
    return(
      rfsrc(
        surv.formula,
        df,
        ntree = n.trees,
        ...
      )
    )
  }
}

survivalBootstrap <- function(
  predict.vars, df, df.test, model.type = 'survreg',
  n.trees = 500, split.rule = 'logrank', n.threads = 1, tod.round = 0.1,
  bootstraps = 200, ...
) {
  
  # Depending on model.type, change the name of the variable for survival time
  if(model.type %in% c('survreg')) {
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
    stop('model.type cph not yet implemented')
  } else if(model.type == 'survreg') {
    return(
      boot(
        formula = surv.formula,
        data = df,
        statistic = bootstrapFitSurvreg,
        R = bootstraps,
        parallel = 'multicore',
        ncpus = n.threads,
        test.data = df.test
      )
    )
  } else if(model.type == 'ranger') {
    stop('model.type ranger not yet implemented')
  } else if(model.type == 'rfsrc') {
    # rfsrc, if you installed it correctly, controls threading by changing an
    # environment variable
    options(rf.cores = n.threads)
    
    return(
      boot(
        formula = surv.formula,
        data = df,
        statistic = bootstrapFitRfsrc,
        R = bootstraps,
        parallel = 'no',
        ncpus = 1, # disable parallelism because rfsrc can be run in parallel
        n.trees = n.trees,
        test.data = df.test,
        
        ...
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
      c.train = cIndex(fit, data),
      c.test = cIndex(fit, test.data),
      calibration.score =
        # Only the area, let's not return the standard error for now...
        calibrationScore(calibrationTable(fit, test.data))$area
    )
  )
}

bootstrapFitSurvregMice <- function(formula, data, indices, test.data) {
  # Wrapper function to pass a survreg fit with c-index calculations to boot.
  
  d <- data[indices,]
  fit <- survreg(formula, data = d, dist = 'exponential')
  # Return fit coefficients, c-index on training data, c-index on test data
  return(
    c(
      fit$qbar, # Bombined coefficients from the imputed runs are stored here
      c.train = cIndex(fit, data, model.type = 'survreg'),
      c.test = cIndex(fit, test.data, model.type = 'survreg')
    )
  )
}

bootstrapFitRfsrc <- function(formula, data, indices, n.trees, test.data, ...) {
  # Wrapper function to pass an rfsrc fit with c-index calculations to boot.
  
  # There used to be an rf-cores declaration here, but keep it as-is and run
  # bootstrapping as a single-threaded process...
  
  d <- data[indices,]
  fit <-  
    rfsrc(
      formula,
      data,
      ntree = n.trees,
      ...
    )
  
  # Check the model calibration on the test set
  calibration.table <- calibrationTable(fit, test.data, ...)
  calibration.score <- calibrationScore(calibration.table, curve = FALSE)
  
  # Return fit coefficients, c-index on training data, c-index on test data
  return(
    c(
      c.test = cIndex(fit, test.data, ...),
      calibration.score = calibration.score$area
    )
  )
}

bootStats <- function(bootfit, uncertainty = 'sd', transform = identity) {
  # Return a data frame with the statistics from a bootstrapped fit
  #
  # Args:
  #     bootfit: A boot object.
  # uncertainty: Function to use for returning uncertainty, defaulting to 'sd'
  #              which returns the standard deviation.
  #   transform: Optional transform for the statistics, defaults to identity, ie
  #              leave the values as they are. Useful if you want the value and
  #              variance of the exp(statistic), etc.
  #
  
  if(uncertainty == 'sd'){
    return(
      data.frame(
        val  = transform(bootfit$t0),
        err  = apply(transform(bootfit$t), 2, sd)
      )
    )
  } else if(uncertainty == '95ci') {
    ci <- apply(transform(bootfit$t), 2, quantile, probs = c(0.025, 0.5, 0.975))
    return(
      data.frame(
        val  = t(ci)[, 2],
        lower = t(ci)[, 1],
        upper = t(ci)[, 3]
      )
    )
  } else {
    stop("Unknown value '", uncertainty, "' for uncertainty parameter.")
  }
}

bootMIStats <- function(boot.mi, uncertainty = '95ci', transform = identity) {
  # Return a data frame with the statistics from a bootstrapped fit
  #
  # Args:
  #     bootfit: A boot object.
  # uncertainty: Function to use for returning uncertainty, defaulting to 'sd'
  #              which returns the standard deviation.
  #   transform: Optional transform for the statistics, defaults to identity, ie
  #              leave the values as they are. Useful if you want the value and
  #              variance of the exp(statistic), etc.
  #
  
  boot.mi.combined <-
    do.call(
      # rbind together...
      rbind,
      # ...a list of matrices of bootstrap estimates extracted from the list of
      # bootstrap fits
      lapply(boot.mi, function(x){x$t})
    )
 
  if(uncertainty == 'sd'){
    return(
      data.frame(
        val  = apply(transform(boot.mi.combined), 2, mean),
        err  = apply(transform(boot.mi.combined), 2, sd),
        row.names = names(boot.mi[[1]]$t0)
      )
    )
  } else if(uncertainty == '95ci') {
    ci <-
      apply(
        transform(boot.mi.combined), 2, quantile, probs = c(0.025, 0.5, 0.975)
      )
    return(
      data.frame(
        val  = t(ci)[, 2],
        lower = t(ci)[, 1],
        upper = t(ci)[, 3],
        row.names = names(boot.mi[[1]]$t0)
      )
    )
  } else {
    stop("Unknown value '", uncertainty, "' for uncertainty parameter.")
  }
}

negExp <- function(x) {
  exp(-x)
}

getRisk <- function(model.fit, df, risk.time = 5, ...) {
  
  # Make predictions for the data df based on the model model.fit
  predictions <- predict(model.fit, df, ...)
  
  # Then, for any model other than cph, they will need to be transformed in some
  # way to get a proxy for risk...
  
  # If we're dealing with a ranger model, then we need to get a proxy for risk
  # TODO: Replace model.type ifs (which now do nothing) with class-based ones
  if(modelType(model.fit) == 'ranger') {
    risk.bin <- which.min(abs(predictions$unique.death.times - risk.time))
    # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
    predictions <- 1 - predictions$survival[, risk.bin]
  } else if(modelType(model.fit) == 'rfsrc') {
    # If we're dealing with a randomForestSRC model, extract the 'predicted' var
    predictions <- predictions$predicted
  } else if(modelType(model.fit) == 'survreg') {
    # survreg type models give larger numbers for longer survival...this is a
    # hack to make this return C-indices which make sense!
    predictions <- max(predictions) - predictions
  }
  
  predictions
}

getRiskAtTime <- function(model.fit, df, risk.time = 5, ...) {
  
  # If we're dealing with a ranger model, then we need to get a proxy for risk
  if(modelType(model.fit) == 'ranger') {
    # Make predictions for the data df based on the model model.fit
    predictions <- predict(model.fit, df, ...)
    
    risk.bin <- which.min(abs(predictions$unique.death.times - risk.time))
    # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
    predictions <- 1 - predictions$survival[, risk.bin]
    
    
  } else if(modelType(model.fit) == 'rfsrc') {
    # Make predictions for the data df based on the model model.fit
    predictions <- predict(model.fit, df, ...)
    
    # If we're dealing with a randomForestSRC model, do the same as ranger but
    # with different variable names
    risk.bin <- which.min(abs(predictions$time.interest - risk.time))
    # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
    predictions <- 1 - predictions$survival[, risk.bin]
    
    
  } else if(modelType(model.fit) == 'survreg') {
    # Make predictions for the data df based on the model
    # 'quantile' returns the quantiles of risk, ie the 0.01 quantile would mean
    # 0.01 ie 1% of patients would be dead by x. Returning the risk of death
    # at a time t requires reverse-engineering this table.
    # It doesn't make sense to go to p = 1 because technically by any model the
    # 100th percentile is at infinity.
    # It's really fast, so do 1000 quantiles for accuracy. Could make this a
    # passable parameter...
    risk.quantiles <- seq(0,0.999, 0.001)
    
    predictions <-
      predict(model.fit, df, type = 'quantile', p = risk.quantiles, ...)
    
    predictions <-
      # Find the risk quantile...
      risk.quantiles[
        # ...by choosing the element corresponding to the matrix output of
        # predict, which has the same number of rows as df and a column per
        # risk.quantiles...
        apply(
          predictions,
          # ...and find the quantile closest to the risk.time being sought
          FUN = function (x) {
            which.min(abs(x - risk.time))
          },
          MARGIN = 1
        )
      ]
  } else if(modelType(model.fit) == 'survfit') {
    # For now, survfit is just a Kaplan-Meier fit, and it only deals with a
    # single variable for KM strata. For multiple strata, this would require a
    # bit or parsing to turn names like 'age=93, gender=Men' into an n-column
    # data frame.
    varname <-  substring(
      names(model.fit$strata)[1], 0,
      # Position of the = sign
      strPos('=', names(model.fit$strata)[1]) - 1
    )
    
    km.df <- data.frame(
      var = rep(
        # Chop off characters before and including = (eg 'age=') and turn into a
        # number (would also need generalising for non-numerics, eg factors)
        as.numeric(
          substring(
            names(model.fit$strata),
            # Position of the = sign
            strPos('=', names(model.fit$strata)[1]) + 1
          )
        ),
        # Repeat each number as many times as there are patients that age
        times = model.fit$strata
      ),
      time = model.fit$time,
      surv = model.fit$surv
    )
    
    risk.by.var <-
      data.frame(
        var = unique(km.df$var),
        risk = NA
      )
    
    for(var in unique(km.df$var)) {
      # If anyone with that variable value lived long enough for us to make a
      # prediction...
      if(max(km.df$time[km.df$var == var]) > risk.time) {
        # Find the first event after that point, which gives us the survival,
        # and do 1 - surv to get risk
        risk.by.var$risk[risk.by.var$var == var] <- 1-
          km.df$surv[
            # The datapoint needs to be for the correct age of patient
            km.df$var == var &
              # And pick the time which is the smallest value greater than the
              # time in which we're interested.
              km.df$time ==
              minGt(km.df$time[km.df$var == var], risk.time)
            ]
      }
    }

    # The predictions are then the risk for a given value of var
    predictions <-
      # join from pylr preserves row order
      join(
        # Slight kludge...make a data frame with one column called 'var' from
        # the var (ie variable, depending on variable!) column of the data
        data.frame(var = df[, varname]),
        risk.by.var[, c('var', 'risk')]
      )$risk
  }
  
  # However obtained, return the predictions
  predictions
}

partialEffectTable <-
  function(
    model.fit, df, variable, n.patients = 1000, max.values = 200,
    risk.time = 5, ...
  ) {
    # The number of values we look at will be either max.values, or the number
    # of unique values if that's lower. Remove NAs because they cause errors.
    n.values <- min(max.values, length(NArm(unique(df[,variable]))))
    
    # Take a sample of df, but repeat each one of those samples n.values times
    df.sample <- df[rep(sample(1:nrow(df), n.patients), each = n.values),]
    # Give each value from the original df an id, so we can keep track
    df.sample$id <- rep(1:n.patients, each = n.values)
    
    # Each individual patient from the original sample is then assigned every
    # value of the variable we're interested in exploring
    df.sample[, variable] <-
      sort(
        samplePlus(df[, variable], n.values, na.rm = TRUE, only.unique = TRUE)
      )
    
    # Use the model to make predictions
    df.sample$risk <- getRisk(model.fit, df.sample, risk.time, ...)
    
    # Use ddply to normalise the risk for each patient by the mean risk for that
    # patient across all values of variable, thus averaging out any risk offsets
    # between patients, and return that data frame.
    as.data.frame(
      df.sample %>%
        group_by(id) %>%
        mutate(risk.normalised = risk/mean(risk))
    )[, c('id', variable, 'risk.normalised')] # discard all unnecessary columns
  }

calibrationTable <- function(
  model.fit, df, risk.time = 5, tod.round = 0.1, ...
  ) {
  
  if(modelType(model.fit) == 'rfsrc') {
    # rfsrc throws an error unless the y-values in the provided data are
    # identical to those used to train the model, so recreate the rounded ones..
    df$surv_time_round <-
      round_any(df$surv_time, tod.round)
    # This means we need to use surv_time_round in the formula
    surv.time <- 'surv_time_round'
  } else {
    # Otherwise, our survival time variable is just surv_time
    surv.time <- 'surv_time'
  }
  
  # Get risk values given this model
  df$risk <- getRiskAtTime(model.fit, df, risk.time, ...)
  
  # Was there an event? Start with NA, because default is unknown (ie censored)
  df$event <- NA
  # Event before risk.time
  df$event[df$surv_event & df$surv_time <= risk.time] <- TRUE
  # Event after, whether censorship or not, means no event by risk.time
  df$event[df$surv_time > risk.time] <- FALSE
  # Otherwise, censored before risk.time, leave as NA
  
  df[, c('risk', 'event')]
}

calibrationPlot <- function(df) {
  # Convert risk to numeric, because ggplot treats logicals like categoricals
  df$event <- as.numeric(df$event)
  
  # Create a dummy data frame of censored values to plot
  censored.df <- df[is.na(df$event),]
  censored.df$event <- 0.5
  
  ggplot(df, aes(x = risk, y = event)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(position = position_jitter(w = 00, h = 0.25), alpha = 0.1) +
    geom_point(
      data = censored.df, colour = 'grey', alpha = 0.1,
      position = position_jitter(w = 00, h = 0.125)
    ) +
    geom_smooth(method = 'loess') +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))
}

calibrationScore <- function(
  calibration.table, risk.breaks = seq(0, 1, 0.001), curve = FALSE,
  extremes = TRUE
  ) {
  #
  # extremes: If set to true, this assumes predictions of 0 below 0.5, and 1
  # above 0.5, providing a worst-case estimate for cases when the prediction
  # model only provides predictions within a narrower range. This allows such
  # models to be fairly compared to others with broader predictive values.
  #
  # * Could rewrite this with the integrate built-in function
  # * Not totally sure about the standard error here...I assume just integrating
  #   the uncertainty region will result in an overestimate?
  
  
  # Fit a LOESS model to the data
  loess.curve <- loess(event ~ risk, calibration.table)
  
  # Get the bin widths, which we'll need in a bit when integrating
  risk.binwidths <- diff(risk.breaks)
  # And the midpoints of the risk bins to calculate predictions at
  risk.mids <- risk.breaks[1:(length(risk.breaks) - 1)] + risk.binwidths / 2
  
  predictions <-
    predict(loess.curve, data.frame(risk = risk.mids), se = TRUE)
  
  pred.curve <- predictions$fit
  
  if(anyNA(pred.curve)) {
    if(extremes) {
      # Get the bins where we don't have a valid prediction
      missing.risks <- risk.mids[is.na(pred.curve)]
      # And predict 0 is < 0.5, 1 if greater, for a worst-case step-function
      missing.risks <- as.numeric(missing.risks > 0.5)
      # Finally, substitute them in
      pred.curve[is.na(pred.curve)] <- missing.risks
    } else {
      # If there are missing values but extremes = FALSE, ie don't extend, then
      # issue a warning to let the user know.
      if(length(is.na(risk.mids) < 10)) {
        warning.examples <- paste(risk.mids[is.na(risk.mids)], collapse = ', ')
      } else {
        warning.examples <-
          paste(
            paste(head(risk.mids[is.na(risk.mids)], 3), collapse = ', '),
            '...',
            paste(tail(risk.mids[is.na(risk.mids)], 3), collapse = ', ')
          )
      }
      warning(
        'Some predictions (for risk bins at ', warning.examples, ') return ',
        'NA. This means calibration is being performed outside the range of ',
        'the data which may mean values are not comparable. Set extremes = ',
        'TRUE to assume worst-case predictions beyond the bounds of the ',
        'actual predictions.'
      )
    }

  }
  
  curve.area <-
    sum(
      abs(pred.curve - risk.mids) * risk.binwidths,
      na.rm = TRUE
    )
  
  # Not sure what if anything to do about NAs in the standard error, so silently
  # ignore them for now...
  curve.area.se <-
    sum(
      abs(predictions$se.fit) * risk.binwidths,
      na.rm = TRUE
    )
  
  # If the curve was requested...
  if(curve) {
    # ...return area between lines and standard error, plus the curve
    list(
      area = curve.area,
      se = curve.area.se,
      curve = predictions$fit, # Not pred.curve, because we want to see NAs
      curve.se = predictions$se.fit
    )
  } else {
    # ...otherwise, just return the summary statistic
    list(
      area = curve.area,
      se = curve.area.se
    )
  }
}

testSetIndices <- function(df, test.fraction = 1/3, random.seed = NA) {
  # Get indices for the test set in a data frame, with a random seed to make the
  # process deterministic if requested.
  
  n.data <- nrow(df)
  if(!is.na(random.seed)) set.seed(random.seed)
  
  sample.int(n.data, round(n.data * test.fraction))
}

summary2 <- function(x) {
  # Practical summary function for summarising medical records data columns
  # depending on number of unique values...
  if('data.frame' %in% class(x)) {
    lapply(x, summary2)
  } else {
    if(length(unique(x)) < 30) {
      if(length(unique(x)) < 10) {
        return(round(c(table(x))/length(x), 3)*100)
      } else {
        summ <- sort(table(x), decreasing = TRUE)
        return(
          round(
            c(
              summ[1:5],
              other = sum(summ[6:length(summ)]),
              missing = sum(is.na(x))
              # divide all by the length and turn into %
            )/length(x), 3)*100
        )
      }
    } else {
      return(
        c(
          min = min(x, na.rm = TRUE),
          max = max(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE),
          missing = round(sum(is.na(x))/length(x), 3)*100
        )
      )
    }
  }
}

lookUpDescriptions <- function(x, bnf.lookup.filename = '../../data/') {
  # Create blank columns for which dictionary a given variable comes from, its
  # code in that dictionary, and a human-readable description looked up from the
  # CALIBER tables
  
  # Make a vector to hold descriptions, fill it with x so it's a) the right
  # length and b) as a fallback
  description <- x
  code <- x
  
  # Look up ICD and OPCS codes
  relevant.rows <- startsWith(x, 'hes.icd.')
  code[relevant.rows] <- textAfter(x, 'hes.icd.')
  for(i in which(relevant.rows)) {
    # Some of these don't work, so add in an if statement to catch the error
    if(
      length(CALIBER_DICT[dict == 'icd10' & code == code[i], term]) > 0
    ){
      description[i] <-
        CALIBER_DICT[dict == 'icd10' & code == code[i], term]
    } else {
      description[i] <- 'ERROR: ICD not matched'
    }
  }
  
  relevant.rows <- startsWith(x, 'hes.opcs.')
  code[relevant.rows] <- textAfter(x, 'hes.opcs.')
  for(i in which(relevant.rows)) {
    description[i] <- CALIBER_DICT[dict == 'opcs' & code == code[i], term]
  }
  
  relevant.rows <- startsWith(x, 'clinical.history.')
  code[relevant.rows] <- textAfter(x, 'clinical.history.')
  for(i in which(relevant.rows)) {
    # Some of these don't work, so add in an if statement to catch the error
    if(
      length(CALIBER_DICT[dict == 'read' & medcode == code[i], term]) > 0
    ){
      description[i] <-
        CALIBER_DICT[dict == 'read' & medcode == code[i], term]
    } else {
      description[i] <- 'ERROR: medcode not matched'
    }
  }
  
  relevant.rows <- startsWith(x, 'clinical.values.')
  code[relevant.rows] <- textAfter(x, 'clinical.values.')
  for(i in which(relevant.rows)) {
    testtype.datatype <- strsplit(code[i], '_', fixed =TRUE)[[1]]
    description[i] <-
      paste0(
        CALIBER_ENTITY[enttype == testtype.datatype[1], description],
        ', ',
        CALIBER_ENTITY[enttype == testtype.datatype[1], testtype.datatype[2], with = FALSE]
      )
  }
  
  relevant.rows <- startsWith(x, 'bnf.')
  code[relevant.rows] <- textAfter(x, 'bnf.')
  for(i in which(relevant.rows)) {
    # Some of these don't work, so add in an if statement to catch the error
    if(
      length(CALIBER_BNFCODES[bnfcode == code[i], bnf]) > 0
    ){
      description[i] <-
        CALIBER_BNFCODES[bnfcode == code[i], bnf]
    } else {
      description[i] <- 'ERROR: BNF code not matched'
    }
  }
  
  relevant.rows <- startsWith(x, 'tests.enttype.data3.')
  code[relevant.rows] <- textAfter(x, 'tests.enttype.data3.')
  for(i in which(relevant.rows)) {
    testtype.datatype <- strsplit(code[i], '_', fixed =TRUE)[[1]]
    description[i] <-
      CALIBER_ENTITY[enttype == testtype.datatype[1], description]
  }
  
  description
}

getVarNums <- function(x, frac = 0.2) {
  # Number of iterations until there's only one variable left
  n <- -ceil(log(x)/log(1 - frac))
  unique(round(x*((1 - frac)^(0:n))))
}

percentMissing <- function(x) {
  sum(is.na(x))/length(x) * 100
}