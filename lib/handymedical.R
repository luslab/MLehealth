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
  'e1071',
  'Rgraphviz',
  'data.table',
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

binByQuantile <- function(x, probs, duplicate.discard = TRUE) {
  # discretises data by binning a vector of values x into quantile-based bins
  # defined by probs
  breaks <- quantile(x, probs, na.rm = TRUE)
  if(duplicate.discard) {
    breaks <- unique(breaks)
  } else if (sum(duplicated(breaks))) {
    stop(
      'Non-unique breaks and discarding of duplicates has been disabled. ',
      'Please choose different quantiles to split at.'
    )
  }
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
  }
  
  # Add event column to predictors to create full column list
  columns <- c(predictors, col.time, col.event)
  
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
        
        df[, col.name] <-
          process.fun(
            df[, col.name],
            process.settings$settings[[j]]
          )
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
  
  # Create a column denoting censorship or otherwise of events
  df$surv_event <- df[, col.event] == event.yes
  
  # Remove the event column so we don't use it as a covariate later
  df <- df[, names(df) != col.event]
  
  # If there's any more preprocessing to do, do it now!
  if(!is.null(extra.fun)) {
    df <- extra.fun(df)
  }
  
  # Return prepped data
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
      Surv(time_death, surv_event) ~ rf.risk.proxy,
      df
    )$concordance
  )
}

calcCoxCIndex <- function(model.fit, df) {
  # Calculate the C-index for a Cox proportional hazards model on data in df
  predictions <- predict(model.fit, df)
  as.numeric(
    survConcordance(
      Surv(time_death, surv_event) ~ predictions,
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
    surv.times = max(df$time_death)*seq(0, 1, length.out = 100)
  ) {
    # return a large, melted data frame of the relevant curves
    data.frame(
      #anonpatid = rep(df$anonpatid, each = length(surv.times)),
      id = rep(1:nrow(df), each = length(surv.times)),
      time_death = rep(df$time_death, each = length(surv.times)),
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
      time_death = rep(df$time_death, each = length(surv.times)),
      surv_event = rep(df$surv_event, each = length(surv.times)),
      t = rep(surv.times, times = nrow(df)),
      s = c(t(predict.rf$survival))
    )
  }