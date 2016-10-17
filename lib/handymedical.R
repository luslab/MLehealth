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
      paste0('Non-unique breaks and discarding of duplicates has been disabled. ',
      'Please choose different quantiles to split at.')
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