#' # Random forest cross-validation
#' 
#' Random forests don't really have any tunable parameters beyond the choice of 
#' splitting criterion. However, our random forests do: because we're looking
#' at continuous data with missing values, we've decided to bin the data and
#' thereby make it categorical, with missing as a final, additional category.
#' It's also necessary to bin the times of death in order to make the problem
#' computationally tractable.
#' 
#' What difference do these parameters make?

#+ setup, message=FALSE

require(reshape2)
require(ggplot2)

rf.cv <- read.csv('../../output/rf-crossvalidation-try2.csv')

#' First, create a data frame of averages.

rf.cv.average <-
  aggregate(
    cbind(
      c.index.rf.train, c.index.rf.val, time.learn.rf, 
      c.index.cph.train, c.index.cph.val
    ) ~ discretise.bins + tod.round,
    data = rf.cv,
    mean
  )

#' ## Overview plot
#' 
#' This heatmap shows the effect of both number of bins used to discretise and
#' accuracy of rounding of time of death. Clearly, using more bins improves
#' C-index on the validation set, but increasing resolution of time of death
#' makes very little difference.

ggplot(
  rf.cv.average,
  aes(x = discretise.bins, y = factor(tod.round), fill = c.index.rf.val)
) +
  geom_tile() +
  geom_text(aes(label = round(c.index.rf.val, 3)))

#' ## Number of bins
#' 
#' Looking more carefully at number of bins, this is a scatterplot of C-index
#' asa  function of number of bins with time of death rounded to a constant 0.1
#' years. The red points represent the score on the training set, and the blue
#' on the validation set.
#' 
#' The random forests certainly overfit to a fair extent. Not only are the
#' C-indices on the training set higher, but they increase more or less
#' monotonically with increasing numbers of bins. By contrast, the C-index on
#' the validation set has probably plateaued by around the five-bin mark.

ggplot(subset(rf.cv, tod.round == 0.1), aes(x = discretise.bins)) +
  geom_point(aes(y = c.index.rf.val), colour = 'blue') +
  geom_point(aes(y = c.index.rf.train), colour = 'red')

#' ## Time of death rounding
#' 
#' Looking at situations where all variables are discretised to 5 bins, it's
#' pretty obvious that rounding of the time of death makes very little
#' difference. Especially when you can see the variability between
#' cross-validation folds...

ggplot(subset(rf.cv, discretise.bins == 5), aes(x = tod.round)) +
  geom_point(aes(y = c.index.rf.val), colour = 'blue') +
  geom_point(aes(y = c.index.rf.train), colour = 'red')



#' ## Cox models
#' 
#' The cox models show much the same trend as the random forests, but less
#' noisy. Time of death rounding makes almost no difference-and besides, this
#' actually isn't necessary for a Cox model anyway!

ggplot(
  rf.cv.average,
  aes(x = discretise.bins, y = factor(tod.round), fill = c.index.cph.val)
) +
  geom_tile() +
  geom_text(aes(label = round(c.index.rf.val, 3)))

#' ### Number of bins
#' 
#' Cox models seem far less prone to overfitting than random forests, and also
#' plateau, but perhaps a little later..?

ggplot(subset(rf.cv, tod.round == 0.1), aes(x = discretise.bins)) +
  geom_point(aes(y = c.index.cph.val), colour = 'blue') +
  geom_point(aes(y = c.index.cph.train), colour = 'red')

#' ### Time of death rounding
#' 
#' Clearly makes almost no difference.

ggplot(subset(rf.cv, discretise.bins == 5), aes(x = tod.round)) +
  geom_point(aes(y = c.index.cph.val), colour = 'blue') +
  geom_point(aes(y = c.index.cph.train), colour = 'red')
