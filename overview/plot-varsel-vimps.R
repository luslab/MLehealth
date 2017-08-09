source('../lib/handymedical.R', chdir = TRUE)
requirePlus('ggplot2', 'cowplot')

plotVarImp <- function(filename, n.max = 20) {
  # Read in variable importances
  var.imp.df <- read.csv(filename)
  
  # Make sure n.max is not greater than the number of rows, or errors will ensue
  n.max <- min(n.max, nrow(var.imp.df))
  
  # Cut off below top 20, we don't need them and they cause errors when duplicate
  # factor levels appear
  var.imp.df <- var.imp.df[order(var.imp.df$imp, decreasing = TRUE)[1:n.max], ]
  
  var.imp.df$description <- lookUpDescriptions(as.character(var.imp.df$var))
  
  # Order the levels in the factor var to make ggplot behave
  var.imp.df$description <- factor(var.imp.df$description, levels = var.imp.df$description[order(var.imp.df$imp)])
  
  ggplot(var.imp.df, aes(x = description, y = imp)) +
    geom_bar(stat='identity') +
    coord_flip()
}

plotVarImp('../../output/rf-bigdata-try11-varselrf-varimp-initial.csv')
plotVarImp('../../output/rf-bigdata-try11-varselrf-varimp-final.csv')
plotVarImp('../../output/rf-bigdata-try11-varselmiss-varimp.csv', 20)
