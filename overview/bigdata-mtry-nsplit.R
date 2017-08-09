source('../lib/handy.R')
requirePlus('ggplot2', 'cowplot')

# Read in two cross-validation data files
cv.performance <-
  read.csv('../../output/rf-bigdata-try7-ALL-cv-largemtry-calibration.csv')
cv.performance <-
  rbind(
    cv.performance,
    read.csv('../../output/rf-bigdata-try7-ALL-cv-smallmtry-calibration.csv')
  )

# Read in overall model performance
models.performance <- readTablePlus('../../output/models-performance.tsv')



cv.performance.avg <-
  aggregate(
    c.index.val ~ n.splits + m.try,
    data = cv.performance,
    mean
  )

ggplot(cv.performance.avg, aes(x = n.splits, y = c.index.val, colour = factor(m.try), group = m.try)) +
  geom_line() +
  geom_point(data = cv.performance) +
  geom_hline(
    yintercept =
      models.performance$c.index[models.performance$model == 'rf-varselrf'],
    colour = 'grey'
  ) +
  geom_hline(
    yintercept =
      models.performance$c.index[models.performance$model == 'rf-varselmiss'],
    colour = 'grey'
  ) +
  coord_cartesian(
    ylim = c(0.65, 0.8)
  )
