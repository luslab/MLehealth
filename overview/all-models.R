source('../lib/handy.R')
requirePlus('ggplot2', 'cowplot')

models.performance <- readTablePlus('../../output/models-performance.tsv')

models.performance$x.labels <-
  paste(
    models.performance$model,
    ifelse(models.performance$imputation, 'imp', ''),
    ifelse(models.performance$discretised, 'disc', '')
  )


plot.c.index <-
  ggplot(models.performance, aes(x = x.labels, y = c.index)) +
  geom_bar(stat='identity', aes(fill = model)) +
  geom_errorbar(
    aes(ymin = c.index.lower, ymax = c.index.upper), width = 0.1
  ) +
  coord_cartesian(
    ylim = c(0.7, 0.8)
  )
  
plot.calibration <-
  ggplot(models.performance, aes(x = x.labels, y = 1 - calibration.score)) +
  geom_bar(stat='identity', aes(fill = model)) +
  geom_errorbar(
    aes(
      ymin = 1 - calibration.score.lower, ymax = 1 - calibration.score.upper
    ),
    width = 0.1
  ) +
  coord_cartesian(
    ylim = c(0.8, 1.0)
  )

plot_grid(
  plot.c.index, plot.calibration,
  labels = c("A", "B"),
  align = "v", ncol = 1
)
