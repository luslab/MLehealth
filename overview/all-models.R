models.include <-
  c(
    'age', 'cox', 'cox disc', 'cox imp', 'cox imp disc', 'rfsrc', 'rfsrc imp',
    'rf-logrank', 'cox-logrank disc', 'cox-elnet disc'
  )

source('../lib/handy.R')
requirePlus('ggplot2', 'cowplot')

models.performance.all <- readTablePlus('../../output/models-performance-manual.tsv')

models.performance.all$x.labels <-
  paste0(
    models.performance.all$model,
    ifelse(models.performance.all$imputation, ' imp', ''),
    ifelse(models.performance.all$discretised, ' disc', '')
  )

# Currently different scripts quote either the area under the curve or a pre
# one-minused the area under the curve...so standardise that
big.calibration.scores <- models.performance.all$calibration.score > 0.5
models.performance.all[big.calibration.scores, c('calibration.score', 'calibration.score.lower', 'calibration.score.upper')] <-
  1 - models.performance.all[big.calibration.scores, c('calibration.score', 'calibration.score.lower', 'calibration.score.upper')]

models.performance <- data.frame()
for(model in models.include) {
  models.performance <-
    rbind(
      models.performance,
      models.performance.all[models.performance.all$x.labels == model, ]
    )
}
# Turn this into a factor with defined levels so ggplot respects the order above
models.performance$x.labels <-
  factor(models.performance$x.labels, levels = models.include)

plot.c.index <-
  ggplot(models.performance, aes(x = x.labels, y = c.index)) +
  geom_bar(stat='identity', aes(fill = model)) +
  geom_errorbar(
    aes(ymin = c.index.lower, ymax = c.index.upper), width = 0.1
  ) +
  coord_cartesian(
    ylim = c(0.75, 0.81)
  ) +
  theme(legend.position = "none")
  
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
  ) +
  theme(legend.position = "none")

plot_grid(
  plot.c.index, plot.calibration,
  labels = c("A", "B"),
  align = "v", ncol = 1
)

ggsave(
  '../../output/all-models-performance.pdf',
  width = 16,
  height = 10,
  units = 'cm',
  useDingbats = FALSE
)
