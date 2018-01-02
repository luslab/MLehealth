source('../lib/handymedical.R', chdir = TRUE)
requirePlus('ggplot2', 'cowplot')

cox.calibration <-
  read.csv('../../output/cox-bigdata-varsellogrank-01-calibration-table.csv')
rf.calibration <-
  read.csv('../../output/rfsrc-cv-nsplit-try3-calibration-table.csv')

cox.calibration.plot <- calibrationPlot(cox.calibration, max.points = 2000)
rf.calibration.plot <- calibrationPlot(rf.calibration, max.points = 2000)

plot_grid(
  cox.calibration.plot, rf.calibration.plot,
  labels = c("C", ""),
  align = "v", ncol = 2
)

ggsave(
  '../../output/models-calibration-eg.pdf',
  width = 16,
  height = 8,
  units = 'cm',
  useDingbats = FALSE
)