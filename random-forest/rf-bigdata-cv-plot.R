require(ggplot2, data.table)

calibration.data <- fread('../../output/rf-bigdata-try4-calibration.csv')

calibration.data <-
  rbind(
    calibration.data,
    fread('../../output/rf-bigdata-try5-calibration.csv')
  )

# Aggregate by mtry and nsplit
calibration.data.avg <-
  calibration.data[, mean(c.index.val), by = c('m.try', 'n.splits')]

names(calibration.data.avg)[3] <- 'c.index'

calibration.data.avg$m.try <- factor(calibration.data.avg$m.try)
calibration.data.avg$n.splits <- factor(calibration.data.avg$n.splits)

ggplot(calibration.data.avg, aes(x = n.splits, y = m.try, fill = c.index)) +
  geom_tile()
