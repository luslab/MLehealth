source('../lib/handymedical.R', chdir = TRUE)
requirePlus('cowplot')


# Load the variable importances from the Cox model
cox.miss <-
  readRDS('../../output/caliber-replicate-with-missing-survreg-6-linear-age-surv-boot.rds')
cox.miss.vars <- bootStats(cox.miss, uncertainty = '95ci')
cox.miss.vars$var <- rownames(cox.miss.vars)
cox.miss.vars <- subset(cox.miss.vars, startsWith(var, 'vimp.c.index'))
cox.miss.vars$var <-
  substring(cox.miss.vars$var, nchar('vimp.c.index.') + 1)

# Load the variable importances from the random forest
rf.boot <- read.csv('../../output/rfsrc-cv-nsplit-try3-boot-all.csv')
rf.boot.vars <- bootStatsDf(rf.boot)
rf.boot.vars$var <- rownames(rf.boot.vars)
rf.boot.vars <- subset(rf.boot.vars, startsWith(var, 'vimp.c.index'))
rf.boot.vars$var <-
  substring(rf.boot.vars$var, nchar('vimp.c.index.') + 1)

var.imp.compare <- merge(cox.miss.vars, rf.boot.vars, by = c('var'))


# Plot a scatterplot of them
rf.vs.cox <-
  ggplot(
    var.imp.compare,
    aes(
      x = val.x, xmin = lower.x, xmax = upper.x,
      y = val.y, ymin = lower.y, ymax = upper.y
    )
  ) +
    geom_point() +
    geom_errorbar() +
    geom_errorbarh() +
    coord_cartesian(xlim = c(0, 0.03), ylim = c(0, 0.03))

print('Spearman correlation coefficient of variable importances:')
print(cor(var.imp.compare$val.x, var.imp.compare$val.y, method = 'spearman'))

# Load the variable importances from the big data model
cox.bigdata <- read.csv('../../output/cox-bigdata-varsellogrank-01-boot-all.csv')
cox.bigdata.vars <- bootStatsDf(cox.bigdata)
cox.bigdata.vars$var <- rownames(cox.bigdata.vars)
cox.bigdata.vars <- subset(cox.bigdata.vars, startsWith(var, 'vimp.c.index'))
cox.bigdata.vars$var <-
  substring(cox.bigdata.vars$var, nchar('vimp.c.index.') + 1)

cox.bigdata.vars <-
  cox.bigdata.vars[order(cox.bigdata.vars$val, decreasing = TRUE)[1:20], ]

cox.bigdata.vars <- cox.bigdata.vars[order(cox.bigdata.vars$val, decreasing = FALSE), ]

cox.bigdata.vars$description <- lookUpDescriptions(cox.bigdata.vars$var)

cat('c(', paste0("'", as.character(cox.bigdata.vars$description), "',"), ')', sep = '\n')

cox.bigdata.vars$description.manual <-
  factorOrderedLevels(
      c(
        'ALT',
        'PVD',
        'Hb',
        'Dementia',
        'Albumin',
        'Cardiac glycosides',
        'LV failure',
        'Home visit',
        'Oestrogens/HRT',
        'Chest pain',
        'Na',
        'WCC',
        'ALP',
        'Lymphocyte count',
        'Diabetes',
        'BMI ',
        'Weight',
        'Loop diuretics',
        'Smoking status',
        'Age'
      )
  )

# Plot a bar graph of them
cox.bigdata.plot <-
  ggplot(
    cox.bigdata.vars,
    aes(x = description.manual, y = val, ymin = lower, ymax = upper)
    ) +
    geom_bar(stat = 'identity') + 
    geom_errorbar(width = 0.25) +
    coord_flip()

# Random forest big data
rf.bigdata <- read.csv('../../output/rf-bigdata-varsellogrank-02-boot-all.csv')
rf.bigdata.vars <- bootStatsDf(rf.bigdata)
rf.bigdata.vars$var <- rownames(rf.bigdata.vars)
rf.bigdata.vars <- subset(rf.bigdata.vars, startsWith(var, 'vimp.c.index'))
rf.bigdata.vars$var <-
  substring(rf.bigdata.vars$var, nchar('vimp.c.index.') + 1)

rf.bigdata.vars <-
  rf.bigdata.vars[order(rf.bigdata.vars$val, decreasing = TRUE)[1:20], ]

rf.bigdata.vars <- rf.bigdata.vars[order(rf.bigdata.vars$val, decreasing = FALSE), ]

cat('c(', paste0("'", as.character(rf.bigdata.vars$description), "',"), ')', sep = '\n')

rf.bigdata.vars$description <- lookUpDescriptions(rf.bigdata.vars$var)

rf.bigdata.vars$description.manual <-
  factorOrderedLevels(
    c(
      'Fit note',
      'Stimulant laxatives',
      'Urea',
      'Hypertension',
      'Cardiac glycosides',
      'Beta2 agonists',
      'Telephone encounter',
      'Feet examination',
      'Creatinine',
      'Screening',
      'Osmotic laxatives',
      'ACE inhibitors',
      'Beta blockers',
      'Home visit',
      'Analgesics',
      'Blood pressure',
      'Chest pain',
      'Loop diuretics',
      'Smoking status',
      'Age'
    )
  )

# Plot a bar graph of them
rf.bigdata.plot <-
  ggplot(
    rf.bigdata.vars,
    aes(x = description.manual, y = val, ymin = lower, ymax = upper)
  ) +
    geom_bar(stat = 'identity') + 
    geom_errorbar(width = 0.25) +
    coord_flip()



# Combine for output
plot_grid(
  rf.bigdata.plot, cox.bigdata.plot,
  labels = c('A', 'B', 'C'),
  align = "h", ncol = 2
)

ggsave(
  '../../output/variable-importances.pdf',
  width = 16,
  height = 9,
  units = 'cm',
  useDingbats = FALSE
)
