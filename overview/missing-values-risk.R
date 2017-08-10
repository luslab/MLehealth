source('../lib/handymedical.R', chdir = TRUE)

requirePlus('survminer', 'cowplot')

models.base <- '../../output'
cox.missing.filename <- 'caliber-replicate-with-missing-model-survreg-bootstrap-1.rds'
cox.missing.riskdists.filename <- 'caliber-replicate-with-missing-survreg-3-risk-violins.csv'
cox.missing.riskcats.filename <- 'caliber-replicate-with-missing-survreg-3-risk-cats.csv'

data.filename <- '../../data/cohort-sanitised.csv'




cox.missing.boot <- readRDS(file.path(models.base, cox.missing.filename))

fit.risks <- bootStats(cox.missing.boot, '95ci')

disc.vs.cont.risks.plot <-
  ggplot(
    fit.risks,
    aes(
      x = val, xmin = lower, xmax = upper,
      y = val, ymin = lower, ymax = upper
    )
  ) +
  geom_point() +
  geom_errorbar() +
  geom_errorbarh()



risk.dist.by.var <-
  read.csv(file.path(models.base, cox.missing.riskdists.filename))
risk.cats <- read.csv(file.path(models.base, cox.missing.riskcats.filename))


# Plot the results
risk.violins.plot <- 
  ggplot() +
    # First, and therefore at the bottom, draw the reference line at risk = 1
    geom_hline(yintercept = 1) +
    # Then, on top of that, draw the violin plot of the risk from the data
    geom_violin(data = risk.dist.by.var, aes(x = quantity, y = risk)) +
    geom_pointrange(
      data = risk.cats,
      aes(x = quantity, y = our_value, ymin = our_lower,
          ymax = our_upper),
      
      position = position_jitter(width = 0.1)
    ) +
    geom_text(
      data = risk.cats,
      aes(
        x = quantity,
        y = our_value,
        label = quantity.level
      )
    )


# Kaplan-Meier survival curves for a few example variables being missing
COHORT <- fread(data.filename)
COHORT <- subset(COHORT, !exclude & time_death > 0)

km.hdl <-
  survfit(
    Surv(time_death, endpoint_death == 'Death') ~ is.na(hdl_6mo),
    data = COHORT
  )

km.hdl.plot <-
  ggsurvplot(
    km.hdl, data = COHORT, conf.int = TRUE, censor = FALSE, legend = 'none'
  )

km.total_chol <-
  survfit(
    Surv(time_death, endpoint_death == 'Death') ~ is.na(total_chol_6mo),
    data = COHORT
  )

km.total_chol.plot <-
  ggsurvplot(km.total_chol, data = COHORT, conf.int = TRUE, censor = FALSE, legend = 'none')

km.crea <-
  survfit(
    Surv(time_death, endpoint_death == 'Death') ~ is.na(crea_6mo),
    data = COHORT
  )

km.crea.plot <-
  ggsurvplot(km.crea_6mo, data = COHORT, conf.int = TRUE, censor = FALSE, legend = 'none')




plot_grid(
  disc.vs.cont.risks.plot, risk.violins.plot,
  #labels = c("A", "B"),
  align = "v", nrow = 1
)


arrange_ggsurvplots(
  list(km.hdl.plot, km.total_chol.plot, km.crea.plot),
  ncol = 1, nrow = 3
  )
