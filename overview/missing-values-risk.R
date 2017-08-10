source('../lib/handymedical.R', chdir = TRUE)

requirePlus('survminer', 'cowplot')

models.base <- '../../output'
cox.missing.filename <- 'caliber-replicate-with-missing-model-survreg-bootstrap-1.rds'
cox.missing.riskdists.filename <- 'caliber-replicate-with-missing-survreg-3-risk-violins.csv'
cox.missing.riskcats.filename <- 'caliber-replicate-with-missing-survreg-3-risk-cats.csv'

data.filename <- '../../data/cohort-sanitised.csv'

survPlots <- function(...) {
  surv.fits <- list(...)
  
  df <- data.frame()
  for(i in 1:length(surv.fits)) {
    df <-
      rbind(
        df,
        data.frame(
          variable = as.character(surv.fits[[i]]$call)[2],
          stratum = '1',
          time = surv.fits[[i]][1]$time,
          surv = surv.fits[[i]][1]$surv,
          lower = surv.fits[[i]][1]$lower,
          upper = surv.fits[[i]][1]$upper
        ),
        data.frame(
          variable = as.character(surv.fits[[i]]$call)[2],
          stratum = '2',
          time = surv.fits[[i]][2]$time,
          surv = surv.fits[[i]][2]$surv,
          lower = surv.fits[[i]][2]$lower,
          upper = surv.fits[[i]][2]$upper
        )
      )
  }
  
  ggplot(
    df,
    aes(x = time, y = surv, ymin = lower, ymax = upper)
  ) +
    geom_line(aes(colour = stratum)) +
    geom_ribbon(aes(fill = stratum), alpha = 0.4) +
    facet_grid(variable ~ .)
}

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
    ) +
  scale_y_continuous(breaks = c(0.75, 1.0, 1.25, 1.5))


# Kaplan-Meier survival curves for a few example variables being missing
COHORT <- fread(data.filename)
COHORT <- subset(COHORT, !exclude & time_death > 0)

# Calculate the curves
km.hdl <-
 survfit(
    Surv(time_death, endpoint_death == 'Death') ~ is.na(hdl_6mo),
    data = COHORT
 )
km.total_chol <-
  survfit(
    Surv(time_death, endpoint_death == 'Death') ~ is.na(total_chol_6mo),
    data = COHORT
  )
km.crea <-
  survfit(
    Surv(time_death, endpoint_death == 'Death') ~ is.na(crea_6mo),
    data = COHORT
  )

km.missingness <-
  survPlots(km.hdl, km.total_chol, km.crea) +
  theme(
    legend.position = "none",
    # Remove grey labels on facets
    strip.background = element_blank(),
    strip.text = element_blank()
  )

theme_set(theme_cowplot(font_size = 10))
plot_grid(
  disc.vs.cont.risks.plot, risk.violins.plot,
  km.missingness,
  labels = c("A", "B", "C"),
  align = "h", nrow = 1
)
ggsave(
  '../../output/missing-values-risk.pdf',
  width = 16,
  height = 5,
  units = 'cm',
  useDingbats = FALSE
)
