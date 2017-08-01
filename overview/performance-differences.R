source('../lib/handymedical.R', chdir = TRUE)

bootstrap.base <- '../../output'

bootstrap.files <-
  c(
    cox.miss = 'caliber-replicate-with-missing-survreg-2-model-bootstrap.rds',
    cox.disc = 'all-cv-survreg-boot-try5-surv-model.rds'
  )

n <- length(bootstrap.files)

bootstraps <- list()

for(i in 1:n) {
  bootstraps[[i]] <- readRDS(file.path(bootstrap.base, bootstrap.files[i]))
}

x1x2 <- combn(1:n, 2)
x1 <- x1x2[1,]
x2 <- x1x2[2,]


bootstrap.differences <- data.frame()
for(i in 1:length(x1)) {
  # C-index
  row.1.c.index <-
    which(names(bootstraps[[x1[i]]]$t0) %in% c('c.test', 'c.index'))
  row.2.c.index <-
    which(names(bootstraps[[x2[i]]]$t0) %in% c('c.test', 'c.index'))
  boot.diff <-
    bootstrapDiff(
      bootstraps[[x1[i]]]$t[, row.1.c.index],
      bootstraps[[x2[i]]]$t[, row.2.c.index]
    )
  
  bootstrap.differences <-
    rbind(
      bootstrap.differences,
      data.frame(
        model.1 = names(bootstrap.files)[x1[i]],
        model.2 = names(bootstrap.files)[x2[i]],
        var = 'c.index',
        diff = boot.diff['val'],
        lower = boot.diff['lower'],
        upper = boot.diff['upper']
      )
    )
  
  # Calibration score
  row.1.calib <- which(names(bootstraps[[x1[i]]]$t0) == 'calibration.score')
  row.2.calib <- which(names(bootstraps[[x2[i]]]$t0) == 'calibration.score')
  boot.diff <-
    bootstrapDiff(
      bootstraps[[x1[i]]]$t[, row.1.calib],
      bootstraps[[x2[i]]]$t[, row.1.calib]
    )
  
  bootstrap.differences <-
    rbind(
      bootstrap.differences,
      data.frame(
        model.1 = names(bootstrap.files)[x1[i]],
        model.2 = names(bootstrap.files)[x2[i]],
        var = 'calibration.score',
        diff = boot.diff['val'],
        lower = boot.diff['lower'],
        upper = boot.diff['upper']
      )
    )
}

# Remove nonsense row names
rownames(bootstrap.differences) <- NULL