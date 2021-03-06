source('../lib/handymedical.R', chdir = TRUE)

bootstrap.base <- '../../output'

bootstrap.files <-
  c(
    cox.miss = 'caliber-replicate-with-missing-survreg-6-linear-age-surv-boot.rds',
    cox.disc = 'all-cv-survreg-boot-try5-surv-model.rds',
    cox.imp = 'caliber-replicate-imputed-survreg-4-surv-boot-imp.rds',
    rf = 'rfsrc-cv-nsplit-try3-boot-all.csv',
    rf.imp = 'rf-imputed-try1-boot.rds',
    rfbig = 'rf-bigdata-varsellogrank-02-boot-all.csv',
    coxbig = 'cox-bigdata-varsellogrank-01-boot-all.csv'
  )

# Helper functions

# Turn a boot object into a data frame
bootstrap2Df <- function(x) {
  df <- data.frame(x$t)
  names(df) <- names(x$t0)
  df
}

# Make sure calibration scores are bigger = better
calibrationFix <- function(x) {
  if(mean(x) < 0.5) {
    x <- 1 - x
  }
  x
}

n <- length(bootstrap.files)

bootstraps <- list()

for(i in 1:n) {
  if(fileExt(bootstrap.files[i]) == 'rds'){
    bootstraps[[i]] <- readRDS(file.path(bootstrap.base, bootstrap.files[i]))
    
    if(class(bootstraps[[i]]) == 'list') {
      # If it's a list, then it's from an imputed dataset with separate bootstraps
      # Turn each of these into a data frame and then combine them together.
      # (data.frame is needed because rbindlist returns a data.table)
      bootstraps[[i]] <-
        data.frame(rbindlist(lapply(bootstraps[[i]], bootstrap2Df)))
    } else {
      bootstraps[[i]] <- bootstrap2Df(bootstraps[[i]] )
    }
  } else{
    bootstraps[[i]] <- read.csv(file.path(bootstrap.base, bootstrap.files[i]))
  }
}

x1x2 <- combn(1:n, 2)
x1 <- x1x2[1,]
x2 <- x1x2[2,]


bootstrap.differences <- data.frame()
for(i in 1:length(x1)) {
  # C-index
  col.1.c.index <-
    which(names(bootstraps[[x1[i]]]) %in% c('c.test', 'c.index'))
  col.2.c.index <-
    which(names(bootstraps[[x2[i]]]) %in% c('c.test', 'c.index'))
  boot.diff <-
    bootstrapDiff(
      bootstraps[[x1[i]]][, col.1.c.index],
      bootstraps[[x2[i]]][, col.2.c.index]
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
  col.1.calib <-
    which(names(bootstraps[[x1[i]]]) == 'calibration.score')
  col.2.calib <-
    which(names(bootstraps[[x2[i]]]) == 'calibration.score')
  boot.diff <-
    bootstrapDiff(
      calibrationFix(bootstraps[[x1[i]]][, col.1.calib]),
      calibrationFix(bootstraps[[x2[i]]][, col.2.calib])
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

print(cbind(bootstrap.differences[, c('model.1', 'model.2', 'var')], round(bootstrap.differences[, 4:6], 3)))

write.csv(bootstrap.differences, '../../output/bootstrap-differences.csv')