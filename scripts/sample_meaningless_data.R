source('../lib/handymedical.R', chdir = TRUE)

source.data <- fread('../data/cohort-sample.csv')

n <- 3

meaningless.data <-
  data.frame(
      lapply(
      source.data,
      function(x, n){
        sample(x, n, replace = TRUE)
      }
    )
  )

write.csv(meaningless.data, '../data/cohort-meaningless.csv')
