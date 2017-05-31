for(i in 1:length(names(COHORT.full))) {
  if(names(COHORT.full)[i] %in% surv.predict) {
    guide.type = 'n'
  } else {
    # If it's not in the list, exclude it
    guide.type = 'x'
  }
  cat(i, names(COHORT.full)[i], guide.type, '\n')
}