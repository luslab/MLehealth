################################################################################
### USER VARIABLES #############################################################
################################################################################

setwd('R:/Pop_Health/Farr_Luscombe/')


### Relating to constructed cohort #############################################
# the R data file containing the patient cohort
cohort.file <- '2_Cohort/COHORT_SCAD_makecohort2.rda'

# columns to remove
kill.cols <- c('anonpatid', 'pracid','year_of_birth')

# columns to make relative to indexdate
date.cols <- c('afterentry_acs','afterentry_mi_nos',
			'afterentry_nstemi','afterentry_stemi','afterentry_ua',
			'crd','date_entry','date_exit','date_mi_endpoint','deathdate','dob','dod',
			'dod_combined','earliest_chd','earliest_hf','earliest_mi','earliest_sa',
			'earliest_ua','frd','hes_end','hes_start','indexdate','praclcd',
			'pracuts','tod','recent_acs','recent_mi_nos','recent_nstemi',
			'recent_stemi','recent_ua','smokdate','endpoint_coronary_date',
			'endpoint_death_date'
)

### Relating to other files ####################################################
data.path <- '1a_ExtractedData'

hes.file <- file.path(data.path, 'hes_diag_epi.csv')
n.top.icd <- 100

procedures.file <- file.path(data.path, 'hes_procedures.csv')
n.top.procedures <- 100

tests.files <-
  file.path(
    data.path,
    paste0('test.part.', 0:3)
  )
# Number of tests to extract
n.top.tests <- 100
# Find tests as close as possible to the first value here, and no more than the
# second timepoint before it
test.timepoints <- c(0, 183) #days

clinical.files <-
  file.path(
    data.path,
    paste0('clinical.part.', 0:3)
  )
n.top.clinical.values <- 100
n.top.clinical.history <- 100

therapy.files <-
  file.path(
    data.path,
    paste0('therapy.part.', 0:7)
  )
n.top.therapy <- 100


################################################################################
### END USER VARIABLES #########################################################
################################################################################


### Load the patient cohort to act as a base ###################################

# load the standard COHORT variable, which is a data table called COHORT
load(cohort.file)

### Do a bunch of silly fixes on the data ######################################
# cabg_6mo is 1s and 0s, which should be TRUE and FALSE
COHORT$cabg_6mo <- as.logical(COHORT$cabg_6mo)
# long_nitrate is 1s and NAs, which should be TRUE and FALSE...
COHORT$long_nitrate <- as.logical(COHORT$long_nitrate)
COHORT$long_nitrate[is.na(COHORT$long_nitrate)] <- FALSE
# pci_6mo is 1s and 0s, which should be TRUE and FALSE
COHORT$pci_6mo <- as.logical(COHORT$pci_6mo)

### Append other data types ####################################################

# We'll be needing handymedical from here on in
source('Andrew/lib/handymedical.R', chdir = TRUE)
require(CALIBERlookups)

percentMissing <- function(x, sf = 3) {
  round(sum(is.na(x))/length(x), digits = sf)*100
}

### Hospital episodes statistics ###############################################

# read in the hospital data (HES)
hes.diag.epi <- readMedicalData(
  hes.file,
  c("anonpatid", "epistart", "epiend", "icd"),
  c("integer", "date", "date", "factor")
)

# remove non alphanumerics (trailing -s on some ICD codes)
hes.diag.epi$icd <- gsub('[^[:alnum:]]', '', hes.diag.epi$icd)

# Now, merge with the indexdate...we only want data before then, because looking
# after it is cheating, and using all data rather than just stuff before may
# distort our choice of variables as some variables may be very common after
# entering the study, but less so before... (This actually isn't much of an
# issue...only 9 variables differ. Still!)
hes.diag.epi <-
  merge(
    hes.diag.epi, COHORT[, c('anonpatid', 'indexdate')],
    by = 'anonpatid', all = TRUE
  )
hes.diag.epi$relativedate <- hes.diag.epi$indexdate - hes.diag.epi$epistart

# Now, remove all the negative relative dates, because they're in the past
hes.diag.epi <- hes.diag.epi[hes.diag.epi$relativedate >= 0, ]

# get aggregate statistics for each ICD code
hes.by.icd <- hes.diag.epi[, length(unique(anonpatid)), by = icd]
names(hes.by.icd) <- c('icd', 'n.pat')

# Take the top n by number of patients with that code
top.icd <-
  hes.by.icd$icd[order(hes.by.icd$n.pat, decreasing = TRUE)[1:n.top.icd]]

# And now, let's put how far in the past the patient was first diagnosed with
# each of these things into the table...

# First, discard all the rows corresponding to non-top ICDs
hes.diag.epi.top <- hes.diag.epi[hes.diag.epi$icd %in% top.icd, ]
# Then, keep only the earliest instance of each per patient, because those are
# the ones we care about
hes.diag.epi.top.earliest <- 
  hes.diag.epi.top[, min(relativedate), by = c('anonpatid', 'icd')]

# Per ICD code, add a new column to the cohort and put in numbers
hes.diag.epi.top.earliest <-
  dcast(
    data = hes.diag.epi.top.earliest,
    formula = anonpatid ~ icd,
    value.var = "V1"
  )

# Add a prefix to the names to keep track
names(hes.diag.epi.top.earliest)[2:(n.top.icd + 1)] <-
  paste0('hes.icd.', names(hes.diag.epi.top.earliest)[2:(n.top.icd + 1)])

# Keep the names for when we need to make them relative to the indexdate when
# anonymising later
names.hes.diag.icd <- names(hes.diag.epi.top.earliest)[2:(n.top.icd + 1)]

# Now, merge with the original cohort
COHORT <-
  merge(
    COHORT,
    hes.diag.epi.top.earliest,
    by = c('anonpatid'),
    all = TRUE
  )

hes.icd.summary <-
  data.frame(
    percent.missing = 
      sort(
        sapply(
          COHORT[, names(COHORT)[startsWith(names(COHORT), 'hes.icd.')], with = FALSE],
          percentMissing
        )
      )
  )

hes.icd.summary$code <- substring(rownames(hes.icd.summary), 9)

hes.icd.summary <- merge(hes.icd.summary, CALIBER_DICT[, c('code', 'term')], by = 'code')

print(hes.icd.summary[order(hes.icd.summary$percent.missing),])

### Hospital procedures data ###################################################

# read in the hospital data (HES)
hes.procedures <- readMedicalData(
  procedures.file,
  c("anonpatid", "opcs", "evdate"),
  c("integer", "factor", "date")
)

# Now, merge with the indexdate...we only want data before then, because looking
# after it is cheating, and using all data rather than just stuff before may
# distort our choice of variables as some variables may be very common after
# entering the study, but less so before... (14 differ...)
hes.procedures <-
  merge(
    hes.procedures, COHORT[, c('anonpatid', 'indexdate')],
    by = 'anonpatid', all = TRUE
  )
hes.procedures$relativedate <- hes.procedures$indexdate - hes.procedures$evdate

# Now, remove all the negative relative dates, because they're in the past
hes.procedures <- hes.procedures[relativedate >= 0]

# get aggregate statistics for each OPCS code
hes.by.opcs <- hes.procedures[, length(unique(anonpatid)), by = opcs]
names(hes.by.opcs) <- c('opcs', 'n.pat')

# Take the top n by number of patients with that code
top.opcs <-
  hes.by.opcs$opcs[order(hes.by.opcs$n.pat, decreasing = TRUE)[1:n.top.procedures]]

# And now, let's put how far in the past the patient was first diagnosed with
# each of these things into the table...

# First, discard all the rows corresponding to non-top OPCS codes
hes.procedures.top <- hes.procedures[hes.procedures$opcs %in% top.opcs, ]
# Then, keep only the earliest instance of each per patient, because those are
# the ones we care about
hes.procedures.top.earliest <-
  hes.procedures.top[, min(relativedate), by = c('anonpatid', 'opcs')]

# Per OPCS code, add a new column to the cohort and put in numbers
hes.procedures.top.earliest <-
  dcast(
    data = hes.procedures.top.earliest,
    formula = anonpatid ~ opcs,
    value.var = "V1"
  )

# Add a prefix to the names to keep track
names(hes.procedures.top.earliest)[2:(n.top.procedures + 1)] <-
  paste0('hes.opcs.', names(hes.procedures.top.earliest)[2:(n.top.procedures + 1)])

# Keep the names for when we need to make them relative to the indexdate when
# anonymising later
names.procedures.opcs <- names(hes.procedures.top.earliest)[2:(n.top.icd + 1)]

# Now, merge with the original cohort
COHORT <-
  merge(
    COHORT,
    hes.procedures.top.earliest,
    by = c('anonpatid'),
    all = TRUE
  )

# Summarise by percent missing
hes.opcs.summary <-
  data.frame(
    percent.missing = 
      sort(
        sapply(
          COHORT[, names(COHORT)[startsWith(names(COHORT), 'hes.opcs.')], with = FALSE],
          percentMissing
        )
      )
  )

hes.opcs.summary$code <- substring(rownames(hes.opcs.summary), 10)

hes.opcs.summary <- merge(hes.opcs.summary, CALIBER_DICT[, c('code', 'term')], by = 'code')

print(hes.opcs.summary[order(hes.opcs.summary$percent.missing),])

### Test results ###############################################################

# read in the test data
tests.data <- readMedicalData(
  tests.files,
  # data1 is operator (eg =, > etc), data2 is value, data 3 is unit of measure
  c("anonpatid", "eventdate", "enttype", "data1", "data2", "data3"),
  c("integer", "date", "integer", "integer", "numeric", "integer")
)

# First, discard those where operator is not =, because > and < etc will
# introduce complexity, and drop data1 since it's now useless
tests.data <- tests.data[data1 == 3]
tests.data$data1 <- NULL

# Now, let's subtract the indexdate from every test so we can choose the ones
# closest to the desired dates...
tests.data <-
  merge(
    tests.data, COHORT[, c('anonpatid', 'indexdate')],
    by = 'anonpatid', all = TRUE
  )
tests.data$relativedate <- tests.data$indexdate - tests.data$eventdate

# Only keep positive values; negative ones are in the future which is cheating!
tests.data <- tests.data[relativedate >= 0]
# Only 89001 have any test results from before the indexdate!

# Discard any test results from times greater than the longest time ago to check
tests.data <- tests.data[relativedate < test.timepoints[2]]
# And this drops to 57972 in the preceding six months

# Per patient and per test, keep the smallest relativedate value only
tests.data <-
  tests.data[, .SD[which.min(relativedate)], by = list(anonpatid, enttype)]

# aggregate by test type (enttype, which covers a range of Read codes which can
# mean the same test) and unit of measure (so we can choose the tests with the
# units with the best coverage)
# We do this after all the preprocessing (ie only looking at the most recent
# value per patient before but not too far before the indexdate) because
# otherwise a lot of the top tests change. Presumably 
tests.by.test <-
  tests.data[, length(unique(anonpatid)), by = c('enttype', 'data3')]
names(tests.by.test) <- c('enttype', 'data3', 'n.pat')

top.tests <- tests.by.test$enttype[order(tests.by.test$n.pat, decreasing = TRUE)[1:n.top.tests]]
top.units <- tests.by.test$data3[order(tests.by.test$n.pat, decreasing = TRUE)[1:n.top.tests]]

# Now, discard all the rows corresponding to non-top tests
tests.data <-
  tests.data[
    # Needs to be exact match of enttype and variable combination
    paste0(enttype, '!!!', data3) %in% paste0(top.tests, '!!!', top.units)
  ]

# Make a column per test
tests.wide <-
  dcast(
    data = tests.data,
    formula = anonpatid ~ enttype + data3,
    value.var = "data2"
  )

# Add a prefix to the names to keep track
# (There may not be n.top.tests columns here as some get lost during the paring
# down processes above...)
names(tests.wide)[-1] <-
  paste0('tests.enttype.data3.', names(tests.wide)[-1])

# Now, merge with the original cohort
COHORT <-
  merge(
    COHORT,
    tests.wide,
    by = c('anonpatid'),
    all = TRUE
  )


### GP diagnosis data ##########################################################

# read in the GP data (clinical)
clinical.data <- readMedicalData(
  clinical.files,
  c("anonpatid", "eventdate", "medcode", "enttype", "data1", "data2", "data3", "data4", "data5", "data6", "data7"),
  c("integer", "date", "integer", "integer", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
)

# Now, merge with the indexdate...we only want data before then, because looking
# after it is cheating, and using all data rather than just stuff before may
# distort our choice of variables as some variables may be very common after
# entering the study, but less so before... (14 differ...)
clinical.data <-
  merge(
    clinical.data, COHORT[, c('anonpatid', 'indexdate')],
    by = 'anonpatid', all = TRUE
  )
clinical.data$relativedate <- clinical.data$indexdate - clinical.data$eventdate

# Now, remove all the negative relative dates, because they're in the past
clinical.data <- clinical.data[relativedate >= 0]

# Find out which enttypes have associated data values, to distinguish purely
# binary variables (like medical history, family history etc) from test results
# and so on which have data1, data2, etc values.
# From a quick scan, anything with any data at all has a data1 value, so we can
# use that as a proxy.
clinical.by.data1 <- clinical.data[, sum(!is.na(data1)), by = enttype]

# Split into two data tables, the ones with associated data, and without
clinical.history <-
  clinical.data[
    enttype %in% clinical.by.data1$enttype[clinical.by.data1$V1 == 0]
  ]
clinical.values <-
  clinical.data[
    enttype %in% clinical.by.data1$enttype[clinical.by.data1$V1 > 0]
    ]

### Clinical history

# get aggregate statistics for each medcode
clinical.history.by.medcode <-
  clinical.history[, length(unique(anonpatid)), by = medcode]
names(clinical.history.by.medcode) <- c('medcode', 'n.pat')

# Print a table of the leading medcodes we've found
print(
  merge(
    clinical.history.by.medcode[order(n.pat, decreasing = TRUE)[1:100],],
    CALIBER_DICT[, c('medcode', 'term')], by = 'medcode'
  )
)

# Take the top n by number of patients with that code
top.history <-
  clinical.history.by.medcode$medcode[
    order(
      clinical.history.by.medcode$n.pat, decreasing = TRUE
    )[1:n.top.clinical.history]
  ]

# And now, let's put how far in the past the patient was first diagnosed with
# each of these things into the table...

# First, discard all the rows corresponding to non-top OPCS codes
clinical.history <- clinical.history[clinical.history$medcode %in% top.history, ]
# Then, keep only the earliest instance of each per patient, because those are
# the ones we care about
clinical.history <-
  clinical.history[, min(relativedate), by = c('anonpatid', 'medcode')]

# Per medcode, add a new column to the cohort and put in numbers
clinical.history <-
  dcast(
    data = clinical.history,
    formula = anonpatid ~ medcode,
    value.var = "V1"
  )

# Add a prefix to the names to keep track
names(clinical.history)[2:(n.top.clinical.history + 1)] <-
  paste0('clinical.history.', names(clinical.history)[2:(n.top.clinical.history + 1)])

# Now, merge with the original cohort
COHORT <-
  merge(
    COHORT,
    clinical.history,
    by = c('anonpatid'),
    all = TRUE
  )

### Clinical values

# Next, let's melt the values data by data1, data2 etc. Each potentially
# contains a separate measurement (eg for entity type 13, which is weight, data1
# is the weight, data2 is the weight centile [always blank in this dataset!] and
# data3 is BMI)
clinical.values <-
  melt(
    clinical.values,
    id.vars = c("anonpatid", "relativedate", "medcode", "enttype"),
    measure.vars = paste0("data", 1:7)
  )
# Because there are lots of data types, the value column gets coerced to double.
# This is fine for now, because it's all numeric, and whilst some values are
# categorical represented as integers, random forests don't care about that.

# Remove all NA values
clinical.values <- clinical.values[!is.na(value)]

# Remove all values measured too long ago
clinical.values <- clinical.values[relativedate <= 183]

# aggregate by enttype and which data column the test came from
clinical.values.by.type <-
  clinical.values[, length(unique(anonpatid)), by = c('enttype', 'variable')]
names(clinical.values.by.type) <- c('enttype', 'variable', 'n.pat')

top.clinical.enttypes <-
  clinical.values.by.type$enttype[order(clinical.values.by.type$n.pat, decreasing = TRUE)[1:n.top.clinical.values]]
top.clinical.dataN <-
  clinical.values.by.type$variable[order(clinical.values.by.type$n.pat, decreasing = TRUE)[1:n.top.clinical.values]]

# Now, discard all the rows corresponding to non-top tests
clinical.values <-
  clinical.values[
    # Needs to be exact match of enttype and variable combination
    paste0(enttype, '!!!', variable) %in%
      paste0(top.clinical.enttypes, '!!!', top.clinical.dataN),
    ]

# Per patient and per test, keep the smallest relativedate value only
clinical.values.most.recent <-
  clinical.values[, .SD[which.min(relativedate)], by = c('anonpatid', 'enttype', 'variable')]

# Make a column per test
clinical.values.most.recent.wide <-
  dcast(
    data = clinical.values.most.recent,
    formula = anonpatid ~ enttype + variable,
    value.var = "value"
  )

# Add a prefix to the names to keep track
# (There may not be n.top.tests columns here as some get lost during the paring
# down processes above...)
names(clinical.values.most.recent.wide)[-1] <-
  paste0('clinical.values.', names(clinical.values.most.recent.wide)[-1])

# Now, merge with the original cohort
COHORT <-
  merge(
    COHORT,
    clinical.values.most.recent.wide,
    by = c('anonpatid'),
    all = TRUE
  )


### Therapy ####################################################################

# read in the therapy data
therapy.data <- readMedicalData(
  therapy.files,
  c("anonpatid", "eventdate", "bnfcode"),
  c("integer", "date", "integer")
)
# The other option than bnfcode is prodcode, which refers to specific products
# rather than BNF categories. There are far more of these so, assuming the BNF
# classification is somewhat rational, I'm going to go with that first to
# reduce data sparsity.

# Now, merge with the indexdate...we only want data before then, because looking
# after it is cheating, and using all data rather than just stuff before may
# distort our choice of variables as some variables may be very common after
# entering the study, but less so before... (14 differ...)
therapy.data <-
  merge(
    therapy.data, COHORT[, c('anonpatid', 'indexdate')],
    by = 'anonpatid', all.x = TRUE
  )

therapy.data$relativedate <- therapy.data$indexdate - therapy.data$eventdate

# Now, remove all the negative relative dates, because they're in the past
therapy.data <- therapy.data[relativedate >= 0]
# And remove all data which is too far into the past
therapy.data <- therapy.data[relativedate < 366]

# get aggregate statistics for each BNF code
therapy.by.bnf <- therapy.data[, length(unique(anonpatid)), by = bnfcode]
names(therapy.by.bnf) <- c('bnfcode', 'n.pat')
# Take the top n by number of patients with that code
top.bnf <-
  therapy.by.bnf$bnfcode[order(therapy.by.bnf$n.pat, decreasing = TRUE)[1:n.top.therapy]]

# Discard all the rows corresponding to non-top BNF codes
therapy.data <- therapy.data[therapy.data$bnfcode %in% top.bnf, ]

# Aggregate by number of prescriptions per patient
therapy.data <- therapy.data[, .N, by = list(anonpatid, bnfcode)]

# Per BNF code, add a new column to the cohort and put in numbers
therapy.wide <-
  dcast(
    data = therapy.data,
    formula = anonpatid ~ bnfcode,
    value.var = "N"
  )

# Add a prefix to the names to keep track
names(therapy.wide)[2:(n.top.therapy + 1)] <-
  paste0('bnf.', names(therapy.wide)[2:(n.top.therapy + 1)])

# And merge into the overall cohort
COHORT <-
  merge(
    COHORT,
    therapy.wide,
    by = c('anonpatid'),
    all = TRUE
  )

### Anonymisation steps ########################################################

# delete the columns with obvious privacy issues
COHORT[, (kill.cols) := NULL]

# make the remaining columns relative to the indexdate
indexdate <- as.Date(COHORT$indexdate)

# date.cols specified
for(date.col in c(date.cols)) {
	COHORT[[date.col]] <-
	  # Do indexdate minus, so this is days in the past
	  as.numeric(indexdate - as.Date(COHORT[[date.col]]))
	# ...and therefore negative values are in the future, hence cheating
	COHORT[[date.col]][COHORT[[date.col]] < 0] <- NA
}
# make age an integer
COHORT$age <- round(COHORT$age)

# make pracregion and ethnicity into categories
lookup_pracregion <-
  sample(unique(COHORT$pracregion),length(unique(COHORT$pracregion)))
lookup_ethnicity <-
  sample(unique(COHORT$hes_ethnicity),length(unique(COHORT$hes_ethnicity)))
write.csv(lookup_pracregion, 'lookup_pracregion.csv')
write.csv(lookup_ethnicity, 'lookup_ethnicity.csv')

COHORT$pracregion <-
  as.integer(factor(COHORT$pracregion, levels = lookup_pracregion))
COHORT$hes_ethnicity <-
  as.integer(factor(COHORT$hes_ethnicity, levels = lookup_ethnicity))

# make IMD score into deciles-ish
COHORT$imd_score <- round(COHORT$imd_score/10)

write.csv(COHORT, 'cohort-datadriven.csv')

fake.df <- data.frame(id = 1:10)
for (colname in names(COHORT)) {
	fake.df[,colname] <- sample(COHORT[!is.na(COHORT[, colname]), colname], 10)
}

write.csv(fake.df, 'cohort-sample.csv')
