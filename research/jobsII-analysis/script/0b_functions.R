



..code_table1 <- function(dat, X.names, X.labels) {

    dat$treat     <- factor(dat$treat,     labels = c("Control",
                                                      "Treatment"))
    dat$edu       <- factor(dat$edu,       labels = c('less than highschool',
                                                      'highschool',
                                                      'some college',
                                                      'college degree',
                                                      'graduate work'))
    dat$marital   <- factor(dat$marital,   labels = c('never married',
                                                      'married',
                                                      'divorced/separated/widowed'))
    dat$hh.income <- factor(dat$hh.income, labels = c('LT15K',
                                                      '15-24K',
                                                      '25-39K',
                                                      '40-49K',
                                                      '50K+'))

    if (!missing(X.names))
        dat <- dat %>% dplyr::rename_with(~ X.labels, all_of(X.names))

    dat
}
