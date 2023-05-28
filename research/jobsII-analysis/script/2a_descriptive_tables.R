


library(tidyverse)
library(PIsens)


jobs.dir   <- here::here('research', 'jobsII-analysis')
script.dir <- file.path(jobs.dir, 'script')
tables.dir <- file.path(jobs.dir, 'outputs', 'tables')


jobs <- readRDS(file.path(jobs.dir, "data", "jobsII.rds"))
names(jobs)

X.names <- c("age", "sex", "race", "edu", "marital", "hh.kids",
             "hh.income", "econ.hard", "occu", "wks.unemp",
             "part.motiv", "seek.motiv", "seek.effi",
             "assertive", "depress1")

X.labels <- c('Age', 'Sex', 'Race', 'Education',
              'Marital status', 'Cohabiting children', 'Household income',
              'Economic hardship',
              'Occupation (last steady job)', 'Weeks unemployed',
              'Motivation to participate',
              'Job-seeking motivation', 'Job-seeking self-efficacy',
              'Assertiveness', 'Depressive symptoms')




########################################
#### a simple table 1

source(file.path(script.dir, "0a_load_function.R"))
..loadFunction(file = file.path(script.dir, "0b_functions.R"),
               function.name = "..code_table1")

tab1.unw <- tableone::CreateTableOne(vars = X.labels,
                                     strata = "treat",
                                     data = ..code_table1(jobs, X.names, X.labels),
                                     addOverall = TRUE,
                                     test = FALSE)

..unloadFunction(function.name = "..code_table1")

tab1.unw

# for latex
tab1.unw.matrix <- print(tab1.unw, printToggle = FALSE, noSpaces = TRUE)
tab1.unw.latex  <- xtable::xtable(tab1.unw.matrix)


########################################
#### an inverse propensity score weighted table 1

X.cont   <- c("age", "econ.hard", "wks.unemp", "part.motiv", "seek.motiv",
              "seek.effi", "assertive", "depress1")
X.square <- paste0(paste0("I("     , X.cont), "^2)")
X.sqroot <- paste0(paste0("I(sqrt(", X.cont), "))" )
X.all    <- c(X.names, X.square, X.sqroot)
rm(X.con, X.square, X.sqroot)

z.prob <- glm(paste("treat ~", paste(X.all, collapse = " + ")),
              family = binomial,
              data = jobs)$fitted.values

z.wt <- 1 / (jobs$treat*z.prob + (1-jobs$treat)*(1-z.prob))

tab1.wtd <- tableone::svyCreateTableOne(vars = X.names,
                                        strata = "treat",
                                        data = survey::svydesign(
                                            ids = ~1,
                                            weights = ~z.wt,
                                            data = cbind(jobs, z.wt = z.wt)),
                                        addOverall = TRUE,
                                        test = FALSE)
tab1.wtd

# for latex
tab1.wtd.matrix <- print(tab1.wtd, printToggle = FALSE, noSpaces = TRUE)
rownames(tab1.wtd.matrix) <- rownames(tab1.unw.matrix)
tab1.wtd.latex  <- xtable::xtable(tab1.wtd.matrix)


########################################
#### treatment group stratified by compliance type

tab1.treat.wtd <-
    tableone::svyCreateTableOne(vars = X.names,
                                strata = "complier",
                                data = survey::svydesign(
                                    ids = ~1,
                                    weights = ~z.wt,
                                    data = cbind(jobs, z.wt = z.wt)[jobs$treat==1,]),
                                test = FALSE)
tab1.treat.wtd

# for latex
tab1.treat.wtd.matrix <- print(tab1.treat.wtd, printToggle = FALSE, noSpaces = TRUE)
rownames(tab1.treat.wtd.matrix) <- rownames(tab1.unw.matrix)
tab1.treat.wtd.latex  <- xtable::xtable(tab1.treat.wtd.matrix)


########################################
#### control group stratified by work3 (one of the three outcomes)

tab1.control.wtd <-
    tableone::svyCreateTableOne(vars = X.names,
                                strata = "work3",
                                data = survey::svydesign(
                                    ids = ~1,
                                    weights = ~z.wt,
                                    data = cbind(jobs, z.wt = z.wt)[jobs$treat==0,]),
                                test = FALSE)
tab1.control.wtd

# for latex
tab1.control.wtd.matrix <- print(tab1.control.wtd, printToggle = FALSE, noSpaces = TRUE)
rownames(tab1.control.wtd.matrix) <- rownames(tab1.unw.matrix)
tab1.control.wtd.latex  <- xtable::xtable(tab1.control.wtd.matrix)

saveRDS(mget(c('tab1.unw.latex', 'tab1.wtd.latex', 'tab1.treat.wtd.latex', 'tab1.control.wtd.latex')),
        file = here::here(tables.dir, 'table1data.rds'))
