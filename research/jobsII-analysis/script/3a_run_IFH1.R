

library(tidyverse)
library(PIsens)

jobsII <- readRDS(here::here("data-raw", "dat.rds"))

X.names <- c("age", "sex", "race", "edu", "marital", "hh.kids",
             "hh.income", "econ.hard", "occu", "wks.unemp",
             "part.motiv", "seek.motiv", "seek.effi",
             "assertive", "depress1")

X.cont <- c("age", "econ.hard", "wks.unemp", "part.motiv", "seek.motiv",
            "seek.effi", "assertive", "depress1")
# squares and square roots of continuous covariates
X.square <- paste0(paste0("I("     , X.cont), "^2)")
X.sqroot <- paste0(paste0("I(sqrt(", X.cont), "))" )
X.all <- c(X.names, X.square, X.sqroot)
rm(X.cont, X.square, X.sqroot)

Z.form <- paste("z ~ ", paste(X.all,   collapse = " + "))
C.form <- paste("c ~ ", paste(X.all,   collapse = " + "))
Y.form <- paste("y ~ ", paste(X.names, collapse = " + "))

sens.dir <- here::here("research", "jobsII-analysis", "outputs", "results", "IFH1-only", "sens")
main.dir <- here::here("research", "jobsII-analysis", "outputs", "results", "IFH1-only", "main")


##################################################
# outcome 1: work for pay


# declare data
dat.w <- declare_data(data   = jobsII,
                      z.var  = "treat",
                      c.var  = "complier",
                      y.var  = "work3",
                      x.vars = X.names,
                      id     = "id")

# main analysis
# TODO: code is currently broken (.clean_nuis_inputs requires sens.type)
work.PI <- PImain(data      = dat.w,
                  estimator = "IFH",
                  Z.form    = Z.form,
                  C.form    = C.form,
                  Y.form    = Y.form,
                  Y.model   = "quasibinomial",
                  boot.seed = 12345)

# OR-based sens analysis
work.OR <- PIsens(data        = dat.w,
                  estimator   = "IFH",
                  Z.form      = Z.form,
                  C.form      = C.form,
                  Y.form      = Y.form,
                  Y.model     = "quasibinomial",
                  sens.type   = "OR",
                  sens.range  = c(1/3, 3),
                  sens.step   = .05,
                  boot.seed   = 12345,
                  double.boot = TRUE)

saveRDS(work.OR,
        file = file.path(sens.dir, "work_OR.rds"))




##################################################
# outcome 2: depressive symptoms


# declare data
dat.d <- declare_data(data   = jobsII,
                      z.var  = "treat",
                      c.var  = "complier",
                      y.var  = "depress3",
                      x.vars = X.names,
                      id     = "id")

# main analysis
# TODO: code is currently broken (.clean_nuis_inputs requires sens.type)
depress.PI <- PImain(data      = dat.d,
                     estimator = "IFH",
                     Z.form    = Z.form,
                     C.form    = C.form,
                     Y.form    = Y.form,
                     Y.model   = "quasibinomial",
                     Y.bounds  = c(1,5),
                     boot.seed = 12345)

# GOR-based sens analysis
depress.GOR <- PIsens(data        = dat.d,
                      estimator   = "IFH",
                      Z.form      = Z.form,
                      C.form      = C.form,
                      Y.form      = Y.form,
                      Y.model     = "quasibinomial",
                      Y.bounds    = c(1,5),
                      sens.type   = "GOR",
                      sens.range  = c(1/3, 1.25),
                      sens.step   = .05,
                      boot.seed   = 12345,
                      double.boot = TRUE)

saveRDS(depress.GOR,
        file = file.path(sens.dir, "depress_GOR.rds"))

# SMDe-based sens analysis
depress.SMDe <- PIsens(data        = dat.d,
                       estimator   = "IFH",
                       Z.form      = Z.form,
                       C.form      = C.form,
                       Y.form      = Y.form,
                       Y.model     = "quasibinomial",
                       Y.bounds    = c(1,5),
                       sens.type   = "SMDe",
                       sens.range  = c(-1,.1),
                       sens.step   = .05,
                       boot.seed   = 12345,
                       double.boot = TRUE)

saveRDS(depress.SMDe,
        file = file.path(sens.dir, "depress_SMDe.rds"))


# GOR-based sens analysis -- symmetric (for finite-sample bias investigation)
depress.GOR.symmetric <- PIsens(data        = dat.d,
                                estimator   = "IFH",
                                Z.form      = Z.form,
                                C.form      = C.form,
                                Y.form      = Y.form,
                                Y.model     = "quasibinomial",
                                Y.bounds    = c(1,5),
                                sens.type   = "GOR",
                                sens.range  = c(1/3, 3),
                                sens.step   = .05,
                                boot.seed   = 12345,
                                double.boot = TRUE)

saveRDS(depress.GOR.symmetric,
        file = file.path(sens.dir, "depress_GOR_symmetric.rds"))

# SMDe-based sens analysis -- symmetric (for finite-sample bias investigation)
depress.SMDe.symmetric <- PIsens(data        = dat.d,
                                 estimator   = "IFH",
                                 Z.form      = Z.form,
                                 C.form      = C.form,
                                 Y.form      = Y.form,
                                 Y.model     = "quasibinomial",
                                 Y.bounds    = c(1,5),
                                 sens.type   = "SMDe",
                                 sens.range  = c(-1,1),
                                 sens.step   = .05,
                                 boot.seed   = 12345,
                                 double.boot = TRUE)

saveRDS(depress.SMDe.symmetric,
        file = file.path(sens.dir, "depress_SMDe_symmetric.rds"))



##################################################
# outcome 3: earnings


# declare data
dat.e <- declare_data(data       = jobsII,
                      z.var      = "treat",
                      c.var      = "complier",
                      y.var      = "earn3",
                      x.vars     = X.names,
                      other.vars = "work3",
                      id         = "id")


# set up function to estimate conditional outcome means
earn3_0cn.1cn_targeted1 <- function(data) {
    earn3_FUN1(data = data, targeted = TRUE, include.mu.1c = TRUE)
}



# main analysis
# TODO: code is currently broken (.clean_nuis_inputs requires sens.type)
earn.PI <- PImain(data      = dat.e,
                  estimator = "IFH",
                  Z.form    = Z.form,
                  C.form    = C.form,
                  Y.FUN     = earn3_0cn.1cn_targeted1,
                  boot.seed = 12345)

# MR-based sens analysis
earn.MR <- PIsens(data        = dat.e,
                  estimator   = "IFH",
                  Z.form      = Z.form,
                  C.form      = C.form,
                  Y.FUN       = earn3_0cn.1cn_targeted1,
                  sens.type   = "MR",
                  sens.range  = c(1/3,3),
                  sens.step   = .05,
                  boot.seed   = 12345,
                  double.boot = TRUE)

saveRDS(earn.MR,
        file = file.path(sens.dir, "earn_MR.rds"))




























