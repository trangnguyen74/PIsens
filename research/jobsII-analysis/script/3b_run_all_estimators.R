

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

sens.dir <- here::here("research", "jobsII-analysis", "outputs", "results", "all-estimators", "sens")


##################################################
# PART 1: work for pay


dat.w <- declare_data(data   = jobsII,
                      z.var  = "treat",
                      c.var  = "complier",
                      y.var  = "work3",
                      x.vars = X.names,
                      id     = "id")

pimu0 <- PIsens(data        = dat.w,
                estimator   = "pimu",
                targeted    = FALSE,
                C.form      = C.form,
                Y.form      = Y.form,
                Y.model     = "quasibinomial",
                sens.type   = "OR",
                sens.range  = c(1/3, 3),
                sens.step   = .05,
                boot.seed   = 12345,
                double.boot = TRUE)

pimu1 <- PIsens(data        = dat.w,
                estimator   = "pimu",
                Z.form      = Z.form,
                C.form      = C.form,
                Y.form      = Y.form,
                Y.model     = "quasibinomial",
                sens.type   = "OR",
                sens.range  = c(1/3, 3),
                sens.step   = .05,
                boot.seed   = 12345,
                double.boot = TRUE)

emu0 <- PIsens(data        = dat.w,
               estimator   = "emu",
               targeted    = FALSE,
               Z.form      = Z.form,
               C.form      = C.form,
               Y.form      = Y.form,
               Y.model     = "quasibinomial",
               sens.type   = "OR",
               sens.range  = c(1/3, 3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

emu1 <- PIsens(data        = dat.w,
               estimator   = "emu",
               Z.form      = Z.form,
               C.form      = C.form,
               Y.form      = Y.form,
               Y.model     = "quasibinomial",
               sens.type   = "OR",
               sens.range  = c(1/3, 3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

MS1 <- PIsens(data        = dat.w,
              estimator   = "MS",
              Z.form      = Z.form,
              C.form      = C.form,
              Y.form      = Y.form,
              Y.model     = "quasibinomial",
              sens.type   = "OR",
              sens.range  = c(1/3, 3),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IF0 <- PIsens(data        = dat.w,
              estimator   = "IF",
              targeted    = FALSE,
              Z.form      = Z.form,
              C.form      = C.form,
              Y.form      = Y.form,
              Y.model     = "quasibinomial",
              sens.type   = "OR",
              sens.range  = c(1/3, 3),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IF1 <- PIsens(data        = dat.w,
              estimator   = "IF",
              Z.form      = Z.form,
              C.form      = C.form,
              Y.form      = Y.form,
              Y.model     = "quasibinomial",
              sens.type   = "OR",
              sens.range  = c(1/3, 3),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IFH0 <- PIsens(data        = dat.w,
               estimator   = "IFH",
               targeted    = FALSE,
               Z.form      = Z.form,
               C.form      = C.form,
               Y.form      = Y.form,
               Y.model     = "quasibinomial",
               sens.type   = "OR",
               sens.range  = c(1/3, 3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

IFH1 <- PIsens(data        = dat.w,
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

saveRDS(mget(c("pimu0", "pimu1", "emu0", "emu1", "MS1", "IF0", "IF1", "IFH0", "IFH1")),
        file = file.path(sens.dir, "work_OR.rds"))

rm(dat.w, pimu0, pimu1, emu0, emu1, MS1, IF0, IF1, IFH0, IFH1)



##################################################
# PART 2: depressive symptoms, GOR-based sens analysis


dat.d <- declare_data(data   = jobsII,
                      z.var  = "treat",
                      c.var  = "complier",
                      y.var  = "depress3",
                      x.vars = X.names,
                      id     = "id")

pimu0 <- PIsens(data        = dat.d,
                estimator   = "pimu",
                targeted    = FALSE,
                C.form      = C.form,
                Y.form      = Y.form,
                Y.model     = "quasibinomial",
                Y.bounds    = c(1,5),
                sens.type   = "GOR",
                sens.range  = c(1/3, 3),
                sens.step   = .05,
                boot.seed   = 12345,
                double.boot = TRUE)

pimu1 <- PIsens(data        = dat.d,
                estimator   = "pimu",
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

emu0 <- PIsens(data        = dat.d,
               estimator   = "emu",
               targeted    = FALSE,
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

emu1 <- PIsens(data        = dat.d,
               estimator   = "emu",
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

MS1 <- PIsens(data        = dat.d,
              estimator   = "MS",
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

IF0 <- PIsens(data        = dat.d,
              estimator   = "IF",
              targeted    = FALSE,
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

IF1 <- PIsens(data        = dat.d,
              estimator   = "IF",
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

IFH0 <- PIsens(data        = dat.d,
               estimator   = "IFH",
               targeted    = FALSE,
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

IFH1 <- PIsens(data        = dat.d,
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

saveRDS(mget(c("pimu0", "pimu1", "emu0", "emu1", "MS1", "IF0", "IF1", "IFH0", "IFH1")),
        file = file.path(sens.dir, "depress_GOR.rds"))

rm(pimu0, pimu1, emu0, emu1, MS1, IF0, IF1, IFH0, IFH1)



##################################################
# PART 3: depressive symptoms, SMDe-based sens analysis


pimu0 <- PIsens(data        = dat.d,
                estimator   = "pimu",
                targeted    = FALSE,
                C.form      = C.form,
                Y.form      = Y.form,
                Y.model     = "quasibinomial",
                Y.bounds    = c(1,5),
                sens.type   = "SMDe",
                sens.range  = c(-1, 1),
                sens.step   = .05,
                boot.seed   = 12345,
                double.boot = TRUE)

pimu1 <- PIsens(data        = dat.d,
                estimator   = "pimu",
                Z.form      = Z.form,
                C.form      = C.form,
                Y.form      = Y.form,
                Y.model     = "quasibinomial",
                Y.bounds    = c(1,5),
                sens.type   = "SMDe",
                sens.range  = c(-1, 1),
                sens.step   = .05,
                boot.seed   = 12345,
                double.boot = TRUE)

emu0 <- PIsens(data        = dat.d,
               estimator   = "emu",
               targeted    = FALSE,
               Z.form      = Z.form,
               C.form      = C.form,
               Y.form      = Y.form,
               Y.model     = "quasibinomial",
               Y.bounds    = c(1,5),
               sens.type   = "SMDe",
               sens.range  = c(-1, 1),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

emu1 <- PIsens(data        = dat.d,
               estimator   = "emu",
               Z.form      = Z.form,
               C.form      = C.form,
               Y.form      = Y.form,
               Y.model     = "quasibinomial",
               Y.bounds    = c(1,5),
               sens.type   = "SMDe",
               sens.range  = c(-1, 1),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

MS1 <- PIsens(data        = dat.d,
              estimator   = "MS",
              Z.form      = Z.form,
              C.form      = C.form,
              Y.form      = Y.form,
              Y.model     = "quasibinomial",
              Y.bounds    = c(1,5),
              sens.type   = "SMDe",
              sens.range  = c(-1, 1),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IF0 <- PIsens(data        = dat.d,
              estimator   = "IF",
              targeted    = FALSE,
              Z.form      = Z.form,
              C.form      = C.form,
              Y.form      = Y.form,
              Y.model     = "quasibinomial",
              Y.bounds    = c(1,5),
              sens.type   = "SMDe",
              sens.range  = c(-1, 1),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IF1 <- PIsens(data        = dat.d,
              estimator   = "IF",
              Z.form      = Z.form,
              C.form      = C.form,
              Y.form      = Y.form,
              Y.model     = "quasibinomial",
              Y.bounds    = c(1,5),
              sens.type   = "SMDe",
              sens.range  = c(-1, 1),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IFH0 <- PIsens(data        = dat.d,
               estimator   = "IFH",
               targeted    = FALSE,
               Z.form      = Z.form,
               C.form      = C.form,
               Y.form      = Y.form,
               Y.model     = "quasibinomial",
               Y.bounds    = c(1,5),
               sens.type   = "SMDe",
               sens.range  = c(-1, 1),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

IFH1 <- PIsens(data        = dat.d,
               estimator   = "IFH",
               Z.form      = Z.form,
               C.form      = C.form,
               Y.form      = Y.form,
               Y.model     = "quasibinomial",
               Y.bounds    = c(1,5),
               sens.type   = "SMDe",
               sens.range  = c(-1, 1),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

saveRDS(mget(c("pimu0", "pimu1", "emu0", "emu1", "MS1", "IF0", "IF1", "IFH0", "IFH1")),
        file = file.path(sens.dir, "depress_SMDe.rds"))

rm(pimu0, pimu1, emu0, emu1, MS1, IF0, IF1, IFH0, IFH1)




##################################################
# PART 4: earnings


dat.e <- declare_data(data       = jobsII,
                      z.var      = "treat",
                      c.var      = "complier",
                      y.var      = "earn3",
                      x.vars     = X.names,
                      other.vars = "work3",
                      id         = "id")

earn3_0cn.1cn_targeted0 <- function(data) {
    earn3_FUN1(data = data, targeted = FALSE, include.mu.1c = TRUE)
}
earn3_0cn.1cn_targeted1 <- function(data) {
    earn3_FUN1(data = data, targeted = TRUE, include.mu.1c = TRUE)
}

earn3_0cn_targeted0 <- function(data) {
    earn3_FUN1(data = data, targeted = FALSE, include.mu.1c = FALSE)
}
earn3_0cn_targeted1 <- function(data) {
    earn3_FUN1(data = data, targeted = TRUE, include.mu.1c = FALSE)
}


# TODO: broken code. epi estimator has error: I cannot find Y.form.

epi0 <- PIsens(data        = dat.e,
               estimator   = "epi",
               targeted    = FALSE,
               Z.form      = Z.form,
               C.form      = C.form,
               sens.type   = "MR",
               sens.range  = c(1/3,3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

epi1 <- PIsens(data        = dat.e,
               estimator   = "epi",
               Z.form      = Z.form,
               C.form      = C.form,
               sens.type   = "MR",
               sens.range  = c(1/3,3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)


pimu0 <- PIsens(data        = dat.e,
                estimator   = "pimu",
                targeted    = FALSE,
                C.form      = C.form,
                Y.FUN       = earn3_0cn.1cn_targeted0,
                sens.type   = "MR",
                sens.range  = c(1/3,3),
                sens.step   = .05,
                boot.seed   = 12345,
                double.boot = TRUE)

pimu1 <- PIsens(data        = dat.e,
                estimator   = "pimu",
                Z.form      = Z.form,
                C.form      = C.form,
                Y.FUN       = earn3_0cn.1cn_targeted1,
                sens.type   = "MR",
                sens.range  = c(1/3,3),
                sens.step   = .05,
                boot.seed   = 12345,
                double.boot = TRUE)

emu0 <- PIsens(data        = dat.e,
               estimator   = "emu",
               targeted    = FALSE,
               Z.form      = Z.form,
               C.form      = C.form,
               Y.FUN       = earn3_0cn_targeted0,
               sens.type   = "MR",
               sens.range  = c(1/3,3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

emu1 <- PIsens(data        = dat.e,
               estimator   = "emu",
               Z.form      = Z.form,
               C.form      = C.form,
               Y.FUN       = earn3_0cn_targeted1,
               sens.type   = "MR",
               sens.range  = c(1/3,3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

MS1 <- PIsens(data        = dat.e,
              estimator   = "MS",
              Z.form      = Z.form,
              C.form      = C.form,
              Y.FUN       = earn3_0cn.1cn_targeted1,
              sens.type   = "MR",
              sens.range  = c(1/3,3),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IF0 <- PIsens(data        = dat.e,
              estimator   = "IF",
              targeted    = FALSE,
              Z.form      = Z.form,
              C.form      = C.form,
              Y.FUN       = earn3_0cn.1cn_targeted0,
              sens.type   = "MR",
              sens.range  = c(1/3,3),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IF1 <- PIsens(data        = dat.e,
              estimator   = "IF",
              Z.form      = Z.form,
              C.form      = C.form,
              Y.FUN       = earn3_0cn.1cn_targeted1,
              sens.type   = "MR",
              sens.range  = c(1/3,3),
              sens.step   = .05,
              boot.seed   = 12345,
              double.boot = TRUE)

IFH0 <- PIsens(data        = dat.e,
               estimator   = "IFH",
               targeted    = FALSE,
               Z.form      = Z.form,
               C.form      = C.form,
               Y.FUN       = earn3_0cn.1cn_targeted0,
               sens.type   = "MR",
               sens.range  = c(1/3,3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

IFH1 <- PIsens(data        = dat.e,
               estimator   = "IFH",
               Z.form      = Z.form,
               C.form      = C.form,
               Y.FUN       = earn3_0cn.1cn_targeted1,
               sens.type   = "MR",
               sens.range  = c(1/3,3),
               sens.step   = .05,
               boot.seed   = 12345,
               double.boot = TRUE)

saveRDS(mget(c("epi0", "epi1", "pimu0", "pimu1", "emu0", "emu1", "MS1", "IF0", "IF1", "IFH0", "IFH1")),
        file = file.path(sens.dir, "earn_MR.rds"))

rm(epi0, epi1, pimu0, pimu1, emu0, emu1, MS1, IF0, IF1, IFH0, IFH1)





