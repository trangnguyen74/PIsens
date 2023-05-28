


library(tidyverse)
library(PIsens)


jobs.dir    <- here::here('research', 'jobsII-analysis')
figs.dir <- file.path(jobs.dir, "outputs", "figures")


dat <- readRDS(file.path(jobs.dir, "data", "jobsII.rds"))

X.names  <- c("age", "sex", "race", "edu", "marital", "hh.kids",
              "hh.income", "econ.hard", "occu", "wks.unemp",
              "part.motiv", "seek.motiv", "seek.effi",
              "assertive", "depress1")

dat <- declare_data(data = dat,
                    z.var = "treat",
                    c.var = "complier",
                    y.var = "earn3",
                    x.vars = X.names,
                    other.vars = "work3")



X.cont   <- c("age", "econ.hard", "wks.unemp", "part.motiv", "seek.motiv",
              "seek.effi", "assertive", "depress1")
X.all    <- c(X.names,
              paste0(paste0("I("     , X.cont), "^2)"),
              paste0(paste0("I(sqrt(", X.cont), "))" ))
rm(X.cont)

Z.form <- paste("z~", paste(X.all, collapse = " + "))
C.form <- paste("z~", paste(X.all, collapse = " + "))
Y.FUN  <- function(data) {
    earn3_FUN1(data = data, targeted = TRUE, include.mu.1c = FALSE)
}



pdf(file.path(figs.dir, 'fig5.pdf'), width = 12, height = 3.5)
visualize_MR_implication(MRs = c(1/4, 1/3, 1/2, 1, 2, 3, 4),
                         MRlabs = c('1/4', '1/3', '1/2', '1', '2', '3', '4'),
                         data = dat,
                         estimator = "IFH",
                         targeted = TRUE,
                         Z.form = Z.form,
                         C.form = C.form,
                         Y.FUN  = Y.FUN)
dev.off()
