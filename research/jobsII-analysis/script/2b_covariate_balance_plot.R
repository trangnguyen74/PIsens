


library(tidyverse)
library(PIsens)


jobs.dir    <- here::here('research', 'jobsII-analysis')
script.dir  <- file.path(jobs.dir, 'script')
figs.dir <- file.path(jobs.dir, "outputs", "figures")


dat <- readRDS(file.path(jobs.dir, "data", "jobsII.rds"))

X.names  <- c("age", "sex", "race", "edu", "marital", "hh.kids",
              "hh.income", "econ.hard", "occu", "wks.unemp",
              "part.motiv", "seek.motiv", "seek.effi",
              "assertive", "depress1")

X.cont   <- c("age", "econ.hard", "wks.unemp", "part.motiv", "seek.motiv",
              "seek.effi", "assertive", "depress1")
X.square <- paste0(paste0("I("     , X.cont), "^2)")
X.sqroot <- paste0(paste0("I(sqrt(", X.cont), "))" )
X.all    <- c(X.names, X.square, X.sqroot)
rm(X.con, X.square, X.sqroot)


# inverse propensity weights
z.prob <- glm(paste("treat ~", paste(X.all, collapse = " + ")),
              family = binomial,
              data = dat)$fitted.values

dat$z.wt <- 1 / (dat$treat*z.prob + (1-dat$treat)*(1-z.prob))


# compliance type weights
c.mod <- glm(paste("complier ~", paste(X.all, collapse = " + ")),
             family = binomial,
             data = dat[dat$treat==1,])

dat$c.wt <- predict(c.mod, newdata = dat, type = "response")
dat$n.wt <- 1 - dat$c.wt


########################################
#### treated vs control balance
#### (inverse propensity score weighted)

z.balplot <-
    cobalt::love.plot(cobalt::bal.tab(x          = dat[X.names],
                                      treat      = dat$treat,
                                      weights    = dat$z.wt,
                                      continuous = "std",
                                      binary     = "raw",
                                      s.d.denom  = "pooled",
                                      un         = TRUE),
                      stars        = "std",
                      grid         = TRUE,
                      wrap         = 20,
                      sample.names = c("unweighted", "prop.score weighted"),
                      position     = "bottom",
                      shapes       = c("circle", "triangle"),
                      colors       = c("red", "blue")) +
    labs(x     = "mean differences (standardized if *)",
         title = "treated vs. control") +
    theme(legend.title       = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_continuous(limits = c(-.4, .4),
                       breaks = c(-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4))


########################################
#### treated compliers vs control balance
#### (inverse propensity score weighted + principal score weighted)

dat.c    <- dat[!(!is.na(dat$complier) & dat$complier==0),]
dat.c$wt <- dat.c$treat*dat.c$z.wt + (1-dat.c$treat)*dat.c$z.wt*dat.c$c.wt

zc.balplot <-
    cobalt::love.plot(cobalt::bal.tab(x          = dat.c[X.names],
                                      treat      = dat.c$treat,
                                      weights    = dat.c$wt,
                                      continuous = "std",
                                      binary     = "raw",
                                      s.d.denom  = "pooled",
                                      un         = TRUE),
                      stars        = "std",
                      grid         = TRUE,
                      wrap         = 20,
                      sample.names = c("unweighted",
                                       "prop.score & prin.score weighted"),
                      position     = "bottom",
                      shapes       = c("circle", "triangle"),
                      colors       = c("red", "blue")) +
    labs(x     = "mean differences (standardized if *)",
         title = "treated compliers vs. control") +
    theme(legend.title       = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_continuous(limits = c(-.4, .4),
                       breaks = c(-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4))




########################################
#### treated noncompliers vs control balance
#### (inverse propensity score weighted + principal score weighted)

dat.n    <- dat[!(!is.na(dat$complier) & dat$complier==1),]
dat.n$wt <- dat.n$treat*dat.n$z.wt + (1-dat.n$treat)*dat.n$z.wt*dat.n$n.wt

zn.balplot <-
    cobalt::love.plot(cobalt::bal.tab(x          = dat.n[X.names],
                                      treat      = dat.n$treat,
                                      weights    = dat.n$wt,
                                      continuous = "std",
                                      binary     = "raw",
                                      s.d.denom  = "pooled",
                                      un         = TRUE),
                      stars        = "std",
                      grid         = TRUE,
                      wrap         = 20,
                      sample.names = c("unweighted",
                                       "prop.score & prin.score weighted"),
                      position     = "bottom",
                      shapes       = c("circle", "triangle"),
                      colors       = c("red", "blue")) +
    labs(x     = "mean differences (standardized if *)",
         title = "treated noncompliers vs. control") +
    theme(legend.title       = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_continuous(limits = c(-.4, .4),
                       breaks = c(-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4))


pdf(file.path(figs.dir, 'fig4.pdf'), width = 14, height = 8)
cowplot::plot_grid(z.balplot, zc.balplot, zn.balplot, nrow = 1)
dev.off()

























