
########################################
# PREPARATIONS

X.names <- c("age", "sex", "race", "edu", "marital", "hh.kids",
             "hh.income", "econ.hard", "occu", "wks.unemp",
             "part.motiv", "seek.motiv", "seek.effi",
             "assertive", "depress1")

X.cont <- c("age", "econ.hard", "wks.unemp", "part.motiv", "seek.motiv",
            "seek.effi", "assertive", "depress1")
X.square <- paste0(paste0("I("     , X.cont), "^2)")
X.sqroot <- paste0(paste0("I(sqrt(", X.cont), "))" )
rm(X.cont)

Z.form <- paste("z ~ ", paste(c(X.names, X.square, X.sqroot), collapse = " + "))
C.form <- paste("c ~ ", paste(c(X.names, X.square, X.sqroot), collapse = " + "))
Y.form <- paste("y ~ ", paste(X.names,                        collapse = " + "))


# JOBS II data
dat <- readRDS(here::here("data-raw", "dat.rds"))

# base data for simulating (X,Z)
base       <- dat[c(X.names, "treat")]
base$treat <- factor(base$treat)

# models for simulating (C,Y)
c.form   <- paste("complier ~", paste(c(X.names, X.square, X.sqroot), collapse = " + "))
y.form.w <- paste("work3 ~",    paste(X.names, collapse = " + "))
y.form.e <- paste("earn3 ~",    paste(X.names, collapse = " + "))



c.mod    <- glm(formula = c.form,   family = "binomial", data = dat[dat$treat==1,])
y0.mod.w  <- glm(formula = y.form.w, family = "binomial", data = dat[dat$treat==0,])
y1c.mod.w <- glm(formula = y.form.w, family = "binomial", data = dat[dat$treat==1 & dat$complier==1,])
y1n.mod.w <- glm(formula = y.form.w, family = "binomial", data = dat[dat$treat==1 & dat$complier==0,])

y0.mod.e  <- glm(formula = y.form.e, family = Gamma(link="log"), data = dat[dat$treat==0 & dat$work3==1,])
y1c.mod.e <- glm(formula = y.form.e, family = Gamma(link="log"), data = dat[dat$treat==1 & dat$work3==1 & dat$complier==1,])
y1n.mod.e <- glm(formula = y.form.e, family = Gamma(link="log"), data = dat[dat$treat==1 & dat$work3==1 & dat$complier==0,])

y0.alpha.e  <- MASS::gamma.shape(y0.mod.e)[[1]]
y1c.alpha.e <- MASS::gamma.shape(y1c.mod.e)[[1]]
y1n.alpha.e <- MASS::gamma.shape(y1n.mod.e)[[1]]

y.max.e <- max(dat$earn3)


my_rtgamma <- function(n, shape, scale, truncation) {

    x <- rgamma(n = n, shape = shape, scale = scale)

    while(sum(x>truncation) > 0) {

        y <- rgamma(n = n, shape = shape, scale = scale)

        x <- ifelse(x>truncation, y, x)
    }

    x
}



sim_earn3 <- function(xz.base     = base,
                      c.model     = c.mod,
                      y0.model.w  = y0.mod.w,
                      y1c.model.w = y1c.mod.w,
                      y1n.model.w = y1n.mod.w,
                      y0.model.e  = y0.mod.e,
                      y1c.model.e = y1c.mod.e,
                      y1n.model.e = y1n.mod.e,
                      y0.shape.e  = y0.alpha.e,
                      y1c.shape.e = y1c.alpha.e,
                      y1n.shape.e = y1n.alpha.e,
                      y.bound.e   = y.max.e,
                      size,
                      seed) {

    sim <- synthpop::syn(xz.base, seed = seed, k = size)$syn
    sim$treat <- as.numeric(sim$treat) - 1

    # simulate complier
    c.sim <- predict(c.model, newdata = sim, type = "response")
    c.sim <- rbinom(n = nrow(sim), size = 1, prob = c.sim)

    sim$complier <- ifelse(sim$treat==0, NA, c.sim)

    # simulate work3 before earn3
    y0.sim  <- predict(y0.model.w,  newdata = sim, type = "response")
    y1c.sim <- predict(y1c.model.w, newdata = sim, type = "response")
    y1n.sim <- predict(y1n.model.w, newdata = sim, type = "response")

    y0.sim  <- rbinom(n = nrow(sim), size = 1, prob = y0.sim)
    y1c.sim <- rbinom(n = nrow(sim), size = 1, prob = y1c.sim)
    y1n.sim <- rbinom(n = nrow(sim), size = 1, prob = y1n.sim)

    sim$work3 <- ifelse(sim$treat==0, y0.sim,
                        ifelse(sim$complier==1, y1c.sim, y1n.sim))

    rm(y0.sim, y1c.sim, y1n.sim)

    # simulate earn3
    y0.scale  <- predict(y0.model.e,  newdata = sim, type = "response")/y0.shape.e
    y1c.scale <- predict(y1c.model.e, newdata = sim, type = "response")/y1c.shape.e
    y1n.scale <- predict(y1n.model.e, newdata = sim, type = "response")/y1n.shape.e

    y0.sim  <- my_rtgamma(n = nrow(sim), shape = y0.shape.e,  scale = y0.scale,  truncation = y.bound.e)
    y1c.sim <- my_rtgamma(n = nrow(sim), shape = y1c.shape.e, scale = y1c.scale, truncation = y.bound.e)
    y1n.sim <- my_rtgamma(n = nrow(sim), shape = y1n.shape.e, scale = y1n.scale, truncation = y.bound.e)

    sim$earn3 <- ifelse(sim$work3==0, 0,
                        ifelse(sim$treat==0, y0.sim,
                               ifelse(sim$complier==1, y1c.sim, y1n.sim)))

    sim
}






########################################
# APPROXIMATE TRUE VALUE

true.num  <- 10
true.seed <- 123
true.size <- 500000


set.seed(true.seed)
true.seeds <- sample(1:.Machine$integer.max, true.num)

true.values <- list()


for (z in 1:true.num) {

    pop <- sim_earn3(size = true.size,
                     seed = true.seeds[z])

    pop <- declare_data(pop,
                        z.var      = "treat",
                        c.var      = "complier",
                        y.var      = "earn3",
                        x.vars     = X.names,
                        other.vars = "work3")

    true.values[[z]] <- tryCatch(PIsens(data       = pop,
                                        estimator  = "IFH",
                                        Z.form     = Z.form,
                                        C.form     = C.form,
                                        Y.FUN      = earn3_FUN1_IFH,
                                        sens.type  = "MR",
                                        sens.range = c(1/3,3),
                                        sens.step  = .05,
                                        boot.num   = 0)$sens,
                                 error = function(e) return(NULL))

    rm(pop)
    gc()
    gc()
    gc()
}


true.values <- true.values[lengths(true.values) != 0]
true.ok     <- length(true.values)
true.values <- simplify2array(true.values)

sens.params <- true.values[,1:2,1]
true.values <- true.values[,-c(1:2),]

true.se     <- apply(true.values, c(1,2), function(z) sd(z, na.rm = TRUE))
true.values <- apply(true.values, c(1,2), function(z) mean(z, na.rm = TRUE))


earn.true <- list(sens.params = sens.params,
                  values      = true.values,
                  se          = true.se,
                  size        = true.size,
                  num         = true.num,
                  seed        = true.seed,
                  ok          = true.ok)

rm(true.num, true.seed, true.size, true.seeds, true.se, true.ok, z)

saveRDS(list(true = earn.true),
        file = here::here("simulation", "earn_sim.rds"))




########################################
# SIMULATION

sim.num  <- 500
sim.seed <- 987
sim.size <- 465 # size of analytic data


set.seed(sim.seed)
sim.seeds <- sample(1:.Machine$integer.max, sim.num)


sens.params <- earn.true$sens.params
true.values <- earn.true$values

error0                <- array(dim = c(dim(true.values), sim.num))
dimnames(error0)[[2]] <- colnames(true.values)
error1 <- error2      <- error0


for (z in 162:500) {

    set.seed(sim.seeds[z])

    tryCatch({

        sim <- sim_earn3(size = sim.size,
                         seed = sim.seeds[z])

        sim <- declare_data(sim,
                            z.var  = "treat",
                            c.var  = "complier",
                            y.var  = "earn3",
                            x.vars = X.names,
                            other.vars = "work3")

        sens <- PIsens(data        = sim,
                       estimator   = "IFH",
                       Z.form      = Z.form,
                       C.form      = C.form,
                       Y.FUN       = earn3_FUN1_IFH,
                       sens.type   = "MR",
                       sens.range  = c(1/3,3),
                       sens.step   = .05,
                       boot.seed   = 12345,
                       BCa.CI      = FALSE,
                       double.boot = TRUE)$sens

        sens <- simplify2array(sens)

        error0[,,z] <- sens[,"point",]     - true.values
        error1[,,z] <- sens[,"point.bc1",] - true.values
        error2[,,z] <- sens[,"point.bc2",] - true.values

        rm(sim, sens)

        saveRDS(list(true = earn.true,
                     sim = list(sens.params = sens.params,
                                errors      = list(bc0 = error0,
                                                   bc1 = error1,
                                                   bc2 = error2),
                                size        = sim.size,
                                num         = sim.num,
                                seed        = sim.seed,
                                ok          = sum(!is.na(error0[1,1,])))),
                file = here::here("simulation", "earn_sim.rds"))
    })

    gc()
}




earn.sim <- list(sens.params = sens.params,
                 errors      = list(bc0 = error0,
                                    bc1 = error1,
                                    bc2 = error2),
                 size        = sim.size,
                 num         = sim.num,
                 seed        = sim.seed,
                 ok          = sum(!is.na(error0[1,1,])))

rm(sim.num, sim.seed, sim.size, sim.seeds, z,
   sens.params, true.values,
   error0, error1, error2)




pdf(here::here("simulation", "earn_sim_figure.pdf"), width = 9, height = 9)
plot_sim_result(sim         = readRDS(here::here("simulation", "earn_sim.rds"))$sim,
                sens.type   = "MR",
                ylimits.std = c(-1,1))
dev.off()



















