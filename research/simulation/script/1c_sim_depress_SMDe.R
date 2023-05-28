
########################################
# PREPARATIONS

base.dir <- here::here("research",
                       "simulation")

script.dir  <- file.path(base.dir, "script")
sim.dir     <- file.path(base.dir, "outputs", "depress", "SMDe")
figs.dir <- file.path(base.dir, "outputs", "figures")

source(file.path(script.dir, "functions.R"))

true.params <- ..get_true_params()
sim.params  <- ..get_sim_params()
model.forms <- ..get_jobsII_formulas()
Y.bounds    <- c(1,5)
sim.basis   <- ..get_sim_basis_depress()




########################################
# APPROXIMATE TRUE VALUE

true.values <- array(dim = c(41, 6, true.params$master$num))



for (z in 1:dim(true.values)[3]) {

    cat(paste("\n\n\niteration", z, "\n\n\n"))

    pop <- do.call("..sim_depress",
                   args = c(sim.basis,
                            size = true.params$master$size,
                            seed = true.params$seeds[z]))

    .cleanMem(5)

    pop <- declare_data(pop,
                        z.var  = "treat",
                        c.var  = "complier",
                        y.var  = "depress3",
                        x.vars = model.forms$X.names)

    if (z==1) {
        tmp <- tryCatch(PIsens(data       = pop,
                               estimator  = "IFH",
                               Z.form     = model.forms$Z.form,
                               C.form     = model.forms$C.form,
                               Y.form     = model.forms$Y.form,
                               Y.model    = "quasibinomial",
                               Y.bounds   = Y.bounds,
                               sens.type  = "SMDe",
                               sens.range = c(-1, 1),
                               sens.step  = .05,
                               boot.num   = 0)$sens,
                        error = function(e) return(NA))

        sens.params <- tmp[,1:2]
        true.values[,,z] <- tmp[,-c(1:2)]
        dimnames(true.values)[[2]] <- colnames(tmp)[-c(1:2)]

        rm(tmp)

    } else {
        true.values[,,z] <- tryCatch(PIsens(data       = pop,
                                            estimator  = "IFH",
                                            Z.form     = model.forms$Z.form,
                                            C.form     = model.forms$C.form,
                                            Y.form     = model.forms$Y.form,
                                            Y.model    = "quasibinomial",
                                            Y.bounds   = Y.bounds,
                                            sens.type  = "SMDe",
                                            sens.range = c(-1, 1),
                                            sens.step  = .05,
                                            boot.num   = 0)$sens[,-c(1:2)],
                                     error = function(e) return(NA))
    }

    rm(pop)

    .cleanMem(5)

    saveRDS(list(master      = true.params$master,
                 sens.params = sens.params,
                 true.values = true.values),
            file = file.path(sim.dir, "true.rds"))
}


true.ok     <- sum(!is.na(true.values[1,1,]))
true.se     <- apply(true.values, c(1,2), function(z) sd(z, na.rm = TRUE)) / sqrt(true.params$master$num)
true.values <- apply(true.values, c(1,2), function(z) mean(z, na.rm = TRUE))


saveRDS(list(master      = true.params$master,
             sens.params = sens.params,
             values      = true.values,
             se          = true.se,
             ok          = true.ok),
        file = file.path(sim.dir, "true.rds"))



rm(true.params, true.se, true.values, true.ok, sens.params, z)

.cleanMem(5)





########################################
# SIMULATION

true.values <- readRDS(file.path(sim.dir, "true.rds"))
sens.params <- true.values$sens.params
true.values <- true.values$values


raw.est                <- array(dim = c(dim(true.values), sim.params$master$num))
dimnames(raw.est)[[2]] <- colnames(true.values)
boot.mean1 <- boot.mean2  <- raw.est


for (z in 1:dim(raw.est)[3]) {

    cat(paste("\n\n\nsimulation", z, "\n\n\n"))

    sim <- tryCatch(do.call("..sim_depress",
                            args = c(sim.basis,
                                     size = sim.params$master$size,
                                     seed = sim.params$seeds[z])),
                    error = function(e) return(NULL))

    .cleanMem(2)

    if (!is.null(sim)) {

        sim <- declare_data(sim,
                            z.var  = "treat",
                            c.var  = "complier",
                            y.var  = "depress3",
                            x.vars = model.forms$X.names)

        sens <- tryCatch(PIsens(data        = sim,
                                estimator   = "IFH",
                                Z.form      = model.forms$Z.form,
                                C.form      = model.forms$C.form,
                                Y.form      = model.forms$Y.form,
                                Y.model     = "quasibinomial",
                                Y.bounds    = Y.bounds,
                                sens.type   = "SMDe",
                                sens.range  = c(-1,1),
                                sens.step   = .05,
                                boot.seed   = 12345,
                                BCa.CI      = FALSE,
                                double.boot = TRUE)$sens,
                         error = function(e) return(NULL))

        rm(sim)

        .cleanMem(5)

        if (!is.null(sens)) {

            sens <- simplify2array(sens)

            raw.est[,,z]    <- sens[,"point",]
            boot.mean1[,,z] <- sens[,"boot.mean",]
            boot.mean2[,,z] <- sens[,"boot2.mean",]

            rm(sens)

            saveRDS(list(master      = sim.params$master,
                         sens.params = sens.params,
                         results     = list(raw = raw.est,
                                            bm1 = boot.mean1,
                                            bm2 = boot.mean2),
                         ok          = sum(!is.na(raw.est[1,1,]))),
                    file = file.path(sim.dir, "sim.rds"))

            .cleanMem(2)
        }
    }
}




rm(sens.params, true.values,
   raw.est, boot.mean1, boot.mean2,
   sim.params, sim.basis, z)

.cleanMem(5)



sim.results <- ..process_sim_results(sim.results  = readRDS(file.path(sim.dir, "sim.rds")),
                                     true.results = readRDS(file.path(sim.dir, "true.rds")))

saveRDS(sim.results, file = file.path(sim.dir, "sim.rds"))

rm(sim.results)

.cleanMem(5)









#

















