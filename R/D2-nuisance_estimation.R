################################################################################
#### .clean_nuis_inputs

#' (internal) Clean inputs used to specify nuisance functions to be estimated
#'
#' @param sensitivity \code{TRUE} (default) means do this for a sensitivity analysis, and \code{FAULT} means for a main analysis. (The only difference this makes is for the \eqn{e\mu} estimator where the main analysis but not the sensitivity analysis requires \code{C.form}.)
#' @inheritParams PIsens
#' @return A list of the specification inputs required for the combination of (\code{estimator}, \code{targeted} and \code{sensitivity}).
#' @keywords internal


.clean_nuis_inputs <- function(estimator,
                               targeted,
                               sensitivity = TRUE,
                               sens.type,
                               y,
                               Z.form = NULL,
                               C.form = NULL,
                               Y.form = NULL,
                               Y.model = NULL,
                               Y.bounds = NULL,
                               Y.FUN = NULL) {


    error1 <- "This estimator requires specifying Z.form."
    error2 <- "This estimator requires specifying C.form."
    error3 <- "All sensitivity analysis methods require specifying C.form."
    error4 <- "I cannot find Y.form."
    error5 <- "Y.form must be of length 1 or length 3."
    error6 <- "Y.model must be specified either as a character string (e.g., \"linear\", \"logit\", \"poisson\", \"gamma\") or as an object of class family. See details on options in function documentation."
    error7 <- "To fit a logit model to a general range outcome, the outcome's lower and upper bounds must be specified through Y.bounds."
    error8 <- "The outcome variable includes values outside of the specified Y.bounds."


    inputs <- list()


    ##########
    # Z.form

    if (targeted) {
        inputs[["Z.form"]] <- .replace_y(Z.form, new.y = "z")

    } else if (estimator != "pimu") {
        if (is.null(Z.form)) stop(error1)
        inputs[["Z.form"]] <- .replace_y(Z.form, new.y = "z")
    }


    ##########
    # C.form

    # for main analysis
    if (!sensitivity) {
        if (targeted) {
            inputs[["C.form"]] <- .replace_y(C.form, new.y = "c")

        } else if (estimator != "emu") {
            if (is.null(C.form)) stop(error2)
            inputs[["C.form"]] <- .replace_y(C.form, new.y = "c")
        }
    }

    # for sensitivity analysis
    if (sensitivity) {
        if (is.null(C.form)) stop(error3)
        inputs[["C.form"]] <- .replace_y(C.form, new.y = "c")
    }


    ##########
    # Y model inputs

    if (estimator!="epi") {

        # Y.FUN
        if (!is.null(Y.FUN)) {
            inputs[["Y.FUN"]] <- Y.FUN


        } else {

            # Y.form or Y.forms
            if (is.null(Y.form)) stop(error4)

            if (estimator=="emu") {
                inputs[["Y0.form"]] <- .replace_y(Y.form, new.y = "y")

            } else if (length(Y.form)==1) {
                inputs[["Y1c.form"]] <-
                    inputs[["Y1n.form"]] <-
                    inputs[["Y0.form"]]  <- .replace_y(Y.form, new.y = "y")

            } else if (length(Y.form)==3) {
                inputs[["Y1c.form"]] <- .replace_y(Y.form[1], new.y = "y")
                inputs[["Y1n.form"]] <- .replace_y(Y.form[2], new.y = "y")
                inputs[["Y0.form"]]  <- .replace_y(Y.form[3], new.y = "y")

            } else {
                stop(error5)
            }


            # Y.family
            if (is.null(Y.model)) {
                Y.family = stats::gaussian()

            } else if (class(Y.model)=="family") {

                Y.family <- Y.model

                if (Y.family$family=="binomial") {
                    Y.family$family <- "quasibinomial"

                } else if (Y.family$family=="poisson") {
                    Y.family$family <- "quasipoisson"
                }

            } else if (class(Y.model)=="character") {

                if (Y.model %in% c("linear", "gaussian", "identity")) {
                    Y.family = stats::gaussian()

                } else if (Y.model %in% c("binomial", "quasibinomial",
                                          "logit", "logistic", "odds ratio")) {
                    Y.family = stats::quasibinomial()

                } else if (Y.model %in% c("poisson", "quasipoisson",
                                          "count", "rate ratio")) {
                    Y.family = stats::quasipoisson()

                } else if (Y.model %in% c("Gamma", "gamma",
                                          "positive outcome")) {
                    Y.family = stats::Gamma(link = "log")
                }
            } else {
                stop(error6)
            }

            inputs[["Y.family"]] <- Y.family


            # Y.bounds
            if (Y.family$family=="quasibinomial" && !is.within.bounds(y, 0:1)) {

                if (is.null(Y.bounds)) stop(error7)

                if (!is.within.bounds(y, Y.bounds)) stop(error8)

                inputs[["Y.bounds"]] <- Y.bounds
            }


            # get.sigma
            if (sens.type=="SMDe")  inputs[["get.sigma"]] <- TRUE

        }

    }








    ##########
    # TODO: This last part needs to be re-considered after paper 2.

    # if (estimator != "epi") {
    #
    #     if (sensitivity && sens.type=="beta.quant") {
    #
    #         if (!is.null(muquantFUN)) { inputs[["Y.FUN"]] <- muquantFUN
    #         } else                    { stop("muquantFUN not found")
    #         }
    #
    #     } else {
    #
    #         if (!is.null(Y.FUN)) { inputs[["Y.FUN"]] <- Y.FUN
    #
    #         } else {
    #             if (is.null(Y0.form)) { stop("Y0.form not found.")
    #             } else { inputs[["Y0.form"]] <- .replace_y(Y0.form, new.y = "y")
    #             }
    #
    #             if (is.null(Y.family)) { stop("Y.family not found.")
    #             } else                 { inputs[["Y.family"]] <- Y.family
    #             }
    #
    #             if (estimator != "emu") {
    #                 if (is.null(Y1c.form)) { stop("Y1c.form not found.")
    #                 } else { inputs[["Y1c.form"]] <- .replace_y(Y1c.form, new.y = "y")
    #                 }
    #                 if (is.null(Y1n.form)) { stop("Y1n.form not found.")
    #                 } else { inputs[["Y1n.form"]] <- .replace_y(Y1n.form, new.y = "y")
    #                 }
    #             }
    #         }
    #     }
    #
    #
    # }

    inputs

}




################################################################################
#### .estimate_nuisance

#' (internal) Estimate nuisance functions needed for effect estimation
#'
#' @inheritParams PIsens
#' @param Y1c.form Formula for Y|X,Z=1,C=1 model, processed by \code{.clean_nuis_inputs()}.
#' @param Y1n.form Formula for Y|X,Z=1,C=0 model, processed by \code{.clean_nuis_inputs()}.
#' @param Y0.form Formula for Y|X,Z=0 model, processed by \code{.clean_nuis_inputs()}.
#' @param Y.family The family to be used for outcome models, processed by \code{.clean_nuis_inputs()}.
#' @importFrom stats glm predict
#'
#' @returns A dataframe with variables \code{s.wt}; grouping variables \code{z} (treatment vs control arm), \code{c} (=1 if complier in treatment arm, 0 otherwise), \code{n} (=1 if noncomplier in treatment arm, 0 otherwise); outcome \code{y}; plus variables holding the requested nuisance functions: \code{z.wt} (for inverse probability of assigned treatment weight), \code{c.wt} and \code{n.wt} (for the two principal score weights), \code{mu.1c}, \code{mu.1n}, \code{mu.0c}, \code{mu.0n} (for the four outcome mean functions).
#'
#' @keywords internal
#'
#' TODO: edit returns

.estimate_nuisance <- function(data,
                               estimator,
                               targeted,
                               get.sigma = FALSE,
                               Z.form    = NULL,
                               C.form    = NULL,
                               Y1c.form  = NULL,
                               Y1n.form  = NULL,
                               Y0.form   = NULL,
                               Y.family  = NULL,
                               Y.bounds  = NULL,
                               Y.FUN     = NULL) {


    dat <- data; rm(data)

    dat$c <- ifelse(dat$z==0, 0, dat$c)
    dat$n <- ifelse(dat$z==0, 0, 1-dat$c)

    nuis <- dat[c("id", "s.wt", "z", "c", "n", "y")]



    ########################################
    # Z weights

    if (!is.null(Z.form)) {

        suppressWarnings({
            z.mod <- do.call("glm",
                             list(formula = Z.form,
                                  data    = dat,
                                  weights = dat$s.wt,
                                  family  = "quasibinomial"))
        })

        if(!z.mod$converged) stop("glm error")

        nuis$z.wt <-
            dat$z / z.mod$fitted.values + (1-dat$z) / (1-z.mod$fitted.values)


        if (targeted) dat$z.wt <- nuis$z.wt
    }



    ########################################
    # C weights

    if (!is.null(C.form)) {

        if (targeted) {
            suppressWarnings({
                c.mod <- do.call("glm",
                                 list(formula = C.form,
                                      data    = dat[dat$z==1, ],
                                      weights = (dat$s.wt*dat$z.wt)[dat$z==1],
                                      family  = "quasibinomial"))
            })

        } else {
            suppressWarnings({
                c.mod <- do.call("glm",
                                 list(formula = C.form,
                                      data    = dat[dat$z==1, ],
                                      weights = dat$s.wt[dat$z==1],
                                      family  = "quasibinomial"))
            })
        }

        nuis$c.wt <- predict(c.mod, newdata = dat, type = "response")
        nuis$n.wt <- 1 - nuis$c.wt


        if (targeted) {
            dat$c.wt <- nuis$c.wt
            dat$n.wt <- nuis$n.wt
        }
    }



    ########################################
    # weights estimation only -- called by estimate_weights()

    if (missing(estimator)) return(nuis)




    ########################################
    # Y means version 2: using Y.FUN

    if (!is.null(Y.FUN)) {

        # TODO: make an external function that gives the user this dataset dat2
        # to use in building their Y.FUN

        dat2 <- dat
        dat2$group <- ifelse(dat$z==0, "0",
                             ifelse(dat$c==1, "1c", "1n"))

        if (targeted)
            dat2 <- .df.rename(dat2,
                               from = c("z.wt", "c.wt", "n.wt"),
                               to   = c("z.ipw", "c.ps", "n.ps"))


        out <- do.call(Y.FUN, args = c(data = list(dat2)))


        if (is.data.frame(out))
            return(merge(nuis, out))

        out$nuis <- merge(nuis, out$nuis)
        return(out)
    }





    ########################################
    # Y means etc. version 1: fitting models

    if (!is.null(Y.bounds))
        dat$y <- .trans_bounds(dat$y, from = Y.bounds, to = 0:1)


    #########################
    # Y means for Z=1

    if (!is.null(Y1c.form)) {

        if (targeted) {
            Y1c.wt <- (dat$s.wt*dat$z.wt)[dat$c==1]
            Y1n.wt <- (dat$s.wt*dat$z.wt)[dat$n==1]
        } else {
            Y1c.wt <- dat$s.wt[dat$c==1]
            Y1n.wt <- dat$s.wt[dat$n==1]
        }

        nuis$mu.1c <- .get_mu_hat(form     = Y1c.form,
                                  data     = dat[dat$c==1, ],
                                  wts      = Y1c.wt,
                                  y.fam    = Y.family,
                                  pred.dat = dat)$mu

        nuis$mu.1n <- .get_mu_hat(form     = Y1n.form,
                                  data     = dat[dat$n==1, ],
                                  wts      = Y1n.wt,
                                  y.fam    = Y.family,
                                  pred.dat = dat)$mu

        if (!is.null(Y.bounds)) {
            nuis$mu.1c <- .trans_bounds(nuis$mu.1c, from = 0:1, to = Y.bounds)
            nuis$mu.1n <- .trans_bounds(nuis$mu.1n, from = 0:1, to = Y.bounds)
        }
    }


    #########################
    # Y means etc. for Z=0

    if (!is.null(Y0.form)) {

        if (targeted) {
            mu.0c <- .get_mu_hat(form           = Y0.form,
                                 data           = dat[dat$z==0, ],
                                 wts            = (dat$s.wt*dat$z.wt*dat$c.wt)[dat$z==0],
                                 y.fam          = Y.family,
                                 pred.dat       = dat,
                                 get.dispersion = get.sigma)

            mu.0n <- .get_mu_hat(form           = Y0.form,
                                 data           = dat[dat$z==0, ],
                                 wts            = (dat$s.wt*dat$z.wt*dat$n.wt)[dat$z==0],
                                 y.fam          = Y.family,
                                 pred.dat       = dat,
                                 get.dispersion = get.sigma)

            nuis$mu.0c <- mu.0c$mu
            nuis$mu.0n <- mu.0n$mu

            if (get.sigma) {
                dispersion <-
                    mean(nuis$c.wt) * mu.0c$dispersion +
                    mean(nuis$n.wt) * mu.0n$dispersion

                nuis$sigma.0c <- sqrt(dispersion * Y.family$variance(nuis$mu.0c))
                nuis$sigma.0n <- sqrt(dispersion * Y.family$variance(nuis$mu.0n))
            }

        } else {
            mu.0 <- .get_mu_hat(form           = Y0.form,
                                data           = dat[dat$z==0, ],
                                wts            = dat$s.wt[dat$z==0],
                                y.fam          = Y.family,
                                pred.dat       = dat,
                                get.dispersion = get.sigma)

            nuis$mu.0n <- nuis$mu.0c <- mu.0$mu

            if (get.sigma) {
                nuis$sigma.0c <- nuis$sigma.0n <-
                    sqrt(mu.0$dispersion * Y.family$variance(mu.0$mu))
            }

        }


        if (!is.null(Y.bounds)) {
            nuis$mu.0c <- .trans_bounds(nuis$mu.0c, from = 0:1, to = Y.bounds)
            nuis$mu.0n <- .trans_bounds(nuis$mu.0n, from = 0:1, to = Y.bounds)

            if (get.sigma) {
                nuis$sigma.0c <- nuis$sigma.0c * abs(diff(Y.bounds))
                nuis$sigma.0n <- nuis$sigma.0n * abs(diff(Y.bounds))
            }

        }

    }









    nuis
}



