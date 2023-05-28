#############################
#### .assemble_GOR

#' (internal) Assemble OR-based sensitivity analysis for a specified estimator based on main analysis estimates and estimated nuisance functions.
#'
#' @inheritParams .assemble_sens
#' @return A matrix with one row (named \code{rho}) holding the values of the sensitivity parameter and two rows for the corresponding CATE and NACE point estimates.
#'
#' @keywords internal

.assemble_GOR <- function(estimator,
                          main,
                          nuis,
                          sens.range,
                          sens.step,
                          Y.bounds = NULL) {



    # the three type 2 estimators

    if (estimator %in% c("pimu", "emu", "MS")) {

        if (estimator=="pimu") {
            cace.wt <- nuis$c.wt
            nace.wt <- nuis$n.wt
        } else if (estimator=="emu") {
            cace.wt <- nuis$c * nuis$z.wt
            nace.wt <- nuis$n * nuis$z.wt
        } else if (estimator=="MS") {
            cace.wt <- nuis$z*nuis$z.wt*(nuis$c - nuis$c.wt) + nuis$c.wt
            nace.wt <- nuis$z*nuis$z.wt*(nuis$n - nuis$n.wt) + nuis$n.wt
        }

        cace.wt <- nuis$s.wt * cace.wt
        nace.wt <- nuis$s.wt * nace.wt


        sens <-
            t(sapply(.get_rho.vec(sens.range,
                                  sens.step),
                     function(r) {
                         if (r==1) {
                             return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                         } else {

                             gor <- .compute_mu_epsilon_GOR(nuis        = nuis,
                                                            r           = r,
                                                            Y.bounds     = Y.bounds,
                                                            get.epsilon = FALSE)

                             tau.0c <- .wtd_mean(x = gor$mu.0c, w = cace.wt)
                             tau.0n <- .wtd_mean(x = gor$mu.0n, w = nace.wt)

                             cace <- main[["tau.1c"]] - tau.0c
                             nace <- main[["tau.1n"]] - tau.0n

                             return(c(cace   = cace,
                                      nace   = nace,
                                      tau.0c = tau.0c,
                                      tau.0n = tau.0n))
                         }
                     }))
    }


    # the two type 3 estimators

    if (estimator %in% c("IF", "IFH")) {

        z <- nuis$z
        c <- nuis$c
        n <- nuis$n
        y <- nuis$y

        z.wt <- nuis$z.wt
        c.wt <- nuis$c.wt
        n.wt <- nuis$n.wt
        mu.0c <- nuis$mu.0c
        mu.0n <- nuis$mu.0n

        if (estimator=="IFH")
            z.wt <- z.wt / (z *.wtd_mean(z*z.wt, nuis$s.wt) +
                                (1-z) *.wtd_mean((1-z)*z.wt, nuis$s.wt))

        phi_delta.c <- z*z.wt*(c-c.wt) + c.wt

        delta.c <- .wtd_mean(x = phi_delta.c, w = nuis$s.wt)
        delta.n <- 1 - delta.c


        sens <-
            t(sapply(.get_rho.vec(sens.range,
                                  sens.step),
                     function(r) {

                         if (r==1) {
                             return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                         } else {

                             gor <- .compute_mu_epsilon_GOR(nuis        = nuis,
                                                            r           = r,
                                                            Y.bounds    = Y.bounds,
                                                            get.epsilon = TRUE)

                             term1.c <- (1-z)*z.wt * gor$epsilon.c.mu * (y-mu.0c)
                             term1.n <- (1-z)*z.wt * gor$epsilon.n.mu * (y-mu.0n)

                             term2.c <- z*z.wt * gor$epsilon.c.pi * (c-c.wt)
                             term2.n <- z*z.wt * gor$epsilon.n.pi * (n-n.wt)

                             term3.c <- c.wt * gor$mu.0c
                             term3.n <- n.wt * gor$mu.0n

                             phi_nu.0c <- term1.c + term2.c + term3.c
                             phi_nu.0n <- term1.n + term2.n + term3.n

                             nu.0c <- .wtd_mean(x = phi_nu.0c, w = nuis$s.wt)
                             nu.0n <- .wtd_mean(x = phi_nu.0n, w = nuis$s.wt)

                             tau.0c <- nu.0c / delta.c
                             tau.0n <- nu.0n / delta.n

                             cace <- main[["tau.1c"]] - tau.0c
                             nace <- main[["tau.1n"]] - tau.0n

                             return(c(cace   = cace,
                                      nace   = nace,
                                      tau.0c = tau.0c,
                                      tau.0n = tau.0n))
                         }
                     }))

    }



    sens <- cbind(sens,
                  tau.1c = main[["tau.1c"]],
                  tau.1n = main[["tau.1n"]])

    sens
}




#############################
#### .compute_mu_epsilon_GOR

#' (internal) Compute OR-based mu.0c and mu.0n functions for a value of the sentivity parameter.
#'
#' This function is called by \code{.assemble_OR_mu()}.
#'
#' @inheritParams .assemble_GOR
#' @param r A value of the sensitivity parameter
#' @keywords internal

.compute_mu_epsilon_GOR <- function(nuis,
                                    r,
                                    Y.bounds = NULL,
                                    get.epsilon = FALSE) {

    c.wt <- nuis$c.wt
    n.wt <- nuis$n.wt
    mu.0c <- nuis$mu.0c
    mu.0n <- nuis$mu.0n

    r.c <- r
    r.n <- 1/r

    if (!is.null(Y.bounds)) {
        mu.0c <- .trans_bounds(mu.0c, from = Y.bounds, to = 0:1)
        mu.0n <- .trans_bounds(mu.0n, from = Y.bounds, to = 0:1)
    }



    alpha.c <- (c.wt + mu.0c) * (r.c-1) + 1
    alpha.n <- (n.wt + mu.0n) * (r.n-1) + 1

    beta.c <- sqrt(alpha.c^2 - 4*c.wt*mu.0c*r.c*(r.c-1))
    beta.n <- sqrt(alpha.n^2 - 4*n.wt*mu.0n*r.n*(r.n-1))

    mu.0c.GOR <- ifelse(c.wt==0, 0, (alpha.c - beta.c) / (2*(r.c-1) * c.wt))
    mu.0n.GOR <- ifelse(n.wt==0, 0, (alpha.n - beta.n) / (2*(r.n-1) * n.wt))

    if (!is.null(Y.bounds)) {
        mu.0c.GOR <- .trans_bounds(mu.0c.GOR, from = 0:1, to = Y.bounds)
        mu.0n.GOR <- .trans_bounds(mu.0n.GOR, from = 0:1, to = Y.bounds)
    }

    out <- list(mu.0c = mu.0c.GOR,
                mu.0n = mu.0n.GOR)

    if (get.epsilon) {
        epsilon.c.mu <- .5 - alpha.c/(2*beta.c) + r.c*c.wt/beta.c
        epsilon.n.mu <- .5 - alpha.n/(2*beta.n) + r.n*n.wt/beta.n

        epsilon.c.pi <- .5 - alpha.c/(2*beta.c) + r.c*mu.0c/beta.c
        epsilon.n.pi <- .5 - alpha.n/(2*beta.n) + r.n*mu.0n/beta.n

        if (!is.null(Y.bounds)) {
            epsilon.c.pi <- .trans_bounds(epsilon.c.pi, from = 0:1, to = Y.bounds)
            epsilon.n.pi <- .trans_bounds(epsilon.n.pi, from = 0:1, to = Y.bounds)
        }

        out <- c(out, mget(c("epsilon.c.mu",
                             "epsilon.n.mu",
                             "epsilon.c.pi",
                             "epsilon.n.pi")))
    }


    out
}





































