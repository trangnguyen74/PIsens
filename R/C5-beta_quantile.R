#############################
#### .assemble_beta.quant

#' (internal) Assemble OR-based sensitivity analysis for a specified estimator based on main analysis estimates and estimated nuisance functions.
#'
#' @inheritParams .assemble_sens
#' @return A matrix with one row (named \code{rho}) holding the values of the sensitivity parameter and two rows for the corresponding CATE and NACE point estimates.
#'
#' @keywords internal



.assemble_beta.quant <- function(estimator,
                                 main,
                                 nuis,
                                 quants,
                                 sens.range = c(-.1, .1),
                                 sens.step) {


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
            t(sapply(.get_signed.kappa.vec(sens.range,
                                           sens.step), function(s.kappa) {

                if (s.kappa==0) {
                    return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                } else {

                    tmp <- .compute_mu.0c_beta.quant(nuis = nuis,
                                                     quants = quants,
                                                     s.kappa = s.kappa)

                    tau.0c <- .wtd_mean(x = tmp$mu.0c, w = cace.wt)
                    tau.0n <- .wtd_mean(x = tmp$mu.0n, w = nace.wt)

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

        if (estimator=="IFH") {
            z.wt <-
                (z==1)*z.wt/.wtd_mean(z    *z.wt, nuis$s.wt) +
                (z==0)*z.wt/.wtd_mean((1-z)*z.wt, nuis$s.wt)
        }

        phi_delta.c <- z*z.wt*(c-c.wt) + c.wt

        delta.c <- .wtd_mean(x = phi_delta.c, w = nuis$s.wt)
        delta.n <- 1 - delta.c


        sens <-
            t(sapply(.get_signed.kappa.vec(sens.range,
                                           sens.step), function(s.kappa) {

                if (s.kappa==0) {
                    return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                } else {

                    tmp <- .compute_mu.0c_beta.quant(nuis = nuis,
                                                     quants = quants,
                                                     s.kappa = s.kappa)

                    phi_nu.0c <- (1-z)*z.wt*c.wt*(y-mu.0c)
                    phi_nu.0n <- (1-z)*z.wt*n.wt*(y-mu.0n)

                    phi_nu.0c <- phi_nu.0c + (z*z.wt*(c-c.wt) + c.wt) *tmp$mu.0c
                    phi_nu.0n <- phi_nu.0n + (z*z.wt*(n-n.wt) + n.wt) *tmp$mu.0n

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

    sens
}




#############################
#### .compute_mu.0c_beta.quant

#' (internal) Compute beta-based mu.0c and mu.0n functions for a value of the sentivity parameter.
#'
#' This function is called by \code{.assemble_beta.quant()}.
#'
#' @inheritParams .assemble_beta.quant
#' @param s.kappa A value of the sensitivity parameter signed-\eqn{\kappa}.
#' @keywords internal

.compute_mu.0c_beta.quant <- function(nuis, quants, s.kappa) {

    m <- ncol(quants[[1]])

    c.points <- (1:m / m) - (1 / (2*m))

    if (s.kappa < 0) c.points <- 1 - c.points

    k <- abs(s.kappa)


    mus <- t(sapply(1:nrow(nuis), function(z) {

        c.wt <- nuis$c.wt[z]
        n.wt <- 1 - c.wt

        alpha <- c.wt * (1-k)/k
        beta  <- n.wt * (1-k)/k

        if (alpha <= beta) {
            c.wt.quants <- stats::qbeta(p = c.points, shape1 = alpha, shape2 = beta)
        } else {
            c.wt.quants <- 1 - stats::qbeta(p = 1-c.points, shape1 = beta, shape2 = alpha)
        }


        c.wt.mean <- mean(c.wt.quants)

        if (c.wt.mean==0 || c.wt.mean==1) {

            out <- c(mu.0c = nuis$mu.0c[z],
                     mu.0n = nuis$mu.0n[z])

        } else {

            # calibrate quantiles to c.wt and n.wt
            if (c.wt.mean >= c.wt) {
                c.wt.quants <- c.wt.quants * c.wt / c.wt.mean
                n.wt.quants <- 1 - c.wt.quants
            } else {
                n.wt.quants <- (1-c.wt.quants) * n.wt / (1-c.wt.mean)
                c.wt.quants <- 1 - n.wt.quants
            }

            out <- c(mu.0c = .wtd_mean(quants$y0c[z, ], c.wt.quants),
                     mu.0n = .wtd_mean(quants$y0n[z, ], n.wt.quants))

        }

        out
    }))

    data.frame(mus)

}
