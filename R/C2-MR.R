#############################
#### .assemble_MR

#' (internal) Assemble MR-based sensitivity analysis for a specified estimator based on main analysis estimates and estimated nuisance functions.
#'
#' @inheritParams .assemble_sens
#' @return A matrix with one row (named \code{rho}) holding the values of the sensitivity parameter and two rows for the corresponding CATE and NACE point estimates.
#'
#' @keywords internal

.assemble_MR <- function(estimator,
                         main,
                         nuis,
                         sens.range,
                         sens.step,
                         include.tau.1c = TRUE) {



    ########################################
    # Type A

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
                                  sens.step), function(r) {

                gamma.n <- 1 / ((r-1) * nuis$c.wt + 1)
                gamma.c <- gamma.n * r

                if (r==1) {
                    return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                } else {

                    tau.0c <- .wtd_mean(x = nuis$mu.0c * gamma.c, w = cace.wt)
                    tau.0n <- .wtd_mean(x = nuis$mu.0n * gamma.n, w = nace.wt)

                    cace <- main[["tau.1c"]] - tau.0c
                    nace <- main[["tau.1n"]] - tau.0n

                    return(c(cace   = cace,
                             nace   = nace,
                             tau.0c = tau.0c,
                             tau.0n = tau.0n))
                }
            }))
    }



    ########################################
    # Type B

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
                                  sens.step), function(r) {

                if (r==1) {
                    return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                } else {

                    gamma.n <- 1 / ((r-1) * nuis$c.wt + 1)
                    gamma.c <- gamma.n * r

                    tmp1.c <- (1-z)*z.wt*c.wt*gamma.c*(y-mu.0c)
                    tmp1.n <- (1-z)*z.wt*n.wt*gamma.n*(y-mu.0n)

                    tmp2.c <- z*z.wt*mu.0c*gamma.c*gamma.n*(c-c.wt)
                    tmp2.n <- z*z.wt*mu.0n*gamma.c*gamma.n*(n-n.wt)

                    tmp3.c <- c.wt*gamma.c*mu.0c
                    tmp3.n <- n.wt*gamma.n*mu.0n

                    phi_nu.0c <- tmp1.c + tmp2.c + tmp3.c
                    phi_nu.0n <- tmp1.n + tmp2.n + tmp3.n

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


    ########################################
    # Type C

    if (estimator=="epi") {

        tau.0c.wt <- nuis$s.wt* ((1-nuis$z)*nuis$z.wt * nuis$c.wt)
        tau.0n.wt <- nuis$s.wt* ((1-nuis$z)*nuis$z.wt * nuis$n.wt)


        sens <- t(sapply(.get_rho.vec(sens.range,
                                      sens.step), function(r) {

                                          gamma.n <- 1 / ((r-1) * nuis$c.wt + 1)
                                          gamma.c <- gamma.n * r

                                          if (r==1) {
                                              return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                                          } else {

                                              tau.0c <- .wtd_mean(x = nuis$y * gamma.c, w = tau.0c.wt)
                                              tau.0n <- .wtd_mean(x = nuis$y * gamma.n, w = tau.0n.wt)

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


