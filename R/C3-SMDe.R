#############################
#### .assemble_SMDe

#' (internal) Assemble SMDe-based sensitivity analysis for a specified estimator based on main analysis estimates and estimated nuisance functions.
#'
#' @inheritParams .assemble_sens
#' @return A matrix with one row (named \code{rho}) holding the values of the sensitivity parameter and two rows for the corresponding CATE and NACE point estimates.
#'
#' @keywords internal

.assemble_SMDe <- function(estimator,
                           main,
                           nuis,
                           sens.range,
                           sens.step) {



    # the three type 2 estimators

    if (estimator %in% c("pimu", "emu")) {

        if (estimator=="pimu") {
            cace.wt <- nuis$c.wt
            nace.wt <- nuis$n.wt
        } else if (estimator=="emu") {
            cace.wt <- nuis$c * nuis$z.wt
            nace.wt <- nuis$n * nuis$z.wt
        }

        cace.wt <- nuis$s.wt * cace.wt
        nace.wt <- nuis$s.wt * nace.wt


        sens <-
            t(sapply(.get_eta.vec(sens.range,
                                  sens.step),
                     function(h) {

                         if (h==0) {
                             return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                         } else {

                             smde <- .compute_lambda_upsilon_SMDe(nuis        = nuis,
                                                                  h           = h,
                                                                  get.upsilon = FALSE)

                             tau.0c <- main[["tau.0c"]] + h * .wtd_mean(x = smde$lambda.c, w = cace.wt)
                             tau.0n <- main[["tau.0n"]] - h * .wtd_mean(x = smde$lambda.n, w = nace.wt)

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

    if (estimator %in% c("IF", "IFH", "MS")) {

        z <- nuis$z
        c <- nuis$c
        n <- nuis$n
        y <- nuis$y

        z.wt <- nuis$z.wt
        c.wt <- nuis$c.wt
        n.wt <- nuis$n.wt
        mu.0c <- nuis$mu.0c
        mu.0n <- nuis$mu.0n
        sigma.0c <- nuis$sigma.0c
        sigma.0n <- nuis$sigma.0n

        if (estimator=="IFH")
            z.wt <- z.wt / (z *.wtd_mean(z*z.wt, nuis$s.wt) +
                                (1-z) *.wtd_mean((1-z)*z.wt, nuis$s.wt))

        phi_delta.c <- z*z.wt*(c-c.wt) + c.wt

        delta.c <- .wtd_mean(x = phi_delta.c, w = nuis$s.wt)
        delta.n <- 1 - delta.c


        sens <-
            t(sapply(.get_eta.vec(sens.range,
                                  sens.step),
                     function(h) {

                         if (h==0) {
                             return(main[c("cace", "nace", "tau.0c", "tau.0n")])

                         } else {

                             smde <- .compute_lambda_upsilon_SMDe(nuis        = nuis,
                                                                  h           = h,
                                                                  get.upsilon = TRUE)

                             term1.c <- (1-z)*z.wt * (smde$upsilon /(2*sigma.0c)) * ((y-mu.0c)^2 - sigma.0c^2)
                             term1.n <- (1-z)*z.wt * (smde$upsilon /(2*sigma.0n)) * ((y-mu.0n)^2 - sigma.0n^2)

                             term2.c <- z*z.wt * smde$upsilon.dot * sigma.0c
                             term2.n <- z*z.wt * smde$upsilon.dot * sigma.0n

                             term3.c <- smde$upsilon * sigma.0c
                             term3.n <- smde$upsilon * sigma.0n

                             phi_vtheta.0c <- term1.c + term2.c + term3.c
                             phi_vtheta.0n <- term1.n + term2.n + term3.n

                             vtheta.0c <- .wtd_mean(x = phi_vtheta.0c, w = nuis$s.wt)
                             vtheta.0n <- .wtd_mean(x = phi_vtheta.0n, w = nuis$s.wt)

                             tau.0c <- main[["tau.0c"]] + h * vtheta.0c / delta.c
                             tau.0n <- main[["tau.0n"]] - h * vtheta.0n / delta.n

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
#### .compute_lambda_upsilon_SMDe

#' (internal) Compute SMDe-based lambda or upsilon functions for a value of the sentivity parameter.
#'
#' This function is called by \code{.assemble_SMDe()}.
#'
#' @inheritParams .assemble_SMDe
#' @param h A value of the sensitivity parameter
#' @keywords internal

.compute_lambda_upsilon_SMDe <- function(nuis,
                                         h,
                                         get.upsilon = FALSE) {

    c.wt <- nuis$c.wt
    n.wt <- nuis$n.wt


    denom <- 1 + h^2 * c.wt * n.wt

    if (get.upsilon)
        return(list(upsilon     = c.wt * n.wt / sqrt(denom),
                    upsilon.dot = (n.wt - c.wt) * (nuis$c - c.wt) *
                                  (denom + 1) / (2 * (denom^1.5))))

    return(list(lambda.c = (n.wt / sqrt(denom)) * nuis$sigma.0c,
                lambda.n = (c.wt / sqrt(denom)) * nuis$sigma.0n))

}





































