####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### PImain ####

#' Main analysis assuming principal ignorability (PI)
#'
#' @inheritParams PIsens
#' @export


PImain <- function(
        data,
        estimator,
        targeted,

        Z.form,
        C.form,

        Y.form,
        Y.model,
        Y.bounds,
        Y.FUN = NULL,

        boot.num = 999,
        boot.seed = 12345,
        boot.cont.wt = TRUE,
        double.boot = TRUE,
        BCa.CI = TRUE) {


    msg1 <- "data needs to go through declare_data() before use with PImain() or PIsens()."
    msg2 <- "estimator must be specified."
    msg3 <- "Untargeted nuisance estimation is not implemented for the MS estimator. Proceeding with targeted = TRUE."
    msg4 <- "For targeted nuisance estimation, must specify both Z.form and C.form."



    ########################################
    # Check if data are declared data

    if (!.is.declared(data)) stop(msg1)



    ########################################
    # Assemble key inputs: estimator and targeted

    if (missing(estimator)) stop(msg2)

    if (estimator=="MS" && targeted==FALSE) {
        message(msg3)
        targeted <- TRUE
    }



    # default to targeted unless (i) Z.form and C.form aren't both provided or (ii) targeted==FALSE
    if (missing(Z.form) || missing(C.form)) {
        if (!missing(targeted) && targeted==TRUE) stop(msg4)
        targeted <- FALSE

    } else if (!missing(targeted) && targeted==FALSE) {
        targeted <- FALSE

    } else {
        targeted <- TRUE
    }

    key.inputs <- list(estimator = estimator,
                       targeted  = targeted)



    ########################################
    # Assemble nuis inputs: Z, C, Y inputs

    nuis.inputs <- mget(c("Z.form", "C.form",
                          "Y.form", "Y.model", "Y.bounds",
                          "Y.FUN"),
                        ifnotfound = list(NULL, NULL,
                                          NULL, NULL, NULL,
                                          NULL))
    nuis.inputs <- c(y = list(data$y),
                     nuis.inputs)

    nuis.inputs <- do.call(".clean_nuis_inputs",
                           args = c(key.inputs,
                                    sensitivity = FALSE,
                                    nuis.inputs))



    ########################################
    # Assemble boot inputs

    boot.inputs <- mget(c("boot.num", "boot.seed", "boot.cont.wt",
                          "double.boot", "BCa.CI"))



    ########################################
    # Construct function do_one() to get point estimates from one dataset

    do_one <- function(data,
                       .key.inputs  = key.inputs,
                       .nuis.inputs = nuis.inputs) {

        nuis <- do.call(".estimate_nuisance",
                        args = c(data = list(data),
                                 .key.inputs,
                                 .nuis.inputs))

        do.call(".assemble_main",
                args = c(estimator = .key.inputs$estimator,
                         nuis      = list(nuis)))

    }


    ########################################
    # Run function do_one() through .do_all_main() to bootstrap and all

    do.call(".do_all_main",
            args = c(pointFUN = do_one,
                     data     = list(data),
                     boot.inputs))
}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .do_all_main ####

#' (internal) Compute confidence intervals and bias-correct point estimates for a main analysis
#'
#' @param pointFUN A function constructed to compute the point estimate, with all options wrapped in, so the only argument of this function is \code{data}. This function should output a numeric vector with two estimates, named "cace" and "nace".
#' @inheritParams PImain
#' @return A table with two rows for CACE and NACE, respectively. The table contains columns for the point estimate, the bias-corrected point estimate, the percentile CI, the BCa CI, and also the bootstrap and double bootstrap means of the estimates.
#' @keywords internal


.do_all_main <- function(pointFUN,
                         data,
                         boot.num = 999,
                         boot.seed = 12345,
                         boot.cont.wt = TRUE,
                         double.boot = TRUE,
                         BCa.CI = TRUE) {

    # get point estimates
    point <- pointFUN(data)


    successful.boot.num <- list()


    # estimate acceleration parameters for BCa confidence interval estimation
    if (BCa.CI)
        ahat <- .get_ahat_for_BCa(FUN = pointFUN, data = data)


    gc()


    # get boot estimates
    boot1 <- .get_boot_estimates(FUN          = pointFUN,
                                 data         = data,
                                 boot.num     = boot.num,
                                 boot.cont.wt = boot.cont.wt,
                                 boot.seed    = boot.seed)

    successful.boot.num[["boot1"]] <- dim(boot1)[2]

    # bias-correct point estimates
    boot1.mean  <- rowMeans(boot1, na.rm = TRUE)
    point.bc1 <- 2*(point - boot1.mean) + boot1.mean



    # bootstrap confidence intervals: BCa and percentile types
    c.intervals <- t(sapply(1:length(point), function(z) {

        stats::quantile(boot1[z,], prob = c(0.025, 0.975), na.rm = TRUE)
    }))

    colnames(c.intervals) <- c("perc.low", "perc.high")


    if (BCa.CI) {

        ci.BCa <- t(sapply(1:length(point), function(z) {

            .get_BCa.interval(theta.boot = boot1[z,],
                              theta.hat  = point[z],
                              ahat       = ahat[z])
        }))
        colnames(ci.BCa) <- c("bca.low", "bca.high")

        c.intervals <- cbind(ci.BCa, c.intervals)
    }

    rm(boot1, ahat)

    gc()



    if (double.boot) {
        # single-draw double-bootstrap estimates
        boot2 <- .get_boot_estimates(FUN          = pointFUN,
                                     data         = data,
                                     boot.num     = boot.num,
                                     boot.cont.wt = boot.cont.wt,
                                     boot.seed    = boot.seed,
                                     double.boot  = TRUE)

        successful.boot.num[["boot2"]] <- dim(boot2)[2]


        # bias-correct point estimates
        boot2.mean <- rowMeans(boot2, na.rm = TRUE)
        point.bc2 <- 3*(point - boot1.mean) + boot2.mean

        rm(boot2)

        gc()
    }



    if (double.boot) {

        main <- cbind(point      = point,
                      point.bc1  = point.bc1,
                      point.bc2  = point.bc2,
                      c.intervals,
                      boot.mean  = boot1.mean,
                      boot2.mean = boot2.mean)
    } else {

        main <- cbind(point      = point,
                      point.bc1  = point.bc1,
                      c.intervals,
                      boot.mean  = boot1.mean)
    }


    mget(c("main", "successful.boot.num"))


}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .assemble_main ####

#' (internal) Assemble PI-based main analysis for a specified estimator based on estimated nuisance functions.
#'
#' @inheritParams PImain
#' @param nuis A dataframe with variables \code{s.wt}; grouping variables \code{z} (treatment vs control arm), \code{c} (=1 if complier in treatment arm, 0 otherwise), \code{n} (=1 if noncomplier in treatment arm, 0 otherwise); outcome \code{y}; plus variables holding the requested nuisance functions: \code{z.wt} (for inverse probability of assigned treatment weight), \code{c.wt} and \code{n.wt} (for the two principal score weights), \code{mu.1c}, \code{mu.1n}, \code{mu.0c}, \code{mu.0n} (for the four outcome mean functions).
#' @return A vector of estimates with elements named \code{cace}, \code{nace}, \code{tau.1c}, \code{tau.1n}, \code{tau.0c}, \code{tau.0n}.
#' @keywords internal

.assemble_main <- function(estimator,
                           nuis) {


    if (estimator=="epi") {

        tau.1c <- .wtd_mean(x = nuis$y, w = nuis$s.wt* nuis$c*nuis$z.wt)
        tau.1n <- .wtd_mean(x = nuis$y, w = nuis$s.wt* nuis$n*nuis$z.wt)
        tau.0c <- .wtd_mean(x = nuis$y, w = nuis$s.wt* (1-nuis$z)*nuis$z.wt * nuis$c.wt)
        tau.0n <- .wtd_mean(x = nuis$y, w = nuis$s.wt* (1-nuis$z)*nuis$z.wt * nuis$n.wt)
    }



    if (estimator=="pimu") {

        cace.wt <- nuis$s.wt * nuis$c.wt
        nace.wt <- nuis$s.wt * nuis$n.wt

        tau.1c <- .wtd_mean(x = nuis$mu.1c, w = cace.wt)
        tau.0c <- .wtd_mean(x = nuis$mu.0c, w = cace.wt)

        tau.1n <- .wtd_mean(x = nuis$mu.1n, w = nace.wt)
        tau.0n <- .wtd_mean(x = nuis$mu.0n, w = nace.wt)
    }



    if (estimator=="emu") {

        cace.wt <- nuis$s.wt* (nuis$c * nuis$z.wt)
        nace.wt <- nuis$s.wt* (nuis$n * nuis$z.wt)

        tau.1c <- .wtd_mean(x = nuis$y,     w = cace.wt)
        tau.0c <- .wtd_mean(x = nuis$mu.0c, w = cace.wt)

        tau.1n <- .wtd_mean(x = nuis$y,     w = nace.wt)
        tau.0n <- .wtd_mean(x = nuis$mu.0n, w = nace.wt)
    }



    if (estimator=="MS") {

        cace.wt <- nuis$s.wt* (nuis$z*nuis$z.wt* (nuis$c - nuis$c.wt) + nuis$c.wt)
        nace.wt <- nuis$s.wt* (nuis$z*nuis$z.wt* (nuis$n - nuis$n.wt) + nuis$n.wt)

        tau.1c <- .wtd_mean(x = nuis$mu.1c, w = cace.wt)
        tau.0c <- .wtd_mean(x = nuis$mu.0c, w = cace.wt)

        tau.1n <- .wtd_mean(x = nuis$mu.1n, w = nace.wt)
        tau.0n <- .wtd_mean(x = nuis$mu.0n, w = nace.wt)
    }



    if (estimator=="IF" | estimator=="IFH") {

        z <- nuis$z
        c <- nuis$c
        n <- nuis$n
        y <- nuis$y

        z.wt <- nuis$z.wt
        c.wt <- nuis$c.wt
        n.wt <- nuis$n.wt

        mu.1c <- nuis$mu.1c
        mu.1n <- nuis$mu.1n
        mu.0c <- nuis$mu.0c
        mu.0n <- nuis$mu.0n

        if (estimator=="IFH") {
            z.wt <- z.wt / (z*.wtd_mean(z*z.wt, nuis$s.wt) +
                                (1-z)*.wtd_mean((1-z)*z.wt, nuis$s.wt))
        }

        phi_nu.1c <- c*z.wt*(y-mu.1c) + z*z.wt*mu.1c*(c-c.wt) + c.wt*mu.1c
        phi_nu.1n <- n*z.wt*(y-mu.1n) + z*z.wt*mu.1n*(n-n.wt) + n.wt*mu.1n

        phi_nu.0c <- (1-z)*z.wt*c.wt*(y-mu.0c) + z*z.wt*mu.0c*(c-c.wt) + c.wt*mu.0c
        phi_nu.0n <- (1-z)*z.wt*n.wt*(y-mu.0n) + z*z.wt*mu.0n*(n-n.wt) + n.wt*mu.0n

        phi_delta.c <- z*z.wt*(c-c.wt) + c.wt # delta.c here means E[pi.c(X)]

        delta.c <- .wtd_mean(x = phi_delta.c, w = nuis$s.wt)
        delta.n <- 1 - delta.c

        nu.1c <- .wtd_mean(x = phi_nu.1c, w = nuis$s.wt)
        nu.1n <- .wtd_mean(x = phi_nu.1n, w = nuis$s.wt)
        nu.0c <- .wtd_mean(x = phi_nu.0c, w = nuis$s.wt)
        nu.0n <- .wtd_mean(x = phi_nu.0n, w = nuis$s.wt)

        tau.1c <- nu.1c / delta.c
        tau.1n <- nu.1n / delta.n
        tau.0c <- nu.0c / delta.c
        tau.0n <- nu.0n / delta.n
    }



    cace <- tau.1c - tau.0c
    nace <- tau.1n - tau.0n

    c(cace = cace,
      nace = nace,
      tau.1c = tau.1c,
      tau.1n = tau.1n,
      tau.0c = tau.0c,
      tau.0n = tau.0n)
}


