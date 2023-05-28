####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .get_boot_estimates ####

#' (internal) Compute bootstrap estimates
#'
#' Implement bootstrap to compute bootstrap estimates.
#'
#' @param FUN The function that will act on \code{data} to implement the estimator. The output of \code{FUN} should be in the form of a scalar, a vector or a matrix, NOT a list.
#' @param data A data frame to be used by \code{FUN}. If sampling weights are used, they should come in a variable named \code{s.wt}.
#' @param boot.num Number of bootstrap samples. Default is 999.
#' @param boot.seed Seed for reproducibility. Default is 12345.
#' @param boot.cont.wt If \code{TRUE} (default), implement continuous weights bootstrap (based on the Dirichlet distribution). If \code{FALSE}, implement the resampling bootstrap.
#' @param double.boot If \code{FALSE} (default), implement the standard bootstrap. If \code{TRUE}, implement a single double-bootstrap. This is used for correction of finite sample bias.
#' @return Bootstrap estimates in the structure of a vector if \code{FUN} outputs a scalar, a matrix if \code{FUN} outputs a vector, a 3D array if \code{FUN} outputs a matrix. In the last two cases, the last dimension of the array is the dimension of the bootstrap samples.
#'
#' @keywords internal

.get_boot_estimates <- function(FUN,
                                data,
                                boot.num = 999,
                                boot.cont.wt = TRUE,
                                double.boot = FALSE,
                                boot.seed = 12345) {


    if (all(names(data) != "s.wt")) data$s.wt <- 1



    # FIRST, get boot weights for specific bootstrap method

    samp.size <- nrow(data)

    set.seed(boot.seed)

    # for standard bootstrap, get the whole matrix of boot weights
    # (each row is a bootstrap draw, all draws are based on original sample)
    if (boot.cont.wt){
        boot.wt <- samp.size * gtools::rdirichlet(n     = boot.num,
                                                  alpha = data$s.wt)
    } else {
        boot.wt <- t(stats::rmultinom(n    = boot.num,
                                      size = samp.size,
                                      prob = data$s.wt))
    }


    if (double.boot) {

        # for double bootstrap, get boot weights row by row
        # (each row is a single draw based on a different first-level bootstrap sample)
        boot.wt <- t(sapply(1:boot.num, function(z) {

            if (boot.cont.wt) {
                c(samp.size * gtools::rdirichlet(n = 1,
                                                 alpha = boot.wt[z,]))
            } else {
                c(stats::rmultinom(n = 1,
                                   size = samp.size,
                                   prob = boot.wt[z,]))
            }
        }))
    }



    # THEN, compute boot estimates

    # (manual simplification to remove NULL elements from failed boot samples)
    ests <- sapply(1:boot.num,
                   function(z) {
                       dat <- data
                       dat$s.wt <- boot.wt[z,]
                       out <- tryCatch(FUN(data = dat),
                                       error = function(e) return(NULL))
                       out
                   },
                   simplify = FALSE)

    ests <- ests[lengths(ests) != 0]

    simplify2array(ests)
}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .get_ahat_for_BCa ####

#' (internal) Estimate acceleration parameter for construction of BCa intervals
#'
#' @inheritParams .get_boot_estimates
#' @return Estimated acceleration params in the same format as the output of \code{FUN}, i.e., one ahat for each parameter estimate from \code{FUN}.
#' @keywords internal

.get_ahat_for_BCa <- function(FUN, data) {

    jacks <- .get_jack_estimates(FUN = FUN,
                                 data = data)

    if (is.vector(jacks)) {
        ahat <- .ahat_from_jacks(jacks)

    } else if (length(dim(jacks))==2) {
        ahat <- apply(jacks, 1, .ahat_from_jacks)

    } else if (length(dim(jacks))==3) {
        ahat <- apply(jacks, c(1, 2), .ahat_from_jacks)
    }

    ahat
}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .ahat_from_jacks ####

#' (internal) Compute acceleration parameter based on jackknife estimates
#'
#' This function is called by function \code{.get_ahat_for_BCa()}. Its purpose is to avoid having to retype the unintuitive formula.
#'
#' @param jacks The set of all jackknife estimates.
#' @return The acceleration parameter.
#' @keywords internal

.ahat_from_jacks <- function(jacks) {

    jacks <- jacks[!is.na(jacks)]

    jack.devs <- mean(jacks) - jacks

    sum(jack.devs^3) / (6 * (sum(jack.devs^2))^(3/2))
}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .get_jack_estimates ####

#' (internal) Obtain jackknife estimates for an estimator
#'
#' Compute jackknife estimates.
#'
#' @inheritParams .get_boot_estimates
#' @return Jackknife estimates in the structure of a vector if \code{FUN} outputs a scalar, a matrix if \code{FUN} outputs a vector, a 3D array if \code{fun} outputs a matrix. In the last two cases, the last dimension of the array is the dimension of the bootstrap samples.
#' @keywords internal


.get_jack_estimates <- function(FUN, data) {

    jacks <- sapply(1:nrow(data),
                    function(z) {
                        out <- tryCatch(FUN(data = data[-z,]),
                                        error = function(e) return(NULL))
                        out
                    },
                    simplify = FALSE)

    jacks <- jacks[lengths(jacks) != 0]

    simplify2array(jacks)

}


####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .get_BCa.interval ####

#' (internal) Compute a BCa interval based on a point estimate, a set of boot estimates and an acceleration parameter
#'
#' @param theta.boot A vector of boot estimates.
#' @param theta.hat The point estimate.
#' @param alpha \code{(1-alpha)} is the confidence level.
#' @param ahat The acceleration parameter.
#' @return The BCa interval.
#' @keywords internal

.get_BCa.interval <- function(theta.boot, theta.hat, alpha = .05, ahat) {

    theta.boot <- theta.boot[!is.na(theta.boot)]

    l <- alpha/2
    h <- 1 - alpha/2

    z0 <- stats::qnorm(mean(theta.boot<theta.hat))

    bca.l <- stats::pnorm(z0 + (z0+stats::qnorm(l))/(1 - ahat*(z0+stats::qnorm(l))))
    bca.h <- stats::pnorm(z0 + (z0+stats::qnorm(h))/(1 - ahat*(z0+stats::qnorm(h))))

    stats::quantile(theta.boot, c(bca.l, bca.h))
}







