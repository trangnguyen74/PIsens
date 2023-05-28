##############################
#### PIsens ----

#' Sensitivity analysis for violation of principal ignorability (PI)
#'
#' @param data A data frame that has been processed by function \code{declare_data()}.
#' @param estimator One of the estimators listed in the paper: "epi", "pimu" or "emu" (for the simple estimator \eqn{e\pi}, \eqn{\pi\mu} or \eqn{e\mu}), or "MS", "IF" or "IFH" (for the multi-step, IF-as-estimating-function, or Hajek-ized estimator).
#' @param targeted If \code{TRUE}, the principal score model and the Y1 models are fit to inverse propensity score weighted treatment arm data; and the Y0 models are fit to control arm data that are inverse probability of control weighted and principal score weighted. If \code{FALSE}, such weighting is not used.
#' @param Z.form A character string specifying the formula for the propensity score model.
#' @param C.form A character string specifying the formula for the principal score model.
#' @param Y.form Either a character string specifying one formula to be used for all outcome model(s), or a character vector of length 3 specifying formulas for the outcome in compliers in the treatment arm (first element), in noncompliers in the treatment arm (second element), and in controls (third element.).
#' @param Y.model The type of model to be used for the outcome. Can be (1) a character string or (2) an object of class family. In case (1), if "linear", "gaussian" or "identity", will fit a linear model; if "binomial", "quasibinomial", "logit", "logistic" or "odds ratio", will fit a quasibinomial model with logit link; if "poisson",  "quasipoisson", "count" or "rate ratio", will fit a quasipoisson model with log link; if "Gamma" or "gamma", will fit a quasi model with log link and variance proportional to mean squared.
#' @param Y.bounds A two element vector for the lower and upper bounds of the outcome. Needed when modeling a bounded outcome by a quasibinomial model and/or when conducting a GOR-based sensitivity analysis.
#' @param Y.FUN Optional function for customized estimation of outcome mean functions. For an example, see function \code{jobsII_earn3_mus()} in this package, and see how it is used in the vignette. If use \code{Y.FUN}, arguments for outcome model forms and family do not need to be specified, and are ignored.
#' @param sens.type Name of the sensitivity parameter. Options are: "OR", "GOR", "MR", "MD", "beta.quant".
#' @param sens.range A numeric vector of length 2, holding the minimum and maximum values for the sensitivity parameter. If not specified, defaults to (1/3, 3) for \code{sens.type} being OR, GOR or MR; (-1, 1) for \code{sens.type} being MD; and (-.05, .05) for \code{sens.type} being "beta.quant".
#' @param sens.step The width of the step by which to increase the sensitivity parameter. If not specified, defaults to 0.05 for mean-based sensitivity analyses and 0.05 for distribution-based sensitivity analyses.
#' @param boot.num Number of bootstrap samples. Default is 999.
#' @param boot.seed Seed for reproducibility. Default is 12345.
#' @param boot.cont.wt If \code{TRUE} (default), implement continuous weights bootstrap (with weights based on the Dirichlet distribution). If \code{FALSE}, implement classic bootstrap (with weights based on the multinomial distribution). If
#' @param double.boot Whether to implement bias correction of point estimate by double bootstrap. If FALSE, implements bias correction of point estimate by single bootstrap only; this substantially reduces computational burden for the beta-quantile sensitivity analysis. Default is FALSE.
#' @param BCa.CI Whether to compute BCa confidence intervals. If FALSE, compute percentile intervals only, which is not recommended. However, BCa CIs are computationally burdensome for the beta-quantile sensitivity analysis. Default is TRUE.
#' @export

# TODO: add bounds for GOR (part of sens.inputs)
# TODO: add options for BCa or not, double.boot or not

PIsens <- function(
        data,
        estimator,
        targeted,

        Z.form,
        C.form,

        Y.form,
        Y.model,
        Y.bounds,
        Y.FUN = NULL,

        sens.type,
        sens.range = NULL,
        sens.step = NULL,

        boot.num = 999,
        boot.seed = 12345,
        boot.cont.wt = TRUE,
        double.boot = FALSE,
        BCa.CI = TRUE) {


    msg1 <- "data needs to go through declare_data() before use with PImain() or PIsens()."
    msg2 <- "estimator must be specified."
    msg3 <- "Untargeted nuisance estimation is not implemented for the MS estimator. Proceeding with targeted = TRUE."
    msg4 <- "For targeted nuisance estimation, must specify both Z.form and C.form."
    msg5 <- "For sens.type GOR, Y.bounds must be specified."



    ########################################
    # Check if data are declared data

    if (!.is.declared(data)) stop(msg1)



    ########################################
    # Assemble key inputs: estimator and targeted

    if (missing(estimator)) stop(msg2)

    if (estimator=="MS") {
        if (missing(targeted)) {
            targeted <- TRUE
        } else if (targeted==FALSE) {
            message(msg3)
            targeted <- TRUE
        }
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
                                    sensitivity = TRUE,
                                    sens.type = sens.type,
                                    nuis.inputs))


    ########################################
    # Assemble sens inputs

    if (sens.type=="SMD") sens.type <- "SMDe"

    if (is.null(sens.range)) {
        if (sens.type %in% c("OR", "GOR", "MR")) { sens.range = c(1/3, 3)
        } else if (sens.type=="SMDe")            { sens.range = c(-1, 1)
        } else                                   { sens.range = c(-.05, .05)
        }
    }

    if (is.null(sens.step)) {
        if (sens.type %in% c("OR", "GOR", "MR")) { sens.step = 0.05
        } else if (sens.type=="SMDe")            { sens.step = 0.05
        } else                                   { sens.step = 0.005
        }
    }

    sens.inputs <- list(sens.type = sens.type,
                        sens.range = range(sens.range),
                        sens.step = sens.step)

    if (sens.type=="GOR") {
        if (missing(Y.bounds)) stop(msg5)
        sens.inputs[["Y.bounds"]] <- Y.bounds
    }




    ########################################
    # Assemble boot inputs

    boot.inputs <- mget(c("boot.num", "boot.seed", "boot.cont.wt",
                          "double.boot", "BCa.CI"))



    ########################################
    # Construct function do_one() to get point estimates from one dataset

    do_one <- function(data,
                       .key.inputs  = key.inputs,
                       .nuis.inputs = nuis.inputs,
                       .sens.inputs = sens.inputs) {

        tmp <- do.call(".estimate_nuisance",
                       args = c(data = list(data),
                                .key.inputs,
                                .nuis.inputs))

        if (is.data.frame(tmp)) {
            nuis <- tmp
        } else                  {
            nuis   <- tmp$nuis
            quants <- tmp$quants
        }
        rm(tmp)



        main <- do.call(".assemble_main",
                        args = c(estimator = .key.inputs$estimator,
                                 nuis      = list(nuis)))



        if (sens.type %in% c("OR", "GOR", "MR", "SMDe")) {

            do.call(".assemble_sens",
                    args = c(estimator = .key.inputs$estimator,
                             main      = list(main),
                             nuis      = list(nuis),
                             .sens.inputs))

        } else if (sens.type=="beta.quant") {

            do.call(".assemble_sens",
                    args = c(estimator = .key.inputs$estimator,
                             main      = list(main),
                             nuis      = list(nuis),
                             quants    = list(quants),
                             .sens.inputs))
        }


    }


    ########################################
    # Run function do_one() through .do_all_main() to bootstrap and all

    results <- do.call(".do_all_sens",
                       args = c(pointFUN = do_one,
                                data     = list(data),
                                boot.inputs))




    ########################################
    # add sensitivity parameter columns

    if (sens.type %in% c("OR", "GOR", "MR")) {
        rho   <- .get_rho.vec(sens.range, sens.step)
        rho.x <- .get_rho.coords(sens.range, sens.step)

        if (boot.num==0) {
            results$sens <- cbind(sens.para = rho,
                                  sens.para.coord = rho.x,
                                  results$sens)
        } else {
            results$sens <- lapply(results$sens,
                                   function(z) cbind(sens.para = rho,
                                                     sens.para.coord = rho.x,
                                                     z))
        }


    }

    if (sens.type=="SMDe") {
        eta <- .get_eta.vec(sens.range, sens.step)

        if (boot.num==0) {
            results$sens <- cbind(sens.para = eta,
                                  sens.para.coord = eta,
                                  results$sens)
        } else {
            results$sens <- lapply(results$sens,
                                   function(z) cbind(sens.para = eta,
                                                     sens.para.coord = eta,
                                                     z))
        }
    }

    if (sens.type=="beta.quant") {
        s.kappa <- .get_signed.kappa.vec(sens.range, sens.step)

        if (boot.num==0) {
            results$sens <- cbind(sens.para = s.kappa,
                                  sens.para.coord = s.kappa,
                                  results$sens)
        } else {
            results$sens <- lapply(results$sens,
                                   function(z) cbind(sens.para = s.kappa,
                                                     sens.para.coord = s.kappa,
                                                     z))
        }
    }


    results[["sens.specs"]] <- list(sens.type  = sens.type,
                                    sens.range = sens.range)


    results
}




##############################
#### .do_all_sens ----

# TODO: update documentation
# TODO (or not): consider starting first with boot1 and double sapply to compute list of CI matrices then boot2 to bias correct point estimates, and finally merge point estimates in the list

#' (internal) Compute confidence intervals and bias-correct point estimates for a sensitivity analysis
#'
#' @param pointFUN A function constructed to compute the point estimate vector (one estimate for each value of the sensitivity parameter), with all options wrapped in, so the only argument of this function is \code{data}. This pointFUN function should output a two-column matrix with CACE and NACE in the first and second columns, respectively.
#' @param long If TRUE, print progress messages to screen. Default is FALSE.
#' @inheritParams PIsens
#' @return A list of two matrices for CACE and NACE, respectively. Each matrix contains columns for the point estimate, the bias-corrected point estimate, the percentile CI, the BCa CI, and also the bootstrap and double bootstrap means of the estimates.
#' @keywords internal


.do_all_sens <- function(pointFUN,
                         data,
                         boot.num,
                         boot.seed = 12345,
                         boot.cont.wt = TRUE,
                         double.boot = TRUE,
                         BCa.CI = TRUE) {

    # get point estimates:
    # this is a matrix with columns for parameters and rows for rho values
    point <- pointFUN(data)

    if (boot.num==0) return(list(sens = point))




    successful.boot.num <- list()


    # the order of what follows is a bit counter-intuitive but
    # its purpose is to minimize memory use at any time point


    # double-bootstrap (for finite-sample bias correction of point estimates)
    if (double.boot) {
            message(paste0("\nStarting the double bootstrap.",
                           "\nThis may take a while.",
                           "\nWhile R runs, please keep your computer from going to sleep.\n"))

        boot2 <- .get_boot_estimates(FUN          = pointFUN,
                                     data         = data,
                                     boot.num     = boot.num,
                                     boot.cont.wt = boot.cont.wt,
                                     boot.seed    = boot.seed,
                                     double.boot  = TRUE)

        boot2.mean <- apply(boot2, c(1,2), function(z) mean(z, na.rm = TRUE))
        successful.boot.num[["boot2"]] <- dim(boot2)[3]
        rm(boot2)
    }

    gc()


    # estimate acceleration parameters to be used in BCa confidence interval estimation
    # (this requires jackknife-ing, but we retain only a matrix of ahat paras)
    if (BCa.CI) {
        message(paste0("\nStarting jackknife to estimate acceleration parameter for BCa confidence interval estimation.",
                       "\nThis may take a while.",
                       "\nWhile R runs, please keep your computer from going to sleep.\n"))

        ahat <- .get_ahat_for_BCa(FUN = pointFUN, data = data)
    }

    gc()



    # now, the usual bootstrap

    if (double.boot) {
        message(paste0("\nStarting the standard bootstrap.",
                       "\nThis will take as long as the double bootstrap step, but after this we will be done.",
                       "\nWhile R runs, please keep your computer from going to sleep.\n"))
    } else {
        message(paste0("\nStarting the bootstrap.",
                       "\nThis may take a while, but after this we will be done.",
                       "\nWhile R runs, please keep your computer from going to sleep.\n"))
    }


    boot1 <- .get_boot_estimates(FUN          = pointFUN,
                                 data         = data,
                                 boot.num     = boot.num,
                                 boot.cont.wt = boot.cont.wt,
                                 boot.seed    = boot.seed)

    boot1.mean <- apply(boot1, c(1,2), function(z) mean(z, na.rm = TRUE))
    successful.boot.num[["boot1"]] <- dim(boot1)[3]


    # point estimate bias correction is easy, so let's do that first
    point.bc1 <- 2*(point - boot1.mean) + boot1.mean

    if (double.boot) point.bc2 <- 3*(point - boot1.mean) + boot2.mean



    # lastly, wrap the computation of confidence intervals in two sapply layers
    # also combine the computed CIs with the point estimates for output
    #
    # this returns a list (note the outer sapply uses simplify = FALSE), where
    # each element corresponds to a target parameter (cace, nace, tau.0c, tau.0n)
    # and each is a matrix with columns for point and interval estimates and
    # rows corresponding to rho values


    sens <- sapply(colnames(point), function(z) {

        # percentile confidence interval
        c.intervals <- t(sapply(1:nrow(point), function(u) {

            stats::quantile(boot1[u, z, ],
                            prob = c(0.025, 0.975),
                            na.rm = TRUE)
        }))
        colnames(c.intervals) <- c("perc.low", "perc.high")


        # BCa confidence interval
        if (BCa.CI) {
            ci.BCa <- t(sapply(1:nrow(point), function(u) {
                .get_BCa.interval(theta.boot = boot1[u, z, ],
                                  theta.hat  = point[u, z],
                                  ahat       = ahat [u, z])
            }))
            colnames(ci.BCa) <- c("bca.low", "bca.high")

            c.intervals <- cbind(ci.BCa, c.intervals)
        }


        # grab all outputs: point estimates (raw and bias-corrected) and CIs
        if (double.boot) {

            sens.each <- cbind(point      = point[, z],
                               point.bc1  = point.bc1[, z],
                               point.bc2  = point.bc2[, z],
                               c.intervals,
                               boot.mean  = boot1.mean[, z],
                               boot2.mean = boot2.mean[, z])
        } else {

            sens.each <- cbind(point      = point[, z],
                               point.bc1  = point.bc1[, z],
                               c.intervals,
                               boot.mean  = boot1.mean[, z])

        }

        sens.each

    },
    simplify = FALSE)


    rm(boot1)

    if (BCa.CI) rm(ahat)

    gc()


    mget(c("sens", "successful.boot.num"))

}




#############################
#### .assemble_sens ----

#' (internal) Assemble sensitivity analysis for PI violation for a specified estimator based on main analysis estimates and estimated nuisance functions under a sensitivity assumption.
#'
#' @inheritParams PIsens
#' @inheritParams .assemble_main
#' @param main The main analysis result from \code{.assemble_main()}
#' @param quants A list of two matrices of conditional outcome quantiles, for estimation of \code{mu.0c} and \code{mu.0n} if \code{sens.type} is "beta.quant".
#' @return A matrix of point estimates with columns named \code{cace}, \code{nace}, \code{tau.0c}, \code{tau.0n} and rows for values of the sensitivity parameter.
#' @keywords internal

.assemble_sens <- function(estimator,
                           main,
                           nuis,
                           quants,
                           sens.type,
                           sens.range,
                           sens.step,
                           Y.bounds = NULL) {

    if (sens.type=="beta.quant") {

        inputs <- mget(c("estimator", "main", "nuis", "quants",
                         "sens.range", "sens.step"))

        sens <- do.call(".assemble_beta.quant", inputs)

    } else {

        inputs <- mget(c("estimator", "main", "nuis",
                         "sens.range",  "sens.step"))

        if (sens.type=="GOR") inputs[["Y.bounds"]] <- Y.bounds

        if (sens.type %in% c("OR", "GOR")) sens <- do.call(".assemble_GOR", inputs)

        if (sens.type=="MR")  sens <- do.call(".assemble_MR", inputs)

        if (sens.type=="SMDe")  sens <- do.call(".assemble_SMDe", inputs)

    }

    sens
}





