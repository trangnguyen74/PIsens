

.cleanMem <- function(n=10) { for (i in 1:n) gc() }




##############################
#### .trans_bounds

#' Transformation of a bounded variable (vector) to a different pair of bounds by shifting and scaling
#'
#' @param x The vector to be transformed
#' @param from The pair of bounds of `x`
#' @param to The destination pair of bounds
#' @keywords internal

.trans_bounds <- function(x, from, to) {
    (x - min(from)) / abs(diff(from)) * abs(diff(to)) + min(to)
}







##############################
#### .is.01

#' Check if a vector is binary and coded as 0/1
#'
#' @param x A vector
#' @keywords internal

.is.01 <- function(x) {

    x <- x[!is.na(x)]

    all(x %in% c(0,1))
}








##############################
#### .replace_y

#' (internal) Replace the response variable in a model formula
#'
#' @param formula A model formula
#' @param new.y Name of new response variable
#' @keywords internal

.replace_y <- function(formula, new.y) {
    stats::reformulate(deparse(formula(formula)[[3]]), response = new.y)
}




##############################
#### .wtd_mean

#' Weighted mean function from the stats package
#'
#' @param x Vector of numeric values to be averaged
#' @param w Vector of weights
#' @keywords internal

.wtd_mean <- function(x, w) {
    stats::weighted.mean(x, w)
}




##############################
#### .get_rho.vec

#' Compute the vector of rho values to use for the sensitivity analysis, for a pair of min and max values
#' @inheritParams .assemble_sens
#' @keywords internal

.get_rho.vec <- function(sens.range, sens.step) {

    rho.min <- min(sens.range)
    rho.max <- max(sens.range)

    if (rho.min < 1) {
        left.vec <- 1/ seq(1, 1/rho.min, sens.step)
        if (!rho.min %in% left.vec) left.vec <- c(left.vec, rho.min)
        left.vec <- left.vec[-1]
        left.vec <- left.vec[length(left.vec):1]

    } else {
        left.vec <- NULL
    }

    if (rho.max > 1) {
        right.vec <- seq(1, rho.max, sens.step)
        if (!rho.max %in% right.vec) right.vec <- c(right.vec, rho.max)
        right.vec <- right.vec[-1]

    } else {
        right.vec <- NULL
    }

    c(left.vec, 1, right.vec)
}




##############################
#### .get_eta.vec

#' Compute the vector of eta values to use for the sensitivity analysis, for a pair of min and max values
#' @inheritParams .assemble_sens
#' @keywords internal

.get_eta.vec <- function(sens.range, sens.step) {

    eta.min <- min(sens.range)
    eta.max <- max(sens.range)

    if (eta.min < 0) {
        left.vec <- seq(0, eta.min, -sens.step)
        if (!eta.min %in% left.vec) left.vec <- c(left.vec, eta.min)
        left.vec <- left.vec[-1]
        left.vec <- left.vec[length(left.vec):1]
    } else {
        left.vec <- NULL
    }

    if (eta.max > 0) {
        right.vec <- seq(0, eta.max, sens.step)
        if (!eta.max %in% right.vec) right.vec <- c(right.vec, eta.max)
        right.vec <- right.vec[-1]
    } else {
        right.vec <- NULL
    }

    c(left.vec, 0, right.vec)
}




##############################
#### .get_rho.coords

#' Compute x-axis plot values corresponding to rho values
#'
#' @inheritParams .get_rho.vec
#' @keywords internal

.get_rho.coords <- function(sens.range, sens.step) {

    rho <- .get_rho.vec(sens.range, sens.step)

    (rho<1)*(-1/rho + 1) + (rho>=1)*(rho - 1)
}




##############################
#### .get_signed.kappa.vec

#' Compute the vector of signed-kappa values to use for the sensitivity analysis, for a pair of min and max values
#' @inheritParams .assemble_sens
#' @keywords internal

.get_signed.kappa.vec <- function(sens.range, sens.step) {

    s.kappa.min <- sens.range[1]
    s.kappa.max <- sens.range[2]

    if (s.kappa.min < 0) {
        left.vec <- seq(0, s.kappa.min, -sens.step)
        if (!s.kappa.min %in% left.vec) left.vec <- c(left.vec, s.kappa.min)
        left.vec <- left.vec[-1]
        left.vec <- left.vec[length(left.vec):1]

    } else {
        left.vec <- NULL
    }

    if (s.kappa.max > 0) {
        right.vec <- seq(0, s.kappa.max, sens.step)
        if (!s.kappa.max %in% right.vec) right.vec <- c(right.vec, s.kappa.max)
        right.vec <- right.vec[-1]

    } else {
        right.vec <- NULL
    }

    c(left.vec, 0, right.vec)
}




##############################
#### .extract_dispersion_wtd_glm

#' Extract the conditional distribution's dispersion param from a weighted GLM
#' Compute dispersion parameter from weighted models ignoring variance inflation due to weighting. Such variance inflation is desired when the dispersion parameter is used to compute SEs of regression coefficients. It is not desired when the dispersion parameter is used to estimate the conditional variance of the response variable, ie `var(Y|X)`.
#' @param mod The fitted linear model object from `lm` or `glm`
#' @keywords internal

.extract_dispersion_wtd_glm <- function(mod) {

    # MANUALLY compute squared Pearson residuals
    # Remember NOT TO USE mod$residuals or residuals(mod, "pearson")
    sq.pearson.resids <-
        (mod$y - mod$fitted.values)^2 / mod$family$variance(mod$fitted.values)

    .wtd_mean(sq.pearson.resids, mod$weights) *
        length(mod$fitted.values) / mod$df.residual
}


# function to check if a numeric vector (with no mising data) is all positive
is.positive <- function(x) {
    mean(x>0) == 1
}

# function to check if a numeric vector (with no missing data) is bounded in (0, 1)
is.within.bounds <- function(x, bounds) {
    sum(x>max(bounds)) + sum(x<min(bounds)) == 0
}




##############################
#### .get_mu_hat

#' Compute mu.hat based on outcome model inputs
#' @param form Model formula
#' @param data Data for model fitting.
#' @param wts A vector of weights for model fitting.
#' @param y.fam Model family for `glm()`.
#' @param pred.dat Data to predict mu.hat on
#' @param get.dispersion Whether to also compute the dispersion parameter based on the quasi family.
#' @keywords internal

.get_mu_hat <- function(form,
                        data,
                        wts,
                        y.fam,
                        pred.dat,
                        get.dispersion = FALSE) {

    suppressWarnings({
        mod <- do.call("glm",
                       list(formula = form,
                            data    = data,
                            weights = wts,
                            family  = y.fam))
    })

    out <- list(mu = predict(mod, newdata = pred.dat, type = "response"))

    if (get.dispersion)
        out[["dispersion"]] <- .extract_dispersion_wtd_glm(mod)

    out
}




##############################
#### .df.rename

#' Rename variables in a data frame
#' @param x The data frame.
#' @param from Names of the variables to be renamed.
#' @param to New names to give the variables.
#' @keywords internal

.df.rename <- function(x, from, to) {

    if (length(from)!=length(to))
        stop("Vector 'from' and vector 'to' must be of the same length.")

    if (length(from)==1) {
        names(x)[names(x)==from] <- to

    } else {
        for(z in 1:length(from)) {
            names(x)[names(x)==from[z]] <- to[z]
        }
    }

    x

}




##############################
#### data_for_Y.FUN

#' Make the dataset that will be handled by \code{Y.FUN} in \code{PImain()} and \code{PIsens()}
#'
#' @inheritParams PIsens
#' @export

data_for_Y.FUN <- function(data,
                           targeted,
                           Z.form,
                           C.form) {

    if (!.is.declared(data)) stop("data needs to go through declare_data().")

    data$group <- ifelse(data$z==0, "0",
                         ifelse(data$c==1, "1c", "1n"))

    if (targeted) {

        nuis <- .estimate_nuisance(data = data,
                                   Z.form = Z.form,
                                   C.form = C.form,
                                   targeted = targeted)

        data$z.ipw <- nuis$z.wt
        data$c.ps  <- nuis$c.ps
        data$n.ps  <- nuis$n.ps
    }


    data
}



















