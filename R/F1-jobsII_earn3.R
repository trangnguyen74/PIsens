
##############################
#### earn3_FUN1

#' JOBSII-specific: estimate \eqn{\mu} functions for the outcome \code{earn3}
#'
#' Customized estimation of \eqn{\mu} functions for the outcome \code{earn3} in the paper data example \code{jobsII}. \code{earn2_FUN1()} is the mother function, from which we derive child functions, eg \code{earn3_FUN1_IFH()}, which are fed into the \code{Y.FUN} argument of \code{PImain()} or \code{PIsens()}.
#'
#' @details Important note: This function satisfies what is required of argument \code{Y.FUN} in function \code{PIsens()} uses the \code{s.wt} variable in \code{data} in the estimation. This is required of any custom function for \eqn{\mu} estimation, even if this variable is all 1 (i.e., no sampling weights), because this variable will be used by the bootstrap to incorporate bootstrap weights.
#'
#' @order 1
#' @inheritParams .estimate_nuisance
#' @param data A data frame that is produced within the calling function (\code{PImain()} or \code{PIsens}). See \code{Y.FUN} in documentation of \code{PImain()} or \code{PIsens()} for description of this type of dataset.
#' @param include.mu.1c If FALSE, only estimate mu.0(X), not mu.1c(X) and mu.1n(X).
#' @return A data frame with variable \code{id} and the mu functions.
#' @export



earn3_FUN1 <- function(data,
                       targeted,
                       include.mu.1c) {
    dat <- data; rm(data)

    dat$c <- ifelse(dat$z==0, 0, dat$c)
    dat$n <- ifelse(dat$z==0, 0, 1 - dat$c)

    Y.form1 <-
        "work3 ~ age + sex + race + edu + marital + hh.kids + hh.income + econ.hard + occu + wks.unemp + part.motiv + seek.motiv + seek.effi + assertive + depress1"

    Y.form2 <-
        "y ~ age + sex + race + edu + marital + hh.kids + hh.income + econ.hard + occu + wks.unemp + part.motiv + seek.motiv + seek.effi + assertive + depress1"

    mu_convenient <- function(dat1,
                              dat2,
                              form1 = Y.form1,
                              form2 = Y.form2,
                              wts1,
                              wts2,
                              pred.dat = dat) {

        # fit models
        suppressWarnings({
            mod1 <- do.call("glm",
                            list(formula = form1,
                                 data    = dat1,
                                 weights = wts1,
                                 family  = "quasibinomial"))

            mod2 <- do.call("glm",
                            list(formula = form2,
                                 data = dat2,
                                 weights = wts2,
                                 family = stats::Gamma(link = "log")))
        })


        # compute the mu function
        pred1 <- predict(mod1, newdata = pred.dat, type = "response")
        pred2 <- predict(mod2, newdata = pred.dat, type = "response") *
            .wtd_mean(mod2$y,             mod2$weights) /
            .wtd_mean(mod2$fitted.values, mod2$weights)

        pred2 <- pmin(6000, pred2)


        pred1 * pred2
    }


    mus <- dat[names(dat)=="id"]


    ########################################
    # mu.1c and mu.1n

    if (include.mu.1c) {

        # mu.1c
        dat1 <- dat[dat$group=="1c", ]
        dat2 <- dat1[dat1$work3==1, ]

        wts1 <- dat1$s.wt
        wts2 <- dat2$s.wt

        if (targeted) {
            wts1 <- wts1 * dat1$z.ipw
            wts2 <- wts2 * dat2$z.ipw
        }

        mus$mu.1c <- mu_convenient(dat1 = dat1,
                                   dat2 = dat2,
                                   wts1 = wts1,
                                   wts2 = wts2)
        rm(dat1, dat2, wts1, wts2)

        # mu.1n
        dat1 <- dat[dat$group=="1n", ]
        dat2 <- dat1[dat1$work3==1, ]

        wts1 <- dat1$s.wt
        wts2 <- dat2$s.wt

        if (targeted) {
            wts1 <- wts1 * dat1$z.ipw
            wts2 <- wts2 * dat2$z.ipw
        }

        mus$mu.1n <- mu_convenient(dat1 = dat1,
                                   dat2 = dat2,
                                   wts1 = wts1,
                                   wts2 = wts2)
        rm(dat1, dat2, wts1, wts2)

    }


    ########################################
    # mu.0

    dat1 <- dat[dat$group=="0", ]
    dat2 <- dat1[dat1$work3==1, ]


    if (targeted) {

        # mu.0c
        wts1 <- dat1$s.wt * dat1$z.ipw * dat1$c.ps
        wts2 <- dat2$s.wt * dat2$z.ipw * dat2$c.ps

        mus$mu.0c <- mu_convenient(dat1 = dat1,
                                   dat2 = dat2,
                                   wts1 = wts1,
                                   wts2 = wts2)
        rm(wts1, wts2)


        # mu.0n
        wts1 <- dat1$s.wt * dat1$z.ipw * dat1$n.ps
        wts2 <- dat2$s.wt * dat2$z.ipw * dat2$n.ps

        mus$mu.0n <- mu_convenient(dat1 = dat1,
                                   dat2 = dat2,
                                   wts1 = wts1,
                                   wts2 = wts2)

    } else {

        # both mu.0c and mu.0n
        mus$mu.0c <- mus$mu.0n <- mu_convenient(dat1 = dat1,
                                                dat2 = dat2,
                                                wts1 = dat1$s.wt,
                                                wts2 = dat2$s.wt)

    }


    mus
}





#' @rdname earn3_FUN1
#' @order 2


earn3_FUN1_IFH <-function(data) {
    earn3_FUN1(data = data, targeted = TRUE, include.mu.1c = TRUE)
}











##############################
#### earn3_FUN2

#' JOBSII-specific: estimate \eqn{\mu} functions (and optionally, quantiles) for the outcome \code{earn3}
#'
#' Customized estimation of \eqn{\mu} functions and conditional \eqn{Y|X,Z=0} quantiles for the outcome \code{earn3} in the paper data example \code{jobsII}
#'
#' @details Important note: This function uses the \code{s.wt} variable in \code{data} in the estimation. This is required of any custom function for \eqn{\mu} estimation, even if this variable is all 1 (i.e., no sampling weights), because this variable will be used by the bootstrap to incorporate bootstrap weights.
#'
#' @order 1
#' @inheritParams .estimate_nuisance
#' @param data A data frame that has been processed by \code{declare_data()} and --- if needed for targeted \eqn{\mu} estimation --- has been enhanced by (1) having inverse probability of assigned treatment weight incorporated (ie multiplied) in variable \code{s.wt} and (2) adding variables \code{c.wt} and \code{n.wt} (holding principal score weights, ie probabilities of being a complier and of being a noncomplier, respectively). The reason this input is in this form is that this function will be incorporated in the internal function that estimates nuisance parameters, where this internal function (if required for targeted \eqn{\mu} estimation) will estimate these weights and put the data in this form for \eqn{\mu} functions estimation.
#' @param quants If TRUE, compute quantiles of the estimated conditional \eqn{Y|X,Z=0} distribution. Default is FALSE.
#' @param num.p.points Number of quantile points. Default is 100.
#' @return If \code{quants} is FALSE, a data frame with variable \code{id} and the mu functions. If \code{quants} is TRUE, a list including this data frame (named \code{mus}) and a matrix named \code{quants} containing p-quantiles for \eqn{X} values in each row of \code{data}.
#' @export


earn3_FUN2 <- function(data,
                      estimator,
                      targeted,
                      quants = FALSE,
                      num.p.points = 100) {

    dat <- data; rm(data)

    dat$c <- ifelse(dat$z==0, 0, dat$c)
    dat$n <- ifelse(dat$z==0, 0, 1 - dat$c)

    Y.form1 <-
        "work3 ~ age + sex + race + edu + marital + hh.kids + hh.income + econ.hard + occu + wks.unemp + part.motiv + seek.motiv + seek.effi + assertive + depress1"

    Y.form2 <-
        "y ~ age + sex + race + edu + marital + hh.kids + hh.income + econ.hard + occu + wks.unemp + part.motiv + seek.motiv + seek.effi + assertive + depress1"


    p.points <- (seq(1:num.p.points) / num.p.points) - (1 / (2*num.p.points))


    mus <- list()


    # convenient function for repeated action
    mu_convenient <- function(dat1,
                              dat2,
                              form1 = Y.form1,
                              form2 = Y.form2,
                              wts1,
                              wts2,
                              pred.dat = dat,
                              quants = FALSE,
                              p.vec = p.points) {

        # fit models
        suppressWarnings({
            mod1 <- do.call("glm",
                            list(formula = form1,
                                 data    = dat1,
                                 weights = wts1,
                                 family  = "quasibinomial"))

            mod2 <- do.call("glm",
                            list(formula = form2,
                                 data = dat2,
                                 weights = wts2,
                                 family = stats::quasi(link = "log",
                                                       variance = "mu^2")))
        })


        # compute the mu function
        pred1 <- predict(mod1, newdata = pred.dat, type = "response")
        pred2 <- predict(mod2, newdata = pred.dat, type = "response")
        pred2.calibrated <- pred2 *
            .wtd_mean(dat2$y,                           wts2) /
            .wtd_mean(predict(mod2, type = "response"), wts2)


        mu <- pred1 * pred2.calibrated

        if (!quants)  return(mu)



        # extract gamma parameters for quantile computing
        shape2 <- tryCatch(MASS::gamma.shape(mod2)$alpha,
                           error = function(e) return(NULL),
                           warning = function(w) return(NULL))

        suppressWarnings({
            if (is.null(shape2)) shape2 <- 1 / summary(mod2)$dispersion
        })

        scale2 <- pred2 / shape2


        # compute quantiles
        quants <- sapply(p.vec, function(p) {

            p2 <- (p + pred1 - 1) / pred1
            p2 <- pmax(0, p2)
            stats::qgamma(p = p2, shape = shape2, scale = scale2)
        })


        # calibrate quantiles to mu
        y.bound <- max(dat2$y) * 2
        m <- length(p.vec)

        quants <- t(sapply(1:nrow(quants), function(z) {

            q <- pmin(y.bound, quants[z,])
            q.mean <- mean(q)

            if (q.mean==0) {
                # for small mean values, the constant shape parameter may be not right
                # so quantiles at integration points may all be 0.
                # in this case, ad hoc fix all quantiles to the mean.
                # this happens only for small means, so the impact is negligible.
                q <- rep(mu[z], m)

            } else if (q.mean >= mu[z]) {
                # scale inward to mu
                q <- q * mu[z] / q.mean
            } else {
                # scale outward to mu while respecting y.bound
                q <- y.bound - (y.bound-q) * (y.bound-mu[z]) / (y.bound-q.mean)
            }

            q
        }))


        list(mu = mu,
             quants = quants)
    }




    # estimate mu.1c and mu.1n
    if (estimator != "emu") {

        # mu.1c
        dat1 <- dat[dat$c==1, ]
        dat2 <- dat[dat$c==1 & dat$work3==1, ]

        mus$mu.1c <- mu_convenient(dat1 = dat1,
                                   dat2 = dat2,
                                   wts1 = dat1$s.wt,
                                   wts2 = dat2$s.wt)
        rm(dat1, dat2)

        # mu.1n
        dat1 <- dat[dat$n==1, ]
        dat2 <- dat[dat$n==1 & dat$work3==1, ]

        mus$mu.1n <- mu_convenient(dat1 = dat1,
                                   dat2 = dat2,
                                   wts1 = dat1$s.wt,
                                   wts2 = dat2$s.wt)
        rm(dat1, dat2)
    }



    # estimate PI-based mu.0c and mu.0n
    # (both are mu.0 but if "targeted" estimated with different weights)
    # and
    # associated quantiles of Y | X,Z=0 for estimating A5-beta-based mu.0c and mu.0n

    dat1 <- dat[dat$z==0, ]
    dat2 <- dat[dat$z==0 & dat$work3==1, ]


    if (targeted || estimator=="MS") {

        # mu.0c and associated quantiles
        tmp <- mu_convenient(dat1 = dat1,
                             dat2 = dat2,
                             wts1 = dat1$s.wt*dat1$c.wt,
                             wts2 = dat2$s.wt*dat2$c.wt,
                             quants = quants)

        if (quants) { mus$mu.0c <- tmp$mu; y0c.quants <- tmp$quants
        } else      { mus$mu.0c <- tmp
        }
        rm(tmp)

        # mu.0n and associated quantiles
        tmp <- mu_convenient(dat1 = dat1,
                             dat2 = dat2,
                             wts1 = dat1$s.wt*dat1$n.wt,
                             wts2 = dat2$s.wt*dat2$n.wt,
                             quants = quants)

        if (quants) { mus$mu.0n <- tmp$mu; y0n.quants <- tmp$quants
        } else      { mus$mu.0n <- tmp
        }
        rm(tmp)

    } else {

        # both mu.0c and mu.0n
        tmp <- mu_convenient(dat1 = dat1,
                             dat2 = dat2,
                             wts1 = dat1$s.wt,
                             wts2 = dat2$s.wt,
                             quants = quants)

        if (quants) {
            mus$mu.0c <- mus$mu.0n <- tmp$mu
            y0c.quants <- y0n.quants <- tmp$quants
        } else      {
            mus$mu.0c <- mus$mu.0n <- tmp
        }
        rm(tmp)

    }

    rm(dat1, dat2)

    mus <- cbind(id = dat$id, data.frame(mus))

    if (!quants) return(mus)


    list(mus = mus,
         quants = list(y0c = y0c.quants,
                       y0n = y0n.quants))
}






#' @rdname earn3_FUN2
#' @order 2

earn3_FUN2_pimu1 <- function(data) {
    earn3_FUN2(data = data, estimator = "pimu", targeted = FALSE)
}


#' @rdname earn3_FUN2
#' @order 3

earn3_FUN2_pimu2 <- function(data) {
    earn3_FUN2(data = data, estimator = "pimu", targeted = TRUE)
}


#' @rdname earn3_FUN2
#' @order 4

earn3_FUN2_emu1 <- function(data) {
    earn3_FUN2(data = data, estimator = "emu", targeted = FALSE)
}


#' @rdname earn3_FUN2
#' @order 5

earn3_FUN2_emu2 <- function(data) {
    earn3_FUN2(data = data, estimator = "emu", targeted = TRUE)
}


#' @rdname earn3_FUN2
#' @order 6

earn3_FUN2_MS1 <- function(data) {
    earn3_FUN2(data = data, estimator = "MS", targeted = FALSE)
}


#' @rdname earn3_FUN2
#' @order 7

earn3_FUN2_MS2 <- function(data) {
    earn3_FUN2(data = data, estimator = "MS", targeted = TRUE)
}


#' @rdname earn3_FUN2
#' @order 8

earn3_FUN2_IF1 <- function(data) {
    earn3_FUN2(data = data, estimator = "IF", targeted = FALSE)
}


#' @rdname earn3_FUN2
#' @order 9

earn3_FUN2_IF2 <- function(data) {
    earn3_FUN2(data = data, estimator = "IF", targeted = TRUE)
}


#' @rdname earn3_FUN2
#' @order 10

earn3_FUN2_IFH1 <- function(data) {
    earn3_FUN2(data = data, estimator = "IFH", targeted = FALSE)
}


#' @rdname earn3_FUN2
#' @order 11

earn3_FUN2_IFH2 <- function(data) {
    earn3_FUN2(data = data, estimator = "IFH", targeted = TRUE)
}








#' @rdname earn3_FUN2
#' @order 12

earn3_FUN2.quants_pimu1 <- function(data) {
    earn3_FUN2(data = data, estimator = "pimu", targeted = FALSE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 13

earn3_FUN2.quants_pimu2 <- function(data) {
    earn3_FUN2(data = data, estimator = "pimu", targeted = TRUE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 14

earn3_FUN2.quants_emu1 <- function(data) {
    earn3_FUN2(data = data, estimator = "emu", targeted = FALSE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 15

earn3_FUN2.quants_emu2 <- function(data) {
    earn3_FUN2(data = data, estimator = "emu", targeted = TRUE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 16

earn3_FUN2.quants_MS1 <- function(data) {
    earn3_FUN2(data = data, estimator = "MS", targeted = FALSE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 17

earn3_FUN2.quants_MS2 <- function(data) {
    earn3_FUN2(data = data, estimator = "MS", targeted = TRUE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 18

earn3_FUN2.quants_IF1 <- function(data) {
    earn3_FUN2(data = data, estimator = "IF", targeted = FALSE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 19

earn3_FUN2.quants_IF2 <- function(data) {
    earn3_FUN2(data = data, estimator = "IF", targeted = TRUE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 20

earn3_FUN2.quants_IFH1 <- function(data) {
    earn3_FUN2(data = data, estimator = "IFH", targeted = FALSE, quants = TRUE)
}


#' @rdname earn3_FUN2
#' @order 21

earn3_FUN2.quants_IFH2 <- function(data) {
    earn3_FUN2(data = data, estimator = "IFH", targeted = TRUE, quants = TRUE)
}
