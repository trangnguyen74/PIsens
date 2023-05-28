
##############################
#### depress3_mus

#' JOBSII-specific: estimate \eqn{\mu} functions (and optionally, quantiles) for the outcome \code{depress3}
#'
#' Customized estimation of \eqn{\mu} functions and conditional \eqn{Y|X,Z=0} quantiles for the outcome \code{depress3} in the paper data example \code{jobsII}
#'
#' @details Important note: This function uses the \code{s.wt} variable in \code{data} in the estimation. This is required of any custom function for \eqn{\mu} estimation, even if this variable is all 1 (i.e., no sampling weights), because this variable will be used by the bootstrap to incorporate bootstrap weights.
#'
#' @order 1
#' @inheritParams .estimate_nuisance
#' @param data A data frame that has been processed by \code{declare_data()} and --- if needed for targeted \eqn{\mu} estimation --- has been enhanced by (1) having inverse probability of assigned treatment weight incorporated (ie multiplied) in variable \code{s.wt} and (2) adding variables \code{c.wt} and \code{n.wt} (holding principal score weights, ie probabilities of being a complier and of being a noncomplier, respectively). The reason this input is in this form is that this function will be incorporated in the internal function that estimates nuisance parameters, where this internal function (if required for targeted \eqn{\mu} estimation) will estimate these weights and put the data in this form for \eqn{\mu} functions estimation.
#' @param quants If TRUE, compute quantiles of the estimated conditional \eqn{Y|X,Z=0} distribution. Default is FALSE.
#' @param num.p.points Number of quantile points. Default is 100.
#' @return If \code{quants} is FALSE, a data frame with variable \code{id} and the mu functions. If \code{quants} is TRUE, a list including this data frame (named \code{mus}) and a matrix named \code{quants} containing p-quantiles for \eqn{X} values in each row of \code{data}.
#' @keywords internal


depress3_mus <- function(data,
                         estimator,
                         targeted,
                         quants = FALSE,
                         num.p.points = 100) {

    dat <- data; rm(data)

    dat$c <- ifelse(dat$z==0, 0, dat$c)
    dat$n <- ifelse(dat$z==0, 0, 1 - dat$c)

    dat$y <- (dat$y - 1) / 4

    Y.form <- "y ~ age + sex + race + edu + marital + hh.kids + hh.income + econ.hard + occu + wks.unemp + part.motiv + seek.motiv + seek.effi + assertive + depress1"


    p.points <- (seq(1:num.p.points) / num.p.points) - (1 / (2*num.p.points))


    mus <- list()


    # convenient function for repeated action
    mu_convenient <- function(form = Y.form,
                              data,
                              wts,
                              pred.dat = dat,
                              quants = FALSE,
                              p.vec = p.points) {

        ### estimate mu
        suppressWarnings({
            mean.mod <- do.call("glm",
                                list(formula = form,
                                     data    = data,
                                     weights = wts,
                                     family  = "quasibinomial"))
        })


        mu <- predict(mean.mod, newdata = pred.dat, type = "response")

        if (!quants)  return(mu)



        ### estimate quantiles
        # adjust outcome data in within the (0,1) interval (not including 0, 1)
        data$y <- (data$y * (nrow(data)-1) + 0.5) / nrow(data)
        mu.fit <- predict(mean.mod, type = "response")
        mu.fit <- (mu.fit * (nrow(data)-1) + 0.5) / nrow(data)
        mu.logit <- log(mu.fit / (1-mu.fit))


        suppressWarnings({
            dist.mod <- do.call(eval(parse(text = "betareg::betareg")),
                                list(formula = "y ~ -1",
                                     offset = mu.logit,
                                     link = "logit",
                                     data = data,
                                     weights = wts,
                                     type = "BR"))

            quants <- betareg::predict(dist.mod, newdata = pred.dat,
                                       type = "quantile",
                                       at = p.vec)
        })

        # calibrate quantiles to mu
        q.means <- rowMeans(quants)
        m <- length(p.vec)

        quants <- t(sapply(1:length(q.means), function(z) {

            if (q.means[z]==0 | q.means[z]==1) { q <- rep(mu[z], m)
            } else if (q.means[z] >= mu[z])    { q <- quants[z,] * mu[z] / q.means[z]
            } else                             { q <- 1 - (1-quants[z,]) * (1-mu[z]) / (1-q.means[z])
            }
            q
        }))

        list(mu = mu,
             quants = quants)
    }




    # estimate mu.1c and mu.1n
    if (estimator != "emu") {

        mus$mu.1c <- mu_convenient(data = dat[dat$c==1, ],
                                   wts  = dat$s.wt[dat$c==1])

        mus$mu.1n <- mu_convenient(data = dat[dat$n==1, ],
                                   wts  = dat$s.wt[dat$n==1])
    }



    # estimate PI-based mu.0c and mu.0n
    # (both are mu.0 but if "targeted" estimated with different weights)
    # and
    # associated quantiles of Y | X,Z=0 for estimating A5-beta-based mu.0c and mu.0n

    if (targeted || estimator=="MS") {

        # mu.0c and associated quantiles
        tmp <- mu_convenient(data   = dat[dat$z==0, ],
                             wts    = (dat$s.wt*dat$c.wt)[dat$z==0],
                             quants = quants)

        if (quants) { mus$mu.0c <- tmp$mu; y0c.quants <- tmp$quants
        } else      { mus$mu.0c <- tmp
        }
        rm(tmp)

        # mu.0n and associated quantiles
        tmp <- mu_convenient(data   = dat[dat$z==0, ],
                             wts    = (dat$s.wt*dat$n.wt)[dat$z==0],
                             quants = quants)

        if (quants) { mus$mu.0n <- tmp$mu; y0n.quants <- tmp$quants
        } else      { mus$mu.0n <- tmp
        }
        rm(tmp)

    } else {

        # both mu.0c and mu.0n
        tmp <- mu_convenient(data   = dat[dat$z==0, ],
                             wts    = dat$s.wt[dat$z==0],
                             quants = quants)

        if (quants) {
            mus$mu.0c <- mus$mu.0n <- tmp$mu
            y0c.quants <- y0n.quants <- tmp$quants
        } else      {
            mus$mu.0c <- mus$mu.0n <- tmp
        }
        rm(tmp)

    }

    mus <- data.frame(mus) * 4 + 1
    mus <- cbind(id = dat$id, mus)


    if (!quants) return(mus)


    quants <- list(y0c = y0c.quants * 4 + 1,
                   y0n = y0n.quants * 4 + 1)

    list(mus = mus,
         quants = quants)
}






#' @rdname depress3_mus
#' @order 2

depress3_mus_pimu1 <- function(data) {
    depress3_mus(data = data, estimator = "pimu", targeted = FALSE)
}


#' @rdname depress3_mus
#' @order 3

depress3_mus_pimu2 <- function(data) {
    depress3_mus(data = data, estimator = "pimu", targeted = TRUE)
}


#' @rdname depress3_mus
#' @order 4

depress3_mus_emu1 <- function(data) {
    depress3_mus(data = data, estimator = "emu", targeted = FALSE)
}


#' @rdname depress3_mus
#' @order 5

depress3_mus_emu2 <- function(data) {
    depress3_mus(data = data, estimator = "emu", targeted = TRUE)
}


#' @rdname depress3_mus
#' @order 6

depress3_mus_MS1 <- function(data) {
    depress3_mus(data = data, estimator = "MS", targeted = FALSE)
}


#' @rdname depress3_mus
#' @order 7

depress3_mus_MS2 <- function(data) {
    depress3_mus(data = data, estimator = "MS", targeted = TRUE)
}


#' @rdname depress3_mus
#' @order 8

depress3_mus_IF1 <- function(data) {
    depress3_mus(data = data, estimator = "IF", targeted = FALSE)
}


#' @rdname depress3_mus
#' @order 9

depress3_mus_IF2 <- function(data) {
    depress3_mus(data = data, estimator = "IF", targeted = TRUE)
}


#' @rdname depress3_mus
#' @order 10

depress3_mus_IFH1 <- function(data) {
    depress3_mus(data = data, estimator = "IFH", targeted = FALSE)
}


#' @rdname depress3_mus
#' @order 11

depress3_mus_IFH2 <- function(data) {
    depress3_mus(data = data, estimator = "IFH", targeted = TRUE)
}








#' @rdname depress3_mus
#' @order 12

depress3_mus.quants_pimu1 <- function(data) {
    depress3_mus(data = data, estimator = "pimu", targeted = FALSE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 13

depress3_mus.quants_pimu2 <- function(data) {
    depress3_mus(data = data, estimator = "pimu", targeted = TRUE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 14

depress3_mus.quants_emu1 <- function(data) {
    depress3_mus(data = data, estimator = "emu", targeted = FALSE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 15

depress3_mus.quants_emu2 <- function(data) {
    depress3_mus(data = data, estimator = "emu", targeted = TRUE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 16

depress3_mus.quants_MS1 <- function(data) {
    depress3_mus(data = data, estimator = "MS", targeted = FALSE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 17

depress3_mus.quants_MS2 <- function(data) {
    depress3_mus(data = data, estimator = "MS", targeted = TRUE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 18

depress3_mus.quants_IF1 <- function(data) {
    depress3_mus(data = data, estimator = "IF", targeted = FALSE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 19

depress3_mus.quants_IF2 <- function(data) {
    depress3_mus(data = data, estimator = "IF", targeted = TRUE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 20

depress3_mus.quants_IFH1 <- function(data) {
    depress3_mus(data = data, estimator = "IFH", targeted = FALSE, quants = TRUE)
}


#' @rdname depress3_mus
#' @order 21

depress3_mus.quants_IFH2 <- function(data) {
    depress3_mus(data = data, estimator = "IFH", targeted = TRUE, quants = TRUE)
}
