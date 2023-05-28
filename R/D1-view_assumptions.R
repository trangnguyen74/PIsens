####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### visualize_GOR ####

#' \loadmathjax
#' Figure 2 in the paper
#'
#' Plot the connection between \mjeqn{\mu_{00}(X)}{ASCII representation} and \mjeqn{\mu_{01}(X)}{ASCII representation} implied by a GOR assumption
#'
#' @param point.mu.0,point.c.wt A pair of \mjeqn{\mu_0(X),\pi_c(X)}{ASCII representation} values to be shown as the example point in the plot
#' @export

visualize_GOR <- function(point.mu.0 = .5, point.c.wt = .25) {


    rho <- c(1/4, 1/3, 1/2, 2/3, 3/2, 2, 3, 4)

    rho.labs <- paste("GOR =",
                      c("1/4", "1/3", "1/2", "2/3", "3/2", "2", "3", "4"))


    # background layer

    anchor.mu.0 <- seq(0,1, .1)
    anchor.c.wt <- .5
    anchor.n.wt <- .5


    b.dat <- sapply(1:length(rho), function(z) {

        r.c <- rho[z]
        r.n <- 1/r.c

        a.c <- (anchor.c.wt + anchor.mu.0) * (r.c-1) + 1
        a.n <- (anchor.n.wt + anchor.mu.0) * (r.n-1) + 1

        alpha.c <- sqrt(a.c^2 - 4*anchor.c.wt*anchor.mu.0*r.c*(r.c-1))
        alpha.n <- sqrt(a.n^2 - 4*anchor.c.wt*anchor.mu.0*r.n*(r.n-1))

        beta.c <- (a.c - alpha.c) / 2
        beta.n <- (a.n - alpha.n) / 2

        dat <- data.frame(mu.0c = beta.c / ((r.c-1)*anchor.c.wt),
                          mu.0n = beta.n / ((r.n-1)*anchor.n.wt),
                          anchor = anchor.mu.0)



        dat <- tidyr::pivot_longer(dat,
                                   cols = c("mu.0c", "mu.0n"),
                                   names_to = "type",
                                   values_to = "y")

        dat$x <- ifelse(dat$type=="mu.0c", 1, 0)

        dat$rho <- rho.labs[z]

        dat

    }, simplify = FALSE)

    b.dat <- do.call(rbind, b.dat)

    b.dat$rho <- factor(b.dat$rho, levels = rho.labs)


    b.legend.title <- "lines show GOR-implied connection between:"

    b.legend.labs <- list(expression(mu["01"](X)), expression(mu["00"](X)))



    x <- y <- anchor <- type <- NULL


    p <-
        ggplot2::ggplot(data = b.dat,
                        ggplot2::aes(x = x,
                                     y = y,
                                     group = factor(anchor))) +

        ggplot2::geom_vline(xintercept = c(0, 1),
                            color = "lightgray") +

        ggplot2::geom_line(color = "gray50",
                           size = .2) +

        ggplot2::geom_point(ggplot2::aes(color = type,
                                         shape = type,
                                         size  = type)) +

        ggplot2::scale_color_manual(name   = b.legend.title,
                                    values = rep("gray40", 2),
                                    labels = b.legend.labs,
                                    guide  = ggplot2::guide_legend(order = 1)) +
        ggplot2::scale_shape_manual(name   = b.legend.title,
                                    values = c(16, 1),
                                    labels = b.legend.labs,
                                    guide  = ggplot2::guide_legend(order = 1)) +
        ggplot2::scale_size_manual (name   = b.legend.title,
                                    values = c(2, 2),
                                    labels = b.legend.labs,
                                    guide  = ggplot2::guide_legend(order = 1)) +
        ggplot2::scale_y_continuous(breaks = seq(0,1,.25),
                                    labels = c("l", "", "", "", "h")) +
        ggplot2::scale_x_continuous(breaks = seq(0,1,.25),
                                    labels = c(0, "", "", "", 1)) +
        ggplot2::labs(x = "",
                      y = "") +
        ggplot2::theme_minimal() +
        ggplot2::facet_wrap(~ rho, ncol = 8, scales = "free_y") +
        ggplot2::theme(legend.position    = "bottom",
                       legend.box         = "vertical",
                       legend.box.spacing = ggplot2::unit(0, "cm"),
                       legend.margin      = ggplot2::margin(-4, 0, 0, 0),
                       panel.spacing      = ggplot2::unit(1, "cm"),
                       panel.grid.minor   = ggplot2::element_blank())



    # foreground layer

    mu.0 <- point.mu.0
    c.wt <- point.c.wt
    n.wt <- 1 - c.wt

    f.coords <- data.frame(x = rep(c.wt, 2),
                           y = rep(mu.0, 2),
                           xend = c(c.wt, 0),
                           yend = c(0, mu.0),
                           rho = rho.labs[1])

    f.coords$rho <- factor(f.coords$rho, levels = rho.labs)


    x <- y <- xend <- yend <- NULL

    p <- p +
        ggplot2::geom_segment(data = f.coords,
                              ggplot2::aes(x = x,
                                           y = y,
                                           xend = xend,
                                           yend = yend),
                              inherit.aes = FALSE,
                              color = "dodgerblue",
                              size = .2)



    f.coord.labs <- data.frame(x = c(-.05, c.wt),
                               y = c(mu.0, -.1),
                               hjust = c(1, .3),
                               vjust = c(.5, 0),
                               label = c("mu[0](X)", "pi[1](X)"),
                               rho = rho.labs[1])

    f.coord.labs$rho <- factor(f.coord.labs$rho, levels = rho.labs)


    x <- y <- label <- hjust <- vjust <- NULL

    p <- p +
        ggplot2::geom_text(data = f.coord.labs,
                           ggplot2::aes(x = x,
                                        y = y,
                                        label = label,
                                        hjust = hjust,
                                        vjust = vjust),
                           inherit.aes = FALSE,
                           size = 3,
                           color = "dodgerblue",
                           parse = TRUE) +
        ggplot2::coord_cartesian(xlim = c(0,1),
                                 ylim = c(0,1),
                                 clip = "off")




    f.dat <- sapply(1:length(rho), function(z) {

        r.c <- rho[z]
        r.n <- 1/r.c

        a.c <- (c.wt+mu.0)*(r.c-1)+1
        a.n <- (n.wt+mu.0)*(r.n-1)+1

        alpha.c <- sqrt(a.c^2 - 4*c.wt*mu.0*r.c*(r.c-1))
        alpha.n <- sqrt(a.n^2 - 4*n.wt*mu.0*r.n*(r.n-1))

        beta.c <- (a.c - alpha.c) / 2
        beta.n <- (a.n - alpha.n) / 2

        dat <- data.frame(type = c("mu.0n", "mu.0", "mu.0c"),
                          y = c(beta.n / ((r.n-1)*n.wt),
                                mu.0,
                                beta.c / ((r.c-1)*c.wt)),
                          x = c(0, c.wt, 1))

        dat$type = factor(dat$type, levels = c("mu.0",
                                               "mu.0c",
                                               "mu.0n"))

        dat$rho <- rho.labs[z]

        dat
    },
    simplify = FALSE)

    f.dat <- do.call(rbind, f.dat)

    f.dat$rho <- factor(f.dat$rho, levels = rho.labs)


    f.legend.title <- "example:"
    f.legend.labs <- list(expression(paste("a pair of ", mu[0](X), " and ", pi[1](X), " values")),
                          expression(paste("implied ", mu["01"](X))),
                          expression(paste("implied ", mu["00"](X))))




    x <- y <- type <- NULL

    p <- p +
        ggnewscale::new_scale_color() +
        ggnewscale::new_scale("shape") +
        ggnewscale::new_scale("size")

    p <- p +
        ggplot2::geom_line(data = f.dat,
                           ggplot2::aes(x = x,
                                        y = y),
                           inherit.aes = FALSE,
                           color = "gray50",
                           size = .3) +
        ggplot2::geom_point(data = f.dat,
                            ggplot2::aes(x = x,
                                         y = y,
                                         color = type,
                                         shape = type,
                                         size = type),
                            inherit.aes = FALSE) +
        ggplot2::scale_color_manual(name   = f.legend.title,
                                    values = c("blue", "magenta", "magenta"),
                                    labels = f.legend.labs,
                                    guide  = ggplot2::guide_legend(order = 2)) +
        ggplot2::scale_shape_manual(name   = f.legend.title,
                                    values = c(15, 16, 1),
                                    labels = f.legend.labs,
                                    guide  = ggplot2::guide_legend(order = 2)) +
        ggplot2::scale_size_manual (name   = f.legend.title,
                                    values = c(2.5, 2.5, 2.5),
                                    labels = f.legend.labs,
                                    guide  = ggplot2::guide_legend(order = 2))

    p

}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### visualize_kappa ####

#' Make the kappa plot in the paper
#'
#' @param pi.x A vector of values for \eqn{\pi_c(X)} to be plotted.
#' @export

visualize_kappa <- function(pi.x = c(.15, .3, .5)) {

    p <- seq(0,1, .001)
    kappa <- c(0.005, .01, .02, .05, .1, .5)

    kdat <- sapply(pi.x, function(pi) {

        tmp <- sapply(kappa, function(k) {
            alpha <- pi     * (1-k)/k
            beta  <- (1-pi) * (1-k)/k

            data.frame(pi.x = pi,
                       kappa = k,
                       pi.xy0 = p,
                       density = stats::dbeta(x = p,
                                              shape1 = alpha,
                                              shape2 = beta))

        },
        simplify = FALSE)

        do.call(rbind, tmp)

    }, simplify = FALSE)

    kdat <- do.call(rbind, kdat)

    kdat <- rbind(kdat,
                  data.frame(pi.x = pi.x,
                             kappa = 0,
                             pi.xy0 = pi.x,
                             density = 0),
                  data.frame(pi.x = pi.x,
                             kappa = 0,
                             pi.xy0 = pi.x,
                             density = max(kdat$density)))

    kdat$pi.x <- factor(kdat$pi.x,
                        levels = pi.x,
                        labels = paste0("pi[c](X)~'='~", as.character(pi.x)))

    # kappa.colors <- RColorBrewer::brewer.pal(n = 8, "Reds")[8:3] #there are 9, I exluded the two lighter hues

    kappa.colors <- c("black", "#99000D", "#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1")

    mdat <- data.frame(pi.x = pi.x,
                       xintercept = pi.x)

    pi.xy0 <- density <- NULL

    ggplot2::ggplot(data = kdat,
                    ggplot2::aes(x = pi.xy0,
                                 y = density,
                                 group = factor(kappa),
                                 color = factor(kappa))) +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(name = expression(kappa~'values'),
                                    values = c("black", kappa.colors)) +
        ggplot2::theme_minimal() +
        ggplot2::facet_grid(. ~ pi.x,
                            labeller = ggplot2::label_parsed) +
        ggplot2::labs(x = expression(tilde(pi)[c](X,Y[0]))) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.spacing    = ggplot2::unit(1, "cm"))

}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### visualize_MR_implication ####

#' \loadmathjax
#' Plot the implication of sensitivity MR values (violin plot)
#'
#' Plot the distributions of \mjeqn{\mu_{00}(X)}{ASCII representation} and of \mjeqn{\mu_{01}(X)}{ASCII representation} that are implied by different values of the sensitivity MR parameter. This plot uses side-by-side violin plots.
#'
#' @inheritParams PIsens
#' @param MRs Sensitivity MR values to be used for the plot.
#' @param MRlabs Labels for MRs to put on the plot. For example, for an MR value of 1/3, one may want to use the label "1/3".
#' @export

# TODO: fix visualize_MR_implication

visualize_MR_implication <- function(MRs,
                                     MRlabs,

                                     data,
                                     estimator,
                                     targeted,

                                     Z.form,
                                     C.form,

                                     Y.form,
                                     Y.model,
                                     Y.bounds,
                                     Y.FUN = NULL) {

    # the nuisance estimation here is borrowed from PIsens
    key.inputs <- list(estimator = estimator,
                       targeted  = targeted)

    nuis.inputs <- mget(c("Z.form", "C.form",
                          "Y.form", "Y.model", "Y.bounds",
                          "Y.FUN"),
                        ifnotfound = list(NULL, NULL,
                                          NULL, NULL, NULL, NULL))

    nuis.inputs <- do.call(".clean_nuis_inputs",
                           args = c(key.inputs,
                                    sensitivity = TRUE,
                                    nuis.inputs))

    nuis <- do.call(".estimate_nuisance",
                    args = c(data = list(data),
                             key.inputs,
                             nuis.inputs))

    mu.dat <- lapply(1:length(MRs), function(z) {

        r <- MRs[z]

        gamma.n <- 1 / ((r-1) * nuis$c.wt + 1)
        gamma.c <- gamma.n * r

        rbind(data.frame(MR = r,
                         MRlab = MRlabs[z],
                         stratum = "compliers",
                         mean.Y0 = nuis$mu.0c * gamma.c,
                         weight = nuis$s.wt * nuis$c.wt),
              data.frame(MR = r,
                         MRlab = MRlabs[z],
                         stratum = "noncompliers",
                         mean.Y0 = nuis$mu.0n * gamma.n,
                         weight = nuis$s.wt * nuis$n.wt))
    })

    mu.dat <- do.call(rbind, mu.dat)

    mu.dat$MRlab <- factor(mu.dat$MRlab, levels = MRlabs)
    mu.dat$stratum <- factor(mu.dat$stratum)

    mu.dat$type <- ifelse(mu.dat$MR==1, "main", "sens")


    MRlab <- mean.Y0 <- weight <- stratum <- type <- NULL

    p <- ggplot2::ggplot(mu.dat,
                         ggplot2::aes(x = MRlab,
                                      y = mean.Y0,
                                      weight = weight,
                                      fill = stratum,
                                      color = type)) +
        .geom_split_violin(nudge = .02) +
        ggplot2::scale_color_manual(values = c("black", "darkgray"),
                                    guide = "none") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "sensitivity MR",
                      y = "") +
        ggplot2::theme(legend.title = ggplot2::element_blank())

    suppressWarnings(print(p))

}





