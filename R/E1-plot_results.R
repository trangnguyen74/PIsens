####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### plot_sens ####

#' Plot the result of a sensitivity analysis
#'
#' @param sens Result from function \code{PIsens}.
#' @param CI Type of confidence interval to be plotted. Options are "BCa" and "percentile". If want "BCa", \code{sens} must have been the result of running \code{PIsens()} with \code{BCa.CI = TRUE} (which is the default of \code{PIsens()}).
#' @param point Type of point estimates to be plotted. Options are "uncorrected", "corrected1" (meaning bias-corrected based on the bootstrap) and "corrected2 (meaning bias-corrected based on the double bootstrap). If want "corrected2", \code{sens} must have been the result of running \code{PIsens()} with \code{double.boot = TRUE} (which is the default of \code{PIsens()}).
#' @param effect.ylims Optional. A two-element numeric vector holding the plot limits for the NACE and CACE.
#' @param outcome.ylims Optional. A two-element numeric vector holding the plot limits for the potential outcome means.
#' @param sens.para.lims Optional. A two-element numeric vector holding the plot limits for the sensitivity parameter. To be used to make publishable plots if analysis was run on a larger range than should be shown, e.g., if the relevant range is narrower.
#' @export


plot_sens <- function(sens,
                      CI = "BCa",
                      point = "corrected1",
                      effect.ylims = NULL,
                      outcome.ylims = NULL,
                      sens.para.lims = NULL) {


    sens.para.coord <- value <- estimand <- estimate <- x <- y <- NULL


    error1 <- "Only CI types BCa and percentile are available."
    error2 <- "Only three options available for `point`: uncorrected, corrected1, corrected2"



    if (CI=="BCa")               { ci.names <- c("bca.low", "bca.high")
    } else if (CI=="percentile") { ci.names <- c("perc.low", "perc.high")
    } else { stop(error1)
    }

    if (point=="uncorrected")       { point.name <- "point"
    } else if (point=="corrected1") { point.name <- "point.bc1"
    } else if (point=="corrected2") { point.name <- "point.bc2"
    } else { stop(error2)
    }

    sens.type <- sens$sens.specs$sens.type
    sens      <- sens$sens

    sens <- sapply(names(sens), function(z) {

        cbind(estimand = z,
              data.frame(sens[[z]]))
    },
    simplify = FALSE)

    sens <- do.call(rbind, sens)

    if (!is.null(sens.para.lims)) {
        sens <- sens[sens$sens.para <= max(sens.para.lims), ]
        sens <- sens[sens$sens.para >= min(sens.para.lims), ]
    }

    sens <- sens[names(sens) != "sens.para"]



    sens <- sens[, c("estimand", "sens.para.coord",
                     point.name, ci.names)]
    names(sens)[3:5] <- c("point", "low", "high")

    sens <- tidyr::pivot_longer(sens,
                                c("point", "low", "high"),
                                names_to = "estimate",
                                values_to = "value")

    sens$stratum <- ifelse(sens$estimand %in% c("cace", "tau.1c", "tau.0c"),
                           "compliers", "noncompliers")

    sens$estimand <- ifelse(sens$estimand %in% c("cace", "nace"),
                            "CACE/NACE",
                            ifelse(sens$estimand %in% c("tau.1c", "tau.1n"),
                                   "mean Y1", "mean Y0"))

    sens$type <- sens$estimand

    sens$type <- ifelse(sens$estimand=="CACE/NACE",
                        "average causal effect",
                        "potential outcome means")

    main <- sens[sens$sens.para.coord==0 & sens$estimate=="point", ]

    hline.dat <- data.frame(stratum = c("compliers", "noncompliers"),
                            type = "average causal effect",
                            threshold = 0)

    if (sens.type %in% c("OR", "GOR", "MR")) {

        x.ticks <- unique(round(sens$sens.para.coord))
        x.labels <- ifelse(x.ticks==0, as.character(1),
                           ifelse(x.ticks>0, as.character(x.ticks+1),
                                  paste0("1/", as.character(abs(x.ticks)+1))))

    } else if (sens.type=="SMDe") {
        x.ticks <- unique(round(sens$sens.para.coord * 2)/2)
        x.labels <- x.ticks

    } else if (sens.type=="beta.quant") {

        x.ticks <- unique(round(sens$sens.para.coord, 2))
        x.labels <- x.ticks
    }


    p <-
        ggplot2::ggplot(sens,
                        ggplot2::aes(x = sens.para.coord,
                                     y = value,
                                     group = interaction(estimand, estimate),
                                     color = estimand,
                                     linetype = estimate)) +
        ggplot2::geom_vline(xintercept = 0, color = "gray") +
        ggplot2::geom_hline(data = hline.dat,
                            ggplot2::aes(yintercept = 0), color = "gray") +
        ggplot2::geom_line() +
        ggplot2::geom_point(data = main) +
        ggplot2::ylab("") +
        ggplot2::scale_x_continuous(breaks = x.ticks, labels = x.labels) +
        ggplot2::scale_linetype_manual(values = c(2, 2, 1)) +
        ggplot2::scale_color_manual(name = "", values = c("black", "blue", "red")) +
        ggplot2::guides(linetype = "none") +
        ggplot2::facet_grid(type ~ factor(stratum, levels = c("noncompliers", "compliers")),
                            scales = "free_y") +
        ggplot2::theme_bw()


    if (sens.type %in% c("OR", "GOR", "MR")) {

        p <- p + ggplot2::xlab(paste("sensitivity", sens.type))

    } else if (sens.type=="SMDe") {

        p <- p + ggplot2::xlab("sensitivity SMD")

    } else if (sens.type=="beta.quant") {

        p <- p +
            ggplot2::xlab(expression(paste("signed ", kappa))) +
            ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(angle = 90,
                                                               vjust = 0.5,
                                                               hjust=1))
    }


    ylim.dat <- NULL

    if (!is.null(effect.ylims))
        ylim.dat <- rbind(ylim.dat,
                          data.frame(stratum = "compliers",
                                     type = "average causal effect",
                                     x = rep(0,2),
                                     y = effect.ylims))

    if (!is.null(outcome.ylims))
        ylim.dat <- rbind(ylim.dat,
                          data.frame(stratum = "compliers",
                                     type = "potential outcome means",
                                     x = rep(0,2),
                                     y = outcome.ylims))

    if (!is.null(ylim.dat))
        p <- p +
        ggplot2::geom_blank(data = ylim.dat,
                            ggplot2::aes(x = x, y = y), inherit.aes = FALSE)

    p
}









####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### plot_finite.sample.bias ####

#' Function that makes figure 5 in the paper
#'
#' @param sens The output of \code{PIsens()}.
#' @param include.bias.correction If TRUE, include the curves for bias-corrected point estimates based on the bootstrap and based on the double bootstrap. Default is FALSE.
#' @param include.CI If TRUE, include the curves for the point-wise 95% BCa confidence intervals (for a sense of the scale of the finite sample bias). Default is FALSE.
#' @param effect.ylims An optional two-element vector for the y-axis limits of the average effect panel.
#' @param outcome.ylims An optional two-element vector for the y-axis limits of the mean Y0 panel.
#' @param sens.para.lims An optional two-element vector for the range of the sensitivity parameter to be plotted. Use this if want to limit the range more than is available in \code{sens}.
#' @export


plot_finite.sample.bias <- function(sens,
                                    sens.type,
                                    include.bias.correction = FALSE,
                                    include.CI = FALSE,
                                    effect.ylims = NULL,
                                    outcome.ylims = NULL,
                                    sens.para.lims = NULL) {


    sens <- sapply(names(sens), function(z) {

        cbind(estimand = z,
              data.frame(sens[[z]]))
    },
    simplify = FALSE)

    sens <- do.call(rbind, sens)

    if (!is.null(sens.para.lims)) {
        sens <- sens[sens$sens.para <= max(sens.para.lims), ]
        sens <- sens[sens$sens.para >= min(sens.para.lims), ]
    }

    base.elements <- c("point", "boot.mean", "boot2.mean")
    bc.elements   <- c("point.bc1", "point.bc2")
    ci.elements   <- c("bca.low", "bca.high")
    all.elements  <- c(base.elements, bc.elements, ci.elements)

    sens <- sens[names(sens) %in% c("estimand", "sens.para.coord",
                                    all.elements)]

    sens <- tidyr::pivot_longer(sens,
                                all.elements,
                                names_to = "key", values_to = "value")

    sens$stratum <- ifelse(sens$estimand %in% c("cace", "tau.0c"),
                           "compliers", "noncompliers")
    sens$stratum <- factor(sens$stratum, levels = c("noncompliers", "compliers"))

    sens$type <- ifelse(sens$estimand %in% c("cace", "nace"),
                        "average effect", "mean Y0")

    if (sens.type %in% c("OR", "GOR", "MR")) {

        x.ticks <- unique(round(sens$sens.para.coord))
        x.labels <- ifelse(x.ticks==0, as.character(1),
                           ifelse(x.ticks>0, as.character(x.ticks+1),
                                  paste0("1/", as.character(abs(x.ticks)+1))))

    } else if (sens.type=="SMDe") {
        x.ticks <- unique(round(sens$sens.para.coord * 2)/2)
        x.labels <- x.ticks

    } else if (sens.type=="beta.quant") {

        x.ticks <- unique(round(sens$sens.para.coord, 2))
        x.labels <- x.ticks
    }


    # x.ticks <- unique(sens$sens.para.coord)
    # x.ticks <- x.ticks[x.ticks==floor(x.ticks)]
    # x.labels <- ifelse(x.ticks==0, as.character(1),
    #                    ifelse(x.ticks>0, as.character(x.ticks+1),
    #                           paste("1/", as.character(abs(x.ticks)+1))))

    base.dat <- sens[sens$key %in% base.elements, ]
    bc.dat   <- sens[sens$key %in% bc.elements  , ]
    ci.dat   <- sens[sens$key %in% ci.elements  , ]

    base.dat$base.key <- factor(base.dat$key,
                                levels = base.elements,
                                labels = c("point estimate (uncorrected)",
                                           "mean bootstrap estimate",
                                           "mean double bootstrap estimate"))
    bc.dat$bc.key     <- factor(bc.dat$key,
                                levels = bc.elements,
                                labels = c("based on bootstrap",
                                           "based on double bootstrap"))
    ci.dat$ci.key     <- "95% BCa confidence intervals"

    base.color <- c("black", "red3", "darkorange")
    bc.color   <- c("blue", "dodgerblue")
    ci.color   <- "gray"

    base.linetype <- c(1, 2, 4)
    bc.linetype   <- c(2,4)
    ci.linetype   <- 1



    legend.counter <- 1

    sens.para.coord <- value <- key <- base.key <- NULL

    p <- ggplot2::ggplot(data = base.dat,
                         ggplot2::aes(x = sens.para.coord,
                                      y = value,
                                      group = key,
                                      color = base.key,
                                      linetype = base.key)) +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(
            name = "Key comparison",
            values = base.color,
            guide  = ggplot2::guide_legend(order = legend.counter)) +
        ggplot2::scale_linetype_manual(
            name = "Key comparison",
            values = base.linetype,
            guide  = ggplot2::guide_legend(order = legend.counter)) +
        ggplot2::scale_x_continuous(breaks = x.ticks,
                                    labels = x.labels) +
        ggplot2::facet_grid(type ~ stratum,
                            scales = "free_y") +
        ggplot2::labs(x = "sensitivity parameter",
                      y = "") +
        ggplot2::theme_bw()


    if (include.bias.correction) {

        legend.counter <- legend.counter + 1

        p <- p +
            ggnewscale::new_scale_color() +
            ggnewscale::new_scale("linetype")

        bc.key <- NULL

        p <- p +
            ggplot2::geom_line(data = bc.dat,
                               inherit.aes = FALSE,
                               ggplot2::aes(x = sens.para.coord,
                                            y = value,
                                            group = key,
                                            color = bc.key,
                                            linetype = bc.key)) +
            ggplot2::scale_color_manual(
                name = "Bias-corrected point estimate",
                values = bc.color,
                guide  = ggplot2::guide_legend(order = legend.counter)) +
            ggplot2::scale_linetype_manual(
                name = "Bias-corrected point estimate",
                values = bc.linetype,
                guide  = ggplot2::guide_legend(order = legend.counter))
    }


    if (include.CI) {

        legend.counter <- legend.counter + 1

        p <- p +
            ggnewscale::new_scale_color() +
            ggnewscale::new_scale("linetype")

        ci.key <- NULL

        p <- p +
            ggplot2::geom_line(data = ci.dat,
                               inherit.aes = FALSE,
                               ggplot2::aes(x = sens.para.coord,
                                            y = value,
                                            group = key,
                                            color = ci.key,
                                            linetype = ci.key)) +
            ggplot2::scale_color_manual(
                name = "For a sense of scale",
                values = ci.color,
                guide  = ggplot2::guide_legend(order = legend.counter)) +
            ggplot2::scale_linetype_manual(
                name = "For a sense of scale",
                values = ci.linetype,
                guide  = ggplot2::guide_legend(order = legend.counter))


    }





    ylim.dat <- NULL

    if (!is.null(effect.ylims))
        ylim.dat <- rbind(ylim.dat,
                          data.frame(stratum = "compliers",
                                     type = "average effect",
                                     x = rep(0,2),
                                     y = effect.ylims))

    if (!is.null(outcome.ylims))
        ylim.dat <- rbind(ylim.dat,
                          data.frame(stratum = "compliers",
                                     type = "mean Y0",
                                     x = rep(0,2),
                                     y = outcome.ylims))

    if (!is.null(ylim.dat)) {
        ylim.dat$stratum <- factor(ylim.dat$stratum,
                                   levels = c("noncompliers", "compliers"))

        x <- y <- NULL

        p <- p +
            ggplot2::geom_blank(data = ylim.dat,
                                ggplot2::aes(x = x, y = y),
                                inherit.aes = FALSE)
    }


    p
}




