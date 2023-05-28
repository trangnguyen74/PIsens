####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .plot_sens_multiple ####

#' (internal) Uber plot of the wholee bunch of estimators
#'
#' @param sens A named list of results from different sensitivity analyses from running \code{PIsens()}. The names can be \code{epi1}, \code{epi2}, \code{pimu1}, \code{pimu2}, \code{emu1}, \code{emu2}, \code{MS1}, \code{MS2}, \code{IF1}, \code{IF2}, \code{IFH1}, \code{IFH2} (where the letter part indicates the estimator, and the number part indicates whether the estimator is based on (1) untargeted or (2) targeted nuisance estimation).
#' @param to.plot The parameter for which results are to be plotted. Options are: "cace", "nace", "tau.0c", "tau.0n".
#' @keywords internal


.plot_sens_multiple <- function(sens, to.plot) {

    sens.para.coord <- value <- key <- bias.corrected <- size <- point <- NULL

    sens <- lapply(sens, function(z) z[[to.plot]])

    pdat <- sapply(1:length(sens), function(z) {
        tmp <- data.frame(sens[[z]])

        tmp <- tidyr::pivot_longer(tmp,
                                   c("point", "point.bc1", "point.bc2",
                                     "bca.low", "bca.high",
                                     "perc.low", "perc.high",
                                     "boot.mean", "boot2.mean"),
                                   names_to = "key", values_to = "value")

        tmp$type <- tmp$key
        tmp$type <- ifelse(tmp$type %in% c("bca.low", "bca.high"), "BCa CI", tmp$type)
        tmp$type <- ifelse(tmp$type %in% c("perc.low", "perc.high"), "percentile CI", tmp$type)

        tmp$point <- ifelse(tmp$key %in% c("point", "point.bc1", "point.bc2", "boot.mean", "boot2.mean"),
                            TRUE, FALSE)
        tmp$point <- factor(tmp$point, labels = c("confidence interval", "point estimate"))
        tmp$bias.corrected <- ifelse(tmp$key %in% c("point.bc1", "point.bc2", "bca.low", "bca.high"),
                                     TRUE, FALSE)

        tmp$size <- 4*(tmp$key %in% c("point.bc1", "point.bc2", "bca.low", "bca.high")) +
            3*(tmp$key %in% c("point", "perc.low", "perc.high")) +
            2*(tmp$key=="boot.mean") +
            1*(tmp$key=="boot2.mean")

        estimator <- names(sens)[z]

        if        (grepl("epi",  estimator)) { tmp$estimator <- "epi"
        } else if (grepl("pimu", estimator)) { tmp$estimator <- "pimu"
        } else if (grepl("emu",  estimator)) { tmp$estimator <- "emu"
        } else if (grepl("MS",   estimator)) { tmp$estimator <- "MS"
        } else if (grepl("IFH",  estimator)) { tmp$estimator <- "IFH"
        } else if (grepl("IF",   estimator)) { tmp$estimator <- "IF"
        }

        tmp$estimator <- factor(tmp$estimator, levels = c("epi", "pimu", "emu",
                                                          "MS", "IF", "IFH"))

        if (grepl("2", estimator)) { tmp$version <- "targeted"
        } else                     { tmp$version <- "untargeted"
        }

        tmp
    },
    simplify = FALSE)

    pdat <- do.call(rbind, pdat)

    x.ticks <- unique(pdat$sens.para.coord)
    x.ticks <- x.ticks[x.ticks==floor(x.ticks)]
    x.labels <- ifelse(x.ticks==0, as.character(1),
                       ifelse(x.ticks>0, as.character(x.ticks+1),
                              paste("1/", as.character(abs(x.ticks)+1))))


    ggplot2::ggplot(pdat,
                    ggplot2::aes(x = sens.para.coord, y = value,
                                 group = key,
                                 color = bias.corrected,
                                 size = size,
                                 linetype = point)) +
        ggplot2::geom_line() +
        ggplot2::scale_size_continuous(range = c(.25, 1)) +
        ggplot2::scale_linetype_manual(values = c(2, 1)) +
        ggplot2::scale_x_continuous(breaks = x.ticks, labels = x.labels) +
        ggplot2::guides(size = "none", linetype = "none", color = "none") +
        ggplot2::facet_grid(version ~ estimator) +
        ggplot2::theme_bw() +
        ggplot2::geom_vline(xintercept = 0, color = "gray") +
        ggplot2::geom_hline(yintercept = 0, color = "gray") +
        ggplot2::labs(title = to.plot,
                      x = "rho", y = "")


}




