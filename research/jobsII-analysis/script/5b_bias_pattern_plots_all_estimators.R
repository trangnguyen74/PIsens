


..plot_finite.sample.bias_light <- function(sens,
                                            sens.type,
                                            include.bias.correction = FALSE,
                                            effect.ylims = NULL,
                                            outcome.ylims = NULL,
                                            sens.para.lims = NULL) {

    sens <- sapply(c("cace", "nace", "tau.0c", "tau.0n"), function(z) {

        cbind(estimand = z,
              data.frame(sens[[z]]))
    },
    simplify = FALSE)

    sens <- do.call(rbind, sens)

    if (!is.null(sens.para.lims)) {
        sens <- sens[sens$sens.para <= max(sens.para.lims), ]
        sens <- sens[sens$sens.para >= min(sens.para.lims), ]
    }

    all.elements <- c("point", "boot.mean", "boot2.mean", "point.bc1", "point.bc2")

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

    sens$key <- factor(sens$key,
                       levels = all.elements,
                       labels = c("point", "BM1", "BM2", "BC1", "BC2"))

    colors <- c("black", "red3", "darkorange", "blue", "dodgerblue")
    linetypes <- c(1, 2, 4, 2, 4)


    if (sens.type %in% c("OR", "GOR", "MR")) {

        x.ticks <- unique(round(sens$sens.para.coord))
        x.tick.labels <- ifelse(x.ticks==0, as.character(1),
                                ifelse(x.ticks>0, as.character(x.ticks+1),
                                       paste0("1/", as.character(abs(x.ticks)+1))))
        x.label <- paste("sensitivity", sens.type)


    } else if (sens.type=="SMDe") {
        x.ticks <- unique(round(sens$sens.para.coord * 2)/2)
        x.tick.labels <- x.ticks
        x.label <- "sensitivity SMD"

    } else if (sens.type=="beta.quant") {

        x.ticks <- unique(round(sens$sens.para.coord, 2))
        x.tick.labels <- x.ticks
    }


    p <- ggplot2::ggplot(data = sens,
                         ggplot2::aes(x = sens.para.coord,
                                      y = value,
                                      color = key,
                                      linetype = key)) +
        ggplot2::geom_vline(xintercept = 0, color = "gray") +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(name = "", values = colors) +
        ggplot2::scale_linetype_manual(name = "", values = linetypes) +
        ggplot2::scale_x_continuous(breaks = x.ticks,
                                    labels = x.tick.labels) +
        ggplot2::facet_grid(type ~ stratum,
                            scales = "free_y") +
        ggplot2::labs(x = x.label,
                      y = "") +
        ggplot2::theme_bw()




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







combine_plots <- function(outcome,
                          sens.type,
                          include.bias.correction,
                          include.CI) {

    sens.file.name <- paste0(outcome, "_", sens.type, ".rds")

    sens <- readRDS(here::here("research", "jobsII-analysis", "outputs", "results", "all-estimators", "sens", sens.file.name))

    sens.names <- names(sens)


    p.list <- list()

    for (z in 1:length(sens.names)) {

        p <- ..plot_finite.sample.bias_light(sens[[z]]$sens,
                                             sens.type = sens.type,
                                             include.bias.correction = include.bias.correction) +
            ggplot2::ggtitle(sens.names[[z]])

        if (z==1) p.legend <- cowplot::get_legend(p)

        p.list[[z]] <- p + ggplot2::theme(legend.position = "none")
    }

    if (outcome=="earn") {
        p.list <- c(p.list[1:6], p.legend = list(p.legend), p.list[7:11])
    } else {
        p.list <- c(p.list[1:4], p.legend = list(p.legend), p.list[5:9])
    }


    cowplot::plot_grid(plotlist = p.list, nrow = 2, byrow = FALSE)
}




sens.dir    <- here::here("research", "jobsII-analysis", "outputs", "results", "all-estimators", "sens")
figs.dir <- here::here("research", "jobsII-analysis", "outputs", "figures")



pdf(file.path(figs.dir, "fig13.pdf"), width = 14, height = 8)
combine_plots(outcome = "work", sens.type = "OR", include.bias.correction = TRUE)
dev.off()

pdf(file.path(figs.dir, "fig14.pdf"), width = 14, height = 8)
combine_plots(outcome = "depress", sens.type = "GOR", include.bias.correction = TRUE)
dev.off()

pdf(file.path(figs.dir, "fig15.pdf"), width = 14, height = 8)
combine_plots(outcome = "depress", sens.type = "SMDe", include.bias.correction = TRUE)
dev.off()

pdf(file.path(figs.dir, "fig16.pdf"), width = 14, height = 8)
combine_plots(outcome = "earn", sens.type = "MR", include.bias.correction = TRUE)
dev.off()
























