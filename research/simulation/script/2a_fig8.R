


..plot_for_fig8 <- function(outcome,
                            sens.type,
                            effect.ylims = NULL,
                            outcome.ylims = NULL) {

    sim.dir <- here::here("research", "simulation", "outputs", "simresults")

    if (outcome=="work") sim.dir <- file.path(sim.dir, "work", "OR")
    if (outcome=="earn") sim.dir <- file.path(sim.dir, "earn", "MR")
    if (outcome=="depress" && sens.type=="GOR")  sim.dir <- file.path(sim.dir, "depress", "GOR")
    if (outcome=="depress" && sens.type=="SMDe") sim.dir <- file.path(sim.dir, "depress", "SMDe")

    estimands <- c("cace", "nace", "tau.0c", "tau.0n")

    tmp <- readRDS(file.path(sim.dir, "true.rds"))
    true <- cbind(tmp$sens.params, tmp$values[,estimands])

    tmp <- readRDS(file.path(sim.dir, "sim.rds"))
    raw <- cbind(tmp$sens.params, apply(tmp$results$raw[,estimands,], c(1,2), function(z) mean(z, na.rm = TRUE)))
    bm1 <- cbind(tmp$sens.params, apply(tmp$results$bm1[,estimands,], c(1,2), function(z) mean(z, na.rm = TRUE)))
    bm2 <- cbind(tmp$sens.params, apply(tmp$results$bm2[,estimands,], c(1,2), function(z) mean(z, na.rm = TRUE)))
    bc1 <- cbind(tmp$sens.params, apply(tmp$results$bc1[,estimands,], c(1,2), function(z) mean(z, na.rm = TRUE)))
    bc2 <- cbind(tmp$sens.params, apply(tmp$results$bc2[,estimands,], c(1,2), function(z) mean(z, na.rm = TRUE)))


    dat <- rbind(data.frame(true, key = "true"),
                 data.frame(raw,  key = "raw"),
                 data.frame(bm1,  key = "bm1"),
                 data.frame(bm2,  key = "bm2"),
                 data.frame(bc1,  key = "bc1"),
                 data.frame(bc2,  key = "bc2"))

    dat <- dat[names(dat) != "sense.para"]

    dat <- tidyr::pivot_longer(dat,
                               cols = all_of(estimands),
                               values_to = "estimate",
                               names_to  = "estimand")

    dat$stratum <- ifelse(dat$estimand %in% c("cace", "tau.0c"),
                          "compliers", "noncompliers")
    dat$stratum <- factor(dat$stratum, levels = c("noncompliers", "compliers"))

    dat$type <- ifelse(dat$estimand %in% c("cace", "nace"),
                       "average effect", "mean Y0")

    dat$key <- factor(dat$key,
                      levels = c("true", "raw", "bm1", "bm2", "bc1", "bc2"),
                      labels = c("true value",
                                 "uncorrected point estimate",
                                 "mean bootstrap estimate",
                                 "mean double bootstrap estimate",
                                 "BC1 point estimate",
                                 "BC2 point estimate"))

    if (sens.type %in% c("OR", "GOR", "MR")) {

        x.ticks <- unique(round(dat$sens.para.coord))
        x.labels <- ifelse(x.ticks==0, as.character(1),
                           ifelse(x.ticks>0, as.character(x.ticks+1),
                                  paste0("1/", as.character(abs(x.ticks)+1))))

    } else if (sens.type=="SMDe") {
        x.ticks <- unique(round(dat$sens.para.coord * 2)/2)
        x.labels <- x.ticks

    } else if (sens.type=="beta.quant") {

        x.ticks <- unique(round(dat$sens.para.coord, 2))
        x.labels <- x.ticks
    }


    colors <- c("purple", "black", "red3", "darkorange", "blue", "dodgerblue")
    linetypes <- c(1, 1, 2, 4, 2,4)



    p <- ggplot2::ggplot(dat,
                         mapping = ggplot2::aes(x = sens.para.coord,
                                                y = estimate,
                                                color = key,
                                                linetype = key)) +
        ggplot2::geom_vline(xintercept = 0, color = "gray") +
        ggplot2::facet_grid(type ~ stratum,
                            scales = "free_y") +
        ggplot2::scale_x_continuous(breaks = x.ticks,
                                    labels = x.labels) +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(name = "", values = colors) +
        ggplot2::scale_linetype_manual(name = "", values = linetypes) +
        ggplot2::theme_bw() +
        ggplot2::labs(y = "")


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

    if (outcome=="work")
        p <- p + ggplot2::labs(x = "sensitivity OR",
                               title = "work for pay simulation:\nOR-based sensitivity analysis")

    if (outcome=="earn")
        p <- p + ggplot2::labs(x = "sensitivity MR",
                               title = "earnings simulation:\nMR-based sensitivity analysis")

    if (outcome=="depress" && sens.type=="GOR")
        p <- p + ggplot2::labs(x = "sensitivity GOR",
                               title = "depressive symptoms simulation:\nGOR-based sensitivity analysis")

    if (outcome=="depress" && sens.type=="SMDe")
        p <- p + ggplot2::labs(x = "sensitivity SMD",
                               title = "depressive symptoms simulation:\nSMDe-based sensitivity analysis")


    p
}


work.OR <- ..plot_for_fig8(outcome = "work", sens.type = "OR")

legnd <- cowplot::get_legend(work.OR)

work.OR <- work.OR + ggplot2::theme(legend.position = "none")

depress.GOR <- ..plot_for_fig8(outcome = "depress", sens.type = "GOR") +
    ggplot2::theme(legend.position = "none")

depress.SMDe <- ..plot_for_fig8(outcome = "depress", sens.type = "SMDe") +
    ggplot2::theme(legend.position = "none")

earn.MR <- ..plot_for_fig8(outcome = "earn", sens.type = "MR") +
    ggplot2::theme(legend.position = "none")


figs.dir <- here::here("research", "simulation", "outputs", "figures")

pdf(file.path(figs.dir, "fig8.pdf"), width = 11, height = 9)
cowplot::plot_grid(cowplot::plot_grid(work.OR,
                                      earn.MR,
                                      depress.GOR,
                                      depress.SMDe,
                                      ncol = 2),
                   legnd,
                   ncol = 2,
                   rel_widths = c(5,2))
dev.off()
