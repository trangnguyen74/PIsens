


..plot_sim_result <- function(sim,
                             # true,
                             sens.type,
                             ylimits,
                             ylimits.std,
                             include.Y1 = FALSE) {


    sim.errors  <- sim$errors
    sens.params <- sim$sens.params
    rm(sim)

    # true.se <- true$se
    # rm(true)


    if (sens.type %in% c("OR", "GOR", "MR")) {
        x.ticks <- unique(round(sens.params[,"sens.para.coord"]))
        x.tick.labels <- ifelse(x.ticks==0, as.character(1),
                                ifelse(x.ticks>0, as.character(x.ticks+1),
                                       paste0("1/", as.character(abs(x.ticks)+1))))

        x.label <- paste("sensitivity", sens.type)
    }

    if (sens.type %in% c("SMDe", "SMD")) {
        x.ticks <- unique(round(sens$sens.para.coord * 2)/2)
        x.tick.labels <- x.ticks

        x.label <- "sensitivity SMD"
    }

    color.values <- c("black", "blue", "dodgerblue")

    all.estimands <- c("cace", "nace", "tau.0c", "tau.0n", "tau.1c", "tau.1n")
    c.estimands   <- c("cace", "tau.0c", "tau.1c")
    effects       <- c("cace", "nace")




    ########################################
    # BIAS PLOT

    bias0 <- apply(sim.errors$bc0, c(1,2), function(z) mean(z, na.rm = TRUE))
    bias1 <- apply(sim.errors$bc1, c(1,2), function(z) mean(z, na.rm = TRUE))
    bias2 <- apply(sim.errors$bc2, c(1,2), function(z) mean(z, na.rm = TRUE))

    dat <- rbind(cbind(data.frame(cbind(sens.params, bias0)), estimate = "uncorrected"),
                 cbind(data.frame(cbind(sens.params, bias1)), estimate = "BC1"),
                 cbind(data.frame(cbind(sens.params, bias2)), estimate = "BC2"))
    dat$estimate <- factor(dat$estimate, levels = c("uncorrected", "BC1", "BC2"))

    dat <- tidyr::pivot_longer(dat,
                               cols = all_of(all.estimands),
                               names_to = "estimand",
                               values_to = "metric")
    dat$stratum <- ifelse(dat$estimand %in% c.estimands, "complier", "noncomplier")
    dat$stratum <- factor(dat$stratum, levels = c("noncomplier", "complier"))
    dat$estimand <- ifelse(dat$estimand %in% effects, "average effect",
                           ifelse(dat$estimand %in% c("tau.0c", "tau.0n"), "mean Y0", "mean Y1"))

    if (!include.Y1) dat <- dat[dat$estimand!="mean Y1",]




    # se <- data.frame(cbind(sens.params, true.se))
    # se <- tidyr::pivot_longer(se,
    #                           cols = all.estimands,
    #                           names_to = "estimand",
    #                           values_to = "metric")
    # se$stratum <- ifelse(se$estimand %in% c.estimands, "complier", "noncomplier")
    # se$stratum <- factor(se$stratum, levels = c("noncomplier", "complier"))
    # se$estimand <- ifelse(se$estimand %in% effects, "average effect",
    #                       ifelse(se$estimand %in% c("tau.0c", "tau.0n"), "mean Y0", "mean Y1"))
    #
    # if (!include.Y1) se <- se[se$estimand!="mean Y1",]
    #
    # se.minus <- se
    # se.minus$metric <- -se.minus$metric
    # se <- rbind(cbind(se,       pos = "above"),
    #             cbind(se.minus, pos = "below"))
    # rm(se.minus)




    p.bias <-
        ggplot2::ggplot(dat,
                        ggplot2::aes(x = sens.para.coord,
                                     y = metric,
                                     color = estimate)) +
        ggplot2::geom_hline(yintercept = 0, color = "lightgray") +
        ggplot2::geom_vline(xintercept = 0, color = "lightgray") +
        # ggplot2::geom_line(data = se, ggplot2::aes(group = pos), color = "gray70") +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = x.ticks, labels = x.tick.labels) +
        ggplot2::scale_color_manual(name = "", values = color.values) +
        ggplot2::labs(x = x.label, y = "", title = "Bias") +
        ggplot2::facet_grid(estimand ~ stratum) +
        ggplot2::theme_bw()

    if (!missing(ylimits))
        p.bias <- p.bias +
        ggplot2::coord_cartesian(ylim = ylimits)






    ########################################
    # STANDARDIZED BIAS PLOT

    se0 <- apply(sim.errors$bc0, c(1,2), function(z) sd(z, na.rm = TRUE))
    se1 <- apply(sim.errors$bc1, c(1,2), function(z) sd(z, na.rm = TRUE))
    se2 <- apply(sim.errors$bc2, c(1,2), function(z) sd(z, na.rm = TRUE))

    bias0 <- bias0 / se0
    bias1 <- bias1 / se0
    bias2 <- bias2 / se0



    dat <- rbind(cbind(data.frame(cbind(sens.params, bias0)), estimate = "uncorrected"),
                 cbind(data.frame(cbind(sens.params, bias1)), estimate = "BC1"),
                 cbind(data.frame(cbind(sens.params, bias2)), estimate = "BC2"))
    dat$estimate <- factor(dat$estimate, levels = c("uncorrected", "BC1", "BC2"))

    dat <- tidyr::pivot_longer(dat,
                               cols = all_of(all.estimands),
                               names_to = "estimand",
                               values_to = "metric")
    dat$stratum <- ifelse(dat$estimand %in% c.estimands, "complier", "noncomplier")
    dat$stratum <- factor(dat$stratum, levels = c("noncomplier", "complier"))
    dat$estimand <- ifelse(dat$estimand %in% effects, "average effect",
                           ifelse(dat$estimand %in% c("tau.0c", "tau.0n"), "mean Y0", "mean Y1"))

    if (!include.Y1) dat <- dat[dat$estimand!="mean Y1",]




    # se <- data.frame(cbind(sens.params, true.se / se0))
    # se <- tidyr::pivot_longer(se,
    #                           cols = all.estimands,
    #                           names_to = "estimand",
    #                           values_to = "metric")
    # se$stratum <- ifelse(se$estimand %in% c.estimands, "complier", "noncomplier")
    # se$stratum <- factor(se$stratum, levels = c("noncomplier", "complier"))
    # se$estimand <- ifelse(se$estimand %in% effects, "average effect",
    #                       ifelse(se$estimand %in% c("tau.0c", "tau.0n"), "mean Y0", "mean Y1"))
    #
    # if (!include.Y1) se <- se[se$estimand!="mean Y1",]
    #
    # se.minus <- se
    # se.minus$metric <- -se.minus$metric
    # se <- rbind(cbind(se,       pos = "above"),
    #             cbind(se.minus, pos = "below"))
    # rm(se.minus)


    p.std.bias <-
        ggplot2::ggplot(dat,
                        ggplot2::aes(x = sens.para.coord,
                                     y = metric,
                                     color = estimate)) +
        ggplot2::geom_hline(yintercept = 0, color = "lightgray") +
        ggplot2::geom_vline(xintercept = 0, color = "lightgray") +
        # ggplot2::geom_line(data = se, ggplot2::aes(group = pos), color = "gray70") +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = x.ticks, labels = x.tick.labels) +
        ggplot2::scale_color_manual(name = "", values = color.values) +
        ggplot2::labs(x = x.label, y = "", title = "Standardized bias") +
        ggplot2::facet_grid(estimand ~ stratum) +
        ggplot2::theme_bw()

    if (!missing(ylimits.std))
        p.std.bias <- p.std.bias  +
        ggplot2::coord_cartesian(ylim = ylimits.std)

    rm(bias0, bias1, bias2)


    # Standardized SE difference plot

    se.diff1 <- (se1 - se0) / se0
    se.diff2 <- (se2 - se0) / se0

    rm(se0, se1, se2)

    dat <- rbind(cbind(data.frame(cbind(sens.params, se.diff1)), estimate = "BC1"),
                 cbind(data.frame(cbind(sens.params, se.diff2)), estimate = "BC2"))

    dat$estimate <- factor(dat$estimate, levels = c("BC1", "BC2"))

    dat <- tidyr::pivot_longer(dat,
                               cols = all_of(all.estimands),
                               names_to = "estimand",
                               values_to = "metric")
    dat$stratum <- ifelse(dat$estimand %in% c.estimands,
                          "complier", "noncomplier")
    dat$stratum <- factor(dat$stratum, levels = c("noncomplier", "complier"))
    dat$estimand <- ifelse(dat$estimand %in% c("cace", "nace"),
                           "average effect",
                           ifelse(dat$estimand %in% c("tau.0c", "tau.0n"),
                                  "mean Y0", "mean Y1"))

    if (!include.Y1) dat <- dat[dat$estimand!="mean Y1",]


    p.std.SE.diff <-
        ggplot2::ggplot(dat,
                        ggplot2::aes(x = sens.para.coord,
                                     y = metric,
                                     color = estimate)) +
        ggplot2::geom_hline(yintercept = 0, color = "lightgray") +
        ggplot2::geom_vline(xintercept = 0, color = "lightgray") +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = x.ticks, labels = x.tick.labels) +
        ggplot2::scale_color_manual(name = "", values = color.values[2:3]) +
        ggplot2::labs(x = x.label, y = "", title = "Standardized SE change") +
        ggplot2::facet_grid(estimand ~ stratum) +
        ggplot2::theme_bw()

    if (!missing(ylimits.std))
        p.std.SE.diff <- p.std.SE.diff +
        ggplot2::coord_cartesian(ylim = ylimits.std)



    # Combine into figure
    p.legend <- cowplot::get_legend(
        p.bias +
            ggplot2::guides(color = ggplot2::guide_legend(nrow = 3)))

    cowplot::plot_grid(
        p.bias        + ggplot2::theme(legend.position = "none"),
        p.legend,
        p.std.bias    + ggplot2::theme(legend.position = "none"),
        p.std.SE.diff + ggplot2::theme(legend.position = "none"),
        nrow = 2,
        ncol = 2,
        rel_widths = 1,
        rel_heights = 1
    )
}








#######################################
########################################
########################################
########################################
########################################
########################################
# MAY DISCARD


process_sim <- function(sens.params,
                        error0, error1, error2) {

    sim.ok <- sum(!is.na(error0[1,1,]))


    bias0 <- apply(error0, c(1,2), function(z) mean(z, na.rm = TRUE))
    bias1 <- apply(error1, c(1,2), function(z) mean(z, na.rm = TRUE))
    bias2 <- apply(error2, c(1,2), function(z) mean(z, na.rm = TRUE))

    se0 <- apply(error0, c(1,2), function(z) sd(z, na.rm = TRUE))
    se1 <- apply(error1, c(1,2), function(z) sd(z, na.rm = TRUE))
    se2 <- apply(error2, c(1,2), function(z) sd(z, na.rm = TRUE))

    std.bias0 <- bias0 / se0
    std.bias1 <- bias1 / se0
    std.bias2 <- bias2 / se0

    std.se.diff1 <- (se1 - se0) / se0
    std.se.diff2 <- (se2 - se0) / se0

    bias0 <- cbind(sens.params, bias0)
    bias1 <- cbind(sens.params, bias1)
    bias2 <- cbind(sens.params, bias2)

    std.bias0 <- cbind(sens.params, std.bias0)
    std.bias1 <- cbind(sens.params, std.bias1)
    std.bias2 <- cbind(sens.params, std.bias2)

    se0 <- cbind(sens.params, se0)
    se1 <- cbind(sens.params, se1)
    se2 <- cbind(sens.params, se2)

    std.se.diff1 <- cbind(sens.params, std.se.diff1)
    std.se.diff2 <- cbind(sens.params, std.se.diff2)

    list(errors = list(bc0 = error0,
                       bc1 = error1,
                       bc2 = error2),
         bias = list(bc0 = bias0,
                     bc1 = bias1,
                     bc2 = bias2),
         se = list(bc0 = se0,
                   bc1 = se1,
                   bc2 = se2),
         std.bias = list(bc0 = std.bias0,
                         bc1 = std.bias1,
                         bc2 = std.bias2),
         std.se.diff = list(bc1 = std.se.diff1,
                            bc2 = std.se.diff2))
}

work.sim <- process_sim(sens.params = sens.params,
                        error0 = error0,
                        error1 = error1,
                        error2 = error2)


plot_simulation_metric <- function(standard, bc1, bc2,
                                   xlabel, ylabel, plot.title, ylimits) {

    dat <- rbind(cbind(data.frame(standard), estimate = "uncorrected"),
                 cbind(data.frame(bc1),      estimate = "BC1"),
                 cbind(data.frame(bc2),      estimate = "BC2"))

    dat$estimate <- factor(dat$estimate, levels = c("uncorrected", "BC1", "BC2"))

    color.values <- c("black", "blue", "dodgerblue")


    dat <- tidyr::pivot_longer(dat,
                               cols = c("cace", "nace", "tau.0c", "tau.0n"),
                               names_to = "estimand",
                               values_to = "metric")
    dat$stratum <- ifelse(dat$estimand %in% c("cace", "tau.0c"),
                          "complier", "noncomplier")
    dat$stratum <- factor(dat$stratum, levels = c("noncomplier", "complier"))
    dat$estimand <- ifelse(dat$estimand %in% c("cace", "nace"),
                           "average effect", "mean Y0")

    x.ticks <- unique(round(dat$sens.para.coord))
    x.tick.labels <- ifelse(x.ticks==0, as.character(1),
                            ifelse(x.ticks>0, as.character(x.ticks+1),
                                   paste0("1/", as.character(abs(x.ticks)+1))))

    p <- ggplot2::ggplot(dat,
                         ggplot2::aes(x = sens.para.coord,
                                      y = metric,
                                      color = estimate)) +
        ggplot2::geom_hline(yintercept = 0, color = "lightgray") +
        ggplot2::geom_vline(xintercept = 0, color = "lightgray") +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = x.ticks, labels = x.tick.labels) +
        ggplot2::scale_color_manual(name = "", values = color.values) +
        ggplot2::facet_grid(estimand ~ stratum) +
        ggplot2::theme_bw()

    if (!missing(xlabel))     p <- p + ggplot2::labs(x = xlabel)
    if (!missing(ylabel))     p <- p + ggplot2::labs(y = ylabel)
    if (!missing(plot.title)) p <- p + ggplot2::labs(title = plot.title)

    if (!missing(ylimits)) p <- p + ggplot2::coord_cartesian(ylim = ylimits)

    p
}



plot_simulation_std.se.increase <- function(bc1, bc2,
                                            xlabel, ylabel, plot.title, ylimits) {

    dat <- rbind(cbind(data.frame(bc1), estimate = "BC1"),
                 cbind(data.frame(bc2), estimate = "BC2"))

    color.values <- c("blue", "dodgerblue")

    dat <- dat[!names(dat) %in% c("tau.1c", "tau.1n")]


    dat <- tidyr::pivot_longer(dat,
                               cols = c("cace", "nace", "tau.0c", "tau.0n"),
                               names_to = "estimand",
                               values_to = "metric")
    dat$stratum <- ifelse(dat$estimand %in% c("cace", "tau.0c"),
                          "complier", "noncomplier")
    dat$stratum <- factor(dat$stratum, levels = c("noncomplier", "complier"))
    dat$estimand <- ifelse(dat$estimand %in% c("cace", "nace"),
                           "average effect", "mean Y0")

    x.ticks <- unique(round(dat$sens.para.coord))
    x.tick.labels <- ifelse(x.ticks==0, as.character(1),
                            ifelse(x.ticks>0, as.character(x.ticks+1),
                                   paste0("1/", as.character(abs(x.ticks)+1))))

    p <- ggplot2::ggplot(dat,
                         ggplot2::aes(x = sens.para.coord,
                                      y = metric,
                                      color = estimate)) +
        ggplot2::geom_hline(yintercept = 0, color = "lightgray") +
        ggplot2::geom_vline(xintercept = 0, color = "lightgray") +
        ggplot2::geom_line() +
        ggplot2::labs(ylab = "standardized SE increase",
                      title = "Standardized SE increase") +
        ggplot2::scale_x_continuous(breaks = x.ticks, labels = x.tick.labels) +
        ggplot2::scale_color_manual(name = "", values = color.values) +
        ggplot2::facet_grid(estimand ~ stratum) +
        ggplot2::theme_bw()


    if (!missing(xlabel))     p <- p + ggplot2::labs(x = xlabel)
    if (!missing(ylabel))     p <- p + ggplot2::labs(y = ylabel)
    if (!missing(plot.title)) p <- p + ggplot2::labs(title = plot.title)

    if (!missing(ylimits)) p <- p + ggplot2::coord_cartesian(ylim = ylimits)

    p
}




p.bias <- plot_simulation_metric(standard   = work.sim$bias$bc0,
                                 bc1        = work.sim$bias$bc1,
                                 bc2        = work.sim$bias$bc2,
                                 xlabel     = "sensitivity OR",
                                 ylabel     = "",
                                 plot.title = "Bias")

p.std.bias <- plot_simulation_metric(standard   = work.sim$std.bias$bc0,
                                     bc1        = work.sim$std.bias$bc1,
                                     bc2        = work.sim$std.bias$bc2,
                                     xlabel     = "sensitivity OR",
                                     ylabel     = "",
                                     plot.title = "Standardized bias",
                                     ylimits    = c(-.6,.6))


p.std.se.increase <-
    plot_simulation_std.se.increase(bc1        = work.sim$std.se.diff$bc1,
                                    bc2        = work.sim$std.se.diff$bc2,
                                    xlabel     = "sensitivity OR",
                                    ylabel     = "",
                                    plot.title = "Standardized SE increase",
                                    ylimits    = c(-.6,.6))

p.legend <- cowplot::get_legend(p.bias + ggplot2::guides(color = ggplot2::guide_legend(nrow = 3)))

p.work.sim <- cowplot::plot_grid(
    p.bias + ggplot2::theme(legend.position = "none"),
    p.legend,
    p.std.bias + ggplot2::theme(legend.position = "none"),
    p.std.se.increase + ggplot2::theme(legend.position = "none"),
    nrow = 2,
    ncol = 2,
    rel_widths = 1,
    rel_heights = 1
)































