

jobs.dir    <- here::here("research", "jobsII-analysis")
sens.dir    <- file.path(jobs.dir, "outputs", "results", "IFH1-only", "sens")
figs.dir <- file.path(jobs.dir, "outputs", "figures")


########################################
#### figure 3

work.OR <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "work_OR.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            sens.type = "OR",
                            include.bias.correction = TRUE,
                            include.CI = TRUE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL)

pdf(file.path(figs.dir, "fig3.pdf"), width = 7, height = 5)
work.OR
dev.off()




########################################
#### figure 6

pattern.legend <- cowplot::get_legend(work.OR)

work.OR <- work.OR +
    ggplot2::labs(x = "sensitivity OR",
                  title = "work for pay:\nOR-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")

depress.GOR <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "depress_GOR_symmetric.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            sens.type = "GOR",
                            include.bias.correction = TRUE,
                            include.CI = TRUE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL) +
    ggplot2::labs(x = "sensitivity GOR",
                  title = "depressive symptoms:\nGOR-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")


depress.SMDe <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "depress_SMDe_symmetric.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            sens.type = "SMDe",
                            include.bias.correction = TRUE,
                            include.CI = TRUE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL) +
    ggplot2::labs(x = "sensitivity SMD",
                  title = "depressive symptoms:\nSMDe-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")

earn.MR <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "earn_MR.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            sens.type = "MR",
                            include.bias.correction = TRUE,
                            include.CI = TRUE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL) +
    ggplot2::labs(x = "sensitivity MR",
                  title = "earnings:\nMR-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")


pdf(file.path(figs.dir, "fig6.pdf"), width = 11, height = 9)
cowplot::plot_grid(cowplot::plot_grid(work.OR,
                                      earn.MR,
                                      depress.GOR,
                                      depress.SMDe,
                                      ncol = 2),
                   pattern.legend,
                   ncol = 2,
                   rel_widths = c(5,2))
dev.off()




########################################
#### figure 7

work.OR.2 <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "work_OR.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            sens.type = "OR",
                            include.bias.correction = TRUE,
                            include.CI = FALSE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL)

pattern.legend.2 <- cowplot::get_legend(work.OR.2)


work.OR.2 <- work.OR.2 +
    ggplot2::labs(x = "sensitivity OR",
                  title = "work for pay:\nOR-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")

depress.GOR.2 <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "depress_GOR_symmetric.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            sens.type = "GOR",
                            include.bias.correction = TRUE,
                            include.CI = FALSE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL) +
    ggplot2::labs(x = "sensitivity GOR",
                  title = "depressive symptoms:\nGOR-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")

depress.SMDe.2 <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "depress_SMDe_symmetric.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            sens.type = "SMDe",
                            include.bias.correction = TRUE,
                            include.CI = FALSE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL) +
    ggplot2::labs(x = "sensitivity SMD",
                  title = "depressive symptoms:\nSMDe-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")

earn.MR.2 <-
    plot_finite.sample.bias(sens = readRDS(file.path(sens.dir, "earn_MR.rds"))
                            $sens[c("cace", "nace", "tau.0c", "tau.0n")],
                            include.bias.correction = TRUE,
                            include.CI = FALSE,
                            effect.ylims = NULL,
                            outcome.ylims = NULL,
                            sens.para.lims = NULL) +
    ggplot2::labs(x = "sensitivity MR",
                  title = "earnings:\nMR-based sensitivity analysis") +
    ggplot2::theme(legend.position = "none")








pdf(file.path(figs.dir, "fig7.pdf"), width = 11, height = 9)
cowplot::plot_grid(cowplot::plot_grid(work.OR.2,
                                      earn.MR.2,
                                      depress.GOR.2,
                                      depress.SMDe.2,
                                      ncol = 2),
                   pattern.legend.2,
                   ncol = 2,
                   rel_widths = c(5,2))
dev.off()













