

library(PIsens)

jobs.dir <- here::here("research", "jobsII-analysis")

IFH1.dir    <- file.path(jobs.dir, "outputs", "results", "IFH1-only", "sens")
figs.dir <- file.path(jobs.dir, "outputs", "figures")


earn.plot <-
    plot_sens(sens          = readRDS(file.path(IFH1.dir, "earn_MR.rds")),
              CI            = "BCa",
              point         = "corrected1",
              effect.ylims  = c(-1500,1500),
              outcome.ylims = c(0,3000)) +
    ggplot2::labs(title = "earnings:\nMR-based sensitivity analysis")

work.plot <-
    plot_sens(sens          = readRDS(file.path(IFH1.dir, "work_OR.rds")),
              CI            = "BCa",
              point         = "corrected1",
              outcome.ylims = c(0,1),
              effect.ylims  = c(-.5,.5)) +
    ggplot2::labs(title = "work for pay:\nOR-based sensitivity analysis")

depress.GOR.plot <-
    plot_sens(sens           = readRDS(file.path(IFH1.dir, "depress_GOR.rds")),
              CI             = "BCa",
              point          = "corrected1",
              effect.ylims = c(-2,2),
              outcome.ylims = c(1,5)) +
    ggplot2::labs(title = "depressive symptoms:\nGOR-based sensitivity analysis")

depress.SMDe.plot <-
    plot_sens(sens           = readRDS(file.path(IFH1.dir, "depress_SMDe.rds")),
              CI             = "BCa",
              point          = "corrected1",
              effect.ylims = c(-2,2),
              outcome.ylims = c(1,5)) +
    ggplot2::labs(title = "depressive symptoms:\nSMDe-based sensitivity analysis")

lgnd <- cowplot::get_legend(
    work.plot +
        guides(color = guide_legend(nrow = 3)) +
        theme(legend.position = "right")
)

all.sens.plot <- cowplot::plot_grid(
    work.plot         + theme(legend.position = "none"),
    earn.plot         + theme(legend.position = "none"),
    depress.GOR.plot  + theme(legend.position = "none"),
    depress.SMDe.plot + theme(legend.position = "none"),
    align = 'vh',
    nrow = 2,
    ncol = 2,
    rel_widths = 1,
    rel_heights = 1
)

pdf(file.path(figs.dir, "fig2.pdf"), width = 9, height = 9)
cowplot::plot_grid(all.sens.plot, lgnd, ncol = 2, rel_widths = c(5, 1))
dev.off()
