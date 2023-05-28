

base.dir    <- here::here("research", "other-stuff")
script.dir  <- file.path(base.dir, "script")
figs.dir <- file.path(base.dir, "outputs", "figures")


pdf(file = file.path(figs.dir, "fig1.pdf"), width = 10, height = 4)
visualize_GOR(point.mu.0 = .5, point.c.wt = .25)
dev.off()
