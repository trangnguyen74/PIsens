
# These two functions come from Brad on this website
# https://stackoverflow.com/questions/23092580/r-source-function-by-name-import-subset-of-functions

loadFunction <- function(file,function.name) {

    eval(parse(text=paste(function.name," <- function(x) {0}",sep="")),envir = .GlobalEnv)
    suppressMessages(insertSource(file, functions=function.name))
    eval(parse(text=paste(function.name," <- ",function.name,"@.Data",sep="")),envir = .GlobalEnv)

}

unloadFunction <- function(function.name) {

    eval(parse(text=paste("rm(",function.name,",envir = .GlobalEnv)",sep="")))

}





base.dir <- here::here("research", "simulation")
sim.dir  <- file.path(base.dir, "outputs", "simresults")
figs.dir <- file.path(base.dir, "outputs", "figures")

loadFunction(file = file.path(base.dir, "script", "0b_sim_functions.R"),
             function.name = "..plot_sim_bias")



# figure 9: same as sim_bias_work_OR.pdf
pdf(file.path(figs.dir, "fig9.pdf"), width = 9, height = 11.5)
..plot_sim_bias(sim = readRDS(file.path(sim.dir, "work", "OR", "sim.rds")),
                sens.type = "OR",
                ylimits = c(-.05,.05),
                ylimits.std = c(-1,1))
dev.off()


# figure 10: same as sim_bias_depress_GOR.pdf
pdf(file.path(figs.dir, "fig10.pdf"), width = 9, height = 11.5)
..plot_sim_bias(sim = readRDS(file.path(sim.dir, "depress", "GOR", "sim.rds")),
                sens.type = "GOR",
                ylimits = c(-.2,.2),
                ylimits.std = c(-1.6,1.6))
dev.off()



# figure 11: same as sim_bias_depress_SMDe.pdf
pdf(file.path(figs.dir, "fig11.pdf"), width = 9, height = 11.5)
..plot_sim_bias(sim = readRDS(file.path(sim.dir, "depress", "SMDe", "sim.rds")),
                sens.type = "SMD",
                ylimits = c(-.2,.2),
                ylimits.std = c(-1,1))
dev.off()




# figure 12: same as sim_bias_earn_MR.pdf
pdf(file.path(figs.dir, "fig12.pdf"), width = 9, height = 11.5)
..plot_sim_bias(sim = readRDS(file.path(sim.dir, "earn", "MR", "sim.rds")),
                sens.type = "MR",
                ylimits = c(-150,150),
                ylimits.std = c(-1,1))
dev.off()
