#!/usr/bin/env Rscript

main <- function()
{
    args <- commandArgs(TRUE)

    csvfile <- read.csv("lst_CWBD_lib.csv", header=T)
    plan <- read.csv("lst_CWBD_lib_plan10.csv", header=T)


    output <- "lst_CWBD_lib_pareto.pdf"
    epimax <- max(plan$episcore)
    epimin <- min(plan$episcore)
    emax <- max(csvfile$energy)
    emin <- min(csvfile$energy)
    pdf(file=output, width=8, height=6)
    par(mar=c(5,5,0.5,0.5))
    wt <- subset(csvfile, csvfile$curve==0)
    pareto <- subset(csvfile, csvfile$curve==1)
    plot(31, 0, col="black", pch=19, xlab="Epitope score", ylab="Potential score", cex.axis=1.5, cex.lab=2, cex=2, xlim=c(epimin, epimax), ylim=c(emin, emax), xaxt="n")
    points(plan$episcore, plan$energy, pch=1, col="black", cex=0.2)
    color = rainbow(9)
    color = c("gray90",
              "gray80",
              "gray70",
              "gray60",
              "gray50",
              "gray40",
              "gray30",
              "gray20",
              "gray10")
    for (m in 1:9)
    {
        new_set = subset(plan, plan$nmut == m)
        points(new_set$episcore, new_set$energy, pch=1, col=color[m], cex=0.2)
    }
    #points(plan$episcore, plan$energy, pch=1, col="black", cex=0.2)
    for (m in 2:20)
    {
        this_set <- subset(csvfile, csvfile$curve==m)
        epi <- this_set$episcore
        e <- this_set$energy
        points(epi, e, pch=1, col="gray70", cex=1.5)
        lines(epi, e, pch=1, col="gray70", cex=2, lty=3)
    }
    points(pareto$episcore, pareto$energy, pch=19, col="red", cex=1)
    lines(pareto$episcore, pareto$energy, pch=19, col="red", cex=1)

    particular <- subset(pareto, pareto$episcore == 19.0)
    points(19.5, -0.75813625, pch=1, col="blue", cex=3, lwd=2)
    axis(1,at=plan$episcore, labels=plan$episcore, cex.axis=0.5)
    #axis.break(1,47.5,style="slash") 

    
    garbage <- dev.off()
}
if (!interactive()) main()

