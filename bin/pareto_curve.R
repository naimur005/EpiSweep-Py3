#!/usr/bin/env Rscript

source(file.path(Sys.getenv('EPISWEEP'), 'bin', 'episweep_plot.R'))

main <- function()
{
    library(getopt)
    spec <- matrix(c(
            'design_spec', 'i', 1, "character", "Design spec file [CSV]",
            'output' , 'o', 1, "character", "Output figure file [PDF]",
            'help', 'h', 0, "logical",   "This message"
    ),ncol=5,byrow=T)
    
    opt = getopt(spec)
    
    if (length(opt) < 2) {
            cat(paste(getopt(spec, usage=T),"\n"));
        quit("no");
    }
    design <- load_spec(opt$design_spec)

    csvfile <- read.csv(design, header=T)
    output <- opt$output

    epimax <- max(csvfile$episcore)
    epimin <- min(csvfile$episcore)
    pot_max <- max(csvfile$potential)
    pot_min <- min(csvfile$potential)
    num_curve <- max(csvfile$curve)

    wt <- subset(csvfile, csvfile$index==0)

    pdf(file=output, width=8, height=6)
    par(mar=c(5,5,0.5,0.5))

    pareto <- subset(csvfile, csvfile$curve==1)

    plot(wt$episcore, wt$potential,
         col="black", pch=19,
         xlab="Epitope score", ylab="Potential score",
         cex.axis=1.5, cex.lab=2, cex=2,
         xlim=c(epimin, epimax), ylim=c(pot_min, pot_max), xaxt="n")

    plot.suboptimal(csvfile, num_curve)
    plot.pareto_optimal(pareto)
    
    points(wt$episcore, wt$potential, pch=19, col="black", cex=2)

    axis(1,at=csvfile$episcore, labels=csvfile$episcore, cex.axis=1.5)
    
    garbage <- dev.off()
}
if (!interactive()) main()
