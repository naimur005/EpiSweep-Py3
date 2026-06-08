#!/usr/bin/env Rscript

source(file.path(Sys.getenv('EPISWEEP'), 'bin', 'episweep_plot.R'))

main <- function()
{
    library(getopt)
    spec <- matrix(c(
            'library_plan', 'i', 1, "character", "Library plan file [CSV]",
            'output' , 'o', 1, "character", "Output figure file [PDF]",
            'help', 'h', 0, "logical",   "This message"
    ),ncol=5,byrow=T)
    
    opt = getopt(spec)
    
    if (length(opt) < 2) {
            cat(paste(getopt(spec, usage=T),"\n"));
        quit("no");
    }

    csvfile <- read.csv(opt$library_plan, header=T)
    output <- opt$output

    epimax <- max(csvfile$episcore)
    epimin <- min(csvfile$episcore)
    pot_max <- max(csvfile$potential)
    pot_min <- min(csvfile$potential)
    num_curve <- max(csvfile$curve)

    original <- subset(csvfile, csvfile$index==0)

    pdf(file=output, width=8, height=6)
    par(mar=c(5,5,0.5,0.5))

    lib_plan <- subset(csvfile, csvfile$curve!=0)

    plot(lib_plan$episcore, lib_plan$potential,
         col="black", pch=1,
         xlab="Epitope score", ylab="Potential score",
         cex.axis=1.5, cex.lab=2, cex=1,
         xlim=c(epimin, epimax), ylim=c(pot_min, pot_max), xaxt="n")
    
    #points(original$episcore, original$potential, pch=1, col="blue", cex=5, lwd=2)

    axis(1,at=csvfile$episcore, labels=csvfile$episcore, cex.axis=1.5)
    
    garbage <- dev.off()
}
if (!interactive()) main()
