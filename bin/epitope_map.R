#!/usr/bin/env Rscript

source(file.path(Sys.getenv('EPISWEEP'), 'bin', 'episweep_plot.R'))

main <- function()
{
    library(getopt)
    spec <- matrix(c(
            'epitope', 'i', 1, "character", "Epitope prediction file [CSV]",
            'output' , 'o', 1, "character", "Output figure file [PDF]",
            'help', 'h', 0, "logical",   "This message"
    ),ncol=5,byrow=T)
    
    opt = getopt(spec)
    
    if (length(opt) < 2) {
            cat(paste(getopt(spec, usage=T),"\n"));
        quit("no");
    }

    epitope_prediction <- read.csv(opt$epitope, header=T)

    ###### Change if you want to change size of figure and margin etc
    pdf(opt$output,width=15,height=7.5)
    par(mar=c(5,5,1,0))
    
    plot.epitope_hits(epitope_prediction)
}

if (!interactive()) main()
