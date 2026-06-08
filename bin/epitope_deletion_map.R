#!/usr/bin/env Rscript

source(file.path(Sys.getenv('EPISWEEP'), 'bin', 'episweep_plot.R'))

main <- function()
{
    library(getopt)
    spec <- matrix(c(
            'wild_type', 'w', 1, "character", "Wild type epitope prediction [CSV]",
            'variant', 'v', 1, "character", "Variant epitope prediction [CSV]",
            'output' , 'o', 1, "character", "Output figure file [PDF]",
            'help', 'h', 0, "logical",   "This message"
    ),ncol=5,byrow=T)
    
    opt = getopt(spec)
    
    if (length(opt) < 2) {
            cat(paste(getopt(spec, usage=T),"\n"));
        quit("no");
    }
	
    wt_epitope <- read.csv(opt$wild_type,row.names=1, header=T)
	var_epitope <- read.csv(opt$variant,row.names=1, header=T)

    pdf(opt$output,width=15,height=10)
    par(mar=c(3,5,4,0.5))
	plot.delta.alleles(wt_epitope, var_epitope)
}

if (!interactive()) main()
