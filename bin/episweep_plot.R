load_spec <- function(spec) {
    spec_file <- read.csv(file=spec, header=F)
    keys <- spec_file[,1]
    vals <- spec_file[,2]
    for (l in 1:length(keys)) {
        if (keys[l] == "design_output") {
            design_file = vals[l]
            return (as.character(design_file))
        }
    }
}

seq_parse <- function(prediction) {
    seq = ""
    for (l in 1:nrow(prediction)) {
        if (l == nrow(prediction)) {
            seq <- paste(seq, prediction$peptide[l], sep="")
        }
        else {
            seq <- paste(seq, substr(prediction$peptide[l], 1, 1), sep="")
        }
    }
    return (unlist(strsplit(seq, split="")))
}

make_title <- function(seq1, seq2, mutated_pos) {
    title <- c()
    for (m in mutated_pos)
    {
        title <- c(title, paste(seq1[m], m, seq2[m], sep=""))
    }
    return(paste(title, collapse=", "))
}


plot.epitope_hits <- function(epitope)
{
    xmax <- max(epitope$start)
    na <- length(epitope) - 4 # excluded start, peptide, end and total
    #na <- 26 #length(epitope) - 4 # excluded start, peptide, end and total
    plot(c(1, xmax), c(0,na), type='n', xlab='', ylab='Epitope Score', xaxt='n', yaxt='n', frame=FALSE, cex.lab=2)
    axis(1, epitope$start, pos=0, labels=FALSE)
    axis(2, 0:na, cex.axis=2)
    lines(epitope$start, epitope$total, col='black', type='h', lwd=7, lend='butt', lty=1)
    for (i in 1:xmax)
    {
        if (i%%5 == 0) {
            text(x=i, y=-0.5, i, cex=2, xpd=TRUE)
        }
    }
    mtext("Residue Number", side=1, line=3, cex=2)
}

plot.suboptimal <- function(designs, ncurve) {
    for (c in 2:ncurve) {
        sub <- subset(designs, designs$curve == c)
        points(sub$episcore, sub$potential, pch=1, col="gray70", cex=1.5)
        lines(sub$episcore, sub$potential, pch=1, col="gray70", lty=3)
    }
}

plot.pareto_optimal <- function(pareto) {
    points(pareto$episcore, pareto$potential, pch=19, col="red", cex=2)
    lines(pareto$episcore, pareto$potential, pch=19, col="red", cex=2)
}

plot.delta.alleles <- function(wt_epitope, var_epitope)
{
    wt_seq <- seq_parse(wt_epitope)
    var_seq <- seq_parse(var_epitope)

    mutated_pos <- which(wt_seq != var_seq)

    title <- make_title(wt_seq, var_seq, mutated_pos)

	na <- ncol(wt_epitope) # number of alleles = number of columns - 1 (the total)
	
    wt_epitope <- wt_epitope[,3:(na-1)]
	var_epitope <- var_epitope[,3:(na-1)]

	plot(c(1,length(wt_seq)), c(1,ncol(wt_epitope)),
         main=title, 
         cex.main=3, cex.lab=2, cex.axis=2,
         type='n', xlab='', ylab='', 
         xaxt='n', yaxt='n', 
         ylim=c(0.8, ncol(wt_epitope)+0.2))

    for (m in mutated_pos) {
		points(m, 0.1, pch=2, col='black', xpd=TRUE, cex=3)
    }
    axis(1, 1:length(wt_seq), rep(c(" "), length(wt_seq)), cex.axis=2.2)
	axis(2, 1:ncol(wt_epitope), colnames(wt_epitope), las=1, cex.axis=1.7)

    for (h_line in 1:(ncol(wt_epitope)-1)) {
        abline(h=h_line+0.5, lty=3, lwd=2, col="black")
    }

	for (a in 1:ncol(wt_epitope)) {
		offset <- 0; till <- 0
		for (p in 1:nrow(wt_epitope)) {
			if (wt_epitope[p,a] || var_epitope[p,a]) {
				if (p < till) offset <- offset + 0.2
				else offset <- -0.2
				till <- p+10
				if (wt_epitope[p,a] && var_epitope[p,a]) { lines(c(p,p+8), c(a+offset, a+offset), col="red", lwd=3) }
				if (wt_epitope[p,a] && !var_epitope[p,a]) { lines(c(p,p+8), c(a+offset, a+offset), col='blue', lwd=5, lend='square')  }
				if (!wt_epitope[p,a] && var_epitope[p,a]) { lines(c(p,p+8), c(a+offset, a+offset), col='red', lwd=5, lend='square', lty="42")  }
			}
		}
	}
}