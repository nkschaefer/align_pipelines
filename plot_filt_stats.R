#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    cat("Please provide stats file (from vcf_depth_filter), and optionally, the \
pops file given to vcf_depth_filter.\n")
    q()
}

library(ggplot2)
library(ggthemes)

statsf <- args[1]
statsfbase <- strsplit(statsf, '.', fixed=T)[[1]][1]

stats <- read.table(statsf)
colnames(stats) <- c("type", "id", "num1", "num2")

pop <- FALSE
if (length(args) > 1){
    popsf <- read.table(args[2])
    colnames(popsf) <- c("id", "pop")
    stats <- merge(stats, popsf, by="id")
    pop <- TRUE
}

maxes <- aggregate(stats[which(stats$type=="DP"),]$num2, by=list(id=stats[which(stats$type=="DP"),]$id), FUN=max)
colnames(maxes)[2] <- "maxy"
maxes <- merge(maxes, stats[which(stats$type=="MED_DP_GQ"),], by="id")

width <- 7
height <- 7
plt <- ""
if (pop){
    width <- 9
    plt <- ggplot(stats[which(stats$type=="DP"),]) + 
        geom_line(aes(x=num1, y=num2, colour=factor(pop), group=factor(id))) + 
        theme_bw() + 
        scale_x_log10("Coverage depth") + 
        scale_y_continuous("Frequency") + 
        annotation_logticks(side="b") + 
        geom_vline(data=stats[which(stats$type=="DP_CUTOFF"),], aes(xintercept=num1), lty="dotted") + 
        geom_vline(data=stats[which(stats$type=="DP_CUTOFF"),], aes(xintercept=num2), lty="dotted") + 
        scale_colour_gdocs("Population") + 
        geom_text(data=maxes, aes(label=id, x=num1, y=maxy, colour=pop))
} else{
    plt <- ggplot(stats[which(stats$type=="DP"),]) + 
        geom_line(aes(x=num1, y=num2, colour=id, group=id), show.legend=FALSE) +
        theme_bw() + 
        scale_colour_gdocs("ID") + 
        scale_x_log10("Coverage depth") + 
        scale_y_continuous("Frequency") + 
        annotation_logticks(side="b") + 
        geom_vline(data=stats[which(stats$type=="DP_CUTOFF"),], aes(xintercept=num1, colour=factor(id)), 
            lty="dotted", show.legend=FALSE) + 
        geom_vline(data=stats[which(stats$type=="DP_CUTOFF"),], aes(xintercept=num2, colour=factor(id)), 
            lty="dotted", show.legend=FALSE) + 
        geom_text(data=maxes, aes(label=id, x=num1, y=maxy, colour=factor(id)), show.legend=FALSE)
}

plt <- plt + theme(panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.text.x=element_text(size=12),
                   axis.text.y=element_text(size=12),
                   axis.title.x=element_text(size=14),
                   axis.title.y=element_text(size=14))

outfname <- paste(statsfbase, '.cov.pdf', sep='')
ggsave(outfname, plt, width=width, height=height)

# Make second plot for heterozygosity & missing data
if (pop){
    # Get stuff sorted by population to make bar plot look better
    cols <- unique(stats[,which(colnames(stats) %in% c("id", "pop"))])
    cols <- cols[order(cols$id),]
    cols <- cols[order(cols$pop),]
    ord <- cols$id
    stats$id <- factor(stats$id, labels=ord, levels=ord)
}

plt <- ggplot(stats[which(stats$type=="HET" | stats$type=="MISS"),])
if (pop){
    plt <- plt + geom_bar(aes(x=id, y=num1/num2, fill=factor(pop)), stat='identity') + 
    scale_fill_gdocs("Population")
} else{
    plt <- plt + geom_bar(aes(x=id, y=num1/num2, fill=factor(id)), stat='identity', show.legend=FALSE) + 
    scale_fill_gdocs("ID")
}

plt <- plt + theme_bw() + 
    facet_grid(type~., scales="free_y") + 
    scale_x_discrete("Sample") + 
    scale_y_continuous("Rate") + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          axis.text.y=element_text(size=12),
          axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1))

outfname <- paste(statsfbase, '.hetmiss.pdf', sep='')
ggsave(outfname, plt, width=width, height=height)

