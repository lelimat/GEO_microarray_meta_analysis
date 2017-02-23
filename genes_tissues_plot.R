# genes_tissues_plot.R
# C: Jan 30, 3017
# M: Feb 23, 3017
# A: Leandro Lima <leandro.lima@gladstone.ucsf.edu>

library('Biobase')
library('ggplot2')

setwd('~/Dropbox (Gladstone)/meta_analysis_PD/data/')

tissues <- c("striatum", "cortex", "medulla", "frontal_cortex", "locus_coeruleus",
             "globus_pallidus_interna", "putamen", "cerebellum", 
             "dorsal_nucleus_of_vagus_nerve", "substantia_nigra", "inferior_olivary_nucleus",
             "whole_blood", "peripheral_blood", "blood")

alltdir <- 'tissue_dirs'
all_genes <- NULL
study_tissues <- NULL

for (tissue in tissues) {
    tdir <- paste(alltdir, '/', tissue, sep='')
    for (csv_file in dir(tdir, pattern='*_all_genes.csv')) {
    #for (csv_file in csv_files) {
        #csv_file <- csv_files[2]
        study_tissues <- c(study_tissues, tissue)
        csv <- read.csv(paste(tdir, csv_file, sep='/'))
        cat('Opening', csv_file, '-', tissue, '\n')
        study_info <- gsub('_', '.', strsplit(csv_file, '_all')[[1]][1])
        if (is.null(all_genes)) {
            all_genes <- data.frame(csv[, c('SYMBOL', 'logFC')], rowMedians(as.matrix(csv[, c('adjP.LIMMA', 'adjP.RP', 'adjP.SAMR')])))
            all_genes <- all_genes[which(!is.na(all_genes$SYMBOL)), ]
            colnames(all_genes)[3] <- 'adjP'
            symbol_logFC <- aggregate(logFC~SYMBOL, data=all_genes, FUN=function(x) {x[which.max(abs(x))]})
            all_genes <- merge(symbol_logFC, all_genes, by.x = c(1,2), by.y = c(1,2)) #, all.x=TRUE)
            
        } else {
            study_genes <- data.frame(csv[, c('SYMBOL', 'logFC')], rowMedians(as.matrix(csv[, c('adjP.LIMMA', 'adjP.RP', 'adjP.SAMR')])))
            colnames(study_genes)[3] <- 'adjP'
            symbol_logFC <- aggregate(logFC~SYMBOL, data=study_genes, FUN=function(x) {x[which.max(abs(x))]})
            study_genes <- merge(symbol_logFC, study_genes, by.x = c(1,2), by.y = c(1,2)) #, all.x=TRUE)
            study_genes <- aggregate(study_genes$logFC, list('SYMBOL'=study_genes$SYMBOL, 'adjP'=study_genes$adjP), FUN=function(x) {x[which.max(abs(x))]})
            colnames(study_genes)[3] <- 'logFC'
            study_genes <- study_genes[, c('SYMBOL', 'logFC', 'adjP')]
            all_genes <- merge(all_genes, study_genes, by=c(1,1)) # MERGE
        }
        
        colnames(all_genes)[ncol(all_genes)-1] <- paste('logFC', study_info, tissue, sep='.')
        colnames(all_genes)[ncol(all_genes)]   <- paste('medianP', study_info, tissue, sep='.')

    }
}

# colnames(all_genes) <- gsub('_GPL', '.GPL', colnames(all_genes))
logFC_table <- all_genes[, grep('logFC', colnames(all_genes))]
Pval_table <- all_genes[, grep('medianP', colnames(all_genes))]
studies <- as.data.frame(strsplit(colnames(all_genes)[grep('medianP', colnames(all_genes))], '\\.'), stringsAsFactors=FALSE)
gse <- studies[2, ]
gpl <- studies[3, ]
studies <- paste(gse, gpl, sep = '.')

n_blood <- length(grep('blood', study_tissues))
n_brain <- length(study_tissues) - n_blood
write.csv(all_genes, 'all_genes_for_barplot.csv', quote=FALSE, row.names = FALSE)

# i = which(all_genes$SYMBOL=='DNAJB1')
plot_all = TRUE

if (plot_all) {
    dest_dir = 'all/'
} else {
    dest_dir = ''
}

for (i in 1:nrow(all_genes)) {
    logFC <- logFC_table[i, ]
    Pval  <- Pval_table[i, ]
    if (plot_all || (max(abs(logFC))>=2 && min(Pval, na.rm=TRUE)<=0.05)) {
        symbol <- all_genes$SYMBOL[i]
        logFC_Pval <- data.frame(studies, t(logFC), t(Pval))
        colnames(logFC_Pval)[2:3] <- c('logFC', 'Pval')
        logFC_Pval$study_tissues <- gsub('_', ' ', study_tissues)
        
        png(paste('barplots/', dest_dir, 'barplot_', symbol, '.png', sep=''), width = 600, height = 600)
        par(las=2, mar=c(18,4.1,4.1,2.1))
        logFC_Pval$colors <- ifelse(logFC_Pval$Pval <= 0.05, ifelse(logFC_Pval$logFC > 0, 'red', 'green'), 'grey')
        end_point = 0.5 + nrow(logFC_Pval) + nrow(logFC_Pval)-1
        barplot(logFC_Pval$logFC, col=logFC_Pval$colors, border=logFC_Pval$colors, ylim=c(-3, 3),
                main=symbol, xlab = "", space=1)
        abline(v=2*n_brain+0.5, col='grey', lty=3)
        text(seq(1.5, end_point, by=2), par("usr")[3]-0.25, 
             srt = 60, adj= 1, xpd = TRUE, col=c(rep('darkblue', n_brain), rep('darkred', n_blood)),
             labels = paste(gse, logFC_Pval$study_tissues, sep=' | '), cex=1)
        legend('topright', legend=c('Up (AdjP<=0.05)', 'Not significant', 'Down (AdjP<=0.05)'), col=c('red', 'grey', 'green'), pch=15)
        par(mar=c(5.1,4.1,4.1,2.1))
        dev.off()
    }
}


panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = (cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * abs(r))
}


#jpeg('~/Dropbox/doutorado/tdah/paper_2014/figures/SNVs_CNVs_correls.jpg', width=5, height=5, units="in", res=500)
jpeg('studies_pairwise_correlation_tissue.jpg', width=50, height=50, units="in", res=500)
pairs(logFC_table, lower.panel=panel.smooth, upper.panel=panel.cor,
      main="pair-wise correlations between studies", labels=gsub('_', '\n', study_tissues))
dev.off()

jpeg('studies_pairwise_correlation_gpl.jpg', width=50, height=50, units="in", res=500)
pairs(logFC_table, lower.panel=panel.smooth, upper.panel=panel.cor,
      main="pair-wise correlations between studies", labels=gpl)
dev.off()

jpeg('studies_pairwise_correlation_gse.jpg', width=50, height=50, units="in", res=500)
pairs(logFC_table, lower.panel=panel.smooth, upper.panel=panel.cor,
      main="pair-wise correlations between studies", labels=gse)
dev.off()
