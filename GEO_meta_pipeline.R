# GEO_meta_pipeline.R
# C: Sep 15, 2016
# M: Feb 16, 2017
# A: Leandro Lima <leandro.lima@gladstone.ucsf.edu>

# Usage: R --slave --file=GEO_meta_pipeline.R --args [GSE] [GPL] [optional: samples_info.csv] [results_dir] [tissue]

# Run in parallel
# rm out*txt
# tail -n +2 Hs_samples_final.csv | cut -d, -f2,3 | sort | uniq | awk -F"," '{gsub(/"/, "", $0); print "R --slave --file=GEO_meta_pipeline.R --args  "$1" "$2" > out1."$1"."$2".txt 2> out2."$1"."$2".txt &"}' | bash

# tail -n +2 tissue_dirs/$tissue/Hs_samples_info_$tissue.csv | cut -d, -f2,3 | sort | uniq | awk -F"," -v tissue=$tissue '{gsub(/"/, "", $0); print "R --slave --file=GEO_meta_pipeline.R --args  "$1" "$2" Hs_samples_final.csv tissue_dirs/"tissue" "tissue" > tissue_dirs/"tissue"/out1."$1"."$2".txt 2> tissue_dirs/"tissue"/out2."$1"."$2".txt &"}'

source('~/Dropbox (Gladstone)/programming/R/GEO_meta_pipeline_functions.R')
con <- dbConnect(SQLite(), '~/databases/GEOmetadb.sqlite')

GSE <- commandArgs(TRUE)[1]
GPL <- commandArgs(TRUE)[2]
if (length(commandArgs(TRUE)) > 2) {
    samples_info_file <- commandArgs(TRUE)[3]
} else {
    samples_info_file <- 'Hs_samples_final.csv'
}

# Splitting tissue results "tissue_info"
# and changing results directory
if (length(commandArgs(TRUE)) > 3) {
    res_dir <- commandArgs(TRUE)[4]
} else {
    res_dir <- ''
}
#res_dir = 'tissue_dirs'


# load series and platform data from GEO - taken from GEO2R
PD_label <- "PD" #"Parkinson's disease"
control_label <- "Control"
pheno_palette <- c("lightblue", "coral1")
Hs_samples_final <- read.csv(samples_info_file)

# plot of samples
#p <- ggplot(Hs_samples_final, aes(factor(gse))) #aes(x=PC1, y=PC2)) +
#p + geom_bar()

# Removing samples without phenotype
Hs_samples_final <- Hs_samples_final[!is.na(Hs_samples_final$pheno_info), ]

# qplot(main='All microarray studies', factor(gse), data=Hs_samples_final, geom="bar", fill=factor(tissue_info)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     labs(fill='Tissues', title='All microarray studies', x='GEO studies')
# ggsave('all_microarray_studies_by_tissue.png')
# 
# ggplot(Hs_samples_final, aes(factor(gse), fill=factor(pheno_info, levels = c('PD', 'Control')))) +
#     geom_bar(position="dodge") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     labs(fill='Phenotypes', title='All microarray studies', x='GEO studies')
# ggsave('all_microarray_studies_by_pheno.png')

#GSE_GPL <- unique(Hs_samples_final[,c(2,3,5)])

calculate_LIMMA <- TRUE
calculate_RP    <- TRUE
calculate_SAMR  <- TRUE
#boxplot_unlogged <- FALSE

GPL_title <- get_gpl_title(GPL)
CDF       <- get_gpl_cdf(GPL)

if (length(commandArgs(TRUE)) > 4) {
    analyzed_tissue <- gsub('_', ' ', commandArgs(TRUE)[5])
    samples_cur_annot <- Hs_samples_final[which(Hs_samples_final$gpl == GPL & Hs_samples_final$tissue_info == analyzed_tissue & Hs_samples_final$gse == GSE), c('gsm', 'gsm_title', 'tissue_info', 'tissue_info_general', 'pheno_info')]
} else {
    analyzed_tissue <- ''
    samples_cur_annot <- Hs_samples_final[which(Hs_samples_final$gpl == GPL & Hs_samples_final$gse == GSE), c('gsm', 'gsm_title', 'tissue_info', 'tissue_info_general', 'pheno_info')]
}

if (sum(is.na(samples_cur_annot$pheno_info)) > 0) {
    cat("The following samples will be removed because of lack of phenotype\n")
    for (gsm in samples_cur_annot$gsm[is.na(samples_cur_annot$pheno_info)]) {
        cat(gsm, '\n')
    }
    cat('\n\n')
    
    samples_cur_annot <- samples_cur_annot[!is.na(samples_cur_annot$pheno_info), ]
    #samples_cur_annot <- samples_cur_annot[-which(samples_cur_annot$gsm %in% Hs_samples_final$gsm[grep('neurological disease control', as.character(Hs_samples_final$gsm_title))]), ]
}

cat('\nDownloading/loading dataset...\n')


# If it's Affymetrix, always download the raw data (why not with other platforms?)
if (length(grep('affy', GPL_title, ignore.case = TRUE)) == 1 && !(GPL %in% c('GPL5175'))) {
    dir_name <- paste('CEL_files', GSE, sep='/')
    #dir_name <- paste('AffyFiles', GSE, sep='/')
    if (!dir.exists(dir_name)) {
        # Download Supp. files from GEO
        try(SupFiles <- getGEOSuppFiles(GSE, baseDir = 'CEL_files'), silent = TRUE)
        # Decompress...
        untar(tarfile=paste(dir_name, '/', GSE, '_RAW.tar', sep=''), exdir=dir_name)
        # ...and remove original RAW files to save space
        file.remove(paste(dir_name, dir(dir_name, pattern = 'RAW'), sep='/'))
    }
    
    # Searching for samples in CEL files
    all_cel_files <- dir(dir_name, pattern = 'CEL.gz')
    cel_files_cur_study <- NULL
    cur_pheno <- NULL
    for (s in samples_cur_annot$gsm) {
        cel_file_name <- paste(s, 'CEL.gz', sep='.')
        if (cel_file_name %in% all_cel_files) {
            cat('File', cel_file_name, 'found.\n')
            cel_files_cur_study <- c(cel_files_cur_study, cel_file_name)
        } else {
            if (length(grep(s, all_cel_files)) == 1) {
                cel_file_name <- grep(s, all_cel_files, value=TRUE)
                cel_files_cur_study <- c(cel_files_cur_study, cel_file_name)
                cat('File', cel_file_name, 'found.\n')
            } else if (length(grep(s, all_cel_files)) == 0) {
                cat('Sample', s, 'was NOT found!\n')
            } else {
                cat('Sample', s, 'was found multiple times!\n')
            }
        }
    }
    if (is.null(CDF)){
        esetOrig <- justRMA(filenames = paste(dir_name, '/', cel_files_cur_study, sep=''), normalize=FALSE, compress=TRUE)
    } else {
        esetOrig <- justRMA(filenames = paste(dir_name, '/', cel_files_cur_study, sep=''), normalize=FALSE, compress=TRUE, cdfname=CDF)
    }
    
    sampleNames(esetOrig) <- samples_cur_annot$gsm

} else {
    try(esetOrig <- getGEO(GSE, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = 'series_matrix'), silent = TRUE)
}

if (is.list(esetOrig)) { # && length(esetOrig)==1) { # ==1 or > 1?
    esetOrig = esetOrig[[1]]
}

# Treating expression tables
ex <- exprs(esetOrig)
ex <- ex[, which(colnames(ex) %in% samples_cur_annot$gsm)]

#esetFData <- data.frame(fData(esetOrig)) # doublecheck if this is working well
#rownames(ex) <- esetFData$Gene.ID
phenotypes <- samples_cur_annot$pheno_info
pheno_table <- table(phenotypes)
n_PD       <- pheno_table[PD_label]
n_controls <- pheno_table[control_label]
#colnames(ex) <- apply(samples_cur_annot[ ,c('pheno_info','gsm')], 1, function(x) {paste(x[1], x[2], sep='_')}) # fix this line
pheno_colors <- rep(NA, nrow(samples_cur_annot))
pheno_colors <- ifelse(samples_cur_annot$pheno_info == control_label, pheno_palette[1], pheno_palette[2])
#t_colors <- rep(NA, nrow(samples_cur_annot))

#tissue_palette <- brewer.pal(nlevels(samples_cur_annot$tissue_info), 'BrBG')
# TO DO: change colors to be fix, not random (to avoid similar colors)

tissue_palette <- distinctColorPalette(nlevels(samples_cur_annot$tissue_info))
#tissue_palette <- c("#DD62C9", "#DD9D51", "#CCA192", "#8478D7", "#D1DDA5", "#DC5E75", "#7DDAD2", "#81E659",
#                    "#7BA9D4", "#78DC98", "#A944E4", "#D6DD5A", "#D6D9DD", "#D5A0D1")
tissue_colors <- tissue_palette[samples_cur_annot$tissue_info]

column_annotation <- data.frame(tissue_colors, pheno_colors, stringsAsFactors = FALSE)
samples_cur_annot <- data.frame(samples_cur_annot, column_annotation)

probe_names <- rownames(ex)


cat('\nGetting annotation...\n')


## ANNOTATION
if (esetOrig@annotation == 'hgu133a') {
    probes_to_genes <- AnnotationDbi::select(hgu133a.db, keys = probe_names, keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL"))
} else if (esetOrig@annotation == 'hgu133a2') {
    probes_to_genes <- AnnotationDbi::select(hgu133a2.db, keys = probe_names, keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL"))
} else if (esetOrig@annotation == 'hgu133plus2') {
    probes_to_genes <- AnnotationDbi::select(hgu133plus2.db, keys = probe_names, keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL"))
} else if (esetOrig@annotation %in% c('hugene10sthsentrezgcdf')) {
    probes_to_genes <- AnnotationDbi::select(hugene10stv1hsentrezg.db, keys = probe_names, keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL"))
} else if (esetOrig@annotation == 'illuminaHumanv3') {
    probes_to_genes <- AnnotationDbi::select(illuminaHumanv3, keys = probe_names, keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL"))
} else if (esetOrig@annotation %in% c('GPL6480', 'GPL6104', 'GPL7884', 'GPL5175', 'GPL21436', 'GPL6947')) {
    # Probes to genes
    featureDf = pData(featureData(esetOrig))
    if (esetOrig@annotation %in% c('GPL6480', 'GPL6104', 'GPL6947')) {
        probes_to_genes <- featureDf[, c('ID', 'Gene ID', 'Gene symbol')]
    } else if (esetOrig@annotation %in% c('GPL7884')) {
        probes_to_genes <- featureDf[, c('ID', 'GENE_ID', 'GENE_SYMBOL')]
    } else if (esetOrig@annotation %in% c('GPL21436')) {
        probes_to_genes <- featureDf[, c('ID', 'GeneID', 'GeneID')]
        #probes_to_genes$ID <- NA
    } else if (esetOrig@annotation %in% c('GPL5175')) {
        SYMBOL <- apply(featureDf, 1, function(x) {strsplit(x[10], ' // ')[[1]][2]})
        ENTREZID <- rep(NA, length(SYMBOL))
        PROBEID <- names(SYMBOL)
        probes_to_genes <- data.frame(PROBEID, ENTREZID, SYMBOL)
    } else if (esetOrig@annotation %in% c('GPL21436')) {
        SYMBOL <- featureDf$GeneID
        ENTREZID <- rep(NA, length(SYMBOL))
        PROBEID <- featureDf$ID
        probes_to_genes <- data.frame(PROBEID, ENTREZID, SYMBOL)
    }
    colnames(probes_to_genes) <- c('PROBEID', 'ENTREZID', 'SYMBOL')
    #idx = sapply(featureDf, is.factor)
    #featureDf[idx] = lapply(featureDf[idx], as.character)
    #probe_names <- rownames(ex)
    #mapping = get_probe_mapping(GPL, featureDf)
}

cat('\nCreating plots and running analyses...\n')

pdf(paste(res_dir, paste('plots/plots_', GSE, '_', GPL, '.pdf', sep=''), sep='/'))

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste('GEO Series: ', GSE, '\n\nGEO Platform: ', GPL,
                             '\n', gsub('\\(','\n\\(', gsub('\\[','\n\\[', GPL_title)), '\n\n',
                             'Samples: ',n_PD,' PD x ',n_controls,' controls\n\n',
                             analyzed_tissue, sep=''),
     cex = 1.6, col = "black")

qx_values <- c(0., 0.25, 0.5, 0.75, 0.99, 1.0)
qx <- as.numeric(quantile(ex, qx_values, na.rm=T))
plot(qx, main='quantiles', xaxt="n", xlab='quantiles', ylab='values')
axis(1, at=1:6, labels=qx_values) #, col.axis="red", las=2)
needs_to_log <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

counts <- table(as.character(samples_cur_annot$pheno_info), as.character(samples_cur_annot$tissue_info))
if (ncol(counts) > 1) {
    par(las=2, mar=c(11,4.1,4.1,2.1))
    barplot(counts, beside = TRUE, legend = row.names(counts), col=pheno_palette, main=paste(GSE, '- phenotypes by tissue'))
    par(mar=c(5.1,4.1,4.1,2.1))
}

# Decide whether log2 transformation has to be applied based on quantiles
if (!needs_to_log) {
    if (ncol(counts) > 1) {
        boxplot(ex, main=paste(GSE, '(raw signal -',n_PD,'PD x',n_controls,'controls)\n', gsub('\\(','\n\\(', GPL_title)),
                col=tissue_colors, xaxt='n', mar=c(13,4.1,4.1,2.1), ylim=c(min(na.omit(ex))-0.5, max(na.omit(ex))+3))
        axis(1, at=1:ncol(ex), labels=colnames(ex), las=2)
        legend("topright", legend=unique(as.character(samples_cur_annot$tissue_info)), cex=1.2, fill=c(unique(samples_cur_annot$tissue_colors)))
    }
    boxplot(ex, main=paste(GSE, '(raw signal -',n_PD,'PD x',n_controls,'controls)\n', gsub('\\(','\n\\(', GPL_title), '\n', analyzed_tissue),
            col=samples_cur_annot$pheno_colors, xaxt='n', mar=c(13,4.1,4.1,2.1), ylim=c(min(na.omit(ex))-0.5, max(na.omit(ex))+3))
    axis(1, at=1:ncol(ex), labels=colnames(ex), las=2)
    legend("topright", legend=unique(samples_cur_annot$pheno_info), cex=1.2, fill=unique(samples_cur_annot$pheno_colors))
} else {
    cat('Logging (base 2) data.\n')
    ex[which(ex <= 0)] <- NaN
    # exprs(esetOrig) <- log2(ex)
    # ex <- exprs(esetOrig)
    ex <- log2(ex)
    qx <- as.numeric(quantile(ex, qx_values, na.rm=T))
    boxplot(ex, main=paste(GSE, '(log2 signal -',n_PD,'PD x',n_controls,'controls)\n', gsub('\\(','\n\\(', GPL_title), '\n', analyzed_tissue),
            col=pheno_colors, xaxt='n')
    axis(1, at=1:ncol(ex), labels=colnames(ex), las=2)
    legend("topright", c(control_label, PD_label), cex=1.2, fill=c("lightblue", "coral1"))
    plot(qx, main='quantiles', xaxt="n", xlab='quantiles', ylab='values')
    axis(1, at=1:6, labels=qx_values)
}

if (GSE == 'GSE20186' && GPL == 'GPL96' && analyzed_tissue == 'substantia nigra') {
    # ComBat (ran just for GSE20186, GPL96, substantia nigra)
    batch <- rep(1, nrow(samples_cur_annot))
    batch[grep('HS', samples_cur_annot$gsm_title)] <- 2
    modcombat = model.matrix(~1, data=samples_cur_annot)
    combat_edata = ComBat(dat=ex, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    
    boxplot(combat_edata, main=paste(GSE, '(correction with ComBat -',n_PD,'PD x',n_controls,'controls)\n',
                                gsub('\\(','\n\\(', GPL_title), '\n', analyzed_tissue), col=pheno_colors, xaxt='n', ylim=c(min(na.omit(ex))-0.5, max(na.omit(ex))+3))
    axis(1, at=1:ncol(combat_edata), labels=colnames(combat_edata), las=2)
    legend("topright", legend=unique(samples_cur_annot$pheno_info), cex=1.2, fill=unique(samples_cur_annot$pheno_colors))
    
    
    # modcancer = model.matrix(~, data=samples_cur_annot)
    # modcancer = model.matrix(~cancer, data=pheno)
    # combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    #ex <- 
    
    # Normalization
    norm_ex <- normalizeBetweenArrays(combat_edata, method="cyclicloess") # the same as normalizeCyclicLoess(ex)
} else {
    norm_ex <- normalizeBetweenArrays(ex, method="cyclicloess") # the same as normalizeCyclicLoess(ex)
}
# if (length(grep('affy', GPL_title, ignore.case = TRUE)) == 1) {
#     norm_ex <- normalizeBetweenArrays(ex, method="cyclicloess", cyclic.method = "affy")
# } else {
#     norm_ex <- normalizeBetweenArrays(ex, method="cyclicloess")
# }

boxplot(norm_ex, main=paste(GSE, '(normalized signal -',n_PD,'PD x',n_controls,'controls)\n',
        gsub('\\(','\n\\(', GPL_title), '\n', analyzed_tissue), col=pheno_colors, xaxt='n', ylim=c(min(na.omit(ex))-0.5, max(na.omit(ex))+3))
axis(1, at=1:ncol(norm_ex), labels=colnames(norm_ex), las=2)
legend("topright", legend=unique(samples_cur_annot$pheno_info), cex=1.2, fill=unique(samples_cur_annot$pheno_colors))

exprsByGene <- norm_ex

# Removing NA rows
na_count <- apply(exprsByGene, 1, function(x) {sum(is.na(x))})
exprsByGene <- exprsByGene[na_count < (ncol(exprsByGene)/2), ]

# Calculate expression by gene
# if (esetOrig@annotation %in% c('GPL6480')) {
#     exprsByGene = calcExprsByGene(norm_ex, mapping)
# } else {
#     exprsByGene <- norm_ex
# }
#colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

plot_hclust(exprsByGene, samples_cur_annot, 'all genes', 'topright')

plot_PCA_tSNE(exprsByGene, samples_cur_annot, 'all genes', 'topleft')

############################
# FOLD CHANGE CALCULATION
#logFC <- apply(exprsByGene, 1, function(x) {mean(x[phenotypes==PD_label] - mean(x[phenotypes==control_label]))})
#logFC <- data.frame(names(logFC), logFC)
#colnames(logFC)[1] <- 'probeId'

# List of all genes
if (file.exists(paste(res_dir, paste(GSE, GPL, 'all_genes_simple.csv', sep='_'), sep='/'))) {
    all.genes.simple <- read.csv(paste(res_dir, paste(GSE, GPL, 'all_genes_simple.csv', sep='_'), sep='/'))
    all.genes.limma <- all.genes.simple[all.genes.simple$Method=='LIMMA', ]
    all.genes.RP    <- all.genes.simple[all.genes.simple$Method=='RankProd', ]
    all.genes.SAMR  <- all.genes.simple[all.genes.simple$Method=='SAMR', ]
} else {
    all.genes.simple <- NULL
    all.genes.limma  <- NULL
    all.genes.RP     <- NULL
    all.genes.SAMR   <- NULL
}

# List of DE genes
# if (file.exists(paste(GSE, GPL, 'dif_exp_simple.csv', sep='_'))) {
#     dif_exp_simple <- read.csv(paste(GSE, GPL, 'dif_exp_simple.csv', sep='_'))
#     dif_exp_RP    <- dif_exp_simple[dif_exp_simple$Method=='RankProd', ]
#     dif_exp_limma <- dif_exp_simple[dif_exp_simple$Method=='LIMMA', ]
# } else {
#     dif_exp_RP <- NULL
#     dif_exp_limma <- NULL
# }

############################
# Dif. expression analysis - LIMMA
if (calculate_LIMMA) {
    
    if (is.null(all.genes.limma)) {
        result_limma <- run_LIMMA(GSE, exprsByGene, samples_cur_annot, probes_to_genes) # this will return a list with 'all.genes' and 'top.genes'
        all.genes.limma <- result_limma$all.genes
        all.genes.limma <- merge(probes_to_genes, all.genes.limma, by=c(1,1), all.y=TRUE)
        all.genes.limma <- all.genes.limma[, c('PROBEID', 'ENTREZID', 'SYMBOL', 'logFC', 'Adj.P.Value')]
        all.genes.limma <- data.frame(all.genes.limma, 'LIMMA')
        colnames(all.genes.limma)[ncol(all.genes.limma)] <- 'Method'
    }

    plot(main=paste("Volcano plot - ", GSE, "\nAll genes (LIMMA)", sep=""), all.genes.limma$logFC, -log10(all.genes.limma$Adj.P.Value),
         xlab='log2(Fold Change)', ylab='-log10(Adj. P-value)', pch=20)
    abline(v = -1, col = 'grey')
    abline(v =  1, col = 'grey')
    abline(h = -log10(0.05), col = 'grey')
    with(subset(all.genes.limma, Adj.P.Value<.05 & logFC >= 1), points(logFC, -log10(Adj.P.Value), pch=20, col="red"))
    with(subset(all.genes.limma, Adj.P.Value<.05 & logFC <=-1), points(logFC, -log10(Adj.P.Value), pch=20, col="green"))
    with(subset(all.genes.limma, Adj.P.Value<.05 & abs(logFC)>1), textxy(logFC, -log10(Adj.P.Value), labs=SYMBOL, cex=.8))
    
}

############################
# Dif. expression analysis - RankProd
if (calculate_RP) {
    
    if (is.null(all.genes.RP)) {
        result_RP <- run_RankProd(GSE, exprsByGene, samples_cur_annot, probes_to_genes) # this will return a list with 'all.genes' and 'top.genes'
        all.genes.RP <- result_RP$all.genes
        all.genes.RP <- merge(probes_to_genes, all.genes.RP, by=c(1,1), all.y=TRUE)
        all.genes.RP <- all.genes.RP[, c('PROBEID', 'ENTREZID', 'SYMBOL', 'logFC', 'Adj.P.Value')]
        all.genes.RP <- data.frame(all.genes.RP, 'RankProd')
        colnames(all.genes.RP)[ncol(all.genes.RP)] <- 'Method'

    }
    
    plot(main=paste("Volcano plot - ", GSE, "\nAll genes (RankProd)", sep=""), all.genes.RP$logFC, -log10(all.genes.RP$Adj.P.Value),
         xlab='log2(Fold Change)', ylab='-log10(Adj. P-value)', pch=20) # Adj.P.Value is the PFP measure for RankProd
    abline(v = -1, col = 'grey')
    abline(v =  1, col = 'grey')
    abline(h = -log10(0.05), col = 'grey')
    with(subset(all.genes.RP, Adj.P.Value<.05 & logFC >= 1), points(logFC, -log10(Adj.P.Value), pch=20, col="red"))
    with(subset(all.genes.RP, Adj.P.Value<.05 & logFC <=-1), points(logFC, -log10(Adj.P.Value), pch=20, col="green"))
    with(subset(all.genes.RP, Adj.P.Value<.05 & abs(logFC)>1), textxy(logFC, -log10(Adj.P.Value), labs=SYMBOL, cex=.8))

}

############################
# Dif. expression analysis - SAMR
if (calculate_SAMR) {
    
    if (is.null(all.genes.SAMR)) {
        result_SAMR <- run_SAMR(GSE, exprsByGene, samples_cur_annot, probes_to_genes) # this will return a list with 'all.genes' and 'top.genes'
        all.genes.SAMR <- result_SAMR$all.genes
        all.genes.SAMR <- merge(probes_to_genes, all.genes.SAMR, by=c(1,1), all.y=TRUE)
        all.genes.SAMR <- all.genes.SAMR[, c('PROBEID', 'ENTREZID', 'SYMBOL', 'logFC', 'Adj.P.Value')]
        all.genes.SAMR <- data.frame(all.genes.SAMR, 'SAMR')
        colnames(all.genes.SAMR)[ncol(all.genes.SAMR)] <- 'Method'

    }
    
    plot(main=paste("Volcano plot - ", GSE, "\nAll genes (SAMR)", sep=""), all.genes.SAMR$logFC, -log10(all.genes.SAMR$Adj.P.Value),
         xlab='log2(Fold Change)', ylab='-log10(Adj. P-value)', pch=20)
    abline(v = -1, col = 'grey')
    abline(v =  1, col = 'grey')
    abline(h = -log10(0.05), col = 'grey')
    with(subset(all.genes.SAMR, Adj.P.Value<.05 & logFC >= 1), points(logFC, -log10(Adj.P.Value), pch=20, col="red"))
    with(subset(all.genes.SAMR, Adj.P.Value<.05 & logFC <=-1), points(logFC, -log10(Adj.P.Value), pch=20, col="green"))
    with(subset(all.genes.SAMR, Adj.P.Value<.05 & abs(logFC)>1), textxy(logFC, -log10(Adj.P.Value), labs=SYMBOL, cex=.8))

}


if (is.null(all.genes.simple)) {
    cat('Table all.genes.simple is NULL. Saving results.\n')
    all.genes.simple <- rbind(all.genes.limma, all.genes.RP)
    all.genes.simple <- rbind(all.genes.simple, all.genes.SAMR)
    write.csv(all.genes.simple, paste(res_dir, paste(GSE, GPL, 'all_genes_simple.csv', sep='_'), sep='/'), row.names = FALSE, quote = FALSE)
} else {
    cat('Table all.genes.simple exists (it was saved in the past).\n')
}

all.genes <- compare_dif_exp_genes(GSE, GPL, all.genes.limma, all.genes.RP, all.genes.SAMR, analyzed_tissue)
write.csv(all.genes, paste(res_dir, paste(GSE, GPL, 'all_genes.csv', sep='_'), sep='/'), row.names = FALSE, quote = FALSE)

#grid.newpage(venn.diagram(main="Down-regulated genes", filename=NULL,
#             x = list("LIMMA" = down_genes_limma,"RankProd" = down_genes_RP),
#             category = c("LIMMA", "RankProd"), cat.pos = c(0,0)))

# library("mygene")
# gene_list <- rownames(up.genes)  # c('DDX26B','CCDC83', 'MAST3')
# a <- queryMany(gene_list, scopes="entrezgene", fields=c("ensembl.gene", "reporter"), species="human")
# geneList <- rep(1, nrow(up.genes))
# names(geneList) <- row.names(up.genes)
# 
# all.genes.entrez_id <- ls(hgu133plus2ACCNUM)
# all.genes.entrez_id <- queryMany(all.genes, scopes="reporter", fields=c("ensembl.gene", "entrezgene"), species="human")
# 
# sampleGOdata <- new("topGOdata",
# description = paste(GSE, "GO Biological process", sep = ' - '),
# ontology = "BP", allGenes = geneList,
# geneSel = topDiffGenes, nodeSize = 10,
# annot = annFUN.GO2genes, affyLib = affyLib)

# plot DE genes
#plotRP(RP.out, cutoff=0.05)

#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

#heatmap.2(DEtable, ColSideColors=p_colors, dendrogram="col", main = paste('Heatmap (dif. exp. genes) -', GSE))

dev.off()
