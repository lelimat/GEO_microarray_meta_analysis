# GEO_meta_pipeline_functions.R
# C: Sep 15, 2016
# M: Feb 16, 2017
# A: Leandro Lima <leandro.lima@gladstone.ucsf.edu>

installation_needed <- FALSE

if (installation_needed) {
    
    ## try http:// if https:// URLs are not supported
    source("https://bioconductor.org/biocLite.R")
    biocLite("GEOmetadb")
    biocLite("GEOquery")
    biocLite("affy")
    biocLite("limma")
    install.packages("gplots")
    biocLite("RankProd")
    biocLite("topGO")
    biocLite("hgu133a.db")
    biocLite("hgu133a2.db")
    biocLite("hgu133plus2.db")
    biocLite("huex10stv2cdf")
    biocLite("hugene10stv1cdf")
    biocLite("hugene10stprobeset.db")
    biocLite("sva")
    biocLite("samr")
    install.packages("devtools")
    install.packages('randomcoloR')
    install.packages("Rtsne")
    install.packages('dendextend')
    install.packages('ggfortify')
    install.packages('calibrate')
    install.packages('VennDiagram')
    
    download.file('http://mbni.org/customcdf/21.0.0/entrezg.download/hugene10sthsentrezgcdf_21.0.0.tar.gz', destfile = 'CDF/hugene10sthsentrezgcdf_21.0.0.tar.gz')
    install.packages('CDF/hugene10sthsentrezgcdf_21.0.0.tar.gz', repos=NULL, type='source')
    
    download.file('http://mbni.org/customcdf/21.0.0/entrezg.download/hgu133ahsentrezgcdf_21.0.0.tar.gz', destfile = 'CDF/hgu133ahsentrezgcdf_21.0.0.tar.gz')
    install.packages('CDF/hgu133ahsentrezgcdf_21.0.0.tar.gz', repos = NULL, type = 'source')
    
    download.file('http://mbni.org/customcdf/21.0.0/entrezg.download/huex10sthsentrezgcdf_21.0.0.tar.gz', destfile = 'CDF/huex10sthsentrezgcdf_21.0.0.tar.gz')
    install.packages('CDF/huex10sthsentrezgcdf_21.0.0.tar.gz', repos = NULL, type = 'source')
    
    download.file('http://mbni.org/customcdf/15.0.0/entrezg.download/hugene10stv1hsentrezg.db_15.0.0.tar.gz', destfile = 'CDF/hugene10stv1hsentrezg.db_15.0.0.tar.gz')
    install.packages('CDF/hugene10stv1hsentrezg.db_15.0.0.tar.gz', repos = NULL, type = 'source')
    
}

library("affy")
library("annotate")
library("Biobase")
library('calibrate')
library('dendextend')
library("devtools")
library("GEOmetadb")
library("GEOquery")
library('ggfortify')
library("gplots")
library("hgu133a.db")
library("hgu133a2.db")
library("hgu133plus2.db")
library("hugene10stv1hsentrezg.db")
library('huex10sthsentrezgcdf')
#library("huex10stv2cdf")
library("limma")
library('org.Hs.eg.db')
library('randomcoloR')
library("RankProd")
library("Rtsne")
library("sva")
library("samr")
#library("topGO")
library('VennDiagram')


# Special heatmap
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

get_gpl_title = function(GPL) {
    gpl_title <- dbGetQuery(con, paste("select title from gpl where gpl = '",GPL,"'", sep=''))
    if (is.na(gpl_title)) {
        return(gpl_title)
    } else {
        return(gpl_title[1,1])
    }
}

# Source: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/21.0.0/entrezg.asp
get_gpl_cdf = function(GPL) {
    bioc_package <- dbGetQuery(con, paste("select bioc_package from gpl where gpl = '",GPL,"'", sep=''))
    #cdf_link <- strsplit(as.character(sup_files), ';')[[1]][1]
    
    # cdf_file <- strsplit(cdf_link, '/')[[1]][length(strsplit(cdf_link, '/')[[1]])]
    # if(!file.exists(paste('CDF', cdf_file, sep='/'))) {
    #     download.file(cdf_link, destfile = paste('CDF', cdf_file, sep='/'))
    # }
    # install.packages('hugene10sthsentrezgcdf')
    
    if (is.na(bioc_package)) {
        if (GPL %in% c('GPL17047')) {
            return('hugene10sthsentrezgcdf') #return('hugene10stv1cdf')
        } else if (GPL %in% c('GPL5175')) {
            return('huex10sthsentrezgcdf')
        }
    } else {
        return(bioc_package[1,1])
    }
}

########################
### Downloading data ###
########################

plot_hclust = function(exprsByGene, samples_info, subtitle, legend_pos = 'topright') {
    # Hierarchical Clustering
    dist_method <- "euclidean"
    exp_transp <- t(na.omit(exprsByGene))
    d <- dist(exp_transp, method = dist_method)
    hc <- hclust(d)
    dend <- as.dendrogram(hc)
    labels_colors(dend) <- samples_info$pheno_colors[order.dendrogram(dend)]
    plot(dend, main=paste(GSE, ' - Hierarchical cluster (dist. method = ',dist_method,')\n', subtitle, ' (by phenotype)', sep=''))
    legend(legend_pos, legend=unique(as.character(samples_info$pheno_info)),
           fill = c(unique(as.character(samples_info$pheno_colors))), border = FALSE)
    
    if (length(unique(samples_info$tissue_info)) > 1) {
        labels_colors(dend) <- samples_info$tissue_colors[order.dendrogram(dend)]
        plot(dend, main=paste(GSE, ' - Hierarchical cluster (dist. method = ',dist_method,')\n', subtitle, 'all genes (by tissue)', sep=''))
        legend(legend_pos, legend=unique(as.character(samples_info$tissue_info)),
               fill = c(unique(as.character(samples_info$tissue_colors))), border = FALSE)
    }
}

plot_PCA_tSNE = function(exprsByGene, samples_info, subtitle, legend_pos = 'topright') {
    # GSM_colors  = column_annotation[,2]
    # GSM_tissues = samples_cur_annot$tissue_info
    # subtitle    = 'all genes'
    
    GSM_list <- colnames(exprsByGene)[-1]
    exp_transp <- t(na.omit(exprsByGene))
    for (p in c(5, 10)) {
        if (nrow(exp_transp) - 1 >= 3 * p) {
            rtsne_out <- Rtsne(exp_transp, perplexity = p)
            # PHENOTYPE
            plot(rtsne_out$Y, t='n', main=paste(GSE, " - BarnesHutSNE\n(",subtitle," - perplexity: ",p,")", sep=''))
            text(rtsne_out$Y, labels=data.frame(strsplit(GSM_list, 'GSM'), stringsAsFactors = FALSE)[2,], col = samples_info$pheno_colors)
            legend(legend_pos, legend=unique(as.character(samples_info$pheno_info)),
                   fill = c(unique(as.character(samples_info$pheno_colors))), border = FALSE)
            # TISSUE
            if (length(unique(samples_info$tissue_info)) > 1) {
                plot(rtsne_out$Y, t='n', main=paste(GSE, " - BarnesHutSNE\n(",subtitle," - perplexity: ",p,")", sep=''))
                text(rtsne_out$Y, labels=data.frame(strsplit(GSM_list, 'GSM'), stringsAsFactors = FALSE)[2,], col = samples_info$tissue_colors)
                legend(legend_pos, legend=unique(as.character(samples_info$tissue_info)),
                       fill = c(unique(as.character(samples_info$tissue_colors))), border = FALSE)
            }
        }
    }
    
    # Plotting PCA
    pca <- prcomp(t(na.omit(exprsByGene)))
    #autoplot(pca, colour=samples_info$pheno_colors)
    plot(pca$x[,1:2], t='n', main=paste(GSE, " - PCA\n(",subtitle,")", sep=''))
    text(pca$x[,1:2], labels=data.frame(strsplit(GSM_list, 'GSM'), stringsAsFactors = FALSE)[2,],
         col=samples_info$pheno_colors, pch=16, main=paste(GSE, " - PCA\n(",subtitle,")", sep=''))
    legend(legend_pos, legend=unique(as.character(samples_info$pheno_info)),
           fill = c(unique(as.character(samples_info$pheno_colors))), border = FALSE)
    
    if (length(unique(samples_info$tissue_info)) > 1) {
        plot(pca$x[,1:2], t='n', main=paste(GSE, " - PCA\n(",subtitle,")", sep=''))
        text(pca$x[,1:2], labels=data.frame(strsplit(GSM_list, 'GSM'), stringsAsFactors = FALSE)[2,],
             col=samples_info$tissue_colors, pch=16, main=paste(GSE, " - PCA\n(",subtitle,")", sep=''))
        legend(legend_pos, legend=unique(as.character(samples_info$tissue_info)),
               fill = c(unique(as.character(samples_info$tissue_colors))), border = FALSE)
    }
    
}

run_LIMMA = function(GSE, exprsByGene, samples_info, probes_to_genes = NULL) {
    # TO DO: return full list of genes instead of only DEG
    # run_LIMMA(GSE, ex_tissue, samples_info)
    # samples_info <- samples_cur_annot
    
    if (length(unique(samples_info$tissue_info)) == 1) {
        tissue <- unique(samples_info$tissue_info)[1]
        GSE <- paste(GSE, tissue, sep=' - ')
    }
    cat("Running LIMMA.\n")
    phenotypes <- samples_info$pheno_info
    if (sum(phenotypes==control_label) >= 5 && sum(phenotypes==PD_label) >= 5) {
        design <- model.matrix(~ phenotypes + 0, as.data.frame(exprsByGene))
        colnames(design) <- levels(phenotypes)
        fit <- lmFit(exprsByGene, design)
        cont.matrix <- makeContrasts(PD-Control, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, proportion=1)
        # The function below will not report B: "log-odds that the gene is differentially expressed (omitted for topTreat)"
        all.genes <- topTable(fit2, adjust="fdr", number=Inf) #, sort.by="B", number=10000, p.value=0.05)
        all.genes <- data.frame(rownames(all.genes), all.genes)
        colnames(all.genes)[1] <- 'probeId'
        colnames(all.genes)[colnames(all.genes)=='adj.P.Val'] <- 'Adj.P.Value'
        #all.genes$Adj.P.Value <- as.numeric(all.genes$Adj.P.Value)
        
        top.genes <- topTable(fit2, adjust="fdr", sort.by="B", number=10000, p.value=0.05)
        top.genes <- data.frame(rownames(top.genes), top.genes)
        colnames(top.genes)[1] <- 'probeId'
        colnames(top.genes)[colnames(top.genes)=='adj.P.Val'] <- 'Adj.P.Value'
        
        #top.genes$Adj.P.Value <- as.numeric(top.genes$Adj.P.Value)
        cat("Saving results.\n")
        # if (nrow(top.genes) > 0) {
        #     write.csv(top.genes,   paste(GSE, 'DE_genes_LIMMA.csv',   sep='_'), row.names=FALSE)
        # } else {
        #     write.csv(top.genes,   paste(GSE, 'zero_DE_genes_LIMMA.txt',   sep='_'), row.names=FALSE)
        # }
        
        if (nrow(top.genes) >= 2) {
            cat('nrow(top.genes)=',nrow(top.genes), '\n')
            cat('ncol(top.genes)=',ncol(top.genes), '\n')
            DEtable <- na.omit(exprsByGene[which(row.names(exprsByGene) %in% top.genes$probeId), ])
            top.genes_na.omit <- top.genes[rownames(DEtable),]
            DE_colors <- matrix('grey', ncol=nrow(DEtable))
            DE_colors[1, which(top.genes_na.omit$logFC > 0)] <- 'red'   # up
            DE_colors[1, which(top.genes_na.omit$logFC < 0)] <- 'green' # down
            #DE_colors[1, which(is.na(top.genes_na.omit$logFC))] <- 'grey'
            cat("Plotting heatmap.\n")
            # https://www.biostars.org/p/18211/
            heatmap.3(DEtable, main = paste('Heatmap (', nrow(DEtable), ' DE genes LIMMA) -', GSE, sep=''),
                      ColSideColors=as.matrix(samples_info[,c('tissue_colors', 'pheno_colors')]),
                      RowSideColors=DE_colors, labCol=samples_info$gsm, ColSideColorsSize=2)
            legend('topright', legend=c(PD_label, control_label, '', unique(as.character(samples_info$tissue_info))),
                   fill = c('coral1', 'lightblue', 'white', unique(samples_info$tissue_colors)), border = FALSE)
            
            plot_PCA_tSNE(DEtable, samples_info, paste(nrow(DEtable), 'DE genes LIMMA'))
            
        }
        
        return(list('top.genes'=top.genes, 'all.genes'=all.genes))
        
    } else {
        
        cat("Not enough samples to run LIMMA (it's necessary at least 5 per group).\n")
        return(NULL)
        
    }
}

run_RankProd = function(GSE, exprsByGene, samples_info, probes_to_genes = NULL) {
    #run_RankProd(GSE, exprsByGene, samples_cur_annot)
    #run_RankProd(GSE, ex_tissue, samples_cur_tissue)
    # samples_info = samples_cur_annot
    # exprsByGene = ex_tissue
    # samples_info = samples_cur_tissue
    GSM_list <- samples_info$gsm
    if (length(unique(samples_info$tissue_info)) == 1) {
        tissue <- unique(samples_info$tissue_info)[1]
        GSE <- paste(GSE, tissue, sep=' - ')
    }
    phenotypes <- samples_info$pheno_info
    numeric_pheno <- rep(NA, length(phenotypes))
    # 0: control, 1: case
    
    if (sum(phenotypes==control_label) >= 5 && sum(phenotypes==PD_label) >= 5) {
        numeric_pheno[phenotypes==PD_label]      <- 0 #1
        numeric_pheno[phenotypes==control_label] <- 1 #0
        
        RP.out <- RankProd::RankProducts(exprsByGene, cl = numeric_pheno, gene.names = rownames(exprsByGene), na.rm = FALSE, logged = TRUE)
        
        all.genes <- data.frame(rownames(RP.out$AveFC), RP.out$AveFC, RP.out$pval)
        
        ## Playing with RankProd results
        # dot_colors <- rep('grey', nrow(RP.out$pfp))
        # dot_colors[RP.out$AveFC > 0.5]  <- 'salmon'     # up
        # dot_colors[RP.out$AveFC < -0.5]  <- 'lightgreen' # down
        # dot_colors[RP.out$AveFC > 1]  <- 'red'   # up
        # dot_colors[RP.out$AveFC < -1] <- 'green' # down
        # 
        # plot(RP.out$pfp[,1], RP.out$pfp[,2], col=dot_colors, pch=20, main='PFP of genes',
        #      xlab=paste('PFP :', colnames(RP.out$pfp))[1], ylab=paste('PFP :', colnames(RP.out$pfp))[2])
        # legend('topright', legend=c('up (logFC>0.5)', 'down (logFC<-0.5)'), col=c('salmon', 'lightgreen'), pch=20)
        # 
        # plot(RP.out$pval[,1], RP.out$pval[,2], col=dot_colors, pch=20, main='P-value of genes',
        #      xlab=paste('P-value :', colnames(RP.out$pval))[1], ylab=paste('P-value :', colnames(RP.out$pval))[2])
        # legend('topright', legend=c('up (logFC>0.5)', 'down (logFC<-0.5)'), col=c('salmon', 'lightgreen'), pch=20)
        # 
        # boxplot(exprsByGene, col=pheno_colors)
        # probeId = '205857_at'
        # boxplot(exprsByGene[probeId, ] ~ phenotypes, main=paste('RankProd - Good PFP, Good FC\n', probeId, '\nlogFC =', round(all.genes[probeId, 'logFC'], 2), '| PFP =', all.genes[probeId, 'Adj.P.Value']))
        # probeId = '212637_s_at'
        # boxplot(exprsByGene[probeId, ] ~ phenotypes, main=paste('RankProd - Good PFP, Bad FC\n', probeId, '\nlogFC =', round(all.genes[probeId, 'logFC'], 2), '| PFP =', all.genes[probeId, 'Adj.P.Value']))
        # probeId = '205769_at'
        # boxplot(exprsByGene[probeId, ] ~ phenotypes, main=paste('RankProd - Bad PFP, Bad FC\n', probeId, '\nlogFC =', round(all.genes[probeId, 'logFC'], 2), '| PFP =', all.genes[probeId, 'Adj.P.Value']))
        # probeId = '205113_at'
        # boxplot(exprsByGene[probeId, ] ~ phenotypes, main=paste('RankProd - Bad PFP?, Good FC?\n', probeId, '\nlogFC =', round(all.genes[probeId, 'logFC'], 2), '| PFP =', all.genes[probeId, 'Adj.P.Value']))
        
        colnames(all.genes) <- c('probeId', 'logFC', 'pval.down', 'pval.up')
        P.Value <- with(all.genes, ifelse(logFC > 0, pval.up, pval.down))
        Adj.P.Value <- p.adjust(P.Value, method = "bonferroni")
        all.genes <- data.frame(all.genes, Adj.P.Value)

        up.genes   <- subset(all.genes, logFC >=  1 & Adj.P.Value <= 0.5)
        down.genes <- subset(all.genes, logFC <= -1 & Adj.P.Value <= 0.5)
        
        # if (!is.null(up.genes) && !is.null(down.genes)) {
        if (nrow(up.genes)>0 && nrow(down.genes)) {
            top.genes <- rbind(up.genes, down.genes)
            # write.csv(top.genes, paste(GSE, 'DE_genes_RankProd.csv', sep='_'), row.names=FALSE)
            
            DEtable <- na.omit(exprsByGene[top.genes$probeId, ])
            top.genes_na.omit <- top.genes[rownames(DEtable), ]
            
            DE_colors <- matrix('grey', ncol=nrow(top.genes_na.omit))
            DE_colors[1, which(top.genes_na.omit$logFC > 0)] <- 'red'   # up
            DE_colors[1, which(top.genes_na.omit$logFC < 0)] <- 'green' # down
            
            # pfp = percentage of false predictions
            
            if (nrow(DEtable) >= 2) {
                
                # https://www.biostars.org/p/18211/
                heatmap.3(DEtable, main = paste('Heatmap (', nrow(DEtable), ' DE genes RankProd) - ', GSE, sep=''),
                          ColSideColors=as.matrix(samples_info[,c('tissue_colors', 'pheno_colors')]),
                          RowSideColors=DE_colors, labCol=samples_info$gsm, ColSideColorsSize=2)
                legend('topright', legend=c(PD_label, control_label, '', unique(as.character(samples_info$tissue_info))),
                       fill = c('coral1', 'lightblue', 'white', unique(samples_info$tissue_colors)), border = FALSE)
                
                plot_PCA_tSNE(DEtable, samples_info, 'DE genes RankProd')
            }
            
            return(list('top.genes'=top.genes, 'all.genes'=all.genes))
            
        } else {
            cat("No gene was found differentially expressed.\n")
            return(list('top.genes'=NULL, 'all.genes'=all.genes))
        }
        
    } else {
        
        cat("Not enough samples to run RankProd (it's necessary at least 5 per group).\n")
        return(NULL)
        
    }
}

run_SAMR = function(GSE, exprsByGene, samples_info, probes_to_genes = NULL) {
    # run_SAMR(GSE, exprsByGene, samples_cur_annot, probes_to_genes)
    GSM_list <- samples_info$gsm
    if (length(unique(samples_info$tissue_info)) == 1) {
        tissue <- unique(samples_info$tissue_info)[1]
        GSE <- paste(GSE, tissue, sep=' - ')
    }
    
    phenotypes <- samples_info$pheno_info
    numeric_pheno <- rep(NA, length(phenotypes))
    if (sum(phenotypes==control_label) >= 5 && sum(phenotypes==PD_label) >= 5) {
        numeric_pheno[phenotypes==control_label] <- 1
        numeric_pheno[phenotypes==PD_label]      <- 2

        data <- list(x=as.matrix(exprsByGene), y=numeric_pheno, geneid=rownames(exprsByGene), genenames=rownames(exprsByGene), logged2=TRUE)
        
        samr.obj <- samr(data, resp.type="Two class unpaired", assay.type='array', nperms=100)
        delta.table <- samr.compute.delta.table(samr.obj) #, min.foldchange=0.1, nvals=200)
        siggenes.table <- samr.compute.siggenes.table(samr.obj, del=0, data, delta.table, all.genes=TRUE)
        Adj.P.Value <- samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
        Adj.P.Value <- data.frame(names(Adj.P.Value), Adj.P.Value)
        all.genes <- rbind(siggenes.table$genes.lo, siggenes.table$genes.up)
        all.genes <- merge(data.frame(all.genes, stringsAsFactors = FALSE)[,-1], Adj.P.Value, by=c(1,1))
        logFC <- log2(as.numeric(all.genes$Fold.Change))
        all.genes <- data.frame(all.genes, logFC)
        rownames(all.genes) <- all.genes$Gene.Name
        
        # DE genes
        up.genes   <- all.genes[all.genes$logFC >=  1 & all.genes$Adj.P.Value <= 0.05, ]
        down.genes <- all.genes[all.genes$logFC <= -1 & all.genes$Adj.P.Value <= 0.05, ]
        
        top.genes <- rbind(up.genes, down.genes)
        colnames(top.genes)[1] <- 'probeId'
        
        #if (!is.null(up.genes) && !is.null(down.genes)) {
        if (nrow(up.genes)>0 && nrow(down.genes)>0) {
            
            #write.csv(top.genes, paste(GSE, 'DE_genes_SAMR.csv', sep='_'), row.names=FALSE)
            
            DEtable <- na.omit(exprsByGene[top.genes$probeId, ])
            top.genes_na.omit <- top.genes[rownames(DEtable), ]
            
            DE_colors <- matrix('grey', ncol=nrow(top.genes_na.omit))
            DE_colors[1, which(top.genes_na.omit$logFC > 0)] <- 'red'   # up
            DE_colors[1, which(top.genes_na.omit$logFC < 0)] <- 'green' # down
            
            # pfp = percentage of false predictions
            
            if (nrow(DEtable) >= 2) {
                
                # https://www.biostars.org/p/18211/
                heatmap.3(DEtable, main = paste('Heatmap (', nrow(DEtable), ' DE genes SAMR) -', GSE, sep=''),
                          ColSideColors=as.matrix(samples_info[,c('tissue_colors', 'pheno_colors')]),
                          RowSideColors=DE_colors, labCol=samples_info$gsm, ColSideColorsSize=2)
                legend('topright', legend=c(PD_label, control_label, '', unique(as.character(samples_info$tissue_info))),
                       fill = c('coral1', 'lightblue', 'white', unique(samples_info$tissue_colors)), border = FALSE)
                
                plot_PCA_tSNE(DEtable, samples_info, 'DE genes SAMR')
            }
            
            return(list('top.genes'=top.genes, 'all.genes'=all.genes))
            
        } else {
            
            cat("No gene was found differentially expressed.\n")
            return(list('top.genes'=NULL, 'all.genes'=all.genes))
        }
        
    } else {
        
        cat("Not enough samples to run SAMR (it's necessary at least 5 per group).\n")
        return(NULL)
        
    }

}

volcanoplot.2 = function(method.name = '', all.genes) {
    
    plot(main=paste("Volcano plot - ", GSE, "\nAll genes (", method.name, ")", sep=""), all.genes$logFC, -log10(all.genes$Adj.P.Value),
         xlab='log2(Fold Change)', ylab='-log10(Adj. P-value)', pch=20)
    abline(h =  1, col = 'grey')
    abline(v = -1, col = 'grey')
    abline(v =  1, col = 'grey')
    with(subset(all.genes, Adj.P.Value <= 0.05 & logFC >= 1),   points(logFC, -log10(Adj.P.Value), pch=20, col="red"))
    with(subset(all.genes, Adj.P.Value <= 0.05 & logFC <=-1),   points(logFC, -log10(Adj.P.Value), pch=20, col="green"))
    with(subset(all.genes, Adj.P.Value <= 0.05 & abs(logFC)>1), textxy(logFC, -log10(Adj.P.Value), labs=SYMBOL, cex=.8))
    
}

compare_dif_exp_genes = function(GSE, GPL, genes_limma, genes_RP, genes_SAMR, tissue = NULL) {
    # genes_limma <- all.genes.limma
    # genes_RP    <- all.genes.RP
    # genes_SAMR  <- all.genes.SAMR
    up_genes_limma   <- genes_limma$SYMBOL[genes_limma$logFC >=  1 & genes_limma$Adj.P.Value <= 0.05]
    down_genes_limma <- genes_limma$SYMBOL[genes_limma$logFC <= -1 & genes_limma$Adj.P.Value <= 0.05]
    up_genes_RP      <- genes_RP$SYMBOL[genes_RP$logFC >=  1 & genes_RP$Adj.P.Value <= 0.05]
    down_genes_RP    <- genes_RP$SYMBOL[genes_RP$logFC <= -1 & genes_RP$Adj.P.Value <= 0.05]
    up_genes_SAMR    <- genes_SAMR$SYMBOL[genes_SAMR$logFC >=  1 & genes_SAMR$Adj.P.Value <= 0.05]
    down_genes_SAMR  <- genes_SAMR$SYMBOL[genes_SAMR$logFC <= -1 & genes_SAMR$Adj.P.Value <= 0.05]
    
    print("colnames(genes_limma)")
    print(colnames(genes_limma))
    print(dim(genes_limma))
    print("colnames(genes_RP)")
    print(colnames(genes_RP))
    print(dim(genes_RP))
    print("colnames(genes_SAMR)")
    print(colnames(genes_SAMR))
    print(dim(genes_SAMR))
    
    genes_limma_simple <- NULL
    genes_RP_simple    <- NULL
    genes_SAMR_simple  <- NULL
    
    # Saving dif. exp. results
    if (!is.null(genes_limma) && !('Method' %in% colnames(genes_limma))) {
        genes_limma_simple <- genes_limma[, c('PROBEID', 'ENTREZID', 'SYMBOL', 'logFC', 'Adj.P.Value')]
        genes_limma_simple <- data.frame(genes_limma_simple, 'LIMMA')
        #colnames(genes_limma_simple)[5] <- 'Adj.P.value'
        colnames(genes_limma_simple)[6] <- 'Method'
    } else {
        genes_limma_simple <- genes_limma
    }
    
    if (!is.null(genes_RP) && !('Method' %in% colnames(genes_RP))) {
        genes_RP_simple    <- genes_RP[, c('PROBEID', 'ENTREZID', 'SYMBOL', 'logFC', 'Adj.P.Value')]
        genes_RP_simple <- data.frame(genes_RP_simple, 'RankProd')
        #colnames(genes_RP_simple)[5] <- 'Adj.P.value'
        colnames(genes_RP_simple)[6] <- 'Method'
    } else {
        genes_RP_simple <- genes_RP
    }
    
    if (!is.null(genes_SAMR) && !('Method' %in% colnames(genes_SAMR))) {
        genes_SAMR_simple <- genes_SAMR[, c('PROBEID', 'ENTREZID', 'SYMBOL', 'logFC', 'Adj.P.Value')]
        genes_SAMR_simple <- data.frame(genes_SAMR_simple, 'SAMR')
        #colnames(genes_SAMR_simple)[5] <- 'Adj.P.value'
        colnames(genes_SAMR_simple)[6] <- 'Method'
    } else {
        genes_SAMR_simple <- genes_SAMR
    }
    
    # Plot top 3 genes in each category
    all.genes <- data.frame(genes_limma_simple[, c('PROBEID', 'ENTREZID', 'SYMBOL', 'logFC', 'Adj.P.Value')],
                            genes_RP_simple[, c('Adj.P.Value')],
                            genes_SAMR_simple[, c('Adj.P.Value')])
    colnames(all.genes)[5:ncol(all.genes)] <- c('adjP.LIMMA', 'adjP.RP', 'adjP.SAMR')
    #rownames(all.genes) <- all.genes$PROBEID # Not possible because of duplicated probe Ids
    
    top.up.genes   <- head(all.genes[order(all.genes$logFC, decreasing = TRUE), ], 3)
    top.down.genes <- head(all.genes[order(all.genes$logFC), ], 3)
    top.genes <- rbind(top.up.genes, top.down.genes)
    
    for (i in 1:nrow(top.genes)) {
        probeId <- top.genes$PROBEID[i]
        symbol  <- top.genes$SYMBOL[i]
        #par(mfrow=c(1,2), oma=c(0,0,4,0))
        #boxplot(exprsByGene, ylim=c(0, max(exprsByGene, na.rm=TRUE)), col=pheno_colors)
        #cat(probeId, '\t', nrow(exprsByGene), '\t', ncol(exprsByGene), '\t')
        #cat(exprsByGene[1,], '\n')
        boxplot(exprsByGene[which(rownames(exprsByGene)==probeId), ] ~ phenotypes,
                main=paste('\n', symbol, '|', probeId,
                            '| logFC =', round(top.genes$logFC[i], 2),
                            '\nAdj.P (LIMMA) =', top.genes$adjP.LIMMA[i],
                            '\nAdj.P (RankProd) =', top.genes$adjP.RP[i],
                            '\nAdj.P (SAM) =', top.genes$adjP.SAMR[i])) #, outer=TRUE)
        #par(mfrow=c(1,1), oma=c(0,0,0,0))
    }
    
    genes_simple <- rbind(genes_limma_simple, genes_RP_simple, genes_SAMR_simple)
    
    up_genes <- unique(c(up_genes_limma, up_genes_RP, up_genes_SAMR))
    
    if (length(up_genes) > 0) {
        Counts <- matrix(0, nrow=length(up_genes), ncol=3)
        for (i in 1:nrow(Counts)) {
            Counts[i,1] <- up_genes[i] %in% up_genes_limma
            Counts[i,2] <- up_genes[i] %in% up_genes_RP
            Counts[i,3] <- up_genes[i] %in% up_genes_SAMR
        }
        
        area1 <- length(up_genes_limma)
        area2 <- length(up_genes_RP)
        area3 <- length(up_genes_SAMR)
        n12 <- sum(up_genes %in% up_genes_limma & up_genes %in% up_genes_RP)
        n13 <- sum(up_genes %in% up_genes_limma & up_genes %in% up_genes_SAMR)
        n23 <- sum(up_genes %in% up_genes_RP & up_genes %in% up_genes_SAMR)
        n123 <- sum(up_genes %in% up_genes_limma & up_genes %in% up_genes_RP & up_genes %in% up_genes_SAMR)
        
        grid.newpage()
        draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, n12 = n12, n13 = n13, n23 = n23, n123 = n123,
                         category = c(paste("Up LIMMA", tissue, sep="\n"), paste("Up RankProd", tissue, sep="\n"), paste("Up SAMR", tissue, sep="\n")), cat.pos = c(0,0,180))
    }
    
    down_genes <- unique(c(down_genes_limma, down_genes_RP, down_genes_SAMR))
        
    if (length(down_genes) > 0) {
        Counts <- matrix(0, nrow=length(down_genes), ncol=3)
        for (i in 1:nrow(Counts)) {
            Counts[i,1] <- down_genes[i] %in% down_genes_limma
            Counts[i,2] <- down_genes[i] %in% down_genes_RP
            Counts[i,3] <- down_genes[i] %in% down_genes_SAMR
        }
        
        area1 <- length(down_genes_limma)
        area2 <- length(down_genes_RP)
        area3 <- length(down_genes_SAMR)
        n12 <- sum(down_genes %in% down_genes_limma & down_genes %in% down_genes_RP)
        n13 <- sum(down_genes %in% down_genes_limma & down_genes %in% down_genes_SAMR)
        n23 <- sum(down_genes %in% down_genes_RP & down_genes %in% down_genes_SAMR)
        n123 <- sum(down_genes %in% down_genes_limma & down_genes %in% down_genes_RP & down_genes %in% down_genes_SAMR)
        
        grid.newpage()
        draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, n12 = n12, n13 = n13, n23 = n23, n123 = n123,
                         category = c(paste("Down LIMMA", tissue, sep="\n"), paste("Down RankProd", tissue, sep="\n"), paste("Down SAMR", tissue, sep="\n")), cat.pos = c(0,0,180))
    }
    
    #return(genes_simple)
    return(all.genes)
}
