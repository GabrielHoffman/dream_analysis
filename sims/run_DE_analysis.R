#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(getopt))

a = matrix(c(
	'prefix', 	'p', 1, "character",
	'folder', 	'o', 1, "character",
	'macau2', 	'm', 0, "logical",
	'kr', 		'k', 0, "logical"
),nrow=4)

opt = getopt( t(a));

# source("helper_functions.R")

suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(lmms))
suppressPackageStartupMessages(library(MACAU2))


# install MACAU2
# module load udunits proj gdal geos
# in R
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# remotes::install_github("jakyzhu/MACAU2")


# read data from simulations
countMatrixOrig = readRDS(paste0(opt$folder, '/data/countMatrix_',opt$prefix,'.RDS'))

info = readRDS(paste0(opt$folder, '/data/info_',opt$prefix,'.RDS'))
rownames(info) = info$Experiment

# filter out genes based on read count
isexpr = rowSums(cpm(countMatrixOrig)>.1) >= 3
# isexpr[] = TRUE
countMatrix = countMatrixOrig[isexpr,]

# read DE list
file = paste0("data/deGeneList_", opt$prefix, ".RDS")
deGenes = readRDS(file)

timeMethods = list()

# voom single replicate
idx = seq(1, nrow(info), by=table(info$Individual)[1])
genes = DGEList( countMatrix[,idx] )
genes = calcNormFactors( genes )
design = model.matrix( ~ Disease + Batch, info[idx,])

timeMethods$fit_lmFit = system.time({
file = paste0(opt$folder, '/figures/voom_', opt$prefix, ".pdf")
pdf( file )
vobj = voom( genes, design, plot=TRUE)
dev.off()
design = model.matrix( ~ Disease + Batch, info[idx,])
fit_lmFit = lmFit(vobj, design)
fit_lmFit = eBayes(fit_lmFit)
})

# DESeq2 one sample
# DESeq2 get full count matrix, no filtering
timeMethods$dds_single = system.time({
dds_single <- DESeqDataSetFromMatrix(countData = countMatrixOrig[,idx],
                              colData = info[idx,],
                              design = ~ Disease + Batch )
dds_single = DESeq(dds_single)
})

# Create VOBJ
timeMethods$fit_lmFit2 = system.time({
genes = DGEList( countMatrix )
genes = calcNormFactors( genes )
design = model.matrix( ~ Disease + Batch, info)
vobj_tmp = voom( genes, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=info$Individual)
vobj = voom( genes, design, plot=FALSE,block=info$Individual,correlation=dupcor$consensus)

# include both replicates, don't account
fit_lmFit2 = lmFit(vobj, design)
fit_lmFit2 = eBayes(fit_lmFit2)
})

# DESeq2 all samples
timeMethods$DESeq2 = system.time({
dds <- DESeqDataSetFromMatrix(countData = countMatrixOrig,
                              colData = info,
                              design= ~ Disease + Batch)
dds <- DESeq(dds)
})

# head(results(dds))
# table(results(dds)$padj < 0.05)



# Sum reads from replicates
###########################

# DESeq2 
# Sum reads by sample
timeMethods$DESeq2_sum = system.time({
countMatrix_sum = lapply( unique(info$Individual), function(ID){
	rowSums(countMatrixOrig[,info$Experiment[info$Individual == ID],drop=FALSE])
	} )
countMatrix_sum2 = do.call("cbind", countMatrix_sum)
colnames(countMatrix_sum2) = info$Experiment[idx]


dds_sum <- DESeqDataSetFromMatrix(countData = countMatrix_sum2,
                              colData = info[idx,],
                              design= ~ Disease )
dds_sum = DESeq(dds_sum)
})

# limma
timeMethods$limma_sum = system.time({
genes_sum = DGEList( countMatrix_sum2 )
genes_sum = calcNormFactors( genes_sum )
design_sum = model.matrix( ~ Disease, info[idx,])
vobj_tmp_sum = voom( genes_sum, design_sum, plot=FALSE)
fit_lmFit_sum = lmFit(vobj_tmp_sum, design_sum)
fit_lmFit_sum = eBayes(fit_lmFit_sum)
})


# dupCor
timeMethods$lmFit_dupCor = system.time({
design = model.matrix( ~ Disease + Batch, info)
dupcor <- duplicateCorrelation(vobj,design,block=info$Individual)
fitDupCor <- lmFit(vobj,design,block=info$Individual,correlation=dupcor$consensus)
fitDupCor <- eBayes(fitDupCor)
})


BPPARAM = SnowParam(6, "SOCK", progressbar=TRUE)
register(BPPARAM)

genes = DGEList( countMatrix )
genes = calcNormFactors( genes )
form <- ~ Disease + (1|Batch) + (1|Individual) 

file = paste0(opt$folder, '/figures/voomWithDreamWeights_', opt$prefix, ".pdf")
pdf( file )
vobjDream = voomWithDreamWeights( genes, form, info, plot=TRUE)
dev.off()

# dream: Kenward-Roger approximation
timeMethods$lmm_KR = system.time({
form <- ~ Disease + (1|Batch) + (1|Individual) 
fit2KR = dream( vobjDream, form, info, ddf='Kenward-Roger')
fit2eKR = eBayes( fit2KR )
})

# variancePartition
form <- ~ (1|Disease) + (1|Batch) + (1|Individual) 
vp = fitExtractVarPartModel(vobjDream, form, info)

# dream: Satterthwaite approximation
timeMethods$lmm_Sat = system.time({
form <- ~ Disease + (1|Batch) + (1|Individual) 
fitSat = dream( vobjDream, form, info)
fitSatEB = eBayes( fitSat )
})

# plot(-log10(fitSat$p.value), -log10(fitSatEB$p.value))
# abline(0,1, col="red")


# fit2KR = dream( vobjDream[2313:2329,], form, info, ddf='Kenward-Roger')


# library(lmerTest)
# i=2313
# y = t(t(vobjDream$E[i,]))
# w = t(t(vobjDream$weights[i,]))
# fit = lmer( y ~ Disease + (1|Batch) + (1|Individual), info, weights= w / max(w), REML=FALSE)
# summary(fit)

# fit = lmer( y ~ Disease + (1|Batch) + (1|Individual), info, weights= w / max(w), REML=TRUE )
# summary(fit, ddf="Ken" )

# where does df come from in dream


# res = variancePartition:::.standardized_t_stat(limma::eBayes(fitSat))
# plot(-log10(fitSat$p.value), -log10(res$p.value))
# abline(0,1, col="red")

# res = limma::eBayes(variancePartition:::.standardized_t_stat(fitSat))
# plot(-log10(fitSat$p.value), -log10(res$p.value))
# abline(0,1, col="red")



# MACAU2
########

if( ! is.null(opt$macau2) && opt$macau2 ){
	timeMethods$macau = system.time({
		# create block diagonal relatedness matrix
		K = matrix(0, nrow(info), nrow(info))
		diag(K) = 1
		rownames(K) = info$Experiment
		colnames(K) = info$Experiment

		for( ID in unique(info$Individual) ){
			expr = info$Experiment[info$Individual==ID]
			i = which(rownames(K) %in% expr)
			K[i,i] = 1
		}

		# K[1:5, 1:5]
		macau_fit <- macau2(countMatrix, info$Disease, data.frame(info$Batch), RelatednessMatrix=K, fit.model="PMM",numCore=6, filtering=FALSE)
	}) 

	macau_fit$p.adj = p.adjust(macau_fit$pvalue, 'fdr')

}else{
	macau_padj = macau_pvalue = rep(NA, nrow(countMatrixOrig))
	timeMethods$macau = NA
}

# ggplot(de_res_p, aes(-log10(lmFit), -log10(macau2))) + geom_point()

# lmms
######

# fit_lmmsDE = lmmsDE( data = t(countMatrix[1:10,]),
# 						time = rep(0, nrow(info)),
# 						sampleID = info$Individual,
# 						group = info$Disease, knots=0 )




###################
# combine results #
###################

# Adjusted pvalue
#----------------
de_res = data.frame( EnsID = rownames(countMatrixOrig), true = rep(0,nrow(countMatrixOrig)), stringsAsFactors=FALSE)
de_res$true[de_res$EnsID %in% deGenes] = 1

# fit_lmFit
df = data.frame(EnsID = rownames(fit_lmFit), 
	lmFit = topTable(fit_lmFit, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# DESeq2_single
df = data.frame(EnsID = rownames(dds_single), 
	DESeq2_single = results(dds_single, contrast =c('Disease', '1', '0'))$padj, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# lmFit_sum
df = data.frame(EnsID = rownames(fit_lmFit_sum), 
	lmFit_sum = topTable(fit_lmFit_sum, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# DESeq2_sum
df = data.frame(EnsID = rownames(dds_sum), 
	DESeq2_sum = results(dds_sum, contrast =c('Disease', '1', '0'))$padj, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# lmFit2
df = data.frame(EnsID = rownames(fit_lmFit2), 
	lmFit2 = topTable(fit_lmFit2, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# DESeq2
df = data.frame(EnsID = rownames(dds), 
	DESeq2 = results(dds, contrast =c('Disease', '1', '0'))$padj, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# macau2
df = data.frame(EnsID = rownames(macau_fit), 
	macau2 = macau_fit$p.adj, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# lmFit_dupCor
df = data.frame(EnsID = rownames(fitDupCor), 
	lmFit_dupCor = topTable(fitDupCor, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# lmm_Sat
df = data.frame(EnsID = rownames(fitSat), 
	lmm_Sat = topTable(fitSat, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )


# lmm_Sat_eBayes
df = data.frame(EnsID = rownames(fitSatEB), 
	lmm_Sat_eBayes = topTable(fitSatEB, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# lmm_KR
df = data.frame(EnsID = rownames(fit2KR), 
	lmm_KR = topTable(fit2KR, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

# lmm_KR_eBayes
df = data.frame(EnsID = rownames(fit2eKR), 
	lmm_KR_eBayes = topTable(fit2eKR, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val, stringsAsFactors=FALSE)
de_res = merge( de_res, df, by="EnsID", all=TRUE )

rownames(de_res) = de_res$EnsID
de_res = de_res[,-1]


# pValue
#-------

de_res_p = data.frame( EnsID = rownames(countMatrixOrig), true = rep(0,nrow(countMatrixOrig)), stringsAsFactors=FALSE)
de_res_p$true[de_res_p$EnsID %in% deGenes] = 1

# fit_lmFit
df = data.frame(EnsID = rownames(fit_lmFit), 
	lmFit = topTable(fit_lmFit, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# DESeq2_single
df = data.frame(EnsID = rownames(dds_single), 
	DESeq2_single = results(dds_single, contrast =c('Disease', '1', '0'))$pvalue, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# lmFit_sum
df = data.frame(EnsID = rownames(fit_lmFit_sum), 
	lmFit_sum = topTable(fit_lmFit_sum, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# DESeq2_sum
df = data.frame(EnsID = rownames(dds_sum), 
	DESeq2_sum = results(dds_sum, contrast =c('Disease', '1', '0'))$pvalue, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# lmFit2
df = data.frame(EnsID = rownames(fit_lmFit2), 
	lmFit2 = topTable(fit_lmFit2, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# DESeq2
df = data.frame(EnsID = rownames(dds), 
	DESeq2 = results(dds, contrast =c('Disease', '1', '0'))$pvalue, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# macau2
df = data.frame(EnsID = rownames(macau_fit), 
	macau2 = macau_fit$pvalue, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# lmFit_dupCor
df = data.frame(EnsID = rownames(fitDupCor), 
	lmFit_dupCor = topTable(fitDupCor, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# lmm_Sat
df = data.frame(EnsID = rownames(fitSat), 
	lmm_Sat = topTable(fitSat, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )


# lmm_Sat_eBayes
df = data.frame(EnsID = rownames(fitSatEB), 
	lmm_Sat_eBayes = topTable(fitSatEB, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# lmm_KR
df = data.frame(EnsID = rownames(fit2KR), 
	lmm_KR = topTable(fit2KR, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

# lmm_KR_eBayes
df = data.frame(EnsID = rownames(fit2eKR), 
	lmm_KR_eBayes = topTable(fit2eKR, coef='Disease1', sort.by="none", number=Inf)$P.Value, stringsAsFactors=FALSE)
de_res_p = merge( de_res_p, df, by="EnsID", all=TRUE )

rownames(de_res_p) = de_res_p$EnsID
de_res_p = de_res_p[,-1]


# save results to file
#######################

file = paste0(opt$folder, '/results/', opt$prefix, "_p.adj.RDS")
saveRDS(de_res, file)

file = paste0(opt$folder, '/results/', opt$prefix, "_p.RDS")
saveRDS(de_res_p, file)

file = paste0(opt$folder, '/results/', opt$prefix, "_vp.RDS")
saveRDS(vp, file)

file = paste0(opt$folder, '/results/', opt$prefix, "_timeMethods.RDS")
saveRDS(timeMethods, file)


# Compare p-values and make plot
################################

t1 = topTable(fitDupCor, coef="Disease1", number=Inf, sort.by="none")$P.Value
t2 = topTable(fitSatEB, number=Inf, sort.by="none")$P.Value

fig2 = plotCompareP( t1, t2, vp$Individual, dupcor$consensus)

file = paste0(opt$folder, '/figures/', opt$prefix, "_prelim.pdf")
pdf( file )
fig2
plotVarPart(vp, main="All")
plotVarPart(vp[rownames(vp) %in% deGenes,], main="DE genes")
dev.off()





