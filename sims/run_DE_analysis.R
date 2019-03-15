#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(getopt))

a = matrix(c(
	'prefix', 	'p', 1, "character",
	'folder', 	'o', 1, "character",
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
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(lmms))
suppressPackageStartupMessages(library(MACAU2))

# read data from simulations
countMatrix = readRDS(paste0(opt$folder, '/data/countMatrix_',opt$prefix,'.RDS'))
info = readRDS(paste0(opt$folder, '/data/info_',opt$prefix,'.RDS'))
rownames(info) = info$Experiment

# retain all genes for this analysis
isexpr = rowSums(cpm(countMatrix)>0.1) >= 0

# voom single replicate
idx = seq(1, nrow(info), by=table(info$Individual)[1])
genes = DGEList( countMatrix[,idx] )
genes = calcNormFactors( genes )
design = model.matrix( ~ Disease, info[idx,])
vobj = voom( genes, design, plot=FALSE)
design = model.matrix( ~ Disease, info[idx,])
fit_lmFit = lmFit(vobj, design)
fit_lmFit = eBayes(fit_lmFit)


# DESeq2 one sample
dds_single <- DESeqDataSetFromMatrix(countData = countMatrix[,idx],
                              colData = info[idx,],
                              design= ~ Disease )
dds_single = DESeq(dds_single)


# Create VOBJ
genes = DGEList( countMatrix )
genes = calcNormFactors( genes )
design = model.matrix( ~ Disease, info)
vobj_tmp = voom( genes, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=info$Individual)
vobj = voom( genes, design, plot=FALSE,block=info$Individual,correlation=dupcor$consensus)

# include both replicates, don't account
fit_lmFit2 = lmFit(vobj, design)
fit_lmFit2 = eBayes(fit_lmFit2)

# DESeq2 all samples
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = info,
                              design= ~ Disease )
dds <- DESeq(dds)

# head(results(dds))
# table(results(dds)$padj < 0.05)

# Sum reads from replicates
###########################

# DESeq2 
# Sum reads by sample
countMatrix_sum = lapply( unique(info$Individual), function(ID){
	rowSums(countMatrix[,info$Experiment[info$Individual == ID],drop=FALSE])
	} )
countMatrix_sum2 = do.call("cbind", countMatrix_sum)
colnames(countMatrix_sum2) = info$Experiment[idx]


dds_sum <- DESeqDataSetFromMatrix(countData = countMatrix_sum2,
                              colData = info[idx,],
                              design= ~ Disease )
dds_sum = DESeq(dds_sum)

# limma
genes_sum = DGEList( countMatrix_sum2 )
genes_sum = calcNormFactors( genes_sum )
design_sum = model.matrix( ~ Disease, info[idx,])
vobj_tmp_sum = voom( genes_sum, design_sum, plot=FALSE)
fit_lmFit_sum = lmFit(vobj_tmp_sum, design_sum)
fit_lmFit_sum = eBayes(fit_lmFit_sum)



# dupCor
design = model.matrix( ~ Disease, info)
dupcor <- duplicateCorrelation(vobj,design,block=info$Individual)
fitDupCor <- lmFit(vobj,design,block=info$Individual,correlation=dupcor$consensus)
fitDupCor <- eBayes(fitDupCor)

cl = makeCluster(5)
registerDoParallel( cl )

# dream: Kenward-Roger approximation
form <- ~ Disease + (1|Individual) 
L = getContrast( vobj, form, info, "Disease1")
fit2KR = dream( vobj, form, info, L, ddf='Kenward-Roger')
fit2eKR = eBayes( fit2KR )

# variancePartition
form <- ~ (1|Disease) + (1|Individual) 
vp = fitExtractVarPartModel(vobj, form, info)

# dream: Satterthwaite approximation
form <- ~ Disease + (1|Individual) 
L = getContrast( vobj, form, info, "Disease1")
fitSat = dream( vobj, form, info, L, ddf='Satterthwaite')
fitSatEB = eBayes( fitSat )

# MACAU2
########

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

macau_fit = macau2(countMatrix, info$Disease, RelatednessMatrix=K, fit.model="PMM",numCore=1)#), fit.maxiter=20)

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
de_res = data.frame( true = rep(0,nrow(countMatrix)))
de_res$true[1:500] = 1
de_res$lmFit = topTable(fit_lmFit, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val
de_res$DESeq2_single = results(dds_single)$padj
de_res$lmFit_sum = topTable(fit_lmFit_sum, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val
de_res$DESeq2_sum = results(dds_sum)$padj
de_res$lmFit2 = topTable(fit_lmFit2, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val
de_res$DESeq2 = results(dds)$padj
de_res$macau2 = p.adjust(macau_fit$pvalue, 'fdr')
de_res$lmFit_dupCor = topTable(fitDupCor, coef='Disease1', sort.by="none", number=Inf)$adj.P.Val
de_res$lmm_Sat = p.adjust(fitSat$pValue, "fdr") 
de_res$lmm_Sat_eBayes = topTable(fitSatEB, sort.by="none", number=Inf)$adj.P.Val
de_res$lmm_KR = p.adjust(fit2KR$pValue, "fdr") 
de_res$lmm_KR_eBayes = topTable(fit2eKR, sort.by="none", number=Inf)$adj.P.Val


# Write results
file = paste0(opt$folder, '/results/', opt$prefix, "_p.adj.RDS")
saveRDS(de_res, file)

# pValue
de_res_p = data.frame( true = rep(0,nrow(countMatrix)))
de_res_p$true[1:500] = 1
de_res_p$lmFit = topTable(fit_lmFit, coef='Disease1', sort.by="none", number=Inf)$P.Value
de_res_p$DESeq2_single = results(dds_single)$pvalue
de_res_p$lmFit_sum = topTable(fit_lmFit_sum, coef='Disease1', sort.by="none", number=Inf)$P.Value
de_res_p$DESeq2_sum = results(dds_sum)$pvalue
de_res_p$lmFit2 = topTable(fit_lmFit2, coef='Disease1', sort.by="none", number=Inf)$P.Value
de_res_p$DESeq2 = results(dds)$pvalue
de_res$macau2 = macau_fit$pvalue
de_res_p$lmFit_dupCor = topTable(fitDupCor, coef='Disease1', sort.by="none", number=Inf)$P.Value
de_res_p$lmm_Sat = fitSat$pValue
de_res_p$lmm_Sat_eBayes = topTable(fitSatEB, sort.by="none", number=Inf)$P.Value
de_res_p$lmm_KR = fit2KR$pValue
de_res_p$lmm_KR_eBayes = topTable(fit2eKR, sort.by="none", number=Inf)$P.Value


file = paste0(opt$folder, '/results/', opt$prefix, "_p.RDS")
saveRDS(de_res_p, file)

file = paste0(opt$folder, '/results/', opt$prefix, "_vp.RDS")
saveRDS(vp, file)


# Compare p-values and make plot
t1 = topTable(fitDupCor, coef="Disease1", number=Inf, sort.by="none")$P.Value
t2 = topTable(fitSatEB, number=Inf, sort.by="none")$P.Value
# t2 = fitSat@pValue

# df2 = data.frame(idx=1:length(t1), dupCor=-log10(t1), fitMixedModelTest=-log10(t2), vp, delta = vp$Individual < dupcor$consensus)

# l1 = lm(fitMixedModelTest ~ dupCor, df2[df2$delta,])
# l2 = lm(fitMixedModelTest ~ dupCor, df2[!df2$delta,])

# df_line = data.frame(rbind(coef(l1), coef(l2)))
# colnames(df_line) = c('a', 'b')
# df_line$type = c('blue', 'red')

# fig2 = ggplot(df2, aes(dupCor, fitMixedModelTest, color = Individual, shape=idx <=500)) + geom_abline() + geom_point(size=2) + theme_bw(12) + scale_colour_gradient2(low="blue", mid="green", high="red", midpoint=dupcor$consensus) + theme(aspect.ratio=1) +
# 	xlab("dupCor (-log10 p)") + ylab("fitMixedModelTest (-log10 p)") +
# 	geom_abline( intercept=df_line$a, slope=df_line$b, color=df_line$type, linetype=2)


fig2 = plotCompareP( t1, t2, vp$Individual, dupcor$consensus)

file = paste0(opt$folder, '/figures/', opt$prefix, "_prelim.pdf")
pdf( file )
fig2
plotVarPart(vp, main="All")
plotVarPart(vp[1:500,], main="DE genes")
dev.off()





