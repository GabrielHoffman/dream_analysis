#! /usr/bin/env Rscript 

suppressPackageStartupMessages(library(getopt))

# opt = list(n_samples = 20, n_reps = 3, n_de_genes = 500, disease_fc = 2, hsq = .4)

a = matrix(c(
	'fasta', 		'a', 1, "character",
	'n_samples', 	's', 1, "integer",
	'n_reps', 		'r', 1, "integer",
	'n_de_genes', 	'd', 1, "integer",
	'disease_fc', 	'f', 1, "integer",
	'nthreads', 	't', 1, "integer",
	'seed', 		'i', 1, "integer",
	'prefix', 		'p', 1, "character",
	'param_ID', 	'k', 1, "character", 
	'param_Disease', 'e', 1, "character",
	'param_Batch', 	'b', 1, "character",
	'out', 			'o', 1, "character"
),nrow=4)

opt = getopt( t(a));

# Process variance fraction parameters
#######################################

# parameters are mean and variance of the variance fraction
opt$array_ID = as.numeric(strsplit(opt$param_ID, ' ')[[1]])
opt$array_Disease = as.numeric(strsplit(opt$param_Disease, ' ')[[1]])
opt$array_Batch = as.numeric(strsplit(opt$param_Batch, ' ')[[1]])

# convert mean and variances in parameters of beta distribution
# https://en.wikipedia.org/wiki/Beta_distribution#Parameter_estimation
# there is a limit on the variance given constraints of the beta distribution.
# If var exceeds this limit, assign maximum allows variance
estimate_beta_parameters = function(mu, v){
	v = min(v, 0.8*mu*(1-mu))

	alpha = mu*(mu*(1-mu)/v - 1)
	beta = (1-mu) *(mu*(1-mu)/v - 1)

	c(alpha=alpha, beta=beta)
}

opt$distr_ID = estimate_beta_parameters( opt$array_ID[1], opt$array_ID[2] )
opt$distr_Disease = estimate_beta_parameters( opt$array_Disease[1], opt$array_Disease[2] )
opt$distr_Batch = estimate_beta_parameters( opt$array_Batch[1], opt$array_Batch[2] )

# Load packages
###############
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(polyester)) 
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(mvtnorm))

# set seed
set.seed( opt$seed )

# read transcripts
fastaTranscripts = readDNAStringSet( opt$fasta )

# Create design matrix
######################

n_samples = opt$n_samples
n_reps = opt$n_reps
n_de_genes = opt$n_de_genes
disease_fc = opt$disease_fc
h_sq = opt$hsq

info = data.frame( Individual = paste("ID", sort(rep(1:n_samples, n_reps)), sep=''),
	Disease = as.character(sort(rep(0:1, n_samples*n_reps / 2))), 
	Experiment = paste0("sample_", gsub(" ", "0", format(1:(n_samples*n_reps), width=3)), sep=''))

info$Batch = factor(sample(0:1, nrow(info), replace=TRUE))

# sampling until design matrix is not singular
idx = seq(1, nrow(info), by=table(info$Individual)[1])

while( min(svd(model.matrix(~Disease + Batch, info))$d) <=0 || min(svd(model.matrix(~Disease + Batch, info[idx,]))$d) <=0 || (nlevels(info$Batch) < 2) ){
	info$Batch = factor(sample(0:1, nrow(info), replace=TRUE))
}

# draw names of differentially expressed genes
##############################################

deGeneIdx = sample(1:length(fastaTranscripts), opt$n_de_genes, replace=FALSE)

# Simulate fold changes from variance components
################################################

design_ID = model.matrix( ~ 0+Individual,info)
design_Disease = model.matrix( ~ 0+Disease,info)
design_Batch = model.matrix( ~ 0+Batch,info)

simParams = foreach(j=1:length(fastaTranscripts) ) %do% {
	cat("\r", j, "        ")

	# Individual
	eta_ID = design_ID %*% rnorm(nlevels(info$Individual))

	# Batch
	eta_batch = design_Batch %*% rnorm(nlevels(info$Batch))

	# draw variance fractions
	sigSq_ID = rbeta(1, opt$distr_ID[1], opt$distr_ID[2])
	sigSq_Batch = rbeta(1, opt$distr_Batch[1], opt$distr_Batch[2])

	# each individual has its own error variance
	v_var = rbeta( nlevels(info$Individual), 1, 1) + 0.5
	resid_var = model.matrix( ~ 0+Individual, info) %*% v_var

	if( j %in% deGeneIdx){
		# Disease
		eta_Disease = design_Disease %*% rnorm(nlevels(info$Disease))

		sigSq_Disease = rbeta(1, opt$distr_Disease[1], opt$distr_Disease[2])
		sigSq_Resid = max(1 - sigSq_ID - sigSq_Disease - sigSq_Batch, .05)

		# combine
		y = scale(eta_ID) * (sigSq_ID-sigSq_Disease) + 
			scale(eta_batch) * sigSq_Batch + 
			scale(eta_Disease) * sigSq_Disease +			 
			rnorm(nrow(info), 0, sd=sqrt(resid_var)) * sigSq_Resid
	}else{
		sigSq_Resid = max(1 - sigSq_ID - sigSq_Batch, .05)

		# combine
		y = scale(eta_ID) * sigSq_ID + 
			scale(eta_batch) * sigSq_Batch + 		 
			rnorm(nrow(info), sd=sqrt(resid_var)) * sigSq_Resid
	}

	# fit = lm( y ~ Disease, info)
	# fit = lmer( y ~ (1|Disease) + (1|Batch) + (1|Individual), info)
	# calcVarPart(fit)

	list( FC = t(y) - min(y) + 1)
}

# store fold changes for each gene and sample
FC = lapply( simParams, function(x) x$FC)
FC = do.call("rbind", FC)
rownames(info) = info$Experiment
colnames(FC) = info$Experiment
rownames(FC) = sapply(strsplit(names(fastaTranscripts), '\\|'), function(x) x[2])

# Just for testing
# August 6, 2019
# library(BiocParallel)
# register(SnowParam(12, "SOCK", progressbar=TRUE))

# info$Batch = factor(info$Batch)
# vp = fitExtractVarPartModel( FC[1:2000,], ~ (1|Individual) + (1|Disease) + (1|Batch), info)

# fig = plotVarPart(vp)
# ggsave("Rplots.png", fig)
# q()


lib_sizes = runif(nrow(info), .5, 1.5) /4

# Simulate read counts
######################

# overwrite this function in polyester to decrease runtime
assignInNamespace('sgseq', function(x,...){1}, "polyester")

# try this
#reads_per_transcript=readspertx[1:length(readspertx)]
# meanmodel=FALSE,

FC_scale = t(apply(FC, 1, function(x){
	x = x / 2
	x - min(x) + 1
	})) 

# reads_per_transcript = 2^runif(length(fastaTranscripts), 2, 14)

# b0 = -3.0158
# b1 = 0.8688
# sigma = 4.152
# logmus = b0 + b1 * log2(width(fastaTranscripts)) + rnorm(length(fastaTranscripts),
#     0, sigma)
# reads_per_transcript = 2^logmus - 1
# reads_per_transcript = pmax(reads_per_transcript, 1)

# size = reads_per_transcript * FC_scale / 3

# reads_per_transcript=reads_per_transcript, size=size,
# Saves count matrix to file, so read it in afterwards
# polyester::simulate_experiment(opt$fasta, transcripts=fastaTranscripts, 
	# meanmodel=TRUE,
 #    num_reps = as.matrix(rep(1, n_samples*n_reps)), fold_changes=FC_scale, lib_sizes=lib_sizes,outdir=paste0(opt$out,'/', opt$prefix), gzip=TRUE, reportCoverage=TRUE, simReads=FALSE)


polyester::simulate_experiment(opt$fasta, transcripts=fastaTranscripts, 
	 meanmodel=TRUE, 
    num_reps = as.matrix(rep(1, n_samples*n_reps)), fold_changes=FC_scale, lib_sizes=lib_sizes,outdir=paste0(opt$out,'/', opt$prefix), gzip=TRUE, reportCoverage=TRUE, simReads=FALSE)

# hist(log2(reads_per_transcript))

# load counts_matrix from file
load(paste0(opt$out,'/', opt$prefix ,'/sim_counts_matrix.rda'))

countMatrix = round(counts_matrix)
colnames(countMatrix) = colnames(FC_scale)
rownames(countMatrix) = sapply(strsplit(rownames(countMatrix), '\\|'), function(x) x[2])


range(colSums(countMatrix))
mean( colSums(countMatrix))

# isexpr = rowSums(cpm(countMatrix)>.1) >= 3
# table(isexpr)

# isexpr[] = TRUE

# # voom single replicate
# genes = DGEList( countMatrix[isexpr,] )
# genes = calcNormFactors( genes )
# design = model.matrix( ~ Disease + Batch, info)

# vobj = voom( genes, design, plot=TRUE)


# vobj = voomWithDreamWeights( genes, ~ Disease + (1|Batch) + (1|Individual), info, plot=TRUE)


# i = which.max(rowSums(cpm(countMatrix)))
# cpm(countMatrix)[i,]
# cpm(countMatrix)[i+1,]


# dds <- DESeqDataSetFromMatrix(countData = countMatrix,
#                               colData = info,
#                               design = ~ Disease + Batch )
# dds = DESeq(dds)

# plotDispEsts(dds)

# Save results
###############

saveRDS(info, paste(opt$out, "/info_", opt$prefix, ".RDS", sep=''))


deGeneList = sapply(strsplit(names(fastaTranscripts)[deGeneIdx], '\\|'), function(x) x[2])
saveRDS(deGeneList, paste(opt$out, "/deGeneList_", opt$prefix, ".RDS", sep=''))

saveRDS(countMatrix, paste(opt$out, "/countMatrix_", opt$prefix, ".RDS", sep=''))
save(list=ls(), file=paste(opt$out, "/infoAll_", opt$prefix, ".RData", sep=''))





# Check variance components
###########################

# info$Batch = factor(info$Batch)

# # # filter out genes based on read count
# isexpr = rowSums(cpm(countMatrix)>.1) >= 3
# countMatrix = countMatrix[isexpr,]
# rownames(countMatrix) = sapply(strsplit(rownames(countMatrix), '\\|'), function(x) x[2])

# timeMethods = list()

# # voom single replicate
# genes = DGEList( countMatrix )
# genes = calcNormFactors( genes )
# design = model.matrix( ~ Disease + Batch, info)

# vobj = voom( genes, design, plot=TRUE)

# form <- ~ (1|Disease) + (1|Batch) + (1|Individual) 
# # form <- ~ Batch
# vp_cpm = fitExtractVarPartModel(vobj[1:2000,], form, info)

# fig = plotVarPart(vp_cpm)
# fig
# ggsave("Rplots.png", fig)




# df = merge(vp, vp_cpm, by="row.names")


# plot(df$Disease.x, df$Disease.y)
# abline(0,1, col="red")

# plot(df$Individual.x, df$Individual.y)
# abline(0,1, col="red")

# plot(df$Batch.x, df$Batch.y)
# abline(0,1, col="red")


# plot(df$Residuals.x, df$Residuals.y)
# abline(0,1, col="red")



