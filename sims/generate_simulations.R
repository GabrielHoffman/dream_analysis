#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(getopt))

a = matrix(c(
	'fasta', 		'a', 1, "character",
	'n_samples', 	's', 1, "integer",
	'n_reps', 		'r', 1, "integer",
	'n_de_genes', 	'd', 1, "integer",
	'disease_fc', 	'f', 1, "integer",
	'hsq', 			'h', 1, "numeric",
	'nthreads', 	't', 1, "integer",
	'seed', 		'i', 1, "integer",
	'prefix', 		'p', 1, "character",
	'out', 			'o', 1, "character"
),nrow=4)

opt = getopt( t(a));

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(polyester)) # v1.7.1
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))

source("~/scripts/varPartSims/helper_functions.R")

# cl = makeCluster(opt$nthreads)
# registerDoParallel( cl )

set.seed( opt$seed )

# fasta_use = '~/work/RNA_seq_sim/transcriptome/transcripts_all_unique.fa'
fastaTranscripts = readDNAStringSet( opt$fasta )

i = sample(1:length(fastaTranscripts), length(fastaTranscripts), replace=FALSE)
fastaTranscripts = fastaTranscripts[i]

n_reads_per = .09

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(n_reads_per * width(fastaTranscripts) )
# sum(readspertx)

# design matrix
n_samples = opt$n_samples
n_reps = opt$n_reps
n_de_genes = opt$n_de_genes
disease_fc = opt$disease_fc
h_sq = opt$hsq

info = data.frame( Individual = paste("ID", sort(rep(1:n_samples, n_reps)), sep=''),
	Disease = as.character(sort(rep(0:1, n_samples*n_reps / 2))), 
	Experiment = paste("sample_", gsub(" ", "0", format(1:(n_samples*n_reps), digits=2)), sep=''), stringsAsFactors=FALSE)

design = model.matrix( ~ Individual + Disease+0,info)


simParams = foreach(j=1:length(fastaTranscripts), .packages=c("lme4", "variancePartition") ) %do% {
	cat("\r", j, "        ")

	# indiv_variance = rnorm(1, .5, 2)
	indiv_variance = 1
	beta = rnorm(n_samples, 0, max(indiv_variance, .1))

	if( j <= n_de_genes){
		# beta[] = 0
		eta = design %*% c(beta, disease_fc)
		error_var = (1-h_sq)/h_sq  * var(eta)
	}else{		
		eta = design %*% c(beta, 0)

		h_sq_other = rbeta(1, 1, 1.6)
		error_var = (1-h_sq_other)/h_sq_other  * var(eta)
	}	

	y = eta + rnorm(nrow(eta), 0, sqrt(error_var))

	fit = lmer( y ~ (1|Individual) + (1|Disease), info, REML=FALSE)
	v = calcVarPart( fit ) 
	v

	list( FC = t(y) - min(y) + 1, modelStats = v[order(names(v))] )
}

FC = matrix(NA, nrow=length(fastaTranscripts), ncol=n_samples*n_reps)
modelStats = matrix(NA, nrow(FC), ncol=3)
colnames(modelStats) = names(simParams[[1]]$modelStats)

for(j in 1:nrow(FC)){
	FC[j,] = simParams[[j]]$FC
	modelStats[j,] = simParams[[j]]$modelStats
}

deGeneList = names(fastaTranscripts)[1:n_de_genes]
# deGeneList = gsub(".*gene=(\\S+)", "\\1", deGeneList)

lib_sizes = runif(nrow(info), .5, 1.5)

#  rm -f bam/* simulated_reads/*
saveRDS(info, paste(opt$out, "/info_", opt$prefix, ".RDS", sep=''))
saveRDS(deGeneList, paste(opt$out, "/deGeneList_", opt$prefix, ".RDS", sep=''))
saveRDS(modelStats, paste(opt$out, "/modelStats_", opt$prefix, ".RDS", sep=''))
save(list=ls(), file=paste(opt$out, "/infoAll_", opt$prefix, ".RData", sep=''))

# simulation call:

# simulate_experiment(fasta_use, reads_per_transcript=readspertx, 
    # num_reps = rep(1, n_samples*n_reps), fold_changes=FC, lib_sizes=lib_sizes,outdir='~/work/RNA_seq_sim/simulated_reads', gzip=TRUE, reportCoverage=TRUE) 

countMatrix = simulate_experiment_justCounts(transcripts=fastaTranscripts, reads_per_transcript=readspertx, 
    num_reps = rep(1, n_samples*n_reps), fold_changes=FC, lib_sizes=lib_sizes,outdir=opt$out, gzip=TRUE, reportCoverage=TRUE, simReads=FALSE) #, meanmodel=FALSE


rownames(countMatrix) = names(fastaTranscripts)
colnames(countMatrix) = info$Experiment

saveRDS(countMatrix, paste(opt$out, "/countMatrix_", opt$prefix, ".RDS", sep=''))







