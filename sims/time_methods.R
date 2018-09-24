#! /usr/bin/env Rscript

# Gabriel Hoffman
# Evaluate methods on simulated data and save the time

suppressPackageStartupMessages(library(getopt))

a = matrix(c(
	'folder', 		'f', 1, "character",
	'prefix', 		'p', 1, "character",
	'out', 			'o', 1, "character",
	'runkr', 		'k', 0, "logical",
	'nthreads', 	't', 1, "integer"
	),nrow=4)

opt = getopt( t(a));

suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(doParallel))

# data(varPartData)
# form <- ~ Batch + (1|Individual)
# idx = (info$Tissue == 'A')
# L = getContrast( geneExpr[,idx], form, info[idx,], "Batch3")

# fit = dream( geneExpr[,idx], form, info[idx,], L)
# fitEB = eBayes( fit )

# fit = dream( geneExpr[,idx], form, info[idx,], L, ddf="Kenward-Roger")
# fitEB = eBayes( fit )


# opt = list(folder="/hpc/users/hoffmg01/work/RNA_seq_sim_v1/", prefix='6_3_500_3_0.4_44')

countMatrix = readRDS(paste0(opt$folder, '/data/countMatrix_',opt$prefix,'.RDS'))
info = readRDS(paste0(opt$folder, '/data/info_',opt$prefix,'.RDS'))
rownames(info) = info$Experiment

isexpr = rowSums(cpm(countMatrix)>0.1) >= 3
# isexpr = rowSums(cpm(countMatrix)>0.5) >= 300

# voom single replicate
idx = seq(1, nrow(info), by=table(info$Individual)[1])
genes = DGEList( countMatrix[,idx] )
genes = calcNormFactors( genes )
design = model.matrix( ~ Disease, info[idx,])
vobj = voom( genes, design, plot=FALSE)
design = model.matrix( ~ Disease, info[idx,])
fit_lmFit = lmFit(vobj, design)
fit_lmFit = eBayes(fit_lmFit)

# Create VOBJ
genes = DGEList( countMatrix )
genes = calcNormFactors( genes )
design = model.matrix( ~ Disease, info)
vobj_tmp = voom( genes, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=info$Individual)
vobj = voom( genes, design, plot=FALSE,block=info$Individual,correlation=dupcor$consensus)

# dupCor
cat("dupCor...\n")
start_dupcor = Sys.time()
design = model.matrix( ~ Disease, info)
dupcor <- duplicateCorrelation(vobj,design,block=info$Individual)
fitDupCor <- lmFit(vobj,design,block=info$Individual,correlation=dupcor$consensus)
fitDupCor <- eBayes(fitDupCor)
end_dupcor = Sys.time()


cl = makeCluster(opt$nthreads)
registerDoParallel( cl )

if( ! is.null(opt$runkr) ){
	cat("dreamKR...\n")
	# dream: Kenward-Roger approximation
	start_dreamkr = Sys.time()
	form <- ~ Disease + (1|Individual) 
	L = getContrast( vobj, form, info, "Disease1")
	fit2KR = dream( vobj, form, info, L, ddf='Kenward-Roger')
	fit2eKR = eBayes( fit2KR )
	end_dreamkr = Sys.time()
}

# dream: Satterthwaite approximation
cat("dream...\n")
start_dream = Sys.time()
form <- ~ Disease + (1|Individual) 
L = getContrast( vobj, form, info, "Disease1")
fitSat = dream( vobj, form, info, L, ddf='Satterthwaite')
fitSatEB = eBayes( fitSat )
end_dream = Sys.time()



if( ! is.null(opt$runkr) ){
timeResults = data.frame( 	dupCor = difftime(end_dupcor, start_dupcor, units="min"),
							dream = difftime(end_dream, start_dream, units="min"),
							dreamkr = difftime(end_dreamkr, start_dreamkr, units="min"))
}else{	
timeResults = data.frame( 	dupCor = difftime(end_dupcor, start_dupcor, units="min"),
							dream = difftime(end_dream, start_dream, units="min"))
}

stopCluster(cl)

file = paste0( opt$folder, '/', opt$out, '.tsv')
write.table(timeResults, file=file, row.names=FALSE, sep="\t", quote=FALSE)








