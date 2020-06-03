#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(getopt))

a = matrix(c(
	'folder', 		'f', 1, "character",
	'noEB', 		'n', 0, "logical",
	'nthreads', 	't', 1, "integer"
	),nrow=4)

opt = getopt( t(a));

folder = opt$folder

# Combine results across simulations
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(binom))
suppressPackageStartupMessages(library(data.table))

allFiles = dir(folder, '_(\\d+)_p.RDS', full.names=FALSE)
prefixes = sort(unique(gsub("_(\\d+)_p.RDS", "", allFiles)))

cat("# files: ", length(allFiles), "\n")
cat("# prefix: ", length(prefixes), "\n")

# registerDoParallel( opt$nthreads )

# define colors
################
# library(RColorBrewer)
# brewer.pal(n = 12, name = "Paired")

colarray = c("Single replicate (limma/voom)", "#A6CEE3",
"Single replicate (DESeq2)", "#1F78B4",
"Sum reads (limma/voom)", "#B2DF8A",
"Sum reads (DESeq2)", "#33A02C",
"Full data, ignore corr (limma/voom)",  "#FB9A99",
"Full data, ignore corr (DESeq2)", "#E31A1C",
"macau2", "gold",
"duplicateCorrelation with limma/voom", "#B15928",
"dream + FMT.vc",  "orange",
"dream + FMT.ws",  "orange3",
"dream (KR)", "#FDBF6F",
"dream",  "#FF7F00",
"dream (KR) [std p]", "#CAB2D6",
"dream [std p]",  "#6A3D9A")
df_plots = data.frame(matrix( colarray, ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
df_plots$X1 = factor(df_plots$X1, df_plots$X1)
colnames(df_plots) = c("method", "color")

if( opt$noEB ){
	df_plots$color[df_plots$method == 'dream'] = '#6A3D9A'
	df_plots$color[df_plots$method == 'dream (KR)'] = '#CAB2D6'
}

## Confidence intervals for AUPR
## https://github.com/kboyd/raucpr/blob/master/precision_recall.r
# aucpr.conf.int.expit <- function(estimate, num.pos, num.neg,conf.level=0.95) {
#   ## Calculates confidence interval for an AUCPR estimate using expit.
  
#   ## convert to logit scale
#   est.logit = log(estimate/(1-estimate))
#   ## standard error (from Kevin Eng)
#   se.logit = sqrt(estimate*(1-estimate)/num.pos)*(1/estimate + 1/(1-estimate))
#   ## confidence interval in logit
#   ci.logit = est.logit+qnorm(c((1-conf.level)/2,(1+conf.level)/2))*se.logit

#   ## back to original scale
#   ci = exp(ci.logit)/(1+exp(ci.logit))
#   attr(ci,"conf.level") = conf.level
#   attr(ci,"method") = "expit"
#   return(ci)
# }

#################
# Plot run time #
#################

resTime = foreach( prefix = prefixes, .combine=rbind ) %dopar% {

	n_donor = as.numeric(gsub("^(\\d+)_(\\d+).*", "\\1", prefix))
	n_reps = as.numeric(gsub("^(\\d+)_(\\d+).*", "\\2", prefix))

	# read results
	files = dir(folder, paste0('^',prefix, '.*_timeMethods.RDS'), full.names=TRUE)

	df_time = foreach( file = files, .combine=rbind) %do% {
		res = data.frame(do.call("rbind", readRDS( file )))

		data.frame(name=rownames(res), time=res$elapsed, 
			n_donor=n_donor, n_reps=n_reps, stringsAsFactors=FALSE)
	}

	df_time$name[df_time$name == 'fit_lmFit']  = "Single replicate (limma/voom)"
	df_time$name[df_time$name == 'dds_single'] = "Single replicate (DESeq2)"
	df_time$name[df_time$name == 'limma_sum'] = "Sum reads (limma/voom)"
	df_time$name[df_time$name == 'DESeq2_sum'] = "Sum reads (DESeq2)"
	df_time$name[df_time$name == 'fit_lmFit2'] = "Full data, ignore corr (limma/voom)"
	df_time$name[df_time$name == 'DESeq2']  = "Full data, ignore corr (DESeq2)"
	df_time$name[df_time$name == 'lmFit_dupCor']  = "duplicateCorrelation with limma/voom"
	df_time$name[df_time$name == 'lmm_Sat'] = "dream"
	df_time$name[df_time$name == 'lmm_KR'] = "dream (KR)"
	df_time$name[df_time$name == 'FMT.vc'] = "dream + FMT.vc"
	df_time$name[df_time$name == 'FMT.ws'] = "dream + FMT.ws"
	df_time$name[df_time$name == 'macau'] = "macau2"
	colnames(df_time)[1] = "method"

	df_time
}

resTime = data.table(resTime)

summTime = resTime[,data.frame(mean=mean(time), sd=sd(time)),by=c("method", "n_donor", "n_reps")]

summTime$method = factor(summTime$method, df_plots$method)
col = df_plots$color[df_plots$method %in% levels(summTime$method)]

file = paste0(folder,'/../figures/combine_time.pdf')
pdf( file )
ggplot(summTime, aes(n_donor,mean/60, color=method)) + geom_line() + geom_errorbar(aes(ymin = (mean-sd)/60, ymax = (mean+sd)/60), width=1) + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, legend.position="bottom")  + scale_color_manual(values=col) + facet_wrap(~n_reps) + xlab("# Donors") + ylab("Run time (minutes)")
dev.off()


################
# Plot results #
################

resList = foreach( prefix = prefixes ) %dopar% {

	message(match(prefix, prefixes))

	n_donor = as.numeric(gsub("^(\\d+)_(\\d+).*", "\\1", prefix))
	n_reps = as.numeric(gsub("^(\\d+)_(\\d+).*", "\\2", prefix))

	# read results
	files = dir(folder, paste0('^',prefix, '.*_p.RDS'), full.names=TRUE)

	de_res = foreach( file = files, .combine=rbind) %do% {
		readRDS( file )
	}
	# for DESeq2, set NA to 1
	de_res[is.na(de_res)] = 1

	if( opt$noEB ){
		de_res[,colnames(de_res) == 'lmm_Sat_eBayes'] = c()
		de_res[,colnames(de_res) == 'lmm_KR_eBayes'] = c()
	}

	colnames(de_res)[colnames(de_res) == 'lmFit' ] = "Single replicate (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2_single' ] = "Single replicate (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit_sum' ] = "Sum reads (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2_sum' ] = "Sum reads (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit2' ] = "Full data, ignore corr (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2' ] = "Full data, ignore corr (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit_dupCor' ] = "duplicateCorrelation with limma/voom"
	colnames(de_res)[colnames(de_res) == 'lmm_Sat_eBayes' ] = "dream"
	colnames(de_res)[colnames(de_res) == 'lmm_KR_eBayes' ] = "dream (KR)"
	colnames(de_res)[colnames(de_res) == 'lmm_Sat' ] = "dream" #dream [std p]"
	colnames(de_res)[colnames(de_res) == 'lmm_KR' ] = "dream (KR)" #"dream (KR) [std p]"
	colnames(de_res)[colnames(de_res) == 'FMT.vc' ] = "dream + FMT.vc"
	colnames(de_res)[colnames(de_res) == 'FMT.ws' ] = "dream + FMT.ws"

	# raw p-values on right
	keepRaw = TRUE
	if( keepRaw ){
		de_res_right = de_res[,colnames(de_res) %in% c("lmm_Sat", "lmm_KR")]
		de_res = de_res[,!(colnames(de_res) %in% c("lmm_Sat", "lmm_KR"))]
		de_res = cbind(de_res, de_res_right)

		# specify colors
		# method_order = colnames(de_res)[-c(1, 12,13)]
		# col = ggColorHue(length(method_order))
		# names(col) = method_order
		# col['lmm_Sat'] = alpha(col['dream'], .5)
		# col['lmm_KR'] = alpha(col['dream (KR)'], .5)
		# method_order = c(method_order, c("lmm_Sat", "lmm_KR"))
	}else{
		de_res = de_res[,!(colnames(de_res) %in% c("lmm_Sat", "lmm_KR"))]
		
		method_order = colnames(de_res)[-1]
		col = ggColorHue(length(method_order))
		names(col) = method_order
	}

	# counts at each cutoff
	#======================
	df = foreach( key = colnames(de_res)[-1], .combine=rbind )  %do% {

		# pv = de_res[[key]]

		# p_cutoff = sort(unique(unlist(pv)))
		# v = max(c(sort(p_cutoff)[1000], 1e-3))
		# p_cutoff = p_cutoff[p_cutoff < v]

		# n_chosen = sapply(p_cutoff, function(x) sum(pv <= x))
		# n_false = sapply(p_cutoff, function(x) sum((pv <= x) & !de_res[['true']]) )

		# much faster version using sorting of p-values
		pv = de_res[[key]]
		true_state = !de_res[['true']]

		idx = order(pv)
		pv = pv[idx]
		true_state = true_state[idx]

		p_cutoff = sort(unique(unlist(pv)))
		v = max(c(sort(p_cutoff)[10000], 1e-2))
		p_cutoff = p_cutoff[p_cutoff < v]

		n_chosen = findInterval(p_cutoff,pv)

		pv_false = pv[true_state]
		n_false = findInterval(p_cutoff,pv_false)

		data.frame(n_donor, n_reps, key, p_cutoff, n_chosen, n_false, stringsAsFactors=FALSE)
	}
	df$key = droplevels(factor(df$key, df_plots$method))
	df$n_chosen = df$n_chosen / length(files)
	df$n_false = df$n_false / length(files)

	col = df_plots$color[match(levels(df$key), as.character(df_plots$method))]

	# plot false positive versus positives
	fig_choose = ggplot(df, aes(n_chosen, n_false, color=key)) + geom_line() + xlim(0, 375) + ylim(0, 50) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=col) + xlab("# genes selected") + ylab("# false positives")

	# Compute PR
	#===========
	prList = foreach( method = colnames(de_res)[-1] ) %do% {
		if( length(unique(de_res[[method]])) == 1 ){
			value = NA
		}else{
			value = pr.curve(scores.class0 = -log10(de_res[de_res$true==1,method]), 
						scores.class1 = -log10(de_res[de_res$true==0,method]), 
						curve=TRUE, rand.compute = TRUE)

			# compute AUPR confidence interval
			# n_pos = sum(de_res$true==1)
			# n_neg = sum(de_res$true==0)
			# value$ci = aucpr.conf.int.expit( value$auc.integral, n_pos, n_neg)

			# compute confidence interval based on asymptotic approximation
			n = length(de_res$true)
			value$ci = binom.confint( value$auc.integral * n, n, methods="asymptotic")[5:6]
		}
		value
	}
	names(prList) = colnames(de_res)[-1] 

	# remove empty entry
	prList[sapply(prList, function(x) all(is.na(x)))] = c()

	aupr = lapply( names(prList), function(method){
		 data.frame(method 	= method,
		 			n_donor = n_donor,
		 			n_reps 	= n_reps,
		 			value 	= prList[[method]]$auc.integral,
		 			low 	= prList[[method]]$ci[1],
		 			high 	= prList[[method]]$ci[2] )
		})
	aupr = do.call(rbind, aupr)

	# aupr = data.frame(method = names(prList) )
	# aupr$value = foreach( method = names(prList), .combine=c ) %do% {
	# 	prList[[method]]$auc.integral
	# }
	aupr$method = droplevels(factor( as.character(aupr$method), df_plots$method))
	# aupr$n_donor = n_donor
	# aupr$n_reps = n_reps
	aupr.rand.score = prList[[method]]$rand$auc.integral	
	col = df_plots$color[match(levels(aupr$method), as.character(df_plots$method))]

	fig_aupr = ggplot(aupr, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_errorbar(aes(method, ymin=lower, ymax=upper), width=.2) + coord_flip() + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("AUPR") + geom_hline(yintercept=aupr.rand.score, linetype=2) 

	# Plot PR
	#========
	
	dfpr = foreach( method = names(prList), .combine=rbind ) %do% {
		pr = data.table(as.data.frame(prList[[method]]$curve))
		colnames(pr) = c( "recall", "precision", "score")
		pr = pr[,score:=NULL]
		pr[,precision:=round(precision,3)]
		res = data.frame(unique(pr), method)
		res[order(res$precision, decreasing=TRUE),]
	}
	dfpr$method = droplevels(factor(dfpr$method, df_plots$method))
	col = df_plots$color[match(levels(dfpr$method), df_plots$method)]

	rnd.value = prList[[1]]$rand$curve[1,2]
	dfpr$rnd.value = rnd.value 
	dfpr$n_donor = n_donor
	dfpr$n_reps = n_reps

	fig_pr = ggplot(dfpr, aes(recall, precision, color=method)) + geom_line() + theme_bw(12) + xlab("Recall") + ylab("Precision") + xlim(0,1) + ylim(0,1) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=rnd.value, color="grey50", linetype=2) + geom_hline(yintercept=aupr.rand.score, linetype=2) + scale_color_manual(values=col)
	
	# Power	at 5% FDR
	#==================
	power_fdr_5 = data.frame(method = names(prList))
	power_fdr_5$value = foreach( method = names(prList), .combine=c ) %do% {
		pr = as.data.frame(prList[[method]]$curve)
		colnames(pr) = c( "recall", "precision", "score")
		# ggplot(pr, aes(recall, precision)) + geom_line()

		i = which.min(abs(pr$precision - 0.95))
		pr$recall[i]
	}
	power_fdr_5$method = droplevels(factor(power_fdr_5$method, df_plots$method))
	power_fdr_5$n_donor = n_donor
	power_fdr_5$n_reps = n_reps
	col = df_plots$color[match(levels(power_fdr_5$method), df_plots$method)]

	fig_power_fdr_5 = ggplot(power_fdr_5, aes(method, value, fill=method)) + geom_bar(stat="identity") + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("Power at FDR 5%") + coord_flip()

	# False positive rate
	#======================
	fpr = apply(de_res[de_res$true ==0,], 2, function(x) sum(x< 0.05)/length(x))
	df_fpr = data.frame(method=names(fpr)[-1])
	df_fpr$value = fpr[-1]
	df_fpr$method = droplevels(factor(df_fpr$method, df_plots$method))
	df_fpr$n_donor = n_donor
	df_fpr$n_reps = n_reps

	# confidence interval for false positive rate
	n = sum(de_res$true==0) # total number of tests
	df_fpr = cbind(df_fpr, binom.confint( df_fpr$value * n, n, method = 'asymptotic' )[,5:6])

	col = df_plots$color[match(levels(df_fpr$method), df_plots$method)]

	fig_fpr = ggplot(df_fpr, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + geom_errorbar(aes(method, ymin=lower, ymax=upper), width=.2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()

	# FP after FDR
	# #============
	# files = dir(folder, paste0('^',prefix, '.*_p.adj.RDS'), full.names=TRUE)
	# de_res = foreach( file = files, .combine=rbind) %do% {
	# 	readRDS( file )
	# }
	# # for DESeq2, set NA to 1
	# de_res[is.na(de_res)] = 1

	# if( opt$noEB ){
	# 	de_res[,colnames(de_res) == 'lmm_Sat_eBayes'] = c()
	# 	de_res[,colnames(de_res) == 'lmm_KR_eBayes'] = c()
	# }

	# colnames(de_res)[colnames(de_res) == 'lmFit' ] = "Single replicate (limma/voom)"
	# colnames(de_res)[colnames(de_res) == 'DESeq2_single' ] = "Single replicate (DESeq2)"
	# colnames(de_res)[colnames(de_res) == 'lmFit_sum' ] = "Sum reads (limma/voom)"
	# colnames(de_res)[colnames(de_res) == 'DESeq2_sum' ] = "Sum reads (DESeq2)"
	# colnames(de_res)[colnames(de_res) == 'lmFit2' ] = "Full data, ignore corr (limma/voom)"
	# colnames(de_res)[colnames(de_res) == 'DESeq2' ] = "Full data, ignore corr (DESeq2)"
	# colnames(de_res)[colnames(de_res) == 'lmFit_dupCor' ] = "duplicateCorrelation with limma/voom"
	# # colnames(de_res)[colnames(de_res) == 'lmm_Sat_eBayes' ] = "dream"
	# # colnames(de_res)[colnames(de_res) == 'lmm_KR_eBayes' ] = "dream (KR)"
	# colnames(de_res)[colnames(de_res) == 'lmm_Sat' ] = "dream"
	# colnames(de_res)[colnames(de_res) == 'lmm_KR' ] = "dream (KR)"

	# de_res = de_res[,!(colnames(de_res) %in% c("lmm_Sat", "lmm_KR", "lmm_KR_eBayes" ))]

	# fd = apply(de_res[de_res$true ==0,], 2, function(x) sum(x< 0.05))
	fd = apply(de_res[de_res$true ==0,], 2, function(x) sum(p.adjust(x, "BH") < 0.05))

	false_discoveries = data.frame(method=names(fd)[-1])
	# divide false discoveries by number of analyses
	false_discoveries$value = as.numeric(fd[-1] / length(files))
	false_discoveries$method = droplevels(factor(false_discoveries$method, df_plots$method))
	false_discoveries$n_donor = n_donor
	false_discoveries$n_reps = n_reps
	col = df_plots$color[match(levels(false_discoveries$method), df_plots$method )]

	fig_fd = ggplot(false_discoveries, aes(method, value, fill=method)) + geom_bar(stat="identity") + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("False discoveries at FDR <0.05") + coord_flip()

	list( df = df, 
		fig_choose = fig_choose, 
		aupr = aupr, 
		fig_aupr = fig_aupr, 
		fig_pr = fig_pr, 
		power_fdr_5 = power_fdr_5,
		fig_power_fdr_5 = fig_power_fdr_5, 
		df_fpr = df_fpr, 
		fig_fpr = fig_fpr, 
		fig_fd = fig_fd, 
		dfpr = dfpr,
		false_discoveries = false_discoveries,
		prefix = prefix)
}
names(resList) = prefixes

donor_array = sort(unique(sapply(strsplit( prefixes, "_"), function(x) as.numeric(x[1]))))


# # Plot combing across datasets
# #=============================
# file = paste0(folder,'/../figures/combine_across_datasets.pdf')
# pdf( file, width=18, height=16 )
# AUPR
aupr = foreach(n_donor = donor_array, .combine=rbind) %do% {
	foreach(n_reps = 2:4, .combine=rbind) %do% {
		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
			resList[[idx]]$aupr 
		}
	}	
}
aupr$method = droplevels(factor(aupr$method, df_plots$method))

# Drop dream-EB
# if( opt$noEB ){
# 	aupr = aupr[!(aupr$method %in% c("dream (KR)", "dream")),]
# 	aupr = droplevels(aupr)	
# 	levels(aupr$method)[levels(aupr$method)=='dream [std p]'] <- "dream"
# 	levels(aupr$method)[levels(aupr$method)=='dream (KR) [std p]'] <- "dream (KR)"
# }

col = df_plots$color[match(levels(aupr$method), as.character(df_plots$method))]

fig = foreach(n_reps = 2:4) %do%{
ggplot(aupr[aupr$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_color_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, 1) + ylab("AUPR")
}
fig_aupr = do.call("grid.arrange", c(fig, ncol=3))

# power_fdr_5
power_fdr_5 = foreach(n_donor = donor_array, .combine=rbind) %do% {
	foreach(n_reps = 2:4, .combine=rbind) %do% {
		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
			resList[[idx]]$power_fdr_5
		}
	}	
}
power_fdr_5$method = droplevels(factor(power_fdr_5$method, df_plots$method))

# Drop dream-EB
# if( opt$noEB ){
# 	power_fdr_5 = power_fdr_5[!(power_fdr_5$method %in% c("dream (KR)", "dream")),]
# 	power_fdr_5 = droplevels(power_fdr_5)	
# 	levels(power_fdr_5$method)[levels(power_fdr_5$method)=='dream [std p]'] <- "dream"
# 	levels(power_fdr_5$method)[levels(power_fdr_5$method)=='dream (KR) [std p]'] <- "dream (KR)"
# }

fig = foreach(n_reps = 2:4) %do%{
ggplot(power_fdr_5[power_fdr_5$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_color_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, 1) + ylab("Power at FDR 5%")
}
fig_power_fdr_5 = do.call("grid.arrange", c(fig, ncol=3))

df_fpr = foreach(n_donor = donor_array, .combine=rbind) %do% {
	foreach(n_reps = 2:4, .combine=rbind) %do% {
		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
			resList[[idx]]$df_fpr	
		}
	}	
}
df_fpr$method = droplevels(factor(df_fpr$method, df_plots$method))

# Drop dream-EB
# if( opt$noEB ){
# 	df_fpr = df_fpr[!(df_fpr$method %in% c("dream (KR)", "dream")),]
# 	df_fpr = droplevels(df_fpr)
# 	levels(df_fpr$method)[levels(df_fpr$method)=='dream [std p]'] <- "dream"
# 	levels(df_fpr$method)[levels(df_fpr$method)=='dream (KR) [std p]'] <- "dream (KR)"
# }

df_fpr = df_fpr[!((df_fpr$method == "macau2")&&(df_fpr$n_donor >14)&&(df_fpr$value == 0)),]

fig = foreach(n_reps = 2:4) %do%{
ggplot(df_fpr[df_fpr$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_hline(yintercept=0.05, linetype=2) + geom_line() + scale_color_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, .2) + ylab("False positive rate") 
}
fig_fpr = do.call("grid.arrange", c(fig, ncol=3))

# False discoveries
df_fd = foreach(n_donor = donor_array, .combine=rbind) %do% {
	foreach(n_reps = 2:4, .combine=rbind) %do% {
		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
			resList[[idx]]$false_discoveries	
		}
	}	
}
df_fd$method = droplevels(factor(df_fd$method, df_plots$method))

# if( opt$noEB ){
# 	df_fd = df_fd[!(df_fd$method %in% c("dream (KR)", "dream")),]
# 	df_fd = droplevels(df_fd)
# 	levels(df_fd$method)[levels(df_fd$method)=='dream [std p]'] <- "dream"
# 	levels(df_fd$method)[levels(df_fd$method)=='dream (KR) [std p]'] <- "dream (KR)"
# }
df_fd = df_fd[!((df_fd$method == "macau2")&&(df_fd$n_donor >14)&&(df_fd$value == 0)),]

fig = foreach(n_reps = 2:4) %do%{
ggplot(df_fd[df_fd$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_color_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, max(df_fd$value)) + ylab("False discoveries at FDR < 5%") 
}
fig_fd = do.call("grid.arrange", c(fig, ncol=3))

fig = foreach(n_reps = 2:4) %do%{
	mvalue = max(df_fd[df_fd$n_reps==n_reps,]$value)
ggplot(df_fd[df_fd$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_color_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, mvalue) + ylab("False discoveries at FDR < 5%") 
}
fig_fd2 = do.call("grid.arrange", c(fig, ncol=3))
# dev.off()

file = paste0(folder,'/../figures/combine_across_datasets.pdf')
pdf( file, width=8, height=12 )
grid.arrange(fig_aupr, fig_power_fdr_5, fig_fpr, fig_fd, fig_fd2,ncol=1)
dev.off()

# Plot PR curves separately due to size
#######################################

# gc()

# file = paste0(folder,'/../figures/','combine_PR', ".pdf")
# pdf( file, width=7, height=25 )
# # PR
# fig_pr = foreach(n_donor = donor_array, .combine=c) %do% {
# 	foreach(n_reps = 2:4, .combine=c) %do% {
# 		cat(n_donor, " ", n_reps, "\n")
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			resList[[idx]]$fig_pr + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + ylim(0, 1) + theme(axis.text.y =  element_text(size=9), axis.text.x =  element_text(size=9))
# 		}
# 	}	
# }
# do.call("grid.arrange", c(fig_pr, ncol=3))
# dev.off()
# dev.off()

# False Discoveries
###################

df = do.call("rbind", lapply(resList, function(x) x$df))
df = data.table(df)

# if( opt$noEB ){
# 	df = df[!(df$key %in% c("dream (KR)", "dream")),]
# 	df = droplevels(df)
# 	levels(df$key)[levels(df$key)=='dream [std p]'] <- "dream"
# 	levels(df$key)[levels(df$key)=='dream (KR) [std p]'] <- "dream (KR)"
# }

# saveRDS(df, file=paste0(opt$folder,'/df.RDS'))
# df = readRDS(paste0(opt$folder,'/df.RDS'))

# col = ggColorHue(length(table(df$key)))
col = df_plots$color[match( levels(df$key), df_plots$method)]

file = paste0(folder,'/../figures/','combine_choose', ".pdf")
pdf( file, width=7, height=20)
fig1 = ggplot(df[(n_donor <= 14),], aes(n_chosen, n_false, color=key)) + geom_line() + facet_grid( n_donor ~ n_reps) + xlim(0, 300) + ylim(0, 50) + theme_bw(10) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size=8), legend.position="bottom") + scale_color_manual(values=col) + xlab("genes selected") + ylab("false discoveries")
if( nrow(df[n_donor > 14,]) > 0 ){
	fig2 = ggplot(df[(n_donor > 14),], aes(n_chosen, n_false, color=key)) + geom_line() + facet_grid( n_donor ~ n_reps) + xlim(0, 550) + ylim(0, 50) + theme_bw(10) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), axis.text.x =  element_text(size=8), legend.position="bottom") + scale_color_manual(values=col) + xlab("genes selected") + ylab("false discoveries")
	grid.arrange(fig1, fig2, ncol=2)
}else{
	fig1
}
dev.off()

# Precision recall
##################

dfpr = do.call("rbind", lapply(resList, function(x) x$dfpr))
dfpr = data.table(dfpr)

# if( opt$noEB ){
# 	dfpr = dfpr[!(dfpr$method %in% c("dream (KR)", "dream")),]
# 	dfpr = droplevels(dfpr)
# 	levels(dfpr$method)[levels(dfpr$method)=='dream [std p]'] <- "dream"
# 	levels(dfpr$method)[levels(dfpr$method)=='dream (KR) [std p]'] <- "dream (KR)"
# }

# saveRDS(dfpr, file=paste0(opt$folder,'/dfpr.RDS'))
# dfpr = readRDS(paste0(opt$folder,'/dfpr.RDS'))

# col = ggColorHue(length(table(df$key)))
col = df_plots$color[match(levels(dfpr$method), df_plots$method)]
randCurve = dfpr[,unique(rnd.value)]

file = paste0(folder,'/../figures/','combine_pr2', ".pdf")
pdf( file, width=8, height=15)
fig1 = ggplot(dfpr[(n_donor <= 14),], aes(recall, precision, color=method)) + geom_line()  + theme_bw(10) + facet_grid( n_donor ~ n_reps) + xlim(0,1) + ylim(0,1) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="bottom", axis.text.x=element_text(size=6)) + scale_color_manual(values=col) + geom_hline(yintercept=randCurve, linetype=2)
if( nrow(df[n_donor > 14,]) > 0 ){
	fig2 = ggplot(dfpr[(n_donor > 14),], aes(recall, precision, color=method)) + geom_line()  + theme_bw(10) + facet_grid( n_donor ~ n_reps) + xlim(0,1) + ylim(0,1) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="bottom", axis.text.x=element_text(size=6)) + scale_color_manual(values=col) + geom_hline(yintercept=randCurve, linetype=2)
	grid.arrange(fig1, fig2, ncol=2)
}else{
	fig1
}
dev.off()

# False positive rate
######################

df_fpr = do.call("rbind", lapply(resList, function(x) x$df_fpr))
df_fpr = data.table(df_fpr)

# if( opt$noEB ){
# 	df_fpr = df_fpr[!(df_fpr$method %in% c("dream (KR)", "dream")),]
# 	df_fpr = droplevels(df_fpr)
# 	levels(df_fpr$method)[levels(df_fpr$method)=='dream [std p]'] <- "dream"
# 	levels(df_fpr$method)[levels(df_fpr$method)=='dream (KR) [std p]'] <- "dream (KR)"
# }

# saveRDS(df_fpr, file=paste0(opt$folder,'/df_fpr.RDS'))
# df_fpr = readRDS(paste0(opt$folder,'/df_fpr.RDS'))


file = paste0(folder,'/../figures/','combine_fpr', ".pdf")
pdf( file, width=15, height=20)
maxValue = max(df_fpr$upper)

df_a = df_fpr[(n_donor <= 14),]
df_a$method = droplevels(df_a$method)
col = df_plots$color[match(levels(df_a$method), df_plots$method)]

fig1 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + geom_errorbar(aes(method, ymin=lower, ymax=upper), width=.2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)

df_a = df_fpr[(n_donor > 14),]
if( nrow(df_a) > 0 ){
# drop macau2 if all empty
# if( df_a[method=="macau2",unique(value)] == 0 ){
# 	df_a = df_a[method!="macau2",]
# }
	df_a$method = droplevels(df_a$method)
	col = df_plots$color[match(levels(df_a$method), df_plots$method)]

	fig2 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2)  + geom_errorbar(aes(method, ymin=lower, ymax=upper), width=.2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)
	grid.arrange(fig1, fig2, ncol=2)
}else{
	fig1
}
dev.off()


# AUPR
######

df_aupr = do.call("rbind", lapply(resList, function(x) x$aupr))
df_aupr = data.table(df_aupr)

# if( opt$noEB ){
# 	df_aupr = df_aupr[!(df_aupr$method %in% c("dream (KR)", "dream")),]
# 	df_aupr = droplevels(df_aupr)
# 	levels(df_aupr$method)[levels(df_aupr$method)=='dream [std p]'] <- "dream"
# 	levels(df_aupr$method)[levels(df_aupr$method)=='dream (KR) [std p]'] <- "dream (KR)"
# }

# saveRDS(df_aupr, file=paste0(opt$folder,'/df_aupr.RDS'))
# df_aupr = readRDS(paste0(opt$folder,'/df_aupr.RDS'))

# col = ggColorHue(length(table(df$key)))

file = paste0(folder,'/../figures/','combine_aupr', ".pdf")
pdf( file, width=15, height=20)
maxValue = max(df_aupr[(n_donor <= 14),upper])

df_a = df_aupr[(n_donor <= 14),]
df_a$method = droplevels(df_a$method)
col = df_plots$color[match(levels(df_a$method), df_plots$method)]
fig1 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + geom_errorbar(aes(method, ymin=lower, ymax=upper), width=.2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("AUPR") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)

df_a = df_aupr[(n_donor > 14),]
if( nrow(df_a) > 0 ){
	df_a$method = droplevels(df_a$method)
	col = df_plots$color[match(levels(df_a$method), df_plots$method)]
	fig2 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + geom_errorbar(aes(method, ymin=lower, ymax=upper), width=.2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("AUPR") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, 1)
	grid.arrange(fig1, fig2, ncol=2)
}else{
	fig1
}
dev.off()








