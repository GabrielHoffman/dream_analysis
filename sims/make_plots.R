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
suppressPackageStartupMessages(library(data.table))

allFiles = dir(folder, '_(\\d+)_p.RDS', full.names=FALSE)
prefixes = sort(unique(gsub("_(\\d+)_p.RDS", "", allFiles)))

cat("# files: ", length(allFiles), "\n")
cat("# prefix: ", length(prefixes), "\n")

registerDoParallel( opt$nthreads )

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
"dream (KR)", "#FDBF6F",
"dream",  "#FF7F00",
"dream (KR) [std p]", "#CAB2D6",
"dream [std p]",  "#6A3D9A")
df_plots = data.frame(matrix( colarray, ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
df_plots$X1 = factor(df_plots$X1, df_plots$X1)
colnames(df_plots) = c("method", "color")

# library(ggplot2)
# df_plots$value = 1
# ggplot(df_plots, aes(method, value, fill=method)) + geom_bar(stat="identity") + coord_flip() + scale_fill_manual(values=df_plots$color)

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

	n_donor = as.numeric(gsub("^(\\d+)_(\\d+).*", "\\1", prefix))
	n_reps = as.numeric(gsub("^(\\d+)_(\\d+).*", "\\2", prefix))

	# read results
	files = dir(folder, paste0('^',prefix, '.*_p.RDS'), full.names=TRUE)

	de_res = foreach( file = files, .combine=rbind) %do% {
		readRDS( file )
	}
	# for DESeq2, set NA to 1
	de_res[is.na(de_res)] = 1

	colnames(de_res)[colnames(de_res) == 'lmFit' ] = "Single replicate (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2_single' ] = "Single replicate (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit_sum' ] = "Sum reads (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2_sum' ] = "Sum reads (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit2' ] = "Full data, ignore corr (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2' ] = "Full data, ignore corr (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit_dupCor' ] = "duplicateCorrelation with limma/voom"
	colnames(de_res)[colnames(de_res) == 'lmm_Sat_eBayes' ] = "dream"
	colnames(de_res)[colnames(de_res) == 'lmm_KR_eBayes' ] = "dream (KR)"
	colnames(de_res)[colnames(de_res) == 'lmm_Sat' ] = "dream [std p]"
	colnames(de_res)[colnames(de_res) == 'lmm_KR' ] = "dream (KR) [std p]"

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
	df$key = factor(df$key, df_plots$method)
	df$n_chosen = df$n_chosen / length(files)
	df$n_false = df$n_false / length(files)

	col = df_plots$color[df_plots$method %in% levels(df$key)]

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
		}
		value
	}
	names(prList) = colnames(de_res)[-1] 

	# remove empty entry
	prList[sapply(prList, function(x) all(is.na(x)))] = c()

	aupr = data.frame(method = names(prList) )
	aupr$value = foreach( method = names(prList), .combine=c ) %do% {
		prList[[method]]$auc.integral
	}
	aupr$method = factor(aupr$method, df_plots$method)
	aupr$n_donor = n_donor
	aupr$n_reps = n_reps
	aupr.rand.score = prList[[method]]$rand$auc.integral	
	col = df_plots$color[df_plots$method %in% levels(aupr$method)]

	fig_aupr = ggplot(aupr, aes(method, value, fill=method)) + geom_bar(stat="identity") + coord_flip() + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("AUPR") + geom_hline(yintercept=aupr.rand.score, linetype=2) 

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
	dfpr$method = factor(dfpr$method, df_plots$method)
	col = df_plots$color[df_plots$method %in% levels(dfpr$method)]

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
	power_fdr_5$method = factor(power_fdr_5$method, df_plots$method)
	power_fdr_5$n_donor = n_donor
	power_fdr_5$n_reps = n_reps
	col = df_plots$color[df_plots$method %in% levels(power_fdr_5$method)]

	fig_power_fdr_5 = ggplot(power_fdr_5, aes(method, value, fill=method)) + geom_bar(stat="identity") + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("Power at FDR 5%") + coord_flip()

	# False positive rate
	#======================
	fpr = apply(de_res[de_res$true ==0,], 2, function(x) sum(x< 0.05)/length(x))
	df_fpr = data.frame(method=names(fpr)[-1])
	df_fpr$value = fpr[-1]
	df_fpr$method = factor(df_fpr$method, df_plots$method)
	df_fpr$n_donor = n_donor
	df_fpr$n_reps = n_reps
	col = df_plots$color[df_plots$method %in% levels(df_fpr$method)]

	fig_fpr = ggplot(df_fpr, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()

	# FP after FDR
	#============
	files = dir(folder, paste0('^',prefix, '.*_p.adj.RDS'), full.names=TRUE)
	de_res = foreach( file = files, .combine=rbind) %do% {
		readRDS( file )
	}
	# for DESeq2, set NA to 1
	de_res[is.na(de_res)] = 1

	colnames(de_res)[colnames(de_res) == 'lmFit' ] = "Single replicate (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2_single' ] = "Single replicate (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit_sum' ] = "Sum reads (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2_sum' ] = "Sum reads (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit2' ] = "Full data, ignore corr (limma/voom)"
	colnames(de_res)[colnames(de_res) == 'DESeq2' ] = "Full data, ignore corr (DESeq2)"
	colnames(de_res)[colnames(de_res) == 'lmFit_dupCor' ] = "duplicateCorrelation with limma/voom"
	colnames(de_res)[colnames(de_res) == 'lmm_Sat_eBayes' ] = "dream"
	colnames(de_res)[colnames(de_res) == 'lmm_KR_eBayes' ] = "dream (KR)"
	colnames(de_res)[colnames(de_res) == 'lmm_Sat' ] = "dream [std p]"
	colnames(de_res)[colnames(de_res) == 'lmm_KR' ] = "dream (KR) [std p]"

	# de_res = de_res[,!(colnames(de_res) %in% c("lmm_Sat", "lmm_KR", "lmm_KR_eBayes" ))]

	fd = apply(de_res[de_res$true ==0,], 2, function(x) sum(x< 0.05))
	false_discoveries = data.frame(method=names(fd)[-1])
	# divide false discoveries by number of analyses
	false_discoveries$value = as.numeric(fd[-1] / length(files))
	false_discoveries$method = factor(false_discoveries$method, df_plots$method)
	false_discoveries$n_donor = n_donor
	false_discoveries$n_reps = n_reps
	col = df_plots$color[df_plots$method %in% levels(false_discoveries$method)]

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

# file = paste0(folder,'../plot_data_full.RDATA')
# save(list=ls(), file=file)

# Write to PDF
#=============
# file = paste0(folder,'/../figures/','combine', ".pdf")
# pdf( file, width=7, height=25 )
# choose
# fig_choose = foreach(n_donor = donor_array, .combine=c) %do% {
# 	foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			resList[[idx]]$fig_choose + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + theme(plot.title = element_text(hjust = 0.5, size=10))
# 		}
# 	}	
# }
# do.call("grid.arrange", c(fig_choose, ncol=3))

# choose
# fig_choose = foreach(n_donor = donor_array, .combine=c) %do% {
# 	foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			df = data.table(resList[[idx]]$df)
# 			resList[[idx]]$fig_choose + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + xlim(0,  df[n_false < 50,max(n_chosen)]) + theme(plot.title = element_text(hjust = 0.5, size=10))
# 		}
# 	}	
# }
# do.call("grid.arrange", c(fig_choose, ncol=3))

# # AUPR
# fig_aupr = foreach(n_donor = donor_array) %do% {
# 	fl = foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			fig = resList[[idx]]$fig_aupr + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + theme(plot.title = element_text(hjust = 0.5, size=10), axis.text.x = element_text(size=8), aspect.ratio=1)
# 			if( n_reps > 2){
# 				fig = fig + theme(axis.title.y=element_blank(), 
# 		        axis.text.y=element_blank(),
# 		        axis.ticks.y=element_blank())
# 			}
# 			if( n_reps != 3){
# 				fig = fig + ylab('')
# 			}
# 			fig 
# 		}
# 	}
# 	cbind(ggplotGrob(fl[[1]]), ggplotGrob(fl[[2]]), ggplotGrob(fl[[3]]), size='last')	
# }
# grid.arrange(grobs = fig_aupr, ncol=1)


# fig_aupr = foreach(n_donor = donor_array) %do% {
# 	fl = foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			fig = resList[[idx]]$fig_aupr + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + ylim(0, 1) + theme(plot.title = element_text(hjust = 0.5, size=10), axis.text.x = element_text(size=8), aspect.ratio=1)
# 			if( n_reps > 2){
# 				fig = fig + theme(axis.title.y=element_blank(), 
# 		        axis.text.y=element_blank(),
# 		        axis.ticks.y=element_blank())
# 			}
# 			if( n_reps != 3){
# 				fig = fig + ylab('')
# 			}
# 			fig 
# 		}
# 	}
# 	cbind(ggplotGrob(fl[[1]]), ggplotGrob(fl[[2]]), ggplotGrob(fl[[3]]), size='last')	
# }
# grid.arrange(grobs = fig_aupr, ncol=1)

# # power
# fig_power_fdr_5 = foreach(n_donor = donor_array) %do% {
# 	fl = foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			fig = resList[[idx]]$fig_power_fdr_5 + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + ylim(0, 1) + theme(plot.title = element_text(hjust = 0.5, size=10), axis.text.x = element_text(size=8), aspect.ratio=1)
# 			if( n_reps > 2){
# 				fig = fig + theme(axis.title.y=element_blank(), 
# 		        axis.text.y=element_blank(),
# 		        axis.ticks.y=element_blank())
# 			}			
# 			if( n_reps != 3){
# 				fig = fig + ylab('')
# 			}
# 			fig + theme(plot.title = element_text(size=10))
# 		}
# 	}	
# 	cbind(ggplotGrob(fl[[1]]), ggplotGrob(fl[[2]]), ggplotGrob(fl[[3]]), size='last')
# }
# grid.arrange(grobs = fig_power_fdr_5, ncol=1)

# # FPR
# maxfpr = max(sapply(resList, function(x) max(x$df_fpr$value)))

# fig_fpr = foreach(n_donor = donor_array) %do% {
# 	fl = foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			fig = resList[[idx]]$fig_fpr + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + ylim(0, maxfpr)  + theme(plot.title = element_text(hjust = 0.5, size=10), axis.text.x = element_text(size=8), aspect.ratio=1)
# 			if( n_reps > 2){
# 				fig = fig + theme(axis.title.y=element_blank(), 
# 		        axis.text.y=element_blank(),
# 		        axis.ticks.y=element_blank())
# 			}
# 			if( n_reps != 3){
# 				fig = fig + ylab('')
# 			}
# 			fig + theme(plot.title = element_text(size=10))
# 		}
# 	}	
# 	cbind(ggplotGrob(fl[[1]]), ggplotGrob(fl[[2]]), ggplotGrob(fl[[3]]), size='last')	
# }
# grid.arrange(grobs = fig_fpr, ncol=1)

# # False discoveries
# maxfd = max(sapply(resList, function(x) max(x$false_discoveries$value)))

# fig_df = foreach(n_donor = donor_array) %do% {
# 	fl = foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			fig = resList[[idx]]$fig_fd + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps))  + theme(plot.title = element_text(hjust = 0.5, size=10), axis.text.x = element_text(size=8), aspect.ratio=1)
# 			if( n_reps > 2){
# 				fig = fig + theme(axis.title.y=element_blank(), 
# 		        axis.text.y=element_blank(),
# 		        axis.ticks.y=element_blank())
# 			}
# 			if( n_reps != 3){
# 				fig = fig + ylab('')
# 			}
# 			fig + theme(plot.title = element_text(size=10))
# 		}
# 	}	
# 	cbind(ggplotGrob(fl[[1]]), ggplotGrob(fl[[2]]), ggplotGrob(fl[[3]]), size='last')
# }
# grid.arrange(grobs = fig_df, ncol=1)

# fig_df = foreach(n_donor = donor_array) %do% {
# 	fl = foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			fig = resList[[idx]]$fig_fd + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + ylim(0, maxfd)  + theme(plot.title = element_text(hjust = 0.5, size=10), axis.text.x = element_text(size=8), aspect.ratio=1)
# 			if( n_reps > 2){
# 				fig = fig + theme(axis.title.y=element_blank(), 
# 		        axis.text.y=element_blank(),
# 		        axis.ticks.y=element_blank())
# 			}
# 			if( n_reps != 3){
# 				fig = fig + ylab('')
# 			}
# 			fig + theme(plot.title = element_text(size=10))
# 		}
# 	}	
# 	cbind(ggplotGrob(fl[[1]]), ggplotGrob(fl[[2]]), ggplotGrob(fl[[3]]), size='last')
# }
# grid.arrange(grobs = fig_df, ncol=1)
# dev.off()



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
aupr$method = factor(aupr$method, df_plots$method)

# Drop dream-EB
if( opt$noEB ){
	aupr = aupr[!(aupr$method %in% c("dream (KR)", "dream")),]
	aupr = droplevels(aupr)
}

col = df_plots$color[df_plots$method %in% levels(aupr$method)]

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
power_fdr_5$method = factor(power_fdr_5$method, df_plots$method)

# Drop dream-EB
if( opt$noEB ){
	power_fdr_5 = power_fdr_5[!(power_fdr_5$method %in% c("dream (KR)", "dream")),]
	power_fdr_5 = droplevels(power_fdr_5)
}

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
df_fpr$method = factor(df_fpr$method, df_plots$method)

# Drop dream-EB
if( opt$noEB ){
	df_fpr = df_fpr[!(df_fpr$method %in% c("dream (KR)", "dream")),]
	df_fpr = droplevels(df_fpr)
}

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
df_fd$method = factor(df_fd$method, df_plots$method)

if( opt$noEB ){
	df_fd = df_fd[!(df_fd$method %in% c("dream (KR)", "dream")),]
	df_fd = droplevels(df_fd)
}
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

if( opt$noEB ){
	df = df[!(df$key %in% c("dream (KR)", "dream")),]
	df = droplevels(df)
}

# saveRDS(df, file=paste0(opt$folder,'/df.RDS'))
# df = readRDS(paste0(opt$folder,'/df.RDS'))

# col = ggColorHue(length(table(df$key)))
col = df_plots$color[df_plots$method %in% levels(df$key)]

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

if( opt$noEB ){
	dfpr = dfpr[!(dfpr$key %in% c("dream (KR)", "dream")),]
	dfpr = droplevels(dfpr)
}

# saveRDS(dfpr, file=paste0(opt$folder,'/dfpr.RDS'))
# dfpr = readRDS(paste0(opt$folder,'/dfpr.RDS'))

# col = ggColorHue(length(table(df$key)))
col = df_plots$color[df_plots$method %in% levels(dfpr$method)]
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

if( opt$noEB ){
	df_fpr = df_fpr[!(df_fpr$key %in% c("dream (KR)", "dream")),]
	df_fpr = droplevels(df_fpr)
}

# saveRDS(df_fpr, file=paste0(opt$folder,'/df_fpr.RDS'))
# df_fpr = readRDS(paste0(opt$folder,'/df_fpr.RDS'))


file = paste0(folder,'/../figures/','combine_fpr', ".pdf")
pdf( file, width=15, height=20)
maxValue = max(df_fpr$value)

df_a = df_fpr[(n_donor <= 14),]
df_a$method = droplevels(df_a$method)
col = df_plots$color[df_plots$method %in% levels(df_a$method)]

fig1 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)

df_a = df_fpr[(n_donor > 14),]
if( nrow(df_a) > 0 ){
# drop macau2 if all empty
# if( df_a[method=="macau2",unique(value)] == 0 ){
# 	df_a = df_a[method!="macau2",]
# }
	df_a$method = droplevels(df_a$method)
	col = df_plots$color[df_plots$method %in% levels(df_a$method)]

	fig2 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)
	grid.arrange(fig1, fig2, ncol=2)
}else{
	fig1
}
dev.off()


# AUPR
######

df_aupr = do.call("rbind", lapply(resList, function(x) x$aupr))
df_aupr = data.table(df_aupr)

if( opt$noEB ){
	df_aupr = df_aupr[!(df_aupr$key %in% c("dream (KR)", "dream")),]
	df_aupr = droplevels(df_aupr)
}

# saveRDS(df_aupr, file=paste0(opt$folder,'/df_aupr.RDS'))
# df_aupr = readRDS(paste0(opt$folder,'/df_aupr.RDS'))

# col = ggColorHue(length(table(df$key)))

file = paste0(folder,'/../figures/','combine_aupr', ".pdf")
pdf( file, width=15, height=20)
maxValue = max(df_aupr[(n_donor <= 14),value])

df_a = df_aupr[(n_donor <= 14),]
df_a$method = droplevels(df_a$method)
col = df_plots$color[df_plots$method %in% levels(df_a$method)]
fig1 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("AUPR") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)

df_a = df_aupr[(n_donor > 14),]
if( nrow(df_a) > 0 ){
	df_a$method = droplevels(df_a$method)
	col = df_plots$color[df_plots$method %in% levels(df_a$method)]
	fig2 = ggplot(df_a, aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("AUPR") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, 1)
	grid.arrange(fig1, fig2, ncol=2)
}else{
	fig1
}
dev.off()








