#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(getopt))

a = matrix(c(
	'folder', 		'f', 1, "character",
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
	# de_res = de_res[,!(colnames(de_res) %in% c("lmm_Sat", "lmm_KR", "lmm_KR_eBayes" ))]

	# raw p-values on right
	keepRaw = TRUE
	if( keepRaw ){
		de_res_right = de_res[,colnames(de_res) %in% c("lmm_Sat", "lmm_KR")]
		de_res = de_res[,!(colnames(de_res) %in% c("lmm_Sat", "lmm_KR"))]
		de_res = cbind(de_res, de_res_right)

		# specify colors
		method_order = colnames(de_res)[-c(1, 12,13)]
		col = ggColorHue(length(method_order))
		names(col) = method_order
		col['lmm_Sat'] = alpha(col['dream'], .5)
		col['lmm_KR'] = alpha(col['dream (KR)'], .5)
		method_order = c(method_order, c("lmm_Sat", "lmm_KR"))
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
	df$key = factor(df$key, method_order)
	df$n_chosen = df$n_chosen / length(files)
	df$n_false = df$n_false / length(files)

	# plot false positive versus positives
	fig_choose = ggplot(df, aes(n_chosen, n_false, color=key)) + geom_line() + xlim(0, 375) + ylim(0, 50) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=col) + xlab("# genes selected") + ylab("# false positives")

	# Compute PR
	#===========
	prList = foreach( method = colnames(de_res)[-1] ) %do% {
		pr.curve(	scores.class0 = -log10(de_res[de_res$true==1,method]), 
					scores.class1 = -log10(de_res[de_res$true==0,method]), 
					curve=TRUE, rand.compute = TRUE)
	}
	names(prList) = colnames(de_res)[-1] 

	aupr = data.frame(method = names(prList) )
	aupr$value = foreach( method = names(prList), .combine=c ) %do% {
		prList[[method]]$auc.integral
	}
	aupr$method = factor(aupr$method, method_order)
	aupr$n_donor = n_donor
	aupr$n_reps = n_reps
	aupr.rand.score = prList[[method]]$rand$auc.integral

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
	dfpr$method = factor(dfpr$method, method_order)

	rnd.value = prList[[1]]$rand$curve[1,2]
	dfpr$rnd.value = rnd.value 
	dfpr$n_donor = n_donor
	dfpr$n_reps = n_reps

	fig_pr = ggplot(dfpr, aes(recall, precision, color=method)) + geom_line() + theme_bw(12) + xlab("Recall") + ylab("Precision") + xlim(0,1) + ylim(0,1) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=rnd.value, color="grey50", linetype=2) + geom_hline(yintercept=aupr.rand.score, linetype=2)
	
	# Power	at 5% FDR
	#===========
	power_fdr_5 = data.frame(method = colnames(de_res)[-1])
	power_fdr_5$value = foreach( method = names(prList), .combine=c ) %do% {
		pr = as.data.frame(prList[[method]]$curve)
		colnames(pr) = c( "recall", "precision", "score")
		# ggplot(pr, aes(recall, precision)) + geom_line()

		i = which.min(abs(pr$precision - 0.95))
		pr$recall[i]
	}
	power_fdr_5$method = factor(power_fdr_5$method, method_order)
	power_fdr_5$n_donor = n_donor
	power_fdr_5$n_reps = n_reps

	fig_power_fdr_5 = ggplot(power_fdr_5, aes(method, value, fill=method)) + geom_bar(stat="identity") + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=col) + ylab("Power at FDR 5%") + coord_flip()

	# False positive rate
	#======================
	fpr = apply(de_res[de_res$true ==0,], 2, function(x) sum(x< 0.05)/length(x))
	df_fpr = data.frame(method=names(fpr)[-1])
	df_fpr$value = fpr[-1]
	df_fpr$method = factor(df_fpr$method, method_order)
	df_fpr$n_donor = n_donor
	df_fpr$n_reps = n_reps

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
	
	# de_res = de_res[,!(colnames(de_res) %in% c("lmm_Sat", "lmm_KR", "lmm_KR_eBayes" ))]

	fd = apply(de_res[de_res$true ==0,], 2, function(x) sum(x< 0.05))
	false_discoveries = data.frame(method=names(fd)[-1])
	# divide false discoveries by number of analyses
	false_discoveries$value = as.numeric(fd[-1] / length(files))
	false_discoveries$method = factor(false_discoveries$method, method_order)
	false_discoveries$n_donor = n_donor
	false_discoveries$n_reps = n_reps

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
# # choose
# fig_choose = foreach(n_donor = donor_array, .combine=c) %do% {
# 	foreach(n_reps = 2:4, .combine=c) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)) ) %do% {
# 			resList[[idx]]$fig_choose + theme(legend.position="none") + ggtitle(paste('# donors:', n_donor, ' # reps:', n_reps)) + theme(plot.title = element_text(hjust = 0.5, size=10))
# 		}
# 	}	
# }
# # do.call("grid.arrange", c(fig_choose, ncol=3))

# # choose
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
# pdf( file, width=7, height=16 )
# # AUPR
# aupr = foreach(n_donor = donor_array, .combine=rbind) %do% {
# 	foreach(n_reps = 2:4, .combine=rbind) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
# 			resList[[idx]]$aupr 
# 		}
# 	}	
# }
# fig = foreach(n_reps = 2:4) %do%{
# ggplot(aupr[aupr$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_fill_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, 1) + ylab("AUPR")
# }
# fig_aupr = do.call("grid.arrange", c(fig, ncol=3))

# # power_fdr_5
# power_fdr_5 = foreach(n_donor = donor_array, .combine=rbind) %do% {
# 	foreach(n_reps = 2:4, .combine=rbind) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
# 			resList[[idx]]$power_fdr_5
# 		}
# 	}	
# }
# fig = foreach(n_reps = 2:4) %do%{
# ggplot(power_fdr_5[power_fdr_5$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_fill_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, 1) + ylab("Power at FDR 5%")
# }
# fig_power_fdr_5 = do.call("grid.arrange", c(fig, ncol=3))

# df_fpr = foreach(n_donor = donor_array, .combine=rbind) %do% {
# 	foreach(n_reps = 2:4, .combine=rbind) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
# 			resList[[idx]]$df_fpr	
# 		}
# 	}	
# }
# fig = foreach(n_reps = 2:4) %do%{
# ggplot(df_fpr[df_fpr$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_hline(yintercept=0.05, linetype=2) + geom_line() + scale_fill_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, .2) + ylab("False positive rate") 
# }
# fig_fpr = do.call("grid.arrange", c(fig, ncol=3))

# # False discoveries
# df_fd = foreach(n_donor = donor_array, .combine=rbind) %do% {
# 	foreach(n_reps = 2:4, .combine=rbind) %do% {
# 		foreach( idx = grep(paste0('^',n_donor, '_',n_reps, '_', ".*"), names(resList)), .combine=rbind ) %do% {
# 			resList[[idx]]$false_discoveries	
# 		}
# 	}	
# }
# fig = foreach(n_reps = 2:4) %do%{
# ggplot(df_fd[df_fd$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_fill_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, max(df_fd$value)) + ylab("False discoveries at FDR < 5%") 
# }
# fig_fd = do.call("grid.arrange", c(fig, ncol=3))

# fig = foreach(n_reps = 2:4) %do%{
# 	mvalue = max(df_fd[df_fd$n_reps==n_reps,]$value)
# ggplot(df_fd[df_fd$n_reps==n_reps,], aes(n_donor, value, color=method)) + geom_line() + scale_fill_manual(values=col) + theme_bw(12) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlim(0, max(donor_array)) + ylim(0, mvalue) + ylab("False discoveries at FDR < 5%") 
# }
# fig_fd2 = do.call("grid.arrange", c(fig, ncol=3))
# dev.off()

# pdf( file, width=8, height=12 )
# grid.arrange(fig_aupr, fig_power_fdr_5, fig_fpr, fig_fd, fig_fd2,ncol=1)
# dev.off()

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

# saveRDS(df, file=paste0(opt$folder,'/df.RDS'))
# df = readRDS(paste0(opt$folder,'/df.RDS'))

col = ggColorHue(length(table(df$key)))

file = paste0(folder,'/../figures/','combine_choose', ".pdf")
pdf( file, width=7, height=20)
fig1 = ggplot(df[(n_donor <= 14),], aes(n_chosen, n_false, color=key)) + geom_line() + facet_grid( n_donor ~ n_reps) + xlim(0, 300) + ylim(0, 50) + theme_bw(10) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size=8), legend.position="bottom") + scale_color_manual(values=col) + xlab("genes selected") + ylab("false discoveries")
fig2 = ggplot(df[(n_donor > 14),], aes(n_chosen, n_false, color=key)) + geom_line() + facet_grid( n_donor ~ n_reps) + xlim(0, 550) + ylim(0, 50) + theme_bw(10) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), axis.text.x =  element_text(size=8), legend.position="bottom") + scale_color_manual(values=col) + xlab("genes selected") + ylab("false discoveries")
grid.arrange(fig1, fig2, ncol=2)
dev.off()

# Precision recall
##################

dfpr = do.call("rbind", lapply(resList, function(x) x$dfpr))
dfpr = data.table(dfpr)

# saveRDS(dfpr, file=paste0(opt$folder,'/dfpr.RDS'))
# dfpr = readRDS(paste0(opt$folder,'/dfpr.RDS'))

col = ggColorHue(length(table(df$key)))
randCurve = dfpr[,unique(rnd.value)]

file = paste0(folder,'/../figures/','combine_pr2', ".pdf")
pdf( file, width=8, height=15)
fig1 = ggplot(dfpr[(n_donor <= 14),], aes(recall, precision, color=method)) + geom_line()  + theme_bw(10) + facet_grid( n_donor ~ n_reps) + xlim(0,1) + ylim(0,1) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="bottom", axis.text.x=element_text(size=6)) + scale_color_manual(values=col) + geom_hline(yintercept=randCurve, linetype=2)
fig2 = ggplot(dfpr[(n_donor > 14),], aes(recall, precision, color=method)) + geom_line()  + theme_bw(10) + facet_grid( n_donor ~ n_reps) + xlim(0,1) + ylim(0,1) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="bottom", axis.text.x=element_text(size=6)) + scale_color_manual(values=col) + geom_hline(yintercept=randCurve, linetype=2)
grid.arrange(fig1, fig2, ncol=2)
dev.off()

# False positive rate
######################

df_fpr = do.call("rbind", lapply(resList, function(x) x$df_fpr))
df_fpr = data.table(df_fpr)

# saveRDS(df_fpr, file=paste0(opt$folder,'/df_fpr.RDS'))
# df_fpr = readRDS(paste0(opt$folder,'/df_fpr.RDS'))

col = ggColorHue(length(table(df$key)))

file = paste0(folder,'/../figures/','combine_fpr', ".pdf")
pdf( file, width=15, height=20)
maxValue = max(df_fpr$value)
fig1 = ggplot(df_fpr[(n_donor <= 14),], aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)
fig2 = ggplot(df_fpr[(n_donor > 14),], aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("False positive rate at p<0.05") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)
grid.arrange(fig1, fig2, ncol=2)
dev.off()


# AUPR
######

df_aupr = do.call("rbind", lapply(resList, function(x) x$aupr))
df_aupr = data.table(df_aupr)

# saveRDS(df_aupr, file=paste0(opt$folder,'/df_aupr.RDS'))
# df_aupr = readRDS(paste0(opt$folder,'/df_aupr.RDS'))

col = ggColorHue(length(table(df$key)))

file = paste0(folder,'/../figures/','combine_aupr', ".pdf")
pdf( file, width=15, height=20)
maxValue = max(df_aupr[(n_donor <= 14),value])
fig1 = ggplot(df_aupr[(n_donor <= 14),], aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("AUPR") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, maxValue)
fig2 = ggplot(df_aupr[(n_donor > 14),], aes(method, value, fill=method)) + geom_bar(stat="identity") + geom_hline(yintercept=0.05, linetype=2) + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x=element_text(size=8)) + scale_fill_manual(values=col) + ylab("AUPR") + coord_flip()  + facet_grid( n_donor ~ n_reps) + ylim(0, 1)
grid.arrange(fig1, fig2, ncol=2)
dev.off()








