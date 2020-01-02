# Gabriel Hoffman
# December 4, 2018

#######################
# Write eQTL R2 table #
#######################

library(RSQLite)
library(biomaRt)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
getAtt = c("ensembl_gene_id","hgnc_symbol")

# identify PrediXcan databases
files = c(CMC = '/sc/hydra/projects/psychgen/resources/predixcan/CMC/CMC_DLPFC-PrediXcan_Models/DLPFC_newMetax.db',
	DGN = '/sc/hydra/projects/psychgen/resources/predixcan/DGN/predb00000059.db'	)

# read CMC R2
#############
con <- dbConnect(RSQLite::SQLite(), files['CMC'])
df_cmc = dbGetQuery(con, "SELECT gene AS ensembl_gene_id, [pred.perf.R2] AS CMC_R2 FROM extra")
dbDisconnect(con)

# get gene identifiers
geneInfoBiomart = getBM(getAtt,filters="ensembl_gene_id", values=df_cmc$ensembl_gene_id ,mart = ensembl)
df_cmc = merge(df_cmc, geneInfoBiomart, by='ensembl_gene_id')


# read DGN R2
#############
con2 <- dbConnect(RSQLite::SQLite(), files['DGN'])
df_dgn = dbGetQuery(con2, "SELECT genename AS hgnc_symbol, R2 AS DGN_R2 FROM extra")
dbDisconnect(con2)

# get gene identifiers
geneInfoBiomart = getBM(getAtt,filters="hgnc_symbol", values=df_dgn$hgnc_symbol ,mart = ensembl)
df_dgn = merge(df_dgn, geneInfoBiomart, by='hgnc_symbol')


# read GTEx R2
##############
# con3 <- dbConnect(RSQLite::SQLite(), files['GTEx'])
# df_gtex = dbGetQuery(con3, "SELECT genename AS hgnc_symbol, R2 AS GTEx_R2 FROM extra")
# dbDisconnect(con3)

# # get gene identifiers
# geneInfoBiomart = getBM(getAtt,filters="hgnc_symbol", values=df_dgn$hgnc_symbol ,mart = ensembl)
# df_gtex = merge(df_gtex, geneInfoBiomart, by='hgnc_symbol')

# merge
df_eqtl = merge(df_dgn, df_cmc, by='ensembl_gene_id', all=TRUE)
df_eqtl$DGN_R2[is.na(df_eqtl$DGN_R2)] = 0
df_eqtl$CMC_R2[is.na(df_eqtl$CMC_R2)] = 0

geneInfoBiomart = getBM(getAtt,filters="ensembl_gene_id", values=df_eqtl$ensembl_gene_id, mart = ensembl)
df_eqtl2 = merge(df_eqtl[,c('ensembl_gene_id', 'DGN_R2', 'CMC_R2')], geneInfoBiomart, by='ensembl_gene_id')

df_eqtl2 = df_eqtl2[,c('ensembl_gene_id', 'hgnc_symbol', colnames(df_eqtl2)[-c(1, ncol(df_eqtl2))]),]

file = "eqtl_r2.tsv"
write.table(df_eqtl2, file=file, row.names=FALSE, sep="\t", quote=FALSE)

system('module load python; synapse add --parentid syn16816470 eqtl_r2.tsv')
file.remove(file)




# cat ~/scripts/dream/df_cor.tsv /hpc/users/hoffmg01/work/dev_dream/dream_analysis/src/df_cor.tsv /hpc/users/hoffmg01/work/dev_dream/dream_analysis/


library(openxlsx)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

# df = read.xlsx('/Users/gabrielhoffman/Dropbox/dream/dream_GenBio_v2/table.xlsx')
df = read.xlsx('/Users/gabrielhoffman/Dropbox/dream/dream_GenBiol_v4/eQTL_table.xlsx')
df$i = rev(1:nrow(df))
df$label = with(df, paste(Dataset, Trait, Tissue))
df$label = factor(df$label, rev(df$label))
i = grep('^<', df$CMC_pValue, invert=TRUE)
df$CMC_pValue[i] = format(as.numeric(df$CMC_pValue[i]), digits=2)
i = grep('^<', df$DGN_pValue, invert=TRUE)
df$DGN_pValue[i] = format(as.numeric(df$DGN_pValue[i]), digits=2)


df2 = melt(df, id.vars=c('Dataset', 'Trait', 'Tissue', 'i', 'label', 'CMC_pValue', 'DGN_pValue'))

fig = ggplot(df2, aes(label, value, fill=variable)) + geom_bar(stat='identity', position='dodge') + theme_bw(20) + theme(aspect.ratio=1, axis.text.y=element_text(size=10)) + coord_flip() + scale_y_continuous(expand=c(0,0), lim=c(0,1)) + scale_fill_manual(values = c("#5580ff", "#ff556c")) 

gr  = tableGrob(df[,c("Dataset", "Trait", "Tissue", "CMC_pValue", 'DGN_pValue')])

fig2 = arrangeGrob(gr, fig, ncol=2,  as.table=TRUE)




file = '/Users/gabrielhoffman/Dropbox/dream/dream_GenBio_v2/Figure_full/orig/barplot.pdf'
ggsave(file, fig2, width=20)











