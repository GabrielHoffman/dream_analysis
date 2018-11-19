# Gabriel Hoffman
# Icahn School Of Medicine at Mount Sinai
# November 16, 2018

# Get genetic Rsq values for every gene from GTEx whole blood
# Get data from: https://gtexportal.org/home/datasets

library(data.table)

# read SNP level eQTL statistics
df_snp = fread('~/Downloads/GTEx_Analysis_v7_eQTL/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz')

# compute r2 for each gene-SNP pair
df_snp[,r2 :=slope^2*(maf*(1-maf))]

# get highest r2 for each gene
df_best = df_snp[,data.frame(r2=max(r2)), by=gene_id]

# get ENSEMBL gene prefix
df_best$gene_id = sapply(strsplit(df_best$gene_id, '\\.'), function(x) x[1])

# use ENSEMBL v75, but can be other ENSEMBL transcript database
library(EnsDb.Hsapiens.v75)
txdb = EnsDb.Hsapiens.v75

cat("Reading gene locations...\n")

# Define map of genes to ENSEMBL ids
target = genes( txdb, columns = c("gene_id", "gene_name")) #

# add gene names to df_best
idx = match(df_best$gene_id, target$gene_id)
df_best$gene_name = target$gene_name[idx]

write.table( df_best, file="./eQTL_R2_GTEx_Whole_Blood.tsv", row.names=FALSE, quote=FALSE)


# get CommonMind eQTL
##############
library(data.table)

# eQTL data
dfeqtl = fread(cmd = 'cat /sc/orga/projects/CommonMind/results/phaseI/eQTL/CMC_EQTL_FROM_SYNAPSE/SVA/Caucasian/Combined/CMC_eqtl_SVA_Caucasian_combined_*_cis.out ')
colnames(dfeqtl)[colnames(dfeqtl) == 't-stat'] = 'tstat'
colnames(dfeqtl)[colnames(dfeqtl) == 'p-value'] = 'pvalue'
dfeqtl$beta = as.numeric(dfeqtl$beta)
dfeqtl$pvalue = as.numeric(dfeqtl$pvalue)
dfeqtl$tstat = as.numeric(dfeqtl$tstat)


# MAF
dmaf = fread("/sc/orga/projects/CommonMind/results/phaseI/eQTL/CMC_EQTL_FROM_SYNAPSE/SNP_positions_A1freq.txt.gz")
# dmaf$SNP = as.character(dmaf$SNP)

# get best SNP per gene 
df_eqtl_best = dfeqtl[,.SD[which.min(pvalue)], by = gene]

# merge eQTL and MAF
dfsub_merge = merge(df_eqtl_best, dmaf, by="SNP")

# compute R2 = beta^2*(maf*(1-maf))
df_final = dfsub_merge[,r2:=beta^2*(Caucasian*(1-Caucasian))]



with(df_final, cor(r2, -log10(pvalue), method="spearman"))
with(df_final, cor(r2, abs(tstat), method="spearman"))
with(df_final, cor(abs(tstat), -log10(pvalue), method="spearman"))











