# Gabriel Hoffman
# August 5, 2018
# 
# Code for simulating repeated measures RNA-seq
# This code is generates some basic files, writes job scripts that
# 	can be submitted to a cluster, and makes some plots of the results
# This is a mix of bash and R-code, so each section must be pasted
# 	into the console were approporiate.


# Given a list of all protein coding transcripts
# extract first transcript per gene
# all transcripts
DIR=/sc/orga/projects/psychencode/gabriel/RNA_seq_sim/RNA_seq_sim_v3/
cd $DIR
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
gunzip gencode.v19.pc_transcripts.fa.gz
FASTALL=$DIR/transcriptome/gencode.v19.pc_transcripts.fa
FASTA=$DIR/transcriptome/gencode.v19.genes.fa
cd $(dirname $FASTA)

# unique genes
cat $FASTALL | grep '>' | cut -f2 -d'|' | cut -f1 -d'.' | sort -u | parallel -P10 "grep -m1 {} $FASTALL" | sed 's/^>//g' > headers.lst

# get genes greater than 400 bp and less thatn 16 kb
cat headers.lst | awk -vFS='|' '{if($7>500 && $7 < 8000) print $0}' > headers_size.lst

# faSomeRecords $FASTALL headers.lst $FASTA
faSomeRecords $FASTALL headers_size.lst $FASTA


# Simulate read counts
#######################

FOLDER=/sc/orga/projects/psychencode/gabriel/RNA_seq_sim/RNA_seq_sim_v3/
cd $FOLDER
# rm -rf data figures jobs logs results
mkdir -p data figures jobs logs results
FASTA=$FOLDER/transcriptome/gencode.v19.genes.fa

# save full name of sim
N_DE=500
FC=3
HSQ=0.4

LOG=$FOLDER/logs/
# for N_SAMPLES in $(echo $(seq 4 2 20) 30 40 50);
for N_SAMPLES in $(echo $(seq 8 2 14));
do
for N_REPS in $(seq 2 4);
do
# for SEED in $(seq 1 50);
for SEED in $(seq 1 5);
do
PFX=${N_SAMPLES}_${N_REPS}_${N_DE}_${FC}_${HSQ}_${SEED}
echo '#!/bin/bash' > jobs/sims_${PFX}.lsf
echo "#BSUB -J ${PFX}
#BSUB -P acc_psychencode
#BSUB -q premium
#BSUB -n 3
#BSUB -R span[hosts=1]
#BSUB -W 12:00
#BSUB -o $LOG/sims_${PFX}_%J.stdout
#BSUB -eo $LOG/sims_${PFX}_%J.stderr
#BSUB -L /bin/bash

module purge
module load R/3.5.3

# load local version of variancePartition
R_LIBS=\$R_LIBS_USER:\$R_LIBS

export OMP_NUM_THREADS=1

echo \$( R -e \"packageVersion('variancePartition')\")

/hpc/users/hoffmg01/work/dev_dream/dream_analysis/sims/generate_simulations.R --fasta $FASTA --n_samples ${N_SAMPLES} --n_reps ${N_REPS} --n_de_genes ${N_DE} --disease_fc ${FC} --nthreads 20  --param_ID '.3 .03' --param_Disease '.2 0.005' --param_Batch '.2 0.01' --seed ${SEED} --out $FOLDER/data --prefix $PFX " >> jobs/sims_${PFX}.lsf
done
done
done

ls jobs/sims_*lsf | parallel -P1 "bsub < {}; sleep .2"


# resubmit crashed jobs
# completed
ls data/countMatrix*.RDS | parallel -P1 basename {} .RDS | sed 's/countMatrix_//g' > complete.lst

# running
bjobs -w | grep -v JOB_NAME | awk '{print $7}' > running.lst

# all
ls jobs/sims_*lsf | parallel -P1 basename {} .lsf | sed 's/sims_//g' > all.lst

# submit remaining jobs
comm -23 <(sort all.lst) <(cat running.lst complete.lst | sort) | parallel -P1 ls jobs/sims_{}.lsf | parallel -P1 "bsub < {}; sleep .2"


# Run differential expression #
###############################

# cd /hpc/users/hoffmg01/work/dev_dream/dream_analysis/sims
# git pull

 # \rm -f figures/* jobs/* logs/* results/* data/*

N_DE=500
FC=3
HSQ=0.4

FOLDER=/sc/orga/projects/psychencode/gabriel/RNA_seq_sim/RNA_seq_sim_v3/
cd $FOLDER

LOG=$FOLDER/logs

# for N_SAMPLES in $(echo $(seq 4 2 20) 30 40 50);
for N_SAMPLES in $(echo $(seq 8 2 14));
do
EXTRA='--macau2'
if [ ${N_SAMPLES} -lt 14 ];
then
 EXTRA='--macau2'
fi
for N_REPS in $(seq 2 4);
do
# for SEED in $(seq 1 50);
for SEED in $(seq 1 5);
do
PFX=${N_SAMPLES}_${N_REPS}_${N_DE}_${FC}_${HSQ}_${SEED}
echo '#!/bin/bash' > jobs/scripts_${PFX}.lsf
echo "#BSUB -J ${PFX}
#BSUB -P acc_psychencode
#BSUB -q premium
#BSUB -n 6
#BSUB -R span[hosts=1]
#BSUB -W 8:00
#BSUB -o $LOG/scripts_${PFX}_%J.stdout
#BSUB -eo $LOG/scripts_${PFX}_%J.stderr
#BSUB -L /bin/bash

module purge
module load R/3.5.3

# load local version of variancePartition
R_LIBS=\$R_LIBS_USER:\$R_LIBS

export OMP_NUM_THREADS=1

echo \$( R -e \"packageVersion('variancePartition')\")

/hpc/users/hoffmg01/work/dev_dream/dream_analysis/sims/run_DE_analysis.R --prefix ${PFX} --folder $FOLDER $EXTRA " >> jobs/scripts_${PFX}.lsf
done
done
done

ls jobs/scripts_*lsf | parallel -P1 "bsub < {}; sleep .2"
	
	
# only sims 1-10
seq 1 10 | parallel -P1 ls jobs/scripts_*_{}.lsf | parallel -P1 "bsub < {}; sleep .2"


# 1) show eBayes vs raw in the supplement
# 2) time course


# resubmit crashed jobs
# completed
ls results/*_p.RDS | parallel -P1 basename {} _p.RDS > complete.lst

# running
bjobs -w | grep -v JOB_NAME | awk '{print $7}' > running.lst

# all
ls jobs/scripts_*lsf | parallel -P1 basename {} .lsf | sed 's/scripts_//g' > all.lst

# submit remaining jobs
comm -23 <(sort all.lst) <(cat running.lst complete.lst | sort) | parallel -P1 ls jobs/scripts_{}.lsf | parallel -P1 "bsub < {}; sleep .2"


# Make plots
###############

FOLDER=/sc/orga/projects/psychencode/gabriel/RNA_seq_sim/RNA_seq_sim_v3/
cd $FOLDER

/hpc/users/hoffmg01/work/dev_dream/dream_analysis/sims/make_plots.R --folder $FOLDER/results --nthreads 32



# Time methods
##############
qFOLDER=/hpc/users/hoffmg01/work/RNA_seq_sim_v1
SCRIPT=/hpc/users/hoffmg01/scripts/varPartSims/time_methods.R
LOG=$FOLDER/logs
for N_SAMPLES in $(echo $(seq 4 2 20) 30 40 50);
# for N_SAMPLES in $(seq 100 100 600);
do
for N_REPS in $(seq 2 4);
do
for SEED in $(seq 1 2);
do
for NTHREADS in $(echo 1 12);
do
MEM=$((8*$NTHREADS))
MEM=$(($MEM > 60 ? 60 : $MEM ))
MEM=$(($MEM/$NTHREADS*1000))
PFX=${N_SAMPLES}_${N_REPS}_${NTHREADS}_${SEED}
echo '#!/bin/bash' > jobs/times_${PFX}.lsf
echo "#BSUB -J ${PFX}
#BSUB -P acc_psychencode
#BSUB -q premium
#BSUB -n $NTHREADS
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=$MEM]
#BSUB -W 96:00
#BSUB -o $LOG/times_${PFX}_%J.stdout
#BSUB -eo $LOG/times_${PFX}_%J.stderr
#BSUB -L /bin/bash

module purge
module load R/3.4.3

$SCRIPT --folder $FOLDER --prefix ${N_SAMPLES}_${N_REPS}_500_3_0.4_1 --runkr --nthreads 1 --out ${N_SAMPLES}_${N_REPS}_${NTHREADS}_${SEED} --folder $FOLDER" >> jobs/times_${PFX}.lsf
done
done
done 
done


ls jobs/times_* | parallel -P1 "bsub < {}; sleep .2"



ls jobs/times_* | grep -v "00" |  parallel -P1 "bsub < {}; sleep .2"










# plot times
#############

library(foreach)
library(data.table)
library(reshape2)
library(ggplot2)
files = dir('/hpc/users/hoffmg01/work/RNA_seq_sim_v1/runTimes/', "*.tsv", full.names=TRUE)

get_idx = function( x, sep, i){
	sapply(strsplit(x, sep), function(x) x[i])
}

df = foreach( file = files, .combine=rbind) %do% {
	df = fread( file )
	df = melt(df)
	colnames(df) = c("method", "time")
	id = gsub(".tsv", "", basename(file))
	df$donors = as.numeric(get_idx( id, '_', 1))
	df$reps = as.numeric(get_idx( id, '_', 2))
	df$threads = as.numeric(get_idx( id, '_', 3))
	# df$seed = get_idx( id, '_', 4)
	df[,prefix:=paste(donors, reps, threads)]
	df
}
df = data.table(df)

df2 = df[,data.frame(time = mean(time),
				     donors, reps, threads), by=c('method', "prefix")]
df2 = unique(df2)

xmax = max(df2$donors)

# pdf("/hpc/users/hoffmg01/work/RNA_seq_sim_v1/run_time.pdf", height=5, width=5)
# ggplot(df2, aes(donors, log2(time), color=method)) + geom_line() + facet_wrap(~ threads+ reps) + theme_bw(10) + theme(aspect.ratio=1) + xlim(0, xmax) + ylab("Run time (minutes)")
# dev.off()

cols = c( '#00b6eb', '#a58aff', "#fb61d7")

pdf("/hpc/users/hoffmg01/work/RNA_seq_sim_v1/run_time.pdf", height=5, width=5)
breaks = -1:5
labels = 10^breaks
ggplot(df, aes(donors, log10(time), color=method)) + geom_point(size=.2) + facet_wrap(~ threads+ reps) + theme_bw(10) + theme(aspect.ratio=1) + xlim(0, xmax) + ylab("Run time (minutes)") + geom_smooth(se=FALSE, size=.6, span=1) + scale_y_continuous( breaks=breaks, labels=labels) + scale_color_manual(values=cols)

ggplot(df, aes(donors, time, color=method)) + geom_point(size=.2) + facet_wrap(~ threads+ reps) + theme_bw(10) + theme(aspect.ratio=1) + xlim(0, xmax) + ylab("Run time (minutes)") + geom_smooth(se=FALSE, size=.6, span=1) + scale_color_manual(values=cols)
dev.off()


df[which.max(time),]
df[which.max(time)-1,]












