dir=/sc1/groups/pls-redbfx/pipeline_runs/iPETE/production/2019-03-12-Expt56_nextseq_testing/analysis/Exp56_PBMC_222_N714_S1/seq/trimmomatic
fq1=${dir}/Exp56_PBMC_222_N714_S1_R1_001_quality_filtered.fastq
fq2=${dir}/Exp56_PBMC_222_N714_S1_R2_001_quality_filtered.fastq

./bin/pair_consensus -1 ${fq1} -2 ${fq2} -p 1 -s 100 -n 100 
