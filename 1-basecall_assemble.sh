#! /bin/bash

WORKDIR=`realpath ./`
RAW_DATA_DIR=`realpath ./0-1-raw_data`
DATA_DIR=`realpath ./0-2-data`
THREADS=64

cd $WORKDIR

# Basecall all pod5 files in the raw-data directory and place the resulting BAM file in the data directory
dorado basecaller $RAW_DATA_DIR/ \
    -v dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v1,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1 \
    --batchsize=32 \
    --reference $DATA_DIR/pSMT3-m.fasta \
    > $DATA_DIR/basecalls.bam

# Index and sort the modBAM file using samtools
samtools sort $DATA_DIR/basecalls.bam -o $DATA_DIR/sorted-calls.bam
samtools index $DATA_DIR/sorted-calls.bam

# Transform basecalled reads into a fastq file and perform read correction through the HERRO algorithm
samtools fastq $DATA_DIR/sorted-calls.bam > $DATA_DIR/raw.fastq
dorado correct $DATA_DIR/raw.fastq > $DATA_DIR/corrected.fasta

# Generate read statistics through assembly-stats and Nanoplot
assembly-stats $DATA_DIR/raw.fastq > raw-reads.stats
NanoPlot -t $THREADS --color green --fastq $DATA_DIR/raw.fastq -o ./nanoplot-raw-reads

assembly-stats $DATA_DIR/corrected.fasta > corrected-reads.stats
NanoPlot -t $THREADS --color green --fasta $DATA_DIR/corrected.fasta -o ./nanoplot-corrected-reads

# Create genome size estimation through jellyfish k-mer counts
jellyfish count -t $THREADS -C -m 21 -s 5G -o $DATA_DIR/21mer_out $DATA_DIR/raw.fastq
jellyfish histo -o $DATA_DIR/21mer_out.histo $DATA_DIR/21mer_out

# TODO: Rscript or online?
# Rscript /data/6.PANEL_EXPERTS/luthfibio/apps/genomescope2.0-2.0.1/genomescope.R -i 15mer_out.histo -o 24_1_15mer -k 15

# Assemble the raw and corrected reads through Flye
mkdir $DATA_DIR/raw-assembly
flye --nano-raw "$DATA_DIR/raw.fastq" \
     --threads $THREADS \
     --out-dir "$DATA_DIR/raw-assembly"
cp $DATA_DIR/raw-assembly/assembly.fasta ./raw-assembly.fasta

mkdir $DATA_DIR/corrected-assembly 
flye --nano-corr "$DATA_DIR/corrected.fasta" \
    --threads $THREADS \
    --out-dir "$DATA_DIR/corrected-assembly"
cp $DATA_DIR/corrected-assembly/assembly.fasta ./corrected-assembly.fasta

# Generate assembly statistics through assembly-stats and BUSCO
assembly-stats raw-assembly.fasta > raw-assembly.stats
assembly-stats corrected-assembly.fasta > corrected-assembly.stats

busco -i raw-assembly.fasta -m genome -l actinopterygii_odb10 -o $DATA_DIR/BUSCO-raw-assembly -c $THREADS
busco -i corrected-assembly.fasta -m genome -l actinopterygii_odb10 -o $DATA_DIR/BUSCO-corrected-assembly -c $THREADS

# Create snail plots through blobtoolkit
mkdir $DATA_DIR/plot-raw-assembly
blobtools create --fasta raw-assembly.fasta --busco $DATA_DIR/BUSCO-raw-assembly/run_actinopterygii_odb10/full_table.tsv $DATA_DIR/plot-raw-assembly
blobtools view --plot --view snail $DATA_DIR/plot-raw-assembly

mkdir $DATA_DIR/plot-corrected-assembly
blobtools create --fasta corrected-assembly.fasta --busco $DATA_DIR/BUSCO-corrected-assembly/run_actinopterygii_odb10/full_table.tsv $DATA_DIR/plot-corrected-assembly
blobtools view --plot --view snail $DATA_DIR/plot-corrected-assembly
