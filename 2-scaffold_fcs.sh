#! /bin/bash

WORKDIR=`realpath ./`
RAW_DATA_DIR=`realpath ./0-1-raw_data`
DATA_DIR=`realpath ./0-2-data`
THREADS=64

RAW_ASSEMBLY="$WORKDIR/raw-assembly.fasta"
CORRECTED_ASSEMBLY="$WORKDIR/corrected-assembly.fasta"
DANRER_REFSEQ="$DATA_DIR/"  # TODO

cd $WORKDIR

# Scaffold assemblies based on the zebrafish reference genome
ragtag.py scaffold $DANRER_REFSEQ $RAW_ASSEMBLY
mv ragtag_output $DATA_DIR/raw-scaffold
cp $DATA_DIR/raw-scaffold/ragtag.scaffold.fasta raw-scaffold.fasta

ragtag.py scaffold $DANRER_REFSEQ $CORRECTED_ASSEMBLY
mv ragtag_output $DATA_DIR/corrected-scaffold
cp $DATA_DIR/corrected-scaffold/ragtag.scaffold.fasta corrected-scaffold.fasta

# Generate assembly statistics through assembly-stats and BUSCO
assembly-stats raw-scaffold.fasta > raw-scaffold.stats
assembly-stats corrected-scaffold.fasta > corrected-scaffold.stats

busco -i raw-scaffold.fasta -m genome -l actinopterygii_odb10 -o $DATA_DIR/BUSCO-raw-scaffold -c $THREADS
busco -i corrected-scaffold.fasta -m genome -l actinopterygii_odb10 -o $DATA_DIR/BUSCO-corrected-scaffold -c $THREADS

# Create snail plots through blobtoolkit
mkdir $DATA_DIR/plot-raw-scaffold
blobtools create --fasta raw-scaffold.fasta --busco $DATA_DIR/BUSCO-raw-scaffold/run_actinopterygii_odb10/full_table.tsv $DATA_DIR/plot-raw-scaffold
blobtools view --plot --view snail $DATA_DIR/plot-raw-scaffold

mkdir $DATA_DIR/plot-corrected-scaffold
blobtools create --fasta corrected-scaffold.fasta --busco $DATA_DIR/BUSCO-corrected-scaffold/run_actinopterygii_odb10/full_table.tsv $DATA_DIR/plot-corrected-scaffold
blobtools view --plot --view snail $DATA_DIR/plot-corrected-scaffold

# Perform FCS through the NCBI-FCS pipeline
# TODO: Scripts or Galaxy?
# cp ? ./clean-scaffold.fasta

# Generate assembly statistics through assembly-stats and BUSCO
assembly-stats clean-scaffold.fasta > clean-scaffold.stats
busco -i clean-scaffold.fasta -m genome -l actinopterygii_odb10 -o $DATA_DIR/BUSCO-clean-scaffold -c $THREADS

# Create a snail plot through blobtoolkit
mkdir $DATA_DIR/plot-clean-scaffold
blobtools create --fasta clean-scaffold.fasta --busco $DATA_DIR/BUSCO-clean-scaffold/run_actinopterygii_odb10/full_table.tsv $DATA_DIR/plot-clean-scaffold
blobtools view --plot --view snail $DATA_DIR/plot-clean-scaffold
