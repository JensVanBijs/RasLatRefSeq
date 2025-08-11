#! /bin/bash

WORKDIR=`realpath ./`
RAW_DATA_DIR=`realpath ./0-1-raw_data`
DATA_DIR=`realpath ./0-2-data`
THREADS=64

DANRER_REFSEQ="$DATA_DIR/"  # TODO
AUGUSTUS_PATH="/home/jens/.augustus"  # TODO

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

# Transform basecalled reads into a fastq file
samtools fastq $DATA_DIR/sorted-calls.bam > $DATA_DIR/raw.fastq

# Generate read statistics through assembly-stats and Nanoplot
assembly-stats $DATA_DIR/raw.fastq > raw-reads.stats
NanoPlot -t $THREADS --color green --fastq $DATA_DIR/raw.fastq -o ./nanoplot-raw-reads

# Create genome size estimation through jellyfish k-mer counts
jellyfish count -t $THREADS -C -m 21 -s 5G -o $DATA_DIR/21mer_out $DATA_DIR/raw.fastq
jellyfish histo -o $DATA_DIR/21mer_out.histo $DATA_DIR/21mer_out

# TODO: Rscript or online?
# Rscript /data/6.PANEL_EXPERTS/luthfibio/apps/genomescope2.0-2.0.1/genomescope.R -i 15mer_out.histo -o 24_1_15mer -k 15

# Assemble the reads through Flye
mkdir $DATA_DIR/raw-assembly
flye --nano-raw "$DATA_DIR/raw.fastq" \
     --threads $THREADS \
     --out-dir "$DATA_DIR/raw-assembly"
cp $DATA_DIR/raw-assembly/assembly.fasta ./raw-assembly.fasta

# Generate assembly statistics through assembly-stats and BUSCO
assembly-stats raw-assembly.fasta > raw-assembly.stats
busco -i raw-assembly.fasta -m genome -l actinopterygii_odb10 -o $DATA_DIR/BUSCO-raw-assembly -c $THREADS

# Create snail plots through blobtoolkit
mkdir $DATA_DIR/plot-raw-assembly
blobtools create --fasta raw-assembly.fasta --busco $DATA_DIR/BUSCO-raw-assembly/run_actinopterygii_odb10/full_table.tsv $DATA_DIR/plot-raw-assembly
blobtools view --plot --view snail $DATA_DIR/plot-raw-assembly

# Scaffold assemblies based on the zebrafish reference genome
ragtag.py scaffold $DANRER_REFSEQ raw-assembly.fasta
mv ragtag_output $DATA_DIR/raw-scaffold
cp $DATA_DIR/raw-scaffold/ragtag.scaffold.fasta raw-scaffold.fasta

# Generate assembly statistics through assembly-stats and BUSCO
assembly-stats raw-scaffold.fasta > raw-scaffold.stats
busco -i raw-scaffold.fasta -m genome -l actinopterygii_odb10 -o $DATA_DIR/BUSCO-raw-scaffold -c $THREADS

# Create snail plots through blobtoolkit
mkdir $DATA_DIR/plot-raw-scaffold
blobtools create --fasta raw-scaffold.fasta --busco $DATA_DIR/BUSCO-raw-scaffold/run_actinopterygii_odb10/full_table.tsv $DATA_DIR/plot-raw-scaffold
blobtools view --plot --view snail $DATA_DIR/plot-raw-scaffold

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

# Create a RepeatModeler BLAST database
mkdir $DATA_DIR/cleanRmDb
BuildDatabase -name "$DATA_DIR/cleanRmDb/rasLatClean" $CLEAN_SCAFFOLD

# Run RepeatModeler for de novo repeat annotation
mkdir $DATA_DIR/RMClean
cd $DATA_DIR/RMClean

RepeatModeler -threads 64 -LTRStruct -database $DATA_DIR/cleanRmDb/rasLatClean

# TODO: Postprocess RM results
cat rasLat-families.fa | seqkit fx2tab | awk '{ print "rasLat"$0 }' | seqkit tab2fx > rasLat-alpha.fasta
cat rasLat-alpha.fasta | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > seqAlphaKnown.fasta
cat rasLat-alpha.fasta | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > seqAlphaUnknown.fasta
grep -c ">" $DATA_DIR/seqAlphaKnown.fasta
grep -c ">" $DATA_DIR/seqAlphaUnknown.fasta

# Run iterative repeat annotation using RepeatMasker
cd $WORKDIR
mkdir repeatAnnotation

# First annotate simple repeats
RepeatMasker -pa 16 -a -dir $DATA_DIR/seqA -noint -xsmall clean-scaffold.fasta
cp "$DATA_DIR/seqA/clean-scaffold.fasta.masked" repeatAnnotation/seqA.fasta

# round 2: annotate/mask Tetrapoda elements sourced from Repbase using output from 1st round of RepeatMasker
RepeatMasker -pa 16 -a -dir $DATA_DIR/seqB -nolow -xsmall \
    -species repeatAnnotation/seqA.fasta
cp $DATA_DIR/seqB/seqA.fasta.masked repeatAnnotation/seqB.fasta

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output from 2nd round of RepeatMasker
RepeatMasker -pa 16 -a -dir $DATA_DIR/seqC -nolow -xsmall \
    -lib $DATA_DIR/seqAlphaKnown.fasta \
    repeatAnnotation/seqB.fasta
cp $DATA_DIR/seqC/SeqB.fasta.masked repeatAnnotation/seqC.fasta

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output froom 3nd round of RepeatMasker
RepeatMasker -pa 16 -a -dir $DATA_DIR/SeqD -nolow -xsmall \
    -lib $DATA_DIR/seqAlphaUnknown.fasta \
    repeatAnnotation/seqC.fasta
cp $DATA_DIR/SeqD/seqC.fasta.masked repeatAnnotation/seqD.fasta

REPDIR=`realpath repeatAnnotation`

cd $REPDIR
OUTDIR="results"
mkdir $OUTDIR

# combine full RepeatMasker result files - .cat.gz
cat seqA/scaffold.fasta.cat.gz \
    seqB/seqA.fasta.cat.gz \
    seqC/seqB.fasta.cat.gz \
    seqD/seqC.fasta.cat.gz \
    > $OUTDIR/rasLat.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat seqA/scaffold.fasta.out \
    <(cat seqB/seqA.fasta.out | tail -n +4) \
    <(cat seqC/seqB.fasta.out | tail -n +4) \
    <(cat seqD/seqC.fasta.out | tail -n +4) \
    > $OUTDIR/rasLat.out

# copy RepeatMasker tabular files for simple repeats - .out
cat seqA/scaffold.fasta.out > $OUTDIR/simple-repeats.out

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat seqB/seqA.fasta.out \
    <(cat seqC/seqB.fasta.out | tail -n +4) \
    <(cat seqD/seqC.fasta.out | tail -n +4) \
    > $OUTDIR/complex-repeats.out

# combine RepeatMasker repeat alignments for all repeats - .align
cat seqA/scaffold.fasta.align \
    seqB/seqA.fasta.align \
    seqC/seqB.fasta.align \
    seqD/seqC.fasta.align \
    > $OUTDIR/rasLat.align

# Resummarize repeat compositions from combined analysis of all RepeatMasker rounds using RepeatMasker's ProcessRepeats
ProcessRepeats -a -species danio $OUTDIR/rasLat.cat.gz

# calculate the length of the genome sequence in the FASTA
allLen=`seqtk comp clean-scaffold.fasta | datamash sum 2`; 

# calculate the length of the N sequence in the FASTA
nLen=`seqtk comp clean-scaffold.fasta | datamash sum 9`; 

# tabulate repeats per subfamily with total bp and proportion of genome masked
cat $OUTDIR/rasLat.out | tail -n +4 | awk -v OFS="\t" '{ print $6, $7, $11 }' | 
    awk -F '[\t/]' -v OFS="\t" '{ if (NF == 3) print $3, "NA", $2 - $1 +1; else print $3, $4, $2 - $1 +1 }' | 
    datamash -sg 1,2 sum 3 | grep -v "\?" | 
    awk -v OFS="\t" -v genomeLen="${allLen}" '{ print $0, $4 / genomeLen }' > $OUTDIR/rasLat.tabulate

#  use Daren's custom script to convert .out to .gff3 for all repeats, simple repeats only, and complex repeats only
rmOutToGFF3custom -o $OUTDIR/rasLat.out > $OUTDIR/rasLatRep.gff3
rmOutToGFF3custom -o $OUTDIR/simple-repeats.out > $OUTDIR/simple-repeats.gff3
rmOutToGFF3custom -o $OUTDIR/complex-repeats.out > $OUTDIR/complex-repeats.gff3

# create simple repeat soft-masked genome
bedtools maskfasta -soft -fi scaffold.fasta -bed $OUTDIR/simple-repeats.gff3 \
    -fo $OUTDIR/rasLat-simpleRepeatMasked.fasta

# create complex repeat soft-masked genome
bedtools maskfasta -soft -fi $OUTDIR/rasLat-simpleRepeatMasked.fasta \
    -bed $OUTDIR/complex-repeats.gff3 \
    -fo $OUTDIR/rasLat-simpleComplexRepeatMasked.fasta

allLen=`seqtk comp scaffold.fasta | datamash sum 2`;
parseRM.pl -v -i $OUTDIR/rasLat.align -p -g ${allLen} -l 50,1

cp $OUTDIR/rasLat-simpleComplexRepeatMasked.fasta $WORKDIR/masked-reference.fasta
cd $WORKDIR

# Run structural annotation on the softmasked genomes
singularity exec braker3.sif braker.pl \
    --genome=masked-reference.fasta \
     --prot_seq=Vertebrata.fa \
    --softmasking \
    --workingdir="$DATA_DIR/rasLatStruct" \
    --threads $THREADS \
    --gff3 \
    --busco_lineage=actinopterygii_odb10 \
    --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_PATH

# TODO: cp BUSCO

# TODO: Runs and funannotate

# TODO: cp BUSCO and final ref and gff

