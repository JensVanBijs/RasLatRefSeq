#! /bin/bash

WORKDIR=`realpath ./`
RAW_DATA_DIR=`realpath ./0-1-raw_data`
DATA_DIR=`realpath ./0-2-data`
THREADS=64

RAW_SCAFFOLD="$WORKDIR/raw-scaffold.fasta"
CORRECTED_SCAFFOLD="$WORKDIR/corrected-scaffold.fasta"
CLEAN_SCAFFOLD="$WORKDIR/clean-scaffold.fasta"

cd $WORKDIR

# Create a RepeatModeler BLAST databases
mkdir $DATA_DIR/rawRmDb
BuildDatabase -name "$DATA_DIR/rawRmDb/rasLatRaw" $RAW_SCAFFOLD

mkdir $DATA_DIR/corrRmDb
BuildDatabase -name "$DATA_DIR/corrRmDb/rasLatCorr" $CORRECTED_SCAFFOLD

mkdir $DATA_DIR/cleanRmDb
BuildDatabase -name "$DATA_DIR/cleanRmDb/rasLatClean" $CLEAN_SCAFFOLD

# Run RepeatModeler for de novo repeat annotation
mkdir $DATA_DIR/RMRaw
mkdir $DATA_DIR/RMCorr
mkdir $DATA_DIR/RMClean

cd $DATA_DIR/RMRaw
RepeatModeler -threads 64 -LTRStruct -database $DATA_DIR/rawRmDb/rasLatRaw

cd $DATA_DIR/RMCorr
RepeatModeler -threads 64 -LTRStruct -database $DATA_DIR/corrRmDb/rasLatCorr

cd $DATA_DIR/RMClean
RepeatModeler -threads 64 -LTRStruct -database $DATA_DIR/cleanRmDb/rasLatClean

cd $WORKDIR

# TODO: Postprocess RM results
cat rasLat-families.fa | seqkit fx2tab | awk '{ print "rasLat"$0 }' | seqkit tab2fx > rasLat-alpha.fasta

cat rasLat-alpha.fasta | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > seqAlphaKnown.fasta
cat rasLat-alpha.fasta | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > seqAlphaUnknown.fasta

grep -c ">" seqAlphaKnown.fasta
grep -c ">" seqAlphaUnknown.fasta

# Run iterative repeat annotation using RepeatMasker
# TODO: For all scaffolds! For loop
# for draftAssembly


# First annotate simple repeats
RepeatMasker -pa 16 -a -dir seqA -noint -xsmall $draftAssembly
cp "seqA/${draftAssembly}.fasta.masked" ./seqA.fasta

# round 2: annotate/mask Tetrapoda elements sourced from Repbase using output from 1st round of RepeatMasker
RepeatMasker -pa 16 -a -dir seqB -nolow -xsmall \
    -species seqA.fasta
cp seqB/seqA.fasta.masked ./seqB.fasta

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 03_known_out -nolow \
-lib round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.known \
02_tetrapoda_out/reference-genome.tetrapoda_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log
# round 3: rename outputs
rename tetrapoda_mask.masked.fasta known_mask 03_known_out/reference-genome*
rename .masked .masked.fasta 03_known_out/reference-genome*

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output froom 3nd round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 04_unknown_out -nolow \
-lib round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.unknown \
03_known_out/reference-genome.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log
# round 4: rename outputs
rename known_mask.masked.fasta unknown_mask 04_unknown_out/reference-genome*
rename .masked .masked.fasta 04_unknown_out/reference-genome*

WORKDIR="./RepMasker"

cd $WORKDIR
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

# TODO: Add processRepeats to PATH
# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species danio $OUTDIR/rasLat.cat.gz

# calculate the length of the genome sequence in the FASTA
allLen=`seqtk comp reference-genome.fasta | datamash sum 2`; 
# calculate the length of the N sequence in the FASTA
nLen=`seqtk comp reference-genome.fasta | datamash sum 9`; 

# tabulate repeats per subfamily with total bp and proportion of genome masked
cat $OUTDIR/rasLat.out | tail -n +4 | awk -v OFS="\t" '{ print $6, $7, $11 }' | 
    awk -F '[\t/]' -v OFS="\t" '{ if (NF == 3) print $3, "NA", $2 - $1 +1; else print $3, $4, $2 - $1 +1 }' | 
    datamash -sg 1,2 sum 3 | grep -v "\?" | 
    awk -v OFS="\t" -v genomeLen="${allLen}" '{ print $0, $4 / genomeLen }' > $OUTDIR/rasLat.tabulate

#  use Daren's custom script to convert .out to .gff3 for all repeats, simple repeats only, and complex repeats only
rmOutToGFF3custom -o $OUTDIR/rasLat.out > $OUTDIR/rasLat.gff3
rmOutToGFF3custom -o $OUTDIR/simple-repeats.out > $OUTDIR/simple-repeats.gff3
rmOutToGFF3custom -o $OUTDIR/complex-repeats.out > $OUTDIR/complex-repeats.gff3

# create masked genome FASTA files

# create simple repeat soft-masked genome
bedtools maskfasta -soft -fi scaffold.fasta -bed $OUTDIR/simple-repeats.gff3 \
    -fo $OUTDIR/rasLat-simpleRepeatMasked.fasta

# create complex repeat soft-masked genome
bedtools maskfasta -soft -fi $OUTDIR/rasLat-simpleRepeatMasked.fasta \
    -bed $OUTDIR/complex-repeats.gff3 \
    -fo $OUTDIR/rasLat-simpleComplexRepeatMasked.fasta

allLen=`seqtk comp scaffold.fasta | datamash sum 2`;
parseRM.pl -v -i $OUTDIR/rasLat.align -p -g ${allLen} -l 50,1
