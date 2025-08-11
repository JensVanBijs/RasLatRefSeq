#! /bin/bash

WORKDIR=`realpath ./`
RAW_DATA_DIR=`realpath ./0-1-raw_data`
DATA_DIR=`realpath ./0-2-data`

RAW_REF="$DATA_DIR/"  # TODO
CORRECTED_REF="$DATA_DIR/"  # TODO
CLEANED_REF="$DATA_DIR/"  # TODO

cd $WORKDIR

# Run structural annotation on the softmasked genomes
# TODO: not galba??
singularity exec braker3.sif braker.pl --genome=Raslat_raw-scaffold_mask_sort.fasta --prot_seq=Vertebrata.fa --softmasking --workingdir=Raslat_raw1 --threads 16 --gff3 --busco_lineage=actinopterygii_odb10 --AUGUSTUS_CONFIG_PATH=/home/nurhidayatl/.augustus/ 

# TODO: BUSCO

# TODO: Runs and funannotate

# TODO: BUSCO
