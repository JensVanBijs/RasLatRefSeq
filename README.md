# *Rasbora lateristriata* Draft Reference Genome Creation and Computational Lineage Estimation

> **NOTE**: This repository contains the code that was used for creating the annotated draft reference genome for the Yogyakartan *Rasbora lateristriata* fish. Code was (mostly) written by Jens van Bijsterveld as part of the master's thesis computer science (bioinformatics) at Leiden University, which was conducted in collaboration with the Universitas Gadjah Mada (UGM) faculty of biology.
<!-- TODO: Check UGM and faculty spelling. -->

## Usage note

All utilized code is currently present in the repository, but is still being processed into a (mostly) automated continuous pipeline. Currently, running the pipeline might take some additional user input. The most up-to-date version of the pipeline can be found under `pipeline.sh`.

## Thesis abstract

*Rasbora lateristriata* is a fish species native to the Indonesian Yogyakarta region that is theorized to be an improved model organism for researching afflictions that require a higher body temperature than the *Danio rerio* model organism can reasonably survive. In
order to develop the *Rasbora lateristriata* as a model organism, a reference genome was assembled from Oxford Nanopore Technologies reads using flye, and annotated using InterProScan, EggNOG, and Pfam. Multiple assemblies were fully processed through the assembly and annotation pipeline to reach the eventually selected reference genome. A synteny analysis and a phylogenetic analysis were then performed on the reference genome, resulting in further evidence that the *Rasbora lateristriata* indeed belongs to the *Danionidae* family. The creation of the reference genome and confirmation of the *Rasbora lateristriata* lineage enables future work in this new model organism like gene knockout studies, the creation of transgenic lines, and methylation analysis.
