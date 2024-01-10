#!usr/bin/env bash

alias bigBedToBed='/mnt/hpc-apps/UCSC/UCSC-2023/bigBedToBed'




bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/mm39/jaspar/JASPAR2022.bb -chrom=chr11 -start=97005200 -end=97006801 output/data/klrb1c.bed