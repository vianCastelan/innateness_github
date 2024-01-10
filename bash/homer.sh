#!/usr/bin/env bash

perl /home/gascui/local/homer/.//configureHomer.pl -install mm10

PEAK=/mnt/BioAdHoc/Groups/KronenbergLab/gascui/innateness_github/new_index.bed
ANNPEAK=/mnt/BioAdHoc/Groups/KronenbergLab/gascui/innateness_github/data/countmatrix/ImmGEN_ATACseq_annotated_peaks.txt
annotatePeaks.pl $PEAK mm10  > $ANNPEAK