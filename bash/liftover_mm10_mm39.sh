#!/bin/bash


cat "Running liftover and bigwigtowig and wig2bed to generated mm39 tracks"
### config
OUTDIR=/mnt/bioadhoc-temp/Groups/KronenbergLab/gascui/innateness_github/data/mm39bed
BW=("GSE100738_NK.27-11b+.Sp.bw"
"GSE100738_NKT.Sp.bw"
"GSE100738_Tgd.Sp.bw"
"GSE100738_T.8.Nve.Sp.bw"
"GSE100738_T.4.Nve.Sp.bw"
)
NAME=("NK27-11b+_Sp"
"NKT_Sp"
"gdT_Sp"
"T8n_Sp"
"T4n_Sp"
)

cd /mnt/bioadhoc-temp/Groups/KronenbergLab/gascui/innateness_github/
## get overchain
#wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
cat "over chain file downloaded!"

## liftOver in .bashrc as alias. 
for i in "${!BW[@]}"
    do
    echo "${BW[$i]}" 
    # make wig file 
    bigWigToWig "data/bigwig/${BW[$i]}" "data/wig/${NAME[$i]}.wig"
    # make bed file 
    wig2bed < "data/wig/${NAME[$i]}.wig" > "data/bed/${NAME[$i]}.bed"
    # lift over 
    ~/liftOver "data/bed/${NAME[$i]}.bed" mm10ToMm39.over.chain.gz "${OUTDIR}/mm39_${NAME[$i]}.bed" "${OUTDIR}/unmapped_${NAME[$i]}.bed"
done



# make wig files
#bigWigToWig bigwig/GSE100738_NK.27-11b+.Sp.bw wig/NK.wig

# make bed files
#wig2bed < wig/NK.wig > bed/NK.bed


## liftover
#liftOver data/bed/NK.bed mm10ToMm39.over.chain.gz data/mm39/mm39_NK.bed data/mm39bed/unmapped_NK.bed


## 
## pip3 install CrossMap

