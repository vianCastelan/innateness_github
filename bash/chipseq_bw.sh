#!/usr/bin/env bash

## get common UCSC apps (new server)
alias bedGraphToBigWig='/mnt/hpc-apps/UCSC/UCSC-2023/bedGraphToBigWig'
alias fetchChromSizes='/mnt/hpc-apps/UCSC/UCSC-2023/fetchChromSizes'
## fetch chromosome sizes
fetchChromSizes mm10 > "data/mm10.chrom.sizes"
fetchChromSizes mm9 > "data/mm9.chrom.sizes"
fetchChromSizes mm39 > "data/mm39.chrom.sizes"
## parameters
READDIR="/mnt/bioadhoc-temp/Groups/KronenbergLab/gascui/innateness_github"
cd $READDIR
# GSE101593 CD4+ T cells SMAD4
# GSE55834 H3K4me3, H3K27me3, H3K4me3 NK cells (this is on mm9)
# GSE50129 NK, CD8+ T cell RUNX3
# GSE135533 CD8+ T cells SMAD4
# GSE58775 CD4+ T cell H3K4me3 and H3K27me3
# GSE75724 CD4+ T cells Blimp1

#===============#  datasets with bed.gz files #===============# 
#data/GSE135533
#data/GSE101593 ## this data is flawed....
#data/GSE50129 ## mm9 only gives 0s ... use bed instead here? 
## datasets with .bed.gz files 
directories=()
# Find directories containing .bed.gz files and add them to the array
while IFS= read -r -d '' dir; do
    if [[ ! " ${directories[@]} " =~ " ${dir} " ]]; then
        directories+=("$dir")
    fi
done < <(find data -type f -name "*.bed.gz" -exec dirname {} \; | sort -u -z)
# Remove duplicate directories using awk
unique_directories=($(printf "%s\n" "${directories[@]}" | awk '!seen[$0]++'))

## Generate bedgraph and Bigwig files 
for dir in "${unique_directories[@]}"; do
    DATASET=$(basename "$dir")
    readarray -t FILES <<< "$(find $DATADIR -type f -name "*.bed.gz")"
    for file in "${FILES[@]}"; do
        echo "${file%.bed.gz}"
        bedfile="${file%.bed.gz}.bedgraph"
        zcat "$file" | awk '{printf "%s\t%d\t%d\t%0.3f\n" , $1,$2,$3,$5}' | sort -k1,1 -k2,2n > "$bedfile"
        if [ $DATASET == 'GSE50129' ]; then 
            ## use mm9 
            bedGraphToBigWig $bedfile "data/mm9.chrom.sizes" "${bedfile%.bedgraph}.bw"
            ## 
            outfile="${file%.bed.gz}_mm10"
            ## run crossmap to convert to mm10
            CrossMap.py wig mm9ToMm10.over.chain.gz $bedfile $outfile
        else
            ## convert to BigWig
            bedGraphToBigWig $bedfile "data/mm10.chrom.sizes" "${bedfile%.bedgraph}.bw"
        fi
    done
    echo "$DATASET all done"
done

#===============#  datasets with .bedgraph.gz files #===============# 
#data/GSE135533 SMAD4 CD*
#data/GSE58775
#data/GSE101593
#data/GSE75724
#data/GSE50129
#data/GSE55834 mm9
directories=()
# Find directories containing .bed.gz files and add them to the array
while IFS= read -r -d '' dir; do
    if [[ ! " ${directories[@]} " =~ " ${dir} " ]]; then
        directories+=("$dir")
    fi
done < <(find data -type f \( -iname "*.bed.gz" -o -iname "*.bedGraph.gz" \) -exec dirname {} \; | sort -u -z)
# Remove duplicate directories using awk
unique_directories=($(printf "%s\n" "${directories[@]}" | awk '!seen[$0]++'))

## Generate  Bigwig files 
for dir in "${unique_directories[@]}"; do
    DATASET=$(basename "$dir")
    readarray -t FILES <<< "$(find $DATADIR -type f -name "*.bedgraph.gz")"
    for file in "${FILES[@]}"; do
        echo "${file%.bedgraph.gz}"
        bedfile="${file%.bed.gz}.bedgraph"
        zcat "$file" | sort -k1,1 -k2,2n > "$bedfile"
        if [ $DATASET == 'GSE55834' ]; then 
            ## use mm9 
            bedGraphToBigWig $bedfile "data/mm9.chrom.sizes" "${bedfile%.bedgraph}.bw"
            outfile="${file%.bedgraph.gz}_mm10"
            ## run crossmap to convert to mm10
            CrossMap.py wig mm9ToMm10.over.chain.gz $bedfile $outfile
        else
            ## convert to BigWig
            bedGraphToBigWig $bedfile "data/mm10.chrom.sizes" "${bedfile%.bedgraph}.bw"
        fi
    done
    echo "$DATASET all done"
done

