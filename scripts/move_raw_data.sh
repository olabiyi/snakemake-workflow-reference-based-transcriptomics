#!/usr/bin/env bash 

# A script to make a seperate directory for each sample and then move the sample's files there

SAMPLES=("0A" "0B" "0C" "0D" "25A" "25B" "25C" "25D")


for sample in ${SAMPLES[*]}; do
    [ -d 01.raw_data/${sample} ] || mkdir 01.raw_data/${sample} && \
    find 01.raw_data/ -type f -name "${sample}*" | xargs -I {} mv {} 01.raw_data/${sample}/ 
done
