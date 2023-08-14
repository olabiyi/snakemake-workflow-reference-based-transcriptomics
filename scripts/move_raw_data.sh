#!/usr/bin/env bash 

# A script to make a seperate directory for each sample and then move the sample's files there

SAMPLES=("Control_A" "Control_B" "Control_C" "Control_D" \
         "Test_1A" "Test_1B" "Test_1C" \
         "Test_10A" "Test_10B" "Test_10C" "Test_10D")


for sample in ${SAMPLES[*]}; do
    [ -d 01.raw_data/${sample} ] || mkdir 01.raw_data/${sample} && \
    find 01.raw_data/ -type f -name "${sample}*" | xargs -I {} mv {} 01.raw_data/${sample}/ 
done
