#!/usr/bin/env bash

# A script to download the reference genome

# Here we download the human genome fro Esemble
FASTA="http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz"
GTF="http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz"

[ -d  01.Download_reference/ ] || mkdir  01.Download_reference/

wget -O 01.Download_reference/reference.fa.gz  ${FASTA}

wget -O 01.Download_reference/reference.gtf  ${GTF}
