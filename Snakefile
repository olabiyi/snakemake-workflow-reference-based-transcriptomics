from os import path, getcwd
import pandas as pd

# ------- A pipeline to perform Reference based RNA-Seq analysis using bowtie and RSEM ----------- #
# Run on your local computer like so:
# snakemake pr --cores 10 --keep-going --rerun-incomplete --restart-times 3
# Rule graph
# snakemake -s Snakefile --rulegraph | dot -Tpng > rulegraph.png

configfile: "config/config.yaml"

#mode=config['mode'] # "single" or "paired"
adaptors=config['parameters']['trimmomatic']['adaptors']
min_length=config['parameters']['trimmomatic']['min_len']

sample_file = config["SAMPLE_FILE"]

CONDA_ACTIVATE = config['conda']['non_model_RNA_Seq']['env']

INDEX="02.Build_rRNA_index/{DB}.done".format(DB=config["DATABASE_NAME"])

#metadata = pd.read_table(sample_file)


# Generate the rule graph on the commadline
# Rule graph
# snakemake -s snakefile --rulegraph | dot -Tpng > rulegraph.png
# Directed Acyclic Graph (DAG)
# snakemake -s snakefile --dag | dot -Tpng > dag.png

# Run snmakemake on the cluster
#snakemake --cluster "qsub -q bioinfo.q -S /bin/bash -cwd -V -N {rule} " --jobs 100


onsuccess:
    print("Workflow completed without any error")

rule all:
    input:
        "01.Download_reference/reference.fa.gz",
        "02.Create_gene2name_file/gene2name.tsv",
        "04.QC/pre_trim/multiqc_report.html",
        "04.QC/post_trim/multiqc_report.html",
        INDEX,
        expand("05.SortSam/{sample}/{sample}.sorted.bam.bai", sample=config["SAMPLES"]),
        expand(["06.remove_rRNA/{sample}/{sample}_R1.fastq.gz", "06.remove_rRNA/{sample}/{sample}_R2.fastq.gz"],
                sample=config['SAMPLES']),
        "07.QC/unmapped_reads/multiqc_report.html",
        expand("09.Plot_lengths/{sample}/{sample}_diagnostic.pdf", sample=config['SAMPLES']),
        "10.Generate_count_matrix/gene_counts_matrix.tsv"

rule Download_reference:
    output:
        fasta = "01.Download_reference/reference.fa.gz",
        gtf =   "01.Download_reference/reference.gtf.gz"
    log: "logs/Download_reference/Download_reference.log"
    params: 
        fasta = config['FASTA'],
        gtf = config['GTF'],
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 10
#    conda: config['CONDA']
    shell:
        """
         {params.CONDA_ACTIVATE}
         # Download the reference gene fasta sequences and GTF abnnotation
         wget -O {output.fasta} {params.fasta}
         wget -O {output.gtf} {params.gtf}
        """


# Create a file that maps gene ids to names to be used when making heatmaps
# for easy of results interpretation
rule Create_gene2name_file:
    input: rules.Download_reference.output.gtf
    output: "02.Create_gene2name_file/gene2name.tsv"
    log: "logs/Create_gene2name_file/Create_gene2name_file.log"
    threads: 10
#    conda: config['CONDA']
    params:
        CONDA_ACTIVATE=CONDA_ACTIVATE
    shell:
        """
        {params.CONDA_ACTIVATE}
        bioawk -c gff  \
           '$feature=="gene"{{printf "%s\\t%s\\n", $seqname, $group}}' \
            {input} > {output}
        """

REF="03.Prepare_reference/ref/{organism}_ref".format(organism=config['ORGANISM'])
rule Prepare_reference:
    input:
        gtf=rules.Download_reference.output.gtf,
        fasta=rules.Download_reference.output.fasta
    output:
        multiext(REF,".grp",".ti", ".transcripts.fa", ".seq",
                 ".chrlist", ".idx.fa", ".n2g.idx.fa")
    log: "logs/Prepare_reference/Prepare_reference.log"
    params:
        aligner=config['ALIGNER'], # "bowtie", "bowtie2", "star",
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 20
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
        ALIGNER={params.aligner}

        [ -e GENOME.fasta ] || zcat {input.fasta} > GENOME.fasta
        [ -e GENOME.GTF ] || zcat {input.gtf} > GENOME.GTF

        if [ ${{ALIGNER}} == "bowtie" ]; then
          rsem-prepare-reference \
              --gtf GENOME.GTF \
              --bowtie \
              GENOME.fasta  {REF}

        elif [ ${{ALIGNER}} == "star" ]; then

          rsem-prepare-reference \
              --gtf GENOME.GTF \
              --star \
              GENOME.fasta {REF}

        else

          rsem-prepare-reference \
              --gtf GENOME.GTF \
              --bowtie2 \
              GENOME.fasta {REF}

        fi

        rm -rf GENOME.GTF
        rm -rf GENOME.fasta
        """

rule QC_pre_trim:
    input:
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz"
    output:
        forward_html="04.QC/pre_trim/{sample}/{sample}_R1_fastqc.html",
        rev_html="04.QC/pre_trim/{sample}/{sample}_R2_fastqc.html"
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=10,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    log: "logs/QC_pre_trim/{sample}/{sample}.log"
    threads: 10
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
          fastqc --outdir {params.out_dir}/ \
             --threads {params.threads} \
             {input.forward} {input.rev} > {log} 2>&1
        """


rule SummarizeQC_pre_trim:
    input:
        forward_html=expand("04.QC/pre_trim/{sample}/{sample}_R1_fastqc.html", sample=config['SAMPLES']),
        rev_html=expand("04.QC/pre_trim/{sample}/{sample}_R2_fastqc.html", sample=config['SAMPLES'])
    output:
        "04.QC/pre_trim/multiqc_report.html"
    log: "logs/SummarizeQC_pre_trim/multiqc.log"
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 10
#    conda: config['CONDA']
    shell:
        """
         {params.CONDA_ACTIVATE}
          multiqc \
              --interactive \
              -f {params.out_dir} \
              -o {params.out_dir} > {log} 2>&1
        """

rule Trim_reads:
    input:
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz"
    output:
        r1="05.Trim_reads/{sample}/{sample}_R1.fastq.gz",
        r2="05.Trim_reads/{sample}/{sample}_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="05.Trim_reads/{sample}/{sample}_R1.unpaired.fastq.gz",
        r2_unpaired="05.Trim_reads/{sample}/{sample}_R2.unpaired.fastq.gz"
    log: "logs/Trim_reads/{sample}/{sample}.log"
    params:
        trimmer="ILLUMINACLIP:{adaptors}:2:30:10"
                " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20"
                " MINLEN:{min_length}".format(adaptors=adaptors,
                                          min_length=min_length),
        threads=20,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 20
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
        trimmomatic PE \
        -threads {params.threads} \
        {input.forward} {input.rev} \
        {output.r1} {output.r1_unpaired} \
        {output.r2} {output.r2_unpaired} \
        {params.trimmer} > {log} 2>&1
        """

# Automatically gues and trime adapters using trim galore
rule Trim_galore_trim_adaptors:
    input: 
        forward=rules.Trim_reads.output.r1,
        rev=rules.Trim_reads.output.r2
    output:
        forward_reads="05.Trim_galore_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        rev_reads="05.Trim_galore_trim_adaptors/{sample}/{sample}_R2.fastq.gz",
    log: "logs/Trim_galore_trim_adaptors/{sample}/{sample}.log"
    threads: 1
#    conda: config['CONDA']
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        CONDA_ACTIVATE=CONDA_ACTIVATE
    shell:
        """ 
        {params.CONDA_ACTIVATE}
        trim_galore \
           -o {params.out_dir} \
           --paired {input.forward} {input.rev} > {log} 2>&1

        #Rename the files
        #Fastq files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1.fq.gz {params.out_dir}/{wildcards.sample}_R1.fastq.gz
        mv {params.out_dir}/{wildcards.sample}_R2_val_2.fq.gz {params.out_dir}/{wildcards.sample}_R2.fastq.gz
        """

rule QC_post_trim:
    input:
        forward=rules.Trim_galore_trim_adaptors.output.forward_reads,
        rev=rules.Trim_galore_trim_adaptors.output.rev_reads  
    output:
        forward_html="04.QC/post_trim/{sample}/{sample}_R1_fastqc.html",
        rev_html="04.QC/post_trim/{sample}/{sample}_R2_fastqc.html"
    log: "logs/QC_post_trim/{sample}/{sample.log}"
    threads: 10
#    conda: config['CONDA']
    params:
        threads=10,
        out_dir=lambda w, output: path.dirname(output[0]),
        CONDA_ACTIVATE=CONDA_ACTIVATE
    shell:
        """
        {params.CONDA_ACTIVATE}
        fastqc --outdir {params.out_dir} \
            --threads {params.threads} \
            {input.forward} {input.rev} > {log}  2>&1
        """

rule SummarizeQC_post_trim:
    input:
        forward_html=expand("04.QC/post_trim/{sample}/{sample}_R1_fastqc.html",
                 sample=config['SAMPLES']),
        rev_html=expand("04.QC/post_trim/{sample}/{sample}_R2_fastqc.html",
                 sample=config['SAMPLES'])
    output:
        "04.QC/post_trim/multiqc_report.html"
    log: "logs/SummarizeQC_post_trim/multiqc.log" 
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 10
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
        multiqc --interactive -f {params.out_dir} \
            -o {params.out_dir} > {log} 2>&1
        """

rule Download_rRNA:
    output:  "01.Download_rRNA/silva.fasta"
    threads: 10
    params:
        LSU=config["LSU_URL"],
        SSU=config["SSU_URL"],
        out_dir=lambda wildcards, output: path.dirname(output[0]),
        CONDA_ACTIVATE=CONDA_ACTIVATE
    log: "logs/Download_rRNA/Download_rRNA.log"
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
        wget -O {params.out_dir}/LSU.fasta.gz {params.LSU} > {log} 2>&1
        wget -O {params.out_dir}/SSU.fasta.gz {params.SSU} >> {log} 2>&1

        zcat {params.out_dir}/LSU.fasta.gz \
                {params.out_dir}/SSU.fasta.gz > {params.out_dir}/temp.rna.fasta
        bioawk -c fastx '{{print ">"$name" "$comment; gsub(/U/,"T",$seq); print $seq}}' \
            {params.out_dir}/temp.rna.fasta  >  {output}
        """

rule Build_rRNA_index:
    input: rules.Download_rRNA.output
    output:
        index="02.Build_rRNA_index/{DB}.1.bt2".format(DB=config["DATABASE_NAME"]),
        done="02.Build_rRNA_index/{DB}.done".format(DB=config["DATABASE_NAME"])
    log: "logs/Build_rRNA_index/Build_rRNA_index.log"
    threads: 20
#    conda: config['CONDA']
    params:
        index="02.Build_rRNA_index/{DB}".format(DB=config["DATABASE_NAME"]),
        threads=20,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    shell:
        """
        {params.CONDA_ACTIVATE}
        bowtie2-build \
            --threads {params.threads} \
                {input} {params.index} > {log} 2>&1 && \
        touch {output.done}
        """

rule Map_reads2rRNA:
    input:
        index_done=rules.Build_rRNA_index.output.done,
        forward=rules.Trim_galore_trim_adaptors.output.forward_reads,
        rev=rules.Trim_galore_trim_adaptors.output.rev_reads
    output:
        stats="04.Map_reads2rRNA/{sample}/{sample}.stats",
        sam="04.Map_reads2rRNA/{sample}/{sample}.sam"
    log: "logs/Map_reads2rRNA/{sample}/{sample}.log"
#    conda: config['CONDA']
    params:
        index="02.Build_rRNA_index/{DB}".format(DB=config["DATABASE_NAME"]),
        threads=20,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 20
    shell:
        """
        {params.CONDA_ACTIVATE}
        bowtie2 \
            -p {params.threads} \
            --rg-id {wildcards.sample} \
            --rg SM:{wildcards.sample} \
            -x {params.index} \
            -1 {input.forward} \
            -2 {input.rev} \
            --met-file {output.stats} \
            -S {output.sam} > {log} 2>&1
        """

rule SortSam:
    input: rules.Map_reads2rRNA.output.sam
    output:
        sorted_bam="05.SortSam/{sample}/{sample}.sorted.bam",
        index="05.SortSam/{sample}/{sample}.sorted.bam.bai"
    log: "logs/SortSam/{sample}/{sample}.log"
#    conda: config['CONDA']
    params:
        threads=20,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 20
    shell:
        """
        {params.CONDA_ACTIVATE}
            [ -e {output.sorted_bam} ] || \
            samtools sort -@ {params.threads} \
            -o {output.sorted_bam} {input} > {log} 2>&1

            samtools index {output.sorted_bam} {output.index} >> {log} 2>&1
        """


rule Remove_rRNA:
    input: rules.SortSam.output.sorted_bam
    output:
        flag_stat="06.remove_rRNA/{sample}/{sample}.flagstat.txt",
        stats="06.remove_rRNA/{sample}/{sample}.stats.txt",
        idxstats="06.remove_rRNA/{sample}/{sample}.idxstats.txt",
        forward="06.remove_rRNA/{sample}/{sample}_R1.fastq.gz",
        rev="06.remove_rRNA/{sample}/{sample}_R2.fastq.gz",
        single="06.remove_rRNA/{sample}/{sample}.S.fastq.gz"
    log: "logs/Remove_rRNA/{sample}/{sample}.log"
#    conda: config['CONDA']
    params:
        flags="-f 12 -t -F 256",
        forward="06.remove_rRNA/{sample}/{sample}_R1.fastq",
        rev="06.remove_rRNA/{sample}/{sample}_R2.fastq",
        threads=20,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 20
    shell:
        """
        {params.CONDA_ACTIVATE}
        samtools flagstat {input} > {output.flag_stat} 2> {log}
        samtools stats --remove-dups {input} > {output.stats} 2> {log}
        samtools idxstats {input} > {output.idxstats} 2> {log}
        samtools fastq  {params.flags} -1 {params.forward} -2 {params.rev} \
                -0 {output.single}  {input} >> {log} 2>&1

        parallel -j 2 'gzip {{}}' ::: {params.forward} {params.rev}

        """

rule QC_unmapped_reads:
    input:
        forward=rules.Remove_rRNA.output.forward,
        rev=rules.Remove_rRNA.output.rev
    output:
        forward_html="07.QC/unmapped_reads/{sample}/{sample}_R1_fastqc.html",
        reverse_html="07.QC/unmapped_reads/{sample}/{sample}_R2_fastqc.html"
#    conda: config['CONDA']
    params:
        out_dir= lambda w, output: path.dirname(output[0]),
        threads=10,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    log: "logs/QC_unmapped_reads/{sample}/{sample}.log"
    threads: 10
    shell:
        """
        {params.CONDA_ACTIVATE}
        fastqc --outdir {params.out_dir} \
               --threads {threads} \
               {input.forward} {input.rev} > {log} 2>&1
        """


rule SummarizeQC_unmapped:
    input:
        expand(["07.QC/unmapped_reads/{sample}/{sample}_R1_fastqc.html",
                "07.QC/unmapped_reads/{sample}/{sample}_R2_fastqc.html"],
               sample=config['SAMPLES'])
    output:
        "07.QC/unmapped_reads/multiqc_report.html"
    threads: 10
    log: "logs/SummarizeQC_unmapped_reads/multiqc.log"
#    conda: config['CONDA']
    params:
        program = config['programs_path']['multiqc'],
        out_dir = lambda w, output: path.dirname(output[0]),
        align_stat_dir="06.remove_rRNA/",
        CONDA_ACTIVATE=CONDA_ACTIVATE
    shell:
        """
        {params.CONDA_ACTIVATE}
        multiqc --interactive \
             -f {params.out_dir} {params.align_stat_dir} \
             -o {params.out_dir} > {log} 2>&1
         """


rule Estimate_abundance:
    input: 
        rules.Prepare_reference.output,
        forward= rules.Remove_rRNA.output.forward,
        rev=rules.Remove_rRNA.output.rev,
    output: 
        "08.Estimate_abundance/{sample}/{sample}.isoforms.results",
        "08.Estimate_abundance/{sample}/{sample}.genes.results",
        #"08.Estimate_abundance/{sample}/{sample}.genome.sorted.bam",
        #"08.Estimate_abundance/{sample}/{sample}.genome.sorted.bam.bai"
    log: "logs/Estimate_abundance/{sample}/{sample}.log"
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        aligner=config['ALIGNER'], # "bowtie", "bowtie2", "star"
        threads=20,
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 20
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
        ALIGNER={params.aligner}
        
	# Unzip the fastq files
        [ -e {params.out_dir}/{wildcards.sample}_R1.fastq ] || \
          zcat {input.forward} > {params.out_dir}/{wildcards.sample}_R1.fastq
        [ -e {params.out_dir}/{wildcards.sample}_R2.fastq ] || \
          zcat {input.rev} > {params.out_dir}/{wildcards.sample}_R2.fastq

        if [ ${{ALIGNER}} == "bowtie" ]; then

          rsem-calculate-expression \
              -p {params.threads} \
              --bowtie \
              --estimate-rspd \
              --append-names \
              --output-genome-bam \
              --paired-end \
              {params.out_dir}/{wildcards.sample}_R1.fastq \
              {params.out_dir}/{wildcards.sample}_R2.fastq \
              {REF} {params.out_dir}/{wildcards.sample} > {log} 2>&1

        elif [ ${{ALIGNER}} == "star" ]; then

          rsem-calculate-expression \
              -p {params.threads} \
              --star \
              --estimate-rspd \
              --append-names \
              --output-genome-bam \
              --paired-end \
              {params.out_dir}/{wildcards.sample}_R1.fastq \
              {params.out_dir}/{wildcards.sample}_R2.fastq \
              {REF} {params.out_dir}/{wildcards.sample} > {log} 2>&1

        else

          rsem-calculate-expression \
              -p {params.threads} \
              --bowtie2 \
              --estimate-rspd \
              --append-names \
              --output-genome-bam \
              --paired-end \
              {params.out_dir}/{wildcards.sample}_R1.fastq \
              {params.out_dir}/{wildcards.sample}_R2.fastq \
              {REF} {params.out_dir}/{wildcards.sample} > {log} 2>&1

        fi
   
        #clean
        rm -rf \
            {params.out_dir}/{wildcards.sample}_R1.fastq \
            {params.out_dir}/{wildcards.sample}_R2.fastq
        """

# Generate Transcripts length plot
rule Plot_lengths:
    input:
        rules.Estimate_abundance.output
    output:
        "09.Plot_lengths/{sample}/{sample}_diagnostic.pdf"
    log: "logs/Plot_lengths/{sample}/{sample}.log"
    params:
        in_dir=lambda w, input: path.dirname(input[0]),
        CONDA_ACTIVATE=CONDA_ACTIVATE
    threads: 10
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
        rsem-plot-model {params.in_dir}/{wildcards.sample} {output} > {log} 2>&1
        """

# Generate genes count matrix
rule Generate_count_matrix:
    input:
        expand("08.Estimate_abundance/{sample}/{sample}.genes.results",
               sample=config['SAMPLES'])
    output:
        "10.Generate_count_matrix/gene_counts_matrix.tsv"
    log: "logs/Generate_count_matrix/count_matrix.log"
    threads: 10
    params:
        CONDA_ACTIVATE=CONDA_ACTIVATE
#    conda: config['CONDA']
    shell:
        """
        {params.CONDA_ACTIVATE}
        rsem-generate-data-matrix {input} > {output} 2> {log}
        """
