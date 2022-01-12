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

metadata = pd.read_table(sample_file)

# Generate the rule graph on the commadline
# Rule graph
# snakemake -s snakefile --rulegraph | dot -Tpng > rulegraph.png
# Directed Acyclic Graph (DAG)
# snakemake -s snakefile --dag | dot -Tpng > dag.png

# Run snmakemake on the cluster
#snakemake --cluster "qsub -q bioinfo.q -S /bin/bash -cwd -V -N {rule} " --jobs 100


onsuccess:
    print("Workflow completed without any error")


RULES = ["Download_reference","Rename_files", "Create_gene2name_file", "QC_pre_trim",
         "SummarizeQC_pre_trim", "Prepare_reference", "Trim_reads", "QC_post_trim",
         "SummarizeQC_post_trim", "Estimate_abundance", "Get_mapping_statistics",
         "Summarize_mapping_statistics", "Plot_lengths", "Generate_count_matrix"]


rule all:
    input:
        expand("01.raw_data/{sample}/{sample}.fastq.gz", sample=config['SAMPLES']),
        "01.Download_reference/reference.fa.gz",
        "02.Create_gene2name_file/gene2name.tsv",
        "04.QC/pre_trim/multiqc_report.html",
        "04.QC/post_trim/multiqc_report.html",
        #"08.Summarize_mapping_statistics/multiqc_report.html",
        expand("09.Plot_lengths/{sample}/{sample}_diagnostic.pdf", sample=config['SAMPLES']),
        "10.Generate_count_matrix/gene_counts_matrix.tsv"




# This rule will Make rule specific log directories
# # in order to easily store the standard input and standard error
# # generated when submiting jobs to the cluster
rule Make_logs_directories:
    output:
        directory("logs/Download_reference/"),
        directory("logs/Rename_files/"),
        directory("logs/Create_gene2name_file/"),
        directory("logs/QC_pre_trim/"),
        directory("logs/SummarizeQC_pre_trim/"),
        directory("logs/Prepare_reference/"),
        directory("logs/Estimate_abundance/"),
        directory("logs/Plot_lengths/"),
        directory("logs/Generate_count_matrix/")
    threads: 1
    log: "logs/Make_logs_directories/Make_logs_directories.log"
    shell:
        """
         [ -d logs/ ] || mkdir -p logs/
         cd logs/
         for RULE in {RULES}; do
          [ -d ${{RULE}}/ ] || mkdir -p ${{RULE}}/
         done
        """


rule Download_reference:
    input:
        log_dirs = rules.Make_logs_directories.output
    output:
        fasta = "01.Download_reference/reference.fa.gz",
        gtf =   "01.Download_reference/reference.gtf"
    log: "logs/Download_reference/Download_reference.log"
    params: 
        fasta = config['FASTA'],
        gtf = config['GTF']
    threads: 5
    shell:
        """
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
    threads: 1
    shell:
        """
        function geneid2name(){{

	     # A function to generate a table of geneid to name from a GTF file
	     # USAGE: geneid2name <path_to_gtf_file>
	     local GTF=$1

         grep -v "#" ${{GTF}} | \
	     awk -F'\t' '{{print $9}}' | \
	     awk -F';' 'BEGIN{{OFS="\t"}} {{print $1,$3}}' | \
	     grep "gene_name" | \
	     sed -E 's/gene_id \"(.+?)\" gene_name \"(.+?)\"/\1\t\2/g'

        }}

         geneid2name {input} > {output}
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
        program=config['programs_path']['rsem']['prepare_reference'],
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        aligner=config['ALIGNER'] # "bowtie", "bowtie2", "star"
    threads: 12
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        ALIGNER={params.aligner}

        [ -e GENOME.fasta ] || zcat {input.fasta} > GENOME.fasta
        [ -e GENOME.GTF ] || zcat {input.gtf} > GENOME.GTF

        if [ ${{ALIGNER}} == "bowtie" ]; then
          {params.program} \
              --gtf GENOME.GTF \
              --bowtie \
              GENOME.fasta  {REF}

        elif [ ${{ALIGNER}} == "star" ]; then

          {params.program} \
              --gtf GENOME.GTF \
              --star \
              GENOME.fasta {REF}

        else

          {params.program} \
              --gtf GENOME.GTF \
              --bowtie2 \
              GENOME.fasta {REF}

        fi

        rm -rf GENOME.GTF
        rm -rf GENOME.fasta
        """

rule Rename_files:
#    input: 
#        log_dirs=rules.Make_logs_directories.output
    output:
        expand("01.raw_data/{sample}/{sample}.fastq.gz", sample=config['SAMPLES'])
    log: "logs/Rename_files/Rename_files.log"
    threads: 5
    run:
        for old,new in zip(metadata.Old_name,metadata.New_name):
            shell("[ -f {new} ] || mv {old} {new}".format(old=old, new=new))



rule QC_pre_trim:
    input:
        read="01.raw_data/{sample}/{sample}.fastq.gz",
#        log_dirs=rules.Make_logs_directories.output
    output:
        "04.QC/pre_trim/{sample}/{sample}_fastqc.html"
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        threads=5
    log: "logs/QC_pre_trim/{sample}/{sample}.log"
    threads: 5
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

          {params.program} --outdir {params.out_dir}/ \
             --threads {params.threads} {input.read}

        """


rule SummarizeQC_pre_trim:
    input:
        expand("04.QC/pre_trim/{sample}/{sample}_fastqc.html", sample=config['SAMPLES'])
    output:
        "04.QC/pre_trim/multiqc_report.html"
    log: "logs/SummarizeQC_pre_trim/multiqc.log"
    params:
        program=config['programs_path']['multiqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

          {params.program} \
              --interactive \
              -f {params.out_dir} \
              -o {params.out_dir}
        """


rule Trim_reads:
    input:
        read="01.raw_data/{sample}/{sample}.fastq.gz",
#        log_dirs=rules.Make_logs_directories.output
    output:
        "05.Trim_reads/{sample}/{sample}.fastq.gz"
    log:
        "logs/Trim_reads/{sample}/{sample}.log"
    params:
        program=config['programs_path']['trimmomatic'],
        trimmer="ILLUMINACLIP:{adaptors}:2:30:10"
                " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20"
                " MINLEN:{min_length}".format(adaptors=adaptors,
                                          min_length=min_length),
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib']
    threads: 5
    resources:
        mem_mb=1024
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} SE -threads {threads} {input.read} {output} {params.trimmer} >  {log} 2>&1

        """

rule QC_post_trim:
    input:
        "05.Trim_reads/{sample}/{sample}.fastq.gz"
    output:
        "04.QC/post_trim/{sample}/{sample}_fastqc.html"
    threads: 1
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib']
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} --outdir {params.out_dir} \
        --threads {threads} {input} 
        """

rule SummarizeQC_post_trim:
    input:
        expand("04.QC/post_trim/{sample}/{sample}_fastqc.html",
                 sample=config['SAMPLES'])
    output:
        "04.QC/post_trim/multiqc_report.html"
    params:
        program=config['programs_path']['multiqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} --interactive -f {params.out_dir} -o {params.out_dir}
        """



rule Estimate_abundance:
    input: 
        rules.Prepare_reference.output,
        forward="05.Trim_reads/{sample}/{sample}.fastq.gz"
    output: 
        "06.Estimate_abundance/{sample}/{sample}.isoforms.results",
        "06.Estimate_abundance/{sample}/{sample}.genes.results",
        #"06.Estimate_abundance/{sample}/{sample}.genome.sorted.bam",
        #"06.Estimate_abundance/{sample}/{sample}.genome.sorted.bam.bai"
    log: "logs/Estimate_abundance/{sample}.log"
    params:
        program=config['programs_path']['rsem']['calculate_expression'],
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        out_dir=lambda w, output: path.dirname(output[0]),
        aligner=config['ALIGNER'] # "bowtie", "bowtie2", "star"
    threads: 10
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        ALIGNER={params.aligner}
        
        [ -e {params.out_dir}/{wildcards.sample}.fasta ] || zcat {input.forward} > {params.out_dir}/{wildcards.sample}.fasta

        if [ ${{ALIGNER}} == "bowtie" ]; then

          {params.program} \
              -p {threads} \
              --bowtie \
              --estimate-rspd \
              --append-names \
              --output-genome-bam \
              {params.out_dir}/{wildcards.sample}.fasta \
              {REF} {params.out_dir}/{wildcards.sample}

        elif [ ${{ALIGNER}} == "star" ]; then

          {params.program} \
              -p {threads} \
              --star \
              --estimate-rspd \
              --append-names \
              --output-genome-bam \
              {params.out_dir}/{wildcards.sample}.fasta \
              {REF} {params.out_dir}/{wildcards.sample}


        else

          {params.program} \
              -p {threads} \
              --bowtie2 \
              --estimate-rspd \
              --append-names \
              --output-genome-bam \
              {params.out_dir}/{wildcards.sample}.fasta \
              {REF} {params.out_dir}/{wildcards.sample}


        fi
   
        #clean
        rm -rf {params.out_dir}/{wildcards.sample}.fasta
        """



#rule Get_mapping_statistics:
#    input:
#        sorted_bam="06.Estimate_abundance/{sample}/{sample}.genome.sorted.bam"
#    output:
#        flag_stat="07.Get_mapping_statistics/{sample}/{sample}.flagstat.txt",
#        stats="07.Get_mapping_statistics/{sample}/{sample}.stats.txt",
#        idxstats="07.Get_mapping_statistics/{sample}/{sample}.idxstats.txt"
#    log: "logs/Get_mapping_statistics/{sample}.log"
#    params:
#        program=config['programs_path']['samtools'],
#        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
#        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib']
#    threads: 2
#    shell:
#        """
#        set +u
#        {params.conda_activate}
#        {params.PERL5LIB}
#        set -u
#        {params.program} flagstat {input.sorted_bam} > {output.flag_stat}
#        {params.program} stats --remove-dups {input.sorted_bam} > {output.stats}
#        {params.program} idxstats {input.sorted_bam} > {output.idxstats}
#       """
#
#
#rule Summarize_mapping_statistics:
#    input: 
#        expand("07.Get_mapping_statistics/{sample}/{sample}.idxstats.txt",
#               sample=config['SAMPLES'])
#    output:
#        "08.Summarize_mapping_statistics/multiqc_report.html"
#    log: "logs/Summarize_mapping_statistics/Summarize_mapping_statistics.log"
#    threads: 2
#    params:
#        program = config['programs_path']['multiqc'],
#        stats_dir = lambda w, input: path.dirname(input[0]).split('/')[0],
#        out_dir = lambda w, output: path.dirname(output[0])
#    shell:
#        "{params.program} --interactive -f {params.stats_dir}  -o {params.out_dir}"
#


# Generate Transcripts length plot
rule Plot_lengths:
    input:
        "06.Estimate_abundance/{sample}/{sample}.isoforms.results",
        "06.Estimate_abundance/{sample}/{sample}.genes.results"
    output:
        "09.Plot_lengths/{sample}/{sample}_diagnostic.pdf"
    params:
        program=config['programs_path']['rsem']['plot_model'],
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        in_dir=lambda w, input: path.dirname(input[0])
    threads: 2
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} {params.in_dir}/{wildcards.sample} {output}
        """

# Generate genes count matrix
rule Generate_count_matrix:
    input:
        expand("06.Estimate_abundance/{sample}/{sample}.genes.results",
               sample=config['SAMPLES'])
    output:
        "10.Generate_count_matrix/gene_counts_matrix.tsv"
    params:
        program=config['programs_path']['rsem']['generate_matrix'],
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib']
    threads: 10
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} {input} > {output}
        """
