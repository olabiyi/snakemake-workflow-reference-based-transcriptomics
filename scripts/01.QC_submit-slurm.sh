#!/bin/bash
#SBATCH --job-name=super_job       #Set the job name to "JobExample2"
#SBATCH --time=23:50:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks=1          #Request 1 tasks/cores per node
#SBATCH --mem=1G                     #Request 1GB per node 
#SBATCH --output=super-script.o.%j      #Send stdout/err to "Example2Out.[jobID]" 
#SBATCH --error=super-script.e.%j    #Send std err to "Example2error.[jobID]"

#module purge; module load iccifort/2020.1.217 impi/2019.7.217 snakemake/5.26.1-Python-3.8.2
module purge; module load GCC/11.3.0  OpenMPI/4.1.4 snakemake/7.22.0
#module load WebProxy Anaconda3/2020.07
# Generate the rule graph on the commadline
# Rule graph
# snakemake -s Snakefile --rulegraph | dot -Tpng > rulegraph.png
# Directed Acyclic Graph (DAG)
# snakemake -s Snakefile --dag | dot -Tpng > dag.png

# Run snmakemake on the cluster
# --jobs 100 # submit a maximum 100 jobs
# --latency-wait 60 # wait for 60 seconds before declaring that a job has failed
# --skip-script-cleanup 

# run the line below first on the command line to generate the log directories for each rule
# before submitting the jobs to the cluster
#snakemake -npr  make_logs_directories --cores 1

snakemake -pr  -s Snakefile \
        --jobs 20 \
        --keep-going \
        --rerun-incomplete \
        --keep-incomplete \
        --latency-wait 60 \
        --cluster-config config/cluster.yaml \
        --local-cores 7 \
        --cluster "sbatch --partition {cluster.queue} --mem={cluster.mem} --time={cluster.time} --ntasks={cluster.threads}"\
        "04.QC/post_trim/multiqc_report.html"

