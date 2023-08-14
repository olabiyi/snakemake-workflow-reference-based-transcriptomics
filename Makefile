.PHONY: run qc create_env

usage:
	@echo "make QC # Quality check and filter reads"
	@echo "make all # Run the complete pipeline"
	@echo "make create_env # Create the conda environment with tools for RNA-seq"

run: qc
	sbatch scripts/02.submit-slurm.sh
qc: 
	sbatch scripts/01.QC_submit-slurm.sh

create_env:
	conda env create -f envs/non_model_RNA_Seq_conda.yaml
