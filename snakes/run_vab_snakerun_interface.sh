#!/bin/bash 

results_dir="/mnt/scratch/scratch"
target_dir="virtual_aging_cohort_interface"
out_dir="OUT_VIRT_COH_IFACE"
err_dir="ERR_VIRT_COH_IFACE"

echo  "Create results folder"

echo  ${results_dir}/${target_dir}
mkdir ${results_dir}/${target_dir}

echo  "Create I/O snakemake folder"

mkdir ${results_dir}/${out_dir}
mkdir ${results_dir}/${err_dir}

echo ${results_dir}/${err_dir}
echo ${results_dir}/${err_dir}

echo "Activate THE ENVIRONMENT"

. env/bin/activate

echo  "SNAKEMAKE"
echo  "Direct RUN"

# If you want to execute a dry-run, uncomment the following line and comment the --cluster line
# snakemake -s snakes/run_vab_snakemake.snake  -n

snakemake -s snakes/run_vab_snakemake.snake --cluster "sbatch --time 48:00:00 -J vab_SBIF --ntasks-per-socket=4 --mem-per-cpu=32G --output ${results_dir}/${out_dir}/virtual_aging%j.log --error ${results_dir}/${err_dir}/virtual_aging%j.log" -j 2000
