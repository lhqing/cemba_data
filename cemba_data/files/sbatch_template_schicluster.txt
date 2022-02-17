#!/bin/bash
#
# Notes from TACC:
#
#   -- Launch this script by executing
#   -- Copy/edit this script as desired.  Launch by executing
#      "sbatch knl.openmp.slurm" on a Stampede2 login node.
#
#   -- OpenMP codes run on a single node (upper case N = 1).
#        OpenMP ignores the value of lower case n,
#        but slurm needs a plausible value to schedule the job.
#
#   -- Default value of OMP_NUM_THREADS is 1; be sure to change it!
#
#   -- Increase thread count gradually while looking for optimal setting.
#        If there is sufficient memory available, the optimal setting
#        is often 68 (1 thread per core) or 136 (2 threads per core).
#
#----------------------------------------------------

#SBATCH -J {job_name}           # Job name
#SBATCH -o {log_dir}/{job_name}.o%j       # Name of stdout output file
#SBATCH -e {log_dir}/{job_name}.e%j       # Name of stderr error file
#SBATCH -p {queue}              # Queue (partition) name
#SBATCH -N 1                    # Total # of nodes (must be 1 for OpenMP)
#SBATCH -n 1                    # Total # of mpi tasks (should be 1 for OpenMP)
#SBATCH -t {time_str}           # Run time (hh:mm:ss)
{email_str}
{email_type_str}

#----------------------------------------------------
# Clone the whole miniconda into /tmp so the snakemake command do not access $WORK
mkdir /tmp/test_{env_dir_random}

# use micromamba
export PATH=/work/05622/lhq/stampede2/bin:$PATH
micromamba shell init -s bash -p /tmp/test_{env_dir_random}
source ~/.bashrc

# activate base environment
micromamba activate

# create schicluster environment
micromamba create -y -n schicluster python=3.8 numpy scipy scikit-learn h5py \
joblib cooler pandas statsmodels rpy2 anndata xarray snakemake pybedtools htslib=1.9 pysam=0.18
micromamba activate schicluster

# export correct PYTHONPATH
export PYTHONPATH=/tmp/test_{env_dir_random}/envs/schicluster/lib/python3.8/site-packages

# install schicluster
pip install schicluster
which hicluster

# Installation finished
#----------------------------------------------------


# ---------------------------------------------------
# actual command

# print some info
date
hostname
pwd
# If you want to profile the job (CPU, MEM usage, etc.)
# load remora with
# "module load remora"
# and change the command to
# "remora {command}"


# Set thread count (default value is 1)...
export OMP_NUM_THREADS=48

for i in `seq 1 5`
do
    {command} --batch summary=${{i}}/5
done

# {command}

# delete everything in /tmp

rm -rf /tmp/test*
# ---------------------------------------------------
