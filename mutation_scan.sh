#!/bin/bash
#BATCH --job-name=process_goodread_job       # Job name
#SBATCH --output=process_goodread_%j.log      # Output log file (%j will be replaced with the job ID)
#SBATCH --error=process_goodread_%j.err       # Error log file
#SBATCH --time=48:00:00                  # Time limit (hh:mm:ss)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=50               # Number of CPU cores per task
#SBATCH --mem=600GB                      # Memory pool for all cores
#SBATCH --partition=highcpu              # Partition (queue) name
# Set TMPDIR
mkdir -p /scratch/$USER/$SLURM_JOB_ID
export TMPDIR=/scratch/$USER/$SLURM_JOB_ID
#SBATCH --tmp=400G
module load python/3.8.3
echo "Python version: $(python --version)"
echo "Running script with TMPDIR=$TMPDIR"
# Run your Python script
python mutational_analysis_fixed.py || echo "Python script failed with exit code $?"

# Clean up TMPDIR
rm -rf /scratch/$USER/$SLURM_JOB_ID
