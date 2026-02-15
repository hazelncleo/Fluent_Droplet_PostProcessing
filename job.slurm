#!/bin/bash -eu
#SBATCH --job-name=EnsightPostProcess
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=90GB
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=hazel.full@monash.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=outPP.out
#SBATCH --error=errPP.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1

export ENSIGHT10_MAX_THREADS=$SLURM_NTASKS
export ENSIGHT10_MAX_CTHREADS=$SLURM_NTASKS
export CEI_EGL_GPU=0

echo "Starting Post-Processing job batch run"
uv run python -u postprocess_PALM.py "/fs04/scratch2/pa16/submodels/simple_vibration_runs/SIMPATH" "parameters.json" -o "options.json"
echo "Completed Post-Processing job successfully."
