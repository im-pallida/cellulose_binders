#!/usr/bin/env bash
#SBATCH --job-name=RFD3_sizes
#SBATCH --partition=gpu_h200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/SCWF00146/karina/1cbh_clear/01_rfd3_symmetry_generation/logs/RFD3_sizes_%j.out
#SBATCH --error=/scratch/SCWF00146/karina/1cbh_clear/01_rfd3_symmetry_generation/logs/RFD3_sizes_%j.err
set -euo pipefail

STAGE="$(cd "$(dirname"$BASH_SOURCE[0]")../../" && pwd)"

EXPERIMENT="$1"
ORDER="$2"
SIZE="$3"
INDEX="$4"

"$STAGE/scripts/helping_scripts/15082026_02_run_one_job.sh" "$EXPERIMENT" "$ORDER" "$SIZE" "$INDEX"
