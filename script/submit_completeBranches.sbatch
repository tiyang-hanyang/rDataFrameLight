#!/bin/bash
#SBATCH -J compBr
#SBATCH -p high_cpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --open-mode=append
#SBATCH --time=02:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-999%50
#SBATCH --output=/home/tiyang/public/br_complete/logs/%x_%A_%a.out
#SBATCH --error=/home/tiyang/public/br_complete/logs/%x_%A_%a.err

export X509_USER_PROXY=$HOME/proxy/x509up

echo "Using proxy at: $X509_USER_PROXY"
/opt/ohpc/pub/voms/voms-proxy-info --all


line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "/home/tiyang/public/br_complete/jobs.txt")
skim=$(echo "$line" | awk '{print $1}')
raw=$(echo "$line" | awk '{print $2}')
out=$(echo "$line" | awk '{print $3}')
brj=$(echo "$line" | awk '{print $4}')

completeBranches_oneFile --skimmed "$skim" --rawFile "$raw" --branches "$brj" --outFile "$out"

