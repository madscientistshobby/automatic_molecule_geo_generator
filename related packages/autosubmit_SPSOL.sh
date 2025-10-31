#!/bin/bash -l
#SBATCH -J autosubmit
#SBATCH -c 1
#SBATCH --time=00:30:00
set -euo pipefail
shopt -s nullglob

# 모든 하위폴더의 .in 파일 탐색
find . -type f -name "*SPSOL.in" | while read -r f; do
  dir=$(dirname "$f")
  base=$(basename "$f" .in)
  echo "[RUN] srun jsub6 $dir/$base"
  (
    cd "$dir"
    printf '8\no\n' | srun --ntasks=1 --cpus-per-task=1 /home/baikgrp/bin/jbatch_orca600 "$base"
  )
done
