#!/bin/bash
set -euo pipefail

#singlet
REF="YS_singlet_SP2.in"
[[ -f "$REF" ]] || { echo "레퍼런스 파일 없음: $REF"; exit 1; }

# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
find . -type f -name "*GO.xyz" -print0 | while IFS= read -r -d '' in_file; do
  dir=$(dirname "$in_file")
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO_singlet_SP2.in"            # 새 파일 이름

  if [[ -f "$xyz_path" ]]; then
    cp -f "$REF" "$target"
    sed -i "s/filename/${base}.xyz/g" "$target"
    echo "[OK] $target 생성 및 치환 완료"
  else
    echo "[SKIP] $xyz_path 없음, 건너뜀"
  fi
done

#Triplet
REF="YS_triplet_SP2.in"
[[ -f "$REF" ]] || { echo "레퍼런스 파일 없음: $REF"; exit 1; }

# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
find . -type f -name "*GO.xyz" -print0 | while IFS= read -r -d '' in_file; do
  dir=$(dirname "$in_file")
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO_triplet_SP2.in"            # 새 파일 이름

  if [[ -f "$xyz_path" ]]; then
    cp -f "$REF" "$target"
    sed -i "s/filename/${base}.xyz/g" "$target"
    echo "[OK] $target 생성 및 치환 완료"
  else
    echo "[SKIP] $xyz_path 없음, 건너뜀"
  fi
done
