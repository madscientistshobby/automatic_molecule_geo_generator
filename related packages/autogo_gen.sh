#!/bin/bash
set -euo pipefail

REF="YS_GO.in"
[[ -f "$REF" ]] || { echo "레퍼런스 파일 없음: $REF"; exit 1; }

# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
find . -type f -name "*.xyz" -print0 | while IFS= read -r -d '' in_file; do
  dir=$(dirname "$in_file")
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO.in"            # 새 파일 이름

  if [[ -f "$xyz_path" ]]; then
    cp -f "$REF" "$target"
    sed -i "s/filename/${base}.xyz/g" "$target"
    echo "[OK] $target 생성 및 치환 완료"
  else
    echo "[SKIP] $xyz_path 없음, 건너뜀"
  fi
done
