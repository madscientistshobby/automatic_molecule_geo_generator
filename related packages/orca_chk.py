#!/usr/bin/env python3
import os

root = "."
target_line = "ORCA TERMINATED NORMALLY"
target_line2 = "OPTIMIZATION RUN DONE"
ok_count = 0
op_count = 0
fail_count = 0
total_checked = 0

for dirpath, dirnames, filenames in os.walk(root):
    for filename in filenames:
        if filename.endswith("GO.out"):  # <-- 여기서 GO.out만 찾음
            total_checked += 1
            out_file = os.path.join(dirpath, filename)
            with open(out_file, "r", errors="ignore") as f:
                content = f.read()
                if target_line in content:
                    print(f"[OK]   {out_file}")
                    ok_count += 1
                if target_line2 in content:
                    print(f"[OK]   {out_file}")
                    op_count += 1

                else:
                    print(f"[FAIL] {out_file}")
                    fail_count += 1

print("\n=== 요약 ===")
print(f"총 검사한 파일 수: {total_checked}")
print(f"정상 종료 (ORCA TERMINATED NORMALLY): {ok_count}")
print(f"Opt 성공 (OPTIMIZATION RUN DONE): {op_count}")
print(f"비정상 종료: {fail_count}")
