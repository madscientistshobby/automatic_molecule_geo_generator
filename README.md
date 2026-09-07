# automatic_molecule_geo_generator

Joins each molecule in an Excel SMILES list to a reference fragment at marked attachment sites and writes starting 3D geometries as XYZ files.

## The problem

Preparing a substituted-molecule screen by hand means joining each candidate to the same fragment, generating coordinates, and saving a separate calculation folder. At an illustrative five minutes per molecule, 50 candidates take about four hours; this is a workflow estimate, not a measured benchmark.

## What it does

`allinone.py` reads the first sheet of an Excel workbook (`.xlsx`) and a reference SMILES or `.xyz` file. It removes one marked atom from each fragment, joins their neighbors with a single bond, adds hydrogens, generates an RDKit ETKDG conformer, and attempts up to 300 UFF optimization iterations. Coordinates are written in angstroms.

| CONFIG field | Value shipped in the script / meaning |
| --- | --- |
| `excel` | `/home/baikgrp/calcs/YS/MRTADF/Code/allinone/Substance_filtered_single_Cl.xlsx` |
| `smiles_col` | `smiles`; Excel header containing host SMILES |
| `ref_smi` | `None`; a nonempty reference SMILES takes precedence over `ref_xyz` |
| `ref_xyz` | `/home/baikgrp/calcs/YS/MRTADF/Code/allinone/ref_truncated.xyz`; converted to SMILES by `obabel` |
| `dummy` | `Cl`; attachment marker on both fragments |
| `out_root` | `out`; writes `YS_MRTADF_001/YS_MRTADF_001.xyz`, etc. |
| `debug_root` | `_debug`; writes `001_combined.xyz`, etc. |
| `manifest` | `out/manifest.tsv`; appends successful, skipped, and dry-run rows |
| `failures` | `logs/failures.tsv`; appends row numbers, host SMILES, and failure reasons |
| `overwrite` | `False`; skip an existing main XYZ |
| `max_rows` | `None`; optional limit on input rows |
| `dry_run` | `False`; parse SMILES without joining or generating coordinates |

Standalone helpers in `related packages/` also generate ORCA `.in` files from XYZ files, submit selected jobs, check `*GO.out` completion messages, and plot Excel data as PNG files. They are not called by `allinone.py`:

- `autogo_gen.sh` reads `REF="YS_GO.in"`; the SP scripts use `YS_singlet_SP1.in`, `YS_triplet_SP1.in`, `YS_singlet_SP2.in`, and `YS_triplet_SP2.in` and search for `*GO.xyz`.
- `autosubmit_allsubfolder.sh` selects `*GO.in`; `autosubmit_SPSOL.sh` selects `*SPSOL.in`. Both call `/home/baikgrp/bin/jbatch_orca600` through `srun`.
- `orca_chk.py` uses `root="."`, `target_line="ORCA TERMINATED NORMALLY"`, and `target_line2="OPTIMIZATION RUN DONE"`.
- Both regression scripts use `file_path='filtered_smiles_dataset_ph_x.xlsx'`, `targets=['Est', 'E_singlet', 'E_triplet']`, and `sigma_para_exp` as x data (named `x_col` in the nonlinear script). Their `output_dir` values are `regression_plots` and `regression_plots_multi`. The latter fits linear, quadratic, cubic, exponential, logarithmic, and sigmoid curves and ranks them by training-data R².

## Requirements

- Core script: Python, pandas, RDKit, and openpyxl for `.xlsx` input. Open Babel's `obabel` executable must be on `PATH` when using `ref_xyz`.
- `environment.yaml` records Python 3.10.19, pandas 2.3.3, RDKit 2025.09.1, openpyxl 3.1.5, and Open Babel 3.1.1 in conda environment `ys_rdkit`. It is a Linux build-specific export with prefix `/home/baikgrp/calcs/KL/miniconda3/envs/crest/ys_rdkit` and includes packages the scripts do not use.
- Regression helpers additionally import NumPy, Matplotlib, and SciPy. SciPy is missing from the supplied environment file.
- Shell helpers need Bash and GNU-style `sed -i`. Submission needs Slurm, ORCA, and the lab's external `jbatch_orca600` wrapper, which is not included.

## Quickstart

Supply your own workbook and reference XYZ; neither is included. The block below uses the real CLI overrides for `CONFIG["excel"]`, `CONFIG["smiles_col"]`, `CONFIG["ref_xyz"]`, and `CONFIG["dummy"]`, so the lab paths in the source can stay as they are. Edit the four shell values before running. Each fragment must have exactly one terminal `Cl` marker for this setup.

```bash
git clone https://github.com/madscientistshobby/automatic_molecule_geo_generator.git
cd automatic_molecule_geo_generator
conda create -n ys_rdkit -c conda-forge python=3.10 pandas rdkit openpyxl openbabel -y
conda activate ys_rdkit

excel="/absolute/path/to/your.xlsx"
smiles_col="smiles"
ref_xyz="/absolute/path/to/your_reference.xyz"
dummy="Cl"
python allinone.py --excel "$excel" --smiles-col "$smiles_col" \
  --ref-xyz "$ref_xyz" --dummy "$dummy"
```

On a compatible Linux system, `conda env create -f environment.yaml -n ys_rdkit` is an alternative to the `conda create` command above. For a reference SMILES, replace `--ref-xyz "$ref_xyz"` with `--ref-smi 'CCCl'` (illustrative fragment). Other overrides are `--out-root`, `--debug-root`, `--manifest`, `--failures`, `--overwrite`, `--max-rows`, and `--dry-run`. Install `numpy matplotlib scipy` with conda before using the optional regression scripts.

## Example

<!-- TODO: add example input + output screenshot -->

## Notes / Limitations

- `dummy` matches an element symbol, not an atom-map label. Each fragment must contain exactly one such atom with exactly one neighbor; another real chlorine atom will also count. The actual default is `Cl`, despite help text mentioning `I`.
- XYZ bond orders are inferred by Open Babel. Reference coordinates and atom ordering are not retained: the product is rebuilt, with host atoms followed by reference atoms and added hydrogens. Check bonding and stereochemistry before calculation.
- These are starting geometries, not a conformer search or a DFT optimization. Embedding return codes and UFF convergence are not checked; UFF exceptions are ignored. Main and debug XYZ files are embedded separately and can differ. No random seed is set.
- Excel needs an exact column header on its first sheet. Blank cells become strings and can fail parsing. Output numbering starts at the first data row as `001`; filenames always use `YS_MRTADF_`. Reordering input can make existing-file skips misleading. Logs append on reruns, and row failures do not make the final process exit nonzero.
- `--dry-run` checks SMILES parsing only, not attachment sites or geometry, and still creates directories and writes logs. Changing `out_root` does not relocate `manifest` automatically.
- Shell helpers search from the working directory and need their templates there; generated `.in` files are overwritten. ORCA templates contain fixed charge, multiplicity, methods, and memory, plus repeated/conflicting `%pal` settings. In particular, `YS_triplet_SP1.in` specifies `0 1`. Review these as-is before use. The checker counts missing optimization text as failure even if normal termination is present.
- Regression inputs need numeric data; energy units are not specified or converted. The linear script does not drop missing values. The nonlinear script silently skips failed fits; R² rankings are not validation of a physical model.

## 한국어

Excel의 SMILES 목록과 참조 조각을 더미 원자 위치에서 연결하고, 수소를 추가한 초기 3D XYZ 구조를 저장합니다. 위 Quickstart에서 `excel`, `smiles_col`, `ref_xyz`, `dummy`를 실제 입력에 맞게 바꾼 뒤 실행하세요.

<details>
<summary>기존 README 원문</summary>

```text
# automatic_molecule_geo_generator

결국 ref.xyz 와 스크리닝하고 싶은 분자들의 SMILES 가 있는 Excel 을 제공하면 알아서 Geometry.xyz 를 만들어줌!
CONFIG 에서 필요한것들, 경로, 엑셀명, Ref.xyz, dummy 원자 를 입력하고 실행하면 끝

이후는 관련 .sh 나 위 .py 들 사용해서 자동화 가능 (related packages 폴더 안에 내용물들 사용하면 됨)
```

</details>

<details>
<summary>기존 한국어 주석 원문 (파일 순서대로 발췌)</summary>

`allinone.py`

```text
# CONFIG (필요시 여기만 수정)
    "excel": "/home/baikgrp/calcs/YS/MRTADF/Code/allinone/Substance_filtered_single_Cl.xlsx",        # 호스트 SMILES가 있는 엑셀 경로
    "smiles_col": "smiles",            # 엑셀 컬럼명
    "ref_smi": None,                   # 참조 SMILES 문자열(더미 포함). 없으면 None
    "ref_xyz": "/home/baikgrp/calcs/YS/MRTADF/Code/allinone/ref_truncated.xyz",         # 참조 XYZ 경로(더미 포함). ref_smi가 None이면 이걸 사용
    "dummy": "Cl",                      # 더미 원자 심볼 (기본 I) — 양쪽 모두 동일해야 함
    "out_root": "out",                 # 개별 xyz 폴더 루트
    "debug_root": "_debug",            # 디버그 xyz 모아두는 폴더(하위폴더 없음)
    "manifest": "out/manifest.tsv",    # 성공/스킵 기록
    "failures": "logs/failures.tsv",   # 실패 기록
    "overwrite": False,                # 기존 파일 있으면 덮어쓸지 여부
    "max_rows": None,                  # 테스트 시 상한 (예: 5)
    "dry_run": False,                  # 드라이런: 구조 검증만
# 유틸 함수들
    """Open Babel로 XYZ→SMILES. obabel이 PATH에 있어야 함."""
# 결합 로직 (직접 AddBond)
    """더미 원자(예: I) 하나와 그 이웃(결합점) 반환."""
        # 케쿨화/방향족성 등 단계적으로 처리
# per-row 처리 함수
    # 1) host: 더미와 결합점 찾고 더미 제거
    # 2) ref: 더미와 결합점 찾고 더미 제거
    # 3) 결합
    # 4) 저장
    # ref 준비: ref_smi 우선, 없으면 ref_xyz로부터 변환
    # Excel 읽기
```

`related packages/autogo_gen.sh`

```text
# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO.in"            # 새 파일 이름
```

`related packages/autosp_gen_sp1.sh`

```text
# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO_singlet_SP1.in"            # 새 파일 이름
# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO_triplet_SP1.in"            # 새 파일 이름
```

`related packages/autosp_gen_sp2.sh`

```text
# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO_singlet_SP2.in"            # 새 파일 이름
# 모든 하위폴더의 xyz 파일 탐색 (예: something.xyz)
  base=$(basename "$in_file" .xyz)          # filename 추출
  xyz_path="$dir/${base}.xyz"              # filename.xyz 경로
  target="$dir/${base}GO_triplet_SP2.in"            # 새 파일 이름
```

`related packages/autosubmit_SPSOL.sh`

```text
# 모든 하위폴더의 .in 파일 탐색
```

`related packages/autosubmit_allsubfolder.sh`

```text
# 모든 하위폴더의 GO.in 파일 탐색
```

`related packages/excel_reg.py`

```text
# 엑셀 파일 불러오기
file_path = 'filtered_smiles_dataset_ph_x.xlsx'  # 파일 이름 수정
# 결과 저장 폴더 생성
# 분석 대상 열 리스트
    # 선형회귀
    # 회귀 결과 출력
    # 그래프 생성
    # 자동 저장
```

`related packages/excel_reg_not_linear.py`

```text
# ========= 사용자 설정 =========
# 데이터 로드
# ----- 모델 정의 -----
    # a * ln(x) + b (x>0에서만)
# (이름, 함수, 도메인 조건)
# ----- 각 target에 대해 피팅 -----
        # 도메인 조건 적용 (예: log)
            # 해당 모델의 유효 구간에서만 곡선 생성
    # R² 기준 정렬
    # ✅ 모든 모델에 대해 개별 플롯 저장
        # 원 데이터 전체 산점도
        # 해당 모델 피팅 곡선
```

`related packages/orca_chk.py`

```text
        if filename.endswith("GO.out"):  # <-- 여기서 GO.out만 찾음
```

</details>
