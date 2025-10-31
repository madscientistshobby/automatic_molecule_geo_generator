#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import argparse
import subprocess
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# =========================
# CONFIG (필요시 여기만 수정)
# =========================
CONFIG = {
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
}

# ==============
# 유틸 함수들
# ==============
def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def write_manifest_row(path: Path, row):
    first = not path.exists()
    with open(path, "a", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        if first:
            w.writerow(["row","status","host_smiles","ref_smiles","combined_smiles","xyz_path","debug_xyz_path","note"])
        w.writerow(row)

def write_failure_row(path: Path, row):
    first = not path.exists()
    with open(path, "a", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        if first:
            w.writerow(["row","host_smiles","reason"])
        w.writerow(row)

def xyz_to_smiles_obabel(xyz_path: Path) -> str:
    """Open Babel로 XYZ→SMILES. obabel이 PATH에 있어야 함."""
    if not xyz_path.exists():
        raise FileNotFoundError(f"ref XYZ not found: {xyz_path}")
    try:
        out = subprocess.run(
            ["obabel", str(xyz_path), "-osmi"],
            check=True, capture_output=True, text=True
        ).stdout.strip()
    except FileNotFoundError:
        raise RuntimeError("Open Babel(obabel) 미설치. conda-forge로 설치: conda install -c conda-forge openbabel")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"obabel 변환 실패: {e.stderr}")
    for line in out.splitlines():
        line = line.strip()
        if line:
            return line.split()[0]
    raise RuntimeError("obabel 출력에서 SMILES 파싱 실패.")

# ======================
# 결합 로직 (직접 AddBond)
# ======================
def find_single_dummy_and_neighbor(mol: Chem.Mol, dummy_sym: str):
    """더미 원자(예: I) 하나와 그 이웃(결합점) 반환."""
    xs = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == dummy_sym]
    if len(xs) != 1:
        raise ValueError(f"더미 {dummy_sym}는 정확히 1개여야 합니다. 현재 {len(xs)}개")
    x = xs[0]
    nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(x).GetNeighbors()]
    if len(nbrs) != 1:
        raise ValueError(f"더미 {dummy_sym}의 degree는 1이어야 합니다. 현재 {len(nbrs)}")
    return x, nbrs[0]

def remove_atom_and_sanitize(mol: Chem.Mol, rm_idx: int) -> Chem.Mol:
    rw = Chem.RWMol(mol)
    rw.BeginBatchEdit()
    rw.RemoveAtom(rm_idx)
    rw.CommitBatchEdit()
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        # 케쿨화/방향족성 등 단계적으로 처리
        Chem.SanitizeMol(out, sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                                   Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                                   Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                                   Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                                   Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION)
    return out

def combine_on_sites(mol1: Chem.Mol, mol2: Chem.Mol, site1: int, site2: int) -> Chem.Mol:
    combo = Chem.CombineMols(mol1, mol2)
    rw = Chem.RWMol(combo)
    offset = mol1.GetNumAtoms()
    rw.AddBond(site1, offset + site2, Chem.BondType.SINGLE)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out

def to_xyz_file(mol: Chem.Mol, xyz_path: Path, comment="Generated"):
    m3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
    try:
        AllChem.UFFOptimizeMolecule(m3d, maxIters=300)
    except Exception:
        pass
    conf = m3d.GetConformer()
    with open(xyz_path, "w") as f:
        f.write(f"{m3d.GetNumAtoms()}\n{comment}\n")
        for i, a in enumerate(m3d.GetAtoms()):
            p = conf.GetAtomPosition(i)
            f.write(f"{a.GetSymbol():2s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}\n")

# =====================
# per-row 처리 함수
# =====================
def process_row(idx1: int, host_smi: str, ref_smi: str, dummy_sym: str,
                out_root: Path, debug_root: Path, overwrite: bool) -> dict:
    tag = f"{idx1:03d}"
    folder = out_root / f"YS_MRTADF_{tag}"
    xyz_path = folder / f"YS_MRTADF_{tag}.xyz"
    debug_xyz_path = debug_root / f"{tag}_combined.xyz"

    if xyz_path.exists() and not overwrite:
        return dict(status="skip_exists", combined_smiles="", xyz_path=str(xyz_path),
                    debug_xyz_path=str(debug_xyz_path), note="exists")

    host = Chem.MolFromSmiles(host_smi)
    if host is None:
        raise ValueError("host SMILES 파싱 실패")
    ref = Chem.MolFromSmiles(ref_smi)
    if ref is None:
        raise ValueError("ref SMILES 파싱 실패")

    # 1) host: 더미와 결합점 찾고 더미 제거
    d_host, site_host = find_single_dummy_and_neighbor(host, dummy_sym)
    host_wo = remove_atom_and_sanitize(host, d_host)
    site_host_adj = site_host - 1 if d_host < site_host else site_host

    # 2) ref: 더미와 결합점 찾고 더미 제거
    d_ref, site_ref = find_single_dummy_and_neighbor(ref, dummy_sym)
    ref_wo = remove_atom_and_sanitize(ref, d_ref)
    site_ref_adj = site_ref - 1 if d_ref < site_ref else site_ref

    # 3) 결합
    merged = combine_on_sites(host_wo, ref_wo, site_host_adj, site_ref_adj)
    combined_smiles = Chem.MolToSmiles(merged)

    # 4) 저장
    ensure_dir(folder)
    ensure_dir(debug_root)
    to_xyz_file(merged, xyz_path, comment=f"YS_MRTADF_{tag}")
    if debug_xyz_path != xyz_path:
        to_xyz_file(merged, debug_xyz_path, comment=f"DEBUG_{tag}")

    return dict(status="ok", combined_smiles=combined_smiles,
                xyz_path=str(xyz_path), debug_xyz_path=str(debug_xyz_path), note="")

# =========
# main
# =========
def main():
    cfg = CONFIG.copy()

    ap = argparse.ArgumentParser(description="MRTADF orchestrator (dummy-direct bonding)")
    ap.add_argument("--excel", help="Excel with host SMILES")
    ap.add_argument("--smiles-col", help="Column name with host SMILES")
    ap.add_argument("--ref-smi", help="Reference SMILES (already contains dummy atom)")
    ap.add_argument("--ref-xyz", help="Reference XYZ path (already contains dummy atom)")
    ap.add_argument("--dummy", help="Dummy atom symbol (default: I)")
    ap.add_argument("--out-root", help="Output root directory")
    ap.add_argument("--debug-root", help="Debug XYZ flat directory")
    ap.add_argument("--manifest", help="Manifest TSV path")
    ap.add_argument("--failures", help="Failures TSV path")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--max-rows", type=int)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # CLI override → cfg
    for k in ["excel","smiles_col","ref_smi","ref_xyz","dummy","out_root","debug_root","manifest","failures"]:
        v = getattr(args, k) if hasattr(args, k) else None
        if v is not None:
            cfg[k] = v
    if args.overwrite: cfg["overwrite"] = True
    if args.max_rows is not None: cfg["max_rows"] = args.max_rows
    if args.dry_run: cfg["dry_run"] = True

    excel = Path(cfg["excel"])
    out_root = Path(cfg["out_root"])
    debug_root = Path(cfg["debug_root"])
    manifest_path = Path(cfg["manifest"])
    failures_path = Path(cfg["failures"])
    dummy_sym = cfg["dummy"]

    ensure_dir(out_root); ensure_dir(debug_root)
    ensure_dir(manifest_path.parent); ensure_dir(failures_path.parent)

    # ref 준비: ref_smi 우선, 없으면 ref_xyz로부터 변환
    if cfg["ref_smi"]:
        ref_smi = cfg["ref_smi"]
    else:
        if not cfg["ref_xyz"]:
            raise RuntimeError("ref 입력 필요: --ref-smi 또는 --ref-xyz 중 하나 지정")
        ref_smi = xyz_to_smiles_obabel(Path(cfg["ref_xyz"]))

    # Excel 읽기
    df = pd.read_excel(excel)
    col = cfg["smiles_col"]
    if col not in df.columns:
        raise KeyError(f"엑셀에 '{col}' 컬럼이 없습니다.")
    hosts = df[col].astype(str).tolist()
    if cfg["max_rows"]:
        hosts = hosts[:cfg["max_rows"]]

    ok = fail = skip = 0
    for i, host_smi in enumerate(hosts, start=1):
        tag = f"{i:03d}"
        try:
            if cfg["dry_run"]:
                assert Chem.MolFromSmiles(host_smi) is not None
                assert Chem.MolFromSmiles(ref_smi) is not None
                status = "dry"; combined_smiles = xyz_path = debug_xyz_path = note = ""
            else:
                res = process_row(i, host_smi, ref_smi, dummy_sym,
                                  out_root, debug_root, cfg["overwrite"])
                status = res["status"]
                combined_smiles = res.get("combined_smiles","")
                xyz_path = res.get("xyz_path","")
                debug_xyz_path = res.get("debug_xyz_path","")
                note = res.get("note","")

            if status == "ok": ok += 1
            elif status == "skip_exists": skip += 1

            write_manifest_row(manifest_path,
                               [i, status, host_smi, ref_smi, combined_smiles, xyz_path, debug_xyz_path, note])
            print(f"[{status.upper()}] Row {tag} → {xyz_path or '-'}")
        except Exception as e:
            fail += 1
            write_failure_row(failures_path, [i, host_smi, str(e)])
            print(f"[FAIL] Row {tag}: {e}")

    print(f"\nDone. OK={ok}, SKIP={skip}, FAIL={fail}")

if __name__ == "__main__":
    main()
