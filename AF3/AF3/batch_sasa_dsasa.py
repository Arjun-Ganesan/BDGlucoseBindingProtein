#!/usr/bin/env python
# batch_sasa_dsasa.py ---------------------------------------------------
# 用法:  python batch_sasa_dsasa.py  [pdb_root_dir]
# 如果未提供目录则默认当前目录

import sys
import argparse
from pathlib import Path

import pandas as pd
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta import core, protocols

# ---- 配置：根据自己数据集修改 ---------------------------------------
PROT_CHAINS = "A"          # 蛋白链，可写 "AB"
LIG_CHAINS  = "B"          # 配体链，可写 "C"
PROBE_RADIUS = 1.4         # 1.4 Å 水探针

# ---- 解析命令行参数 --------------------------------------------------
parser = argparse.ArgumentParser(
    description="Batch compute protein SASA and interface dSASA for every .pdb file found (recursively) under a directory.")
parser.add_argument("pdb_root_dir", nargs="?", default=".",
                    help="Root directory containing PDB files (default: current directory)")
args = parser.parse_args()

pdb_root = Path(args.pdb_root_dir).resolve()
out_csv  = pdb_root / "sasa_summary.csv"
err_log  = pdb_root / "errors.txt"

# ---- 初始化 PyRosetta -----------------------------------------------
init("-mute all")          # 静默模式

from pyrosetta.rosetta.core.scoring import sasa as _sasa

_sasa_calc = _sasa.SasaCalc()
_sasa_calc.set_probe_radius(PROBE_RADIUS)

def calc_sasa_and_dsasa(pdb_file: Path):
    """返回 (protein_SASA, interface_dSASA)"""
    pose = pose_from_pdb(str(pdb_file))

    # 选链
    sel_prot = core.select.residue_selector.ChainSelector(PROT_CHAINS)
    prot_mask = sel_prot.apply(pose)

    # ---- total protein SASA ----
    _sasa_calc.calculate(pose)
    sasa_vec = _sasa_calc.get_residue_sasa()  # vector1<double>
    protein_sasa = sum(sasa_vec[i] for i in range(1, pose.size()+1)
                       if prot_mask[i])

    # ---- interface dSASA ----
    intf = f"{PROT_CHAINS}_{LIG_CHAINS}"
    iam = protocols.analysis.InterfaceAnalyzerMover(intf)
    iam.apply(pose)
    iam.add_score_info_to_pose(pose)
    interface_dsasa = pose.scores.get("dSASA_int", float("nan"))

    return protein_sasa, interface_dsasa


# ---- 主循环 -----------------------------------------------------------
results, errors = [], []

pdb_files = list(pdb_root.rglob("*.pdb"))
if not pdb_files:
    sys.exit(f"No .pdb files found under {pdb_root}")

for pdb_path in pdb_files:
    try:
        prot_sasa, d_sasa = calc_sasa_and_dsasa(pdb_path)
        results.append({
            "pdb_file"      : str(pdb_path.relative_to(pdb_root)),
            "protein_sasa"  : prot_sasa,
            "interface_dsasa": d_sasa
        })
        print(f"✓ {pdb_path.name:35s}  protSASA={prot_sasa:8.1f}  dSASA={d_sasa:8.1f}")
    except Exception as e:
        msg = f"{pdb_path}: {e}"
        errors.append(msg)
        print("✗", msg)

# ---- 写 CSV / 错误日志 -----------------------------------------------
pd.DataFrame(results).to_csv(out_csv, index=False)
if errors:
    err_log.write_text("\n".join(errors))

print(f"\n完成：{len(results)} 条成功，{len(errors)} 条失败")
print(f"结果已写入 {out_csv}")
if errors:
    print(f"错误日志写入 {err_log}")