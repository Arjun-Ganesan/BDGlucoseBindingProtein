#!/usr/bin/env python
# af3_rmsd_cif2pdb_lig.py

import os, sys, pathlib, pandas as pd, gemmi, numpy as np
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.protocols.simple_moves import SuperimposeMover
init("-mute all")

# ───── 用户只改这 3 行 ─────────────────────────────────
native_path = "/home/qzhong/LigandMPNN/inputs/06302025/6471_no_constraint/6471_0001.pdb"
af3_dir     = "/home/qzhong/AlphaFold3/Outputs/07032025/6471_no_constraint"
output_dir  = "/home/qzhong/AlphaFold3/Outputs/07032025/6471_no_constraint/RMSD_Calc"
# ──────────────────────────────────────────────────────

os.makedirs(output_dir, exist_ok=True)
pdb_cache_root = os.path.join(output_dir, "pdb_cache")
os.makedirs(pdb_cache_root, exist_ok=True)

# 1 ── 参考蛋白区间
native_pose = pose_from_file(native_path)
def prot_range(p):
    pro = [i for i in range(1, p.size()+1) if p.residue(i).is_protein()]
    return (pro[0], pro[-1]) if pro else (None, None)
ref_start, ref_end = prot_range(native_pose)
if None in (ref_start, ref_end):
    sys.exit("❌ 参考结构中未检测到蛋白质残基")

# 2 ── 计算 ligand-only RMSD
def calc_ligand_rmsd(p1, p2):
    xyz1, xyz2 = [], []
    for i in range(1, p1.size()+1):
        r1 = p1.residue(i)
        if r1.is_protein() or r1.is_water():
            continue                          # 只留配体
        r2 = p2.residue(i)                    # 假设序号匹配
        for ai in range(1, r1.nheavyatoms()+1):
            name = r1.atom_name(ai)
            try:
                aj = r2.atom_index(name)
            except KeyError:
                continue
            xyz1.append(r1.xyz(ai))
            xyz2.append(r2.xyz(aj))
    if not xyz1:
        return None                           # 没有配体原子
    a1 = np.array([[c.x, c.y, c.z] for c in xyz1])
    a2 = np.array([[c.x, c.y, c.z] for c in xyz2])
    return np.sqrt(((a1 - a2) ** 2).sum() / len(a1))

# 3 ── 主循环
results, errors = [], []
for root, _, files in os.walk(af3_dir):
    for fn in files:
        if not fn.endswith(".cif"):
            continue
        cif_path = os.path.join(root, fn)
        try:
            rel      = pathlib.Path(cif_path).relative_to(af3_dir)
            pdb_path = pathlib.Path(pdb_cache_root, rel).with_suffix(".pdb")
            if not pdb_path.exists():
                pdb_path.parent.mkdir(parents=True, exist_ok=True)
                st = gemmi.read_structure(cif_path)
                st.write_pdb(str(pdb_path))

            pose = pose_from_file(str(pdb_path))
            mov_start, mov_end = prot_range(pose)
            if None in (mov_start, mov_end):
                raise ValueError("当前结构中未找到蛋白质残基")

            SuperimposeMover(native_pose, ref_start, ref_end,
                             mov_start, mov_end, True).apply(pose)

            bb_rmsd  = CA_rmsd(native_pose, pose)
            lig_rmsd = calc_ligand_rmsd(native_pose, pose)

            results.append({
                "Structure"        : os.path.relpath(cif_path, af3_dir),
                "Backbone_CA_RMSD" : bb_rmsd,
                "Ligand_RMSD"      : lig_rmsd
            })

        except Exception as e:
            errors.append(f"Error processing {cif_path}: {e}")
            print(errors[-1])

pd.DataFrame(results).to_csv(os.path.join(output_dir, "rmsd_results.csv"), index=False)
if errors:
    with open(os.path.join(output_dir, "rmsd_errors.txt"), "w") as f:
        f.writelines(e + "\n" for e in errors)

print(f"\n✅ 计算完成；成功 {len(results)} 条，失败 {len(errors)} 条")