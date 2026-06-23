#!/usr/bin/env python
"""
make_jsons.py
用法：
    python make_jsons.py  path/to/group.fasta   ATP
前者 = FASTA 路径，后者 = 配体 CCD 码（可省；默认 ATP）
"""
import os, sys, re, json, random
from pathlib import Path
from Bio import SeqIO

# ── 参数解析 ───────────────────────────────────────────
if len(sys.argv) < 2:
    sys.exit("Usage: python make_jsons.py group.fasta [ATP]")
fasta_path = Path(sys.argv[1]).resolve()
ligand_id  = sys.argv[2] if len(sys.argv) > 2 else "ATP"

out_dir = fasta_path.parent
print(f"→ 输出 JSON 保存目录: {out_dir}")

id_pat = re.compile(r"id\s*=\s*(\d+)")        # 捕获 id=123

# ── 逐条生成 JSON ─────────────────────────────────────
skipped, written = 0, 0
for rec in SeqIO.parse(str(fasta_path), "fasta"):
    header = rec.description
    m = id_pat.search(header)
    if not m:
        skipped += 1
        continue                              # 无 id= → 跳过原序列

    short_id  = m.group(1)                   # 仅数字
    seed_val = random.randint(0, 2**32 - 1)

    payload = {
        "name": header.split()[0],           # 长名首段
        "modelSeeds": [seed_val],            # 只用一个 seed
        "sequences": [
            {"protein": {"id": "A", "sequence": str(rec.seq)}},
            {"ligand" : {"id": "B", "ccdCodes": [ligand_id]}}
        ],
        "dialect": "alphafold3",
        "version": 3
    }

    json_file = out_dir / f"{short_id}.json"
    with open(json_file, "w") as fp:
        json.dump(payload, fp, indent=2)
    written += 1
    print(f"✓ 生成 {json_file.name}   seed={seed_val}")

print(f"完成：写入 {written} 个 JSON，跳过 {skipped} 条无 id= 序列")