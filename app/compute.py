"""
Protein Structure Analysis Engine

研究室のノートブック(DSA_Cis_250317.ipynb)のロジックを完全移植
- cis結合検出
- Cα距離行列計算
- ヒートマップ生成
"""

from __future__ import annotations
import os
import json
import math
from datetime import datetime
from typing import List, Dict, Any, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # GUIなし環境
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import MMCIFParser, PPBuilder, is_aa


# ─── 保存先を絶対パスに固定 ─────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]  # app/compute.py -> protein-service/
ART_ROOT = BASE_DIR / "data" / "artifacts"
ART_ROOT.mkdir(parents=True, exist_ok=True)


def _ensure_dir(p: Path):
    """ディレクトリ作成"""
    p.mkdir(parents=True, exist_ok=True)


# ─── mmCIF読み込み ─────────────────────────────────────
def _load_structures(cif_paths: List[str]):
    """
    mmCIFファイルを読み込んで構造オブジェクトのリストを返す
    """
    parser = MMCIFParser(QUIET=True)
    structs = []
    for p in cif_paths:
        try:
            s = parser.get_structure(Path(p).stem, p)
            structs.append(s)
        except Exception as e:
            print(f"Warning: Failed to parse {p}: {e}")
            continue
    return structs


def _collect_residues(struct) -> List:
    """
    構造から標準アミノ酸残基のリストを取得（1モデルのみ）
    """
    residues = []
    for model in struct:
        for chain in model:
            for res in chain:
                if is_aa(res, standard=True):
                    residues.append(res)
        break  # 1モデルのみ
    return residues


def _resname(res) -> str:
    """残基名を取得"""
    try:
        return res.get_resname().strip()
    except Exception:
        return "UNK"


def _atom(res, name):
    """原子を取得"""
    try:
        return res[name]
    except Exception:
        return None


def _dihedral(a, b, c, d) -> float | None:
    """
    4原子の二面角を計算（度数法）
    ノートブックのロジックを移植
    """
    if a is None or b is None or c is None or d is None:
        return None
    try:
        v1 = a.get_vector()
        v2 = b.get_vector()
        v3 = c.get_vector()
        v4 = d.get_vector()

        # BioPythonのPPBuilderの二面角計算を使用
        ang = math.degrees(PPBuilder()._calc_dihedral(v1, v2, v3, v4))

        # -180～180の範囲に正規化
        if ang < -180:
            ang += 360
        if ang > 180:
            ang -= 360
        return ang
    except Exception:
        return None


# ─── 解析：cis結合検出 & Cα距離 ─────────────────────────────
def _detect_cis(residues: List) -> List[Tuple[int, int, float, str]]:
    """
    ペプチド結合のcis形式を検出

    ω角（omega）が-30°～30°の範囲をcisと判定
    ノートブックの`run_DSA`関数のロジックを移植

    Returns:
        List of (res_i, res_j, omega_deg, resname_j)
    """
    out = []
    for i in range(len(residues) - 1):
        r_i = residues[i]
        r_ip1 = residues[i + 1]

        # ω角の計算に必要な原子: C(i) - N(i+1) - CA(i+1) - C(i+1)
        C = _atom(r_i, "C")
        N = _atom(r_ip1, "N")
        CA = _atom(r_ip1, "CA")
        C2 = _atom(r_ip1, "C")

        if C is None or N is None or CA is None or C2 is None:
            continue

        omega = _dihedral(C, N, CA, C2)
        if omega is None:
            continue

        # cis判定: -30° <= ω <= 30°
        if -30.0 <= omega <= 30.0:
            out.append((i, i + 1, float(omega), _resname(r_ip1)))

    return out


def _ca_coords(residues: List) -> np.ndarray:
    """
    Cα座標を抽出してNumPy配列で返す

    Returns:
        Shape (n_residues, 3) の配列
    """
    coords = []
    for r in residues:
        ca = _atom(r, "CA")
        if ca is None:
            coords.append([np.nan, np.nan, np.nan])
        else:
            v = ca.get_coord()
            coords.append([float(v[0]), float(v[1]), float(v[2])])
    return np.asarray(coords, dtype=float)


def _pairwise_dist(C: np.ndarray) -> np.ndarray:
    """
    Cα座標からペアワイズ距離行列を計算

    Args:
        C: Shape (n, 3) の座標配列

    Returns:
        Shape (n, n) の距離行列（対角はNaN）
    """
    n = C.shape[0]
    D = np.full((n, n), np.nan, dtype=float)

    # NaNでない座標のみを抽出
    mask = ~np.isnan(C).any(axis=1)
    idx = np.where(mask)[0]

    if len(idx) == 0:
        return D

    X = C[idx]

    # ユークリッド距離の高速計算
    aa = np.sum(X * X, axis=1, keepdims=True)
    bb = aa.T
    ab = X @ X.T
    d2 = np.maximum(aa + bb - 2 * ab, 0.0)
    d = np.sqrt(d2, dtype=float)

    # 結果を元の行列に戻す
    for ii, i in enumerate(idx):
        for jj, j in enumerate(idx):
            D[i, j] = np.nan if i == j else d[ii, jj]

    return D


def _pad_square(D: np.ndarray, nmax: int) -> np.ndarray:
    """
    行列を正方形にパディング
    """
    R = np.full((nmax, nmax), np.nan, dtype=float)
    n = D.shape[0]
    R[:n, :n] = D
    return R


def _mean_matrix(mats: List[np.ndarray]) -> np.ndarray:
    """
    複数の距離行列の平均を計算（ノートブックの統合処理）
    """
    if not mats:
        return np.empty((0, 0))

    nmax = max(m.shape[0] for m in mats)
    padded = [_pad_square(m, nmax) for m in mats]
    stack = np.stack(padded, axis=0)

    with np.errstate(invalid="ignore"):
        mean = np.nanmean(stack, axis=0)

    return mean


# ─── ヒートマップ生成（ノートブックスタイル） ───────────────
def _generate_heatmap(
    distance_matrix: np.ndarray,
    uniprot_id: str,
    output_path: Path,
    title: str = "Mean Cα Distance",
):
    """
    距離行列からヒートマップを生成

    ノートブックのseaborn + matplotlibスタイルを再現
    """
    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)

    # NaNをマスク
    masked = np.ma.masked_invalid(distance_matrix)

    # seabornでヒートマップ作成
    sns.heatmap(
        masked,
        vmax=130,
        vmin=20,
        square=True,
        center=75,
        cmap="rainbow_r",
        cbar=True,
        cbar_kws={"label": "Cα distance (Å)"},
        ax=ax,
    )

    ax.set_title(f"{uniprot_id} - {title}", fontsize=14, fontweight="bold")
    ax.set_xlabel("Residue index", fontsize=12)
    ax.set_ylabel("Residue index", fontsize=12)

    plt.tight_layout()
    fig.savefig(output_path, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ─── 表面API（FastAPIから呼ぶ） ──────────────────────────
def analyze_structures(
    uniprot_id: str,
    method: str,
    seq_ratio: int,
    cif_paths: List[str],
    alphafold: List[Dict],
) -> Dict[str, Any]:
    """
    タンパク質構造解析のメイン関数

    研究室のノートブック処理を完全自動化

    Args:
        uniprot_id: UniProt ID
        method: 実験手法 (X-ray, EM, NMR)
        seq_ratio: シーケンス比率パラメータ
        cif_paths: mmCIFファイルのパスリスト
        alphafold: AlphaFoldモデルのメタデータリスト

    Returns:
        解析結果の辞書
    """
    # 出力ディレクトリ作成
    outdir = ART_ROOT / uniprot_id
    _ensure_dir(outdir)

    # 構造ファイル読み込み
    structs = _load_structures(cif_paths)

    if not structs:
        return {
            "uniprot_id": uniprot_id,
            "method": method,
            "seq_ratio": seq_ratio,
            "kpi": {"cis_count": 0, "midrange_dist_fraction": None},
            "artifacts": {},
            "note": "No valid mmCIF parsed",
            "ts": datetime.utcnow().isoformat() + "Z",
        }

    # 各構造から残基リストを抽出
    reslists = []
    for s in structs:
        residues = _collect_residues(s)

        # 計算量対策: 1200残基まで
        if len(residues) > 1200:
            residues = residues[:1200]

        if residues:
            reslists.append(residues)

    if not reslists:
        return {
            "uniprot_id": uniprot_id,
            "method": method,
            "seq_ratio": seq_ratio,
            "kpi": {"cis_count": 0, "midrange_dist_fraction": None},
            "artifacts": {},
            "note": "No valid residues found",
            "ts": datetime.utcnow().isoformat() + "Z",
        }

    # cis結合検出（最初の構造のみ）
    cis_list = _detect_cis(reslists[0])

    # 各構造の距離行列を計算
    mats = []
    for residues in reslists:
        C = _ca_coords(residues)
        D = _pairwise_dist(C)
        mats.append(D)

    # 平均距離行列を計算
    meanD = _mean_matrix(mats)

    # ─── アーティファクト生成 ───────────────────────

    # 1. cis結合のCSV
    cis_df = pd.DataFrame(
        cis_list, columns=["res_i", "res_j", "omega_deg", "resname_j"]
    )
    cis_csv = outdir / "cis.csv"
    cis_df.to_csv(cis_csv, index=False)

    # 2. 距離行列のCSV
    n = meanD.shape[0]
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            val = float(meanD[i, j]) if not np.isnan(meanD[i, j]) else np.nan
            rows.append((i + 1, j + 1, val))

    dist_df = pd.DataFrame(rows, columns=["i", "j", "mean_ca_dist"])
    dist_csv = outdir / "distance.csv"
    dist_df.to_csv(dist_csv, index=False)

    # 3. ヒートマップ（ノートブックスタイル）
    heat_png = outdir / "heatmap.png"
    _generate_heatmap(meanD, uniprot_id, heat_png)

    # ─── KPI計算 ───────────────────────────────────

    # 中距離割合（10-20Åの割合）
    with np.errstate(invalid="ignore"):
        finite = meanD[np.isfinite(meanD)]
        pct_mid = (
            float(np.mean((finite >= 10.0) & (finite <= 20.0)))
            if finite.size > 0
            else None
        )

    # cis結合の平均距離とスコア（研究室の指標）
    cis_distances = []
    if cis_list:
        for res_i, res_j, omega, resname in cis_list:
            if res_i < meanD.shape[0] and res_j < meanD.shape[0]:
                d = meanD[res_i, res_j]
                if not np.isnan(d):
                    cis_distances.append(d)

    mean_cis_dist = float(np.mean(cis_distances)) if cis_distances else None
    std_cis_dist = float(np.std(cis_distances)) if cis_distances else None

    # ─── 結果のJSON保存 ───────────────────────────
    results_json = outdir / "results.json"
    summary = {
        "uniprot_id": uniprot_id,
        "method": method,
        "seq_ratio": seq_ratio,
        "inputs": {
            "pdb_cifs": [str(Path(p).name) for p in cif_paths],
            "alphafold_models": len(alphafold),
        },
        "kpi": {
            "cis_count": int(len(cis_list)),
            "midrange_dist_fraction": pct_mid,
            "mean_cis_distance": mean_cis_dist,
            "std_cis_distance": std_cis_dist,
            "total_structures": len(structs),
            "total_residues": reslists[0] and len(reslists[0]) or 0,
        },
        "artifacts": {
            "results_json": str(results_json),
            "cis_csv": str(cis_csv),
            "distance_csv": str(dist_csv),
            "heatmap_png": str(heat_png),
        },
        "ts": datetime.utcnow().isoformat() + "Z",
    }

    with results_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    return summary
