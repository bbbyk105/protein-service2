"""
FastAPI Main Application

研究室のJupyterノートブックを完全自動化するRESTful API
"""

import re
from typing import List, Optional, Union
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from app.services.uniprot import fetch_uniprot_core, fetch_pdb_ids
from app.services.alphafold import fetch_alphafold_with_cif
from app.services.pdb import download_mmCIF_batch
from app.compute import analyze_structures
from app.util import uniq


# ─── FastAPI App ───────────────────────────────────────
app = FastAPI(
    title="Protein Structure Analysis Service",
    description="研究室のタンパク質構造解析を完全自動化",
    version="1.0.0",
)

# CORS設定（Next.jsなどフロントエンドからのアクセスを許可）
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 本番環境では制限すること
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ─── Models ────────────────────────────────────────────
class AnalyzeReq(BaseModel):
    """解析リクエスト"""

    uniprot_ids: Union[List[str], str] = Field(
        ...,
        description="UniProt IDs (array or comma/space/newline separated)",
        examples=[["P17456", "C6H0Y9"], "P17456,C6H0Y9", "P17456 C6H0Y9"],
    )
    method: Union[List[str], str] = Field(
        "X-ray",
        description="Experimental method(s): X-ray | EM | NMR (can specify multiple)",
        examples=["X-ray", ["X-ray", "EM"], "X-ray,EM"],
    )
    seq_ratio: int = Field(
        20, ge=1, le=100, description="Sequence ratio parameter (研究室パラメータ)"
    )
    negative_pdbids: Optional[Union[List[str], str]] = Field(
        None, description="PDB IDs to exclude from analysis"
    )
    forced_pdbids: Optional[Union[List[str], str]] = Field(
        None,
        description="Force specific PDB IDs (e.g., '1ROB,7RSA' for RNase A cis-Pro validation)",
    )


class AnalyzeRes(BaseModel):
    """解析レスポンス"""

    results: List[dict]
    warnings: List[str] = []


# ─── Helper Functions ──────────────────────────────────
def _to_list(v: Union[List[str], str, None]) -> List[str]:
    """
    文字列またはリストを正規化してリストに変換

    Args:
        v: 入力値

    Returns:
        正規化されたリスト
    """
    if v is None:
        return []
    if isinstance(v, list):
        return [s.strip() for s in v if s and s.strip()]

    # カンマ、スペース、改行で分割
    parts = re.split(r"[, \n]+", str(v))
    return [s.strip() for s in parts if s and s.strip()]


# ─── Endpoints ─────────────────────────────────────────
@app.get("/health")
def health():
    """ヘルスチェック"""
    return {"status": "ok", "service": "Protein Structure Analysis", "version": "1.0.0"}


@app.post("/analyze", response_model=AnalyzeRes)
def analyze(req: AnalyzeReq):
    """
    タンパク質構造解析のメインエンドポイント

    研究室のノートブック処理を完全自動化：
    1. UniProtから基本情報取得
    2. PDB IDリスト取得
    3. mmCIFファイルダウンロード
    4. AlphaFoldモデル取得（PDBが不足の場合）
    5. cis結合・距離行列解析
    6. ヒートマップ生成
    7. 結果をCSV/JSON/PNG保存

    Args:
        req: 解析リクエスト

    Returns:
        解析結果とワーニング
    """
    warnings = []
    out = []

    # UniProt IDsの正規化
    ids = [s for s in _to_list(req.uniprot_ids) if s.lower() != "undefined"]
    ids = uniq(ids)

    if not ids:
        raise HTTPException(status_code=400, detail="No valid UniProt IDs provided")

    # 除外するPDB IDsのセット
    negset = set(s.upper() for s in _to_list(req.negative_pdbids))

    # 各UniProt IDを処理
    for uid in ids:
        print(f"\n{'='*60}")
        print(f"Processing: {uid}")
        print(f"{'='*60}")

        # 1. UniProt基本情報取得
        core = fetch_uniprot_core(uid)
        if not core:
            warnings.append(f"{uid}: UniProt entry not found")
            out.append(
                {
                    "uniprot_id": uid,
                    "status": "not_found",
                    "error": "UniProt entry not found",
                }
            )
            continue

        print(f"✓ UniProt: {core['name']} ({core['organism']})")
        print(f"  Length: {core['length']} aa")

        # methodsを正規化（文字列またはリスト）
        methods = _to_list(req.method) if req.method else ["X-ray"]
        if not methods:
            methods = ["X-ray"]

        # 2. PDB IDs取得
        # forced_pdbidsが指定されている場合は、それを優先
        forced = _to_list(req.forced_pdbids) if req.forced_pdbids else []

        if forced:
            # 強制指定されたPDB IDを使用
            all_pdb_ids = [p.upper() for p in forced if p.upper() not in negset]
            print(f"✓ PDB IDs (forced): {all_pdb_ids} ({len(all_pdb_ids)} entries)")
        else:
            # 通常の自動取得（全てのmethodを試す）
            all_pdb_ids = []
            for m in methods:
                pdb_ids = fetch_pdb_ids(uid, method=m, max_pdb=8, max_resolution=3.5)
                all_pdb_ids.extend(pdb_ids)

            # 重複削除 & 除外リスト適用
            all_pdb_ids = list(set(all_pdb_ids))
            all_pdb_ids = [p for p in all_pdb_ids if p.upper() not in negset]
            print(
                f"✓ PDB IDs ({', '.join(methods)}): {all_pdb_ids} ({len(all_pdb_ids)} entries)"
            )

        # 3. mmCIFダウンロード
        cif_paths = download_mmCIF_batch(all_pdb_ids)
        print(f"✓ Downloaded: {len(cif_paths)} mmCIF files")

        # 4. AlphaFold取得（PDBが不足の場合）
        af_models = []
        af_cif_paths = []
        if len(cif_paths) < 2:
            print(f"⚠ PDB不足 → AlphaFold取得中...")
            try:
                af_models, af_cif_paths = fetch_alphafold_with_cif(uid)
                print(f"✓ AlphaFold: {len(af_models)} models, {len(af_cif_paths)} CIFs")
            except Exception as e:
                print(f"✗ AlphaFold取得失敗: {e}")
                af_models = []
                af_cif_paths = []

        # 全CIFパスを統合
        all_cif_paths = cif_paths + af_cif_paths

        if not all_cif_paths:
            warnings.append(f"{uid}: no PDB/AlphaFold models available")

        # 5. 構造解析実行
        print(f"▶ 解析開始...")
        try:
            analysis = analyze_structures(
                uniprot_id=uid,
                method=", ".join(methods),  # 複数methodを結合
                seq_ratio=req.seq_ratio,
                cif_paths=all_cif_paths,
                alphafold=af_models,
            )
            print(f"✓ 解析完了")
            print(f"  - cis結合: {analysis['kpi']['cis_count']} 個")
            print(f"  - 中距離割合: {analysis['kpi'].get('midrange_dist_fraction')}")

        except Exception as e:
            print(f"✗ 解析エラー: {e}")
            warnings.append(f"{uid}: Analysis failed - {str(e)}")
            analysis = {
                "error": str(e),
                "kpi": {"cis_count": 0, "midrange_dist_fraction": None},
            }

        # 6. 結果格納
        out.append(
            {
                "uniprot_id": uid,
                "core": core,
                "pdb_ids": all_pdb_ids,
                "methods": methods,
                "alphafold_models": len(af_models),
                "analysis": analysis,
                "status": "ok" if "error" not in analysis else "error",
            }
        )

    print(f"\n{'='*60}")
    print(f"✓ 全処理完了: {len(out)} entries")
    print(f"{'='*60}\n")

    return {"results": out, "warnings": warnings}


@app.get("/")
def root():
    """ルートエンドポイント"""
    return {
        "service": "Protein Structure Analysis Service",
        "version": "1.0.0",
        "description": "研究室のタンパク質構造解析を完全自動化",
        "endpoints": {
            "/health": "ヘルスチェック",
            "/analyze": "構造解析（POST）",
            "/docs": "API仕様書（Swagger UI）",
        },
        "usage": {
            "example": {
                "method": "POST",
                "url": "/analyze",
                "body": {
                    "uniprot_ids": ["P17456", "C6H0Y9"],
                    "method": "X-ray",
                    "seq_ratio": 20,
                },
            }
        },
    }


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
