# 🧬 Protein Structure Analysis Service

研究室のJupyterノートブック（`DSA_Cis_250317.ipynb`）を**完全自動化**するFastAPIサービス

## 🎯 概要

このシステムは、岡田研究室で手動実行していたタンパク質構造解析を完全に自動化します：

### 従来の方法
- ✋ UniProt IDを1つずつ手入力
- ✋ AlphaFold/PDBから手動ダウンロード
- ✋ Jupyterノートブックで1つずつ計算
- ✋ 結果を手動で保存

### このシステム
- 🚀 API `/analyze` にIDを投げるだけ
- 🚀 すべて自動処理
- 🚀 成果物（CSV, heatmap, JSON）が自動保存
- 🚀 複数IDを並列処理

---

## 📁 プロジェクト構造

```
protein-service/
├── app/
│   ├── __init__.py
│   ├── main.py              # FastAPIエントリーポイント
│   ├── compute.py           # 構造解析エンジン（ノートブックロジック）
│   ├── util.py              # ユーティリティ関数
│   └── services/
│       ├── __init__.py
│       ├── uniprot.py       # UniProt API
│       ├── pdb.py           # PDB mmCIFダウンロード
│       └── alphafold.py     # AlphaFold API + CIFダウンロード
├── data/
│   ├── artifacts/           # 解析結果（UniProt ID別）
│   │   └── P17456/
│   │       ├── results.json
│   │       ├── cis.csv
│   │       ├── distance.csv
│   │       └── heatmap.png
│   ├── pdb/                 # PDB mmCIFキャッシュ
│   └── alphafold/           # AlphaFold CIFキャッシュ
├── requirements.txt
├── start.sh                 # 起動スクリプト
├── test_api.py              # テストスクリプト
└── README.md
```

---

## 🚀 セットアップ

### 1. 環境構築

```bash
# リポジトリクローン（または作業ディレクトリに移動）
cd protein-service

# 仮想環境作成
python -m venv .venv

# 仮想環境アクティベート
source .venv/bin/activate  # macOS/Linux
# または
.venv\Scripts\activate     # Windows

# 依存関係インストール
pip install -r requirements.txt
```

### 2. サーバー起動

```bash
# 方法1: 起動スクリプト
./start.sh

# 方法2: 直接起動
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

起動後、以下にアクセス：
- **API**: http://localhost:8000
- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

---

## 📖 使い方

### API エンドポイント

#### 1. ヘルスチェック
```bash
curl http://localhost:8000/health
```

#### 2. 構造解析

**単一ID**
```bash
curl -X POST http://localhost:8000/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "uniprot_ids": "P17456",
    "method": "X-ray",
    "seq_ratio": 20
  }'
```

**複数ID**
```bash
curl -X POST http://localhost:8000/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "uniprot_ids": ["P17456", "C6H0Y9"],
    "method": "X-ray",
    "seq_ratio": 20
  }'
```

**カンマ区切り**
```bash
curl -X POST http://localhost:8000/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "uniprot_ids": "P17456,C6H0Y9",
    "method": "X-ray",
    "seq_ratio": 20
  }'
```

### パラメータ

| パラメータ | 型 | 説明 | デフォルト |
|-----------|-----|------|-----------|
| `uniprot_ids` | string or array | UniProt ID（複数可） | **必須** |
| `method` | string | 実験手法: `X-ray`, `EM`, `NMR` | `X-ray` |
| `seq_ratio` | int | シーケンス比率パラメータ (1-100) | `20` |
| `negative_pdbids` | string or array | 除外するPDB ID | `null` |

---

## 🧪 テスト

### 自動テスト実行
```bash
# FastAPIが起動していることを確認してから
python test_api.py
```

### 期待される出力例
```
🧪 FastAPI Test Suite
============================================================

TEST 1: Health Check
Status: 200
Response: {"status": "ok", ...}

TEST 2: Single UniProt ID Analysis
Request: {"uniprot_ids": "P17456", ...}
Status: 200
Time: 5.23s

Results:
  UniProt: P17456
  Name: Capsid protein
  Organism: Cymbidium ringspot virus
  cis count: 3
  midrange fraction: 0.42

============================================================
✅ All tests passed!
============================================================
```

---

## 📊 出力ファイル

解析結果は `data/artifacts/<UniProtID>/` に自動保存されます：

### 1. `results.json`
```json
{
  "uniprot_id": "P17456",
  "method": "X-ray",
  "seq_ratio": 20,
  "kpi": {
    "cis_count": 3,
    "midrange_dist_fraction": 0.42,
    "mean_cis_distance": 12.5,
    "std_cis_distance": 2.1,
    "total_structures": 5,
    "total_residues": 380
  },
  "artifacts": {
    "cis_csv": "...",
    "distance_csv": "...",
    "heatmap_png": "..."
  },
  "ts": "2025-10-26T06:00:00Z"
}
```

### 2. `cis.csv`
cis結合のリスト
```csv
res_i,res_j,omega_deg,resname_j
45,46,-12.3,PRO
123,124,8.7,PRO
```

### 3. `distance.csv`
Cα間距離行列（上三角のみ）
```csv
i,j,mean_ca_dist
1,2,3.8
1,3,6.2
2,3,3.9
...
```

### 4. `heatmap.png`
距離行列のヒートマップ（研究室スタイル）

---

## 🔮 次のステップ

### 1. n8n連携
```
Next.js フォーム
    ↓
Webhook → n8n
    ↓
FastAPI /analyze
    ↓
結果を Slack/LINE/Gmail に通知
```

### 2. 自動実行フロー
- Notionデータベースと連携
- Google Drive自動アップロード
- 定期実行（cron）

### 3. Next.js UI
- UniProt ID入力フォーム
- リアルタイム進捗表示
- 結果のビジュアライゼーション

---

## 🛠️ トラブルシューティング

### エラー: "No valid mmCIF parsed"
- **原因**: PDB IDが見つからない、またはダウンロード失敗
- **解決**: AlphaFoldデータが自動で使用されます

### エラー: "Connection Error"
- **原因**: FastAPIが起動していない
- **解決**: `./start.sh` を実行

### エラー: "Module not found"
- **原因**: 依存関係がインストールされていない
- **解決**: `pip install -r requirements.txt`

---

## 📝 技術スタック

- **FastAPI**: RESTful API
- **BioPython**: タンパク質構造解析
- **NumPy/Pandas**: データ処理
- **Matplotlib/Seaborn**: ヒートマップ生成
- **Requests**: 外部API連携

---

## 📧 サポート

問題が発生した場合は、以下を確認してください：
1. `./start.sh` でサーバーが起動しているか
2. `python test_api.py` でテストが通るか
3. `data/artifacts/` に結果が保存されているか

---

## 🎓 研究室での利用

このシステムにより：
- ✅ 学部生でも「IDを打つだけ」で全構造解析が完了
- ✅ 研究データの再現性と効率化
- ✅ AlphaFold/PDBデータベースの更新に自動追従
- ✅ 今後の統合解析への拡張が容易

**Happy Analyzing! 🧬🔬**
