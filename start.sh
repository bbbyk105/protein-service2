#!/bin/bash
# FastAPI起動スクリプト

echo "🚀 Starting Protein Analysis Service..."
echo ""

# 仮想環境の確認
if [ ! -d ".venv" ]; then
    echo "⚠️  仮想環境が見つかりません"
    echo "   以下を実行してください:"
    echo "   python -m venv .venv"
    echo "   source .venv/bin/activate"
    echo "   pip install -r requirements.txt"
    exit 1
fi

# 起動
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
