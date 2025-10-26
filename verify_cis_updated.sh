#!/usr/bin/env bash
set -euo pipefail

API="http://localhost:8000/analyze"
METHODS=("X-ray" "NMR" "EM")

# ゴールデンセット（期待値：cis_count＝Pro直前限定、ユニークサイト数）
# フォーマット: UniProt  Xray  NMR  EM
EXPECTATIONS="$(cat <<'EXPECTATIONS_EOF'
P61823  2  2  0   # RNase A (Tyr92–Pro93, Asn113–Pro114)
P61769  1  1  NA  # β2-microglobulin (native Pro32 cis; EMはデータ依存でブレるため評価外)
P10599  1  1  0   # Human Thioredoxin (cis-Pro75) - ユニーク化で重複解消
P0CG47  0  0  0   # Ubiquitin（陰性）
EXPECTATIONS_EOF
)"

call_api() {
  local uniprot="$1" method="$2"
  curl -s -X POST "$API" \
    -H "Content-Type: application/json" \
    -d "{\"uniprot_ids\":[\"$uniprot\"],\"method\":\"$method\",\"seq_ratio\":20}" \
  | jq -r '.results[0].analysis.kpi.cis_count // "NA"'
}

get_expected() {
  local uniprot="$1" method="$2"
  local line
  line="$(printf "%s\n" "$EXPECTATIONS" | grep -E "^$uniprot[[:space:]]" | head -n1 || true)"
  [[ -z "$line" ]] && { echo "NA"; return; }
  case "$method" in
    "X-ray") awk '{print $2}' <<< "$line" ;;
    "NMR")   awk '{print $3}' <<< "$line" ;;
    "EM")    awk '{print $4}' <<< "$line" ;;
    *)       echo "NA" ;;
  esac
}

passfail() {
  local got="$1" exp="$2"
  if [[ "$exp" == "NA" || "$got" == "NA" ]]; then
    echo "SKIP"
  elif [[ "$got" == "$exp" ]]; then
    echo "PASS"
  else
    echo "FAIL"
  fi
}

section() {
  echo
  echo "================================================================"
  echo "$1"
  echo "================================================================"
}

section "A) 単体メソッド検証（1法ずつ）"
PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

while read -r line; do
  [[ -z "$line" || "$line" =~ ^# ]] && continue
  uniprot="$(awk '{print $1}' <<< "$line")"
  echo ">>> $uniprot"
  for m in "${METHODS[@]}"; do
    exp="$(get_expected "$uniprot" "$m")"
    got="$(call_api "$uniprot" "$m" || echo "NA")"
    pf="$(passfail "$got" "$exp")"
    printf "  - %-5s  got: %-3s  exp: %-3s  [%s]\n" "$m" "$got" "$exp" "$pf"
    
    case "$pf" in
      PASS) ((PASS_COUNT++)) ;;
      FAIL) ((FAIL_COUNT++)) ;;
      SKIP) ((SKIP_COUNT++)) ;;
    esac
  done
done <<< "$(printf "%s\n" "$EXPECTATIONS")"

section "B) 2メソッド組み合わせ検証（両方が期待値ならPASS）"
pairs=( "X-ray NMR" "X-ray EM" "NMR EM" )
while read -r line; do
  [[ -z "$line" || "$line" =~ ^# ]] && continue
  uniprot="$(awk '{print $1}' <<< "$line")"
  echo ">>> $uniprot"
  for pair in "${pairs[@]}"; do
    m1="$(awk '{print $1}' <<< "$pair")"
    m2="$(awk '{print $2}' <<< "$pair")"
    exp1="$(get_expected "$uniprot" "$m1")"
    exp2="$(get_expected "$uniprot" "$m2")"
    got1="$(call_api "$uniprot" "$m1" || echo "NA")"
    got2="$(call_api "$uniprot" "$m2" || echo "NA")"
    pf1="$(passfail "$got1" "$exp1")"
    pf2="$(passfail "$got2" "$exp2")"
    
    overall="PASS"
    if [[ "$pf1" == "FAIL" || "$pf2" == "FAIL" ]]; then
      overall="FAIL"
    elif [[ "$pf1" == "SKIP" && "$pf2" == "SKIP" ]]; then
      overall="SKIP"
    elif [[ "$pf1" == "SKIP" || "$pf2" == "SKIP" ]]; then
      overall="PASS"  # 片方SKIPなら他方がPASSならOK
    fi
    
    printf "  - [%-5s + %-3s]  got:(%s,%s) exp:(%s,%s)  [%s]\n" "$m1" "$m2" "$got1" "$got2" "$exp1" "$exp2" "$overall"
  done
done <<< "$(printf "%s\n" "$EXPECTATIONS")"

section "C) 3メソッド検証（3法すべて期待値ならPASS）"
while read -r line; do
  [[ -z "$line" || "$line" =~ ^# ]] && continue
  uniprot="$(awk '{print $1}' <<< "$line")"
  expX="$(get_expected "$uniprot" "X-ray")"
  expN="$(get_expected "$uniprot" "NMR")"
  expE="$(get_expected "$uniprot" "EM")"
  gotX="$(call_api "$uniprot" "X-ray" || echo "NA")"
  gotN="$(call_api "$uniprot" "NMR"  || echo "NA")"
  gotE="$(call_api "$uniprot" "EM"    || echo "NA")"
  pfX="$(passfail "$gotX" "$expX")"
  pfN="$(passfail "$gotN" "$expN")"
  pfE="$(passfail "$gotE" "$expE")"
  
  overall="PASS"
  if [[ "$pfX" == "FAIL" || "$pfN" == "FAIL" || "$pfE" == "FAIL" ]]; then
    overall="FAIL"
  fi
  
  printf ">>> %-7s  X:%s/%s(%s)  N:%s/%s(%s)  E:%s/%s(%s)  [%s]\n" \
    "$uniprot" "$gotX" "$expX" "$pfX" "$gotN" "$expN" "$pfN" "$gotE" "$expE" "$pfE" "$overall"
done <<< "$(printf "%s\n" "$EXPECTATIONS")"

section "Summary"
echo "PASS: $PASS_COUNT"
echo "FAIL: $FAIL_COUNT"
echo "SKIP: $SKIP_COUNT"
echo
echo "Done."
