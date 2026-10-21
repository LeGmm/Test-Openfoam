#!/usr/bin/env bash
# ============================================================================
#  v8_9~v8_14 六结构 Step1 CHT 预处理：体检 → 转换 → 耦合 → checkMesh → preflight
#  用法：bash tools/_prep_six.sh           （在 OpenFOAM/ 目录下执行）
#  日志：log.prep_<v>.txt（每结构一个，gitignore 已排除 log.*）
# ============================================================================
set -u
PY="C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
ROOT="F:/B_Simulation/fangzhen-3"
OF="$ROOT/OpenFOAM"
TEMPLATE="$OF/caseCHT_v88"

declare -A CAS=(
  [v8_9]="$ROOT/Qishi_v8_9_20261012/qishi_v8_9_files/dp0/FLTG/Fluent/FLTG.1-3.cas.h5"
  [v8_10]="$ROOT/Qishi_v8_10_20261014/qishi-v8-10_files/dp0/FLTG/Fluent/FLTG-2.cas.h5"
  [v8_11]="$ROOT/Qishi_v8_11_20261014/qishi-v8-11_files/dp0/FLTG/Fluent/FLTG-2.cas.h5"
  [v8_12]="$ROOT/Qishi_v8_12_20261014/qishi-v8-12_files/dp0/FLTG/Fluent/FLTG-2.cas.h5"
  [v8_13]="$ROOT/Qishi_v8_13_20261014/qishi-v8-13_files/dp0/FLTG/Fluent/FLTG-4.cas.h5"
  [v8_14]="$ROOT/Qishi_v8_14_20261014/qishi-v8-14_files/dp0/FLTG/Fluent/FLTG-7.cas.h5"
)

fail=0
for v in v8_9 v8_10 v8_11 v8_12 v8_13 v8_14; do
  case_dir="$OF/caseCHT_${v}"
  cas="${CAS[$v]}"
  log="$OF/log.prep_${v}.txt"
  {
    echo "##################################################"
    echo "# $v   $(date)"
    echo "# cas: $cas"
    echo "##################################################"

    if [ ! -f "$cas" ]; then echo "FATAL: 源文件不存在 $cas"; fail=1; continue; fi

    echo "=== [1/6] inspect_cas_h5 体检 ==="
    "$PY" "$OF/tools/inspect_cas_h5.py" "$cas" --json "$OF/_inspect_${v}.json" 2>&1 | tail -25
    echo "inspect_rc=$?"

    echo "=== [2/6] 转换（双区域）==="
    if [ -d "$case_dir" ]; then echo "警告：$case_dir 已存在，先移走旧目录"; mv "$case_dir" "${case_dir}_bak_$(date +%H%M%S)"; fi
    "$PY" "$OF/tools/fluent_cas_to_foam_cht.py" "$cas" "$case_dir" 2>&1 | tail -12
    rc=$?; echo "convert_rc=$rc"
    if [ $rc -ne 0 ] || [ ! -d "$case_dir/constant/fluid/polyMesh" ]; then
      echo "FATAL: 转换失败，跳过 $v"; fail=1; continue
    fi

    echo "=== [3/6] 拷 system/ + regionProperties + 0/（模板 caseCHT_v88）==="
    cp -r "$TEMPLATE/system" "$case_dir/system"
    cp "$TEMPLATE/constant/regionProperties" "$case_dir/constant/regionProperties"
    cp -r "$TEMPLATE/0" "$case_dir/0"
    echo "copied."

    echo "=== [4/6] setup_cht_coupling（界面统一 interface + 共轭 T）==="
    "$PY" "$OF/tools/setup_cht_coupling.py" "$case_dir" 2>&1
    echo "coupling_rc=$?"

    echo "=== [5/6] checkMesh 两区 ==="
    "$PY" "$OF/tools/of8.py" "$case_dir" "checkMesh -region fluid" 2>&1 | grep -E "Mesh OK|FAILED|FATAL|errors|cells:|bounds" | head -8
    "$PY" "$OF/tools/of8.py" "$case_dir" "checkMesh -region solid" 2>&1 | grep -E "Mesh OK|FAILED|FATAL|errors|cells:" | head -8

    echo "=== [6/6] preflight 门禁 ==="
    "$PY" "$OF/tools/preflight.py" "$case_dir" 2>&1 | tail -30
    echo "preflight_rc=$?"

    echo "=== $v 预处理完成 $(date) ==="
  } >> "$log" 2>&1
  tail -3 "$log"
done

echo "ALL_DONE fail=$fail"
