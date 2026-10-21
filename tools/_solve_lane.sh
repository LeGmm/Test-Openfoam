#!/usr/bin/env bash
# ============================================================================
#  单车道求解驱动：依次跑完本车道负责的结构
#  用法：bash tools/_solve_lane.sh <车道号> <结构...>
#    例：bash tools/_solve_lane.sh 1 v8_9 v8_10 v8_11
#  每结构：decomposePar → chtMultiRegionFoam(-np 8, 15000 外迭代)
#          → reconstructPar → condensation_metrics
#  前提：该结构已完成 _prep_six.sh 且 preflight 无 FAIL
# ============================================================================
set -u
PY="C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
OF="F:/B_Simulation/fangzhen-3/OpenFOAM"
cd "$OF"

lane=$1; shift
cases="$*"
echo "[lane $lane] start $(date) cases: $cases"

for v in $cases; do
  case_dir="caseCHT_${v}"
  plog="log.prep_${v}.txt"

  # ---- 门禁复核：预处理完成且 preflight 无 FAIL ----
  if [ ! -d "$case_dir/0" ] || [ ! -d "$case_dir/constant/fluid/polyMesh" ]; then
    echo "[lane $lane] $v SKIP：预处理产物不全（无 0/ 或 mesh）"; continue
  fi
  if grep -q "FAIL" "$plog" 2>/dev/null; then
    echo "[lane $lane] $v SKIP：prep 日志有 FAIL，人工检查后再跑"; continue
  fi

  echo "[lane $lane] ===== $v decomposePar $(date) ====="
  "$PY" tools/run_of_cmd.py "$case_dir" "decomposePar -allRegions -force" "log.decomp_${v}.txt"
  drc=$?
  if [ $drc -ne 0 ]; then echo "[lane $lane] $v decomposePar rc=$drc，跳过求解"; continue; fi

  echo "[lane $lane] ===== $v chtMultiRegionFoam 开跑 $(date) ====="
  "$PY" tools/run_of_cmd.py "$case_dir" "@_cmd_cht_v88.txt" "log.cht_${v}.txt"
  src=$?
  last_t=$(ls "$case_dir" 2>/dev/null | grep -E '^[0-9]+$' | sort -n | tail -1)
  echo "[lane $lane] $v 求解结束 rc=$src 最后时间步=${last_t:-无} $(date)"

  echo "[lane $lane] ===== $v reconstructPar $(date) ====="
  "$PY" tools/run_of_cmd.py "$case_dir" "reconstructPar -allRegions" "log.recon_${v}.txt"
  rrc=$?

  if [ -n "${last_t:-}" ]; then
    echo "[lane $lane] ===== $v condensation_metrics @${last_t} ====="
    "$PY" tools/condensation_metrics.py "$case_dir" "$last_t" --json "_m_${v}.json" 2>&1 | tail -5
  else
    echo "[lane $lane] $v 无时间步，跳过指标"
  fi
  echo "[lane $lane] $v 全流程结束 rc(solve=$src recon=$rrc) $(date)"
done

echo "[lane $lane] ALL_DONE $(date)"
