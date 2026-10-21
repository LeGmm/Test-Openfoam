#!/usr/bin/env bash
# ============================================================================
#  v8_9~v8_14 六结构 Step1 CHT 战役监视器
#  职责：① 等 _prep_six.sh 六结构预处理全部完成（带超时）
#        ② 按预处理结果分两车道（各 -np 8，双路并行）开跑求解
#        ③ 守候两车道结束，输出战役总结
#  日志：log.campaign_lane1.txt / log.campaign_lane2.txt / log.campaign.txt
# ============================================================================
set -u
PY="C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
OF="F:/B_Simulation/fangzhen-3/OpenFOAM"
cd "$OF"

ALL="v8_9 v8_10 v8_11 v8_12 v8_13 v8_14"
MAX_WAIT_MIN=200          # 预处理最长等待（分钟）

log() { echo "[campaign $(date +%H:%M)] $*"; }

# ---------- ① 等预处理 ----------
t0=$(date +%s)
while true; do
  done_cnt=0
  for v in $ALL; do
    grep -q "预处理完成" "log.prep_${v}.txt" 2>/dev/null && done_cnt=$((done_cnt+1))
  done
  waited=$(( ($(date +%s) - t0) / 60 ))
  log "预处理完成 $done_cnt/6（已等 ${waited} 分钟）"
  [ "$done_cnt" -eq 6 ] && break
  [ "$waited" -ge "$MAX_WAIT_MIN" ] && { log "等待超时，按已完成的开跑"; break; }
  sleep 60
done

# ---------- ② 分级：通过 / 未过 ----------
passed=""; blocked=""
for v in $ALL; do
  if grep -q "预处理完成" "log.prep_${v}.txt" 2>/dev/null && ! grep -q "FAIL" "log.prep_${v}.txt" 2>/dev/null; then
    passed="$passed $v"
  else
    blocked="$blocked $v"
    log "$v 未过门禁或无产物，本轮不开跑（见 log.prep_${v}.txt）"
  fi
done
log "开跑清单:$passed ；阻断:$blocked"

# ---------- ③ 分车道（交错的负载/风险搭配） ----------
lane1=""; lane2=""
i=0
for v in $passed; do
  if [ $((i % 2)) -eq 0 ]; then lane1="$lane1 $v"; else lane2="$lane2 $v"; fi
  i=$((i+1))
done
log "lane1:$lane1"
log "lane2:$lane2"

# ---------- ④ 发射两车道 ----------
pids=""
[ -n "$(echo $lane1)" ] && { bash tools/_solve_lane.sh 1 $lane1 > log.campaign_lane1.txt 2>&1 & pids="$pids $!"; }
[ -n "$(echo $lane2)" ] && { bash tools/_solve_lane.sh 2 $lane2 > log.campaign_lane2.txt 2>&1 & pids="$pids $!"; }
log "车道已发射 pids:$pids"

# ---------- ⑤ 守候 ----------
for p in $pids; do wait "$p" 2>/dev/null; done
log "两车道全部结束"

echo "================ 战役总结 $(date) ================"
for v in $ALL; do
  last_t=$(ls "caseCHT_${v}" 2>/dev/null | grep -E '^[0-9]+$' | sort -n | tail -1)
  if [ -f "_m_${v}.json" ]; then
    line=$("$PY" - "$v" <<'PYEOF' 2>/dev/null
import json, sys
v = sys.argv[1]
try:
    d = json.load(open(f"_m_{v}.json", encoding="utf-8"))
    z = d.get("zoned", {})
    print(f"UA={d.get('UA_W_per_K', float('nan')):.5f} W/K  冷凝={z.get('m_cond_mg_h', float('nan')):.3f} mg/h  "
          f"w_out={z.get('w_out_g_per_kg', float('nan')):.3f} g/kg  T_out_mw={d.get('T_out_mw', float('nan')):.3f} K")
except Exception as e:
    print(f"指标读取失败: {e}")
PYEOF
)
    echo "$v: 时间步=${last_t:-无}  $line"
  else
    echo "$v: 无指标（时间步=${last_t:-无}）"
  fi
done
echo "=================================================="
log "CAMPAIGN_ALL_DONE"
