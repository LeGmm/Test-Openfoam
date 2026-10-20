#!/bin/bash
# ==========================================================================
#  一键跑完所有结构 + 出对比排名
#  用法：在 BlueCFD-Core 终端里执行
#     cd /f/B_Simulation/fangzhen-3/OpenFOAM
#     ./run_all.sh                # 跑全部
#     ./run_all.sh v8_1 v8_8      # 只跑指定结构
# ==========================================================================
set -u

ROOT="$(cd "$(dirname "$0")" && pwd)"
CASES="$ROOT/cases"
LOGDIR="$ROOT/logs"
mkdir -p "$LOGDIR"

# 结构列表：未指定参数则跑 cases/ 下全部
if [ $# -gt 0 ]; then
    STRUCTS="$@"
else
    STRUCTS=$(ls "$CASES")
fi

echo "================================================================"
echo " 结构列表: $STRUCTS"
echo "================================================================"

for s in $STRUCTS; do
    CASE="$CASES/$s"
    if [ ! -d "$CASE" ]; then
        echo "  [跳过] $s 不存在"
        continue
    fi
    echo ""
    echo "----------------------------------------------------------------"
    echo "### 结构 $s"
    echo "----------------------------------------------------------------"

    cd "$CASE" || continue

    # 1. 网格校验
    echo "  [1/2] checkMesh ..."
    checkMesh > "$LOGDIR/${s}_checkMesh.log" 2>&1
    if grep -qiE "Failed|FOAM FATAL" "$LOGDIR/${s}_checkMesh.log"; then
        echo "    !! checkMesh 失败，见 $LOGDIR/${s}_checkMesh.log"
        continue
    fi
    echo "    OK"

    # 2. 求解
    echo "  [2/2] multiphaseEulerFoam ..."
    multiphaseEulerFoam > "$LOGDIR/${s}_run.log" 2>&1
    if grep -qiE "FOAM FATAL ERROR|FOAM FATAL IO ERROR" "$LOGDIR/${s}_run.log"; then
        echo "    !! 求解失败，见 $LOGDIR/${s}_run.log"
        tail -20 "$LOGDIR/${s}_run.log"
        continue
    fi
    echo "    完成 -> $LOGDIR/${s}_run.log"

    # 3. 后处理指标
    echo "  提取指标 ..."
    postProcess -func "surfaceFieldValue(name=outlet,operation=sum,fields=(alphaPhi.liquid))" \
        > "$LOGDIR/${s}_pp.log" 2>&1 || true
done

cd "$ROOT"
echo ""
echo "================================================================"
echo " 全部结束，开始对比 ..."
echo "================================================================"
python tools/compare.py $STRUCTS || \
    "C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe" tools/compare.py $STRUCTS
