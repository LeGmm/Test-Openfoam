# TODO · caseCHT_v89（v8_9 结构，未完成，勿直接开跑）

状态：**网格已转换并体检通过，但缺 `0/` 字典，preflight 会报 layout FAIL（缺 0/fluid、0/solid）**。
本目录放在独立分支 `v8_9-setup`，不影响 master。

## 已完成（2026-10-21）
- 源文件体检 `tools/inspect_cas_h5.py`：122.34 MB / Fluent 23.2，节点 3,068,872、
  面 4,092,408（内部 3,927,457 / 边界 164,951）、单元 705,138；
  **0 负体积 / 0 近零体积 / 0 非法面 / 0 孤立单元**；总体积 106.07 cm³。
- 双区域转换 `tools/fluent_cas_to_foam_cht.py`：fluid **489,728** 单元（4 补丁）、
  solid **215,410** 单元（3 补丁）。
  补丁名完整识别：fluid `inlet/outlet/fluid_to_solid/wall`，solid `solid_to_fluid/coldWall/outer`
  —— 源 cas.h5 的 11 个面区里 10 个名字是 `zone<id>` 占位符，转换器靠关键字仍正确归类，
  此前「边界识别不可靠」的判定解除。
- `checkMesh -region fluid` / `-region solid` 均 **Mesh OK**。
- 几何同口径核对（与 v8_1/v8_8 一致）：coldWall **16.000 cm²**、outer 128.628 cm²、
  inlet/outlet 均 **1.130 cm²**；流固界面 **136.187 cm²**（比 v8_8 大 10.2%）。
- `system/` 与 `constant/regionProperties` 从 caseCHT_v88 复制（求解设置按同口径共享）。

## 待办
1. **缺 `0/`**：需要 CHT 字典生成器（待办 T6）或从 caseCHT_v88 手工拷 0/ 后把
   `interface`/`fluid_to_solid`/`solid_to_fluid` 的 nonuniform 列表改成 uniform 占位
   （70,511 面，与 v8_8 的 62,158 不一致，直接拷会 FATAL）。
2. 0/ 就绪后跑 `python tools/preflight.py caseCHT_v89`，须全 PASS。
3. 求解按 README §3.7（decomposePar → mpiexec → reconstructPar → condensation_metrics）。
4. 跑完后与 v8_8 对比：两者网格尺度只差 **1.10×**，是当前最干净的一对可比结构；
   v8_1 与它们差 7.6×，需先加密（待办 T2/T18）。
