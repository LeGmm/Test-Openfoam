# 湿空气冷凝仿真 —— Fluent → OpenFOAM 移植

把 `routeA_condensation_setup.py`（Fluent 路线 A：Mixture 多相 + 组分输运 + Lee 相变）移植到
**OpenFOAM（multiphaseEulerFoam + 界面组分相变）**，用于多结构「冷凝效率」对比实验。

> 当前可跑的是 v8_1 两步共轭链路：`caseCHT_qishi/`（Step1，9000 次外迭代）+ `case_v81_1tau/`
> （Step2，1τ 诊断，冷壁 9459 面、T_int 均值 278.279 K）。v8_8 就绪算例为 `caseCHT_v88/`、`case_v88/`。
> 旧的 `case/`、`caseCHT/` 已于 2026-10-20 归档移出（见 `_trash_2026-10-20/`）。

---

## 一、目录结构

```
OpenFOAM/
├── README.md                  # 本文件
├── 安装指南.md                 # WSL2 + OpenFOAM 安装步骤
├── tools/
│   ├── fluent_cas_to_foam.py  # Fluent .cas.h5 → 单流体区 polyMesh（Step2 用，已实测验证）
│   ├── fluent_cas_to_foam_cht.py # Fluent .cas.h5 → fluid+solid 双区域 polyMesh（Step1 用）
│   ├── setup_cht_coupling.py  # 两侧界面统一改名 interface + T 写 coupledBaffleMixed（Step1 必做）
│   ├── align_step2_thermo.py  # Step2 热物性对齐（液相 multiComponent+Henry、气相 eConst）※生成算例后必跑
│   ├── ofmesh.py              # polyMesh 读写公共模块（新脚本统一从这里 import）
│   ├── preflight.py           # 开跑前门禁：语法/依赖/网格/拓扑/补丁名/controlDict/耦合/热物性
│   ├── map_interface_T.py     # Step1→Step2 界面温度映射（extract/extract-vtk/map/verify）
│   ├── mesh_topology_report.py # 判"一对单元多张面"是【真重复】还是【铰链】；--fix 只删真重复
│   ├── merge_dup_cells.py     # ⚠️ 已弃用：迭代并单元在 polyhedral 网格上会级联毁网格（64729→28378）
│   ├── fix_pinched_cells.py   # ⚠️ 已弃用：单对并细胞会让缺陷扩散（1 对 → 7 对 → …）
│   ├── drop_duplicate_faces.py # ⚠️ 已弃用：删面会让铰链单元开口（open cells 2）
│   ├── plot_z0_slice.py       # 读 surfaces FO 的切片 VTK，画 z=0 温度/流场（scipy griddata）
│   ├── fix_zero_volume_cells.py  # 零体积单元节点注入工具
│   ├── check_face_orientation.py # 网格朝向/排序自检（纯 Python，秒级）
│   ├── of8.py                  # 在本机 blueCFD(OpenFOAM 8) 环境里执行任意 OF 命令
│   ├── generate_case.py       # 参数化生成全套 0/ constant/ system/ 字典（含六项指标 FO）
│   ├── prepare_all_structures.py # 批量转换 8 个结构 + 生成字典
│   ├── ablate.py             # 单因素对照探针（求解发散/崩溃归因）
│   ├── mesh_stats.py / compare.py # 网格面积统计 / 跨结构指标对比
│   └── diagnostics/           # 排错期诊断脚本归档（audit_sources/dup_faces/liq_probe 等，见其 README）
├── caseCHT_qishi/             # Step1：v8_1 共轭传热（0…9000，能量比 1.080，待续算）
├── case_v81_1tau/             # Step2：v8_1 冷凝 1τ（T_int 已映射，当前诊断运行）
├── caseCHT_v88/               # Step1：v8_8 共轭（就绪，未跑）
├── case_v88/                  # Step2：v8_8（字典待重生成 + 映射后再跑）
└── tests/                     # 管线验证留存（bcDissolving_map = 映射端到端验证算例）
```

---

## 二、物理映射（Fluent 路线 A → OpenFOAM）

| Fluent（2023 R2） | OpenFOAM（multiphaseEulerFoam） | 说明 |
|---|---|---|
| Mixture 多相（2 相） | Eulerian 多相（gas + liquid） | 各相独立速度/温度，可重力驱动凝液排水 |
| 主相 humid-air（air+h2o） | `gas` 相，`multiComponentPhaseModel` | 组分 air（不凝气体）+ H2O（可凝） |
| 次相 water-liquid | `liquid` 相，`purePhaseModel`（纯水） | |
| Lee 蒸发-凝结（温度驱动） | `interfaceCompositionPhaseChange`（界面饱和 + Stefan 扩散） | 见下方「关键差异」 |
| wall_cool 定壁温 T_c | `coldWall` 定壁温（`fixedValue`） | 固/流界面简化为冷壁 |
| inlet 质量流量 + 含湿 | `U.gas` 定速 + `H2O.gas` 质量分数 | 流量→流速换算见参数表 |
| outlet 压力出口 | `p_rgh` 定压 + `inletOutlet` 防倒流 | |
| 重力 -Y（排水方向） | `constant/g = (0 -9.81 0)` | 已从 case 实测 `gravity/y=-9.81` |

**关键差异（务必理解）**：Fluent 的 Lee 模型是**温度驱动**的经验式
`ṁ = β·ρ·(T_sat−T)/T_sat`，β 是经验系数。OpenFOAM 这里用的是**界面组分相变**模型：
气液界面取饱和（ArdenBuck 饱和压力），水蒸气经**不凝气体（空气）扩散**到界面凝结
（Stefan 流）。后者是含不凝气体冷凝的**物理正确**做法，直接对应你 SOP 里强调的
「不凝气体扩散阻力」；它不需要标定 Lee 频率 β，但结果会受气相组分扩散与界面传热控制。

---

## 三、运行步骤

### 1. 安装 OpenFOAM（见《安装指南.md》）
本机实际使用 **blueCFD-Core OpenFOAM 8**（Windows 原生，`D:\APP\openfoam`）；WSL2 方案保留作备选。

### 2. 把网格转换进算例（v8_1 已完成，其余结构重复此步）
```bash
# 在本机 Windows（用工程内置 Python）
C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe \
    tools/fluent_cas_to_foam.py "<结构>.cas.h5" "case_qishi"
```
输出 `case_qishi/constant/polyMesh/`，并把边界区命名固定为 `inlet / outlet / coldWall / wall`。
不同结构的 `.cas.h5` 都能用同一个脚本（自动识别流体/固体 cell zone 与固流耦合界面）。

> ⚠️ **8 结构批量链路尚未重建**：`cases/` 目录与 `prepare_all_structures.py` 的批量流程已停用，
> 现在按单结构手工走两步共轭 SOP（见 §6）。批量化是待办 P1-3（先写 `generate_case_cht.py`）。

### 3. 校核网格

先用本机脚本做快速检查（无需进 WSL/OpenFOAM）：
```bash
C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe \
    tools/check_face_orientation.py case_qishi/constant/polyMesh
```
检查项：A owner<neighbour、B (owner,neighbour) 字典序、C 边界分区连续性、
D1/D2 面法向朝向、E 未使用节点。要求最后一行 `RESULT: CONFORMS`。

再用 OpenFOAM 自带的权威检查（blueCFD 环境）：
```bash
C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe \
    tools/of8.py case_qishi "checkMesh -constant"
```
要求输出 `Mesh OK`（`Failed 0 mesh checks`）。
注意：只检查 `checkMesh` 的「非上三角排序」是不够的——**同owner 的面之间 neighbour 也必须递增**，
否则会报 `Faces not in upper triangular order`（旧网格实测 127895 个）。

> 求解若在第 2 个 PIMPLE 外循环崩溃（`janafThermo::limit` FATAL / T 掉到 200 K 以下），
> 用消融探针定位根因，不要靠猜：
> ```bash
> python tools/ablate.py --list
> python tools/ablate.py baseline every_step --endtime 0.02
> ```
> `every_step` 会在崩溃前抢先把场统计落盘——求解器崩在写盘之前，这是唯一可靠的数据来源。

### 4. 校核字典（可选，确认字典能解析）
```bash
foamDictionary constant/phaseProperties    # 快速核对字典能解析
```

### 5. 运行
```bash
multiphaseEulerFoam > log 2>&1 &
tail -f log
```

### 6. 两步共轭全流程 SOP（当前正式路线，2026-10-20 定稿）

> 为什么两步：wall_cool（5℃ 面）是**固体外表面**而非流固界面，multiphaseEulerFoam 无法带固体；
> OF8 没有「多相相变 + 固体导热」一步求解器。故 Step1 用 chtMultiRegionFoam 算共轭传热得
> 真实界面温度分布 T_int(x)，Step2 把 T_int 作为 nonuniform 边界跑冷凝。
> **前提（2026-10-20 更新，已不再是阻塞）**：checkMesh 报的 `multiple in between faces`
> 经三天实测定性为**铰链 pinch**（同一对单元绕一条棱折起来的两个面瓣，面心距 339 µm、
> 面积差 11 倍、共用 2 个节点），是**合法拓扑**：每个单元都闭合（`Max cell openness` ~1e-16），
> 有限体积逐面独立积分，不影响守恒。**不要修**——删面会让单元开口，并单元会让缺陷扩散。
> 判定用 `tools/mesh_topology_report.py <polyMesh>`（它只会对"真重复面"动手）。
> 真要彻底消掉才回 CAD/网格器让两体面接触或留 ~0.05 mm 间隙。

```bash
PY=C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe
CASE=<修复后的>.cas.h5

# Step 0  开跑门禁（2 秒；FAIL 时退出码=1，先处理后开跑）
$PY tools/preflight.py caseCHT_qishi            # 语法/依赖/网格计数/拓扑/补丁名/controlDict/耦合
$PY tools/preflight.py caseCHT_qishi --checkmesh  # 追加 checkMesh（慢约 10 s）

# Step 0.5  双区域转换（fluid + solid 各自 polyMesh + mapped 耦合补丁）
$PY tools/fluent_cas_to_foam_cht.py "$CASE" caseCHT_qishi
$PY tools/of8.py caseCHT_qishi "checkMesh -region fluid"     # 要求 Mesh OK
$PY tools/of8.py caseCHT_qishi "checkMesh -region solid"

# Step 1  共轭传热（仅 coldWall=278.15K，其余外表面绝热）
$PY tools/of8.py caseCHT_qishi "chtMultiRegionFoam > log.cht 2>&1"
# 守恒核对：postProcessing/fluid/{m_in,m_out,T_in,T_out,Q_interface}
#   m_in ≈ m_out（ρConst），Q_interface·面积 ≈ m·cp·(T_out−T_in)

# Step 1.5  界面温度提取（取收敛的最后时刻场）
$PY tools/map_interface_T.py extract caseCHT_qishi/constant/fluid/polyMesh \
        <最后时刻>/fluid/T fluid_to_solid _interface_T.csv

# Step 1.6  同网格一致性校验（正式管线必须 PASS；跨几何会 FAIL 并给出偏差报告）
$PY tools/map_interface_T.py verify _interface_T.csv case_v81_1tau/constant/polyMesh coldWall

# Step 2  冷凝（映射为 coldWall nonuniform fixedValue 后运行）
$PY tools/align_step2_thermo.py case_v81_1tau  # 必须：液相 multiComponent+Henry / 气相 eConst
                                               # 否则 janaf 反解温度不收敛（Maximum number of iterations 100）
$PY tools/map_interface_T.py map _interface_T.csv case_v81_1tau/constant/polyMesh coldWall case_v81_1tau/0/T.gas
$PY tools/map_interface_T.py map _interface_T.csv case_v81_1tau/constant/polyMesh coldWall case_v81_1tau/0/T.liquid
$PY tools/of8.py case_v81_1tau "multiphaseEulerFoam > log 2>&1"

# Step 3  六项指标对比（mdot_cond/Y_H2O_out/T_out/Q_wall/dp 自动汇总）
$PY tools/compare.py
```

管线状态（2026-10-20 23:30）：压力口径已修正（`p_rgh` 绝压，preflight 第 11 项把关）；
v8_1 的 T_int 已用 9000 步重映射（9459 面、均值 278.279 K、verify PASS）。
当前 P0 = **液相温度发散**：1τ 在 t=0.0467 s 报 `Negative initial temperature T0`（液相先到 218–443 K）；
Step1 能量比 1.080（阈值 <1.05），需续算后再取最终 T_int。
旧结论"网格源头修复 / 薄片单元"已由 v8_8 审计与压力修复推翻，不要再按它排产。
先跑满 1τ（约 27 h / 8 进程），观察 `alpha.liquid > 0`、沿 -Y 下沉与出口液相。

---

## 四、参数在哪里改

都在 `tools/generate_case.py` 顶部（改完重跑脚本再运行算例）：

| 参数 | 默认值 | 对应 Fluent |
|---|---|---|
| `T_C` | 5.0 ℃ | wall_cool 定壁温 |
| `T_IN` | 35.0 ℃ | inlet 总温 |
| `MDOT` | 2e-5 kg/s | inlet 质量流量 |
| `Y_H2O` | 0.0289 | 入口水蒸气质量分数（80%RH） |
| `P_OP` | 101325 Pa | 操作压力 |
| `G` | (0, -9.81, 0) | 重力/排水方向 |

> 扫描 T_c 对比结构：只改 `T_C` 重跑 `generate_case.py`，再用不同 `OUT_DIR` 并行。

---

## 五、考核指标后处理（对应 Fluent Report Definition）

在算例目录里用 `postProcess`（或加入 `controlDict` 的 `functions`）：

```bash
# ① 出口凝液质量流率 ṁ_cond（液相经 outlet 排出）
postProcess -func "phasePatchFlowRate(name=outlet, phase=liquid)"

# ② 出口水蒸气质量分数（H2O.gas 的质量加权平均）
postProcess -func "patchExpression(name=outlet, expression='alpha.gas*thermo:rho.gas*H2O.gas/alpha.gas*thermo:rho.gas', accumulations='weightedAverage')"

# ③ 冷壁换热量 Q_wall（对 coldWall 面做热流积分）
postProcess -func "patchHeatFlux(name=coldWall)"   # 或 wallHeatFlux

# ④ 压降 Δp = p@inlet − p@outlet
postProcess -func "patchAverage(name=inlet,p)" -func "patchAverage(name=outlet,p)"
```

**守恒校核（必做，同 SOP §8）**：
- 质量：入口 H₂O 流量 ≈ 出口气相 H₂O 流量 + ṁ_cond（误差 < 2%）
- 能量：Q_wall ≈ ṁ_air·cp·(T_out−T_in) + ṁ_cond·h_fg

---

## 六、已知假设与局限（务必阅读）

1. **固/流共轭传热被简化**：`coldWall` 直接设 T_c，等价于把铝结构当作理想导体。
   **（2026-10-19 晚间修正）**该做法被否决：wall_cool（5℃ 面）是**固体外表面**（如 v8_12
   顶板约 16cm²），不是流固界面；省略固体后冷源与流体断开，且界面全部被人为压到 5℃
   （过度冷却）。改为**两步共轭方案**：
   - Step 1 `chtMultiRegionFoam`：流体（单相湿空气）+ 固体（铝，含 v8_12~14 的空腔），
     仅 wall_cool 定 278.15 K，稳态求得流固界面的温度分布 T_int(x)；
   - Step 2 `multiphaseEulerFoam`：把 T_int 映射为 coldWall 的 nonuniform fixedValue 后跑冷凝。
   局限：单向耦合（未计入凝液潜热对固体的反馈），必要时用 Step 2 的 Q_wall 迭代修正。
   转换器需扩展：导出固体 cell zone（cellZones + splitMeshRegions 分区）。
2. **界面组分相变 ≠ Lee 温度驱动**：两者机理不同，绝对凝液量会有差异；但结构**相对排序**
   受「接触面积 + 停留时间 + 流场组织」主导，结论方向一致。若一定要复刻 Lee，需改用
   `thermalPhaseChangeSystem`（温度驱动），可自行对照 `bubbleColumnEvaporatingDissolving` 教程。
3. **Eulerian 相模型**用气泡/液滴直径描述分散相，对「壁面膜状冷凝」是近似；
   若关注液膜形态，需转 VOF（`interCondensatingEvaporatingFoam`）或 wall-film 模型。
4. **入口流速方向**按 inlet 面法向设为 -Y（流入），运行前建议用 ParaView 看一次流场确认。
5. **初始无液相**（`alpha.liquid=0`），凝液完全由相变生成；稳态是否成立取决于凝液能否排走。

---

## 七、下一步建议（对应 SOP §5 的改进路线）

- **【已撤回，不必再做】** ~~v8_1 网格源头修复（修捏点）~~、~~给进出口面下薄片单元加密~~。
  捏点已定性为铰链拓扑（合法）；"薄片单元"已被 v8_8 审计用 `q=V/A^1.5` 与 `checkMesh` 反驳
  （两套网格都没有薄片，只有 67–96 µm 的等轴小胞）。真正的启动阻塞是压力口径（已修），
  当前的 1τ 阻塞是液相温度发散（P0）。
- 当前真正的进展项：Step1 已到 9000 步（能量比 1.080），T_int 已映射；
  Step2 待解决"液相能量方程在 α_l→1e-6 时的病态"（候选：T 限幅 / 液相 hConst / 共享温度场）。
- 备选扩展（不阻塞主线）：
  - 单向耦合的潜热反馈补偿：用 Step2 的 Q_wall 反代 Step1 冷壁热流迭代一次；
  - `tools/diagnostics/audit_sources.py` 升级为正式工具并接入 prepare_all_structures 前置检查；
  - 网格无关性（GCI）：冷壁与相变区加密，两套网格结果差 <1%。
  - 收敛判据统一：残差 + 指标量（ṁ_cond、出口 H₂O 分数）在 N 步内变化 < 阈值。
