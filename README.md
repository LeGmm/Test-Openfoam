# 湿空气冷凝仿真 —— Fluent → OpenFOAM 移植

把 `routeA_condensation_setup.py`（Fluent 路线 A：Mixture 多相 + 组分输运 + Lee 相变）移植到
**OpenFOAM（multiphaseEulerFoam + 界面组分相变）**，用于多结构「冷凝效率」对比实验。

> 当前主线是**稳态结构对比**（见 §3.7）：v8_1 `caseCHT_qishi/`（15000 外迭代）与 v8_8 `caseCHT_v88/`
> （15000 外迭代，2026-10-21 续算至收敛），用 `tools/condensation_metrics.py` 出四项对比指标。
> 修正后结论（**分区集总 + 活塞流口径**）：v8_1 冷凝量 1.645 mg/h、v8_8 1.626 mg/h（−1.2%）——
> **两者在现有精度下不可区分**（模型形式不确定度 ±11% 已大于结构差异 1.2%）；UA 差 +23.1% 是稳健结论。
> 详见《报告_修复与v8_9体检_2026-10-21.html》；《报告_三任务交付_2026-10-21.html》§3 的冷凝量/UA/NTU_m
> /出口含湿量数字**已作废**（旧报告保留了作废横幅，供追溯）。
> 已废弃的 Step2 瞬态路线（`case_v81_1tau/`）保留作 P0-1 诊断用，但**不再是出数的前提**。
> 旧的 `case/`、`caseCHT/` 已于 2026-10-20 归档移出，回收站本体已于 2026-10-21 删除。
> 版本控制：2026-10-21 起本目录是 git 仓库，远端
> `https://github.com/LeGmm/Test-Openfoam`（PRIVATE，分支 `master`，tag `baseline-2026-10-21`；
> `.gitignore` 排除网格/结果/日志/缓存，仅入库配置与脚本）。

---

## 一、目录结构

```
OpenFOAM/
├── README.md                  # 本文件
├── 安装指南.md                 # ⚠ 方案 A 写的是 WSL2 + OpenFOAM.org v11，本机实际用 blueCFD-Core 8
├── .gitignore / .gitattributes
├── tools/
│   ├── fluent_cas_to_foam.py  # Fluent .cas.h5 → 单流体区 polyMesh（Step2 用，已实测验证）
│   ├── fluent_cas_to_foam_cht.py # Fluent .cas.h5 → fluid+solid 双区域 polyMesh（Step1 用）
│   ├── inspect_cas_h5.py      # Fluent .cas.h5 体检器：完整性/规模/区域名/负体积/非法单元
│   ├── setup_cht_coupling.py  # 两侧界面统一改名 interface + T 写 coupledBaffleMixed（Step1 必做）
│   ├── align_step2_thermo.py  # Step2 热物性对齐（液相 multiComponent+Henry、气相 eConst）※生成算例后必跑
│   ├── warmstart_step2.py     # Step2 初场/warmstart（拷 p_rgh +101325 等绝压口径处理）
│   ├── ofmesh.py              # polyMesh 读写公共模块（新脚本统一从这里 import）
│   ├── preflight.py           # 开跑前门禁：语法/依赖/网格/拓扑/补丁名/controlDict/耦合/热物性/压力
│   ├── map_interface_T.py     # Step1→Step2 界面温度映射（extract/extract-vtk/map/verify）
│   ├── condensation_metrics.py # ★ 稳态结构对比的四项指标（含分区集总冷凝模型，见 §3.7）
│   ├── mesh_topology_report.py # 判"一对单元多张面"是【真重复】还是【铰链】；--fix 只删真重复
│   ├── audit_cleanup.py       # 清理候选三项核对器（只读，绝不自删）
│   ├── _do_cleanup.py         # 按 audit_cleanup 的机读结果执行删除（一次性）
│   ├── merge_dup_cells.py     # ⚠️ 已弃用：迭代并单元在 polyhedral 网格上会级联毁网格（64729→28378）
│   ├── fix_pinched_cells.py   # ⚠️ 已弃用：单对并细胞会让缺陷扩散（1 对 → 7 对 → …）
│   ├── drop_duplicate_faces.py # ⚠️ 已弃用：删面会让铰链单元开口（open cells 2）
│   ├── plot_z0_slice.py       # 读 surfaces FO 的切片 VTK，画 z=0 温度/流场（scipy griddata）
│   ├── fix_zero_volume_cells.py  # 零体积单元节点注入工具
│   ├── check_face_orientation.py # 网格朝向/排序自检（纯 Python，秒级）
│   ├── of8.py / run_of_cmd.py # 在本机 blueCFD(OpenFOAM 8) 环境执行 OF 命令（后者实时写日志，适合长跑）
│   ├── generate_case.py       # 参数化生成全套 0/ constant/ system/ 字典
│   ├── prepare_all_structures.py # 批量转换 8 个结构（⚠ 依赖不存在的 cases/，当前停用）
│   ├── ablate.py              # 单因素对照探针（求解发散/崩溃归因）
│   ├── mesh_stats.py / compare.py # ⚠ 均硬编码依赖 cases/（不存在），当前停用
│   ├── fetch_installer.py / _diag_field_probe.py # 安装器抓取 / 场极值单元定位（排错用）
│   └── diagnostics/           # 排错期诊断脚本归档（audit_sources/dup_faces/liq_probe 等，见其 README）
├── caseCHT_qishi/             # Step1：v8_1 共轭传热（已跑到 15000 外迭代，收敛）
├── caseCHT_v88/               # Step1：v8_8 共轭传热（已跑到 15000 外迭代，收敛）
├── case_v81_1tau/             # Step2：v8_1 冷凝 1τ（T_int 已映射；崩溃诊断场 0.046 留存）
├── case_v88/                  # Step2：v8_8（字典仍是启动测试版，未映射、未跑）
├── case_qishi/                # 旧单区域算例（保留 0/ 与 4e-06 爆温证据）
└── tests/                     # 管线验证留存（bcDissolving_map = 映射端到端验证算例）
```

根目录还有若干报告与记录文件：`报告_修复与v8_9体检_2026-10-21.html`（★最新结论）、
`报告_三任务交付_2026-10-21.html`（§3 已作废）、`待办清单_2026-10-21.md`、`清理记录_2026-10-21.md`、
`代码审计报告_2026-10-20.html`、`网格修复指南.md`。

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

> ⚠️ **8 结构批量链路尚未重建**：`cases/` 目录与 `prepare_all_structures.py` / `compare.py` /
> `mesh_stats.py` / `run_all.sh` 的批量流程已停用（前三者硬编码依赖不存在的 `cases/`），
> 现在按单结构手工走两步共轭 SOP（见 §6）。批量化对应待办 T6（先写 `generate_case_cht.py`）。
> 另：转换完立刻跑 `checkMesh` 时若缺 `system/controlDict` 会 FATAL——先确保字典齐备再体检。

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

# Step 3  指标提取（稳态结构对比用 §3.7 的 condensation_metrics.py；
#         ⚠ tools/compare.py 依赖不存在的 cases/，已停用，不要跑）
$PY tools/condensation_metrics.py caseCHT_qishi 15000 --json _m_v81.json
```

管线状态（2026-10-21 更新）：压力口径已修正（`p_rgh` 绝压，preflight 第 11 项把关）；
v8_1、v8_8 的 Step1 均已跑到 **15000 外迭代并收敛**（出口漂移 v8_1 ≈1.5e-6 K/步）。
Step2 瞬态（1τ）已**降级为可选验证项**：唯一残留阻塞是液相温度发散
（t=0.0467 s 报 `Negative initial temperature T0`，液相先到 218–443 K），
当前主线不依赖它出数——结构对比走 §3.7 的稳态路线。
旧结论"网格源头修复 / 薄片单元"已由 v8_8 审计与压力修复推翻，不要再按它排产。

---

### 7. 结构对比：稳态路线（2026-10-21 定稿，当前主线）

> 为什么不跑瞬态：Step2（`multiphaseEulerFoam`）实测 0.949 s/步、dt≈2.93e-5，
> 跑满 1τ(2.97 s) ≈ **27 h/结构**，且当前崩在 t=0.0467 s（极小 α_l 单元液相温度发散）。
> 本机负荷不允许，因此结构对比按「稳态为主」改走下面这条链路。

```bash
PY=C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe

# ① 双区域转换（两结构各一次，已有产物可跳过）
$PY tools/fluent_cas_to_foam_cht.py "<结构>.cas.h5" caseCHT_<结构>
# ② 分区 + 稳态共轭 + 重构（8 进程；v8_8 实测 15000 外迭代 ≈ 5300 s，0.58 s/外迭代）
$PY tools/run_of_cmd.py caseCHT_X "decomposePar -allRegions -force"  log.decomp.txt
$PY tools/run_of_cmd.py caseCHT_X "@_cmd_cht_v88.txt"                log.cht.txt
$PY tools/run_of_cmd.py caseCHT_X "reconstructPar -allRegions"       log.recon.txt
# ③ 四项指标（壁面温度 / 换热系数 / 冷凝量 / 出口含湿量）
$PY tools/condensation_metrics.py caseCHT_X <最后外迭代数> --json _m_X.json
```

指标口径要点（**对比时务必保持一致**，2026-10-21 复核后定稿）：
- 热负荷用**焓平衡** `ṁ·cp·ΔT`（`Q_flux_W`），**不要**用 `wallHeatFlux` 场积分（`Q_W`）——
  后者在现网格上与焓平衡差 7.6%（v8_1）/ 14.4%（v8_8），两者之比存于 `energy_ratio`。
  已排除「未收敛残差」这一解释：流体热容 C=0.061 J/K，残差 0.0449 W 若属暂态应给出
  0.73 K/步，而实测出口漂移仅 1.5e-6 K/步，差 5×10⁵ 倍。
- 总热导 `UA = Q_flux/ΔT_lm`；**所有**传热与传质量都由 UA 派生：
  `h = UA/A`、`h_m = h/(ρ·cp·Le^{2/3})`、`NTU_m = UA/(cp·ṁ·Le^{2/3})`（代码里有自洽断言）。
  逐面 `h_global` 只在 q<0 且 ΔT>0.5 K 的面定义，仅作空间分布参考。
- 冷凝量主口径 = **分区集总**（按壁温等面积分带 + 活塞流对数平均推动力），见 `zoned` 字段。
  单一平均壁温对壁温分布宽的结构会失效：v8_8 壁温中位 283.6 K 高于干工况出口 280.4 K，
  取中位算 ΔT_lm 得到负温差、数学上无定义。
- **不再使用「出口取平衡饱和」**：分区模型实测出口仅 38% 饱和，该假设不成立。
  `m_cond_cstr_*`（全混口径，下界）与活塞流之差即模型形式不确定度（两结构上约 ±11%，
  已大于 v8_1/v8_8 之间 1.2% 的结构差异——这是当前对比的主要限制）。
- 假设与适用条件见 `tools/condensation_metrics.py` 文件头与
  《报告_三任务交付_2026-10-21.html》§3.6–3.7。
- 换结构进管线前先体检源文件：`python tools/inspect_cas_h5.py "<结构>.cas.h5"`
  （需 h5py；检查完整性、规模、区域名、负体积与非法单元）。

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
| `A_IN` | 1.1298e-4 m² | 入口面积（**硬编码，无校验**；换结构必须核对，见待办 T12） |

> 扫描 T_c 对比结构：只改 `T_C` 重跑 `generate_case.py`，再用不同 `OUT_DIR` 并行。

---

## 五、考核指标后处理（OF8 可用命令，2026-10-21 实测校订）

`postProcess -func X` 是按名字在 `$FOAM_ETC/caseDicts/postProcessing/` 找模板文件。
blueCFD OF8 自带的模板实测只有：`surfaceFieldValue/{patch,patchAverage,patchIntegrate}`、
`flowRate/flowRatePatch`、`fields/wallHeatFlux`、`pressure/pressureDifferencePatch`、
`minMax/*`、`numerical/residuals`、`visualization/surfaces` 等。
**旧版 README 里的 `phasePatchFlowRate`、`patchExpression`、`patchHeatFlux`、`surfaceFieldValue`
四个名字在本机 OF8 中均不存在**（实测报 `Cannot find functionObject file`），
它们是 OpenFOAM.com / 更新版本的接口。可用写法如下：

```bash
# ① 出口质量流率（flowRatePatch 只给混合物通量）
postProcess -latestTime -func "flowRatePatch(name=outlet)"

# ② 出口 H2O.gas 质量分数（面积平均；质量加权见 controlDict 的 weightedAverage FO）
postProcess -latestTime -func "patchAverage(name=outlet, H2O.gas)"

# ③ 冷壁热流：wallHeatFlux 只在【有 compressible momentumTransport 的区】可用
postProcess -latestTime -func wallHeatFlux        # Step1 fluid 区可用；Step2 多相会 FATAL

# ④ 压差（patchAverage 逐个取，或 pressureDifferencePatch 一次算）
postProcess -latestTime -func "patchAverage(name=inlet,p)"
postProcess -latestTime -func "patchAverage(name=outlet,p)"
```

> ⚠️ **实测根因（2026-10-21）**：`wallHeatFlux` 在 Step2 会报
> `Unable to find compressible turbulence model in the database`（wallHeatFlux.C:226）——
> 它要查 compressible momentumTransport，Step1（单相可压缩）有、Step2（multiphaseEulerFoam）没有。
> 所以 **Q_wall 只能从 Step1 拿，Step2 的能量闭合只能用气相焓平衡间接核对**。
>
> ⚠️ **逐相质量流率没有现成 FO**：凝液贴壁、当前模型无排液/液膜，`alphaRhoPhi.liquid`
> 在 `coldWall` 是 `zeroGradient + noSlip`，即使取到也是 0，**读出 0 不代表没冷凝**。
> 出口含湿量请用 §3.7 的 `condensation_metrics.py`（气相收支法）。
> 逐相流率的可行写法是像 `caseCHT_*/system/controlDict` 那样用
> `surfaceFieldValue` + `alphaRhoPhi.liquid`（solver 内置类，不受上面模板目录限制）。
> 另：④ 里 Step1 的压降用 `p`，Step2 用 `p_rgh`，两者口径不同。

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
- 当前真正的进展项：Step1 v8_1/v8_8 均已到 15000 步（能量比 1.078/1.102，口径为热流积分/焓平衡）；
  Step2 待解决"液相能量方程在 α_l→1e-6 时的病态"（候选：T 限幅 / 液相 hConst / 共享温度场）。
- 备选扩展（不阻塞主线）：
  - 单向耦合的潜热反馈补偿：用 Step2 的 Q_wall 反代 Step1 冷壁热流迭代一次；
  - `tools/diagnostics/audit_sources.py` 升级为正式工具并接入 prepare_all_structures 前置检查；
  - 网格无关性（GCI）：冷壁与相变区加密，两套网格结果差 <1%。
  - 收敛判据统一：残差 + 指标量（ṁ_cond、出口 H₂O 分数）在 N 步内变化 < 阈值。

## 八、版本控制与工具环境

- 仓库根 = 本目录，远端 `https://github.com/LeGmm/Test-Openfoam`（PRIVATE），分支 `master`，
  tag `baseline-2026-10-21`（清理前基线）。提交信息三段式：改动内容 / 影响范围 / 验证方式。
- 推送前先 `gh auth setup-git`（把 gh 的凭据接给 git）；已验证 `git push -u origin master` 可用。
- Python：脚本统一用 `C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe`
  （该 venv 已装 numpy 2.5.2 + h5py 3.16.0）。用不带 h5py 的解释器跑 `preflight.py`
  只会让 `deps` 一项 FAIL，不影响求解。
- 清理走 `tools/audit_cleanup.py`（只读三项核对）→ `tools/_do_cleanup.py`（执行），删除前先打 tag。
  已删 5.98 GB，工程 7.4→2.0 GB。
- 入库内容自检：`git ls-files` 共 272 个文件 / 20 MB，polyMesh、`processor*/`、
  `postProcessing*/`、`*.h5`、`log.*`、`__pycache__`、时间目录追踪数均为 0。

---
