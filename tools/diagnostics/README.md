# tools/diagnostics — 诊断脚本归档区

本目录存放 2026-10-19~20 排错期间使用、仍有复用价值的诊断脚本。
**它们不属于正式管线**（正式管线见 `../README.md`），使用前请阅读各文件头注释。

| 脚本 | 用途 | 状态 |
|---|---|---|
| `audit_sources.py` | 批量解析 8 个 FLTG*-Setup-Output.cas.h5 的模型/边界/材料设置，输出一致性审计表（原 `audit_sources.py（原 _audit_sources.py）`） | 可用；建议后续升级为正式 `tools/audit_sources.py` 并接入 prepare_all_structures 前置检查 |
| `dup_faces.py` | 检测 polyMesh 中「同一对单元被多个内部面连接」的非流形缺陷（原 `dup_faces.py（原 _dup_faces.py）`）。注意 numpy 高级索引副本坑：`mask[order][1:] \|= x` 是无效操作 | 可用；v8_1 实测发现面 113533/113534 双连接 |
| `liq_probe.py` | multiphaseEulerFoam 教程对齐配置生成器：`make_liq_multicomponent`（液相 multiComponent+Henry）+ `GAS_ECONST`（气相 eConst）。可消除 janaf FATAL（原 `liq_probe.py（原 _liq_probe.py）`） | 可用 |
| `pc_probe.py` | 相变隔离探针：`disable_phase_change` 置空 interfaceComposition（原 `pc_probe.py（原 _pc_probe.py）`）。正则用非贪婪匹配 `interfaceComposition.gas\n(...)`，勿改回被 `species ( H2O );` 内联截坏的形式 | 可用 |
| `num_probe.py` | 教程 fvSchemes/fvSolution 整套对照 + fieldMinMax 定位探针（原 `num_probe.py（原 _num_probe.py）`） | 可用 |

已删除的一次性探针（历史结论已沉淀到项目日志与技能库，无需保留代码）：
`_ablate*.py`（迭代版，终版在 `tools/ablate.py`）、`_catch/_fo/_probe_T/_runof/_xval.py`。

相关结论速查：
- janaf 物性正确（构造时自动乘 R 换算），勿再「修单位」；
- v8_1 网格缺陷 = 捏点（枢纽薄片 22767，14 面 14 邻居），OpenFOAM 侧不可修，需 Workbench/Fluent 源头修复；
- 修复后的正式流程 = `fluent_cas_to_foam_cht.py` → `caseCHT`(chtMultiRegionFoam) → `map_interface_T.py` → `case`(multiphaseEulerFoam)。
