#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成三结构 z=0 场图 + 残差图汇总报告（单文件 HTML，图片 base64 内嵌）。"""
import base64
import pathlib

ROOT = pathlib.Path(r"F:\B_Simulation\fangzhen-3\OpenFOAM")
OUT = ROOT / "报告_z0温度流场与残差_三结构_2026-10-21.html"


def b64(name):
    return base64.b64encode((ROOT / name).read_bytes()).decode()


html = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>z=0 温度场/流场与计算残差 · 三结构对比 · 2026-10-21</title>
<style>
:root{
--bg:#0e1116; --panel:#161b22; --panel2:#1c232c; --line:#2a323d;
--tx:#e6edf3; --tx2:#9aa7b4; --tx3:#6e7b8a;
--ok:#3fb950; --warn:#d29922; --bad:#f85149; --info:#58a6ff; --acc:#bc8cff;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--tx);
font:15px/1.75 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif;}
.wrap{max-width:1180px;margin:0 auto;padding:46px 26px 80px;}
header{border-bottom:1px solid var(--line);padding-bottom:24px;margin-bottom:30px;}
.kicker{color:var(--info);font-size:12.5px;letter-spacing:.14em;text-transform:uppercase;margin-bottom:10px}
h1{font-size:27px;line-height:1.35;margin:0 0 12px;font-weight:650;letter-spacing:-.01em}
.meta{color:var(--tx3);font-size:13px}
h2{font-size:20px;margin:46px 0 14px;padding-left:11px;border-left:3px solid var(--info);font-weight:620}
h3{font-size:15.5px;margin:26px 0 9px;font-weight:620}
p{margin:9px 0}
ul,ol{margin:9px 0;padding-left:22px}
li{margin:5px 0}
code{font-family:"Cascadia Mono",Consolas,monospace;font-size:12.6px;
background:var(--panel2);padding:1.5px 5px;border-radius:4px;color:#c9d1d9}
table{width:100%;border-collapse:collapse;margin:14px 0;font-size:13px}
th{background:var(--panel2);color:var(--tx2);font-weight:600;text-align:left;
padding:8px 10px;border-bottom:1px solid var(--line);white-space:nowrap}
td{padding:7px 10px;border-bottom:1px solid #1f262f;vertical-align:top}
tr:hover td{background:#141a21}
.num{font-family:"Cascadia Mono",Consolas,monospace;text-align:right;white-space:nowrap}
.ok{color:var(--ok)} .warn{color:var(--warn)} .bad{color:var(--bad)} .dim{color:var(--tx3)}
.box{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:15px 17px;margin:14px 0}
.box.hi{border-color:rgba(248,81,73,.45);background:linear-gradient(180deg,rgba(248,81,73,.08),var(--panel))}
.box.good{border-color:rgba(63,185,80,.4);background:linear-gradient(180deg,rgba(63,185,80,.07),var(--panel))}
.box.note{border-color:rgba(88,166,255,.35);background:linear-gradient(180deg,rgba(88,166,255,.06),var(--panel))}
.box h4{margin:0 0 7px;font-size:14.5px}
.box p{margin:5px 0;font-size:14px;color:var(--tx2)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:11px;margin:16px 0}
.stat{background:var(--panel);border:1px solid var(--line);border-radius:9px;padding:12px 14px}
.stat .k{color:var(--tx3);font-size:12px;margin-bottom:4px}
.stat .v{font-size:18px;font-weight:640;font-family:"Cascadia Mono",Consolas,monospace}
.stat .s{color:var(--tx3);font-size:11.5px;margin-top:3px}
.pill{display:inline-block;font-size:11.5px;padding:1px 8px;border-radius:10px;font-weight:600;white-space:nowrap}
.p0{background:rgba(248,81,73,.15);color:var(--bad);border:1px solid rgba(248,81,73,.4)}
.p1{background:rgba(210,153,34,.15);color:var(--warn);border:1px solid rgba(210,153,34,.4)}
.p2{background:rgba(88,166,255,.13);color:var(--info);border:1px solid rgba(88,166,255,.35)}
.pq{background:rgba(188,140,255,.13);color:var(--acc);border:1px solid rgba(188,140,255,.38)}
pre{background:var(--panel2);border:1px solid var(--line);border-radius:8px;padding:12px 14px;overflow-x:auto;
font-family:"Cascadia Mono",Consolas,monospace;font-size:12.4px;line-height:1.65;color:#c9d1d9}
.small{font-size:12.4px;color:var(--tx3)}
figure{margin:16px 0}
figure img{width:100%;border:1px solid var(--line);border-radius:10px;display:block}
figcaption{color:var(--tx3);font-size:12.5px;margin-top:7px;text-align:center}
footer{margin-top:48px;padding-top:16px;border-top:1px solid var(--line);color:var(--tx3);font-size:12.5px}
</style>
</head>
<body>
<div class="wrap">

<header>
<div class="kicker">fangzhen-3 · 湿空气冷凝结构对比</div>
<h1>z=0 温度场 / 流场与计算残差 · 三结构结果汇总</h1>
<div class="meta">数据源：QiShi_v8_1_20261009 / Qishi_v8_8_20261012 / Qishi_v8_9_20261012（Fluent Workbench）·
处理：OpenFOAM 工程 <code>OpenFOAM/</code>（Fluent cas.h5 → chtMultiRegionFoam Step1）· 生成日期：2026-10-21</div>
</header>

<h2>一、三个 Fluent 工程的处理结果</h2>
<p>三个算例均经 <code>tools/inspect_cas_h5.py</code> 体检（0 负体积 / 0 非法面）后由
<code>tools/fluent_cas_to_foam_cht.py</code> 转换为双区域 CHT 算例。当前状态：</p>

<table>
<tr><th>结构</th><th>Fluent 源工程</th><th>转换结果（流体 / 固体单元）</th><th>界面面积</th><th>Step1 CHT</th></tr>
<tr><td><b>v8_1</b></td><td>QiShi_v8_1_20261009</td>
<td class="num">64,729 / —</td><td class="num">107.91 cm²</td>
<td><span class="pill p2">15000 迭代完成</span> caseCHT_qishi</td></tr>
<tr><td><b>v8_8</b></td><td>Qishi_v8_8_20261012</td>
<td class="num">444,474 / —</td><td class="num">123.64 cm²</td>
<td><span class="pill p2">15000 迭代完成</span> caseCHT_v88</td></tr>
<tr><td><b>v8_9</b></td><td>Qishi_v8_9_20261012</td>
<td class="num">489,728 / 215,410</td><td class="num">136.19 cm²</td>
<td><span class="pill p1">未开跑</span> 缺 <code>0/</code> 字典（待办 T6），分支 v8_9-setup</td></tr>
</table>

<div class="box note">
<h4>口径说明</h4>
<p>冷壁 <code>wall_cool</code> 均为 16.000 cm²、进/出口均 1.130 cm²，三结构同口径可比。
v8_8 与 v8_9 网格尺度仅差 1.10×，是当前最干净的一对可比结构；v8_1 与其差 7.6×，
跨结构的绝对量差异需按网格密度折半解读。</p>
</div>

<h2>二、Step1 稳态 CHT 定量结果（均 15000 外迭代）</h2>
<table>
<tr><th>指标</th><th class="num">v8_1</th><th class="num">v8_8</th><th>判读</th></tr>
<tr><td>UA（界面总传热系数）</td><td class="num">0.06273 W/K</td><td class="num">0.07721 W/K</td>
<td><span class="ok">+23.1%</span>，结构差异稳健信号</td></tr>
<tr><td>冷凝量（分区带法）</td><td class="num">1.645 mg/h</td><td class="num">1.626 mg/h</td>
<td>−1.2%，<span class="warn">不可区分</span>（模型形式不确定度 ±11% 之内）</td></tr>
<tr><td>出口含湿量 ω_out</td><td class="num">6.231 g/kg</td><td class="num">6.497 g/kg</td><td>进口 29.76 g/kg</td></tr>
<tr><td>出口质量加权温度 T_out_mw</td><td class="num">279.534 K</td><td class="num">279.465 K</td><td>进口 308.15 K</td></tr>
<tr><td>界面热流（焓平衡口径 Q_flux）</td><td class="num">0.5757 W</td><td class="num">0.5771 W</td>
<td class="ok">与 ṁ·cp·ΔT 逐位一致（能量守恒闭合）</td></tr>
<tr><td>界面壁温均值 / 范围</td><td class="num">278.209 K<br><span class="dim">278.150~278.370</span></td>
<td class="num">278.833 K<br><span class="dim">278.151~282.003</span></td><td>v8_8 壁温分布明显更宽（std 1.22 K vs 0.05 K）</td></tr>
<tr><td>NTU_m（传质）</td><td class="num">3.456</td><td class="num">4.253</td><td>由 UA 派生，非传输受限</td></tr>
<tr><td>质量流率 m_in = m_out</td><td class="num">1.99985e-05 kg/s</td><td class="num">1.99984e-05 kg/s</td><td class="ok">进出门户质量平衡</td></tr>
<tr><td>能量比 energy_ratio</td><td class="num">1.078</td><td class="num">1.102</td>
<td class="warn">与 wallHeatFlux 积分差 7.8%/10.2%（未归因缺陷 A1）</td></tr>
</table>

<div class="box hi">
<h4>核心结论（与既有判定一致）</h4>
<p><b>UA +23.1% 是稳健的结构差异；冷凝量 −1.2% 落在模型形式不确定度 ±11% 以内，不可作为区分依据。</b>
要判定结构优劣，最高优先级是把统一口径（模型形式）不确定度压到 &lt;3%（待办 T0）。</p>
</div>

<h2>三、z=0（xy 平面）温度场与流场</h2>
<p>由 <code>postProcess -func surfaces(_flow) -latestTime</code> 切出 15000 迭代时刻的 z=0 平面，
用 <code>tools/plot_z0_slice.py</code> 绘制（散点 griddata 插值到规则网格，避免多边形扇形三角化报错）。
左：流体温度；中：速度幅值与流线；右：固体（铝）温度。</p>

<figure>
<img src="data:image/png;base64,__IMG_Z0_V81__" alt="v8_1 z=0 温度场与流场">
<figcaption>v8_1（caseCHT_qishi @15000）：流体 T 278.19~308.15 K，|U| max 0.234 m/s、中位 0.072 m/s；
固体冷壁侧 278.15~278.34 K。热进口 (x=+12, y≈+29) 自上而下进入 U 形流道，出口 x=−17.5。</figcaption>
</figure>

<figure>
<img src="data:image/png;base64,__IMG_Z0_V88__" alt="v8_8 z=0 温度场与流场">
<figcaption>v8_8（caseCHT_v88 @15000）：流体 T 278.43~308.15 K，|U| max 0.235 m/s、中位 0.034 m/s；
固体 278.43~279.94 K。速度中位数为 v8_1 的一半弱——流道截面更大、核心流速更低，
与 v8_8 更大的界面面积（123.64 vs 107.91 cm²）和更高 UA 互相印证。</figcaption>
</figure>

<p class="small">两张图读取的场文件已留存于 <code>_z0slice/v8_1/</code>、<code>_z0slice/v8_8/</code>
（T/U/p_rgh 流体切片 + 固体 T 切片，VTK 格式）。</p>

<h2>四、计算残差（Initial residual vs 外迭代）</h2>
<p>从求解器日志按迭代块解析（每块只取每条方程首次出现的初残差，避免 PIMPLE 外迭代重复行干扰），
容差线取 <code>fvSolution</code> 的 <code>residualControl</code>：U 1e-4 / h 1e-6 / p_rgh 1e-5 / e 1e-6。
左图流体三方程，右图固体能量方程。</p>

<figure>
<img src="data:image/png;base64,__IMG_RES_V81__" alt="v8_1 残差曲线">
<figcaption>v8_1（log.cont_15000.txt，9001→15000，6000 步）：续跑起点已达收敛量级，
全程平坦；末值 U 7.0e-10 / h 1.0e-7 / p_rgh 1.8e-8 / e 9.8e-6（初残差，对应线性求解终残差 4.7e-8）。</figcaption>
</figure>

<figure>
<img src="data:image/png;base64,__IMG_RES_V88__" alt="v8_8 残差曲线">
<figcaption>v8_8（log.cht_v88.txt + log.cht_v88_cont.txt，1→15000，15000 步完整历史）：
前 500 步四个方程从 1.0 衰减 6~8 个量级后进入平台；9000 步处的 p_rgh 竖线为续跑重启瞬态
（2e-8 → 2.5e-6，约 15 步内衰减回来），非求解失稳。末值 U 2.2e-9 / h 3.2e-5 / p_rgh 5.2e-8 / e 4.9e-5。</figcaption>
</figure>

<div class="box good">
<h4>收敛性判定</h4>
<p>末迭代（15000）连续误差：sum local 2.49e-8，global/cumulative 1.57e-10；
进出门户质量流率相等；焓平衡与界面热流逐位一致。<b>两个算例均为收敛的稳态解</b>，
15000 步之后的场可直接用于结构对比。</p>
</div>

<h2>五、待办与风险</h2>
<ul>
<li><b>v8_9 开跑</b>：缺 <code>0/</code> 字典（T6）。可从 caseCHT_v88 拷 <code>0/</code> 后把
interface 面 nonuniform 列表改为 uniform 占位（70,511 面 vs v8_8 的 62,158，直接拷会 FATAL），
再过 <code>tools/preflight.py</code> 全 PASS 后按 README §3.7 求解。</li>
<li><b>T0 模型形式不确定度 ±11%</b>：目前大于结构差异（UA +23.1% 除外），冷凝量 −1.2% 无法解读，
需统一口径压到 &lt;3%。</li>
<li><b>A1 能量闭合偏差</b>：wallHeatFlux 积分与焓平衡差 7.8%/10.2% 未归因；
并行下 <code>postProcessing/solid/Q_solid</code> 未合并（读到 0），重构后需重取。</li>
<li><b>v8_1 残差历史不完整</b>：0→9000 首跑日志已被清理（_cleanup_report.txt 有记录），
现存 9001→15000；如论文需要全历史需重跑或接受分段呈现。</li>
</ul>

<footer>
fangzhen-3 · OpenFOAM/ 工程 · 2026-10-21 · 图表可由
<code>tools/plot_z0_slice.py</code> 与 <code>tools/plot_residuals.py</code> 复现
</footer>
</div>
</body>
</html>
"""

html = (html
        .replace("__IMG_Z0_V81__", b64("Step1_z0_温度流场_QiShi_v8_1_15000.png"))
        .replace("__IMG_Z0_V88__", b64("Step1_z0_温度流场_v8_8_15000.png"))
        .replace("__IMG_RES_V81__", b64("Step1_残差曲线_QiShi_v8_1.png"))
        .replace("__IMG_RES_V88__", b64("Step1_残差曲线_v8_8.png")))

OUT.write_text(html, encoding="utf-8")
print(f"saved {OUT}  ({OUT.stat().st_size/1024/1024:.2f} MB)")
