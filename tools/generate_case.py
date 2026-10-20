#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 湿空气冷凝 OpenFOAM 算例生成器（multiphaseEulerFoam + 界面组分相变）
================================================================================
 物理对应（Fluent 路线 A -> OpenFOAM）：
   * Mixture 多相 + 组分输运 + Lee 相变
     -> multiphaseEulerFoam + interfaceCompositionPhaseChangeMultiphaseSystem
   * 气相 = 湿空气（air 不凝气体 + H2O 水蒸气，多组分）
   * 液相 = 纯水（water）
   * 相变驱动 = 界面饱和（ArdenBuck 饱和压力）+ 组分扩散（Stefan 流）
     —— 这是含不凝气体冷凝的物理正确做法（Fluent 的 Lee 模型是其温度驱动近似）

 运行（用工程内置 Python）：
   C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe \
       tools/generate_case.py [输出算例目录，默认 ../case]

 生成后仍需（见 README.md）：
   1. 用 tools/fluent_cas_to_foam.py 生成 constant/polyMesh（本脚本不生成网格）
   2. 在 OpenFOAM 里跑 foamDictionary 校核、checkMesh、setFields（如需初始积液）
================================================================================
"""

from __future__ import annotations

import os
import sys

# ============================================================================
# ① 用户参数 —— 与 routeA_condensation_setup.py 一一对应
# ============================================================================
T_C   = 5.0          # wall_cool 定壁温 = TEC 冷面温度 [℃]
T_IN  = 35.0         # 入口温度 [℃]
MDOT  = 2.0e-5       # 入口质量流量 [kg/s]
Y_H2O = 0.0289       # 入口水蒸气质量分数（35℃ / 80%RH）
T_SAT = 31.4         # 露点（仅作参考；实际饱和温度由 ArdenBuck 按当地分压给出）
P_OP  = 101325.0     # 操作压力 [Pa]

RHO_IN = 1.145       # 35℃ 湿空气密度 [kg/m³]（用于流量->流速换算）
A_IN   = 1.1298e-4   # 入口面积 [m²]（由网格换算，见 fluent_cas_to_foam 输出）
U_IN   = MDOT / (RHO_IN * A_IN)   # 入口流速 [m/s]

# 液相种子体积分数：不能让 alpha.liquid 初始为 0（会让液相方程奇异、产生 NaN）。
# 10^-6 相当于每立方米空气里有极微量的雾滴，对凝液总量影响可忽略，
# 但足以保证界面面积与液相方程非奇异。
ALPHA_LIQ_SEED = 1.0e-6

G      = (0.0, -9.81, 0.0)        # 重力：-Y = 排水方向（与 Fluent gravity/y=-9.81 一致）

# 边界分区名（由 fluent_cas_to_foam.py 生成，勿改）
P_INLET    = "inlet"
P_OUTLET   = "outlet"
P_COLD     = "coldWall"   # 固/流界面 = 冷壁（定壁温 T_C）
P_WALL     = "wall"       # 其余绝热壁

# 相名（与 phaseProperties 一致）
GAS   = "gas"
LIQ   = "liquid"


def c2k(t_c: float) -> float:
    return t_c + 273.15


def header(cls: str, obj: str, loc: str) -> str:
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | Version:  8                                     |\n"
        "|   \\\\  /    A nd           | Website:  https://openfoam.org                  |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       " + cls + ";\n"
        "    location    \"" + loc + "\";\n"
        "    object      " + obj + ";\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    )


def write(path: str, content: str) -> None:
    d = os.path.dirname(path)
    os.makedirs(d, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  wrote {os.path.relpath(path, start=case_dir)}")


# ============================================================================
# 各字典内容
# ============================================================================

def controlDict() -> str:
    return header("dictionary", "controlDict", "system") + f"""
application     multiphaseEulerFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         20;          // 先瞬态跑 20 s 观察凝液形成；稳态难收敛则继续加长
deltaT          0.001;       // 小时间步（相变刚性强）

writeControl    adjustableRunTime;
writeInterval   0.5;
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
adjustTimeStep  yes;         // 用 Courant 限制自动调步长
maxCo           0.3;
maxDeltaT       0.01;

functions
{{
    // ---- ① 出口凝液质量流率 ṁ_cond  kg/s ----
    // alphaRhoPhi.liquid = 液相【质量】通量 (alpha*rho*phi)，对 outlet 求和直接就是 kg/s。
    // （注意不可用 alphaPhi.liquid —— 那是体积通量，量纲是 m³/s）
    mdot_cond
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        writeFields     no;
        regionType      patch;
        name            outlet;
        operation       sum;
        fields          ( alphaRhoPhi.liquid );
    }}

    // ---- ② 出口水蒸气质量分数（面积平均，气相侧）----
    YH2O_out
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        writeFields     no;
        regionType      patch;
        name            outlet;
        operation       areaAverage;
        fields          ( H2O.gas );
    }}

    // ---- ③ 出口温度（用于能量平衡反算 Q_wall）----
    T_out
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        writeFields     no;
        regionType      patch;
        name            outlet;
        operation       areaAverage;
        fields          ( T.gas );
    }}

    // ---- ④ 压降：入口 / 出口 面积平均压力 ----
    p_in
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        writeFields     no;
        regionType      patch;
        name            inlet;
        operation       areaAverage;
        fields          ( p );
    }}
    p_out
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        writeFields     no;
        regionType      patch;
        name            outlet;
        operation       areaAverage;
        fields          ( p );
    }}
}}
"""


def fvSchemes() -> str:
    return header("dictionary", "fvSchemes", "system") + """
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    // 兜底格式：upwind 稳定。真实项按 OpenFOAM 8 教程逐条给出。
    default                         Gauss upwind;

    // 注意：引号内的 key 是正则，括号必须转义为 \\( \\)，否则匹配不到真实项名
    // （如 "div(phi,alpha.gas)"），会全部回落到 default 导致解析失败。
    "div\\(phi,alpha.*\\)"           Gauss vanLeer;
    "div\\(phir,alpha.*\\)"          Gauss vanLeer;

    "div\\(alphaRhoPhi.*,U.*\\)"     Gauss limitedLinearV 1;
    "div\\(phi.*,U.*\\)"             Gauss limitedLinearV 1;

    "div\\(alphaRhoPhi.*,Yi\\)"      Gauss limitedLinear 1;
    "div\\(alphaRhoPhi.*,(h|e).*\\)" Gauss limitedLinear 1;
    "div\\(alphaRhoPhi.*,K.*\\)"     Gauss limitedLinear 1;
    "div\\(alphaPhi.*,p\\)"          Gauss limitedLinear 1;

    "div\\(\\(\\(\\(alpha.*\\*thermo:rho.*\\)*nuEff.*\\)\\*dev2\\(T\\(grad\\(U.*\\)\\)\\)\\)\\)" Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear limited corrected 0.5;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         limited corrected 0.5;
}
"""


def fvSolution() -> str:
    return header("dictionary", "fvSolution", "system") + """
solvers
{
    "alpha.*"
    {
        nAlphaCorr      2;
        nAlphaSubCycles 1;
        cAlpha          1;
    }

    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.01;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|T).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }

    "(h|e).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }

    "Y.*"
    {
        // 组分输运由 MULES 限幅求解，必须给出 residualAlpha 与 cAlpha
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
        residualAlpha   1e-6;
        cAlpha          1;
    }
}

PIMPLE
{
    nOuterCorrectors    2;
    nCorrectors         2;
    nNonOrthogonalCorrectors 1;
    correctPhi          yes;
    ddtCorr             yes;
}

relaxationFactors
{
    equations
    {
        ".*"            0.7;
    }
}
"""


def gravity() -> str:
    return header("uniformDimensionedVectorField", "g", "constant") + f"""
dimensions      [0 1 -2 0 0 0 0];
value           ({G[0]} {G[1]} {G[2]});
"""


def momentumTransport(phase: str) -> str:
    return header("dictionary", f"momentumTransport.{phase}", "constant") + f"""
simulationType  laminar;
"""


def phaseProperties() -> str:
    # OpenFOAM 8 语法（BlueCFD-Core 2020）：相间用 (A in B) / (A and B) 写法，逐相 heatTransfer 等
    return header("dictionary", "phaseProperties", "constant") + f"""
type    interfaceCompositionPhaseChangeMultiphaseSystem;

phases (gas liquid);

gas
{{
    type            multiComponentPhaseModel;
    diameterModel   isothermal;
    isothermalCoeffs
    {{
        d0              3e-3;
        p0              1e5;
    }}

    residualAlpha   1e-6;
}}

liquid
{{
    type            purePhaseModel;
    diameterModel   constant;
    constantCoeffs
    {{
        d               1e-4;
    }}

    residualAlpha   1e-6;
}}

blending
{{
    default
    {{
        type            linear;
        minPartlyContinuousAlpha.gas 0.3;
        minFullyContinuousAlpha.gas 0.7;
        minPartlyContinuousAlpha.liquid 0.3;
        minFullyContinuousAlpha.liquid 0.7;
    }}

    heatTransfer
    {{
        type            linear;
        minPartlyContinuousAlpha.gas 0;
        minFullyContinuousAlpha.gas 1;
        minPartlyContinuousAlpha.liquid 0;
        minFullyContinuousAlpha.liquid 1;
    }}

    diffusiveMassTransfer
    {{
        $heatTransfer;
    }}
}}

surfaceTension
(
    (gas and liquid)
    {{
        type            constant;
        sigma           0.07;
    }}
);

aspectRatio
(
    (gas in liquid)
    {{
        type            constant;
        E0              1.0;
    }}

    (liquid in gas)
    {{
        type            constant;
        E0              1.0;
    }}
);

drag
(
    (gas in liquid)
    {{
        type            SchillerNaumann;
        residualRe      1e-3;
        swarmCorrection
        {{
            type        none;
        }}
    }}

    (liquid in gas)
    {{
        type            SchillerNaumann;
        residualRe      1e-3;
        swarmCorrection
        {{
            type        none;
        }}
    }}
);

virtualMass
(
    (gas in liquid)
    {{
        type            constantCoefficient;
        Cvm             0.5;
    }}

    (liquid in gas)
    {{
        type            constantCoefficient;
        Cvm             0.5;
    }}
);

// 界面组分相变（逐相声明）：气相侧取饱和，液相侧无（纯水）
interfaceComposition.gas
(
    (gas and liquid)
    {{
        type saturated;
        species ( H2O );
        Le 1.0;
        saturationPressure
        {{
            type ArdenBuck;
        }}
    }}
);

interfaceComposition.liquid
();

heatTransfer.gas
(
    (gas in liquid)
    {{
        type spherical;
        residualAlpha 1e-4;
    }}

    (liquid in gas)
    {{
        type RanzMarshall;
        residualAlpha 1e-4;
    }}
);

heatTransfer.liquid
(
    (gas in liquid)
    {{
        type RanzMarshall;
        residualAlpha 1e-4;
    }}

    (liquid in gas)
    {{
        type spherical;
        residualAlpha 1e-4;
    }}
);

diffusiveMassTransfer.gas
(
    (gas in liquid)
    {{
        type spherical;
        Le 1.0;
    }}

    (liquid in gas)
    {{
        type Frossling;
        Le 1.0;
    }}
);

diffusiveMassTransfer.liquid
(
    (gas in liquid)
    {{
        type Frossling;
        Le 1.0;
    }}

    (liquid in gas)
    {{
        type spherical;
        Le 1.0;
    }}
);

phaseTransfer
();

lift
();

wallLubrication
();

turbulentDispersion
();

interfaceCompression
();

// Minimum allowable pressure
pMin            10000;
"""


def physicalProperties_gas() -> str:
    # air + H2O（janaf 多项式系数取自 OpenFOAM 官方 bubbleColumnEvaporating 教程）
    return header("dictionary", "thermophysicalProperties.gas", "constant") + """
thermoType
{
    type            heRhoThermo;
    mixture         multiComponentMixture;
    transport       sutherland;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}

// multiComponentMixture 用 inertSpecie（不是 defaultSpecie！）
// defaultSpecie 属于纯组分/反应体系，会让框架改用 Ydefault 代表 air，
// 与显式给出的 air.gas 冲突 → correctSpecies() 里 ΣY 失衡、崩溃。
inertSpecie air;

species
(
    air
    H2O
);

H2O
{
    specie
    {
        molWeight       18.0153;
    }
    thermodynamics
    {
        Tlow            200;
        Thigh           3500;
        Tcommon         1000;
        highCpCoeffs    ( 3.03399 0.00217692 -1.64073e-07 -9.7042e-11 1.68201e-14 -30004.3 4.96677 );
        lowCpCoeffs     ( 4.19864 -0.00203643 6.5204e-06 -5.48797e-09 1.77198e-12 -30293.7 -0.849032 );
    }
    transport
    {
        As              1.67212e-06;
        Ts              170.672;
    }
}

air
{
    specie
    {
        molWeight       28.9596;
    }
    thermodynamics
    {
        Tlow            200;
        Thigh           3500;
        Tcommon         1000;
        highCpCoeffs    ( 3.57304 -7.24383e-04 1.67022e-06 -1.26501e-10 -4.20580e-13 -1047.41 3.12431 );
        lowCpCoeffs     ( 3.09589 1.22835e-03 -4.14267e-07 6.56910e-11 -3.87021e-15 -983.191 5.34161 );
    }
    transport
    {
        As              1.67212e-06;
        Ts              170.672;
    }
}
"""


def physicalProperties_liquid() -> str:
    # 纯水（rhoConst 近似；密度/Cv/粘度取常温值）
    return header("dictionary", "thermophysicalProperties.liquid", "constant") + """
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          eConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleInternalEnergy;
}

defaultSpecie H2O;

mixture
{
    specie
    {
        molWeight       18.0153;
    }
    equationOfState
    {
        rho             998.2;      // kg/m³ 液态水
    }
    thermodynamics
    {
        Cv              4182;       // J/(kg·K)
        Hf              -1.5879e+07;  // J/kg 标准生成焓（与水蒸气之差即潜热 h_fg）
    }
    transport
    {
        mu              1.003e-3;   // Pa·s
        Pr              7.0;        // k = mu*cp/Pr ≈ 0.6 W/(m·K)
    }
}
"""


# ---- 0/ 场文件 ----
def field_header(cls: str, obj: str, dim: str, internal: str) -> str:
    return (header(cls, obj, "0") +
            f"dimensions      {dim};\n"
            f"internalField   {internal};\n"
            "boundaryField\n{\n")


def bnd_inlet_value(v: str) -> str:
    return f'        {P_INLET}\n        {{\n            type            fixedValue;\n            value           {v};\n        }}\n'


def bnd_outlet_inletOutlet(v: str) -> str:
    return (f'        {P_OUTLET}\n        {{\n            type            inletOutlet;\n'
            f'            inletValue     {v};\n            value           {v};\n        }}\n')


def bnd_zeroGradient(name: str) -> str:
    return f'        {name}\n        {{\n            type            zeroGradient;\n        }}\n'


def bnd_cold_T() -> str:
    return f'        {P_COLD}\n        {{\n            type            fixedValue;\n            value           uniform {c2k(T_C):.4f};\n        }}\n'


def bnd_wall_T() -> str:
    return f'        {P_WALL}\n        {{\n            type            zeroGradient;\n        }}\n'


def close() -> str:
    return "}\n"


def field_alpha(phase: str, val_inlet: float, val_other_wall: str) -> str:
    # 关键：分散相（液相）不能初始为严格的 0。
    # alpha.liquid ≡ 0 会让液相的能量/组分方程退化为奇异矩阵，
    # 求解出 NaN 后经界面源项（0*NaN）污染气相温度场，
    # 表现为第一步 e.gas 残差上千、T 冲到 -10^4 K。
    # 因此给液相一个极小的种子体积分数（相当于空气中已有的微细雾滴），
    # 同时保证 condense 模型有可用的界面面积。10^-6 量级对结果影响可忽略。
    if phase == LIQ:
        internal = f"uniform {ALPHA_LIQ_SEED:g}"
    else:
        internal = f"uniform {1.0 - ALPHA_LIQ_SEED:g}"
    return (field_header("volScalarField", f"alpha.{phase}",
                         "[0 0 0 0 0 0 0]", internal) +
            bnd_inlet_value(f"uniform {val_inlet}") +
            bnd_outlet_inletOutlet(f"uniform {val_inlet}") +
            bnd_zeroGradient(P_COLD) +
            bnd_zeroGradient(P_WALL) +
            close())


def field_T(phase: str) -> str:
    return (field_header("volScalarField", f"T.{phase}",
                         "[0 0 0 1 0 0 0]", f"uniform {c2k(T_IN):.4f}") +
            bnd_inlet_value(f"uniform {c2k(T_IN):.4f}") +
            bnd_outlet_inletOutlet(f"uniform {c2k(T_IN):.4f}") +
            bnd_cold_T() +
            bnd_wall_T() +
            close())


def field_U(phase: str, inlet_vec: str) -> str:
    return (field_header("volVectorField", f"U.{phase}",
                         "[0 1 -1 0 0 0 0]", "uniform (0 0 0)") +
            bnd_inlet_value(f"uniform {inlet_vec}") +
            bnd_outlet_inletOutlet("uniform (0 0 0)") +
            f'        {P_COLD}\n        {{\n            type            noSlip;\n        }}\n'
            f'        {P_WALL}\n        {{\n            type            noSlip;\n        }}\n' +
            close())


def field_species(phase: str, sp: str, val_inlet: float) -> str:
    # 物种出口用 inletOutlet 时必须显式指定该相的 face flux（phi.<phase>），
    # 否则多相下无法判断流入/流出，回流出错——官方教程同写法。
    outlet = (f'        {P_OUTLET}\n        {{\n'
              f'            type            inletOutlet;\n'
              f'            phi             phi.{phase};\n'
              f'            inletValue      uniform {val_inlet};\n'
              f'            value           uniform {val_inlet};\n'
              f'        }}\n')
    return (field_header("volScalarField", f"{sp}.{phase}",
                         "[0 0 0 0 0 0 0]", f"uniform {val_inlet}") +
            bnd_inlet_value(f"uniform {val_inlet}") +
            outlet +
            bnd_zeroGradient(P_COLD) +
            bnd_zeroGradient(P_WALL) +
            close())


def bnd_fixedFluxPressure(name: str) -> str:
    """p_rgh 在壁面必须用 fixedFluxPressure：它把重力/体积力放进压力梯度，
    用 zeroGradient 会在有重力时产生壁面伪速度（凝液排水算例尤其明显）。"""
    return f'        {name}\n        {{\n            type            fixedFluxPressure;\n        }}\n'


def field_p_rgh() -> str:
    """p_rgh 在 multiphaseEulerFoam 里是【绝压】，不是表压。

    依据 OF8 源码 pU/pEqn.H:391：每个压力修正步无条件执行
        p = max(p_rgh + rho*gh, pMin);   p_rgh = p - rho*gh;
    气相为 perfectGas，密度由 p 算出。若这里按表压写 uniform 0，第一步 p 就被
    pMin=10000 Pa 夹住 → 密度只剩常压的 1/10 → 出口(0)与内部(1atm)差一个大气压
    → 速度一步冲到 3000 m/s → FOAM_SIGFPE。
    （2026-10-20 实测：改成本函数当前的绝压口径后，同一算例 Co 由 143 降到 0.008，
      连续跑到 2 ms 零 FATAL。）

    口径与官方教程 multiphaseEulerFoam/laminar/bubbleColumnEvaporating/0/p_rgh 一致。
    """
    return (field_header("volScalarField", "p_rgh",
                         "[1 -1 -2 0 0 0 0]", f"uniform {P_OP:g}") +
            bnd_fixedFluxPressure(P_INLET) +
            f'        {P_OUTLET}\n        {{\n'
            f'            type            prghPressure;\n'
            f'            p               uniform {P_OP:g};\n'
            f'            value           uniform {P_OP:g};\n'
            f'        }}\n' +
            bnd_fixedFluxPressure(P_COLD) +
            bnd_fixedFluxPressure(P_WALL) +
            close())


def field_Ydefault(val_inlet: float) -> str:
    """Ydefault 是【单字段、不带相后缀】的默认组分质量分数场。
    文件名 0/Ydefault 必须对应 object Ydefault（写成 Ydefault.gas 会 FATAL 报对象名不匹配）。

    注意：本算例的默认组分（air）已经由 `air.gas` 显式给出，因此
    Ydefault 必须为 0；若也填 0.9711，会与 air.gas 重复计数导致
    ΣY=1.94，correctSpecies() 归一化时除零/出现负值 → FOAM_SIGFPE 崩溃。
    OpenFOAM 8 官方教程 bubbleColumnEvaporating 的 Ydefault 同样是 uniform 0。
    """
    return (field_header("volScalarField", "Ydefault",
                         "[0 0 0 0 0 0 0]", "uniform 0") +
            bnd_inlet_value("uniform 0") +
            bnd_outlet_inletOutlet("uniform 0") +
            bnd_zeroGradient(P_COLD) +
            bnd_zeroGradient(P_WALL) +
            close())


def field_p() -> str:
    return (field_header("volScalarField", "p",
                         "[1 -1 -2 0 0 0 0]", f"uniform {P_OP}") +
            bnd_zeroGradient(P_INLET) +
            bnd_zeroGradient(P_OUTLET) +
            bnd_zeroGradient(P_COLD) +
            bnd_zeroGradient(P_WALL) +
            close())


# ============================================================================
def main(out_dir: str) -> None:
    global case_dir
    case_dir = os.path.abspath(out_dir)
    print(f"生成算例到：{case_dir}")
    print(f"入口流速 U_in = {U_IN:.4f} m/s（{MDOT:.2e} kg/s ÷ {RHO_IN} ÷ {A_IN:.3e} m²）")
    print(f"T_c = {c2k(T_C):.2f} K   T_in = {c2k(T_IN):.2f} K   Y_H2O = {Y_H2O}")

    Uin = f"(0 {-U_IN:.5f} 0)"   # 入口外向法向 +Y -> 流入 -Y

    files = {
        "system/controlDict": controlDict(),
        "system/fvSchemes": fvSchemes(),
        "system/fvSolution": fvSolution(),
        "constant/g": gravity(),
        "constant/momentumTransport.gas": momentumTransport("gas"),
        "constant/momentumTransport.liquid": momentumTransport("liquid"),
        "constant/phaseProperties": phaseProperties(),
        "constant/thermophysicalProperties.gas": physicalProperties_gas(),
        "constant/thermophysicalProperties.liquid": physicalProperties_liquid(),
        "0/alpha.gas": field_alpha("gas", 1.0, "uniform 1"),
        "0/alpha.liquid": field_alpha("liquid", 0.0, "uniform 0"),
        "0/T.gas": field_T("gas"),
        "0/T.liquid": field_T("liquid"),
        "0/U.gas": field_U("gas", Uin),
        "0/U.liquid": field_U("liquid", "(0 0 0)"),
        "0/air.gas": field_species("gas", "air", 1.0 - Y_H2O),
        "0/H2O.gas": field_species("gas", "H2O", Y_H2O),
        "0/air.liquid": field_species("liquid", "air", 0.0),
        "0/H2O.liquid": field_species("liquid", "H2O", 1.0),
        "0/Ydefault": field_Ydefault(1.0 - Y_H2O),
        "0/p_rgh": field_p_rgh(),
        "0/p": field_p(),
    }
    for rel, content in files.items():
        write(os.path.join(case_dir, rel), content)

    print("\n完成。接下来：")
    print("  1) 生成网格： fluent_cas_to_foam.py <case.cas.h5> " + case_dir)
    print("  2) 在 OpenFOAM 里： checkMesh && foamDictionary 校核 && 运行")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python tools/generate_case.py <输出算例目录>")
        print("（默认 ../case 已停用：避免把字典写进一个空壳目录）")
        raise SystemExit(2)
    main(sys.argv[1])
