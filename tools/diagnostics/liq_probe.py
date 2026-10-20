#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
对照探针：把液相相模型按官方教程 bubbleColumnEvaporatingDissolving 对齐。
  liq_mc          : 液相 multiComponentPhaseModel + Henry 空气溶解 + 多组分热物性
  liq_mc_gasconst : 在 liq_mc 基础上再把气相换成教程的 eConst 热物性（完全教程化）
用法:
  python _liq_probe.py liq_mc liq_mc_gasconst --endtime 0.05
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tools'))
from ablate import fresh, solve, summarize, print_row  # noqa: E402

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
|  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|   \\\\    /   O peration    | Version:  8                                     |
|    \\\\  /    A nd           | Website:  https://openfoam.org                  |
|     \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

LIQ_THERMO = HEADER.format(obj='thermophysicalProperties.liquid') + """
thermoType
{
    type            heRhoThermo;
    mixture         multiComponentMixture;
    transport       const;
    thermo          eConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleInternalEnergy;
}

species
(
    air
    H2O
);

inertSpecie H2O;

"(mixture|H2O)"
{
    specie
    {
        molWeight   18.0153;
    }
    equationOfState
    {
        rho         998.2;
    }
    thermodynamics
    {
        Cv          4182;
        Hf          -1.5879e+07;
    }
    transport
    {
        mu          1.003e-3;
        Pr          7.0;
    }
}

air
{
    specie
    {
        molWeight   28.9596;
    }
    equationOfState
    {
        rho         998.2;
    }
    thermodynamics
    {
        Cv          724.8;
        Hf          0;
    }
    transport
    {
        mu          1.003e-3;
        Pr          7.0;
    }
}

// ************************************************************************* //
"""

GAS_ECONST = HEADER.format(obj='thermophysicalProperties.gas') + """
thermoType
{
    type            heRhoThermo;
    mixture         multiComponentMixture;
    transport       const;
    thermo          eConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}

species
(
    air
    H2O
);

inertSpecie air;

"(mixture|air)"
{
    specie
    {
        molWeight   28.9596;
    }
    thermodynamics
    {
        Hf          0;
        Cv          724.8;
    }
    transport
    {
        mu          1.84e-05;
        Pr          0.7;
    }
}

H2O
{
    specie
    {
        molWeight   18.0153;
    }
    thermodynamics
    {
        Hf          -1.3435e+07;
        Cv          1396.3;
    }
    transport
    {
        mu          1.84e-05;
        Pr          0.7;
    }
}

// ************************************************************************* //
"""


def make_liq_multicomponent(d):
    pp = os.path.join(d, 'constant', 'phaseProperties')
    t = open(pp, encoding='utf-8').read()

    old_type = 'type            purePhaseModel;'
    assert old_type in t, 'purePhaseModel not found'
    t = t.replace(old_type, 'type            multiComponentPhaseModel;')

    old_ic = 'interfaceComposition.liquid\n();'
    assert old_ic in t, 'interfaceComposition.liquid empty block not found'
    new_ic = ('interfaceComposition.liquid\n'
              '(\n'
              '    (liquid and gas)\n'
              '    {\n'
              '        type Henry;\n'
              '        species ( air );\n'
              '        k ( 1.492e-2 );\n'
              '        Le 1.0;\n'
              '    }\n'
              ');')
    t = t.replace(old_ic, new_ic)
    open(pp, 'w', encoding='utf-8').write(t)

    open(os.path.join(d, 'constant', 'thermophysicalProperties.liquid'),
         'w', encoding='utf-8').write(LIQ_THERMO)


def probe_liq_mc(d):
    make_liq_multicomponent(d)


def probe_liq_mc_gasconst(d):
    make_liq_multicomponent(d)
    open(os.path.join(d, 'constant', 'thermophysicalProperties.gas'),
         'w', encoding='utf-8').write(GAS_ECONST)


PROBES = {
    'liq_mc': probe_liq_mc,
    'liq_mc_gasconst': probe_liq_mc_gasconst,
}


def main():
    argv = sys.argv[1:]
    names = [a for a in argv if not a.startswith('--')]
    endtime = '0.05'
    for i, a in enumerate(argv):
        if a.startswith('--endtime='):
            endtime = a.split('=', 1)[1]
        elif a == '--endtime' and i + 1 < len(argv):
            endtime = argv[i + 1]
            names = [n for n in names if n != argv[i + 1]]
    if not names:
        print('用法: python _liq_probe.py liq_mc liq_mc_gasconst [--endtime=0.05]')
        return 1
    for name in names:
        fn = PROBES[name]
        d = fresh(name, endtime=endtime)
        fn(d)
        rc, out, log = solve(d, timeout=900)
        s = summarize(rc, out)
        print_row(name, s)
        print(f'{"":22s}    log: {log}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
