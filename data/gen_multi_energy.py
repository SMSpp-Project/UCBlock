#!/usr/bin/env python3
"""
Writes the multi-energy UCBlock instances, in which a plant produces
electricity and heat: an electric node 0 with demand 60, 120, 150, 60, a
heat node 1 with demand 100, 50, 150, 100, a peaker (100 per MWh, 0 to 300)
on node 0, a boiler (60 per MWh of heat, 0 to 300) on node 1, a lossless
heat store of 100 MWh and +-50 MW on node 1 and a line 1 -> 0 of efficiency
0 that dumps the heat. The plant (100 per instant on, 500 per start-up, on
before the horizon) is written in two ways: as a ConversionUnitBlock with
its electricity on node 0 and its heat on node 1 (me_<case>_unit.nc4), and
as a ThermalUnitBlock on a private node 2 whose output the hyperarcs of the
network split into electricity and heat (me_<case>_hyperarc.nc4). The cases:

- bp:    back-pressure, heat = 1.25 electricity, electricity in [ 30 , 100 ];
- bpres: bp with a primary reserve of 5 on the electric node, PrimaryRho 0.1
         for the plant and 0.5 for the peaker, the heat following the
         electricity in the deployment [ReserveDirection];
- bppol: bp with a budget of 177 t of CO2 over all the nodes, 0.5 per MWh of
         electricity of the plant, 0.2 for the peaker and 0.25 per MWh of
         heat of the boiler;
- ec:    extraction-condensing, Q = e + 0.15 h in [ 30 , 100 ], e >= 0.8 h,
         h <= 80, at 40 per unit of Q (the hyperarc form has a condensing
         line and a back-pressure hyperarc out of node 2).

Then me_mesh_unit.nc4: a triangle of electric nodes 0, 1, 2 with lines of
susceptance 5, the line 0 - 2 limited to 40, a heat node 3 with no line, the
back-pressure plant with its electricity on node 0 and its heat on node 3, a
peaker on node 2 and a boiler on node 3, demands 100 on nodes 2 and 3 over
4 instants: Kirchhoff's voltage law caps the plant at 60.

usage: gen_multi_energy.py <output directory>

The instances are written by ncgen from their CDL, and the reference
objectives (by hand and by an independent MILP) are printed in the form the
batteries read them.
"""

import subprocess
import sys
from pathlib import Path

T = 4
DE = [60, 120, 150, 60]     # electric demand
DH = [100, 50, 150, 100]    # heat demand
CM, CV, HMX = 0.8, 0.15, 80.0
INF = "Infinity"

REF = {"bp": 21700, "bpres": 21700, "bppol": 22350, "ec": 28900,
       "mesh": 31600}


def fmt(v):
    return ", ".join(x if isinstance(x, str) else repr(float(x))
                     if not isinstance(x, int) else str(x) for x in v)


def thermal(i, pmin, pmax, lin, cst, suc, prho=0.0):
    o = ["group: UnitBlock_%d {" % i, "variables:", " double MinPower ;",
         " double MaxPower ;", " double LinearTerm ;", " double ConstTerm ;",
         " double StartUpCost ;", " double InitialPower ;",
         " int InitUpDownTime ;", " uint MinUpTime ;", " uint MinDownTime ;"]
    if prho > 0:
        o += [" double PrimaryRho ;"]
    o += [':type = "ThermalUnitBlock" ;', "data:",
          " MinPower = %r ;" % float(pmin), " MaxPower = %r ;" % float(pmax),
          " LinearTerm = %r ;" % float(lin), " ConstTerm = %r ;" % float(cst),
          " StartUpCost = %r ;" % float(suc),
          " InitialPower = %r ;" % float(pmin), " InitUpDownTime = 5 ;",
          " MinUpTime = 1 ;", " MinDownTime = 1 ;"]
    if prho > 0:
        o += [" PrimaryRho = %r ;" % prho]
    return o + ["}"]


def battery(i):
    return ["group: UnitBlock_%d {" % i, "variables:", " double MinPower ;",
            " double MaxPower ;", " double MaxStorage ;",
            " double MinStorage ;", " double InitialStorage ;",
            ':type = "BatteryUnitBlock" ;', "data:", " MinPower = -50.0 ;",
            " MaxPower = 50.0 ;", " MaxStorage = 100.0 ;",
            " MinStorage = 0.0 ;", " InitialStorage = 0.0 ;", "}"]


def conversion(i, ec, res, limits=True):
    """The plant as a ConversionUnitBlock, electricity first, heat second;
    with limits, the back-pressure plant starts up and shuts down at its
    minimum, as the ThermalUnitBlock of the hyperarc form does by default
    (the extraction-condensing one would need the limit on Q = e + 0.15 h,
    a combination, and has none)."""
    limits = limits and not ec
    rep = lambda v: v * T
    if ec:
        mn, mx, M = [0, 0], [100, HMX], 2
        A, lhs, rhs = [1, CV, 1, -CM], [30, 0], [100, INF]
        lin, init = [40, 40 * CV], [30, 0]
    else:
        mn, mx, M = [30, 37.5], [100, 125], 1
        A, lhs, rhs = [-1.25, 1], [0], [0]
        lin, init = [40, 0], [30, 37.5]
    o = ["group: UnitBlock_%d {" % i, "dimensions:", " NumberGenerators = 2 ;",
         " NumberOperatingRows = %d ;" % M, "variables:",
         " double MinPower(TimeHorizon, NumberGenerators) ;",
         " double MaxPower(TimeHorizon, NumberGenerators) ;",
         " double OperatingMatrix(NumberOperatingRows, NumberGenerators) ;",
         " double OperatingLHS(TimeHorizon, NumberOperatingRows) ;",
         " double OperatingRHS(TimeHorizon, NumberOperatingRows) ;",
         " double LinearTerm(TimeHorizon, NumberGenerators) ;",
         " double ConstTerm ;", " double StartUpCost ;",
         " double InitialPower(NumberGenerators) ;", " int InitUpDownTime ;",
         " uint MinUpTime ;", " uint MinDownTime ;"]
    if res:
        o += [" double PrimaryRho(TimeHorizon, NumberGenerators) ;",
              " double ReserveDirection(NumberGenerators, NumberGenerators) ;"]
    if limits:
        o += [" double StartUpLimit(TimeHorizon, NumberGenerators) ;",
              " double ShutDownLimit(TimeHorizon, NumberGenerators) ;"]
    o += [':type = "ConversionUnitBlock" ;', "data:",
          " MinPower = %s ;" % fmt(rep(mn)), " MaxPower = %s ;" % fmt(rep(mx)),
          " OperatingMatrix = %s ;" % fmt(A),
          " OperatingLHS = %s ;" % fmt(rep(lhs)),
          " OperatingRHS = %s ;" % fmt(rep(rhs)),
          " LinearTerm = %s ;" % fmt(rep(lin)), " ConstTerm = 100.0 ;",
          " StartUpCost = 500.0 ;", " InitialPower = %s ;" % fmt(init),
          " InitUpDownTime = 5 ;", " MinUpTime = 1 ;", " MinDownTime = 1 ;"]
    if res:
        o += [" PrimaryRho = %s ;" % fmt(rep([0.1, 0.0])),
              " ReserveDirection = 1.0, 0.0, 1.25, 1.0 ;"]
    if limits:
        o += [" StartUpLimit = %s ;" % fmt(rep(mn)),
              " ShutDownLimit = %s ;" % fmt(rep(mn))]
    return o + ["}"]


def cdl(name, unit):
    ec, res, pol = name == "ec", name == "bpres", name == "bppol"
    nn = 2 if unit else 3
    ng = 5 if unit else 4
    # the lines: ( start , [ ( end , eta ) ] , max )
    L = []
    if not unit:
        if ec:
            L.append((2, [(0, 1.0)], 100))
            L.append((2, [(0, CM / (CM + CV)), (1, 1 / (CM + CV))],
                      (CM + CV) * HMX))
        else:
            L.append((2, [(0, 1.0), (1, 1.25)], 100))
    L.append((1, [(0, 0.0)], 1000))
    nb = sum(len(l[1]) for l in L)
    o = ["netcdf me_%s {" % name, ":SMS++_file_type = 1 ;",
         "group: Block_0 {", "dimensions:", " TimeHorizon = %d ;" % T,
         " NumberUnits = 4 ;", " NumberElectricalGenerators = %d ;" % ng,
         " NumberNodes = %d ;" % nn, " NumberLines = %d ;" % len(L)]
    if not unit:
        o += [" NumberBranches = %d ;" % nb]
    if res:
        o += [" NumberPrimaryZones = 1 ;"]
    if pol:
        o += [" NumberPollutants = 1 ;", " PolOne = 1 ;"]
    br = "NumberLines" if unit else "NumberBranches"
    o += ["variables:",
          " double ActivePowerDemand(NumberNodes, TimeHorizon) ;",
          " uint GeneratorNode(NumberElectricalGenerators) ;",
          " uint StartLine(%s) ;" % br, " uint EndLine(%s) ;" % br,
          " double Efficiency(%s) ;" % br,
          " double MaxPowerFlow(NumberLines) ;",
          " double MinPowerFlow(NumberLines) ;"]
    if not unit:
        o += [" uint HyperArcID(NumberBranches) ;"]
    if res:
        o += [" uint PrimaryZones(NumberNodes) ;",
              " double PrimaryDemand(NumberPrimaryZones, TimeHorizon) ;"]
    if pol:
        o += [" uint NumberPollutantZones(NumberPollutants) ;",
              " uint PollutantZones(NumberPollutants, NumberNodes) ;",
              " double PollutantBudget(PolOne) ;",
              " double PollutantRho(PolOne, NumberPollutants, "
              "NumberElectricalGenerators) ;"]
    o += [':type = "UCBlock" ;', "data:",
          " ActivePowerDemand = %s ;" % fmt(DE + DH + ([] if unit
                                                        else [0] * T)),
          " GeneratorNode = %s ;" % fmt([0, 1, 0, 1, 1] if unit
                                        else [2, 0, 1, 1])]
    st, en, eta, hid = [], [], [], []
    for k, (s, b, _) in enumerate(L):
        for (e, h) in b:
            st.append(s); en.append(e); eta.append(h); hid.append(k)
    o += [" StartLine = %s ;" % fmt(st), " EndLine = %s ;" % fmt(en),
          " Efficiency = %s ;" % fmt(eta),
          " MaxPowerFlow = %s ;" % fmt([l[2] for l in L]),
          " MinPowerFlow = %s ;" % fmt([0] * len(L))]
    if not unit:
        o += [" HyperArcID = %s ;" % fmt(hid)]
    if res:
        o += [" PrimaryZones = %s ;" % fmt([0, 1] if unit else [0, 1, 0]),
              " PrimaryDemand = 5, 5, 5, 5 ;"]
    if pol:
        o += [" NumberPollutantZones = 1 ;",
              " PollutantZones = %s ;" % fmt([0] * nn),
              " PollutantBudget = 177.0 ;",
              " PollutantRho = %s ;" % fmt([0.5, 0.0, 0.2, 0.25, 0.0] if unit
                                           else [0.5, 0.2, 0.25, 0.0])]
    if unit:
        o += conversion(0, ec, res)
    else:
        o += thermal(0, 30, 100, 40, 100, 500, 0.1 if res else 0.0)
    o += thermal(1, 0, 300, 100, 0, 0, 0.5 if res else 0.0)
    o += thermal(2, 0, 300, 60, 0, 0)
    o += battery(3)
    return "\n".join(o + ["}", "}"]) + "\n"


def cdl_mesh():
    o = ["netcdf me_mesh {", ":SMS++_file_type = 1 ;", "group: Block_0 {",
         "dimensions:", " TimeHorizon = %d ;" % T, " NumberUnits = 3 ;",
         " NumberElectricalGenerators = 4 ;", " NumberNodes = 4 ;",
         " NumberLines = 3 ;", "variables:",
         " double ActivePowerDemand(NumberNodes, TimeHorizon) ;",
         " uint GeneratorNode(NumberElectricalGenerators) ;",
         " uint StartLine(NumberLines) ;", " uint EndLine(NumberLines) ;",
         " double LineSusceptance(NumberLines) ;",
         " double MaxPowerFlow(NumberLines) ;",
         " double MinPowerFlow(NumberLines) ;", ':type = "UCBlock" ;',
         "data:", " ActivePowerDemand = %s ;" % fmt([0] * 2 * T +
                                                    [100] * 2 * T),
         " GeneratorNode = 0, 3, 2, 3 ;", " StartLine = 0, 1, 0 ;",
         " EndLine = 1, 2, 2 ;", " LineSusceptance = 5.0, 5.0, 5.0 ;",
         " MaxPowerFlow = 1000.0, 1000.0, 40.0 ;",
         " MinPowerFlow = -1000.0, -1000.0, -40.0 ;"]
    c = conversion(0, False, False, False)
    # no start-up and no fixed cost: the plant runs at 60 at every instant
    c = [l.replace("ConstTerm = 100.0", "ConstTerm = 0.0")
          .replace("StartUpCost = 500.0", "StartUpCost = 0.0") for l in c]
    o += c
    o += thermal(1, 0, 300, 100, 0, 0)
    o += thermal(2, 0, 300, 60, 0, 0)
    return "\n".join(o + ["}", "}"]) + "\n"


def write(out, name, text):
    c = out / (name + ".cdl")
    c.write_text(text)
    subprocess.run(["ncgen", "-k", "nc4", "-o", str(out / (name + ".nc4")),
                    str(c)], check=True)
    c.unlink()


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: gen_multi_energy.py <output directory>")
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for name in ("bp", "bpres", "bppol", "ec"):
        for unit in (True, False):
            nm = "me_%s_%s" % (name, "unit" if unit else "hyperarc")
            write(out, nm, cdl(name, unit))
            print("REF_OBJ[%s.nc4]=%s" % (nm, REF[name]))
    write(out, "me_mesh_unit", cdl_mesh())
    print("REF_OBJ[me_mesh_unit.nc4]=%s" % REF["mesh"])


if __name__ == "__main__":
    main()
