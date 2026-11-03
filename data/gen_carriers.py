#!/usr/bin/env python3
"""
Writes the UCBlock instances with more than one energy carrier from small
PyPSA networks, each solved by PyPSA with HiGHS, whose optimum is the
reference:

- h2:   an electric, a hydrogen and a gas bus: wind and an electric load, an
        electrolyzer (electricity to hydrogen, 0.7), a fuel cell (hydrogen
        to electricity, 0.5) and a store of hydrogen, a gas turbine (gas to
        electricity, 0.4) and the supply of the gas;
- heat: an electric, a heat and a gas bus: a CHP whose link has two outputs
        (gas to electricity 0.4 and to heat 0.45), a heat pump whose
        efficiency, the COP, changes with the snapshot, a boiler, the supply
        of the gas, an electric backup at 100 per MWh and the two loads (the
        heat of the CHP is bound to its electricity, and with no way to
        dump it the backup covers what the CHP cannot);
- net:  the h2 case over two regions A and B: an electric line between them,
        and a hydrogen pipeline, a link with efficiency 1 that can go both
        ways, between their two hydrogen buses; the wind and the electrolyzer
        are in A, the load and the fuel cell in B.

In the UCBlock each carrier of PyPSA (AC, H2, gas, heat) is a carrier, the
AC one being the electricity, carrier 0; each bus is a node of its carrier.
A link whose two buses have the same carrier and efficiency 1 is a line of
the network of that carrier; any other link is an IntermittentUnitBlock
whose power is the p0 of the link, in [ p_min_pu , p_max_pu ] * p_nom, with
the outputs -1 on the carrier of bus0 and efficiency, efficiency2, ... on
those of bus1, bus2, ... (see UnitBlock::get_output_efficiency()), and the
marginal cost of the link on p0, as in PyPSA. A generator is an
IntermittentUnitBlock with its output on the carrier of its bus, a store a
BatteryUnitBlock with its output on the carrier of its bus, a load the
demand of its node. The lines are written as transport lines, i.e., with
their limits and without susceptance, which is enough since the networks
here are radial.

usage: gen_carriers.py <output directory>

It runs in an environment with PyPSA, HiGHS and netCDF4; the reference
objectives are printed in the form the batteries read them.
"""

import sys
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import pypsa

T = 4
WIND = [0.9, 0.2, 0.6, 0.1]
ELOAD = [60.0, 80.0, 40.0, 90.0]
HLOAD = [50.0, 70.0, 30.0, 60.0]
COP = [3.0, 2.5, 3.5, 2.0]


def network_h2(regions):
    n = pypsa.Network()
    n.set_snapshots(range(T))
    for c in ("AC", "H2", "gas"):
        n.add("Carrier", c)
    for r in regions:
        n.add("Bus", "el " + r, carrier="AC")
        n.add("Bus", "h2 " + r, carrier="H2")
    n.add("Bus", "gas", carrier="gas")
    a, b = regions[0], regions[-1]
    n.add("Generator", "wind", bus="el " + a, p_nom=150,
          p_max_pu=pd.Series(WIND, n.snapshots), marginal_cost=0.0)
    n.add("Generator", "gas supply", bus="gas", p_nom=1000,
          marginal_cost=20.0)
    n.add("Load", "load", bus="el " + b, p_set=pd.Series(ELOAD, n.snapshots))
    n.add("Link", "electrolyzer", bus0="el " + a, bus1="h2 " + a,
          efficiency=0.7, p_nom=80, marginal_cost=1.0)
    n.add("Link", "fuel cell", bus0="h2 " + b, bus1="el " + b,
          efficiency=0.5, p_nom=50, marginal_cost=2.0)
    n.add("Link", "gas turbine", bus0="gas", bus1="el " + b,
          efficiency=0.4, p_nom=200, marginal_cost=3.0)
    n.add("Store", "h2 store", bus="h2 " + b, e_nom=200)
    if len(regions) > 1:
        n.add("Line", "AB", bus0="el " + a, bus1="el " + b, x=0.1,
              s_nom=60)
        n.add("Link", "pipeline", bus0="h2 " + a, bus1="h2 " + b,
              efficiency=1.0, p_nom=40, p_min_pu=-1)
    return n


def network_heat():
    n = pypsa.Network()
    n.set_snapshots(range(T))
    for c in ("AC", "heat", "gas"):
        n.add("Carrier", c)
    n.add("Bus", "el", carrier="AC")
    n.add("Bus", "heat", carrier="heat")
    n.add("Bus", "gas", carrier="gas")
    n.add("Generator", "wind", bus="el", p_nom=60,
          p_max_pu=pd.Series(WIND, n.snapshots), marginal_cost=0.0)
    n.add("Generator", "gas supply", bus="gas", p_nom=1000,
          marginal_cost=20.0)
    n.add("Load", "load", bus="el", p_set=pd.Series(ELOAD, n.snapshots))
    n.add("Generator", "backup", bus="el", p_nom=100, marginal_cost=100.0)
    n.add("Load", "heat load", bus="heat",
          p_set=pd.Series(HLOAD, n.snapshots))
    n.add("Link", "chp", bus0="gas", bus1="el", bus2="heat", efficiency=0.4,
          efficiency2=0.45, p_nom=300, marginal_cost=2.0)
    n.add("Link", "heat pump", bus0="el", bus1="heat",
          efficiency=pd.Series(COP, n.snapshots), p_nom=20,
          marginal_cost=0.5)
    n.add("Link", "boiler", bus0="gas", bus1="heat", efficiency=0.9,
          p_nom=100, marginal_cost=1.0)
    return n


def series(n, comp, attr, name):
    """The time series of attribute attr of component name, the static
    value if it does not change."""
    tv = getattr(n, comp + "_t")[attr]
    if name in tv.columns:
        return np.asarray(tv[name], dtype=float)
    return np.full(T, float(getattr(n, comp).at[name, attr]))


def write(n, path):
    carriers = ["AC"] + sorted(set(n.buses.carrier) - {"AC"})
    bus_node = {}
    for c in carriers:
        buses = list(n.buses.index[n.buses.carrier == c])
        for i, b in enumerate(buses):
            bus_node[b] = (carriers.index(c), i)
    nodes = [sum(1 for b in bus_node if bus_node[b][0] == k)
             for k in range(len(carriers))]
    name = {"AC": "electricity", "H2": "hydrogen"}

    # the links between two buses of one carrier with efficiency 1 are lines
    is_line = [n.buses.carrier[l.bus0] == n.buses.carrier[l.bus1]
               and np.all(series(n, "links", "efficiency", i) == 1)
               for i, l in n.links.iterrows()]
    links = [i for i, li in zip(n.links.index, is_line) if not li]
    pipes = [i for i, li in zip(n.links.index, is_line) if li]

    units = []  # ( kind , name , outputs [ ( bus , eta series ) ] )
    for g in n.generators.index:
        units.append(("gen", g, [(n.generators.at[g, "bus"], np.ones(T))]))
    for lk in links:
        outs = [(n.links.at[lk, "bus0"], -np.ones(T)),
                (n.links.at[lk, "bus1"], series(n, "links", "efficiency",
                                                lk))]
        k = 2
        while True:
            bus = n.links.at[lk, "bus%d" % k] \
                if ("bus%d" % k) in n.links.columns else ""
            if not (isinstance(bus, str) and bus):
                break
            outs.append((bus, series(n, "links", "efficiency%d" % k, lk)))
            k += 1
        units.append(("link", lk, outs))
    for s in n.stores.index:
        units.append(("store", s, [(n.stores.at[s, "bus"], np.ones(T))]))

    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    ds.setncattr("SMS++_file_type", 1)
    g = ds.createGroup("Block_0")
    g.setncattr("type", "UCBlock")
    g.createDimension("TimeHorizon", T)
    g.createDimension("NumberUnits", len(units))
    g.createDimension("NumberGenerators", len(units))
    g.createDimension("NumberCarriers", len(carriers))
    v = g.createVariable("CarrierName", str, ("NumberCarriers",))
    for k, c in enumerate(carriers):
        v[k] = name.get(c, c)

    def load(bus):
        d = np.zeros(T)
        for ld in n.loads.index[n.loads.bus == bus]:
            d += series(n, "loads", "p_set", ld)
        return d

    def network(grp, k, lines):
        """The nodes of carrier k, its lines ( bus0 , bus1 , max , min ,
        susceptance or None ) and its demand in group grp."""
        buses = [b for b in bus_node if bus_node[b][0] == k]
        grp.createDimension("NumberNodes", len(buses))
        if lines:
            grp.createDimension("NumberLines", len(lines))
            for vname, col, tp in (("StartLine", 0, "u4"),
                                   ("EndLine", 1, "u4"),
                                   ("MaxPowerFlow", 2, "f8"),
                                   ("MinPowerFlow", 3, "f8")):
                var = grp.createVariable(vname, tp, ("NumberLines",))
                var[:] = [bus_node[ln[col]][1] if col < 2 else ln[col]
                          for ln in lines]
            if all(ln[4] is not None for ln in lines):
                var = grp.createVariable("LineSusceptance", "f8",
                                         ("NumberLines",))
                var[:] = [ln[4] for ln in lines]
        dname = "ActivePowerDemand" if k == 0 else "EnergyDemand"
        var = grp.createVariable(dname, "f8", ("NumberNodes", "TimeHorizon"))
        var[:] = np.array([load(b) for b in buses])

    for k in range(len(carriers)):
        lines = []
        if k == 0:
            for ln in n.lines.index:
                lines.append((n.lines.at[ln, "bus0"], n.lines.at[ln, "bus1"],
                              n.lines.at[ln, "s_nom"],
                              -n.lines.at[ln, "s_nom"], None))
        for p in pipes:
            if bus_node[n.links.at[p, "bus0"]][0] == k:
                pn = n.links.at[p, "p_nom"]
                lines.append((n.links.at[p, "bus0"], n.links.at[p, "bus1"],
                              pn * n.links.at[p, "p_max_pu"],
                              pn * n.links.at[p, "p_min_pu"], None))
        grp = g if k == 0 else g.createGroup("Carrier_%d" % k)
        network(grp, k, lines)

    gn = g.createVariable("GeneratorNode", "i4",
                          ("NumberGenerators", "NumberCarriers"))
    table = -np.ones((len(units), len(carriers)), dtype=int)
    for i, (_, _, outs) in enumerate(units):
        for bus, _ in outs:
            k, node = bus_node[bus]
            table[i, k] = node
    gn[:] = table

    for i, (kind, nm, outs) in enumerate(units):
        u = g.createGroup("UnitBlock_%d" % i)
        u.setncattr("name", nm)
        if kind == "store":
            u.setncattr("type", "BatteryUnitBlock")
            for vn, val in (("MinPower", -1e4), ("MaxPower", 1e4),
                            ("MinStorage", 0.0),
                            ("MaxStorage", n.stores.at[nm, "e_nom"]),
                            ("InitialStorage",
                             n.stores.at[nm, "e_initial"])):
                u.createVariable(vn, "f8")[...] = val
        else:
            comp = "generators" if kind == "gen" else "links"
            pnom = getattr(n, comp).at[nm, "p_nom"]
            u.setncattr("type", "IntermittentUnitBlock")
            u.createVariable("MaxPower", "f8", ("TimeHorizon",))[:] = \
                pnom * series(n, comp, "p_max_pu", nm)
            u.createVariable("MinPower", "f8", ("TimeHorizon",))[:] = \
                pnom * series(n, comp, "p_min_pu", nm)
            u.createVariable("ActivePowerCost", "f8", ("TimeHorizon",))[:] = \
                series(n, comp, "marginal_cost", nm)
        u.createDimension("NumberOutputs", len(outs))
        oc = u.createVariable("OutputCarrier", str, ("NumberOutputs",))
        for k, (bus, _) in enumerate(outs):
            oc[k] = name.get(n.buses.at[bus, "carrier"],
                             n.buses.at[bus, "carrier"])
        oe = u.createVariable("OutputEfficiency", "f8",
                              ("TimeHorizon", "NumberOutputs"))
        oe[:] = np.array([eta for _, eta in outs]).T
    ds.close()


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: gen_carriers.py <output directory>")
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    cases = {"carriers_h2": network_h2(["A"]),
             "carriers_heat": network_heat(),
             "carriers_net": network_h2(["A", "B"])}
    for nm, n in cases.items():
        n.optimize(solver_name="highs", solver_options={"output_flag": False})
        write(n, out / (nm + ".nc4"))
        print("REF_OBJ[%s.nc4]=%.9e" % (nm, float(n.objective)))


if __name__ == "__main__":
    main()
