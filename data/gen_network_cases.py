#!/usr/bin/env python3
"""
Writes the UCBlock instances of the network edge cases, with the optimum of
PyPSA as their reference objective: small dispatch networks whose DC lines
have a susceptance, so that the three formulations of DCNetworkBlock (PTDF,
CYCLE, KIRCHHOFF) are each exercised, on the topologies they handle apart:

- par:   two DC lines between the same two nodes, with different reactance
         and capacity, so that Kirchhoff decides how the flow splits and the
         capacity of the weaker line binds;
- opp:   the same with the second line in the opposite direction;
- tri:   a triangle with a parallel line on one side, the loop flow of the
         cheap generator meeting the limit of a line;
- link:  two DC components joined by a lossless link, and a third node
         reached from them only by a link with losses.

usage: gen_network_cases.py <output directory>

The instances are written as smspp_net_<case>.nc, and the reference
objectives are printed in the form batch-pypsa reads them.
"""

import sys
from pathlib import Path

import numpy as np
import pypsa

from pypsa2smspp.transformation import Transformation

T = 4  # time steps
LOAD = np.array([60.0, 90.0, 120.0, 80.0])


def base(nodes):
    n = pypsa.Network()
    n.set_snapshots(range(T))
    for b in nodes:
        n.add("Bus", b, carrier="AC")
    n.add("Carrier", "AC")
    n.add("Carrier", "gas")
    n.add("Carrier", "load")
    return n


def gen(n, name, bus, p_nom, cost):
    n.add("Generator", name, bus=bus, p_nom=p_nom, marginal_cost=cost,
          carrier="gas")


def shed(n, bus):
    # load shedding at a high cost, so that every instance is feasible
    n.add("Generator", f"shed {bus}", bus=bus, p_nom=1e4, marginal_cost=1e4,
          carrier="load")


def case_par(opposite):
    n = base(["a", "b"])
    gen(n, "cheap", "a", 200, 10)
    gen(n, "dear", "b", 200, 50)
    n.add("Load", "load", bus="b", p_set=LOAD)
    n.add("Line", "l0", bus0="a", bus1="b", x=0.1, r=0.01, s_nom=40)
    if opposite:
        n.add("Line", "l1", bus0="b", bus1="a", x=0.3, r=0.01, s_nom=100)
    else:
        n.add("Line", "l1", bus0="a", bus1="b", x=0.3, r=0.01, s_nom=100)
    shed(n, "b")
    return n


def case_tri():
    n = base(["a", "b", "c"])
    gen(n, "cheap", "a", 300, 10)
    gen(n, "mid", "b", 300, 30)
    n.add("Load", "load", bus="c", p_set=LOAD)
    n.add("Line", "ab", bus0="a", bus1="b", x=0.1, r=0.01, s_nom=200)
    n.add("Line", "bc", bus0="b", bus1="c", x=0.1, r=0.01, s_nom=200)
    n.add("Line", "ac", bus0="a", bus1="c", x=0.1, r=0.01, s_nom=50)
    n.add("Line", "ac2", bus0="c", bus1="a", x=0.2, r=0.01, s_nom=20)
    shed(n, "c")
    return n


def case_link():
    n = base(["a", "b", "c", "d", "e"])
    gen(n, "cheap", "a", 300, 10)
    gen(n, "mid", "c", 300, 30)
    n.add("Load", "load b", bus="b", p_set=LOAD)
    n.add("Load", "load e", bus="e", p_set=LOAD / 3)
    n.add("Line", "ab", bus0="a", bus1="b", x=0.1, r=0.01, s_nom=70)
    n.add("Line", "ab2", bus0="b", bus1="a", x=0.2, r=0.01, s_nom=70)
    n.add("Line", "cd", bus0="c", bus1="d", x=0.1, r=0.01, s_nom=200)
    # the two DC components, and the node reached only by a lossy link
    n.add("Link", "bc", bus0="b", bus1="c", p_nom=50, p_min_pu=-1,
          efficiency=1.0)
    n.add("Link", "ae", bus0="a", bus1="e", p_nom=100, efficiency=0.5)
    n.add("Link", "de", bus0="d", bus1="e", p_nom=100, efficiency=0.8)
    for b in "bce":
        shed(n, b)
    return n


def case_comp():
    n = base(["a", "b", "c", "d"])
    gen(n, "cheap", "a", 300, 10)
    gen(n, "mid", "c", 300, 30)
    n.add("Load", "load b", bus="b", p_set=LOAD)
    n.add("Load", "load d", bus="d", p_set=LOAD / 2)
    n.add("Line", "ab", bus0="a", bus1="b", x=0.1, r=0.01, s_nom=200)
    n.add("Line", "cd", bus0="c", bus1="d", x=0.1, r=0.01, s_nom=200)
    n.add("Link", "bc", bus0="b", bus1="c", p_nom=50, p_min_pu=-1,
          efficiency=1.0)
    for b in "bd":
        shed(n, b)
    return n


def case_loss():
    n = base(["a", "b", "e"])
    gen(n, "cheap", "a", 300, 10)
    n.add("Load", "load b", bus="b", p_set=LOAD)
    n.add("Load", "load e", bus="e", p_set=LOAD / 3)
    n.add("Line", "ab", bus0="a", bus1="b", x=0.1, r=0.01, s_nom=200)
    n.add("Link", "ae", bus0="a", bus1="e", p_nom=100, efficiency=0.5)
    for b in "be":
        shed(n, b)
    return n


CASES = {
    "comp": case_comp,
    "loss": case_loss,
    "par": lambda: case_par(False),
    "opp": lambda: case_par(True),
    "tri": case_tri,
    "link": case_link,
}


def main(outdir):
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / "work"
    work.mkdir(exist_ok=True)
    for case, make in CASES.items():
        n = make()
        ref = n.copy()
        ref.optimize(solver_name="gurobi")
        obj = float(ref.objective + getattr(ref, "objective_constant", 0.0))

        name = f"net_{case}"
        tr = Transformation(
            capacity_expansion_ucblock=True,
            workdir=work,
            name=name,
            overwrite=True,
            fp_temp="smspp_{name}.nc",
            fp_log="smspp_{name}_log.txt",
            fp_solution="smspp_{name}_solution.nc",
            configfile="auto",
            pysmspp_options={},
        )
        tr.run(n, verbose=False)
        got = float(tr.result.objective_value)
        (work / f"smspp_{name}.nc").replace(outdir / f"smspp_{name}.nc")
        print(f"REF_OBJ[smspp_{name}.nc]={obj:.9e}  # SMS++ {got:.9e}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
