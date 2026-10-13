#!/usr/bin/env python3
"""
Writes the bus UCBlock instances whose own Objective has a constant term:
two thermal units over 4 instants, demand 150, 220, 260, 120, the first
unit (0 to 200, 10 per MWh, 50 per instant on, 100 per start-up) on before
the horizon at 100, the second (30 to 150, 20 per MWh, 30 per instant on,
80 per start-up) off, minimum up and down times of 1, ramps as wide as the
maximum powers. With a single node the constant of the network
("NetworkConstantTerms") is kept in the Objective of the UCBlock, +1000 in
bus_const_pos.nc4 and -1000 in bus_const_neg.nc4.

The optimum, by hand, is 9070 plus the constant: the second unit is needed
at the two peaks, it is started at instant 1 and, its shut-down limit
being its minimum power, it can only be switched off from 30, so it is kept
on until the end at 30 rather than brought down at instant 2.

usage: gen_constant.py <output directory>

The instances are written by ncgen from their CDL, and the reference
objectives are printed in the form the batteries read them.
"""

import subprocess
import sys
from pathlib import Path

DEMAND = [150, 220, 260, 120]
BASE = 9070
CONST = {"pos": 1000, "neg": -1000}


def thermal(i, pmin, pmax, lin, cst, suc, p0, ud):
    return ["group: UnitBlock_%d {" % i, "variables:", " double MinPower ;",
            " double MaxPower ;", " double DeltaRampUp ;",
            " double DeltaRampDown ;", " double StartUpCost ;",
            " double LinearTerm ;", " double ConstTerm ;",
            " double InitialPower ;", " int InitUpDownTime ;",
            " uint MinUpTime ;", " uint MinDownTime ;",
            ':type = "ThermalUnitBlock" ;', "data:",
            " MinPower = %r ;" % float(pmin),
            " MaxPower = %r ;" % float(pmax),
            " DeltaRampUp = %r ;" % float(pmax),
            " DeltaRampDown = %r ;" % float(pmax),
            " StartUpCost = %r ;" % float(suc),
            " LinearTerm = %r ;" % float(lin),
            " ConstTerm = %r ;" % float(cst),
            " InitialPower = %r ;" % float(p0),
            " InitUpDownTime = %d ;" % ud, " MinUpTime = 1 ;",
            " MinDownTime = 1 ;", "}"]


def cdl(const):
    o = ["netcdf bus_const {", ":SMS++_file_type = 1 ;", "group: Block_0 {",
         "dimensions:", " TimeHorizon = %d ;" % len(DEMAND),
         " NumberUnits = 2 ;", " NumberNetworks = 1 ;", "variables:",
         " double ActivePowerDemand(TimeHorizon) ;",
         " double NetworkConstantTerms(NumberNetworks) ;",
         ':type = "UCBlock" ;', "data:",
         " ActivePowerDemand = %s ;" % ", ".join(str(d) for d in DEMAND),
         " NetworkConstantTerms = %r ;" % float(const),
         "group: NetworkData {", "dimensions:", " NumberNodes = 1 ;", "}"]
    o += thermal(0, 50, 200, 10, 50, 100, 100, 1)
    o += thermal(1, 30, 150, 20, 30, 80, 0, -1)
    return "\n".join(o + ["}", "}"]) + "\n"


def write(out, name, text):
    c = out / (name + ".cdl")
    c.write_text(text)
    subprocess.run(["ncgen", "-k", "nc4", "-o", str(out / (name + ".nc4")),
                    str(c)], check=True)
    c.unlink()


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: gen_constant.py <output directory>")
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for name, const in CONST.items():
        nm = "bus_const_%s" % name
        write(out, nm, cdl(const))
        print("CONST_OBJ[%s.nc4]=%s" % (nm, BASE + const))


if __name__ == "__main__":
    main()
