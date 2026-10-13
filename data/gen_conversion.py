#!/usr/bin/env python3
"""
Writes the UCBlock instances with a ConversionUnitBlock that pypsa2smspp
makes of the committable links of its test networks, i.e., the xlsx cases
chp_*, chpx_* and chpy_* of test/configs/data/test of pypsa2smspp that its
sweep test_ucblock solves in the UCBlock form: a back-pressure CHP between a
gas, an electric and a heat node, alone, with no line, with ramps and an
empty start-up range, with a heat store and the CO2 of the fuel as third
output, in a meshed AC network with a hyperarc, with an extendable boiler or
a solar generator, with a heat pump whose COP changes over time (also
delayed), with a delayed heat output, with modules, with the modules built
(also with one on before the horizon) and with a continuous design (also
with ramps and start-up and shut-down limits, and on before the horizon).

Each network is built as test_ucblock builds it (the slack units included),
solved by PyPSA with HiGHS, whose optimum is the reference, and written by
pypsa2smspp as the SMS++ file of its UCBlock, <case>.nc4.

usage: gen_conversion.py <pypsa2smspp directory> <output directory>

The pypsa2smspp directory is a checkout with ConversionUnitBlock (option
enable_conversion_units) and with its test cases, in an environment with
PyPSA, HiGHS and pySMSpp; the reference objectives are printed in the form
the batteries read them.
"""

import os
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: gen_conversion.py <pypsa2smspp directory> "
                 "<output directory>")
    root = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)

    # the helpers of the tests of pypsa2smspp read their data relative to
    # the test directory
    sys.path.insert(0, str(root))
    sys.path.insert(0, str(root / "test"))
    os.chdir(root / "test")

    from conftest import create_test_config
    from network_definition import NetworkDefinition
    from pypsa2smspp.network_correction import (add_slack_unit,
                                                clean_ciclicity_storage)
    from pypsa2smspp.transformation import Transformation

    cases = sorted(p for p in (root / "test/configs/data/test").glob("*.xlsx")
                   if p.name.startswith(("chp_", "chpx_", "chpy_")))

    for xlsx in cases:
        n = clean_ciclicity_storage(NetworkDefinition(
            create_test_config(xlsx)).n)
        n = add_slack_unit(n)

        reference = n.copy()
        reference.optimize(solver_name="highs")
        objective = float(reference.objective
                          + getattr(reference, "objective_constant", 0.0))

        transformation = Transformation(
            capacity_expansion_ucblock=True,
            enable_thermal_units=bool(n.generators.committable.any()),
            enable_conversion_units=bool(n.links.committable.any()),
            workdir=out, name=xlsx.stem, overwrite=True,
        )
        transformation.create_model(n, verbose=False)
        transformation.sms_network.to_netcdf(str(out / f"{xlsx.stem}.nc4"),
                                             force=True)

        print("REF_OBJ[%s.nc4]=%.9e" % (xlsx.stem, objective))


if __name__ == "__main__":
    main()
