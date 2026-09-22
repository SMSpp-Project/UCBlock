# Writes a unit commitment instance with K times the thermal units of a given
# one: each unit is copied K times, the costs of every copy after the first
# multiplied by a factor drawn uniformly in [ 1 - p , 1 + p ] (so that the
# copies are not identical), and the demand multiplied by K, so that the
# instance keeps the same shape with K times the components.
# Usage: python replicate_units.py in.nc4 out.nc4 K [ p [ seed ] ]
import sys
import numpy as np
from netCDF4 import Dataset

src, dst, K = sys.argv[1], sys.argv[2], int(sys.argv[3])
p = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 1)
COSTS = ("QuadTerm", "LinearTerm", "ConstTerm", "StartUpCost")

with Dataset(src) as fi, Dataset(dst, "w", format="NETCDF4") as fo:
    fo.setncatts({a: fi.getncattr(a) for a in fi.ncattrs()})
    bi = fi.groups["Block_0"]
    bo = fo.createGroup("Block_0")
    bo.setncatts({a: bi.getncattr(a) for a in bi.ncattrs()})
    n = len(bi.dimensions["NumberUnits"])
    for name, dim in bi.dimensions.items():
        bo.createDimension(name, n * K if name == "NumberUnits" else len(dim))
    for name, var in bi.variables.items():
        v = bo.createVariable(name, var.datatype, var.dimensions)
        v[:] = var[:] * K if name == "ActivePowerDemand" else var[:]
    units = [g for g in bi.groups if g.startswith("UnitBlock_")]
    others = [g for g in bi.groups if not g.startswith("UnitBlock_")]
    for g in others:  # e.g., NetworkData, copied as it is
        gi, go = bi.groups[g], bo.createGroup(g)
        go.setncatts({a: gi.getncattr(a) for a in gi.ncattrs()})
        for name, dim in gi.dimensions.items():
            go.createDimension(name, len(dim))
        for name, var in gi.variables.items():
            go.createVariable(name, var.datatype, var.dimensions)[:] = var[:]
    assert len(units) == n
    for k in range(K):
        for i in range(n):
            gi = bi.groups[f"UnitBlock_{i}"]
            go = bo.createGroup(f"UnitBlock_{k * n + i}")
            go.setncatts({a: gi.getncattr(a) for a in gi.ncattrs()})
            for name, dim in gi.dimensions.items():
                go.createDimension(name, len(dim))
            f = 1.0 if k == 0 else rng.uniform(1 - p, 1 + p)
            for name, var in gi.variables.items():
                val = var[:]
                if name in COSTS:
                    val = val * f
                go.createVariable(name, var.datatype, var.dimensions)[:] = val
print(dst, n * K, "units")
