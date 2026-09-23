# Writes a unit commitment instance with K times the units of a given one:
# each unit is copied K times, the costs of every copy after the first
# multiplied by a factor drawn uniformly in [ 1 - p , 1 + p ] (so that the
# copies are not identical), and the quantities that count the units, the
# generators and the storages multiplied by K with them, so that the instance
# keeps the same shape with K times the components.
# Usage: python replicate_units.py in.nc4 out.nc4 K [ p [ seed ] ]
import sys
import numpy as np
from netCDF4 import Dataset

src, dst, K = sys.argv[1], sys.argv[2], int(sys.argv[3])
p = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 1)

COSTS = ("QuadTerm", "LinearTerm", "ConstTerm", "StartUpCost")

# the dimensions that count the replicated components: a variable indexed by
# one of them is repeated K times along that axis
TILED = ("NumberUnits", "NumberElectricalGenerators", "NumberStorages")

# what the K times the units ask for, i.e. K times the demand and K times the
# budget: leaving them as they are would make the instance a different problem
# and not the same one with K times the components
SCALED = ("ActivePowerDemand", "ReactivePowerDemand", "PrimaryDemand",
          "SecondaryDemand", "InertiaDemand", "PollutantBudget",
          "PollutantMinBudget")


def tiled_axis(var, dims):
    """the axis of var that is replicated, if any, refusing more than one"""
    ax = [i for i, d in enumerate(var.dimensions) if d in dims]
    if len(ax) > 1:
        sys.exit(f"{var.name} is indexed by {len(ax)} replicated dimensions, "
                 "which this tool does not know how to replicate")
    return ax[0] if ax else None


with Dataset(src) as fi, Dataset(dst, "w", format="NETCDF4") as fo:
    fo.setncatts({a: fi.getncattr(a) for a in fi.ncattrs()})
    bi = fi.groups["Block_0"]
    bo = fo.createGroup("Block_0")
    attrs = {a: bi.getncattr(a) for a in bi.ncattrs()}
    if "n_gen" in attrs:                     # informational, kept truthful
        attrs["n_gen"] = int(attrs["n_gen"]) * K
    bo.setncatts(attrs)

    n = len(bi.dimensions["NumberUnits"])
    for name, dim in bi.dimensions.items():
        bo.createDimension(name, len(dim) * K if name in TILED else len(dim))

    for name, var in bi.variables.items():
        v = bo.createVariable(name, var.datatype, var.dimensions)
        ax = tiled_axis(var, TILED)
        if ax is not None:                   # one entry per unit or generator
            if name in SCALED:
                sys.exit(f"{name} is both indexed by a replicated dimension "
                         "and a quantity to scale: which one it is has to be "
                         "decided before this instance can be replicated")
            v[:] = np.concatenate([var[:]] * K, axis=ax)
        elif name in SCALED:
            v[:] = var[:] * K
        else:
            v[:] = var[:]

    units = [g for g in bi.groups if g.startswith("UnitBlock_")]
    others = [g for g in bi.groups if not g.startswith("UnitBlock_")]

    for g in others:  # e.g., NetworkData, which the copies share
        gi, go = bi.groups[g], bo.createGroup(g)
        go.setncatts({a: gi.getncattr(a) for a in gi.ncattrs()})
        for name, dim in gi.dimensions.items():
            if name in TILED:
                sys.exit(f"group {g} is indexed by {name}, which this tool "
                         "replicates only in Block_0")
            go.createDimension(name, len(dim))
        for name, var in gi.variables.items():
            if var.dimensions and any(d in TILED for d in var.dimensions):
                sys.exit(f"{g}/{name} is indexed by a replicated dimension, "
                         "which this tool does not replicate outside Block_0")
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
