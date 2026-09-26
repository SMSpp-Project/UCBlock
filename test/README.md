# test

A tester for the `UCBlock` module that needs nothing but the core SMS++
library, all of its instances being built in memory.

The dynamic programming solvers of the module, i.e., `ThermalUnitDPSolver` and
`ThermalUnitExtDPSolver` on a `ThermalUnitBlock` and `NuclearUnitExtDPSolver`
on a `NuclearUnitBlock`, are compared with each other and with a brute force
that enumerates the commitment and runs a dynamic programming over the integer
power levels, which is exact when the costs are linear and the data integer
(the economic dispatch of a fixed commitment has integer vertices then). The
instances, some of them with an optimum known in closed form and the others
drawn with a fixed seed, cover a horizon of one instant, the minimum up and
down times that the initial state or the end of the horizon cut, the ramps
with the limits of the start-up and of the shut-down, the commitment fixed ON
or OFF in some instants (an infeasible fixing included), the minimum times set
before the Variable are generated, the scale factor and the modulations of a
nuclear unit. The schedule each solver returns has to be feasible for the unit
and, when the abstract representation is there, to satisfy its Constraint and
to cost what the solver says.

Every `UnitBlock` and `NetworkBlock` of the module, and a `UCBlock` with no
unit or with one unit of each kind on two nodes, is then written in a netCDF
group, read back, written and read again: the two groups written by the Block
have to coincide, and the data read have to be those written, a single
instant and the optional data left out included.

Finally the setters of the costs and of the scale factor of a
`ThermalUnitBlock`, and the one of the demand of a `UCBlock`, are called with
a `FakeSolver` attached: it has to receive the physical Modification, the
Objective (or the node balance) has to follow the change, with the scale
factor, and the DP solver attached to the unit has to find the new optimum.

The exit code is 0 when every check passes, printing `All tests passed!!`, and
1 otherwise, each failed check printing a line. The `makefile` builds the
executable including the `UCBlock` module and the core SMS++ library.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
