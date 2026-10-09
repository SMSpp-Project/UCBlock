/*--------------------------------------------------------------------------*/
/*------------------------- File UCBlockModel.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Documentation-only header: the mathematical model that a UCBlock and its
 * sub-Blocks represent, as a whole, is described in the page
 * \ref ucblock_model. The file contains no code and is included by no other
 * file; each class header of the module describes the constraints of its
 * own Block and refers to that page for the composition.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __UCBlockModel
 #define __UCBlockModel  ///< self-identification: endif at the end of the file

/*--------------------------------------------------------------------------*/
/*---------------------------- THE MODEL PAGE ------------------------------*/
/*--------------------------------------------------------------------------*/
/** \page ucblock_model The unit commitment model of UCBlock

We are concerned with the Unit Commitment (UC) problem that a UCBlock
represents together with its sub-Blocks. A set of generating units, storages
and flexibilities (each one described by a UnitBlock) has to satisfy, at
minimum total cost, the demand of active power at the nodes of a transmission
network (described by one or more NetworkBlock) at each instant of a
(typically short) time horizon. These units also have to meet requirements of
reserve and inertia on zones of nodes and limits on the emissions of
pollutants over the horizon. While each class of the module describes the
constraints of its Block, this page states the model in one place, i.e., the
notation shared by all the classes, the conventions on time, units of measure
and signs, and the composition of the sub-Blocks into one problem (cf.
UCBlock). We also describe how capacities and investments enter the problem,
its Lagrangian dual and the meaning of the dual values of the linking
constraints, and its use as the problem of one stage of a multistage model,
and we list what is not modeled. After the symbols (\ref ucbm_notation) and
the conventions (\ref ucbm_conv), we state the problem (\ref ucbm_problem) and
the contribution of each kind of sub-Block (\ref ucbm_units). Then we describe
the capacities and the investment (\ref ucbm_cap), the Lagrangian dual
(\ref ucbm_dual) and the multistage use (\ref ucbm_multi), and we close with
what is not modeled (\ref ucbm_not).

\tableofcontents

\section ucbm_notation Notation

We use the symbols below in the documentation of all the classes of the
module, and the third column of each table gives the name of the netCDF
variable (or dimension), or the method, that provides each datum or variable.
Time is the first index of each time-dependent quantity, although the netCDF
variables may order their dimensions differently (as the documentation of each
deserialize() states, e.g., "ActivePowerDemand" is indexed [ node , instant
]). Upper-case Latin letters with a superscript denote data, lower-case Latin
letters variables or cost coefficients; \f$ c \f$ with a superscript is always
a cost, an upper-case \f$ C \f$ with a superscript a dimensionless rate or
scaling constant, calligraphic letters denote sets.

<table>
<caption>Sets and indices</caption>
<tr><th>Symbol</th><th>Meaning</th><th>Data or method</th></tr>
<tr><td>\f$ t \in \mathcal{T} = \{ 0 , \ldots , T - 1 \} \f$</td>
    <td>the instants of the horizon, the state before the horizon having
        index \f$ -1 \f$</td><td>"TimeHorizon"</td></tr>
<tr><td>\f$ n \in \mathcal{N} \f$, \f$ N = | \mathcal{N} | \f$</td>
    <td>the nodes of the network</td><td>"NumberNodes"</td></tr>
<tr><td>\f$ l \in \mathcal{L} \f$; \f$ s( l ) \f$, \f$ e( l ) \f$
        (\f$ e_j( l ) \f$ for a line with several ends)</td>
    <td>the lines, their start and end nodes</td>
    <td>"NumberLines", "StartLine", "EndLine"</td></tr>
<tr><td>\f$ \mathcal{L}^{S} \f$, \f$ \mathcal{L}^{H} \f$</td>
    <td>the lines with nonzero susceptance, and those with zero susceptance
        (HVDC lines, whose flow is free within its bounds)</td>
    <td>DCNetworkBlock</td></tr>
<tr><td>\f$ i \in \mathcal{I} \f$</td>
    <td>the UnitBlock of the UCBlock</td>
    <td>"NumberUnits", "UnitBlock_i"</td></tr>
<tr><td>\f$ g \in \mathcal{G} \f$, \f$ \mathcal{G}_n \f$, \f$ i( g ) \f$
        </td>
    <td>the electrical generators, those at node \f$ n \f$, the UnitBlock
        of \f$ g \f$</td>
    <td>"NumberElectricalGenerators", "GeneratorNode"</td></tr>
<tr><td>\f$ k \in \mathcal{K}_n \f$</td>
    <td>the storages of the units at node \f$ n \f$ (those of a unit being
        at the node of its first generator)</td>
    <td>"NumberStorages", UnitBlock::get_number_storages()</td></tr>
<tr><td>\f$ z \in \mathcal{Z}^{pr} \f$, \f$ \mathcal{Z}^{sc} \f$,
        \f$ \mathcal{Z}^{in} \f$</td>
    <td>the zones of the primary reserve, of the secondary reserve and of
        the inertia, each zone a set of nodes</td>
    <td>"NumberPrimaryZones", "PrimaryZones", and so on</td></tr>
<tr><td>\f$ p \in \mathcal{P} \f$, \f$ z \in \mathcal{Z}^{p} \f$</td>
    <td>the pollutants, the zones of pollutant \f$ p \f$</td>
    <td>"NumberPollutants", "NumberPollutantZones", "PollutantZones"</td>
    </tr>
<tr><td>\f$ n \in \mathcal{N}^{hy} \f$, \f$ l \in \mathcal{L}^{hy} \f$,
        \f$ j \in \mathcal{J}_l \f$</td>
    <td>the reservoirs and the arcs of a hydro valley, the pieces of the
        power curve of arc \f$ l \f$</td>
    <td>HydroUnitBlock</td></tr>
<tr><td>\f$ j \in \mathcal{J} \f$</td>
    <td>the rows (cuts) of a polyhedral function</td>
    <td>HydroSystemUnitBlock, PolyhedralFunction</td></tr>
</table>

<table>
<caption>Data and variables of the linking constraints</caption>
<tr><th>Symbol</th><th>Meaning</th><th>Data or method</th></tr>
<tr><td>\f$ \sigma_i \f$, \f$ \sigma_g := \sigma_{i(g)} \f$</td>
    <td>the scale of a UnitBlock (number of identical copies)</td>
    <td>"Scale", UnitBlock::get_scale()</td></tr>
<tr><td>\f$ p^{ac}_{t,g} \f$, \f$ p^{rc}_{t,g} \f$</td>
    <td>active and reactive power of a generator</td>
    <td>UnitBlock::get_active_power(), UnitBlock::get_reactive_power()</td>
    </tr>
<tr><td>\f$ p^{pr}_{t,g} \f$, \f$ p^{sc}_{t,g} \f$</td>
    <td>primary and secondary reserve of a generator</td>
    <td>UnitBlock::get_primary_spinning_reserve(),
        UnitBlock::get_secondary_spinning_reserve()</td></tr>
<tr><td>\f$ u_{t,g} \in \{ 0 , 1 \} \f$</td>
    <td>commitment of a generator (1 if on)</td>
    <td>UnitBlock::get_commitment()</td></tr>
<tr><td>\f$ P^{au}_{t,g} \geq 0 \f$</td>
    <td>power absorbed by a generator when it is off</td>
    <td>"FixedConsumption", UnitBlock::get_fixed_consumption()</td></tr>
<tr><td>\f$ h^{u}_{t,g} \f$, \f$ h^{p}_{t,g} \f$</td>
    <td>inertia given by a generator per unit of commitment and per unit of
        active power</td>
    <td>UnitBlock::get_inertia_commitment(),
        UnitBlock::get_inertia_power()</td></tr>
<tr><td>\f$ v^{sl}_{t,k} \f$</td>
    <td>level of storage \f$ k \f$ at the end of instant \f$ t \f$</td>
    <td>UnitBlock::get_storage_level()</td></tr>
<tr><td>\f$ S_{t,n} \f$, \f$ R_{t,n} \f$</td>
    <td>active and reactive injection at a node</td>
    <td>NetworkBlock::get_node_injection(),
        NetworkBlock::get_reactive_node_injection()</td></tr>
<tr><td>\f$ D^{ac}_{t,n} \f$, \f$ D^{rc}_{t,n} \f$</td>
    <td>active and reactive demand of a node</td>
    <td>"ActivePowerDemand", "ReactivePowerDemand", "ActiveDemand"</td>
    </tr>
<tr><td>\f$ D^{pr}_{t,z} \f$, \f$ D^{sc}_{t,z} \f$, \f$ D^{in}_{t,z} \f$
        </td>
    <td>primary reserve, secondary reserve and inertia required in a zone
        </td>
    <td>"PrimaryDemand", "SecondaryDemand", "InertiaDemand"</td></tr>
<tr><td>\f$ \rho_{t,p,g} \f$, \f$ \rho^{v}_{t,p,k} \f$</td>
    <td>emission of pollutant \f$ p \f$ per unit of active power over one
        instant, and factor of a storage level</td>
    <td>"PollutantRho", "PollutantStorageRho"</td></tr>
<tr><td>\f$ O_{z,p} \f$, \f$ O^{mn}_{z,p} \f$</td>
    <td>upper and lower bound on the emission of a zone over the horizon
        </td>
    <td>"PollutantBudget", "PollutantMinBudget"</td></tr>
<tr><td>\f$ c^{0}_{k} \f$, \f$ c^{U} \f$, \f$ c^{net}_{l} \f$,
        \f$ V_{k,l} \f$</td>
    <td>constant term of NetworkBlock \f$ k \f$, constant term of the
        UCBlock (that of the NetworkBlock in the bus case), cost of the flow
        of a line, auxiliary variable bounding the absolute value of the
        flow</td>
    <td>"ConstantTerm", "NetworkConstantTerms", "NetworkCost"</td></tr>
<tr><td>\f$ y^{ac}_{t,n} \f$, \f$ y^{pr}_{t,z} \f$, \f$ y^{sc}_{t,z} \f$,
        \f$ y^{in}_{t,z} \f$, \f$ y^{p}_{z} \f$</td>
    <td>dual values (multipliers) of the linking constraints (1)-(5)</td>
    <td>"ActivePowerDuals", "PrimaryDuals", "SecondaryDuals",
        "InertiaDuals", "PollutantDuals"</td></tr>
<tr><td>\f$ y^{rc}_{t,n} \f$</td>
    <td>dual value of the reactive node injection constraint (1')</td>
    <td>RowConstraint::get_dual() (not written by UCBlockSolution)</td></tr>
<tr><td>\f$ \Theta( y ) \f$, \f$ \Theta_i \f$, \f$ \Theta^{net}_k \f$</td>
    <td>the dual function and its components</td>
    <td>LagrangianDualSolver</td></tr>
</table>

<table>
<caption>Network (DCNetworkBlock)</caption>
<tr><th>Symbol</th><th>Meaning</th><th>Data or method</th></tr>
<tr><td>\f$ F_l \f$</td><td>flow of a line</td><td>"FlowValue"</td></tr>
<tr><td>\f$ \mathfrak{S}_l \f$, \f$ \theta_n \f$</td>
    <td>susceptance of a line, voltage angle of a node</td>
    <td>"LineSusceptance"</td></tr>
<tr><td>\f$ \eta_{l} \f$, \f$ \eta_{l,j} \f$</td>
    <td>efficiency of an HVDC line (of its \f$ j \f$-th end)</td>
    <td>"Efficiency"</td></tr>
<tr><td>\f$ P^{mn}_l \f$, \f$ P^{mx}_l \f$, \f$ \kappa_l \f$,
        \f$ C^{v} \f$</td>
    <td>flow bounds, capacity multiplier of a line (1 unless changed, it
        is not read from the data), scaling of the flow bounds</td>
    <td>"MinPowerFlow", "MaxPowerFlow", DCNetworkBlock::set_kappa(), the
        Configuration of the static constraints</td></tr>
<tr><td>\f$ \Psi \f$, \f$ \Psi_{ln} \f$</td>
    <td>the matrix of the Power Transfer Distribution Factors (PTDF), with
        no column for the reference nodes \f$ n \in R \f$</td>
    <td>DCNetworkBlock::DCNetworkData::get_PTDF()</td></tr>
<tr><td>\f$ y^{ov}_t \f$, \f$ y^{b}_{t,n} \f$, \f$ \mu_{t,l} \f$</td>
    <td>dual values of the overall balance, of the node balance of node
        \f$ n \f$ and of the flow bounds of line \f$ l \f$ in the network
        of instant \f$ t \f$ (\f$ \mu_{t,l} \geq 0 \f$ if the upper bound is
        active, \f$ \leq 0 \f$ if the lower one is), and
        \f$ \mu^{mx}_{t,l} , \mu^{mn}_{t,l} \geq 0 \f$ its positive and
        negative parts</td>
    <td>DCNetworkBlock::get_dual_prices()</td></tr>
</table>

<table>
<caption>Units (see each UnitBlock for the complete list)</caption>
<tr><th>Symbol</th><th>Meaning</th><th>Data or method</th></tr>
<tr><td>\f$ P^{mn}_t \f$, \f$ P^{mx}_t \f$</td>
    <td>minimum and maximum power of a unit in its rows (for a
        ThermalUnitBlock the operational ones, with the availability
        \f$ \chi_t \f$)</td>
    <td>"MinPower", "MaxPower", "Availability"</td></tr>
<tr><td>\f$ \hat P^{mn}_t \f$, \f$ \hat P^{mx}_t \f$</td>
    <td>nominal minimum and maximum power of a unit, without the
        availability</td>
    <td>"MinPower", "MaxPower", UnitBlock::get_min_power(),
        UnitBlock::get_max_power()</td></tr>
<tr><td>\f$ P^{su}_t \f$, \f$ P^{sd}_t \f$</td>
    <td>start-up and shut-down limits</td>
    <td>"StartUpLimit", "ShutDownLimit"</td></tr>
<tr><td>\f$ \Delta^{+}_t \f$, \f$ \Delta^{-}_t \f$</td>
    <td>largest increase and decrease of the power from instant
        \f$ t - 1 \f$ to instant \f$ t \f$</td>
    <td>"DeltaRampUp", "DeltaRampDown"</td></tr>
<tr><td>\f$ \tau^{+} \f$, \f$ \tau^{-} \f$, \f$ \tau_0 \f$</td>
    <td>minimum up and down times, and the time the unit has been up
        (\f$ > 0 \f$) or down (\f$ \leq 0 \f$) before the horizon</td>
    <td>"MinUpTime", "MinDownTime", "InitUpDownTime"</td></tr>
<tr><td>\f$ p_{-1} \f$</td><td>power before the horizon</td>
    <td>"InitialPower"</td></tr>
<tr><td>\f$ v_t \f$, \f$ w_t \f$</td>
    <td>start-up and shut-down of a thermal unit</td>
    <td>ThermalUnitBlock</td></tr>
<tr><td>\f$ a_t \f$, \f$ b_t \f$, \f$ c_t \f$, \f$ c^{su}_t \f$,
        \f$ c^{sd}_t \f$</td>
    <td>quadratic, linear and fixed cost, start-up and shut-down cost</td>
    <td>"QuadTerm", "LinearTerm", "ConstTerm", "StartUpCost",
        "ShutDownCost"</td></tr>
<tr><td>\f$ c^{pr}_t \f$, \f$ c^{sc}_t \f$, \f$ \rho^{pr}_t \f$,
        \f$ \rho^{sc}_t \f$</td>
    <td>costs of the reserves (0 if not given), and the largest reserves
        as fractions of the active power</td>
    <td>"PrimarySpinningReserveCost", "SecondarySpinningReserveCost",
        "PrimaryRho", "SecondaryRho"</td></tr>
<tr><td>\f$ x \f$, \f$ x_b \f$, \f$ x_c \f$, \f$ x_l \f$,
        \f$ c^{inv} \f$</td>
    <td>design variable of a unit (of a battery and of its converter, of a
        line) and its investment cost</td>
    <td>"InvestmentCost"</td></tr>
<tr><td>\f$ \underline{X} \f$, \f$ \bar X \f$</td>
    <td>bounds of a design variable, integer in
        \f$ \{ \underline{X} , \ldots , | \bar X | \} \f$ if
        \f$ \bar X < 0 \f$</td>
    <td>"MinCapacityDesign", "MaxCapacityDesign" (with the prefix
        "Battery" or "Converter" for a battery)</td></tr>
<tr><td>\f$ H_t \f$</td>
    <td>inertia constant of a synchronous machine, folded in
        \f$ h^u_{t,g} \f$ or \f$ h^p_{t,g} \f$</td>
    <td>"InertiaCommitment", "InertiaPower"</td></tr>
<tr><td>\f$ \tau^{M} \f$, \f$ L^{M} \f$, \f$ \Delta^{M+}_t \f$,
        \f$ \Delta^{M-}_t \f$, \f$ \tau^{v} \f$, \f$ N^{su} \f$,
        \f$ N^{dd} \f$</td>
    <td>the rules of a modulation of a nuclear unit</td>
    <td>NuclearUnitBlock</td></tr>
<tr><td>\f$ p^{+}_t \f$, \f$ p^{-}_t \f$, \f$ v^{ba}_t \f$,
        \f$ u^{ch}_t \f$</td>
    <td>power delivered and drawn by a battery, its level at the end of
        instant \f$ t \f$, the charging indicator</td>
    <td>BatteryUnitBlock</td></tr>
<tr><td>\f$ \rho^{in}_t \f$, \f$ \rho^{out}_t \f$, \f$ \rho^{st}_t \f$
        </td>
    <td>efficiencies of charge and discharge, standing loss</td>
    <td>"StoringBatteryRho", "ExtractingBatteryRho",
        "StandingBatteryRho"</td></tr>
<tr><td>\f$ v^{hy}_{n,t} \f$, \f$ f_{t,l} \f$, \f$ A_{n,t} \f$</td>
    <td>volume of a reservoir, flow of an arc, inflow, all per instant
        </td>
    <td>HydroUnitBlock</td></tr>
<tr><td>\f$ \kappa \f$, \f$ \gamma \f$</td>
    <td>capacity multiplier and reserve margin of an intermittent unit
        (\f$ \kappa \f$ also of a battery)</td>
    <td>"Kappa", "Gamma"</td></tr>
</table>

<table>
<caption>Investment (InvestmentBlock)</caption>
<tr><th>Symbol</th><th>Meaning</th><th>Data or method</th></tr>
<tr><td>\f$ x_i \f$, \f$ \bar x_i \f$, \f$ i = 1 , \ldots , m \f$</td>
    <td>the investment in asset \f$ i \f$ and its installed quantity</td>
    <td>the ColVariable of InvestmentBlock, "InstalledQuantity"</td></tr>
<tr><td>\f$ c_i \f$, \f$ d_i \f$</td>
    <td>cost of investment and of disinvestment</td>
    <td>"Cost", "DisinvestmentCost"</td></tr>
<tr><td>\f$ \mathcal{F}( x ) \f$, \f$ Q( x ) \f$</td>
    <td>the investment function and the operational cost as a function of
        the investment</td>
    <td>InvestmentFunction</td></tr>
<tr><td>\f$ f^1_i \f$, \f$ g^1_i \f$, \f$ \bar z_i \f$</td>
    <td>the cost of one copy of a scaled unit, its contributions to the
        linking constraints, and a solution of it</td>
    <td>UnitBlock::scale()</td></tr>
<tr><td>\f$ \upsilon \in \Upsilon \f$</td>
    <td>alternative models of the uncertainty (meta-scenarios), not
        modeled</td>
    <td>none</td></tr>
</table>

<table>
<caption>Multistage models (SDDPBlock)</caption>
<tr><th>Symbol</th><th>Meaning</th><th>Data or method</th></tr>
<tr><td>\f$ v^{f} \f$</td>
    <td>the vector of the final volumes \f$ v^{hy}_{n,T-1} \f$ of the
        reservoirs of a HydroSystemUnitBlock</td>
    <td>HydroUnitBlock::get_volume()</td></tr>
<tr><td>\f$ \check\nu( v^{f} ) \f$, \f$ \alpha_j \f$, \f$ \beta_j \f$</td>
    <td>the future cost of the water and the constant and the coefficients
        of its row (cut) \f$ j \in \mathcal{J} \f$</td>
    <td>PolyhedralFunction</td></tr>
<tr><td>\f$ \mathfrak{A} \f$</td>
    <td>a matrix aggregating volumes, not modeled</td>
    <td>none</td></tr>
</table>

\section ucbm_conv Conventions

\subsection ucbm_conv_time Instants and data per instant

We denote by \f$ \mathcal{T} = \{ 0 , \ldots , T - 1 \} \f$ the set of
instants of the time horizon, where \f$ T \f$ is "TimeHorizon", and we index
by \f$ -1 \f$ the state of the system before the horizon (e.g., \f$ p_{-1} \f$
is the power of a thermal unit at the instant before the horizon, and
\f$ v^{hy}_{n,-1} \f$ the volume of a reservoir). Note that the instants need
not be equally spaced, since no length of the time step appears in any
constraint or in the objective: each datum is expressed per instant, and the
conversion of data given per hour is the task of whoever writes the data. With
instants of \f$ \Delta t \f$ hours, a ramp is the largest change of a power,
or of a flow, between two consecutive instants, and therefore a gradient of
\f$ G \f$ MW/h becomes \f$ \Delta^{+} = G \Delta t \f$. A minimum up or down
time, a delay and a stability period are numbers of instants (a duration of
\f$ d \f$ hours becomes \f$ \lceil d / \Delta t \rceil \f$ instants), while
the flows and the inflows of a hydro valley are volumes per instant, i.e.,
\f$ 3600 \, \Delta t \f$ times a value in m\f$ ^3 \f$/s. Since each cost
coefficient is the cost of one instant, a cost per MWh (or per hour) is
multiplied by \f$ \Delta t \f$, and so is an emission factor given per MWh.
The factor \f$ \rho_{t,p,g} \f$ is thus in tonnes per MW over one instant,
possibly divided by the efficiency of the generator if the emission is that of
the fuel it burns. A power, a reserve and a demand are instead powers (e.g.,
MW) at the instant, and the level of a storage is an energy (e.g., MWh), or a
volume of water for a reservoir. Hence, the coefficients that turn a power
into a change of the level over one instant (e.g., the efficiencies of a
battery, see BatteryUnitBlock) include \f$ \Delta t \f$, and with all of them
equal to 1 the level is measured in units of power times one instant. However,
the inertia is not multiplied by \f$ \Delta t \f$, since it is an energy
(e.g., MWs) available at the instant. Also, the dual values of the constraints
are values per instant (see \ref ucbm_dual_sign).

\subsection ucbm_conv_sign Signs

We take the active power \f$ p^{ac}_{t,g} \f$ of a generator positive when the
generator injects power into the network and negative when it absorbs it (as a
pump, a charging battery or a load does). A generator that is off absorbs
instead the fixed consumption \f$ P^{au}_{t,g} \geq 0 \f$, which therefore
enters the balance of its node with a minus sign. Both reserves
\f$ p^{pr}_{t,g} \f$ and \f$ p^{sc}_{t,g} \f$ are nonnegative and symmetric,
i.e., a unit offering them must be able to move its power by that amount both
upward and downward, which each UnitBlock enforces through its constraints. A
positive flow \f$ F_l \f$ goes from \f$ s( l ) \f$ to \f$ e( l ) \f$. By the
node injection \f$ S_{t,n} \f$ we mean the net power that the units of node
\f$ n \f$ give to the network, positive when they inject. Finally, the problem
is a minimization.

\subsection ucbm_conv_scale Scaling factors

Four kinds of factors multiply the data or the variables of the sub-Blocks.
First, the scale \f$ \sigma_i \geq 0 \f$ of a UnitBlock multiplies its cost
and its terms in the linking constraints (the Variable of the UnitBlock are
those of one copy of the unit). Second, the capacity multipliers
\f$ \kappa \f$ of an IntermittentUnitBlock or of a BatteryUnitBlock, and
\f$ \kappa_l \f$ of a line, multiply the power, storage and flow limits (see
\ref ucbm_cap). Third, DCNetworkBlock applies a factor \f$ C^{v} \f$ (1 by
default) to all the flow bounds. Finally, some factors are folded in the data,
as the length of the instant in the emission factors and the factor
\f$ 1.2 \f$ in the inertia (see \ref ucbm_link_in).

\section ucbm_problem The problem

\subsection ucbm_sets Sets and sub-Blocks

In a UCBlock, the units are partitioned into the UnitBlock sub-Blocks, each of
which is a set of interdependent units (most often one; a HydroSystemUnitBlock
groups the HydroUnitBlock of a valley) whose technical constraints and cost
the UnitBlock handles entirely. Each UnitBlock \f$ i \in \mathcal{I} \f$ has
one or more electrical generators \f$ g \in \mathcal{G} \f$, numbered unit
after unit, each one connected to one node \f$ n \in \mathcal{N} \f$
("GeneratorNode"), and therefore the generators of one unit may be at
different nodes (e.g., the plants of a valley). On the other hand, the network
constraints of each instant (definition of the flows, flow limits, node
balances) belong to the NetworkBlock sub-Blocks, one for each instant or group
of consecutive instants, whose only interface with the units is the vector of
the node injections \f$ S_{t,n} \f$. When the network has a single node (a
bus) there is no NetworkBlock (see UCBlock::deserialize()).

\subsection ucbm_obj The objective

We minimize the sum of the objectives of the sub-Blocks plus that of the
UCBlock, which is a constant (cf. UCBlock::generate_objective()):
\f[
  \min \; \sum_{ i \in \mathcal{I} } f_i
     + \sum_{ k } \Bigl( c^{0}_k + \sum_{ l \in \mathcal{L} } c^{net}_l
       V_{k,l} \Bigr) + c^{U} \tag{0}
\f]
where \f$ f_i \f$ is the cost of UnitBlock \f$ i \f$, already multiplied by
\f$ \sigma_i \f$, which for a HydroSystemUnitBlock includes the future cost of
the water left in its reservoirs (a convex polyhedral function of the final
volumes, see HydroSystemUnitBlock). In the second sum, \f$ k \f$ runs over the
NetworkBlock, \f$ c^0_k \f$ is its constant term and \f$ c^{net}_l V_{k,l} \f$
is the optional cost of the flow on line \f$ l \f$
(\f$ V_{k,l} \geq | F_l | \f$, see DCNetworkBlock). The constant \f$ c^{U} \f$
is the sum of the constant terms of the NetworkBlock in the bus case, where
the NetworkBlock are not built (zero otherwise). For each NetworkBlock, the
constant term is its "ConstantTerm" if its group has one, otherwise the entry
of "NetworkConstantTerms" of the UCBlock, if given, otherwise zero.

\subsection ucbm_link The linking constraints

The feasible set is the product of the feasible sets of the sub-Blocks
intersected with the linking constraints of
UCBlock::generate_abstract_constraints(), i.e., the node injection rows (1)
(and the reactive ones (1') when the network handles reactive power), the
primary and secondary reserve requirements (2) and (3), the inertia
requirements (4) and the pollutant budgets (5). These are the only constraints
of the UCBlock, which has no Variable. In all of them, each term of a
generator, or of a storage, is multiplied by the scale of its unit, and a term
is absent when the unit does not have the corresponding Variable or datum
(e.g., a unit without commitment variables, which is always on, has no fixed
consumption).

\subsubsection ucbm_link_ac Node injection and demand

For \f$ t \in \mathcal{T} \f$ and \f$ n \in \mathcal{N} \f$, with more
than one node,
\f[
  \sum_{ g \in \mathcal{G}_n } \sigma_g \bigl( p^{ac}_{t,g}
    - P^{au}_{t,g} ( 1 - u_{t,g} ) \bigr) - S_{t,n} = 0 \; ,
  \tag{1}
\f]
while in the bus case, where there is no node injection,
\f[
  \sum_{ g \in \mathcal{G} } \sigma_g \bigl( p^{ac}_{t,g}
    - P^{au}_{t,g} ( 1 - u_{t,g} ) \bigr) = D^{ac}_{t} \tag{1b}
\f]
is directly the balance between the total production and the demand. With
more than one node, the node balance (B) of the NetworkBlock of
\f$ t \f$ in DCNetworkBlock reads
\f[
  - S_{t,n} + \sum_{ l : s( l ) = n } F_l
    - \sum_{ l , j : e_j( l ) = n } \eta_{l,j} F_l = - D^{ac}_{t,n}
\f]
where \f$ \eta_{l,j} = 1 \f$ for the lines with susceptance (the rows actually
written depend on the formulation chosen in DCNetworkBlock, although all of
them have the same set of feasible injections). Composed with (1), it gives
the usual balance of active power at node \f$ n \f$:
\f[
  \sum_{ g \in \mathcal{G}_n } \sigma_g \bigl( p^{ac}_{t,g}
    - P^{au}_{t,g} ( 1 - u_{t,g} ) \bigr) - D^{ac}_{t,n}
  = \sum_{ l : s( l ) = n } F_l - \sum_{ l , j : e_j( l ) = n }
    \eta_{l,j} F_l \; , \tag{6}
\f]
i.e., the net production at the node minus its demand equals the net flow
leaving it. Note that the only differences from the lossless balance are the
fixed consumption of the units that are off, the scale of the units and the
efficiencies of the HVDC lines. We intend the demand of a node to include the
reference consumption profile of any load-shifting flexibility connected to
it, whose deviation from that profile is the active power of the
BatteryUnitBlock that represents it. Nodes may also represent zones of a
larger grid whose lines have been aggregated (see DCNetworkBlock), or a
transmission node together with the distribution grid below it, which makes no
difference to the UCBlock. When the network handles reactive power the same
rows are written for it,
\f[
  \sum_{ g \in \mathcal{G}_n } \sigma_g \, p^{rc}_{t,g} - R_{t,n} = 0
  \; , \qquad \text{or} \qquad
  \sum_{ g \in \mathcal{G} } \sigma_g \, p^{rc}_{t,g} = D^{rc}_{t}
  \text{ in the bus case,} \tag{1'}
\f]
without fixed consumption, which is an active power only.

\subsubsection ucbm_link_zones Zones

Three families of zones are given, i.e., \f$ \mathcal{Z}^{pr} \f$ for the
primary reserve, \f$ \mathcal{Z}^{sc} \f$ for the secondary reserve and
\f$ \mathcal{Z}^{in} \f$ for the inertia; each zone is a set of nodes, the
zones of a family are pairwise disjoint, and each zone \f$ z \f$ carries one
requirement per instant. The families may coincide, and they need not cover
\f$ \mathcal{N} \f$; indeed, a node whose zone index (say, "PrimaryZones"[ n
]) is not smaller than the number of zones of the family is in no zone of it,
and its generators do not enter the corresponding constraints. This holds also
when the family has a single zone; without the vector of the zones, which is
then optional, all the nodes are in the single zone, while more zones require
the vector. In the same way, for each pollutant \f$ p \in \mathcal{P} \f$ a
family \f$ \mathcal{Z}^{p} \f$ of disjoint zones carries the emission budgets.
Active power needs no zones, since its demand is given node by node. In what
follows, we denote by \f$ z^{pr}( g ) \f$ (and so on) the zone of the node of
generator \f$ g \f$, if any.

\subsubsection ucbm_link_res Reserves

For \f$ t \in \mathcal{T} \f$, \f$ z \in \mathcal{Z}^{pr} \f$ and
\f$ z' \in \mathcal{Z}^{sc} \f$,
\f[
  \sum_{ n \in z } \sum_{ g \in \mathcal{G}_n } \sigma_g \, p^{pr}_{t,g}
    \geq D^{pr}_{t,z} \tag{2}
\f]
\f[
  \sum_{ n \in z' } \sum_{ g \in \mathcal{G}_n } \sigma_g \, p^{sc}_{t,g}
    \geq D^{sc}_{t,z'} \tag{3}
\f]
where only the generators whose UnitBlock has reserve variables appear.
Indeed, the UCBlock requires from each UnitBlock the reserve variables of the
kinds that have at least one zone (see UnitBlock::set_reserve_vars()), and
each UnitBlock bounds them by its technical constraints.

\subsubsection ucbm_link_in Inertia

The inertia that generator \f$ g \f$ provides at instant \f$ t \f$ is the
affine function \f$ h^u_{t,g} u_{t,g} + h^p_{t,g} p^{ac}_{t,g} \f$, where each
term is absent if the corresponding datum or Variable is. For a synchronous
machine of inertia constant \f$ H_t \f$ and rated power \f$ \hat P^{mx}_t \f$
the contribution is usually \f$ 1.2 H_t \hat P^{mx}_t \f$ whenever it is on,
i.e., \f$ h^u_{t,g} = 1.2 H_t \hat P^{mx}_t \f$ and \f$ h^p_{t,g} = 0 \f$
(ThermalUnitBlock "InertiaCommitment", SlackUnitBlock "MaxInertia"). For a
hydro turbine or an intermittent unit it is instead taken proportional to the
power produced, i.e., \f$ h^u_{t,g} = 0 \f$ and \f$ h^p_{t,g} = 1.2 H_t \f$
("InertiaPower"). Note that the factor 1.2, which converts active into
apparent power under a constant phase angle, belongs to the data, which the
UCBlock uses as they are, and that a BatteryUnitBlock provides no inertia. For
\f$ t \in \mathcal{T} \f$ and \f$ z \in \mathcal{Z}^{in} \f$
\f[
  \sum_{ n \in z } \sum_{ g \in \mathcal{G}_n } \sigma_g \bigl(
    h^u_{t,g} u_{t,g} + h^p_{t,g} p^{ac}_{t,g} \bigr)
    \geq D^{in}_{t,z} \; . \tag{4}
\f]

\subsubsection ucbm_link_pol Pollutant budgets

For \f$ p \in \mathcal{P} \f$ and \f$ z \in \mathcal{Z}^{p} \f$
\f[
  O^{mn}_{z,p} \leq \sum_{ t \in \mathcal{T} } \sum_{ n \in z } \Bigl(
    \sum_{ g \in \mathcal{G}_n } \sigma_g \, \rho_{t,p,g} \, p^{ac}_{t,g}
    + \sum_{ k \in \mathcal{K}_n } \sigma_k \, \rho^{v}_{t,p,k} \,
      v^{sl}_{t,k} \Bigr) \leq O_{z,p} \tag{5}
\f]
where \f$ \sigma_k \f$ is the scale of the unit of storage \f$ k \f$,
\f$ O^{mn}_{z,p} = - \infty \f$ if "PollutantMinBudget" is not given, and
\f$ O^{mn}_{z,p} = O_{z,p} \f$ is allowed (a limit the emission has to match).
With the term of the storage levels we account for the change over the horizon
of the level of a storage that holds a pollutant, or a fuel that emits it,
typically with a nonzero factor at the last instant only; the contribution of
its initial level, a constant, is then to be subtracted from both bounds by
whoever writes the data. Since the module has no units producing heat, no
emission due to heat is accounted for.

\subsubsection ucbm_link_bounds Bounds on the injections

The UCBlock gives to each NetworkBlock the bounds
\f[
  \sum_{ g \in \mathcal{G}_n } \sigma_g \kappa_g X_g
    \min\{ \hat P^{mn}_{t,g} , - P^{au}_{t,g} \} \leq S_{t,n} \leq
  \sum_{ g \in \mathcal{G}_n } \sigma_g \kappa_g X_g
    \max\{ 0 , \hat P^{mx}_{t,g} \}
\f]
on the injections, with the nominal bounds of the units, i.e., those of
UnitBlock::get_min_power() and UnitBlock::get_max_power() without the
availability of a ThermalUnitBlock (which makes them valid but possibly weaker
than those of the rows). Here \f$ X_g \f$ is the largest value of the design
variable of the unit of \f$ g \f$ that bounds its power, i.e., \f$ | \bar X |
\f$ ("MaxCapacityDesign", or "BatteryMaxCapacityDesign" for a battery; 1 if the
unit has none), and \f$ \kappa_g \f$ is the kappa of an IntermittentUnitBlock
or of a BatteryUnitBlock, 1 for the other units (see UnitBlock::get_design_ub()
and UnitBlock::get_kappa()). The same holds on the reactive injections
\f$ R_{t,n} \f$ with the smallest and the largest reactive power of the unit,
on or off (for a ThermalUnitBlock, \f$ Q^{mn}_t + \min\{ 0 , Q^{mn,on}_t \} \f$
and \f$ Q^{mx}_t + \max\{ 0 , Q^{mx,on}_t \} \f$). All of these bounds are
implied by the constraints of the units, and they are computed when the
UCBlock is read and again whenever the scale or the kappa of a unit, or the
maximum power of a ThermalUnitBlock or of an IntermittentUnitBlock, changes
[see UCBlock::set_node_injection_bounds()]. Only ECNetworkBlock imposes the
active ones as constraints, and ACNetworkBlock the reactive ones, while
DCNetworkBlock bounds the injections through the flow limits.

\section ucbm_units The sub-Blocks

Each kind of sub-Block contributes to (0)-(5) as follows; the complete
model of each one is in the documentation of its class, to which we refer.

- ThermalUnitBlock: one generator with commitment \f$ u_t \f$, start-up
  \f$ v_t \f$ and shut-down \f$ w_t \f$, power \f$ p^{ac}_t \in \{ 0 \} \cup [
  P^{mn}_t , P^{mx}_t ] \f$, ramps, start-up and shut-down limits, minimum up
  and down times and the initial state. It gives to (1) its power and its
  fixed consumption \f$ P^{au}_t \f$ when off, to (2)-(3) its primary and
  secondary reserves, to (4) the term \f$ h^u_t u_t \f$, and to (5) its
  emissions. Its reserves are bounded by fractions of its power and by the
  room that the power bounds and the ramps leave, at instant 0 the ramps from
  the initial power of a unit on before the horizon (see ThermalUnitBlock for
  the exact rows). Its cost is the sum of a convex quadratic function of the
  power, of the fixed, start-up and shut-down costs, of the costs of the
  reserves (0 when the data do not give them), possibly of an investment term,
  of the deviation from a reference schedule and of a price of the reactive
  power, all multiplied by \f$ \sigma \f$.
- NuclearUnitBlock: a ThermalUnitBlock whose power changes by modulations
  between bands, with the rules on their length, stability and number per
  day, and the deep decreases of the power.
- HydroUnitBlock: a valley of reservoirs joined by arcs, one generator per
  arc, whose power is bounded by a concave piecewise-linear function of the
  flow (equal to it in the cases that HydroUnitBlock states, and to a linear
  one for a pump, whose power is negative). The volumes \f$ v^{hy}_{n,t} \f$
  are its storages for (5), and its arcs give power, reserves and \f$ h^p_t
  p^{ac}_t \f$ to (1)-(4).
- HydroSystemUnitBlock: a set of HydroUnitBlock solved together, plus the
  future cost of the water left in all their reservoirs at the end of the
  horizon, a convex polyhedral function of the final volumes.
- BatteryUnitBlock: a storage whose active power \f$ p^{ac}_t = p^{+}_t -
  p^{-}_t \f$ is the power delivered minus the power drawn, with efficiencies
  and standing loss, reserves bounded by data and by the room that the power
  limits leave, the level \f$ v^{ba}_t \f$ at the end of each instant as its
  storage for (5), and no inertia. It also represents load-shifting
  flexibilities.
- IntermittentUnitBlock: a unit whose power is bounded by
  \f$ \kappa P^{mx}_t \f$, where "MaxPower" is the production per unit of
  installable capacity times that capacity, possibly with reserves bounded
  through the margin \f$ \gamma \f$, and an inertia \f$ h^p_t p^{ac}_t \f$
  proportional to the dispatched power.
- SlackUnitBlock: a unit of large cost covering an imbalance, whose power lies
  between 0 and "MaxPower" (which may be negative), whose default is 0 (no
  slack), possibly with reserves and with an inertia given through a
  commitment variable.
- DCNetworkBlock (and its variants): the flows of the lines, their bounds
  \f$ \kappa_l C^{v} P^{mn}_l \leq F_l \leq \kappa_l C^{v} P^{mx}_l \f$
  and the node balances, written in one of three equivalent formulations
  (PTDF, CYCLE, KIRCHHOFF); ECNetworkBlock models instead an energy
  community and ACNetworkBlock an alternating-current network with
  reactive power.

\section ucbm_cap Capacities and investment

\subsection ucbm_cap_param Scale and capacity multipliers

In the model, the capacity of the assets enters through two kinds of
parameters of the sub-Blocks. First, the scale \f$ \sigma_i \f$ of a UnitBlock
(ThermalUnitBlock, IntermittentUnitBlock, BatteryUnitBlock; see
UnitBlock::scale()) represents a number of identical copies of the unit,
committed and dispatched together. The Variable are those of one copy, and the
UnitBlock contributes \f$ \sigma_i f^1_i \f$ to (0) and \f$ \sigma_i g^1_i \f$
to (1)-(5), where \f$ f^1_i \f$ and \f$ g^1_i \f$ are the cost and the
contributions of one copy. Second, the capacity multipliers \f$ \kappa \f$ of
IntermittentUnitBlock and BatteryUnitBlock, and \f$ \kappa_l \f$ of the lines
of DCNetworkBlock, multiply the power, storage and flow limits of a unit or of
a line. For an intermittent unit, with "MaxPower" the production per unit of
installable capacity times that capacity, \f$ \kappa \in [ 0 , 1 ] \f$ is the
installed fraction; the initial level of a battery and the energy bounds of an
intermittent unit are, however, absolute quantities, not multiplied by
\f$ \kappa \f$. Both kinds of parameters can be changed after the Block is
built (UnitBlock::scale(), IntermittentUnitBlock::set_kappa(),
BatteryUnitBlock::set_kappa(), DCNetworkBlock::set_kappa()), and the UCBlock
rewrites the coefficients of (1)-(5) when a scale changes. Therefore, one can
evaluate the optimal cost as a function of the capacities.

\subsection ucbm_cap_design Design variables

Alternatively, the capacities may be variables of the model. A
ThermalUnitBlock with "InvestmentCost" has a binary design variable \f$ x \f$
with \f$ u_t \leq x \f$, while an IntermittentUnitBlock with "InvestmentCost"
has a design variable \f$ x \f$ (continuous or integer between
"MinCapacityDesign" and "MaxCapacityDesign") that multiplies its power bounds.
A BatteryUnitBlock has design variables \f$ x_b \f$ and \f$ x_c \f$ for the
storage and the converter, and a DesignNetworkBlock adds a design variable
\f$ x_l \f$ that multiplies the flow bounds of a line. In this case the
investment cost \f$ \sigma c^{inv} x \f$ (and its analogues) is part of the
cost of the sub-Block, and the problem (investment included) is a single
mixed-integer program.

\subsection ucbm_cap_outer The investment problem

In the decomposed approach the investment is the variable of an outer
problem, an InvestmentBlock whose objective is an InvestmentFunction. Let
\f$ x \in \mathbb{R}^m \f$ be the vector of its Variable, one per asset
\f$ i = 1 , \ldots , m \f$, \f$ \bar x_i \f$ the installed quantities and
\f$ c_i \f$, \f$ d_i \f$ the investment and disinvestment costs; the
function is
\f[
  \mathcal{F}( x ) = \sum_{ i = 1 }^m \bigl( c_i ( x_i - \bar x_i )^+
    + d_i ( \bar x_i - x_i )^+ \bigr) + Q( x ) \tag{7}
\f]
where \f$ Q( x ) \f$ is the operational cost computed by the Solver of the
inner Block (a UCBlock, or a stochastic Block whose scenarios are UCBlock)
once the investment has been written into it. For an asset represented by a
scale, \f$ x_i \f$ is the scale \f$ \sigma \f$ of the unit, i.e., its number
of copies; for one represented by a capacity multiplier, \f$ x_i = \kappa \f$
(or \f$ \kappa_l \f$, on the NetworkBlock of each instant). Note that the
investment term is convex if \f$ c_i + d_i \geq 0 \f$ for all \f$ i \f$. The
optimal value of a convex model of the problem (its continuous relaxation, or
its Lagrangian dual, see \ref ucbm_dual) is a convex function of the
capacities. Indeed, \f$ \sigma \f$ multiplies the complete contribution of a
unit (see UnitBlock::scale()), while the \f$ \kappa \f$ of an intermittent
unit and the \f$ \kappa_l \f$ of a line only appear in right-hand sides of
linear rows, and the optimal value of a convex problem is a convex function of
its right-hand sides. The \f$ \kappa \f$ of a battery also multiplies the
binary variables \f$ u^{ch}_t \f$ in the rows (8) of BatteryUnitBlock;
however, in their continuous relaxation the change of variable
\f$ w_t = \kappa u^{ch}_t \f$ gives rows that are linear in
\f$ ( p , w , \kappa ) \f$ jointly (see InvestmentFunction), whose optimal
value is again convex in \f$ \kappa \f$. None of this holds when a unit or a
line is also in design mode, where \f$ \kappa \f$ multiplies the design
variable. For the problem with integer commitments neither convexity nor the
subgradients below hold, and a fractional \f$ \sigma \f$ is the continuous
relaxation of the number of copies. With a stochastic inner Block
\f$ Q( x ) \f$ is the expected operational cost (with an SDDPBlock, the mean
cost of the simulated policy, see InvestmentFunction), and the linearization
is averaged in the same way. Since the bounds and linear constraints of the
InvestmentBlock give the domain of \f$ \mathcal{F} \f$, with the integrality
of the Variable relaxed minimizing \f$ \mathcal{F} \f$ over it is a convex
nondifferentiable problem, which a bundle method can solve using the values
and the subgradients below.

\subsection ucbm_cap_sens Sensitivities

Let the inner problem be solved as a convex problem, with optimal dual
values \f$ \bar y \f$ of its constraints, which follow the convention of
RowConstraint (see \ref ucbm_dual_sign): the derivative of the optimal
value with respect to the finite side \f$ b \f$ of a row is
\f$ - \bar y \f$. For a unit represented by a scale, the Lagrangian
contribution of the unit is \f$ \sigma_i \f$ times that of one copy, and
therefore
\f[
  f^1_i( \bar z_i ) + \sum_{ r } \bar y_r \, g^1_{i,r}( \bar z_i )
  \in \partial_{ \sigma_i } Q \tag{8}
\f]
where \f$ r \f$ runs over the rows (1)-(5) and (1'), \f$ g^1_{i,r} \f$ is the
left-hand side contribution of one copy to row \f$ r \f$ (e.g., \f$ \sum_{ g }
( p^{ac}_{t,g} - P^{au}_{t,g} ( 1 - u_{t,g} ) ) \f$ for the row (1) of node
\f$ n \f$ and instant \f$ t \f$, the sum running over the generators of
\f$ i \f$ at \f$ n \f$), and \f$ \bar z_i \f$ is a minimizer of the Lagrangian
subproblem of one copy at \f$ \bar y \f$. The value \f$ f^1_i( \bar z_i ) +
\bar y^\top g^1_i( \bar z_i ) \f$ of this minimizer is the slope of the
function of \f$ \sigma_i \f$ that is active in the dual function. Note that
the solution of the unit given by the Solver of a convex problem is such a
minimizer (see UnitBlock::scale() for the argument). A capacity multiplier,
instead, multiplies the data \f$ b^1_r \f$ in the finite side
\f$ \kappa b^1_r \f$ of a set \f$ \mathcal{R}_\kappa \f$ of rows of the
sub-Block (for a two-sided row, the side that the sign of its dual value
indicates as active), and therefore
\f[
  - \sum_{ r \in \mathcal{R}_\kappa } \bar y_r \, b^1_r
  \in \partial_{ \kappa } Q \; , \tag{9}
\f]
which for an intermittent unit includes, besides the power bounds \f$ \kappa
P^{mn}_t \leq p^{ac}_t \leq \kappa P^{mx}_t \f$, the rows of its reserves,
whose right-hand sides also contain \f$ \kappa P^{mx}_t \f$ or \f$ \kappa
P^{mn}_t \f$ (see IntermittentUnitBlock::get_kappa_linearization()), while for
a battery with binary variables the rows (8) of BatteryUnitBlock add the terms
in which \f$ \kappa \f$ multiplies \f$ u^{ch}_t \f$ (see
BatteryUnitBlock::get_kappa_linearization()). For a line \f$ l \f$
\f[
  - \sum_{ t \in \mathcal{T} } \bar\mu_{t,l} \, C^{v} P^{\star}_{t,l}
  \in \partial_{ \kappa_l } Q \tag{10}
\f]
where \f$ \bar\mu_{t,l} \f$ is the dual value of the flow bounds of line
\f$ l \f$ in the NetworkBlock of instant \f$ t \f$ and \f$ P^{\star}_{t,l} \f$
is the maximum flow of the line at that instant if the upper bound is the
active one (\f$ \bar\mu_{t,l} > 0 \f$), its minimum flow otherwise (see
DCNetworkBlock::set_kappa()); a switchable line of an OTSNetworkBlock adds
the dual value of the limit \f$ \min\{ 1 , \kappa_l \} \f$ of its switching
while \f$ \kappa_l < 1 \f$ (see OTSNetworkBlock::set_kappa()). In (9) and
(10) the minus sign makes the
subgradient nonpositive when more capacity relaxes an active upper bound, as
it must for a minimization. A capacity multiplier of a unit or of a line that
is also in design mode is not covered by these formulae, since there
\f$ \kappa \f$ multiplies the design variable and the optimal value need not
be convex in it. Hence, IntermittentUnitBlock::get_kappa_linearization(),
BatteryUnitBlock::get_kappa_linearization() and, for a line with a design
variable, DCNetworkBlock::get_resize_linearization() throw in that case. A
worst case of \f$ Q \f$ over a
set \f$ \Upsilon \f$ of alternative futures, i.e., the problem
\f[
  \min_x \Bigl\{ \sum_{ i = 1 }^m \bigl( c_i ( x_i - \bar x_i )^+
    + d_i ( \bar x_i - x_i )^+ \bigr) + \max_{ \upsilon \in \Upsilon }
    Q_\upsilon( x ) \Bigr\} \; ,
\f]
is not modeled by InvestmentFunction, whose value is a single (possibly
expected) operational cost. However, one can build it from one
InvestmentFunction \f$ \mathcal{F}_\upsilon \f$ per future as the maximum of
them, and a convex combination of the linearizations of those attaining the
maximum is then a subgradient of it (see InvestmentFunction).

\section ucbm_dual Lagrangian dual

\subsection ucbm_dual_fn The dual function

The rows (1)-(5) and (1'), which link the sub-Blocks with each other, are the
only constraints of the UCBlock; since the UCBlock also has no Variable and
only a constant Objective, it has the structure that LagrangianDualSolver
requires. This Solver relaxes all of these rows at once, and it builds one
Lagrangian subproblem for each UnitBlock and one for each NetworkBlock. Each
row \f$ r \f$ with left-hand side \f$ g_r( x ) \f$, written as in (1)-(5) and
(1') with all the variables on the left, and finite side \f$ b_r \f$ is
relaxed by adding \f$ y_r ( g_r( x ) - b_r ) \f$ to the objective. Let
\f$ y^{ac}_{t,n} \f$ (\f$ y^{ac}_t \f$ in the bus case), \f$ y^{rc}_{t,n} \f$
(\f$ y^{rc}_t \f$), \f$ y^{pr}_{t,z} \f$, \f$ y^{sc}_{t,z} \f$,
\f$ y^{in}_{t,z} \f$ and \f$ y^{p}_{z} \f$ be the multipliers of (1), (1'),
(2), (3), (4) and (5), and \f$ \bar O_{z,p} \f$ the finite side of (5) that is
relaxed; the dual function is
\f[
  \Theta( y ) = \sum_{ i \in \mathcal{I} } \Theta_i( y )
    + \sum_{ k } \Theta^{net}_k( y ) + c^{U}
    - \sum_{ t \in \mathcal{T} } \Bigl(
      \sum_{ z \in \mathcal{Z}^{pr} } y^{pr}_{t,z} D^{pr}_{t,z}
    + \sum_{ z \in \mathcal{Z}^{sc} } y^{sc}_{t,z} D^{sc}_{t,z}
    + \sum_{ z \in \mathcal{Z}^{in} } y^{in}_{t,z} D^{in}_{t,z} \Bigr)
    - \sum_{ p \in \mathcal{P} } \sum_{ z \in \mathcal{Z}^{p} }
      y^{p}_{z} \bar O_{z,p} \tag{11}
\f]
(with the further terms \f$ - \sum_t ( y^{ac}_t D^{ac}_t + y^{rc}_t
D^{rc}_t ) \f$ in the bus case), where \f$ \Theta_i( y ) \f$ is the minimum
over the feasible set of UnitBlock \f$ i \f$ of its cost \f$ f_i \f$ plus,
for each \f$ t \f$ and each generator \f$ g \f$ of \f$ i \f$,
\f[
  \sigma_i \Bigl( y^{ac}_{t,n(g)} \bigl( p^{ac}_{t,g} - P^{au}_{t,g}
    ( 1 - u_{t,g} ) \bigr) + y^{rc}_{t,n(g)} p^{rc}_{t,g}
    + y^{pr}_{t,z^{pr}(g)} p^{pr}_{t,g}
    + y^{sc}_{t,z^{sc}(g)} p^{sc}_{t,g} + y^{in}_{t,z^{in}(g)} \bigl(
      h^u_{t,g} u_{t,g} + h^p_{t,g} p^{ac}_{t,g} \bigr)
    + \sum_{ p \in \mathcal{P} } y^{p}_{z^{p}(g)} \rho_{t,p,g}
      p^{ac}_{t,g} \Bigr) \tag{12}
\f]
(a node in no zone gives no term, and the reactive term is present only when
the network handles reactive power), plus the terms \f$ \sigma_i
y^{p}_{z^{p}(k)} \rho^{v}_{t,p,k} v^{sl}_{t,k} \f$ of its storages. In turn,
\f$ \Theta^{net}_k( y ) \f$ is the minimum of \f$ - \sum_{ t } \sum_{ n } (
y^{ac}_{t,n} S_{t,n} + y^{rc}_{t,n} R_{t,n} ) \f$, plus the cost of the
network if any, over the feasible set of NetworkBlock \f$ k \f$, where the sum
runs over the instants that \f$ k \f$ covers. The constant \f$ - y^{ac}_{t,n}
\sum_g \sigma_g P^{au}_{t,g} \f$ of (12) is in \f$ \Theta_i \f$ above, while
UCBlock writes the row (1) with \f$ \sum_g \sigma_g P^{au}_{t,g} \f$ in its
right-hand side, and therefore LagrangianDualSolver counts it in the constant
of the dual function; of course, the sum \f$ \Theta \f$ is the same. Since the
scale multiplies the complete contribution of a unit, \f$ \Theta_i \f$ is
\f$ \sigma_i \f$ times the Lagrangian value of one copy, which is how the dual
handles \f$ \sigma_i \f$ identical units without approximation. The dual
problem maximizes \f$ \Theta \f$ with \f$ y^{ac} \f$ and \f$ y^{rc} \f$ free,
\f$ y^{pr} , y^{sc} , y^{in} \leq 0 \f$, and \f$ y^{p}_{z} \geq 0 \f$ when
only the upper budget is finite, \f$ y^{p}_{z} \leq 0 \f$ when only the lower
one is, free when the two coincide. Note that a pollutant row with two
distinct finite bounds cannot be relaxed (see LagrangianDualSolver). The
optimal value of the dual problem is a lower bound on the optimum of (0)-(5).
When the costs are linear, it equals the optimum of the problem in which the
feasible set of each sub-Block is replaced by its convex hull. When they are
not (e.g., the quadratic cost of a ThermalUnitBlock), it equals the optimum of
the problem in which each sub-Block is replaced by the closed convex hull of
the epigraph of its cost over its feasible set, i.e., its cost by the convex
envelope of the cost on the convex hull of the feasible set. If the Lagrangian
subproblems are in turn solved as continuous relaxations, the bound is that of
the continuous relaxation of the complete problem.

\subsection ucbm_dual_net The network subproblem and the flow limits

Since the node injections belong to the NetworkBlock and the demand is a datum
of the latter, \f$ S_{t,n} - D^{ac}_{t,n} \f$ is the net balance of node
\f$ n \f$, the network is a separate subproblem, and the multiplier that the
generators at node \f$ n \f$ see is the nodal one \f$ y^{ac}_{t,n} \f$ (in the
bus case, a single one per instant). \f$ \Theta^{net}_k \f$ is the value of a
linear program that only depends on the set of feasible injections, hence it
does not depend on the formulation chosen by DCNetworkBlock. This is the dual
obtained by relaxing the node injection rows and keeping the network
constraints in the network subproblem; we call it the dual of version 2. In
the dual of version 1, instead, the network constraints are relaxed as well.
Consider a DCNetworkBlock without cost of the flows, written in its PTDF
formulation (see DCNetworkBlock::generate_abstract_constraints()) with the
flows of the lines in \f$ \mathcal{L}^{S} \f$ replaced by their expressions
(6) there. Then the network of instant \f$ t \f$ is reduced to the overall
balance (7), the node balances (8) at the nodes of the set
\f$ \mathcal{N}^{b} \f$ defined there and the bounds of all the flows. Version
1 relaxes (7), (8) and the bounds of the lines in \f$ \mathcal{L}^{S} \f$,
with multipliers \f$ y^{ov}_t \f$ and \f$ y^{b}_{t,n} \f$ free and
\f$ \mu^{mn}_{t,l} , \mu^{mx}_{t,l} \geq 0 \f$ (the latter two enter the
Lagrangian as \f$ \mu^{mx}_{t,l} ( F_l - \kappa_l C^v P^{mx}_l ) +
\mu^{mn}_{t,l} ( \kappa_l C^v P^{mn}_l - F_l ) \f$). In the network subproblem
it keeps only the node injections \f$ S_{t,n} \f$, which are free, and the
flows \f$ F_h \f$ of the lines \f$ h \in \mathcal{L}^{H} \f$ within their
bounds. Let \f$ \hat A \f$ be the incidence matrix of the lines in
\f$ \mathcal{L}^{S} \f$ (\f$ +1 \f$ at the start node, \f$ -1 \f$ at the end
node), \f$ A^{H} \f$ that of the HVDC lines, with \f$ - \eta_{h,j} \f$ at
their end nodes, and, for \f$ l \in \mathcal{L}^{S} \f$,
\f[
  \pi_{t,l} = \mu^{mx}_{t,l} - \mu^{mn}_{t,l}
    + \sum_{ m \in \mathcal{N}^{b} } \hat A_{lm} \, y^{b}_{t,m} \; ;
\f]
the coefficient of \f$ S_{t,n} \f$ in the Lagrangian, which must be zero for
the minimum over the free \f$ S_{t,n} \f$ to be finite, is then \f$ -
y^{ac}_{t,n} + y^{ov}_t - [ n \in \mathcal{N}^{b} ] y^{b}_{t,n} + [ n \notin R
] \sum_{ l } \Psi_{ln} \pi_{t,l} \f$, where \f$ [ \cdot ] \f$ is 1 if the
condition holds and 0 otherwise and \f$ R \f$ is the set of the reference
nodes. Hence, at each multiplier with a finite dual value,
\f[
  y^{ac}_{t,n} = y^{ov}_t - [ n \in \mathcal{N}^{b} ] \, y^{b}_{t,n}
    + [ n \notin R ] \sum_{ l \in \mathcal{L}^{S} } \Psi_{ln}
      \pi_{t,l} \; , \tag{13}
\f]
and the nodal price \f$ - y^{ac}_{t,n} \f$ splits into the price
\f$ - y^{ov}_t \f$ common to all the nodes, the price of the node balance at
the nodes of \f$ \mathcal{N}^{b} \f$, and a congestion component that sums
over the lines the effect of an injection at \f$ n \f$ on the flow of the line
times the marginal value \f$ \pi_{t,l} \f$ of the line. For a network with one
component and no HVDC line \f$ \mathcal{N}^{b} \f$ is empty, \f$ \pi_{t,l} =
\mu^{mx}_{t,l} - \mu^{mn}_{t,l} \f$, and (13) reads \f$ - y^{ac}_{t,n} = -
y^{ov}_t - \sum_{ l } \Psi_{ln} ( \mu^{mx}_{t,l} - \mu^{mn}_{t,l} ) \f$, and
the price at the reference is then \f$ - y^{ov}_t \f$. Substituting (13) into
(12), and the remaining terms of the network into \f$ \Theta^{net}_k \f$,
gives the dual function of version 1: it is (11) with \f$ y^{ac}_{t,n} \f$
given by (13) in each \f$ \Theta_i \f$ (where the coefficient of
\f$ p^{ac}_{t,g} \f$ in (12) becomes \f$ \sigma_g \bigl( y^{ov}_t - [ n(g) \in
\mathcal{N}^{b} ] y^{b}_{t,n(g)} + [ n(g) \notin R ] \sum_l \Psi_{l n(g)}
\pi_{t,l} + \sum_p y^{p}_{z^p(g)} \rho_{t,p,g} \bigr) + \sigma_g
y^{in}_{t,z^{in}(g)} h^p_{t,g} \f$), and with \f$ \Theta^{net}_k \f$ replaced
by the sum over the instants of
\f[
  - \sum_{ n \in \mathcal{N} } y^{ac}_{t,n} D^{ac}_{t,n}
  - \sum_{ l \in \mathcal{L}^{S} } \kappa_l C^v \bigl( \mu^{mx}_{t,l}
    P^{mx}_l - \mu^{mn}_{t,l} P^{mn}_l \bigr)
  + \sum_{ h \in \mathcal{L}^{H} } \kappa_h C^v \min \bigl\{
    \phi_{t,h} P^{mn}_h \, , \, \phi_{t,h} P^{mx}_h \bigr\} \; ,
\f]
where \f$ \phi_{t,h} = y^{ov}_t ( \sum_j \eta_{h,j} - 1 ) + \sum_{ n \in
\mathcal{N}^{b} } A^{H}_{hn} y^{b}_{t,n} + \sum_{ l \in \mathcal{L}^{S} }
\pi_{t,l} \mathrm{DCDF}_{lh} \f$ is the coefficient of \f$ F_h \f$ in the
Lagrangian, \f$ y^{ac}_{t,n} \f$ is again (13), and the bounds are those of
the instant. We remark that the two versions give the same bound. Indeed,
\f$ \Theta^{net}_k \f$ of version 2 is the optimal value of a linear program
whose dual, by linear programming duality, is the maximum of the expression
above over the multipliers \f$ ( y^{ov} , y^{b} , \mu^{mn} , \mu^{mx} ) \f$
that satisfy (13) for the given \f$ y^{ac} \f$ (since the network is bounded,
the dual is feasible whenever the network is). Maximizing over \f$ y^{ac} \f$
as well, the dual of version 2 becomes the maximum of the dual function of
version 1 over all its multipliers, i.e., the dual of version 1. Hence, an
optimal solution of either gives one of the other. In fact, the multipliers of
version 1 are recovered from the dual values of the network subproblem in the
PTDF formulation, since \f$ y^{ov}_t \f$, \f$ y^{b}_{t,n} \f$ and
\f$ \mu_{t,l} = \mu^{mx}_{t,l} - \mu^{mn}_{t,l} \f$ are the dual values of its
rows (7) and (8) and of the bounds of the lines (see
DCNetworkBlock::get_dual_prices()). Conversely, (13) gives the nodal
multipliers of version 2 from those of version 1. With a single node the two
duals coincide.

\subsection ucbm_dual_sign Sign and meaning of the dual values

The multipliers are the dual values of the rows, as written by the Solver (see
UCBlock::get_Solution() and UCBlockSolution). By the convention of
RowConstraint, for this minimization problem the dual value of a row is the
coefficient of its left-hand side in the Lagrangian, i.e., minus the
derivative of the optimal value with respect to the finite side that is
active. We check the following consequences in the unit test of the module, on
instances whose duals are known. First, the opposite \f$ - y^{ac}_{t,n} \f$ is
the marginal cost of a unit increase of the demand of node \f$ n \f$ at
instant \f$ t \f$, i.e., the locational marginal price, and
\f$ y^{ac}_{t,n} \f$ is therefore nonpositive whenever a unit increase of the
demand costs something. Also, the duals \f$ y^{pr}_{t,z} \f$,
\f$ y^{sc}_{t,z} \f$ and \f$ y^{in}_{t,z} \f$ are nonpositive, and their
opposites are the marginal costs of the requirements. Finally,
\f$ y^{p}_{z} \f$ is nonnegative when the budget \f$ O_{z,p} \f$ is active, in
which case it is the saving of a unit more of budget, and nonpositive when the
lower bound \f$ O^{mn}_{z,p} \f$ is active, in which case it is minus the cost
of a unit more of lower bound. Since the costs are per instant (see
\ref ucbm_conv_time), so are these values: with instants of \f$ \Delta t \f$
hours a price per MWh is \f$ - y^{ac}_{t,n} / \Delta t \f$, and the same holds
for the reserves. Instead, \f$ y^{in}_{t,z} \f$ is a value per unit of inertia
at that instant and \f$ y^{p}_{z} \f$ a value per unit of emission over the
horizon, and these two need no conversion. A row whose requirement is zero may
be degenerate, and its dual value is then not unique.

\subsection ucbm_dual_hydro Future cost of water

We do not relax the future cost of the water of a HydroSystemUnitBlock, which
is part of its cost \f$ f_i \f$. Hence, all the reservoirs whose volumes are
arguments of one such function belong to the same subproblem, and no
multiplier is associated with the rows (cuts) of the function. At the final
volumes, the coefficients \f$ \beta_j \f$ of a row active there form a
subgradient of the future cost, i.e., the marginal future cost of a unit more
of water left in each reservoir, which is nonpositive when the water has a
value. Hence, the values of the water (\f$ \omega = - \beta_j \f$ in
HydroSystemUnitBlock) are minus a subgradient of the future cost, which the
Lagrangian subproblem of the HydroSystemUnitBlock uses implicitly.

\subsection ucbm_dual_cfg Configuration

As any Solver, a LagrangianDualSolver is attached to the UCBlock through a
BlockSolverConfig; for instance, the configuration "BSPar.txt" of the tests of
UCBlock attaches a :MILPSolver, a LagrangianDualSolver configured by
"LDCfg-easy.txt" and a primal heuristic. In "LDCfg.txt" the dual is maximized
by BundleSolver ("str_LDSlv_ISName"). Each Lagrangian subproblem is solved by
the Solver of the BlockSolverConfig "str_LagBF_BSCfg", which may be
"LPBSCfg.txt" (a :MILPSolver that solves the continuous relaxation of the
subproblem) or a configuration that maps the classname of each sub-Block to
its Solver. For instance, "InnerBSCfg-DP.txt" gives the dynamic programming
ThermalUnitDPSolver, which solves the integer subproblem, to each
ThermalUnitBlock, and the continuous relaxation to the others. With
"intDoEasy" set, the subproblems that are linear programs (the NetworkBlock
among them, the classnames in "vstrNoEasy" excepted) are put in the master
problem of the bundle method, and therefore they are not represented by
cutting planes; this changes the way the dual is solved, while its value stays
the same. Then, LagrangianDualSolver::get_dual_solution() writes the
multipliers in the rows (1)-(5).

\section ucbm_multi Multistage models

A UCBlock can describe the problem of one stage of a multistage model (see
SDDPBlock), in which a long horizon is split into consecutive periods, one per
stage. The uncertainty (demands, inflows, renewable production) is revealed
stage after stage, and the decisions of a stage depend only on what has been
revealed so far, as the nested form (1) of SDDPBlock expresses. The
operational cost of an investment is then the expected cost of the stages,
i.e., the value of SDDPBlock, and one can estimate it in two ways. First, the
cuts that SDDPSolver computes on convex (relaxed) stage problems give a lower
model of the expected cost of the following stages. Second, SDDPGreedySolver
simulates on the scenarios the policy that these cuts define (possibly on
stage problems with their integer variables), and the mean cost of the
simulations is an estimate from above of the expected cost of that policy (see
InvestmentFunction). With a LagrangianDualSolver in each stage, the value of a
stage is the value of its Lagrangian dual, and the value of a scenario is the
sum of these values over the stages. The state passed from one stage to the
next is the vector of the volumes of the reservoirs, whose initial values
enter only the right-hand sides of the water balances of the first instant;
hence, the final volumes are an affine function of the initial volumes, of the
flows and of the inflows (see HydroUnitBlock). As a function of the final
volumes, the expected cost of the following stages is represented by the
polyhedral function of each HydroSystemUnitBlock,
\f[
  \check\nu( v^{f} ) = \max_{ j \in \mathcal{J} } \bigl\{ \alpha_j +
    \beta_j^\top v^{f} \bigr\} \; , \tag{14}
\f]
where \f$ v^{f} \f$ is the vector of the volumes \f$ v^{hy}_{n,T-1} \f$ of its
reservoirs (the state \f$ z_t \f$ of SDDPBlock). The cuts that build it are
computed from the dual values of the water balances of the first instant and
of the rows of the function (see HydroSystemUnitBlock for the derivation of a
cut). Note that the cuts are valid only if the problem of the stage is convex,
i.e., if it is solved as a continuous relaxation or through its Lagrangian
dual, and if its optimal value is convex in the initial volumes (see
SDDPBlock). When the stage is solved by a LagrangianDualSolver, the water
balances are rows of the subproblem of a HydroSystemUnitBlock, which is not
relaxed, and their dual values are those of the last solution of that
subproblem. Since the value of the Lagrangian dual is
\f$ \max_y \Theta( y ) \f$, by Danskin's theorem its derivative with respect
to the initial volumes, which enter only that subproblem, is the one of
\f$ \Theta_i \f$ at an optimal multiplier \f$ \bar y \f$. Hence, the dual
values of the water balances give a subgradient only if that last solution was
computed at \f$ \bar y \f$, and only if the Solver of the subproblem is a
CDASolver, which provides them (see
LagrangianDualSolver::get_dual_solution()). We do not model the aggregation of
the volumes of several reservoirs into fewer state variables, which would
reduce the dimension of the state; one can still give a cut computed on
aggregated volumes \f$ \mathfrak{A} v^{f} \f$ as the cut \f$ \alpha_j + (
\mathfrak{A}^\top \beta_j )^\top v^{f} \f$ on the individual volumes. Any
other initial condition (initial power and up or down time of a thermal unit,
initial flow of a plant, initial level of a battery) is a datum of the stage
and is not part of the state, which amounts to fixing it to a prescribed
value in the cuts. In a simulation of the policy one can still give each stage
the initial power and the initial up or down time reached by the previous one,
through the callback of SDDPGreedySolver (whose end function restores the
data); ThermalUnitBlock::set_init_updown_time() takes such a value after the
generation only if it leaves unchanged whether the unit is on and the first
instant at which it may switch, and ThermalUnitBlock::set_initial_power()
refuses an initial power below the minimum power of a unit that is on.

\section ucbm_not What is not modeled

We do not represent the following features, either because the data can
express them otherwise or because they would break the decomposition of the
model:
- units producing heat, heat-only units and their emissions;
- a synchronous condenser, i.e., a unit whose power is nonpositive while it
  provides inertia, with the energy it absorbs at the start-up (its inertia
  alone can be described by a ThermalUnitBlock with zero power bounds, see
  ThermalUnitBlock);
- a start-up cost that depends on how long the unit has been off, a cost
  of the active power given as a maximum of affine functions, and a
  quadratic cost coupling different instants;
- spinning reserves offered in fixed amounts by discrete reserve states,
  with a hierarchy between the two reserves, a minimum duration, a change
  of the power forced by the start of a reserve, and specific rules during
  a modulation: the reserves are continuous quantities bounded by the
  constraints of each UnitBlock;
- the rules of modulation, the limits on the number of start-ups per day
  and the deep decreases for thermal units that are not nuclear, which can
  be represented by a NuclearUnitBlock, and the rule that keeps a unit on
  during the stability that follows a start-up or a modulation;
- hydro units with discrete operating points, a power curve that is not
  concave or that changes over time, a dependence of the power or of the
  efficiency of a pump on the head, the rules that tie spillage to full
  turbining, that forbid simultaneous pumping and turbining in a reversible
  plant, and that forbid two changes of the operating point in three
  instants, and a ramp on the sum of the turbined and spilled flows (see
  HydroUnitBlock);
- an inertia proportional to the available (rather than dispatched) power
  of an intermittent unit, and an inertia given by a battery;
- susceptances that change over time within one NetworkBlock (one
  NetworkBlock per instant, each with its own NetworkData, gives them),
  a PTDF matrix given as data, the reduction of a network by clustering its
  nodes, and the offsets of the flows of the lines of an aggregated network
  (see DCNetworkBlock);
- the investment in lines of the distribution grid and in the connections
  between transmission and distribution nodes, and a worst case of the
  operational cost over several alternative futures in the investment
  problem;
- offers (bids) of producers and consumers on a market, their acceptance,
  the energy requirement of a consumer over a period, and the bidding zones
  in which offers are cleared: the demand is a datum of each node and
  instant, and the costs are those of the units;
- the relaxation of the flow limits in the Lagrangian dual, which gives
  the same bound (see \ref ucbm_dual_net), the relaxation of the future
  cost of water, with multipliers on its cuts or on copies of the
  (aggregated) final volumes, a joint future cost of reservoirs
  belonging to different HydroSystemUnitBlock, and the aggregation of
  volumes into fewer state variables of a multistage model;
- a future value of the energy left in a battery or in a load curtailment
  contract, and load curtailment as a seasonal storage;
- the initial conditions of the units other than the volumes of the
  reservoirs as part of the state of a multistage model, or their
  optimization, and the outages of the units among its random data;
- a slack unit of unbounded capacity, or with a nonlinear cost;
- the length of the time step (see \ref ucbm_conv_time).
*/

/*--------------------------------------------------------------------------*/

#endif /* __UCBlockModel */

/*--------------------------------------------------------------------------*/
/*----------------------- End File UCBlockModel.h --------------------------*/
/*--------------------------------------------------------------------------*/
