/*--------------------------------------------------------------------------*/
/*------------------------- File ThermalUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class ThermalUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" thermal unit
 * of a Unit Commitment Problem. A ThermalUnitBlock corresponds to a single
 * electrical generator.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Tiziano Bacci \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli, Tiziano Bacci
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ThermalUnitBlock
 #define __ThermalUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

#include "DQuadFunction.h"

#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a thermal unit
/** The ThermalUnitBlock class derives from UnitBlock and implements a
 * "reasonably standard" thermal unit of a Unit Commitment problem, i.e., one
 * electrical generator (a nuclear, coal, lignite, gas, combined cycle or oil
 * plant, as well as a plant burning biomass, whose pollutant factors in the
 * UCBlock are then zero). Its operation over the instants \f$ \mathcal{T} =
 * \{ 0 , \ldots , T - 1 \} \f$ is described by a binary commitment
 * \f$ u_t \f$, binary start-up and shut-down indicators \f$ v_t \f$ and
 * \f$ w_t \f$, the active power \f$ p^{ac}_t \geq 0 \f$ and, if the UCBlock
 * requires them, the primary and secondary spinning reserves
 * \f$ p^{pr}_t \geq 0 \f$ and \f$ p^{sc}_t \geq 0 \f$. Its rows are (i) the
 * minimum up and down times, with the initial state, (ii) the minimum and
 * maximum power, which depend on the availability and on the limits at the
 * start-up and at the shut-down, (iii) the ramp-up and ramp-down limits, (iv)
 * the bounds on the reserves, and (v) a separable convex quadratic cost. They
 * are given in generate_abstract_constraints() and generate_objective(), in
 * each of the seven formulations that generate_abstract_variables() selects.
 * These are also the rows that the dynamic programming Solvers
 * ThermalUnitDPSolver and ThermalUnitExtDPSolver solve exactly, with the
 * exceptions listed in ThermalUnitDPSolverBase.h (a ReferenceSchedule and
 * the fixings of Variable other than the commitment and the design, which
 * they refuse, and the parts of their recursion that they compute by
 * interpolation). How the unit is composed with the other
 * Blocks of a UCBlock (balance of the injections, reserve, inertia and
 * pollutant requirements, scale and design of the unit) is described in
 * \ref ucblock_model.
 *
 * The rows come from the literature, one family of formulations at a time.
 * The variables of the 3bin formulation and the logical row (1) are those of
 *
 *  L.L. Garver "Power Generation Scheduling by Integer Programming -
 *  Development of Theory" Transactions of the AIEE, Part III 81(3), 730 -
 *  734, 1962
 *
 * the minimum up and down time rows (2), (3) are those of
 *
 *  D. Rajan, S. Takriti "Minimum Up/Down Polytopes of the Unit Commitment
 *  Problem with Start-Up Costs" IBM Research Report RC23628, 2005
 *
 * the ramp rows (11) of the 3bin formulation are those of
 *
 *  J.M. Arroyo, A.J. Conejo "Optimal Response of a Thermal Unit to an
 *  Electricity Spot Market" IEEE Transactions on Power Systems 15(3), 1098 -
 *  1104, 2000
 *
 *  M. Carrion, J.M. Arroyo "A Computationally Efficient Mixed-Integer
 *  Linear Formulation for the Thermal Unit Commitment Problem" IEEE
 *  Transactions on Power Systems 21(3), 1371 - 1378, 2006
 *
 * the maximum power rows (20)-(23) are those of
 *
 *  G. Morales-Espana, J.M. Latorre, A. Ramos "Tight and Compact MILP
 *  Formulation for the Thermal Unit Commitment Problem" IEEE Transactions
 *  on Power Systems 28(4), 4897 - 4908, 2013
 *
 *  C. Gentile, G. Morales-Espana, A. Ramos "A Tight MIP Formulation of the
 *  Unit Commitment Problem with Start-Up and Shut-Down Constraints" EURO
 *  Journal on Computational Optimization 5(1-2), 177 - 201, 2017
 *
 * the ramp rows (12), (13) of the T formulation are those of
 *
 *  P. Damci-Kurt, S. Kucukyavuz, D. Rajan, A. Atamturk "A Polyhedral Study
 *  of Production Ramping" Mathematical Programming 158(1-2), 175 - 205, 2016
 *
 * its maximum power rows (24)-(26) are those of
 *
 *  J. Ostrowski, M.F. Anjos, A. Vannelli "Tight Mixed Integer Linear
 *  Programming Formulations for the Unit Commitment Problem" IEEE
 *  Transactions on Power Systems 27(1), 39 - 46, 2012
 *
 *  K. Pan, Y. Guan "A Polyhedral Study of the Integrated Minimum-Up/-Down
 *  Time and Ramping Polytope" arXiv:1604.02184, 2016
 *
 * and the T formulation as a whole is the one of
 *
 *  B. Knueven, J. Ostrowski, J.-P. Watson "On Mixed-Integer Programming
 *  Formulations for the Unit Commitment Problem" INFORMS Journal on
 *  Computing 32(4), 857 - 876, 2020, doi:10.1287/ijoc.2019.0944
 *
 * whose rows (38), (40), (41) are (24), (25), (26) with data constant in
 * time. The state-space graph of the runs on and off of the unit, on which
 * the dynamic programming Solvers work, is that of
 *
 *  A. Frangioni, C. Gentile "Solving Nonlinear Single-Unit Commitment
 *  Problems with Ramping Constraints" Operations Research 54(4), 767 - 775,
 *  2006
 *
 * and the pt, DP and SU formulations, which describe a schedule as a path
 * in it (with the rows (4)-(10), (14)-(16) and (27)), are those of
 *
 *  T. Bacci, A. Frangioni, C. Gentile, K. Tavlaridis-Gyparakis "New
 *  Mixed-Integer Nonlinear Programming Formulations for the Unit Commitment
 *  Problems with Ramping Constraints" Operations Research 72(5), 2153 -
 *  2167, 2024, doi:10.1287/opre.2023.2435
 *
 * while the SD and SUSD formulations, the symmetric counterpart of the SU
 * one and the combination of the two, are those of
 *
 *  T. Bacci, A. Frangioni, C. Gentile "Start-Up/Shut-Down MINLP
 *  Formulations for the Unit Commitment with Ramp Constraints" Technical
 *  Report R. 20-01, IASI-CNR, Rome, 2020
 *
 * the SUSD one with, in addition, the ramp rows over several steps, which
 * extend the single-step ones of the report. In the pt, SU, SD and SUSD
 * formulations the maximum power rows (27) cap the first and the last
 * instant of a run by the ramps as well as by the start-up and shut-down
 * limits, and are thus stronger than those of the two papers [see
 * generate_abstract_constraints()]. Finally, the perspective cuts (33) are
 * those of
 *
 *  A. Frangioni, C. Gentile "Perspective Cuts for a Class of Convex 0-1
 *  Mixed Integer Programs" Mathematical Programming 106(2), 225 - 236, 2006
 *
 *  A. Frangioni, C. Gentile, F. Lacalandra "Tighter Approximated MILP
 *  Formulations for Unit Commitment Problems" IEEE Transactions on Power
 *  Systems 24(1), 105 - 113, 2009
 *
 * All the data of the unit are given per instant, and no length of the time
 * step appears in any row. Thus, a ramp is the largest change of the active
 * power between two consecutive instants (a gradient in MW/h times the length
 * of the step, in hours) and a minimum up or down time is a number of
 * instants (a duration in hours is rounded up to an integer number of
 * instants). In particular, a unit started at \f$ s \f$ satisfies a minimum
 * up time \f$ \tau^+ \f$ if it is on at \f$ s , \ldots , s + \tau^+ - 1 \f$.
 * Similarly, every cost is the cost of one instant: a cost per MWh, per hour
 * or per MW\f$ {}^2 \f$h is multiplied by the length of the step beforehand,
 * and the quadratic cost is \f$ a_t ( p^{ac}_t )^2 \f$, with no factor
 * \f$ 1/2 \f$. With a step of \f$ \Delta t \f$ hours, a linear cost
 * \f$ C_t \f$ per MWh, a fixed cost \f$ C^{fx}_t \f$ per hour on, a quadratic
 * cost \f$ \frac{1}{2} Q_t ( p^{ac}_t )^2 \f$ per hour and a start-up cost
 * \f$ C^{st}_t \f$ per start-up are therefore given as
 * \f[
 *   b_t = C_t \Delta t \; , \qquad c_t = C^{fx}_t \Delta t \; , \qquad
 *   a_t = \tfrac{1}{2} Q_t \Delta t \; , \qquad c^{su}_t = C^{st}_t \; ,
 * \f]
 * i.e., LinearTerm, ConstTerm, QuadTerm and StartUpCost (see deserialize()).
 * The fixed consumption \f$ P^{au}_t \f$ of an off unit enters only the
 * balance of the UCBlock, where the unit injects
 * \f$ p^{ac}_t - P^{au}_t ( 1 - u_t ) \f$, and no cost of the unit is charged
 * on it. In the same way, the inertia \f$ h^u_t u_t \f$ of the unit enters
 * only the inertia requirement of the UCBlock (see "InertiaCommitment" in
 * deserialize()).
 *
 * Several features of a thermal unit are not represented, and we list them
 * here once. (i) Only the separable convex quadratic function of
 * generate_objective() is available as the cost of the power: a cost given as
 * the maximum of affine functions, and a quadratic form coupling different
 * instants, are not supported. (ii) The start-up cost depends on the instant
 * of the start-up but not on how long the unit has been off before it (e.g.,
 * logarithmically). (iii) The spinning reserves are continuous quantities
 * bounded by a fraction of the active power and by the ramps (cf.
 * generate_abstract_constraints()). Hence, reserves offered in fixed amounts
 * in discrete reserve states, a secondary reserve that requires the primary
 * one, a minimum duration of a reserve once it is offered, a change of the
 * output forced by the start of a reserve, and different reserve amounts in
 * different power bands are not represented. (iv) A limit on the number of
 * start-ups, the modulation and stability rules and the deep decreases are
 * not rows of a ThermalUnitBlock; a thermal unit subject to them (nuclear or
 * not) is a NuclearUnitBlock (with "DayLength" equal to the horizon for a
 * limit over the whole horizon). (v) A unit whose power is nonpositive, such
 * as a synchronous condenser that only consumes its rotating losses (and the
 * energy of a start-up), cannot be described, since MinPower, MaxPower,
 * InitialPower and \f$ p^{ac}_t \f$ are nonnegative. Yet, its inertia (but
 * not its consumption) can be described by a unit with MinPower = MaxPower =
 * 0, InertiaCommitment set to its contribution, ConstTerm set to its rotating
 * cost and StartUpCost to the cost of a start-up. Of course, this omits the
 * energy that the unit draws from the grid, both at the instants in which it
 * is on and at a start-up. (vi) A fictitious unit covering any imbalance at a
 * high cost is a SlackUnitBlock rather than a ThermalUnitBlock. */

class ThermalUnitBlock : public UnitBlock
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*----------------------------- CONSTANTS ----------------------------------*/

 /// mask for the first three bits of AR, i.e., the formulation code
 static constexpr unsigned char FormMsk = 7;

 /// mask for the 4th bit of AR, == 1 if the perspective cuts are used
 static constexpr unsigned char PCuts = 8;

 /// mask for the 5th bit of AR, == 1 if v_t and w_t are continuous
 static constexpr unsigned char ZWCont = 16;

 /// the "three binaries" (3bin) formulation is used
 static constexpr unsigned char tbinForm = 0;

 /// the T formulation is used
 static constexpr unsigned char TForm = 1;

 /// the p_t formulation is used
 static constexpr unsigned char ptForm = 2;

 /// the "dynamic programming" (DP) formulation is used
 static constexpr unsigned char DPForm = 3;

 /// the "start-up" (SU) formulation is used
 static constexpr unsigned char SUForm = 4;

 /// the "shut-down" (SD) formulation is used
 static constexpr unsigned char SDForm = 5;

 /// the "start-up shut-down" (SUSD) formulation is used
 static constexpr unsigned char SUSDForm = 6;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father block
 /** Constructor of ThermalUnitBlock, taking possibly a pointer of its father
  * Block. */

 explicit ThermalUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of ThermalUnitBlock

 virtual ~ThermalUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the ThermalUnitBlock. Besides the mandatory "type" attribute of any
  * :Block, the group must contain all the data required by the base
  * UnitBlock, as described in the comments to UnitBlock::deserialize(
  * netCDF::NcGroup ), to which we refer for the dimensions "TimeHorizon"
  * (\f$ T \f$), "NumberIntervals" and "ChangeIntervals".
  *
  * Every time-dependent datum below is a netCDF::NcDouble variable that is
  * either of size 1, and then it has the same value at every instant, or
  * indexed over "NumberIntervals" (or over "TimeHorizon" if
  * "NumberIntervals" is not given), and then its entry \f$ i \f$ is the
  * value at every instant \f$ t \f$ with ChangeIntervals[ \f$ i - 1 \f$ ]
  * \f$ < t \leq \f$ ChangeIntervals[ \f$ i \f$ ], with ChangeIntervals[
  * -1 ] = -1 and the last entry taken as \f$ T - 1 \f$ ("ChangeIntervals"
  * is not read if NumberIntervals \f$ \leq 1 \f$
  * or NumberIntervals \f$ \geq T \f$). Each such datum is expanded to one
  * value per instant when it is read, and we denote by the index \f$ t \f$
  * the value at the instant \f$ t \in \mathcal{T} \f$. All the data are
  * per instant (see the class documentation). The group contains:
  *
  * - The variable "MaxPower", the nominal maximum active power
  *   \f$ \hat P^{mx}_t \f$ of the unit. It is mandatory.
  *
  * - The variable "MinPower", the nominal minimum active power
  *   \f$ \hat P^{mn}_t \f$ of the unit when it is on, with
  *   \f$ 0 \leq \hat P^{mn}_t \leq \hat P^{mx}_t \f$. It is optional with
  *   default 0.
  *
  * - The variable "Availability", the availability \f$ \chi_t \in [ 0 , 1 ]
  *   \f$ of the unit (0 meaning, e.g., an outage or a maintenance). It is
  *   optional with default 1. The operational bounds of the active
  *   power, which are the ones all the rows use, are
  *   \f[
  *     P^{mx}_t = \chi_t \hat P^{mx}_t \; , \qquad
  *     P^{mn}_t = \hat P^{mn}_t \;\text{ if } \chi_t > 0 \; , \quad
  *     P^{mn}_t = 0 \;\text{ if } \chi_t = 0 \; ;
  *   \f]
  *   in particular, at an instant with \f$ \chi_t = 0 \f$ the unit may be
  *   on (if it cannot shut down, e.g.) with \f$ p^{ac}_t = 0 \f$, and it
  *   then pays \f$ c_t \f$ and gives its inertia \f$ h^u_t \f$.
  *
  * - The variables "DeltaRampUp" and "DeltaRampDown", the ramp-up limit
  *   \f$ \Delta^+_t \geq 0 \f$ and the ramp-down limit
  *   \f$ \Delta^-_t \geq 0 \f$, i.e., the largest increase and decrease of
  *   the active power from \f$ t - 1 \f$ to \f$ t \f$ (from InitialPower
  *   if \f$ t = 0 \f$); a ramp 0 keeps the power of an on unit constant.
  *   Each is optional: if it is not given there are no ramp rows of that
  *   kind (get_delta_ramp_up() and get_delta_ramp_down() then return an
  *   empty vector, and their versions for one instant \f$ t \f$ return
  *   \f$ \hat P^{mx}_t \f$; the dynamic programming Solvers take a
  *   ramp no move exceeds, see ThermalUnitDPSolverBase).
  *
  * - The variable "StartUpLimit", the start-up limit \f$ P^{su}_t \f$,
  *   i.e., the largest active power at the instant \f$ t \f$ when the unit
  *   starts up at \f$ t \f$, and the variable "ShutDownLimit", the
  *   shut-down limit \f$ P^{sd}_t \f$, i.e., the largest active power at
  *   the instant \f$ t - 1 \f$ when the unit shuts down at \f$ t \f$ (is on
  *   at \f$ t - 1 \f$ and off at \f$ t \f$). Both are optional with default
  *   \f$ P^{mn}_t \f$ (the unit starts up and shuts down at its
  *   minimum power), which follows the availability when it changes [see
  *   set_availability()], and both must lie in
  *   \f$ [ P^{mn}_t , P^{mx}_t ] \f$.
  *
  * - The variables "PrimaryRho" and "SecondaryRho", the largest fractions
  *   \f$ \rho^{pr}_t \f$ and \f$ \rho^{sc}_t \f$ of the active power that
  *   the unit can offer as primary and secondary spinning reserve. Each is
  *   optional; if it is not given, or it is zero at every instant, the unit
  *   offers no reserve of that kind and has no variable for it. Both are
  *   ignored if ignore_reserve() has been called.
  *
  * - The variables "PrimarySpinningReserveCost" and
  *   "SecondarySpinningReserveCost", the costs \f$ c^{pr}_t \f$ and
  *   \f$ c^{sc}_t \f$ of one unit of reserve (typically the opposite of a
  *   multiplier of the reserve requirement of the UCBlock, hence possibly
  *   negative; no sign is checked). Each is optional with default 0, the
  *   fractions \f$ \rho^{pr}_t \f$ and \f$ \rho^{sc}_t \f$ being no
  *   price; whether a reserve whose cost is 0 is in the Objective is
  *   decided by the Configuration of generate_objective().
  *
  * - The variables "QuadTerm", "LinearTerm" and "ConstTerm", the
  *   coefficients \f$ a_t \geq 0 \f$, \f$ b_t \f$ and \f$ c_t \f$ of the
  *   cost \f$ a_t ( p^{ac}_t )^2 + b_t p^{ac}_t + c_t u_t \f$ of the unit at
  *   the instant \f$ t \f$; each is optional with default 0.
  *
  * - The variable "StartUpCost", the cost \f$ c^{su}_t \f$ of a start-up
  *   at the instant \f$ t \f$, and the variable "ShutDownCost", the cost
  *   \f$ c^{sd}_t \f$ of a shut-down at \f$ t \f$; both are optional with
  *   default 0, and if "ShutDownCost" is absent or zero at every
  *   instant the Objective has no shut-down term.
  *
  * - The variable "ReactiveLinearTerm", the coefficient \f$ b^q_t \f$ of
  *   the reactive power in the Objective (see set_reactive_linear_term());
  *   optional with default 0.
  *
  * - The scalar variable "InitialPower", the active power \f$ p_{-1} \f$
  *   of the unit at the instant \f$ -1 \f$, before the horizon; optional with
  *   default 0, and nonnegative. It is used only if the unit is
  *   on before the horizon (InitUpDownTime \f$ > 0 \f$), and then it is
  *   meant to lie in \f$ [ \hat P^{mn}_0 , \hat P^{mx}_0 ] \f$: a smaller
  *   value is raised to \f$ \hat P^{mn}_0 \f$ with a warning (also when
  *   \f$ \chi_0 = 0 \f$), and a larger one is kept, with a warning. With the
  *   ramps, generate_abstract_constraints() also requires
  *   \f$ p_{-1} + \Delta^+_0 \geq P^{mn}_0 \f$ and
  *   \f$ p_{-1} - \Delta^-_0 \leq P^{mx}_0 \f$.
  *
  * - The scalar variable "InitUpDownTime", of type netCDF::NcInt, the
  *   number \f$ \tau_0 \f$ of instants the unit has been on
  *   (\f$ \tau_0 > 0 \f$) or off (\f$ \tau_0 \leq 0 \f$, for
  *   \f$ - \tau_0 \f$ instants) before the instant 0; \f$ \tau_0 = 0 \f$
  *   means that the unit has shut down at the beginning of the instant 0.
  *   It is optional with default \f$ - \tau^- \f$ if InitialPower is
  *   0 and \f$ \tau^+ \f$ otherwise. A unit with a nonzero InvestmentCost
  *   must have \f$ \tau_0 < 0 \f$.
  *
  * - The scalar variables "MinUpTime" and "MinDownTime", of type
  *   netCDF::NcUint, the minimum up time \f$ \tau^+ \f$ and the minimum
  *   down time \f$ \tau^- \f$, i.e., the least number of consecutive
  *   instants the unit stays on after a start-up and off after a shut-down.
  *   Both are optional with default 1, and a value 0 is taken as 1
  *   (starting up and shutting down at the same instant never pays). The
  *   minimum up time is taken at most \f$ T + \max\{ 1 , \tau_0 \} \f$ if
  *   the unit is on before the horizon and at most \f$ T + 1 \f$ otherwise,
  *   the minimum down time at most \f$ T + \max\{ 1 , - \tau_0 \} \f$ if
  *   the unit is off before the horizon and at most \f$ T + 1 \f$
  *   otherwise: the bound already means that the unit never switches within
  *   the horizon, and a larger value means the same.
  *
  * - The variable "FixedConsumption", the power \f$ P^{au}_t \geq 0 \f$
  *   that the unit draws when it is off at the instant \f$ t \f$: the
  *   UCBlock counts the injection \f$ p^{ac}_t - P^{au}_t ( 1 - u_t ) \f$
  *   of the unit, and no cost of the unit is charged on it. Optional with
  *   default 0; a negative value is rejected.
  *
  * - The variable "InertiaCommitment", the inertia \f$ h^u_t \f$ that the
  *   unit gives when it is on at the instant \f$ t \f$, i.e., the
  *   coefficient of \f$ u_t \f$ in the inertia requirement of the UCBlock;
  *   optional with default 0. For a synchronous machine with inertia
  *   constant \f$ H \f$ (in s) and rated power \f$ \hat P^{mx}_t \f$ a
  *   usual value is \f$ h^u_t = 1.2 \, H \hat P^{mx}_t \f$, the factor
  *   1.2 converting the active into the apparent power under a constant
  *   phase angle (see \ref ucblock_model); the contribution does not
  *   depend on the availability of the unit nor on its active power.
  *
  * - The variables "MaxReactivePower", "MinReactivePower",
  *   "MaxReactivePowerOn" and "MinReactivePowerOn", used only when the
  *   UCBlock has the reactive power (an AC network): the reactive power
  *   \f$ q_t \f$ of the unit, negative when it is absorbed, satisfies
  *   \f$ Q^{mn}_t + Q^{mn,on}_t u_t \leq q_t \leq Q^{mx}_t +
  *   Q^{mx,on}_t u_t \f$ (see generate_abstract_constraints()). All are
  *   optional with default 0 (a vector that is zero at every instant
  *   is the same as an absent one). The variable "InitialReactivePower" is
  *   accepted and not read.
  *
  * - The variable "ReferenceSchedule", an active power \f$ \hat{p}_t \f$
  *   the unit is asked to follow; if it is given, the Objective has the
  *   term \f$ \sigma \sum_{ t \in \mathcal{T} } | p^{ac}_t - \hat{p}_t | \f$
  *   (see generate_objective()). Optional.
  *
  * - The scalar variable "FixToMaximum", of type netCDF::NcInt: if it is
  *   positive, the unit produces its operational maximum power at every
  *   instant, \f$ p^{ac}_t \geq P^{mx}_t \f$ (a negative value is taken as
  *   0). Optional with default 0.
  *
  * - The variables "MaxRampUpSteps" and "MaxRampDownSteps", the numbers
  *   \f$ J^\pm_t \f$ of ramp steps that the SUSD formulation bounds from
  *   each instant (see generate_abstract_constraints()), are computed by
  *   the ThermalUnitBlock (see compute_ramp_steps()) and cannot be given:
  *   the formulation reads \f$ T + 1 \f$ entries of them, one more than a
  *   datum indexed over the instants holds, hence a file that contains
  *   either of them is refused with an exception.
  *
  * - The scalar variable "InvestmentCost", the cost \f$ c^{inv} \f$ of
  *   building the unit: if it is nonzero, the unit has a binary design
  *   variable \f$ x \f$, it can be on only if it is built
  *   (\f$ u_t \leq x \f$), and the Objective has the term
  *   \f$ \sigma c^{inv} x \f$. Optional with default 0. A continuous
  *   design, "MaxCapacityDesign" or "MinCapacityDesign", is rejected with an
  *   exception, since it would multiply the binary commitment in the
  *   maximum power rows: a unit that may be built is described by
  *   "InvestmentCost", a fleet of identical modules committed together by
  *   "Scale", an integer number of modules committed independently by as
  *   many ThermalUnitBlock, each with "InvestmentCost".
  *
  * - The scalar variable "Capacity", the capacity the user may install,
  *   which is stored and returned by get_capacity() and is not used in any
  *   row. Optional with default 0.
  *
  * - The scalar variable "Scale", the scale factor \f$ \sigma \f$ of the
  *   unit (see UnitBlock::scale()); optional with default 1. A unit
  *   with \f$ \sigma = N \f$ is a fleet of \f$ N \f$ identical modules,
  *   whose data (powers, ramps, limits, costs, InvestmentCost) are those of
  *   one module: every coefficient of the Objective is multiplied by
  *   \f$ \sigma \f$, and the UCBlock multiplies by \f$ \sigma \f$ the terms
  *   of the unit in its linking rows (see \ref ucblock_model). The
  *   commitment, start-up and shut-down variables are those of one module
  *   and are shared by the \f$ N \f$ modules, which are therefore all on or
  *   all off at every instant.
  *
  * A ThermalUnitBlock has no capacity multiplier \f$ \kappa \f$:
  * get_kappa() returns 1, and a variable "Kappa" in the group is not read.
  */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /* extends UnitBlock::expected_dims()
  * not necessary, no new dimensions

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends UnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the ThermalUnitBlock
 /** Generates the abstract Variable of the ThermalUnitBlock and decides
  * which formulation is its abstract representation. The formulation is
  * given by an int \f$ wf \f$, read from \p stvv or, if \p stvv is
  * nullptr, from f_BlockConfig->f_static_variables_Configuration (if
  * f_BlockConfig exists): \f$ wf \f$ is the f_value of that Configuration
  * if it is a SimpleConfiguration< int >, and \f$ wf = 1 \f$ (the T
  * formulation with binary start-up and shut-down variables and no
  * perspective cuts) otherwise, also when a non-null \p stvv is of another
  * type.
  * The int is coded bit-wise: \f$ wf \,\&\, 7 \f$ (FormMsk) chooses the
  * formulation, among
  *
  * - 0 (tbinForm), the "three binaries" (3bin) formulation,
  *
  * - 1 (TForm), the T formulation,
  *
  * - 2 (ptForm), the pt formulation,
  *
  * - 3 (DPForm), the "dynamic programming" (DP) formulation,
  *
  * - 4 (SUForm), the "start-up" (SU) formulation,
  *
  * - 5 (SDForm), the "shut-down" (SD) formulation,
  *
  * - 6 (SUSDForm), the "start-up shut-down" (SUSD) formulation,
  *
  * a value 7 being rejected with an exception; the bit 8 (PCuts) adds the
  * perspective reformulation of the quadratic cost, with its variables;
  * the bit 16 (ZWCont) declares the start-up and shut-down variables
  * continuous in \f$ [ 0 , 1 ] \f$ instead of binary, which does not change
  * the optimum (the rows force them to be integer when the commitment is)
  * but may change the behavior of a Solver. The rows of each formulation
  * are given in generate_abstract_constraints().
  *
  * Let \f$ \tau_0 \f$, \f$ \tau^+ \f$ and \f$ \tau^- \f$ be InitUpDownTime,
  * MinUpTime and MinDownTime (see deserialize()). The first instant
  * \f$ t_0 \f$ (init_t) at which the commitment is free is
  * \f[
  *   t_0 = \min\{ T , \max\{ 0 , \tau^+ - \tau_0 \} \}
  *   \;\text{ if } \tau_0 > 0 \; , \qquad
  *   t_0 = \min\{ T , \max\{ 0 , \tau^- + \tau_0 \} \}
  *   \;\text{ if } \tau_0 \leq 0 \; ,
  * \f]
  * i.e., the number of instants that the minimum up (down) time still
  * imposes to a unit that is on (off) before the horizon; in particular
  * \f$ \tau_0 = 0 \f$, a unit that has shut down at the beginning of the
  * instant 0, gives \f$ t_0 = \min\{ T , \tau^- \} \f$. We write
  * \f$ u_{-1} = 1 \f$ if \f$ \tau_0 > 0 \f$ and \f$ u_{-1} = 0 \f$
  * otherwise. The Variable of every formulation are
  *
  * - the binary commitment \f$ u_t \f$, \f$ t \in \mathcal{T} \f$
  *   ("u_thermal"), 1 if the unit is on at \f$ t \f$;
  *
  * - the active power \f$ p^{ac}_t \geq 0 \f$, \f$ t \in \mathcal{T} \f$
  *   ("p_thermal");
  *
  * - the start-up and shut-down indicators \f$ v_t \f$ ("v_thermal") and
  *   \f$ w_t \f$ ("w_thermal"), for \f$ t = t_0 , \ldots , T - 1 \f$ only
  *   (none if \f$ t_0 = T \f$), \f$ v_t = 1 \f$ if the unit is off at
  *   \f$ t - 1 \f$ and on at \f$ t \f$, \f$ w_t = 1 \f$ if it is on at
  *   \f$ t - 1 \f$ and off at \f$ t \f$; we take \f$ v_t = w_t = 0 \f$
  *   for \f$ t < t_0 \f$ wherever they appear in a row;
  *
  * - the primary reserve \f$ p^{pr}_t \geq 0 \f$ ("pr_thermal") if the
  *   UCBlock has a primary reserve requirement (bit 0 of the reserve_vars
  *   set by the UCBlock) and \f$ \rho^{pr} \neq 0 \f$, and the secondary
  *   reserve \f$ p^{sc}_t \geq 0 \f$ ("sc_thermal") under the same
  *   conditions with bit 1 and \f$ \rho^{sc} \f$; a reserve with no
  *   variable is zero;
  *
  * - the reactive power \f$ q_t \f$ ("q_thermal"), free in sign, its
  *   bounds being rows, if the UCBlock has the reactive power;
  *
  * - the binary design variable \f$ x \f$ ("x_thermal") if InvestmentCost
  *   is nonzero;
  *
  * - the variables \f$ a^{ref}_t \geq 0 \f$ ("v_abs_refschd") of the deviation
  *   \f$ | p^{ac}_t - \hat{p}_t | \f$ if a ReferenceSchedule is given;
  *
  * - with PCuts, the variables \f$ z_t \geq 0 \f$ ("z_thermal") that
  *   replace \f$ ( p^{ac}_t )^2 \f$ in the Objective.
  *
  * The state before the horizon fixes some of them: if \f$ \tau_0 > 0 \f$,
  * \f$ u_t = 1 \f$ for \f$ t < t_0 \f$ and \f$ v_t = 0 \f$ for
  * \f$ t_0 \leq t < \min\{ t_0 + \tau^- , T \} \f$, since a unit that
  * cannot shut down before \f$ t_0 \f$ cannot start up again before
  * \f$ t_0 + \tau^- \f$; if \f$ \tau_0 \leq 0 \f$,
  * \f$ u_t = p^{ac}_t = p^{pr}_t = p^{sc}_t = 0 \f$ for \f$ t < t_0 \f$ and
  * \f$ w_t = 0 \f$ for \f$ t_0 \leq t < \min\{ t_0 + \tau^+ , T \} \f$.
  *
  * The 3bin and T formulations have no other Variable. The pt, DP, SU, SD
  * and SUSD formulations describe a schedule as a path in a state-space
  * graph, whose arcs are the runs of the unit on and off. With the indices
  * of the code, which count the instants from 1, an on-arc \f$ (h,k) \f$
  * is a run of the unit at the instants \f$ t \f$ with
  * \f$ h \leq t + 1 \leq k \f$, i.e., started at the instant \f$ h - 1 \f$
  * (or before the horizon if \f$ h = 0 \f$) and shut down at the instant
  * \f$ k \f$ (beyond the horizon if \f$ k \geq T \f$, where the arc
  * \f$ ( h , T + 1 ) \f$ is a run that the end of the horizon cuts).
  * Symmetrically, an off-arc \f$ (k,h) \f$ is a period off at the instants
  * \f$ k , \ldots , h - 2 \f$ (to the end of the horizon if
  * \f$ h = T + 1 \f$, since the instant 0 if \f$ k = 0 \f$). The arcs are
  * those that respect the minimum up and down times and the state before
  * the horizon:
  *
  * - if \f$ \tau_0 > 0 \f$, the on-arcs \f$ ( 0 , k ) \f$ for
  *   \f$ k = t_0 , \ldots , T + 1 \f$ (the run in progress, with
  *   \f$ ( 0 , 0 ) \f$, a shut-down at 0, whenever \f$ t_0 = 0 \f$,
  *   whatever \f$ p_{-1} \f$, which a row fixes to 0 if
  *   \f$ p_{-1} > P^{sd}_0 \f$; see generate_abstract_constraints()), the
  *   off-arcs \f$ ( k , h ) \f$ for \f$ k = t_0 , \ldots , T \f$ and
  *   \f$ h = k + \tau^- + 1 , \ldots , T \f$, plus \f$ ( k , T + 1 ) \f$,
  *   and the on-arcs \f$ ( h , k ) \f$ for
  *   \f$ h = t_0 + \tau^- + 1 , \ldots , T \f$ and
  *   \f$ k = h + \tau^+ - 1 , \ldots , T \f$, plus \f$ ( h , T + 1 ) \f$;
  *
  * - if \f$ \tau_0 \leq 0 \f$, the off-arcs \f$ ( 0 , h ) \f$ for
  *   \f$ h = t_0 + 1 , \ldots , T + 1 \f$, the on-arcs \f$ ( h , k ) \f$
  *   for \f$ h = t_0 + 1 , \ldots , T \f$ and
  *   \f$ k = h + \tau^+ - 1 , \ldots , T \f$, plus \f$ ( h , T + 1 ) \f$,
  *   and the off-arcs \f$ ( k , h ) \f$ for
  *   \f$ k = t_0 + \tau^+ , \ldots , T \f$ and
  *   \f$ h = k + \tau^- + 1 , \ldots , T \f$, plus \f$ ( k , T + 1 ) \f$.
  *
  * A run that reaches the end of the horizon and lasts at least
  * \f$ \tau^+ \f$ instants is therefore represented both by the arc
  * \f$ ( h , T + 1 ) \f$ and by the arc \f$ ( h , T ) \f$ followed by the
  * off-arc \f$ ( T , T + 1 ) \f$, which describe the same schedule. The
  * Variable of these formulations are, besides those listed above,
  *
  * - the binary arc variables \f$ y^{hk}_+ \f$ ("y_plus_thermal") of the
  *   on-arcs and \f$ y^{kh}_- \f$ ("y_minus_thermal") of the off-arcs, in
  *   all five;
  *
  * - in the DP formulation, the power \f$ p^{hk}_t \geq 0 \f$
  *   ("p_h_k_thermal") of the run \f$ (h,k) \f$ at each of its instants,
  *   and with PCuts the variables \f$ z^{hk}_t \geq 0 \f$
  *   ("z_h_k_thermal");
  *
  * - in the SU and SUSD formulations, the power \f$ p^h_t \geq 0 \f$
  *   ("p_h_thermal") of the runs started at \f$ h - 1 \f$, for every
  *   start \f$ h \f$ of an on-arc and \f$ t \geq h - 1 \f$, and with PCuts
  *   the variables \f$ z^h_t \geq 0 \f$ ("z_h_thermal");
  *
  * - in the SD and SUSD formulations, the power
  *   \f$ \tilde{p}^k_t \geq 0 \f$ ("p_k_thermal") of the runs shut down at
  *   \f$ k \f$, for every end \f$ k \f$ of an on-arc and \f$ t \leq k - 1
  *   \f$, and with PCuts the variables \f$ \tilde{z}^k_t \geq 0
  *   \f$ ("z_k_thermal");
  *
  * - in the SUSD formulation with PCuts, the variables
  *   \f$ \vartheta_t \geq 0 \f$ ("teta_thermal"), \f$ t \in \mathcal{T} \f$,
  *   that replace \f$ ( p^{ac}_t )^2 \f$ in the Objective. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the ThermalUnitBlock
 /** Generates the abstract Constraint of the ThermalUnitBlock in the
  * formulation chosen by generate_abstract_variables(), whose notation we
  * use: \f$ t_0 \f$ is init_t, \f$ u_{-1} \f$ the state before the horizon,
  * \f$ p_{-1} \f$ InitialPower if \f$ \tau_0 > 0 \f$ and 0 otherwise,
  * \f$ P^{mn}_t \f$ and \f$ P^{mx}_t \f$ the operational bounds of the
  * active power, \f$ \Delta^\pm_t \f$ the ramps of the step from
  * \f$ t - 1 \f$ to \f$ t \f$, \f$ P^{su}_t \f$ the largest power at a
  * start-up at \f$ t \f$ and \f$ P^{sd}_t \f$ the largest power at
  * \f$ t - 1 \f$ before a shut-down at \f$ t \f$ (see deserialize()); the
  * variables \f$ v_t \f$ and \f$ w_t \f$ are 0 for \f$ t < t_0 \f$, and a
  * term in \f$ w_{t+1} \f$ is absent for \f$ t = T - 1 \f$. We also write
  * \f[
  *   \bar{P}^{sd}_t = \min\{ P^{sd}_{t+1} , P^{mx}_t \} \;\text{ for }
  *   t < T - 1 \; , \qquad \bar{P}^{sd}_{T-1} = P^{mx}_{T-1} \; ,
  * \f]
  * the largest power at \f$ t \f$ of a unit that is off at \f$ t + 1 \f$,
  * there being no shut-down limit beyond the horizon. Each group of rows is
  * given with the name under which it is added to the Block; a group whose
  * data are absent (e.g., the ramps) is not generated. If \p stcc (or,
  * if \p stcc is nullptr, f_BlockConfig->f_static_constraints_Configuration)
  * is a SimpleConfiguration< int > with a nonzero value, the bounds
  * \f$ 0 \leq u_t \leq 1 \f$ ("Commitment_bound_Thermal"),
  * \f$ 0 \leq v_t \leq 1 \f$ ("StartUp_binary_bound_Thermal") and
  * \f$ 0 \leq w_t \leq 1 \f$ ("ShoutDown_binary_bound_Thermal") are added
  * as explicit ZOConstraint as well.
  *
  * <b>Commitment, 3bin and T formulations.</b> The rows
  * \f{align*}{
  *   & u_t - u_{t-1} = v_t - w_t
  *     && t = t_0 , \ldots , T - 1 \tag{1} \\
  *   & \sum_{ s = \max\{ t_0 , t - \tau^+ + 1 \} }^{ t } v_s \leq u_t
  *     && t = t^+ , \ldots , T - 1 \tag{2} \\
  *   & \sum_{ s = \max\{ t_0 , t - \tau^- + 1 \} }^{ t } w_s \leq 1 - u_t
  *     && t = t^- , \ldots , T - 1 \tag{3}
  * \f}
  * ("StartUp_ShutDown_Variables_Const_Thermal",
  * "StartUp_Commitment_Const_Thermal",
  * "ShutDown_Commitment_Const_Thermal") hold, with \f$ u_{t_0 - 1} =
  * u_{-1} \f$ in (1) at \f$ t = t_0 \f$, where \f$ t^\pm = \min\{ t_0 +
  * \tau^\pm - 1 , T - 1 \} \f$ (there are no rows (1)-(3) if
  * \f$ t_0 = T \f$). Row (2) requires that a unit started in the last
  * \f$ \tau^+ \f$ instants is on: for binary variables it is equivalent to
  * the conditions \f$ u_t \geq v_s \f$ for each \f$ s \f$ of its window,
  * since it implies each of them and, conversely, two start-ups in a window
  * of \f$ \tau^+ \f$ instants require a shut-down between them, and at the
  * instant \f$ r \f$ of that shut-down the condition \f$ u_r \geq v_s \f$
  * of the first start-up \f$ s \f$ fails; row (3) is the same for the
  * shut-downs. Hence a full window (\f$ t \geq t_0 + \tau^+ - 1 \f$) is
  * written wherever it fits in the horizon, and the windows cut by
  * \f$ t_0 \f$ are not written then, since they are implied by the row at
  * \f$ t_0 + \tau^+ - 1 \f$, whose window contains theirs, and by the
  * fixings of generate_abstract_variables(), also in the continuous
  * relaxation: by (1), \f$ u_{t'} = u_t + \sum_{ t < s \leq t' } ( v_s -
  * w_s ) \f$ for \f$ t' = t_0 + \tau^+ - 1 \f$, and hence the full row at
  * \f$ t' \f$ and \f$ w \geq 0 \f$ give \f$ \sum_{ s \leq t } v_s \leq
  * u_t \f$. When the horizon is too short for a
  * full window, \f$ T < t_0 + \tau^+ \f$, the only row (2) is the one at
  * \f$ T - 1 \f$ with the window cut to \f$ [ t_0 , T - 1 ] \f$: a unit
  * that starts up in the horizon then stays on to its end, and it may not
  * shut down after fewer than \f$ \tau^+ \f$ instants as it could with no
  * row at all; the same for (3) when \f$ T < t_0 + \tau^- \f$. Row (1)
  * replaces the usual \f$ u_t - u_{t-1} \leq v_t \f$: the two describe the
  * same commitments, \f$ v_t = w_t = 1 \f$ being excluded by (2) and (3). If
  * \f$ \tau_0 > 0 \f$ and \f$ \tau_0 < \tau^+ \f$, the fixings \f$ u_t = 1
  * \f$, \f$ t < t_0 \f$, are also added as BoxConstraint
  * ("Commitment_fixed_to_one_Thermal").
  *
  * <b>Commitment, pt, DP, SU, SD and SUSD formulations.</b> The schedule is
  * a path from a source to a sink in the state-space graph described in
  * generate_abstract_variables(), whose nodes are the source, the starts
  * \f$ h \f$ of the on-arcs (OFF nodes), the ends \f$ k \f$ of the on-arcs
  * (ON nodes) and the sink \f$ T + 1 \f$; with the sums over the arcs of
  * the graph only, the rows ("Network_Const_Thermal") are
  * \f{align*}{
  *   & \sum_{ k } y^{0k}_+ = 1 \;\text{ if } \tau_0 > 0 \; , \qquad
  *     \sum_{ h } y^{0h}_- = 1 \;\text{ if } \tau_0 \leq 0 \; ,
  *     \tag{4} \\
  *   & \sum_{ k } y^{kh}_- = \sum_{ k } y^{hk}_+
  *     \qquad \text{for every OFF node } h \; , \tag{5} \\
  *   & \sum_{ h } y^{hk}_+ = \sum_{ h } y^{kh}_-
  *     \qquad \text{for every ON node } k \; , \tag{6} \\
  *   & \sum_{ h } y^{h,T+1}_+ + \sum_{ k } y^{k,T+1}_- = 1 \; , \tag{7}
  * \f}
  * and the variables of the 3bin formulation are linked to the arcs by
  * \f{align*}{
  *   & u_t = \sum_{ (h,k) \,:\, h \leq t + 1 \leq k } y^{hk}_+
  *     && t \in \mathcal{T} \; , \tag{8} \\
  *   & v_t = \sum_{ k \geq t + 1 } y^{t+1,k}_+ \; , \qquad
  *     w_t = \sum_{ h \leq t } y^{ht}_+
  *     && t = t_0 , \ldots , T - 1 \tag{9}
  * \f}
  * ("Eq_Commitment_Const_Thermal", "Eq_StartUp_Const_Thermal",
  * "Eq_ShutDown_Const_Thermal"); the arcs enforce the minimum up and down
  * times, and (1)-(3) are not written. The active power is linked to the
  * powers of the runs ("Eq_ActivePower_Const_Thermal") by
  * \f[
  *   p^{ac}_t = \sum_{ (h,k) \,:\, h \leq t + 1 \leq k } p^{hk}_t \;\;
  *   \text{(DP)} \; , \quad
  *   p^{ac}_t = \sum_{ h \leq t + 1 } p^h_t \;\; \text{(SU, SUSD)}
  *   \; , \quad
  *   p^{ac}_t = \sum_{ k \geq t + 1 } \tilde{p}^k_t \;\;
  *   \text{(SD, SUSD)} \tag{10}
  * \f]
  * for \f$ t \in \mathcal{T} \f$, the pt formulation writing its rows on
  * \f$ p^{ac}_t \f$ itself. Below, \f$ Y_t = \sum_{ (h,k) : h \leq t + 1
  * \leq k } y^{hk}_+ \f$ is the commitment at \f$ t \f$ in terms of the
  * arcs, \f$ Y^h_t = \sum_{ k \geq t + 1 } y^{hk}_+ \f$ that of the runs
  * started at \f$ h - 1 \f$ and \f$ \tilde{Y}^k_t = \sum_{ h \leq t + 1 }
  * y^{hk}_+ \f$ that of the runs shut down at \f$ k \f$.
  *
  * <b>Ramps</b> ("RampUp_Const_Thermal", "RampDown_Const_Thermal",
  * generated if DeltaRampUp, respectively DeltaRampDown, is given).
  *
  * - 3bin formulation, for \f$ t \in \mathcal{T} \f$:
  *   \f[
  *     p^{ac}_t - p^{ac}_{t-1} \leq \Delta^+_t u_{t-1} + P^{su}_t v_t
  *     \; , \qquad
  *     p^{ac}_{t-1} - p^{ac}_t \leq \Delta^-_t u_t + P^{sd}_t w_t \; ,
  *     \tag{11}
  *   \f]
  *   where at \f$ t = 0 \f$ the power and the commitment at \f$ t - 1 \f$
  *   are \f$ p_{-1} \f$ and \f$ u_{-1} \f$. For binary commitments and
  *   the default limits \f$ P^{su}_t = P^{sd}_t = P^{mn}_t \f$ these are
  *   \f$ p^{ac}_t \leq p^{ac}_{t-1} + \Delta^+_t u_{t-1} + P^{mn}_t
  *   ( 1 - u_{t-1} ) \f$ and \f$ p^{ac}_{t-1} \leq p^{ac}_t + \Delta^-_t
  *   u_t + P^{mn}_t ( 1 - u_t ) \f$, i.e., a unit starts up and shuts down
  *   at its minimum power.
  *
  * - T formulation, for \f$ t \in \mathcal{T} \f$:
  *   \f{align*}{
  *     p^{ac}_t - p^{ac}_{t-1} & \leq ( \Delta^+_t + P^{mn}_t ) u_t
  *       - P^{mn}_t u_{t-1} + ( P^{su}_t - P^{mn}_t - \Delta^+_t ) v_t
  *       - ( P^{mn}_{t-1} - P^{mn}_t ) w_t \; , \tag{12} \\
  *     p^{ac}_{t-1} - p^{ac}_t & \leq ( \Delta^-_t + P^{mn}_t ) u_{t-1}
  *       - P^{mn}_t u_t + ( P^{sd}_t - P^{mn}_t - \Delta^-_t ) w_t \; ,
  *       \tag{13}
  *   \f}
  *   with \f$ p_{-1} \f$ and \f$ u_{-1} \f$ at \f$ t = 0 \f$; the last term
  *   of (12), at \f$ t \geq \max\{ 1 , t_0 \} \f$, makes the row read
  *   \f$ p^{ac}_{t-1} \geq P^{mn}_{t-1} \f$ at a shut-down, and it is in
  *   the row also when the minimum power does not change from \f$ t - 1
  *   \f$ to \f$ t \f$, with coefficient 0, so that the Variable of the
  *   row do not depend on the data. For binary
  *   commitments (12) is \f$ p^{ac}_t - p^{ac}_{t-1} \leq \Delta^+_t \f$ on
  *   at \f$ t - 1 \f$ and \f$ t \f$, \f$ p^{ac}_t \leq P^{su}_t \f$ at a
  *   start-up, and the same holds for (13).
  *
  * - pt formulation, for \f$ t \geq 1 \f$:
  *   \f{align*}{
  *     p^{ac}_t - p^{ac}_{t-1} & \leq \Delta^+_t \sum_{ h \leq t < k }
  *       y^{hk}_+ - P^{mn}_{t-1} \sum_{ h \leq t } y^{ht}_+
  *       + P^{su}_t \sum_{ k \geq t + 1 } y^{t+1,k}_+ \; , \\
  *     p^{ac}_{t-1} - p^{ac}_t & \leq \Delta^-_t \sum_{ h \leq t < k }
  *       y^{hk}_+ + P^{sd}_t \sum_{ h \leq t } y^{ht}_+
  *       - P^{mn}_t \sum_{ k \geq t + 1 } y^{t+1,k}_+ \; , \tag{14}
  *   \f}
  *   and, if \f$ \tau_0 > 0 \f$, at \f$ t = 0 \f$
  *   \f$ ( p_{-1} - \Delta^-_0 ) \sum_{ k \geq 1 } y^{0k}_+ \leq
  *   p^{ac}_0 \leq ( p_{-1} + \Delta^+_0 ) \sum_{ k \geq 1 } y^{0k}_+ \f$.
  *
  * - DP formulation, for every run \f$ (h,k) \f$ on at \f$ t - 1 \f$ and
  *   \f$ t \geq 1 \f$ (\f$ h \leq t \f$, \f$ k \geq t + 1 \f$):
  *   \f[
  *     - \Delta^-_t y^{hk}_+ \leq p^{hk}_t - p^{hk}_{t-1} \leq
  *     \Delta^+_t y^{hk}_+ \; , \tag{15}
  *   \f]
  *   and, if \f$ \tau_0 > 0 \f$, at \f$ t = 0 \f$ for the runs \f$ (0,k) \f$
  *   \f$ ( p_{-1} - \Delta^-_0 ) y^{0k}_+ \leq p^{0k}_0 \leq
  *   ( p_{-1} + \Delta^+_0 ) y^{0k}_+ \f$.
  *
  * - SU formulation, for every start \f$ h \leq t \f$ and \f$ t \geq 1 \f$:
  *   \f{align*}{
  *     p^h_t - p^h_{t-1} & \leq \Delta^+_t Y^h_t - P^{mn}_{t-1} y^{ht}_+
  *       \; , \\
  *     p^h_{t-1} - p^h_t & \leq \Delta^-_t Y^h_t + P^{sd}_t y^{ht}_+
  *       \; , \tag{16}
  *   \f}
  *   and, if \f$ \tau_0 > 0 \f$, at \f$ t = 0 \f$
  *   \f$ ( p_{-1} - \Delta^-_0 ) Y^0_0 \leq p^0_0 \leq
  *   ( p_{-1} + \Delta^+_0 ) Y^0_0 \f$.
  *
  * - SD formulation, for every end \f$ k \geq t + 1 \f$ and \f$ t \geq 1
  *   \f$:
  *   \f{align*}{
  *     \tilde{p}^k_t - \tilde{p}^k_{t-1} & \leq \Delta^+_t \tilde{Y}^k_{t-1}
  *       + P^{su}_t y^{t+1,k}_+ \; , \\
  *     \tilde{p}^k_{t-1} - \tilde{p}^k_t & \leq \Delta^-_t
  *       \tilde{Y}^k_{t-1} - P^{mn}_t y^{t+1,k}_+ \; , \tag{17}
  *   \f}
  *   and, if \f$ \tau_0 > 0 \f$, at \f$ t = 0 \f$
  *   \f$ ( p_{-1} - \Delta^-_0 ) y^{0k}_+ \leq \tilde{p}^k_0 \leq
  *   ( p_{-1} + \Delta^+_0 ) y^{0k}_+ \f$.
  *
  * - SUSD formulation: the rows (16) downwards on \f$ p^h_t \f$ and the
  *   rows (17) upwards on \f$ \tilde{p}^k_t \f$, with their rows at
  *   \f$ t = 0 \f$ if \f$ \tau_0 > 0 \f$, and, in place of the rows (16)
  *   upwards and (17) downwards, the ramps over several steps: for every
  *   start \f$ h \f$, every \f$ t \geq h - 1 \f$ and
  *   \f$ j = 1 , \ldots , J^+_{t+1} \f$, and for every end \f$ k \f$, every
  *   \f$ t \leq k - 1 - j \f$ and \f$ j = 1 , \ldots , J^-_{t+1} \f$,
  *   \f{align*}{
  *     p^h_{t+j} - p^h_t & \leq \Bigl( \sum_{ r = t + 1 }^{ t + j }
  *       \Delta^+_r \Bigr) \sum_{ k \geq t + j + 1 } y^{hk}_+
  *       - P^{mn}_t \sum_{ k = t + 1 }^{ t + j } y^{hk}_+ \; , \\
  *     \tilde{p}^k_t - \tilde{p}^k_{t+j} & \leq \Bigl( \sum_{ r = t + 1
  *       }^{ t + j } \Delta^-_r \Bigr) \sum_{ h \leq t + 1 } y^{hk}_+
  *       - P^{mn}_{t+j} \sum_{ h = t + 2 }^{ t + j + 1 } y^{hk}_+ \; ,
  *       \tag{18}
  *   \f}
  *   ("RampUp_SUSD_Const_Thermal", "RampDown_SUSD_Const_Thermal", the
  *   groups "RampUp_Const_Thermal" and "RampDown_Const_Thermal" holding
  *   the rows (17) upwards and (16) downwards), where \f$ J^+_{t+1} \f$
  *   (MaxRampUpSteps) is the largest \f$ j \leq T - t - 1 \f$ with
  *   \f$ \sum_{ s = t + 1 }^{ t + j } \Delta^+_s \leq \max_{ s \geq t }
  *   P^{mx}_s - P^{mn}_t \f$, and
  *   \f$ J^-_{t+1} \f$ (MaxRampDownSteps) the largest one with
  *   \f$ \sum_{ s = t + 1 }^{ t + j } \Delta^-_s \leq P^{mx}_t -
  *   \min_{ s \geq t } P^{mn}_s \f$: beyond them the rows are implied by
  *   the bounds, since a run on at \f$ t \f$ and \f$ t + j \f$ has
  *   \f$ p^h_{t+j} - p^h_t \leq P^{mx}_{t+j} - P^{mn}_t \f$ (and a run that
  *   ends in between leaves the row \f$ - p^h_t \leq - P^{mn}_t \f$), also
  *   when the bounds change in time, e.g., at an instant with
  *   \f$ \chi_t = 0 \f$. If \f$ \tau_0 > 0 \f$ there are also the rows
  *   \f$ p^0_j \leq ( p_{-1} + \sum_{ s = 0 }^{ j } \Delta^+_s )
  *   \sum_{ k \geq j + 1 } y^{0k}_+ \f$ for \f$ j = 0 , \ldots , J^+_0 - 1
  *   \f$, where \f$ J^+_0 \f$ is one more than the last \f$ j \f$ at which
  *   the bound is at most \f$ P^{mx}_j \f$ (0 if there is none), and
  *   \f$ ( p_{-1} - \sum_{ s = 0 }^{ j } \Delta^-_s ) y^{0k}_+ \leq
  *   \tilde{p}^k_j \f$ for \f$ j = 0 , \ldots , J^-_0 - 1 \f$, with
  *   \f$ J^-_0 \f$ one more than the last \f$ j \f$ at which the bound is
  *   at least \f$ P^{mn}_j \f$; at the other instants the bound is not
  *   tighter than \f$ P^{mx}_j \f$, respectively \f$ P^{mn}_j \f$. With
  *   \f$ j = 1 \f$ the rows (18) are the rows (16) upwards and (17)
  *   downwards, and with \f$ j = 0 \f$ those at \f$ 0 \f$ are their rows
  *   at \f$ t = 0 \f$: thus, the SUSD formulation contains all the ramp
  *   rows of the SU and SD formulations, as in the technical report cited
  *   in the class description, together with those over more than one
  *   step. Both families matter for the continuous relaxation: the powers
  *   \f$ p^h_t \f$ and \f$ \tilde{p}^k_t \f$ of a fractional point are
  *   linked by the active power alone [see (10)], so that the upward rows
  *   on \f$ p^h_t \f$ do not imply those on \f$ \tilde{p}^k_t \f$, nor
  *   the downward rows on \f$ \tilde{p}^k_t \f$ those on \f$ p^h_t \f$.
  *
  * <b>Shut-down at the instant 0.</b> A unit on before the horizon with
  * \f$ t_0 = 0 \f$ shuts down at 0 with the power \f$ p_{-1} \f$ at
  * \f$ -1 \f$, which, as the power before any shut-down, has to be at most
  * the shut-down limit: the shut-down at 0 is allowed only if
  * \f$ p_{-1} \leq P^{sd}_0 \f$, whatever the ramps, in every formulation
  * and in the dynamic programming Solver. In the 3bin and T formulations
  * the ramp-down row (11), respectively (13), at \f$ t = 0 \f$ says so if
  * DeltaRampDown is given, and otherwise the row
  * ("ShutDownZero_Const_Thermal")
  * \f[
  *   w_0 \leq [ \, p_{-1} \leq P^{sd}_0 \, ]
  * \f]
  * is there, where \f$ [ \cdot ] \f$ is 1 if the condition holds and 0
  * otherwise; in the pt, DP, SU, SD and SUSD formulations the same row is
  * written on \f$ y^{00}_+ \f$, the interval of the shut-down at 0,
  * whether or not the ramps are given. A unit that cannot shut down at 0
  * can at 1 if the ramp-down lets it reach \f$ P^{sd}_1 \f$ at 0, which
  * is always the case without DeltaRampDown.
  *
  * <b>Minimum power</b> ("MinPower_Const_Thermal"): for
  * \f$ t \in \mathcal{T} \f$,
  * \f[
  *   P^{mn}_t u_t \leq p^{ac}_t \;\text{ (3bin, T)} \; , \quad
  *   P^{mn}_t Y_t \leq p^{ac}_t \;\text{ (pt)} \; , \quad
  *   P^{mn}_t y^{hk}_+ \leq p^{hk}_t \;\text{ (DP)} \; , \quad
  *   P^{mn}_t Y^h_t \leq p^h_t \; , \quad
  *   P^{mn}_t \tilde{Y}^k_t \leq \tilde{p}^k_t \; , \tag{19}
  * \f]
  * the last two in the SU and SD formulations respectively, both in SUSD.
  *
  * <b>Maximum power</b> ("MaxPower_Const_Thermal").
  *
  * - 3bin formulation: \f$ p^{ac}_t \leq P^{mx}_t u_t \f$ for
  *   \f$ t < t_0 - 1 \f$, \f$ p^{ac}_{t_0-1} \leq P^{mx}_{t_0-1}
  *   u_{t_0-1} - ( P^{mx}_{t_0-1} - \bar{P}^{sd}_{t_0-1} ) w_{t_0} \f$ if
  *   \f$ 0 < t_0 < T \f$ (the shut-down at \f$ t_0 \f$ ends the run that the
  *   state before the horizon imposes), \f$ p^{ac}_{T-1} \leq P^{mx}_{T-1}
  *   u_{T-1} - ( P^{mx}_{T-1} - P^{su}_{T-1} ) v_{T-1} \f$ if
  *   \f$ T - 1 \geq t_0 \f$, and for \f$ t_0 \leq t < T - 1 \f$
  *   \f[
  *     p^{ac}_t \leq P^{mx}_t u_t - ( P^{mx}_t - \bar{P}^{sd}_t ) w_{t+1}
  *     - ( P^{mx}_t - P^{su}_t ) v_t \tag{20}
  *   \f]
  *   if \f$ \tau^+ \geq 2 \f$, while if \f$ \tau^+ = 1 \f$, when a run of
  *   one instant is possible, the two rows
  *   \f{align*}{
  *     p^{ac}_t & \leq P^{mx}_t u_t - ( P^{mx}_t - \bar{P}^{sd}_t ) w_{t+1}
  *       - ( \bar{P}^{sd}_t - \min\{ P^{su}_t , P^{mx}_t \} )^+ v_t
  *       \; , \\
  *     p^{ac}_t & \leq P^{mx}_t u_t - ( \min\{ P^{su}_t , P^{mx}_t \}
  *       - \bar{P}^{sd}_t )^+ w_{t+1} - ( P^{mx}_t - P^{su}_t ) v_t \tag{21}
  *   \f}
  *   (the same as (23) below) give the cap
  *   \f$ \min\{ P^{su}_t , \bar{P}^{sd}_t \} \f$ to such a run,
  *   \f$ P^{su}_t \f$ to a start-up, \f$ \bar{P}^{sd}_t \f$ to the last
  *   instant before a shut-down and \f$ P^{mx}_t \f$ otherwise.
  *
  * - T formulation: \f$ p^{ac}_t \leq P^{mx}_t u_t \f$ for
  *   \f$ t \in \mathcal{T} \f$, with the term
  *   \f$ - ( P^{mx}_{t_0-1} - \bar{P}^{sd}_{t_0-1} ) w_{t_0} \f$ at
  *   \f$ t = t_0 - 1 \f$ if \f$ 0 < t_0 < T \f$, as in the 3bin formulation,
  *   and for \f$ t = t_0 , \ldots , T - 1 \f$
  *   \f[
  *     p^{ac}_t \leq P^{mx}_t u_t - ( P^{mx}_t - P^{su}_t ) v_t
  *     - ( P^{mx}_t - \bar{P}^{sd}_t ) w_{t+1} \tag{22}
  *   \f]
  *   if \f$ \tau^+ \geq 2 \f$, while if \f$ \tau^+ = 1 \f$ the two rows
  *   \f{align*}{
  *     p^{ac}_t & \leq P^{mx}_t u_t - ( P^{mx}_t - P^{su}_t ) v_t
  *       - ( \min\{ P^{su}_t , P^{mx}_t \} - \bar{P}^{sd}_t )^+ w_{t+1}
  *       \; , \\
  *     p^{ac}_t & \leq P^{mx}_t u_t - ( P^{mx}_t - \bar{P}^{sd}_t ) w_{t+1}
  *       - ( \bar{P}^{sd}_t - \min\{ P^{su}_t , P^{mx}_t \} )^+ v_t
  *       \; , \tag{23}
  *   \f}
  *   the second with no term at all at \f$ t = T - 1 \f$, the terms
  *   \f$ ( \cdot )^+ \f$ being in the rows, with coefficient 0, also when
  *   they vanish. With the ramps,
  *   let \f$ \nu^+_t = \lfloor ( P^{mx}_t - P^{su}_t ) / \Delta^+_t \rfloor
  *   \f$ and \f$ \nu^-_t = \lfloor ( P^{mx}_t - P^{sd}_t ) / \Delta^-_t
  *   \rfloor \f$, the numbers of whole ramps between the limits and the
  *   maximum power (taken as \f$ T + \tau^+ \f$, larger than every bound
  *   they are compared with, if the ramp is 0 or the quotient is larger
  *   than that), \f$ K^d_t = \min\{ \tau^+ - 1 , T - t - 2 , \nu^-_t \} \f$
  *   and \f$ K^u_t = \min\{ \tau^+ - 2 - \max\{ 0 , K^d_t \} , \nu^+_t ,
  *   t - t_0 \} \f$, and let
  *   \f[
  *     \lambda^+_{t,s} = \Bigl( P^{mx}_t - P^{su}_{t-s} - \sum_{ r = t - s
  *       + 1 }^{ t } \Delta^+_r \Bigr)^+ , \qquad
  *     \lambda^-_{t,s} = \Bigl( P^{mx}_t - P^{sd}_{t+1+s} - \sum_{ r = t + 1
  *       }^{ t + s } \Delta^-_r \Bigr)^+ ;
  *   \f]
  *   then, for \f$ t = t_0 , \ldots , T - 1 \f$ and with the terms in
  *   \f$ v_{t-s} \f$ only for \f$ t - s \geq t_0 \f$ and those in
  *   \f$ w_{t+1+s} \f$ only for \f$ t + 1 + s \leq T - 1 \f$,
  *   \f{align*}{
  *     p^{ac}_t & \leq P^{mx}_t u_t - ( P^{mx}_t - \bar{P}^{sd}_t ) w_{t+1}
  *       - \sum_{ s = 0 }^{ \tau^+ - 2 } \lambda^+_{t,s} v_{t-s}
  *       \; , \tag{24} \\
  *     p^{ac}_t & \leq P^{mx}_t u_t - \sum_{ s = 0 }^{ \min\{ \tau^+ - 1 ,
  *       \nu^+_t \} } \lambda^+_{t,s} v_{t-s}
  *       \qquad \text{if } 2 \leq \tau^+ < \nu^+_t + 2 \; , \tag{25} \\
  *     p^{ac}_t & \leq P^{mx}_t u_t - \sum_{ s = 0 }^{ K^d_t }
  *       \lambda^-_{t,s} w_{t+1+s} - \sum_{ s = 0 }^{ K^u_t }
  *       \lambda^+_{t,s} v_{t-s}
  *       \qquad \text{if } K^d_t > 0 \; , \tag{26}
  *   \f}
  *   (24) and (25) if DeltaRampUp is given, (26) if both ramps are. The
  *   rows (24) have the terms in \f$ v_{t-s} \f$ for every
  *   \f$ s \leq \tau^+ - 2 \f$, those that cannot bind with coefficient 0
  *   (where \f$ \lambda^+_{t,s} \f$ vanishes), so that their Variable do
  *   not depend on the data; the
  *   rows (25) and (26), which exist only at the instants where the
  *   condition on \f$ \nu^+_t \f$, respectively \f$ K^d_t \f$, holds, are
  *   dynamic ("MaxPower5_Const_Thermal", "MaxPower6_Const_Thermal"),
  *   since a change of the data adds or removes some of them. A unit
  *   started at \f$ t - s \f$ reaches at most \f$ P^{su}_{t-s} + \sum_{
  *   r = t - s + 1 }^{ t } \Delta^+_r \f$ at \f$ t \f$, and one that
  *   shuts down at \f$ t + 1 + s \f$ at most \f$ P^{sd}_{t+1+s} + \sum_{
  *   r = t + 1 }^{ t + s } \Delta^-_r \f$, which are the bounds that
  *   \f$ \lambda^+_{t,s} \f$ and \f$ \lambda^-_{t,s} \f$ express; with
  *   data constant in time, (24), (25) and (26) are the rows (38), (40) and
  *   (41) of the T formulation of Knueven, Ostrowski and Watson (see the
  *   class description), whose (38) also contains (22). A
  *   row is valid because at most one of its terms in \f$ v \f$ and
  *   \f$ w \f$ can be 1 for an integer commitment: two start-ups (or two
  *   shut-downs) of a row are fewer than \f$ \tau^+ \f$ instants apart, and
  *   a start-up and a shut-down of a row would make a run of fewer than
  *   \f$ \tau^+ \f$ instants (\f$ K^u_t + K^d_t \leq \tau^+ - 2 \f$), both
  *   of which the rows (2) and (3) forbid.
  *
  * - DP formulation: \f$ p^{hk}_t \leq P y^{hk}_+ \f$ for every run
  *   \f$ (h,k) \f$ and instant \f$ t \f$ of it, where
  *   \f$ P = \min\{ P^{su}_t , \bar{P}^{sd}_t \} \f$ for a run of one
  *   instant (\f$ h = k = t + 1 \f$), \f$ P = P^{su}_t \f$ at its first
  *   instant (\f$ h = t + 1 < k \f$), \f$ P = \bar{P}^{sd}_t \f$ at its
  *   last one (\f$ h < t + 1 = k \f$), and \f$ P = P^{mx}_t \f$ otherwise.
  *
  * - pt, SU, SD and SUSD formulations: for every instant \f$ t \f$ of a
  *   run \f$ (h,k) \f$, i.e., \f$ h \leq t + 1 \leq k \f$, the first and
  *   the last one included, let
  *   \f[
  *     \psi^{hk}_t = \min\{ P^{mx}_t , \Xi^h_t , \Theta^k_t \} \; ,
  *     \qquad
  *     \Xi^h_t = P^{su}_{h-1} + \sum_{ r = h }^{ t } \Delta^+_r \; ,
  *     \quad \Theta^k_t = P^{sd}_k + \sum_{ r = t + 1 }^{ k - 1 }
  *     \Delta^-_r \; ,
  *   \f]
  *   and \f$ \Xi^0_t = p_{-1} + \sum_{ r = 0 }^{ t } \Delta^+_r \f$ for the
  *   run from before the horizon, the largest power that the start-up and
  *   shut-down limits and the ramps allow at \f$ t \f$ in that run: the sums
  *   are empty at the first instant, \f$ \Xi^h_{h-1} = P^{su}_{h-1} \f$, and
  *   at the last one, \f$ \Theta^k_{k-1} = P^{sd}_k \f$. Without
  *   DeltaRampUp, \f$ \Xi^h_t \f$ is there only at the first instant and
  *   \f$ \Xi^0_t \f$ is dropped (as it is if \f$ \tau_0 \leq 0 \f$);
  *   without DeltaRampDown, \f$ \Theta^k_t \f$ is there only at the last
  *   instant; and \f$ \Theta^k_t \f$ is dropped if \f$ k \geq T \f$, there
  *   being no shut-down within the horizon (the run from before the horizon
  *   keeps it, as every other run). Then the pt formulation has, for
  *   \f$ t \in \mathcal{T} \f$,
  *   \f[
  *     p^{ac}_t \leq \sum_{ (h,k) \,:\, h \leq t + 1 \leq k }
  *     \psi^{hk}_t y^{hk}_+ \; , \tag{27}
  *   \f]
  *   the SU (and SUSD) formulation has \f$ p^h_t \leq \sum_{ k \geq t + 1 }
  *   \psi^{hk}_t y^{hk}_+ \f$ for every start \f$ h \leq t + 1 \f$, and the
  *   SD (and SUSD) formulation has \f$ \tilde{p}^k_t \leq \sum_{ h \leq t +
  *   1 } \psi^{hk}_t y^{hk}_+ \f$ for every end \f$ k \geq t + 1 \f$. If
  *   \f$ \psi^{hk}_t < P^{mn}_t \f$ the run cannot be on at \f$ t \f$, and
  *   the rows (27) and (19) together exclude it; this happens, e.g., when
  *   the shut-down limit \f$ P^{sd}_k \f$, which deserialize() checks
  *   against the bounds of \f$ k \f$, is below \f$ P^{mn}_{k-1} \f$, as at
  *   an instant \f$ k \f$ at which the unit is unavailable. Also,
  *   \f$ \psi^{hk}_t \f$ is at most \f$ P^{su}_t \f$ at the first instant of
  *   a run, at most \f$ \bar{P}^{sd}_t \f$ at its last one, and at most
  *   \f$ \min\{ P^{su}_t , \bar{P}^{sd}_t \} \f$ for a run of one instant
  *   (\f$ h = k = t + 1 \f$, which exists only if \f$ \tau^+ = 1 \f$). The
  *   rows of the papers cited in the class description, which cap the first
  *   and the last instant of a run by these limits alone, are therefore
  *   implied by (27), which is stronger
  *   at the first instant of a run that shuts down soon after it and at the
  *   last instant of a run that started shortly before; also, the arcs
  *   \f$ ( h , T ) \f$ and \f$ ( h , T + 1 ) \f$, which describe the same
  *   schedule, have the same coefficients in (27).
  *
  * <b>Spinning reserves.</b> If the reserve variables exist, then
  * ("PrimaryRho_Const_Thermal", "SecondaryRho_Const_Thermal")
  * \f[
  *   p^{pr}_t \leq \rho^{pr}_t p^{ac}_t \; , \qquad
  *   p^{sc}_t \leq \rho^{sc}_t p^{ac}_t \qquad t \in \mathcal{T} \; ,
  *   \tag{28}
  * \f]
  * and, with \f$ r_t = p^{pr}_t + p^{sc}_t \f$ (the reserve that has no
  * variable being 0), for every formulation ("Reserve_Const_Thermal")
  * \f{align*}{
  *   & P^{mn}_t u_t \leq p^{ac}_t - r_t
  *     && t \in \mathcal{T} \; , \tag{29} \\
  *   & p^{ac}_t + r_t \leq P^{mx}_t u_t + ( P^{su}_t - P^{mx}_t ) v_t \; ,
  *     \quad p^{ac}_t + r_t \leq P^{mx}_t u_t + ( P^{sd}_{t+1} - P^{mx}_t )
  *     w_{t+1}
  *     && t \in \mathcal{T} \; , \tag{30} \\
  *   & p^{ac}_t + r_t - p^{ac}_{t-1} \leq \Delta^+_t u_{t-1}
  *     + P^{mx}_t ( 1 - u_{t-1} )
  *     && t = 1 , \ldots , T - 1 \; , \tag{31} \\
  *   & p^{ac}_{t-1} - p^{ac}_t + r_t \leq \Delta^-_t u_t
  *     + P^{mx}_{t-1} ( 1 - u_t )
  *     && t = 1 , \ldots , T - 1 \; , \tag{32}
  * \f}
  * (31) if DeltaRampUp is given and (32) if DeltaRampDown is; if the unit
  * is on before the horizon (\f$ \tau_0 > 0 \f$) they are also written at
  * \f$ t = 0 \f$, with \f$ u_{-1} = 1 \f$ and the initial power
  * \f$ p_{-1} \f$ in place of \f$ p^{ac}_{-1} \f$ and of
  * \f$ P^{mx}_{-1} \f$, i.e.,
  * \f[
  *   p^{ac}_0 + r_0 \leq p_{-1} + \Delta^+_0 \; , \qquad
  *   p_{-1} - p^{ac}_0 + r_0 \leq \Delta^-_0 u_0 + p_{-1} ( 1 - u_0 )
  *   \; , \tag{31'}
  * \f]
  * the second of which a shut-down at 0 leaves slack. The rows (29)
  * and (30) keep the reserve within the band of the unit, the cap of (30)
  * being \f$ P^{su}_t \f$ at a start-up and \f$ \min\{ P^{sd}_{t+1} ,
  * P^{mx}_t \} = \bar{P}^{sd}_t \f$ at the last instant before a shut-down
  * (the first row capping at \f$ P^{mx}_t \f$ when \f$ v_t = 0 \f$); the rows
  * (31) and (32) make the reserve deliverable within the ramp left after the
  * move of the output, and bind only when the unit is on at \f$ t - 1 \f$ and
  * \f$ t \f$, so that no reserve is held at an instant where the output moves
  * by the full ramp, the instant 0 included when the output moves from the
  * initial power.
  * With the default start-up and shut-down limits and a constant minimum
  * power, (29) and (30) describe the same set as the plain band
  * \f$ P^{mn}_t u_t \leq p^{ac}_t - r_t \f$, \f$ p^{ac}_t + r_t \leq
  * P^{mx}_t u_t \f$ together with the maximum power rows (20)-(27), which
  * exist also without ramp data and already cap \f$ p^{ac}_t \f$ at
  * \f$ P^{su}_t = P^{mn}_t \f$ at a start-up and at \f$ \bar{P}^{sd}_t =
  * P^{mn}_t \f$ at the last instant on, hence force \f$ p^{ac}_t =
  * P^{mn}_t \f$ and \f$ r_t = 0 \f$ there; (31) and (32) are a restriction
  * that the plain band does not have. These rows use only \f$ u \f$, \f$ v
  * \f$, \f$ w \f$ and \f$ p^{ac} \f$, hence the reserve model is the same in
  * all the formulations, and it is the one the DP Solvers implement.
  *
  * <b>Perspective reformulation</b> ("Init_PC_Const_Thermal",
  * "Max_SUSD_PC_Const_Thermal", "Eq_PC_Const_Thermal"), with PCuts. For
  * \f$ P \in \{ P^{mn}_t , P^{mx}_t \} \f$ the two tangent cuts of the
  * perspective of \f$ ( p^{ac}_t )^2 \f$ at the bounds,
  * \f[
  *   z_t \geq 2 P p^{ac}_t - P^2 u_t \;\text{ (3bin, T)} \; , \quad
  *   z_t \geq 2 P p^{ac}_t - P^2 Y_t \;\text{ (pt)} \; , \quad
  *   z^{hk}_t \geq 2 P p^{hk}_t - P^2 y^{hk}_+ \;\text{ (DP)} \; ,
  *   \tag{33}
  * \f]
  * \f$ z^h_t \geq 2 P p^h_t - P^2 Y^h_t \f$ (SU, SUSD) and
  * \f$ \tilde{z}^k_t \geq 2 P \tilde{p}^k_t - P^2 \tilde{Y}^k_t \f$ (SD,
  * SUSD), hold for every \f$ t \f$ and run; the variables of the runs are
  * linked to \f$ z_t \f$ by \f$ z_t = \sum_{ (h,k) } z^{hk}_t \f$ (DP),
  * \f$ z_t = \sum_{ h } z^h_t \f$ (SU, SUSD) and
  * \f$ z_t = \sum_{ k } \tilde{z}^k_t \f$ (SD, SUSD), and in the SUSD
  * formulation \f$ \vartheta_t \geq \sum_{ h } z^h_t \f$ and
  * \f$ \vartheta_t \geq \sum_{ k } \tilde{z}^k_t \f$. The further cuts are
  * separated by generate_dynamic_constraints().
  *
  * <b>Other rows.</b> If InvestmentCost is nonzero,
  * \f$ u_t \leq x \f$ for \f$ t \in \mathcal{T} \f$
  * ("CommitmentDesign_Const_Thermal"). If FixToMaximum is positive,
  * \f$ p^{ac}_t \geq P^{mx}_t \f$ ("FixedGeneration"). If a
  * ReferenceSchedule \f$ \hat{p} \f$ is given,
  * \f$ p^{ac}_t - a^{ref}_t \leq \hat{p}_t \f$ and
  * \f$ - p^{ac}_t - a^{ref}_t \leq - \hat{p}_t \f$
  * ("Norm1_Reference_Schedule"). If the UCBlock has the reactive power,
  * an absent bound being 0 as for the getters and the DP Solvers,
  * \f$ Q^{mn}_t \leq q_t \leq Q^{mx}_t \f$ is added as
  * a BoxConstraint ("ReactivePowerBound_thermal") if neither
  * MaxReactivePowerOn nor MinReactivePowerOn is given and there is no
  * design variable, and otherwise the rows ("ReactivePowerMax_thermal",
  * "ReactivePowerMin_thermal")
  * \f[
  *   q_t \leq Q^{mx,on}_t u_t + Q^{mx}_t \xi \; , \qquad
  *   q_t \geq Q^{mn,on}_t u_t + Q^{mn}_t \xi \; ,
  * \f]
  * where \f$ \xi = x \f$ if the unit has a design variable (and the bound is
  * finite) and \f$ \xi = 1 \f$ otherwise. In every formulation,
  * \f$ 0 \leq p^{ac}_t \leq P^{mx}_t \f$ is added as a BoxConstraint
  * ("ActivePowerBound_thermal"): the rows that tie the power to the
  * commitment imply it, but a Solver that only reads the bounds of the
  * Variable (e.g., a BoxSolver bounding the Objective) needs it to see the
  * power bounded; is_feasible() checks it with the other rows.
  *
  * <b>Changes after the generation.</b> The rows are written by
  * build_rows(), from the data as they are, both here and when a setter
  * changes MaxPower, Availability, InitialPower or (within the values that
  * set_init_updown_time() accepts) InitUpDownTime after the generation
  * [see update_rows()]: the Variable do not depend on these data (the
  * interval \f$ ( 0 , 0 ) \f$ is there whatever \f$ p_{-1} \f$ and the
  * limits, see generate_abstract_variables()), nor do the Variable of a
  * static row, hence a change enters the coefficients and the sides of the
  * static rows, which get the new values, and the number of the dynamic
  * rows (18), (25) and (26), which are added and removed; the abstract
  * representation is then the one a fresh load of the changed data would
  * generate, apart from the order of the dynamic rows.
  *
  * Before the ramp rows are built, the method throws std::logic_error if the
  * unit is on before the horizon and DeltaRampUp is given with
  * \f$ p_{-1} + \Delta^+_0 < P^{mn}_0 \f$, or DeltaRampDown is given with
  * \f$ p_{-1} - \Delta^-_0 > P^{mx}_0 \f$ (also if a shut-down at 0 would
  * make the schedule feasible). */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the dynamic constraint of the ThermalUnitBlock
 /** With PCuts (see generate_abstract_variables()), separates the
  * perspective cuts of the quadratic cost at the current values of the
  * Variable and adds the violated ones to the dynamic group "PC_cuts_Thermal";
  * otherwise it does nothing. Two parameters drive the separation: a
  * tolerance \f$ \varepsilon \f$ (1e-6 by default) below which the
  * commitment of a run is taken as 0, and a relative threshold
  * \f$ \delta \f$ (1e-7 by default) on the violation. They are read from
  * \p dycc if it is a SimpleConfiguration< double > (\f$ \delta \f$) or a
  * SimpleConfiguration< std::pair< double , double > > (\f$ \delta \f$,
  * \f$ \varepsilon \f$), and otherwise in the same way from
  * f_BlockConfig->f_dynamic_constraints_Configuration.
  *
  * Let \f$ ( \hat{p} , \hat{z} , \hat{y} ) \f$ be the current values of a
  * power, of its perspective variable and of its commitment, i.e.,
  * \f$ ( p^{ac}_t , z_t , u_t ) \f$ in the 3bin and T formulations,
  * \f$ ( p^{ac}_t , z_t , Y_t ) \f$ in the pt one,
  * \f$ ( p^{hk}_t , z^{hk}_t , y^{hk}_+ ) \f$ in the DP one,
  * \f$ ( p^h_t , z^h_t , Y^h_t ) \f$ in the SU and SUSD ones and
  * \f$ ( \tilde{p}^k_t , \tilde{z}^k_t , \tilde{Y}^k_t ) \f$ in the SD and
  * SUSD ones, in the notation of generate_abstract_constraints(). If
  * \f$ \hat{y} > \varepsilon \f$, \f$ \check{p} = \hat{p} / \hat{y} \f$ and
  * \f$ m = \max\{ 1 , \hat{z} + \check{p}^2 \hat{y} \} \f$ satisfy
  * \f$ \hat{p}^2 / \hat{y} - \hat{z} \geq \delta m \f$ (the violation of
  * the perspective constraint \f$ z \, y \geq p^2 \f$ at the point,
  * relative to its size, and absolute when the size is below 1), the
  * cut
  * \f[
  *   z \geq 2 \check{p} \, p - \check{p}^2 \, y \; ,
  * \f]
  * tangent to the perspective of \f$ p^2 \f$ at \f$ \check{p} \f$, is added
  * with the corresponding variables, the commitment \f$ y \f$ being the
  * sum of arc variables that defines \f$ Y_t \f$, \f$ Y^h_t \f$ or
  * \f$ \tilde{Y}^k_t \f$ in the formulations where it is one. */

 void generate_dynamic_constraints( Configuration * dycc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the ThermalUnitBlock
 /** Generates the Objective of the ThermalUnitBlock, a DQuadFunction to be
  * minimized, which in the notation of generate_abstract_variables() is
  * \f{align*}{
  *   \sigma \Bigl( & \sum_{ t = t_0 }^{ T - 1 } ( c^{su}_t v_t
  *     + c^{sd}_t w_t ) + \sum_{ t \in \mathcal{T} } \bigl( a_t
  *     ( p^{ac}_t )^2 + b_t p^{ac}_t + c_t u_t \bigr) \\
  *   & + \sum_{ t \in \mathcal{T} } \bigl( c^{pr}_t p^{pr}_t
  *     + c^{sc}_t p^{sc}_t + a^{ref}_t + b^q_t q_t \bigr)
  *     + c^{inv} x \Bigr) \; , \tag{34}
  * \f}
  * where \f$ \sigma \f$ is the scale factor, \f$ a_t \f$, \f$ b_t \f$,
  * \f$ c_t \f$, \f$ c^{su}_t \f$, \f$ c^{sd}_t \f$, \f$ c^{pr}_t \f$,
  * \f$ c^{sc}_t \f$, \f$ b^q_t \f$ and \f$ c^{inv} \f$ are the data of
  * deserialize(), and every term whose variables do not exist is absent:
  * the shut-down term if ShutDownCost is absent or zero, the reserve terms
  * if the reserve variables do not exist, the term
  * \f$ a^{ref}_t \geq | p^{ac}_t - \hat{p}_t | \f$ (the variable
  * "v_abs_refschd", with coefficient \f$ \sigma \f$) if there is no
  * ReferenceSchedule, the reactive term if the UCBlock has no reactive
  * power, and the investment term if InvestmentCost is zero. The costs
  * \f$ c^{pr}_t \f$ and \f$ c^{sc}_t \f$ are 0 when they are not given,
  * and the term of a reserve whose Variable exist is in the Objective if
  * its cost is nonzero at some instant or if the Configuration of the
  * Objective asks for it: \p objc, or, if it is not a
  * SimpleConfiguration< int >, the f_objective_Configuration of the
  * BlockConfig, whose value has bit 0 set for the primary reserve and bit
  * 1 for the secondary one (default 0). A term with cost 0 is there for a
  * Solver that changes the coefficients of the Objective, say a
  * LagrangianDualSolver pricing the reserve requirements of the UCBlock;
  * without it, set_primary_spinning_reserve_cost() and
  * set_secondary_spinning_reserve_cost() refuse a nonzero cost once the
  * Objective is generated, as set_shutdown_costs() does for a unit whose
  * Objective has no shut-down term.
  *
  * Each term is the cost of one instant (see the class documentation),
  * and the start-up cost \f$ c^{su}_t \f$ depends on the instant of the
  * start-up only, not on how long the unit has been off before it. All the
  * formulations have the same Objective, with one exception: with PCuts the
  * quadratic term \f$ a_t ( p^{ac}_t )^2 \f$ is replaced by \f$ a_t z_t \f$
  * in the 3bin, T and pt formulations, by
  * \f$ a_t \sum_{ (h,k) } z^{hk}_t \f$ in the DP one,
  * \f$ a_t \sum_{ h } z^h_t \f$ in the SU one,
  * \f$ a_t \sum_{ k } \tilde{z}^k_t \f$ in the SD one and
  * \f$ a_t \vartheta_t \f$ in the SUSD one (see
  * generate_abstract_constraints()), the variables \f$ p^{ac}_t \f$ then
  * having no quadratic coefficient. The order of the Variable in the
  * DQuadFunction is that of (34): \f$ v \f$, \f$ w \f$ (if any),
  * \f$ p^{ac} \f$, \f$ u \f$, \f$ a^{ref} \f$ (if any), \f$ p^{pr} \f$,
  * \f$ p^{sc} \f$ (if any), the perspective variables (if any), \f$ q \f$
  * (if any), and \f$ x \f$ (if any) last. */

 void generate_objective( Configuration * objc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// ignore reserve netCDF variables when a ThermalUnitBlock is deserialized
 /** This function instructs the ThermalUnitBlock to ignore the reserve netCDF
  * variables, namely "PrimaryRho" and "SecondaryRho", when it is
  * deserialized. */

 static void ignore_reserve( void ) {
  f_ignore_netcdf_vars |= 1;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------- Methods for checking the ThermalUnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the ThermalUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this ThermalUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this ThermalUnitBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this ThermalUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  * - If \p fsbc is not a nullptr and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr and it is a
  *   SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *   fsbc->f_value.first and the type of violation is determined by
  *   fsbc->f_value.second (any nonzero number for relative violation and
  *   zero for absolute violation);
  *
  * - Otherwise, if both f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of the
  *   parameters are obtained analogously as above;
  *
  * - Otherwise, by default, the tolerance is Block::DefaultFeasTol and the
  *   relative violation is considered.
  *
  * This function considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed (Constraint::compute()).
  *
  * @param useabstract This parameter is ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns true if the schedule in the Solution is feasible
 /** Returns true if the schedule that the given :UnitBlockSolution holds,
  * i.e. the commitment, the active power and the spinning reserves of the
  * unit, is feasible: the values are read out of it and checked against the
  * data of the unit, i.e. the operational bounds of the power, the ramps
  * with the limits of the start-up and of the shut-down, the minimum up and
  * down times and the state the unit comes from, so that the Variable of the
  * ThermalUnitBlock are neither needed nor touched
  * [see Block::is_sol_feasible()]. Those are the constraints of the unit for
  * an integral commitment, which is what the formulations of it encode; a
  * commitment that is not integral is therefore not declared feasible, as a
  * Solution that says it holds a direction is not, the feasible region of a
  * unit being bounded. A unit that carries something the schedule does not
  * answer for is left to the check of the base class
  * [see is_sol_feasible_physical()]. */

 bool is_sol_feasible( Solution * sol ,
		       Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// true if is_sol_feasible() reads the Solution and leaves the Variable be
 /** Returns true if is_sol_feasible() checks the schedule against the data of
  * the unit, which it does unless the unit carries something that the
  * schedule does not answer for: a dimensioning variable, the reactive
  * power, a reference schedule, a scale of its own or a fixed Variable of a
  * formulation that the schedule only implies. In those cases the check goes
  * through the abstract representation, and therefore through the Variable,
  * as the base class does it. A :ThermalUnitBlock with constraints of its
  * own that the schedule does not answer for says so here. */

 [[nodiscard]] bool is_sol_feasible_physical( void ) const override;


/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE ThermalUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ThermalUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of electrical generation units, i.e.:
 *
 * - fixed consumption when the unit is off;
 *
 * - the contribution to the inertia depending on the commitment status.
 * @{ */

 /// returns the initial power value
 double get_initial_power( void ) const { return( f_InitialPower ); }

 /// returns the init up and down time value
 int get_init_up_down_time( void ) const { return( f_InitUpDownTime ); }

 /// returns the minimum allowed up time value
 Index get_min_up_time( void ) const { return( f_MinUpTime ); }

 /// returns the minimum allowed down time value
 Index get_min_down_time( void ) const { return( f_MinDownTime ); }

 /// returns the investment cost
 double get_investment_cost( void ) const { return( f_InvestmentCost ); }

 /// returns the current cost of the design (investment) variable
 /** Returns the current cost of the design (investment) variable, in the same
  * unscaled units as get_investment_cost(): i.e. the variable's Objective
  * coefficient divided by the scale factor. This equals f_InvestmentCost when
  * the Objective is in its "original" state, but it may differ if a dualizing
  * Solver has changed the coefficient (e.g. the non-anticipativity multiplier
  * in a nested Lagrangian); Solvers that consume the structural data (the DP
  * Solvers) must use this, rather than get_investment_cost(), as the actual
  * cost of building the unit. Returns 0 if the unit has no design variable,
  * and f_InvestmentCost if the Objective has not been generated yet. */
 double get_design_cost( void ) const;

 /// returns the installable capacity by the user
 double get_capacity( void ) const { return( f_Capacity ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nominal minimum active power output
 /** This method returns (a const reference to) the vector containing the
  * nominal minimum active power output of the unit for all time steps. When
  * the unit is available, get_min_power()[ t ] gives the minimum active power
  * output of the unit at time t in { 0, ..., get_time_horizon() - 1}. */

 const std::vector< double > & get_min_power( void ) const {
  return( v_MinPower );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the minimum power of the \p generator at time \p t

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( ( v_MinPower.size() > t ) ? v_MinPower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nominal maximum active power output
 /** This method returns (a const reference to) the vector containing the
  * nominal maximum active power output of the unit for all time steps. When
  * the unit is available, get_max_power()[ t ] gives the maximum active power
  * output of the unit at time t in { 0 , ..., get_time_horizon() - 1 }. */

 const std::vector< double > & get_max_power( void ) const {
  return( v_MaxPower );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the maximum power of \p generator at time \p t

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( ( v_MaxPower.size() > t ) ? v_MaxPower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum reactive power of \p generator at time t

 double get_min_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MinReactivePower.size() > t ) ? v_MinReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of \p generator at time t

 double get_max_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MaxReactivePower.size() > t ) ? v_MaxReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the commitment-gated minimum reactive power coefficient
 /** The reactive power bound depends on the commitment \f$ u_t \f$:
  * \f[
  *   Q^{mn}_t + Q^{mn,on}_t u_t \leq q_t \leq Q^{mx}_t + Q^{mx,on}_t u_t \; ,
  * \f]
  * where get_min_reactive_power() and get_max_reactive_power() return
  * \f$ Q^{mn}_t \f$ and \f$ Q^{mx}_t \f$ (the bounds when the unit is
  * off) and this method returns \f$ Q^{mn,on}_t \f$. Zero (the default)
  * gives the plain box \f$ [ Q^{mn}_t , Q^{mx}_t ] \f$ (see
  * generate_abstract_constraints()). */

 double get_min_reactive_power_on( Index t , Index generator = 0 ) const {
  return( ( v_MinReactivePowerOn.size() > t ) ?
	  v_MinReactivePowerOn[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the commitment-gated maximum reactive power coefficient
 /** The coefficient \f$ Q^{mx,on}_t \f$ of \f$ u_t \f$ in the reactive
  * bound; see get_min_reactive_power_on(). */

 double get_max_reactive_power_on( Index t , Index generator = 0 ) const {
  return( ( v_MaxReactivePowerOn.size() > t ) ?
	  v_MaxReactivePowerOn[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the operational minimum active power output at the time t
 /** This method returns the operational minimum active power output of the
  * unit at time \p t. See get_availability() for the definition of
  * operational minimum power.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The operational minimum active power output of the unit at the
  *         given time. */

 double get_operational_min_power( Index t ) const {
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_operational_min_power: "
			    "invalid time index " + std::to_string( t ) ) );
  return( compute_operational_min_power( v_MinPower[ t ] ,
                                         get_availability( t ) ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the operational maximum active power output at time \p t
 /** This method returns the operational maximum active power output of the
  * unit at time \p t. See get_availability() for the definition of
  * operational maximum power.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The operational maximum active power output of the unit at the
  *         given time. */

 double get_operational_max_power( Index t ) const {
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_operational_max_power: "
			    "invalid time index " + std::to_string( t ) ) );
  return( compute_operational_max_power( v_MaxPower[ t ] ,
                                         get_availability( t ) ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the availability of the unit at all time instants
 /** This method returns (a const reference to) the vector containing the
  * availability of the unit at all time instants. For each t in { 0 , ... ,
  * get_time_horizon() - 1 }, get_availability()[ t ] is the availability of
  * the unit at time t, which is a number between 0 and 1. When the
  * availability of the unit is zero, the unit is not under operation (for
  * instance, due to an outage or maintenance). When the availability of the
  * unit is 1, it is fully available and operating at maximum capacity.
  *
  * The availability of the unit determines its operational minimum and
  * maximum active power output, i.e., the effective bounds on the active
  * power output under which the unit operates. For each t in {0, ...,
  * get_time_horizon() - 1}, the operational maximum and minimum power are
  * \f$ P^{mx}_t = \chi_t \, \mathrm{MaxPower}[ t ] \f$ and
  * \f$ P^{mn}_t = \mathrm{MinPower}[ t ] \f$ if \f$ \chi_t > 0 \f$,
  * \f$ P^{mn}_t = 0 \f$ if \f$ \chi_t = 0 \f$, where \f$ \chi_t \f$ is the
  * availability and MinPower[ t ], MaxPower[ t ] are the nominal bounds
  * returned by get_min_power() and get_max_power() (see deserialize()). */

 const std::vector< double > & get_availability( void ) const {
  return( v_Availability );
  }

/*--------------------------------------------------------------------------*/
 /// returns the availability of the unit at time \p t

 double get_availability( Index t ) const {
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_availability: invalid "
                            "time index " + std::to_string( t ) ) );
  if( v_Availability.empty() )
   return( 1.0 );
  return( v_Availability[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary rho
 /** The returned vector contains the primary rho at each time; it is
  * empty if the datum is not given (or, for a rho, if it is zero at every
  * instant), and otherwise its size is get_time_horizon(). */

 const std::vector< double > & get_primary_rho( void ) const {
  return( v_PrimaryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary rho
 /** The returned vector contains the secondary rho at each time; it is
  * empty if the datum is not given (or, for a rho, if it is zero at every
  * instant), and otherwise its size is get_time_horizon(). */

 const std::vector< double > & get_secondary_rho( void ) const {
  return( v_SecondaryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary spinning-reserve cost
 /** The cost \f$ c^{pr}_t \f$ of the primary reserve, which, unlike the
  * participation factor (get_primary_rho()), may be modified (e.g., by a
  * Lagrangian price). Empty if it is zero at every instant, which is its
  * default. */

 const std::vector< double > & get_primary_spinning_reserve_cost( void )
  const {
  return( v_PrimarySpinningReserveCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary spinning-reserve cost
 /** The objective coefficient on the secondary spinning-reserve variables;
  * see get_primary_spinning_reserve_cost(). */

 const std::vector< double > & get_secondary_spinning_reserve_cost( void )
  const {
  return( v_SecondarySpinningReserveCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of delta ramp-up
 /** The returned vector contains the delta ramp-up at each time; it is
  * empty if the datum is not given (or, for a rho, if it is zero at every
  * instant), and otherwise its size is get_time_horizon(). */

 const std::vector< double > & get_delta_ramp_up( void ) const {
  return( v_DeltaRampUp );
  }

/*--------------------------------------------------------------------------*/
 /// returns the delta ramp-up at the time \p t
 /** This function return the delta ramp-up at time \p t, which is assumed
  * to be between 0 and get_time_horizon() - 1, i.e., the maximum increase
  * of the active power from t - 1 to t (from InitialPower if t == 0). */

 double get_delta_ramp_up( Index t ) const {
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_delta_ramp_up: invalid "
                            "time index " + std::to_string( t ) ) );
  if( v_DeltaRampUp.empty() )
   return( get_max_power( t ) );
  return( v_DeltaRampUp[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of delta ramp-down
 /** The returned vector contains the delta ramp-down at each time; it is
  * empty if the datum is not given (or, for a rho, if it is zero at every
  * instant), and otherwise its size is get_time_horizon(). */

 const std::vector< double > & get_delta_ramp_down( void ) const {
  return( v_DeltaRampDown );
  }

/*--------------------------------------------------------------------------*/
 /// returns the delta ramp-down at the time \p t
 /** This function return the delta ramp-down at the time \p t, which is
  * assumed to be between 0 and get_time_horizon() - 1, i.e., the maximum
  * decrease of the active power from t - 1 to t (from InitialPower if
  * t == 0). */

 double get_delta_ramp_down( Index t ) const {
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_delta_ramp_down: invalid "
                            "time index " + std::to_string( t ) ) );
  if( v_DeltaRampDown.empty() )
   return( get_max_power( t ) );
  return( v_DeltaRampDown[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the increase the ramp-up limits allow from \p from - 1 to \p to
 /** Returns \f$ \sum_{r = from}^{to} \Delta^+_r \f$, i.e., the
  * maximum increase of the active power over the steps from \p from - 1
  * (InitialPower if \p from == 0) to \p to, 0 if \p from > \p to. */

 double ramp_up_sum( Index from , Index to ) const {
  double sum = 0;
  for( Index t = from ; t <= to ; ++t )
   sum += get_delta_ramp_up( t );
  return( sum );
  }

/*--------------------------------------------------------------------------*/
 /// returns the decrease the ramp-down limits allow from \p from - 1 to \p to
 /** Returns \f$ \sum_{r = from}^{to} \Delta^-_r \f$, i.e., the
  * maximum decrease of the active power over the steps from \p from - 1
  * (InitialPower if \p from == 0) to \p to, 0 if \p from > \p to. */

 double ramp_down_sum( Index from , Index to ) const {
  double sum = 0;
  for( Index t = from ; t <= to ; ++t )
   sum += get_delta_ramp_down( t );
  return( sum );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of quadratic term
 /** Returns the coefficients \f$ a_t \f$ of the quadratic term of the cost,
  * of size get_time_horizon(), one value per instant (deserialize() expands
  * a datum given by intervals, and fills an absent one with zeros). */

 const std::vector< double > & get_quad_term( void ) const {
  return( v_QuadTerm );
  }

/*--------------------------------------------------------------------------*/
/// returns the coefficient of the quadratic term of the power cost function
/** This function returns the coefficient of the quadratic term of the
 * quadratic function that represents the cost of the power produced by the
 * unit at the time \p t.
 *
 * @param t A time instant between 0 and get_time_horizon() - 1.
 *
 * @return The coefficient of the quadratic term of the quadratic function
 *         that represents the cost of the power produced by the unit at the
 *         given time instant. */

 double get_quad_term( Index t ) const {
  if( v_QuadTerm.empty() )
   return( 0 );
  if( v_QuadTerm.size() == 1 )
   return( v_QuadTerm.front() );
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_quad_term: invalid "
                            "time index " + std::to_string( t ) ) );
  return( v_QuadTerm[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of linear term
 /** Returns the coefficients \f$ b_t \f$ of the linear term of the cost,
  * of size get_time_horizon(), one value per instant (deserialize() expands
  * a datum given by intervals, and fills an absent one with zeros). */

 const std::vector< double > & get_linear_term( void ) const {
  return( v_LinearTerm );
 }

/*--------------------------------------------------------------------------*/
/// returns the coefficient of the linear term of the power cost function
/** This function returns the coefficient of the linear term of the quadratic
 * function that represents the cost of the power produced by the unit at the
 * time \p t.
 *
 * @param t A time instant between 0 and get_time_horizon() - 1.
 *
 * @return The coefficient of the linear term of the quadratic function that
 *         represents the cost of the power produced by the unit at the given
 *         time instant. */

 double get_linear_term( Index t ) const {
  if( v_LinearTerm.empty() )
   return( 0 );
  if( v_LinearTerm.size() == 1 )
   return( v_LinearTerm.front() );
  assert( v_LinearTerm.size() == f_time_horizon );
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_linear_term: invalid "
                            "time index " + std::to_string( t ) ) );
  return( v_LinearTerm[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the (dualized) linear cost coefficient of the reactive power
 /** The reactive power variable q[t] has no cost in the unit's own objective,
  * but a Solver that dualizes a constraint involving it (e.g. a
  * LagrangianDualSolver dualizing the reactive-power node balance) injects a
  * linear coefficient on q[t] into the objective; set_reactive_linear_term()
  * mirrors it here. The vector follows the same empty/size-1/size-horizon
  * convention as get_linear_term(). */

 const std::vector< double > & get_reactive_linear_term( void ) const {
  return( v_ReactiveLinearTerm );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of constant term
 /** Returns the coefficients \f$ c_t \f$ of the commitment in the cost,
  * of size get_time_horizon(), one value per instant (deserialize() expands
  * a datum given by intervals, and fills an absent one with zeros). */

 const std::vector< double > & get_const_term( void ) const {
  return( v_ConstTerm );
  }

/*--------------------------------------------------------------------------*/
/// returns the constant term of the power cost function
/** This function returns the constant term of the function that represents
 * the cost of the power produced by the unit at the given time instant. This
 * is the fixed cost incurred when the unit is committed at time \p t.
 *
 * @param t A time instant between 0 and get_time_horizon() - 1.
 *
 * @return The fixed cost when the unit is committed at the given time
 *         instant. */

 double get_const_term( Index t ) const {
  if( v_ConstTerm.empty() )
   return( 0 );
  if( v_ConstTerm.size() == 1 )
   return( v_ConstTerm.front() );
  assert( v_ConstTerm.size() == f_time_horizon );
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_const_term: invalid "
                            "time index " + std::to_string( t ) ) );
  return( v_ConstTerm[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start-up costs
 /** Returns the start-up costs \f$ c^{su}_t \f$,
  * of size get_time_horizon(), one value per instant (deserialize() expands
  * a datum given by intervals, and fills an absent one with zeros). */

 const std::vector< double > & get_start_up_cost( void ) const {
  return( v_StartUpCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the reference schedule of the unit, if it has one
 /** Returns the vector of the power the unit is asked to follow, empty if it
  * has none; the deviation from it is a term of the Objective, weighed with
  * the scale factor [see generate_objective()]. */

 const std::vector< double > & get_reference_schedule( void ) const {
  return( v_RefSchedule );
  }

/*--------------------------------------------------------------------------*/
 /// tells whether the unit produces its maximum power at every instant
 /** Returns true if "FixToMaximum" is positive [see deserialize()], i.e., if
  * the rows \f$ p^{ac}_t \geq P^{mx}_t \f$ ("FixedGeneration") are there
  * [see generate_abstract_constraints()]. */

 bool is_fixed_to_maximum( void ) const { return( f_fixToMax > 0 ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of shut-down costs
 /** The returned vector contains the shut-down cost \f$ c^{sd}_t \f$, empty
  * if the unit pays nothing to shut down and of get_time_horizon() elements
  * otherwise. The cost is paid at
  * the instant the unit goes off. */

 const std::vector< double > & get_shut_down_cost( void ) const {
  return( v_ShutDownCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of fixed consumption
 /** Returns a pointer to the fixed consumption \f$ P^{au}_t \f$ of the
  * unit when it is off (see deserialize()), one value per instant, or
  * nullptr if the unit has none. */

 const double * get_fixed_consumption( Index generator ) const override {
  if( v_FixedConsumption.empty() )
   return( nullptr );
  return( &( v_FixedConsumption.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inertia commitment
 /** Returns a pointer to the inertia \f$ h^u_t \f$ the unit gives when it
  * is on (see deserialize()), one value per instant, or nullptr if the
  * unit gives none. */

 const double * get_inertia_commitment( Index generator ) const override {
  if( v_InertiaCommitment.empty() )
   return( nullptr );
  return( &( v_InertiaCommitment.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the start-up limit

 const std::vector< double > & get_start_up_limit( void ) const {
  return( v_StartUpLimit );
  }

/*--------------------------------------------------------------------------*/
 /// returns the shut-down limit

 const std::vector< double > & get_shut_down_limit( void ) const {
  return( v_ShutDownLimit );
  }

/*--------------------------------------------------------------------------*/
 /// returns the scale factor

 double get_scale( void ) const override { return( f_scale ); }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Variable OF THE ThermalUnitBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the ThermalUnitBlock
 *
 * These methods allow to read the each group of Variable that any
 * ThermalUnitBlock in principle has (although some may not):
 *
 * - commitment variables;
 *
 * - active and reactive power variables;
 *
 * - primary spinning reserve variables;
 *
 * - secondary spinning reserve variables;
 *
 * - start up variables;
 *
 * - shut down variables.
 * @{ */

 /// returns the vector of commitment variables

 bool has_commitment( void ) const override { return( true ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ColVariable * get_commitment( Index generator ) override {
  if( v_commitment.empty() )
   return( nullptr );
  return( &( v_commitment.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active_power variables

 ColVariable * get_active_power( Index generator ) override {
  if( v_active_power.empty() )
   return( nullptr );
  return( &( v_active_power.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of reactive power variables

 ColVariable * get_reactive_power( Index generator ) override {
  if( v_reactive_power.empty() )
   return( nullptr );
  return( &( v_reactive_power.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary_spinning_reserve variables

 /// the unit provides the reserve only if it has any to give
 /** The enclosing UCBlock asking for the reserve is not enough, the unit has
  * to have some to give: this is the very condition with which the Variable
  * are generated, said in terms of the data alone. */

 bool has_primary_reserve( void ) const override {
  return( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_secondary_reserve( void ) const override {
  return( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( v_primary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_primary_spinning_reserve.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary_spinning_reserve variables
 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( v_secondary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_secondary_spinning_reserve.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start_up variables, or nullptr if not defined
 ColVariable * get_start_up( void ) {
  if( v_start_up.empty() )
   return( nullptr );
  return( &( v_start_up.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// how many start-up (and shut-down) Variable this unit has
 /** The start-up and shut-down Variable are as many as the time instants in
  * which the unit can change state, which is fewer than the time horizon
  * whenever the initial state forces the first ones [see init_t]. */

 [[nodiscard]] Index get_number_start_up( void ) const {
  return( v_start_up.size() );
 }

/*--------------------------------------------------------------------------*/
 /// the start-up and shut-down indicators implied by a commitment profile
 /** Fills su and sd, both sized get_number_start_up(), with the indicators
  * that the given commitment profile implies: a start-up wherever the unit
  * goes off to on, a shut-down wherever it goes the other way, the first
  * instant being decided by the state the unit was in before the horizon.
  *
  * This is the same rule set_solution() applies to the Variable, in the
  * form a Solver filling a Solution needs: the indicators have to be saved
  * rather than derived later, because deriving them from a commitment that
  * is the average of several schedules gives the start-ups of the average,
  * which are fewer than the average of the start-ups. */

 void derive_start_up( const std::vector< double > & u ,
                       std::vector< double > & su ,
                       std::vector< double > & sd ) const;

/*--------------------------------------------------------------------------*/
 /// returns the vector of shut_down variables, or nullptr if not defined
 ColVariable * get_shut_down( void ) {
  if( v_shut_down.empty() )
   return( nullptr );
  return( &( v_shut_down.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the shut_down variable for time t, or nullptr if not defined
 ColVariable * get_shut_down( Index t ) {
  if( v_shut_down.empty() || ( t < init_t ) )
   return( nullptr );
  return( &( v_shut_down[ t - init_t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the design binary variable

 ColVariable & get_design( void ) { return( design ); }

/*--------------------------------------------------------------------------*/
 /// returns the const design binary variable

 const ColVariable & get_const_design( void ) const { return( design ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of commitment_plus variables (DP / SU / SD / pt /
 /// SUSD formulations), or nullptr if not defined

 ColVariable * get_commitment_plus( void ) {
  if( v_commitment_plus.empty() )
   return( nullptr );
  return( &( v_commitment_plus.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of the variables \f$ z_t \f$ of the perspective
 /// cuts, which every formulation with PCuts has, or nullptr

 ColVariable * get_cut( void ) {
  if( v_cut.empty() )
   return( nullptr );
  return( &( v_cut.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of (h, k) index pairs of the y^+ commitment
 /// variables for the DP, pt, SU, SD and SUSD formulations

 const std::vector< std::pair< Index , Index > > & get_Y_plus( void ) const {
  return( v_Y_plus );
  }

/*--------------------------------------------------------------------------*/
 /// returns the formulation code currently in use; the value is one of
 /// tbinForm, TForm, ptForm, DPForm, SUForm, SDForm, SUSDForm

 unsigned char get_formulation( void ) const { return( AR & FormMsk ); }

/*--------------------------------------------------------------------------*/
 /// returns true iff the formulation uses perspective cuts (PCuts bit)

 bool has_perspective_cuts( void ) const { return( AR & PCuts ); }

/*--------------------------------------------------------------------------*/
 /// completes a schedule given by the commitment and the active power
 /** Given the commitment \f$ u_t \f$ and the active power \f$ p^{ac}_t \f$
  * already written in their Variable (and the reserves, if any), sets every
  * other Variable of the formulation consistently with them, so that the
  * point is feasible whenever the schedule is and costs what the schedule
  * costs:
  *
  * - the start-up and shut-down indicators \f$ v_t \f$ and \f$ w_t \f$,
  *   from the changes of \f$ u \f$ and the state \f$ u_{-1} \f$ before the
  *   horizon;
  *
  * - the deviation \f$ a^{ref}_t = | p^{ac}_t - \hat{p}_t | \f$ from the
  *   reference schedule, if any;
  *
  * - with PCuts, \f$ z_t = ( p^{ac}_t )^2 / u_t \f$ (0 if \f$ u_t = 0 \f$),
  *   the perspective of the square, which satisfies every cut (33) and
  *   gives \f$ a_t z_t = a_t ( p^{ac}_t )^2 \f$ for an integer commitment
  *   and, for a fractional one (a convex combination of schedules), at most
  *   the combination of their quadratic costs;
  *
  * - in the pt, DP, SU, SD and SUSD formulations, if the commitment is
  *   integer, the arc variables of the path of the schedule (1 on its runs
  *   on and off, the run that the end of the horizon cuts being the arc
  *   to \f$ T + 1 \f$, and 0 on the other arcs), the power of each run on
  *   its own Variable (\f$ p^{hk}_t \f$, \f$ p^h_t \f$, \f$ \tilde{p}^k_t
  *   \f$ equal to \f$ p^{ac}_t \f$ at the instants of the run and 0
  *   elsewhere) and, with PCuts, its square on the variables of the cuts of
  *   the run and on \f$ \vartheta_t \f$ [see set_path_solution()]; a
  *   fractional commitment does not identify a path, and these Variable
  *   are then left as they are.
  *
  * It is used by the Solvers of a ThermalUnitBlock that compute only
  * \f$ ( p^{ac} , u ) \f$, such as the dynamic programming ones, to write
  * a complete solution at the end of compute(), and by
  * ThermalUnitBlockSolution::write(), which then writes back the saved
  * start-up and shut-down indicators; hence a Solution saved and written
  * back gives a feasible point of the same value in every formulation,
  * except that with PCuts the value is that of the quadratic cost, which
  * is at least the one of the cuts that a :MILPSolver minimizes. A
  * schedule that has no path (e.g., a run shorter than the minimum up
  * time) leaves its arc variables all 0, and one that has a path the rows
  * forbid (e.g., a shut-down at 0 with \f$ p_{-1} > P^{sd}_0 \f$, see
  * generate_abstract_constraints()) is written as it is. */

 virtual void set_solution( void );

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution storing the state of this ThermalUnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this ThermalUnitBlock. This
  * is a ThermalUnitBlockSolution extending UnitBlockSolution with the
  * specific extra solution information of ThermalUnitBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - the first four bits (bit 0 to bit 3) are "taken" by the base
  *   UnitBlock[Solution]
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 15 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [Thermal]UnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE ThermalUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the ThermalUnitBlock
 * @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * ThermalUnitBlock. See ThermalUnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR INITIALIZING THE ThermalUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the ThermalUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "ThermalUnitBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
 /** Methods for changing the data of the ThermalUnitBlock.
  * @{ */

 /** This method has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the ThermalUnitBlock, and
  * "translate" them into both changes of the actual data structures and
  * corresponding "physical Modification". These Modification are those
  * for which Modification::concerns_Block() is true.
  *
  *     THE IMPLEMENTATION OF THIS METHOD IS BOTH PARTIAL AND HORRIBLE,
  *     ONE SINGLE ABSTRACT MODIFICATION CAN GIVE RISE TO MANY MANY MANY
  *     PHYSICAL ONES, IT SHOULD BE COMPLETELY OVERHAULED!!! */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// update the availability of the unit
 /** This method updates the availability of the unit. The \p subset parameter
  * contains a list of time instants and \p values contains the availability
  * of the unit at those time instants. The availability of the unit at time
  * subset[ i ] is given by std::next( values , i ) for each i in {0, ...,
  * subset.size() - 1}, also if \p subset is not ordered (\p ordered being
  * false), and an instant given more than once takes the last of its
  * values.
  *
  * The new values are checked before anything changes, as deserialize()
  * checks the data [see check_power_limits()]: each has to be between 0
  * and 1, the operational bounds it gives have to satisfy \f$ P^{mn}_t
  * \leq P^{mx}_t \f$, and a StartUpLimit or ShutDownLimit given in the
  * data has to lie between them; std::logic_error is thrown otherwise. A
  * StartUpLimit or ShutDownLimit not given in the data is the operational
  * minimum power, which it follows [see set_default_limits()], so that an
  * unavailable unit has 0 for both and cannot produce; the numbers of ramp
  * steps of the SUSD formulation are computed again [see
  * compute_ramp_steps()].
  *
  * Once the Constraint are generated, and unless \p issueAMod is eDryRun,
  * the abstract representation of every formulation follows [see
  * update_rows()]: the coefficients and the sides of the rows in which the
  * operational bounds and the limits are (the ramps, the minimum and the
  * maximum power, the reserves, the Perspective Cuts at the bounds, the
  * box of the active power, FixToMaximum, the shut-down at 0 and, for a
  * NuclearUnitBlock, the rows of its operating rules that contain them)
  * change, and the rows (18), (25) and (26) of
  * generate_abstract_constraints() whose existence depends on the bounds
  * are added or removed, all the Modification in a single GroupModification
  * (nested in the channel of \p issueAMod, if any). If this cannot be done
  * (e.g., the new availability makes \f$ p_{-1} - \Delta^-_0 > P^{mx}_0
  * \f$, which generate_abstract_constraints() refuses), std::logic_error
  * is thrown and the data are restored as they were. The physical
  * Modification is issued last, according to \p issuePMod. */

 void set_availability( MF_dbl_it values ,
                        Subset && subset ,
                        const bool ordered = false ,
                        ModParam issuePMod = eNoBlck ,
                        ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the availability of the unit
 /** This method updates the availability of the unit. The \p rng parameter
  * contains a range of time instants and \p values contains the availability
  * of the unit at those time instants. The availability of the unit at time
  * rng.first + i is given by std::next( values , i ) for each i in {0, ...,
  * ( std::min( rng.second, get_time_horizon() ) - rng.first - 1 )}. The
  * checks, the default limits and the changes of the abstract
  * representation are those of the Subset version. */

 void set_availability( MF_dbl_it values , Range rng = INFRange ,
                        ModParam issuePMod = eNoBlck ,
                        ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the maximum power of the unit
 /** Sets MaxPower[ t ] for the instants t in \p subset to the values
  * pointed by \p values, in the order of \p subset also if it is not
  * ordered (\p ordered being false). The new values are checked before
  * anything changes [see check_power_limits()]: each has to be at least
  * MinPower[ t ], and a StartUpLimit or ShutDownLimit given in the data has
  * to lie between the operational bounds it gives; std::logic_error is
  * thrown otherwise. The numbers of ramp steps of the SUSD formulation are
  * computed again [see compute_ramp_steps()], and once the Constraint are
  * generated the abstract representation of every formulation follows, as
  * for the availability [see set_availability()], the data being restored
  * if it cannot. */

 void set_maximum_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the maximum power of the unit
 /** As the Subset version, for the instants in \p rng. */

 void set_maximum_power( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_startup_costs( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_startup_costs( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the shut-down cost at the instants of \p subset
 /** Sets \f$ c^{sd}_t \f$ for \f$ t \in \f$ \p subset to the values from
  * \p values on; the instants before \f$ t_0 \f$, which have no shut-down
  * variable, cannot be in \p subset. Once the Objective is generated
  * without the shut-down term (ShutDownCost absent or zero, see
  * generate_objective()) only zero costs are accepted, and a nonzero one is
  * refused with std::logic_error before anything changes; otherwise the
  * coefficients of the Objective, \f$ \sigma c^{sd}_t \f$, follow. The same
  * rules hold for the Range version. */

 void set_shutdown_costs( MF_dbl_it values ,
                          Subset && subset ,
                          const bool ordered = false ,
                          ModParam issuePMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_shutdown_costs( MF_dbl_it values , Range rng = INFRange ,
                          ModParam issuePMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_const_term( MF_dbl_it values ,
                      Subset && subset ,
                      const bool ordered = false ,
                      ModParam issuePMod = eNoBlck ,
                      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_const_term( MF_dbl_it values , Range rng = INFRange ,
                      ModParam issuePMod = eNoBlck ,
                      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_linear_term( MF_dbl_it values ,
                       Subset && subset ,
                       const bool ordered = false ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_linear_term( MF_dbl_it values , Range rng = INFRange ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_quad_term( MF_dbl_it values ,
                     Subset && subset ,
                     const bool ordered = false ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_quad_term( MF_dbl_it values , Range rng = INFRange ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/// set the linear cost coefficient of the reactive power variables
/** The reactive power variables q[t] carry no cost in the unit's own
 * objective, but they appear in the objective when a Solver dualizes a
 * constraint they belong to (e.g. the reactive-power node balance in a
 * LagrangianDualSolver): the dual term is a linear coefficient on q[t]. These
 * setters store that coefficient into v_ReactiveLinearTerm so that the DP
 * solvers can price q[t] over its [Qmin,Qmax] box; they are the reactive
 * counterpart of set_linear_term(). */

 void set_reactive_linear_term( MF_dbl_it values ,
                                Subset && subset ,
                                const bool ordered = false ,
                                ModParam issuePMod = eNoBlck ,
                                ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_reactive_linear_term( MF_dbl_it values , Range rng = INFRange ,
                                ModParam issuePMod = eNoBlck ,
                                ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the cost of the primary reserve at the instants of \p subset
 /** Sets \f$ c^{pr}_t \f$ for \f$ t \in \f$ \p subset to the values from
  * \p values on. Before the Variable are generated the cost is stored
  * whatever the reserve Variable will be (the enclosing UCBlock may ask for
  * them later, see UnitBlock::set_reserve_vars()); once they are generated,
  * a unit without primary reserve Variable ignores the call. Once the
  * Objective is generated without the primary reserve term [see
  * generate_objective()] only zero costs are accepted, and a nonzero one is
  * refused with std::logic_error before anything changes; otherwise the
  * coefficients of the Objective, \f$ \sigma c^{pr}_t \f$, follow. The same
  * rules hold for the Range version and for the secondary reserve. */

 void set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                         Subset && subset ,
                                         const bool ordered = false ,
                                         ModParam issuePMod = eNoBlck ,
                                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                         Range rng = INFRange ,
                                         ModParam issuePMod = eNoBlck ,
                                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                           Subset && subset ,
                                           const bool ordered = false ,
                                           ModParam issuePMod = eNoBlck ,
                                           ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                           Range rng = INFRange ,
                                           ModParam issuePMod = eNoBlck ,
                                           ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given \p subset contains the 0 index, this function sets the
  * initial power. If the given \p subset does not contain the index 0, this
  * function does nothing. Since \p subset can have multiple zeros, only the
  * last one is considered, which means that the value for the initial power
  * will be that in the vector pointed by \p it associated with this last
  * zero.
  *
  * The initial power enters the abstract representation only if the unit
  * is on before the horizon, and then the rows that contain it follow
  * [see update_initial_power_in_cnstrs()]: the ramp rows at 0 and the
  * deliverability rows of the reserves at 0 (31'), the row of the
  * shut-down at 0, which is allowed iff \f$ p_{-1} \leq P^{sd}_0 \f$,
  * and, in the pt, DP, SU, SD and SUSD formulations, the coefficients of
  * \f$ \psi \f$ and the ramp rows of several steps from the initial
  * power of the SUSD formulation, which are added or removed. If this
  * cannot be done (e.g., the new value makes \f$ p_{-1} + \Delta^+_0 <
  * P^{mn}_0 \f$, which generate_abstract_constraints() refuses, or it
  * moves the initial power of a NuclearUnitBlock to another band), the
  * method throws std::logic_error and leaves the initial power as it was.
  * The numbers of ramp steps from the initial power that the SUSD
  * formulation uses (MaxRampUpSteps[ 0 ], MaxRampDownSteps[ 0 ]) are
  * computed again. Before any change, the method throws std::logic_error
  * if the value is negative or if the unit is on before the horizon and
  * the value is below \f$ \hat P^{mn}_0 \f$, the two cases that
  * check_data_consistency() refuses (where deserialize() raises such a
  * value to \f$ \hat P^{mn}_0 \f$, the method does not), since with
  * \f$ p_{-1} < \hat P^{mn}_0 \f$ the formulations and the dynamic
  * programming Solvers would not agree on the shut-down at 0. */

 void set_initial_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given Range \p rng contains 0, this function sets the initial
  * power. In this case, if the first element of \p rng is 0, the initial
  * power will be set to the value pointed by the given iterator. In general,
  * the initial power will be the one found at position -rng.first in the
  * vector pointed by \p it if this Range contains the 0 index. If the given
  * Range \p rng does not contain the 0 index, this function does nothing.
  *
  * The initial power enters the abstract representation only if the unit
  * is on before the horizon, and then the rows that contain it follow
  * [see update_initial_power_in_cnstrs()]: the ramp rows at 0 and the
  * deliverability rows of the reserves at 0 (31'), the row of the
  * shut-down at 0, which is allowed iff \f$ p_{-1} \leq P^{sd}_0 \f$,
  * and, in the pt, DP, SU, SD and SUSD formulations, the coefficients of
  * \f$ \psi \f$ and the ramp rows of several steps from the initial
  * power of the SUSD formulation, which are added or removed. If this
  * cannot be done (e.g., the new value makes \f$ p_{-1} + \Delta^+_0 <
  * P^{mn}_0 \f$, which generate_abstract_constraints() refuses, or it
  * moves the initial power of a NuclearUnitBlock to another band), the
  * method throws std::logic_error and leaves the initial power as it was.
  * The numbers of ramp steps from the initial power that the SUSD
  * formulation uses (MaxRampUpSteps[ 0 ], MaxRampDownSteps[ 0 ]) are
  * computed again. Before any change, the method throws std::logic_error
  * if the value is negative or if the unit is on before the horizon and
  * the value is below \f$ \hat P^{mn}_0 \f$, the two cases that
  * check_data_consistency() refuses (where deserialize() raises such a
  * value to \f$ \hat P^{mn}_0 \f$, the method does not), since with
  * \f$ p_{-1} < \hat P^{mn}_0 \f$ the formulations and the dynamic
  * programming Solvers would not agree on the shut-down at 0. */

 void set_initial_power( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial up/down time
 /** Sets InitUpDownTime \f$ \tau_0 \f$ to the value associated with the
  * index 0 of \p subset (the last one, if 0 appears more than once), and
  * does nothing if 0 is not in \p subset. The initial state enters the
  * abstract representation only through whether the unit is on before the
  * horizon, \f$ \tau_0 > 0 \f$, and through the first instant
  * \f$ t_0 \f$ at which it may switch [see first_free_instant()]: these
  * two decide which Variable exist (the start-up and shut-down Variable
  * of the instants \f$ t \geq t_0 \f$, the intervals of the pt, DP, SU, SD
  * and SUSD formulations) and which are fixed, and the rows are written
  * from them. Once the Variable are generated (and the abstract
  * representation is to be changed), the value can therefore change only
  * within the values that give the same state and the same \f$ t_0 \f$,
  * i.e., any two values not smaller than MinUpTime, or any two not larger
  * than -MinDownTime, and only if the bound that deserialize() puts on
  * MinUpTime and MinDownTime does not change. A unit kept on or off for
  * the whole horizon by a minimum time cut to that bound is therefore
  * refused, since the minimum time of the data is unknown. The rows,
  * written anew [see update_rows()], are then those of a unit read with
  * the new value, and only the dynamic programming Solvers see the
  * change. Any other value
  * would make Variable appear or disappear: the method then throws
  * std::logic_error and changes nothing. Generating those Variable always
  * and fixing them, as is done for the interval \f$ ( 0 , 0 ) \f$ of the
  * shut-down at 0, would change the model of every unit with
  * \f$ t_0 > 0 \f$ (by \f$ 2 t_0 \f$ start-up and shut-down Variable
  * fixed to 0, and by the intervals of the other initial state), which is
  * not done. The numbers of ramp steps from the initial power of the SUSD
  * formulation are computed again. Before the Variable are generated any
  * value is taken; the bounds on MinUpTime and MinDownTime computed by
  * deserialize() are not computed again. A dry run of \p issuePMod changes
  * nothing, one of \p issueAMod the data only. */

 void set_init_updown_time( MF_int_it values ,
                            Subset && subset ,
                            const bool ordered = false ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial up/down time
 /** As the Subset version, the value being the one at the position
  * -rng.first of \p values if \p rng contains 0. */

 void set_init_updown_time( MF_int_it values , Range rng = INFRange ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the minimum up and down times, before the Variable are generated
 /** Sets the minimum up time and the minimum down time of the unit. Unlike
  * the other data, these two decide the structure of the model (how many
  * instants of the commitment the initial state fixes, and therefore how many
  * start-up and shut-down Variable there are), hence they can only be set
  * before the Variable are generated, and the method throws otherwise. Each
  * of them is taken between 1 and the time horizon plus the number of
  * instants the unit has been on (off) before it, or plus one if that is
  * smaller or the unit has been off (on), the bound saying that the unit,
  * being on (off) before the horizon, never switches within it [see
  * deserialize()]. The bound is that of the initial up/down time the unit
  * has when the method is called, which set_init_updown_time() does not
  * change afterwards. */

 void set_min_up_down_time( Index min_up_time , Index min_down_time ,
                            ModParam issuePMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the scale factor
 /** This method sets the scale factor.
  *
  * @param values An iterator to a vector containing the scale factor.
  *
  * @param subset If non-empty, the scale factor is set to the value pointed
  *               by \p values. If empty, no operation is performed.
  *
  * @param ordered This parameter is ignored.
  *
  * @param issuePMod Controls how physical Modification are issued.
  *
  * @param issueAMod Controls how abstract Modification are issued. */

 void scale( MF_dbl_it values ,
             Subset && subset ,
             const bool ordered = false ,
             c_ModParam issuePMod = eNoBlck ,
             c_ModParam issueAMod = eNoBlck ) override;

/*--------------------------------------------------------------------------*/

 // For the Range version, use the default implementation defined in UnitBlock
 using UnitBlock::scale;

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/
 /// a dynamic group of rows whose number depends on the data
 /** The rows that build_rows() writes with put_dyn_row(), each with the key
  * that identifies it among the rows of the group (the indices of the loops
  * that build it), in the same order as the rows. */
 struct DynRows {
  std::list< FRowConstraint > rows;          ///< the rows
  std::list< std::vector< int > > keys;      ///< the key of each row
  };

 /// computes the numbers of ramp steps of the SUSD formulation
 /** Fills v_MaxRampSteps and v_MaxRampDownSteps, i.e., \f$ J^+_t \f$ and
  * \f$ J^-_t \f$ of (18) in generate_abstract_constraints(), from the
  * operational bounds and the ramps, and then calls
  * compute_initial_ramp_steps() for the entries of the initial power; it is
  * called by deserialize() and whenever the availability or the maximum
  * power changes, before the rows are updated [see update_rows()]. */

 void compute_ramp_steps( void );

/*--------------------------------------------------------------------------*/
 /// checks the maximum power and the availability of an instant
 /** Throws std::logic_error, the message beginning with \p who, if the
  * availability \p availability is not between 0 and 1, if the maximum
  * power \p max_power is below MinPower[ \p t ], if the operational
  * minimum power that they give is above the operational maximum one, or
  * if a StartUpLimit or ShutDownLimit given in the data is not between the
  * two (a default one follows them, see set_default_limits()). */

 void check_power_limits( const std::string & who , Index t ,
                          double max_power , double availability ) const;

/*--------------------------------------------------------------------------*/
 /// the default StartUpLimit and ShutDownLimit
 /** Sets StartUpLimit[ t ] (ShutDownLimit[ t ]) to the operational minimum
  * power at t, at every instant, if it is not given in the data; called by
  * deserialize() and whenever the availability changes. */

 void set_default_limits( void );

/*--------------------------------------------------------------------------*/
 /// what set_maximum_power() and set_availability() have in common
 /** Sets MaxPower (Availability if \p availability) at the instants \p ts
  * to the values from \p values on, after having checked them with
  * check_power_limits(); then sets the default limits and the numbers of
  * ramp steps, and changes the abstract representation with update_rows()
  * if the Constraint are generated, restoring the data if it throws. The
  * changes are made according to \p issuePMod and \p issueAMod as in the
  * setters, which issue the physical Modification. Returns false if
  * nothing changes, true otherwise. */

 bool guts_of_set_power_limits( const std::string & who , bool availability ,
                                const Subset & ts , MF_dbl_it values ,
                                ModParam issuePMod , ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// the first instant at which the unit may switch, if InitUpDownTime were
 /// \p init
 /** The instant \f$ t_0 \f$ of generate_abstract_variables() for
  * InitUpDownTime \p init and the minimum times of the unit:
  * MinUpTime - \p init if the unit is on before the horizon for less than
  * MinUpTime instants, MinDownTime + \p init if it is off for less than
  * MinDownTime instants, 0 otherwise, and at most the time horizon. */

 Index first_free_instant( int init ) const;

/*--------------------------------------------------------------------------*/
 /// what the two set_init_updown_time() have in common
 /** Sets InitUpDownTime to \p value, with \p issuePMod and \p issueAMod as
  * in set_init_updown_time() [see there]. */

 void guts_of_set_init_updown_time( int value , ModParam issuePMod ,
                                    ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// writes the rows of the abstract representation from the current data
 /** Writes every row of the formulation, from the data as they are, through
  * size_rows(), put_row(), push_row(), put_dyn_row() and put_box(), and
  * registers the groups through add_rows() and add_dyn_rows(); the
  * ZOConstraints only if \p generate_ZOConstraints. Called with the rows
  * being generated by generate_abstract_constraints(), and with the rows
  * being compared by update_rows(), which then finds what differs from the
  * rows there are: the rows are written in one place only, and the update
  * equals a generation from the changed data. A derived class that has
  * rows of its own calls the method of the base class first and then
  * writes them with the same methods. */

 virtual void build_rows( bool generate_ZOConstraints );

/*--------------------------------------------------------------------------*/
 /// updates the rows of the abstract representation to the current data
 /** Writes the rows anew with build_rows() and compares them with those
  * that the abstract representation has: a static row whose coefficients
  * or sides differ gets the new ones, a dynamic row [see put_dyn_row()]
  * that is no longer there is removed and one that is new is added, all
  * the Modification being issued according to \p issueAMod in a single
  * GroupModification. A static row with other Variable than those it has, a
  * static group with another number of rows (also a group written by
  * size_rows() or push_row() to which no row would be written at all), or
  * a row of a dynamic group that the abstract representation does not
  * have, cannot be changed in place: the method then throws
  * std::logic_error before anything is changed, and so it does if
  * build_rows() throws, e.g., because the data violate a condition of the
  * generation. Nothing is done if the Constraint are not generated. */

 void update_rows( c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// true if build_rows() is generating the rows, false if comparing them

 bool generating_rows( void ) const { return( f_row_cmp == nullptr ); }

/*--------------------------------------------------------------------------*/
 /// the number of rows of a static group: \p n
 /** Resizes \p rows to \p n when generating; when comparing, throws
  * std::logic_error if \p rows does not have \p n rows. */

 void size_rows( std::vector< FRowConstraint > & rows , Index n );

/*--------------------------------------------------------------------------*/
 /// records \p rows among the groups that update_rows() checks

 void register_row_group( std::vector< FRowConstraint > & rows );

/*--------------------------------------------------------------------------*/
 /// the row \p i of a static group: \p lhs <= \p vars <= \p rhs
 /** Sets the row when generating; when comparing, records the coefficients
  * and the sides that differ from those of the row, after having checked
  * that the row has the Variable of \p vars in the same order (and thrown
  * std::logic_error if not). */

 void put_row( std::vector< FRowConstraint > & rows , Index i ,
               LinearFunction::v_coeff_pair && vars , double lhs ,
               double rhs );

/*--------------------------------------------------------------------------*/
 /// the next row of a static group: as put_row(), at the end of \p rows
 /** Appends the row when generating; when comparing, it is put_row() of
  * the next row of \p rows, and the number of rows written for each group
  * has then to be the one the group has, also when it is none. */

 void push_row( std::vector< FRowConstraint > & rows ,
                LinearFunction::v_coeff_pair && vars , double lhs ,
                double rhs );

/*--------------------------------------------------------------------------*/
 /// the row with the given \p key of a dynamic group
 /** Appends the row to \p rows when generating; when comparing, the rows
  * written are matched with those there are by the key: the row of a key
  * that is in both is compared as by put_row() (and removed and added
  * anew if it has other Variable), the row of a key that is only there is
  * removed, and the row of a key that is new is added. */

 void put_dyn_row( DynRows & rows , std::vector< int > && key ,
                   LinearFunction::v_coeff_pair && vars , double lhs ,
                   double rhs );

/*--------------------------------------------------------------------------*/
 /// the box \p lhs <= \p var <= \p rhs
 /** Sets the box when generating, records the sides that differ when
  * comparing. */

 void put_box( BoxConstraint & box , ColVariable * var , double lhs ,
               double rhs );

/*--------------------------------------------------------------------------*/
 /// registers a static group when generating, does nothing when comparing

 template< class C >
 void add_rows( C & rows , std::string && name ) {
  if( generating_rows() )
   add_static_constraint( rows , std::move( name ) );
  }

/*--------------------------------------------------------------------------*/
 /// registers a dynamic group when generating, does nothing when comparing

 void add_dyn_rows( DynRows & rows , std::string && name ) {
  if( generating_rows() ) {
   add_dynamic_constraint( rows.rows , std::move( name ) );
   f_dyn_groups.push_back( & rows );
   }
  }

/*--------------------------------------------------------------------------*/
 /// writes the path of the schedule in the path formulations
 /** Called by set_solution() in the pt, DP, SU, SD and SUSD formulations:
  * from the commitment \f$ u \f$ and the active power \f$ p^{ac} \f$, if
  * \f$ u \f$ is integer, sets to 1 the arc variables of the runs of the
  * schedule and to 0 the others, the power of each run to \f$ p^{ac}_t \f$
  * on its Variable and to 0 on those of the other runs, and the variables
  * of the perspective cuts of a run to \f$ ( p^{ac}_t )^2 \f$ (see
  * set_solution()). */

 void set_path_solution( void );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current initial power
 /** This function updates the abstract representation for the current
  * initial power, which it contains only if the unit is on before the
  * horizon: no Variable depends on it (the interval \f$ ( 0 , 0 ) \f$ of
  * the shut-down at 0 being there whatever its value), and the rows that
  * contain it, those of a derived class included, are written anew by
  * update_rows(), which throws std::logic_error, changing nothing, if
  * they cannot be updated in place. A derived class that has rows that
  * update_rows() cannot change (e.g., the bands at 0 of a
  * NuclearUnitBlock) checks them before calling it. */

 virtual void update_initial_power_in_cnstrs( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// computes v_MaxRampSteps[ 0 ] and v_MaxRampDownSteps[ 0 ]
 /** The number of ramp steps from the initial power depends on it and on
  * the initial up/down time, hence it is computed again whenever either of
  * them changes. */

 void compute_initial_ramp_steps( void );

/*--------------------------------------------------------------------------*/
 /// throws if set_initial_power() may not give \p value [see there]

 void check_initial_power( double value ) const;

/*--------------------------------------------------------------------------*/
 /// number of Variable a derived class appends to the Objective
 /** generate_objective() puts the Variable of the ThermalUnitBlock in the
  * Objective in a fixed order, the design variable (if any) being the last;
  * a derived class may append further Variable of its own after it, and it
  * must then return their number here, so that the sections of the
  * ThermalUnitBlock are found at the same indices. A change of the
  * coefficients of those Variable via the abstract representation is
  * routed to objective_tail_change(). */

 virtual Index objective_tail( void ) const { return( 0 ); }

/*--------------------------------------------------------------------------*/
 /// handles a change of the coefficients of the appended Variable
 /** Called by handle_objective_change() when a change of the coefficients
  * of the Objective reaches the Variable appended by a derived class [see
  * objective_tail()], [ \p first , \p last ) being the positions, among
  * them, of the changed ones (a Subset change may not change all of those
  * in between). The default throws, as the ThermalUnitBlock has no such
  * Variable. */

 virtual void objective_tail_change( const DQuadFunction * qf , Index first ,
                                     Index last ) {
  throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: the "
   "coefficients of the Variable appended to the Objective cannot change" ) );
  }

/*--------------------------------------------------------------------------*/
 /// the index of the first reactive power Variable in the Objective
 /** The reactive power Variable are the last section of the Objective before
  * the design variable (if any) and the Variable appended by a derived
  * class [see objective_tail()]. */

 Index reactive_objective_start( const DQuadFunction * qf ) const {
  return( qf->get_num_active_var() - objective_tail() -
          ( ( f_InvestmentCost != 0 ) ? 1 : 0 ) - f_time_horizon );
  }


/*--------------------------------------------------------------------------*/
 /// returns the operational minimum power
 /** This method computes the operational minimum power for the given nominal
  * minimum power and availability.
  *
  * @param nominal_min_power The nominal minimum power.
  *
  * @param availability A number between 0 and 1.
  *
  * @return The operational minimum power. */

 double compute_operational_min_power( double nominal_min_power ,
                                       double availability ) const {
  return( availability > 0.0 ? nominal_min_power : 0.0 );
 }

/*--------------------------------------------------------------------------*/
 /// returns the operational maximum power
 /** This method computes the operational maximum power for the given nominal
  * maximum power and availability.
  *
  * @param nominal_max_power The nominal maximum power.
  *
  * @param availability A number between 0 and 1.
  *
  * @return The operational maximum power. */

 double compute_operational_max_power( double nominal_max_power ,
                                       double availability ) const {
  return( nominal_max_power * availability );
 }

/*--------------------------------------------------------------------------*/
 /// updates the terms of the Objective associated with the start-up cost
 /** This method updates the terms of the Objective that are associated with
  * the start-up cost.
  *
  * @param subset A set of time instants at which the start-up costs must be
  *        updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_start_up( const Subset & subset ,
                                 c_ModParam issueAMod ) const;

/*--------------------------------------------------------------------------*/
 /// updates the terms of the Objective associated with the active power cost
 /** This method updates the terms of the Objective that are associated with
  * the active power cost.
  *
  * @param subset A set of time instants at which the active power costs must
  *        be updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_active_power( const Subset & subset ,
                                     c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the terms of the Objective associated with the fixed cost
 /** This method updates the terms of the Objective that are associated with
  * the fixed cost.
  *
  * @param subset A set of time instants at which the fixed costs must be
  *        updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_commitment( const Subset & subset ,
                                   c_ModParam issueAMod ) const;

/*--------------------------------------------------------------------------*/
 /// updates the remaining terms of the Objective that carry the scale
 /** This method updates the terms of the Objective associated with the
  * shut-down, the primary and secondary spinning reserves, the perspective
  * cuts and the reactive power, each being \f$ \sigma \f$ times the cost of
  * one copy of the unit (see UnitBlock::scale()).
  *
  * @param subset A set of time instants at which the terms must be updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_other_terms( const Subset & subset ,
                                    c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the coefficients of the Variable appended by a derived class
 /** Called by update_objective() after the scale factor has changed, so
  * that the coefficients of the Variable appended to the Objective by a
  * derived class [see objective_tail()] follow it at the time instants in
  * \p subset. The default does nothing, as the ThermalUnitBlock has no such
  * Variable. */

 virtual void update_objective_tail( const Subset & subset ,
                                     c_ModParam issueAMod ) { }

/*--------------------------------------------------------------------------*/
 /// updates the term of the Objective associated with the investment cost
 /** This method updates the coefficient of the design variable in the
  * Objective, i.e., the term \f$ \sigma c^{inv} x \f$, where \f$ \sigma \f$
  * is the current scale factor (see UnitBlock::scale()), \f$ c^{inv} \f$ the
  * investment cost of one module and \f$ x \f$ the design variable. It
  * is meant to be called after f_scale has been modified, so that the
  * Objective coefficient stays in sync with the scaled cost convention
  * used by the other update_objective_* helpers.
  *
  * Besides changing the abstract representation (when issueAMod allows it), it
  * issues a eSetInvCost ThermalUnitBlockMod (when issuePMod allows it) so that
  * Solvers consuming the structural data (the DP Solvers) refresh their copy
  * of the design cost via get_design_cost().
  *
  * @param issuePMod controls how physical Modification are issued.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_investment( ModParam issuePMod , ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the coefficients of the Objective
 /** This method updates the coefficients of the Objective.
  *
  * @param subset A set of time instants at which the coefficients must be
  *        updated.
  *
  * @param issuePMod controls how physical Modification are issued.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective( const Subset & subset , ModParam issuePMod ,
                        c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the coefficients of the Objective
 /** This method updates the coefficients of the Objective.
  *
  * @param rng The Range of time instants at which the coefficients must be
  *        updated.
  *
  * @param issuePMod controls how physical Modification are issued.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective( Range rng , ModParam issuePMod ,
                        c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this ThermalUnitBlock is consistent
 /** This function checks whether the data in this ThermalUnitBlock is
  * consistent. The data is consistent if all the following conditions are met.
  *
  * - The minimum power is not greater than the maximum power.
  *
  * - The availability is between 0 and 1.
  *
  * - The delta ramp-up and ramp-down are nonnegative.
  *
  * - The quadratic term of the objective function is nonnegative.
  *
  * - The minimum power, FixedConsumption and InitialPower are nonnegative.
  *
  * - InvestmentCost is 0 unless InitUpDownTime is negative.
  *
  * - If the unit is on before the horizon, InitialPower is not below
  *   MinPower[ 0 ].
  *
  * - StartUpLimit and ShutDownLimit lie between the operational minimum and
  *   maximum power at each instant.
  *
  * If any of the above conditions is not met, std::logic_error is
  * thrown. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );

 void handle_objective_change( FunctionMod * mod , ChnlName chnl );

 /// returns the index, among the active Variable of the Objective, of the
 /// first perspective-cut variable (see the layout in generate_objective())

 Index cut_section_start( void ) const;

/*--------------------------------------------------------------------------*/
 /// the shift that the shut-down section adds to the later Objective sections
 /** The shut-down variables sit right after the start-up ones when the unit
  * pays anything to shut down [see generate_objective()], so every section
  * after them starts this much further on. */

 Index shut_down_offset( void ) const {
  return( f_shut_down_in_obj ? f_time_horizon - init_t : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// the position in the Objective of the first primary (secondary) reserve
 /** The reserve sections follow the start-up, shut-down, active power,
  * commitment and schedule-deviation ones, the primary before the
  * secondary, each only if generate_objective() has put it there. */

 Index reserve_section_start( bool secondary ) const {
  return( 3 * f_time_horizon - init_t + shut_down_offset() +
          ( v_RefSchedule.empty() ? 0 : f_time_horizon ) +
          ( ( secondary && f_primary_in_obj ) ? f_time_horizon : 0 ) );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the vector of MinPower
 std::vector< double > v_MinPower;

 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the vector of Availability
 std::vector< double > v_Availability;

 /// the vector of PrimaryRho
 /** Participation factor: the maximum fraction of active power that can be
  * used as primary reserve, used as the cap in the PrimaryRho_Const
  * constraints. Kept separate from the primary spinning-reserve cost
  * (v_PrimarySpinningReserveCost), the objective coefficient, which may be
  * changed (e.g. by a Lagrangian price) without altering this cap. */
 std::vector< double > v_PrimaryRho;

 /// the vector of SecondaryRho (participation factor / cap, see v_PrimaryRho)
 std::vector< double > v_SecondaryRho;

 /// primary spinning-reserve cost (objective coefficient on the primary
 /// reserve variables), empty if it is zero at every instant
 std::vector< double > v_PrimarySpinningReserveCost;

 /// secondary spinning-reserve cost (objective coefficient on the secondary
 /// reserve variables), empty if it is zero at every instant
 std::vector< double > v_SecondarySpinningReserveCost;

 /// the vector of RampUp
 std::vector< double > v_DeltaRampUp;

 /// the vector of RampDown
 std::vector< double > v_DeltaRampDown;

 /// the vector of QuadTerm
 std::vector< double > v_QuadTerm;

 /// the vector of LinearTerm
 std::vector< double > v_LinearTerm;

 /// dualized linear cost coefficient on the reactive power variables q[t]
 /// (zero in the unit's own objective; set by a dualizing Solver via
 /// set_reactive_linear_term()). Empty means "all zero".
 std::vector< double > v_ReactiveLinearTerm;

 /// the vector of ConstTerm
 std::vector< double > v_ConstTerm;

 /// the vector of StartUpCost
 std::vector< double > v_StartUpCost;

 /// the vector of ShutDownCost, empty if the unit pays nothing to shut down
 std::vector< double > v_ShutDownCost;

 /// true if the Objective has the shut-down variables
 /** The shut-down variables enter the Objective only if the unit pays
  * anything to shut down when generate_objective() runs; without them, a
  * nonzero shut-down cost set afterwards is refused [see
  * set_shutdown_costs()]. */
 bool f_shut_down_in_obj = false;

 /// true if the Objective has the primary reserve variables
 /** They enter it if, when generate_objective() runs, the Variable exist and
  * either the cost is nonzero at some instant or the Configuration of the
  * Objective asks for them [see generate_objective()]. */
 bool f_primary_in_obj = false;

 /// true if the Objective has the secondary reserve variables
 /** See f_primary_in_obj. */
 bool f_secondary_in_obj = false;

 /// the vector of fixed consumption of generator
 std::vector< double > v_FixedConsumption;

 /// the vector of inertia commitment of generator
 std::vector< double > v_InertiaCommitment;

 /// the vector of start-up limits
 std::vector< double > v_StartUpLimit;

 /// the vector of shut-down limits
 std::vector< double > v_ShutDownLimit;

 /// the vector of max ramps steps (SUSD formulation)
 /// entry t > 0 is the largest j <= T - t with the sum of DeltaRampUp over
 /// t, ..., t + j - 1 at most the largest operational maximum power at
 /// t - 1 or later minus the operational minimum power at t - 1; entry 0
 /// is one more than the last j at which InitialPower plus the ramps up to
 /// j is at most the operational maximum power at j (0 if there is none,
 /// -1 if the unit is off before the horizon)
 std::vector< int > v_MaxRampSteps;

 /// the vector of max ramps steps (SUSD formulation)
 /// as v_MaxRampSteps with DeltaRampDown, the smallest operational minimum
 /// power at t - 1 or later and the operational maximum power at t - 1;
 /// entry 0 is one more than the last j at which InitialPower minus the
 /// ramps down to j is at least the operational minimum power at j
 std::vector< int > v_MaxRampDownSteps;

 /// true once compute_ramp_steps() has filled v_MaxRampSteps (always, the
 /// data do not give it)
 bool f_derived_ramp_up_steps = false;

 /// true once compute_ramp_steps() has filled v_MaxRampDownSteps (always,
 /// the data do not give it)
 bool f_derived_ramp_down_steps = false;

 /// true if StartUpLimit is not in the data, and then it is the
 /// operational minimum power [see set_default_limits()]
 bool f_default_start_up_limit = false;

 /// true if ShutDownLimit is not in the data, and then it is the
 /// operational minimum power [see set_default_limits()]
 bool f_default_shut_down_limit = false;

 struct RowCmp;

 /// what update_rows() collects while build_rows() compares the rows,
 /// nullptr while it generates them [see generating_rows()]
 RowCmp * f_row_cmp = nullptr;

 /// the static groups whose rows build_rows() writes with size_rows() or
 /// push_row() when the Constraint are generated [see update_rows()]
 std::vector< std::vector< FRowConstraint > * > f_row_groups;

 /// the dynamic groups [see put_dyn_row()] that add_dyn_rows() registers
 std::vector< DynRows * > f_dyn_groups;

 /// the vector of MinReactivePower
 std::vector< double > v_MinReactivePower;

 /// the vector of MaxReactivePower
 std::vector< double > v_MaxReactivePower;

 /// the vector of MinReactivePowerOn (commitment coefficient; empty = all 0)
 std::vector< double > v_MinReactivePowerOn;

 /// the vector of MaxReactivePowerOn (commitment coefficient; empty = all 0)
 std::vector< double > v_MaxReactivePowerOn;

 /// the reference Schedule to deviate minimally from if there
 std::vector< double > v_RefSchedule;

 // the vector for separating PC-cuts
 std::vector< double > prevpbar;

 /// the vector of index of the variables \f$p_t^{hk}\f$ of the DP formulation.
 /// In particular, v_P_h_k.first = t, v_P_h_k.second.first = h,
 /// v_P_h_k.second.second = k
 std::vector< std::pair< Index , std::pair< Index , Index > > > v_P_h_k;

 /// the vector of index of the variables \f$z_t^{hk}\f$ of the DP formulation.
 /// In particular, v_Z_h_k.first = t, v_Z_h_k.second.first = h,
 /// v_Z_h_k.second.second = k
 std::vector< std::pair< Index , std::pair< Index , Index > > > v_Z_h_k;

 /// the vector of index of the variables \f$p_t^{h}\f$ of the SU formulation.
 /// In particular, v_P_h.first = t, v_P_h.second = h
 std::vector< std::pair< Index , Index > > v_P_h;

 /// the vector of index of the variables \f$z_t^{h}\f$ of the SU formulation.
 /// In particular, v_Z_h.first = t, v_Z_h.second = h
 std::vector< std::pair< Index , Index > > v_Z_h;

 /// the vector of index of the variables \f$\tilde p_t^{k}\f$ of the SD
 /// formulation. In particular, v_P_k.first = t, v_P_k.second = k
 std::vector< std::pair< Index , Index > > v_P_k;

 /// the vector of index of the variables \f$\tilde z_t^{k}\f$ of the SD
 /// formulation. In particular, v_Z_k.first = t, v_Z_k.second = k
 std::vector< std::pair< Index , Index > > v_Z_k;

 /// the vector of index of the ON nodes in the state-space graph
 std::vector< Index > v_nodes_plus;

  /// the vector of index of the OFF nodes in the state-space graph
 std::vector< Index > v_nodes_minus;

 /// the vector of index of the variables \f$y_+^{hk}\f$ of the DP, pt, SU,
 /// SD and SUSD formulations. In particular, v_Y_plus.first = h,
 /// v_Y_plus.second = k
 std::vector< std::pair< Index , Index > > v_Y_plus;

 /// the vector of index of the variables \f$y_-^{hk}\f$ of the DP, pt, SU,
 /// SD and SUSD formulations. In particular, v_Y_minus.first = h,
 /// v_Y_minus.second = k
 std::vector< std::pair< Index , Index > > v_Y_minus;


 /// the investment cost
 double f_InvestmentCost{};

 /// last design cost for which a eSetInvCost Modification was issued
 /** Cache of the design-variable Objective coefficient for which the last
  * eSetInvCost ThermalUnitBlockMod was issued by
  * update_objective_investment().
  * A dualizing Solver rewrites the whole coefficient vector (design included)
  * at each of its iterations; issuing a eSetInvCost every time would reset the
  * cached state of the (DP) Solvers and, more importantly, force the dualizing
  * Bundle to invalidate the linearizations of this component at every
  * iteration (which prevents convergence). Hence eSetInvCost is issued only
  * when the design cost actually changed. Initialised in
  * generate_objective(). */
 double f_last_design_cost = std::numeric_limits< double >::quiet_NaN();

 /// the installable capacity by the user
 double f_Capacity{};

 // total MVA base of this machine
 double f_MBase{};

 /// the InitialPower value
 double f_InitialPower{};

 /// the InitUpDownTime value
 int f_InitUpDownTime{};

 /// the MinUpTime value
 Index f_MinUpTime = 1;

 /// the MinDownTime value
 Index f_MinDownTime = 1;

 /// variable denoting the time-steps unit is subjected to initial conditions
 Index init_t{};

 /// the scale factor
 double f_scale = 1;

 /// "FixToMaximum": if positive, p_t >= the operational maximum power
 int f_fixToMax = 0;

 // this variable indicates which netCDF variables must be ignored
 inline static bool f_ignore_netcdf_vars;

/*-------------------------------- variables -------------------------------*/

 /// the design binary variable
 ColVariable design;

 /// the start-up binary variables
 std::vector< ColVariable > v_start_up;

 /// the shut-down binary variables
 std::vector< ColVariable > v_shut_down;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

 /// the commitment binary variables, in all formulations
 std::vector< ColVariable > v_commitment;

 /// the y^+ arc variables of the pt, DP, SU, SD and SUSD formulations
 std::vector< ColVariable > v_commitment_plus;

 /// the y^- arc variables of the pt, DP, SU, SD and SUSD formulations
 std::vector< ColVariable > v_commitment_minus;

 /// the active power variables, in all formulations
 std::vector< ColVariable > v_active_power;

 /// the reactive power variables
 std::vector< ColVariable > v_reactive_power;

 /// the active power variables for DP model
 std::vector< ColVariable > v_active_power_h_k;

 /// the active power variables for SU model
 std::vector< ColVariable > v_active_power_h;

 /// the active power variables for SD model
 std::vector< ColVariable > v_active_power_k;

 /// the perspective cuts variables, in all formulations with PCuts
 std::vector< ColVariable > v_cut;

 /// the perspective cuts variables for DP model
 std::vector< ColVariable > v_cut_h_k;

 /// the perspective cuts variables for SU model
 std::vector< ColVariable > v_cut_h;

 /// the perspective cuts variables for SD model
 std::vector< ColVariable > v_cut_k;

 /// the perspective cuts variables for SUSD model
 std::vector< ColVariable > v_cut_teta;

 /// the variables for deviation to reference schedule
 std::vector< ColVariable > v_abs_ref_schedule;

/*------------------------------- constraints ------------------------------*/

 /// the reference schedule constraints
 std::vector< FRowConstraint > Reference_Schedule_Const;

 /// the rows p_t >= the operational maximum power (FixToMaximum)
 std::vector< FRowConstraint > fixed_to_max_Power_Const;

 /// the commitment design constraints
 std::vector< FRowConstraint > CommitmentDesign_Const;

 /// the connection min up and down time constraints
 std::vector< FRowConstraint > StartUp_ShutDown_Variables_Const;

 /// the turn on min up and down time constraints
 std::vector< FRowConstraint > StartUp_Const;

 /// the shut-down min up and down time constraints
 std::vector< FRowConstraint > ShutDown_Const;

 /// the RampUp time constraints
 std::vector< FRowConstraint > RampUp_Const;

 /// the RampDown time constraints
 std::vector< FRowConstraint > RampDown_Const;

 /// the minimum power rows (19)
 std::vector< FRowConstraint > MinPower_Const;

 /// the maximum power rows (20)-(27), the dynamic (25), (26) apart
 std::vector< FRowConstraint > MaxPower_Const;

 /// the PrimaryRho fraction constraints
 std::vector< FRowConstraint > PrimaryRho_Const;

 /// the SecondaryRho fraction constraints
 std::vector< FRowConstraint > SecondaryRho_Const;

 /// the spinning-reserve band (capacity + ramp-deliverability) constraints
 /** The same in all formulations: rows (29)-(32) of
  * generate_abstract_constraints(), built only when the unit offers a
  * reserve. */
 std::vector< FRowConstraint > Reserve_Const;

 /* Constraints connecting power variables of 3bin, T and pt
  * formulations with those of DP, SU, SD and SUSD formulations */
 std::vector< FRowConstraint > Eq_ActivePower_Const;

 /** Constraints connecting commitment variables of 3bin and T
  * formulations with those of pt, DP, SU, SD and SUSD formulations */
 std::vector< FRowConstraint > Eq_Commitment_Const;

 /** Constraints connecting start-up variables of 3bin and T
  * formulations with those of pt, DP, SU, SD and SUSD formulations */
 std::vector< FRowConstraint > Eq_StartUp_Const;

 /** Constraints connecting shut-down variables of 3bin and T
  * formulations with those of pt, DP, SU, SD and SUSD formulations */
 std::vector< FRowConstraint > Eq_ShutDown_Const;

 /// the network constraints of the pt, DP, SU, SD and SUSD formulations
 std::vector< FRowConstraint > Network_Const;

 /// the initial perspective cuts constraints
 std::vector< FRowConstraint > Init_PC_Const;

 /** Constraints connecting perspective cuts variables of 3bin, T and pt
  * formulations with those of DP, SU, SD and SUSD formulations */
 std::vector< FRowConstraint > Eq_PC_Const;

 /** Constraints connecting variables of the SUSD formulations with the
  * maximum of the perspective function of the SU and the SD formulations */
 std::vector< FRowConstraint > Max_SUSD_PC_Const;

 /// the perspective dynamic cuts constraints
 std::list< FRowConstraint > PC_cuts;

 /// the row w_0 (or y^{00}) <= [ InitialPower <= ShutDownLimit[ 0 ] ] of
 /// the shut-down at 0 [see build_rows()]
 std::vector< FRowConstraint > ShutDownZero_Const;

 /// the bound constraints 5 of the maximum power of the T formulation
 DynRows MaxPower5_Const;

 /// the bound constraints 6 of the maximum power of the T formulation
 DynRows MaxPower6_Const;

 /// the ramp-up rows of k steps of the SUSD formulation
 DynRows RampUpSUSD_Const;

 /// the ramp-down rows of k steps of the SUSD formulation
 DynRows RampDownSUSD_Const;

 /// the commitment bound constraints
 std::vector< ZOConstraint > Commitment_bound_Const;

 /// the start-up binary bound constraints
 std::vector< ZOConstraint > StartUp_Binary_bound_Const;

 /// the shut-down binary bound constraints
 std::vector< ZOConstraint > ShutDown_Binary_bound_Const;

 /// the commitment fixed to one BoxConstraints
 std::vector< BoxConstraint > Commitment_fixed_to_One_Const;

 /// the reactive power bound constraints (plain box; no commitment gating)
 std::vector< BoxConstraint > ReactivePower_Bound_Const;

 /// the active power bound constraints, 0 <= p[ t ] <= the operational
 /// maximum power: redundant with the rows of the commitment, but readable
 /// by a Solver that only reads the boxes (say, a BoxSolver)
 std::vector< BoxConstraint > ActivePower_Bound_Const;

 /// the reactive power upper bound constraints (commitment-gated variant):
 /// q[t] - Qmax_on[t] u[t] <= Qmax_off[t]
 std::vector< FRowConstraint > ReactivePowerMax_Const;

 /// the reactive power lower bound constraints (commitment-gated variant):
 /// q[t] - Qmin_on[t] u[t] >= Qmin_off[t]
 std::vector< FRowConstraint > ReactivePowerMin_Const;

 //!! Q <= P
 //!! std::vector< FRowConstraint > Reactive_2_Active_Const;

 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 static void static_initialization( void )
 {
  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_availability" ,
   & ThermalUnitBlock::set_availability );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_availability" ,
   & ThermalUnitBlock::set_availability );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_maximum_power" ,
   & ThermalUnitBlock::set_maximum_power );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_maximum_power" ,
   & ThermalUnitBlock::set_maximum_power );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_startup_costs" ,
   & ThermalUnitBlock::set_startup_costs );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_startup_costs" ,
   & ThermalUnitBlock::set_startup_costs );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_shutdown_costs" ,
   & ThermalUnitBlock::set_shutdown_costs );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_shutdown_costs" ,
   & ThermalUnitBlock::set_shutdown_costs );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_const_term" ,
   & ThermalUnitBlock::set_const_term );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_const_term" ,
   & ThermalUnitBlock::set_const_term );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_linear_term" ,
   & ThermalUnitBlock::set_linear_term );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_linear_term" ,
   & ThermalUnitBlock::set_linear_term );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_quad_term" ,
   & ThermalUnitBlock::set_quad_term );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_quad_term" ,
   & ThermalUnitBlock::set_quad_term );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_primary_spinning_reserve_cost" ,
   & ThermalUnitBlock::set_primary_spinning_reserve_cost );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_primary_spinning_reserve_cost" ,
   & ThermalUnitBlock::set_primary_spinning_reserve_cost );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_secondary_spinning_reserve_cost" ,
   & ThermalUnitBlock::set_secondary_spinning_reserve_cost );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_secondary_spinning_reserve_cost" ,
   & ThermalUnitBlock::set_secondary_spinning_reserve_cost );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_initial_power" ,
   & ThermalUnitBlock::set_initial_power );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_initial_power" ,
   & ThermalUnitBlock::set_initial_power );

  register_method< ThermalUnitBlock , MF_int_it , Subset && , bool >(
   "ThermalUnitBlock::set_init_updown_time" ,
   & ThermalUnitBlock::set_init_updown_time );

  register_method< ThermalUnitBlock , MF_int_it , Range >(
   "ThermalUnitBlock::set_init_updown_time" ,
   & ThermalUnitBlock::set_init_updown_time );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::scale" , & ThermalUnitBlock::scale );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::scale" , & ThermalUnitBlock::scale );
 }

/*--------------------------------------------------------------------------*/

 };  // end( class( ThermalUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a ThermalUnitBlock
class ThermalUnitBlockMod : public UnitBlockMod
{
 public:

 /// public enum for the types of ThermalUnitBlockMod
 enum TUB_mod_type
 {
  eSetMaxP = eUBModLastParam , ///< set max power values
  eSetInitP ,                  ///< set initial power values
  eSetInitUD ,                 ///< set initial up/down times
  eSetAv ,                     ///< set availability
  eSetSUC ,                    ///< set start-up costs
  eSetSDC ,                    ///< set shut-down costs
  eSetLinT ,                   ///< set linear term
  eSetQuadT ,                  ///< set quad term
  eSetConstT ,                 ///< set constant term
  eSetPrSpResCost ,            ///< set primary spinning reserve costs
  eSetSecSpResCost ,           ///< set secondary spinning reserve costs
  eSetReactiveLinT ,           ///< set reactive power linear term
  eSetInvCost ,                ///< set design (investment) cost
  eFixVars ,                   ///< the fixed status of some Variable changed
  eTUBModLastParam   ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of ThermalUnitBlockMod. */
  };

 /// constructor, takes the ThermalUnitBlock and the type
 ThermalUnitBlockMod( ThermalUnitBlock * const fblock , const int type )
  : UnitBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the ThermalUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "ThermalUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetMaxP ):
    output << "set max power values";
    break;
   case( eSetInitP ):
    output << "set initial power values";
    break;
   case( eSetInitUD ):
    output << "Set initial up/down times";
    break;
   case( eSetAv ):
    output << "Set availability";
    break;
   case( eSetSUC ):
    output << "Set start-up costs";
    break;
   case( eSetSDC ):
    output << "Set shut-down costs";
    break;
   case( eSetLinT ):
    output << "Set linear term";
    break;
   case( eSetQuadT ):
    output << "Set quad term";
    break;
   case( eSetConstT ):
    output << "Set constant term";
    break;
   case( eSetPrSpResCost ):
    output << "Set primary spinning reserve costs";
    break;
   case( eSetSecSpResCost ):
    output << "Set secondary spinning reserve costs";
    break;
   case( eSetReactiveLinT ):
    output << "Set reactive power linear term";
    break;
   case( eSetInvCost ):
    output << "Set design (investment) cost";
    break;
   case( eFixVars ):
    output << "Changed fixed status of some Variable";
    break;
   default:;
   }
  }
 };  // end( class( ThermalUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS ThermalUnitBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ThermalUnitBlockMod for "ranged" modifications
class ThermalUnitBlockRngdMod : public ThermalUnitBlockMod
{
 public:

 /// constructor: takes the ThermalUnitBlock, the type, and the range
 ThermalUnitBlockRngdMod( ThermalUnitBlock * const fblock ,
                          const int type ,
                          Block::Range rng )
  : ThermalUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the ThermalUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  ThermalUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

 };  // end( class( ThermalUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ThermalUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ThermalUnitBlockMod for "subset" modifications
class ThermalUnitBlockSbstMod : public ThermalUnitBlockMod
{

 public:

 /// constructor: takes the ThermalUnitBlock, the type, and the subset
 ThermalUnitBlockSbstMod( ThermalUnitBlock * const fblock ,
                          const int type ,
                          Block::Subset && nms )
  : ThermalUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the ThermalUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  ThermalUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

 };  // end( class( ThermalUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS ThermalUnitBlockSolution ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [UnitBlock]Solution of a ThermalUnitBlock
/** The ThermalUnitBlockSolution class derives from UnitBlockSolution and
 * adds to the "standard" information stored in there (active power, possibly
 * commitment and primary/secondary reserve) the other information that is
 * typical of the ThermalUnitBlock, i.e.,
 *
 * - if defined, the value of the Thermal Design Variable */

class ThermalUnitBlockSolution : public UnitBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*----------------------------- CONSTANTS ----------------------------------*/

 static constexpr double dNaN = std::numeric_limits< double >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

/*------------------------------- FRIENDS ----------------------------------*/

 friend ThermalUnitBlock;  ///< make ThermalUnitBlock friend

/*--------- CONSTRUCTING AND DESTRUCTING ThermalUnitBlockSolution ----------*/

 /// constructor, it has nothing to do
 explicit ThermalUnitBlockSolution( void ) : f_design( dNaN ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~ThermalUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*----- METHODS DESCRIBING THE BEHAVIOR OF A ThermalUnitBlockSolution -----*/

 void read( const Block * block ) override;

 void write( Block * block ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a ThermalUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a ThermalUnitBlockSolution into a netCDF::NcGroup.
  * The format is the one of UnitBlockSolution
  * [cf. UnitBlockSolution::serialize()], plus:
  *
  * - The scalar variable "ThermalDesign", of type netCDF::NcDouble,
  *   that represent the value of the design variable; the variable is
  *   optional in that the thermal unit may not have any design
  *   variable. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ThermalUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

/*----------- METHODS FOR READING AND WRITING THE SOLUTION -----------------*/

 /// returns the value of the dimensioning variable saved here

 [[nodiscard]] double get_design( void ) const { return( f_design ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the value of the dimensioning variable saved here
 /** Sets the value of the dimensioning variable saved in this
  * ThermalUnitBlockSolution, which is what a Solver filling the Solution out
  * of its own data structures uses [see set_active_power() and the like in
  * UnitBlockSolution]. */

 void set_design( double design ) { f_design = design; }

/*--------------------------------------------------------------------------*/
 /// the start-up indicators, one per time instant, empty if not saved
 /** The Objective pays the start-up through its own Variable, so a Solution
  * that does not carry it can only have it derived from the commitment when
  * it is written back. That is right for one schedule and wrong for a convex
  * combination of several: the start-ups of an averaged commitment are fewer
  * than the average of the start-ups, hence the combination would come out
  * cheaper than it is. Whoever fills a Solution that may be combined has to
  * set these. */

 [[nodiscard]] const std::vector< double > & get_start_up( void ) const {
  return( v_start_up );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the shut-down indicators, one per time instant, empty if not saved

 [[nodiscard]] const std::vector< double > & get_shut_down( void ) const {
  return( v_shut_down );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// set the start-up indicators [see get_start_up()]

 void set_start_up( std::vector< double > && su ) {
  v_start_up = std::move( su );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// set the shut-down indicators [see get_start_up()]

 void set_shut_down( std::vector< double > && sd ) {
  v_shut_down = std::move( sd );
  }

 ThermalUnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "ThermalUnitBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 double f_design;    ///< the value of the dimensioning variable

 std::vector< double > v_start_up;   ///< the start-up indicators, if saved

 std::vector< double > v_shut_down;  ///< the shut-down indicators, if saved

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( ThermalUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __ThermalUnitBlock */

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
