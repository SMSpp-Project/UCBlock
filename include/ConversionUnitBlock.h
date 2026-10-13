/*--------------------------------------------------------------------------*/
/*----------------------- File ConversionUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class ConversionUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define a unit with one
 * commitment and several coupled generators, each of which may sit on its
 * own node of the UCBlock: a plant that converts energy among several
 * carriers (a combined heat and power plant, an electrolyser with its heat,
 * a heat pump, a power-to-gas unit), whose inputs are generators of
 * nonpositive power and whose outputs are generators of nonnegative power,
 * coupled by an operating region given in H-form by the data. The file also
 * defines ConversionUnitBlockMod, the physical Modification of the class,
 * and ConversionUnitBlockSolution, its Solution.
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ConversionUnitBlock
 #define __ConversionUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <limits>
#include <map>

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
/*---------------------- CLASS ConversionUnitBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a multi-output conversion unit
/** The ConversionUnitBlock class derives from UnitBlock and implements a unit
 * of a Unit Commitment problem that is switched on and off as a whole and
 * has \f$ G \geq 1 \f$ generators \f$ \mathcal{G} = \{ 0 , \ldots , G - 1
 * \} \f$, each of which is placed by the enclosing UCBlock on its own node
 * through "GeneratorNode", as the generators of any other unit. The
 * generators are the flows of energy that cross the boundary of the plant,
 * whatever the carrier: an output (electricity, heat, gas, hydrogen) is a
 * generator whose active power \f$ p_{t,g} \f$ is nonnegative, an input
 * (the fuel drawn from a gas node, the electricity drawn by an electrolyser
 * or a heat pump) is a generator whose active power is nonpositive, the
 * power being, as for every UnitBlock, positive when it is injected into the
 * node. The name follows that of the other units, which say what the unit
 * does rather than which carrier it serves: the unit converts the energy of
 * its inputs into that of its outputs (the "converter" of the energy-hub
 * literature, the "process" of the energy-system models, the "multi-link" of
 * the network ones), and nothing in it depends on the carriers, which are
 * only the nodes its generators sit on.
 *
 * The operation of the unit over the instants \f$ \mathcal{T} = \{ 0 ,
 * \ldots , T - 1 \} \f$ is described by a binary commitment \f$ u_t \f$, by
 * binary start-up and shut-down indicators \f$ v_t \f$ and \f$ w_t \f$, by
 * the active powers \f$ p_{t,g} \f$ and, if the UCBlock requires them and the
 * data allow them, by the primary and secondary spinning reserves
 * \f$ p^{pr}_{t,g} \geq 0 \f$ and \f$ p^{sc}_{t,g} \geq 0 \f$ of some of
 * the generators. When the unit is on at the instant \f$ t \f$ its powers
 * lie in the bounded polytope
 * \f[
 *   \mathcal{P}_t = \bigl\{ \, p \in \mathbb{R}^G \,:\, P^{mn}_{t,g} \leq
 *   p_g \leq P^{mx}_{t,g} \; g \in \mathcal{G} \; , \;\;
 *   \underline{b}_{t,m} \leq \sum_{ g \in \mathcal{G} } a_{t,m,g} p_g \leq
 *   \bar{b}_{t,m} \; m \in \mathcal{M} \, \bigr\} \; ,
 * \f]
 * the operating region of the unit, given by the bounds of each generator
 * and by \f$ M \f$ operating rows \f$ \mathcal{M} = \{ 0 , \ldots , M - 1
 * \} \f$ whose coefficients \f$ a_{t,m,g} \f$ may change with the instant
 * (say, the coefficient of performance of a heat pump, which follows the
 * temperature of its source); when it is off all its powers are 0. The
 * unit is thus described by the perspective of \f$ \mathcal{P}_t \f$ (the rows
 * (4), (5) and (8) of generate_abstract_constraints()), which for one instant
 * is the tightest possible description of the set \f$ \{ ( 0 , 0 ) \} \cup (
 * \mathcal{P}_t \times \{ 1 \} ) \f$ of the pairs \f$ ( p , u ) \f$, whatever
 * the shape of the polytope; the perspective formulation of a disjunction of
 * polytopes is that of
 *
 *  E. Balas "Disjunctive Programming: Properties of the Convex Hull of
 *  Feasible Points" Discrete Applied Mathematics 89(1-3), 3 - 44, 1998
 *
 * (see also J.P. Vielma "Mixed Integer Linear Programming Formulation
 * Techniques" SIAM Review 57(1), 3 - 57, 2015). A row may be an equality
 * (\f$ \underline{b}_{t,m} = \bar{b}_{t,m} \f$), so that the region covers,
 * e.g., a back-pressure combined heat and power plant, whose heat is a fixed
 * multiple of its electricity, the polygon of an extraction-condensing one,
 * bounded by the fuel it burns and by its back-pressure line,
 *
 *  M. Zugno, J.M. Morales, H. Madsen "Commitment and Dispatch of Heat and
 *  Power Units via Affinely Adjustable Robust Optimization" Computers &
 *  Operations Research 75, 191 - 201, 2016
 *
 *  M. Koller, R. Hofmann "Mixed-Integer Linear Programming Formulation of
 *  Combined Heat and Power Units for the Unit Commitment Problem" Journal of
 *  Sustainable Development of Energy, Water and Environment Systems 6(4),
 *  755 - 769, 2018
 *
 * the coupling matrix of a converter with several inputs and outputs,
 *
 *  M. Geidl, G. Andersson "Optimal Power Flow of Multiple Energy Carriers"
 *  IEEE Transactions on Power Systems 22(1), 145 - 155, 2007
 *
 * and a fuel consumption that is a convex piecewise-linear function of the
 * outputs (one row per piece, the fuel being an input generator). A region
 * given by its vertices, as in
 *
 *  R. Lahdelma, H. Hakonen "An Efficient Linear Programming Algorithm for
 *  Combined Heat and Power Production" European Journal of Operational
 *  Research 148(1), 141 - 151, 2003
 *
 * is converted to its H-form beforehand; a region that is not convex (a
 * plant with disjoint operating zones, or a combined cycle with its
 * configurations) needs one binary variable per zone and is not a
 * ConversionUnitBlock.
 *
 * The commitment rows are those of the 3bin formulation of ThermalUnitBlock:
 * the logical row (1) of
 *
 *  L.L. Garver "Power Generation Scheduling by Integer Programming -
 *  Development of Theory" Transactions of the AIEE, Part III 81(3), 730 -
 *  734, 1962
 *
 * and the minimum up and down time rows (2), (3) of
 *
 *  D. Rajan, S. Takriti "Minimum Up/Down Polytopes of the Unit Commitment
 *  Problem with Start-Up Costs" IBM Research Report RC23628, 2005
 *
 * which describe the convex hull of the feasible commitments. The rows of
 * the powers that involve the start-up and shut-down indicators, i.e., the
 * limits at the start-up and at the shut-down and the ramps, are written on
 * the model of the T formulation of ThermalUnitBlock, generalised to powers
 * of either sign and to ramps on any combination of the generators; the
 * spinning reserves are offered along directions given by the data, so that
 * the deployment of the reserve of one generator moves the others as the
 * region requires (see generate_abstract_constraints()).
 *
 * <b>Tightness.</b> Without ramps and without start-up and shut-down limits
 * (the rows (6), (7), (9) and (13) of generate_abstract_constraints()), and
 * with or without the reserves, the rows of the unit describe the convex
 * hull of its feasible points with integer commitment, for any \f$ T \f$,
 * any number of generators and any bounded region. The argument extends
 * the one of
 *
 *  B. Hua, R. Baldick "A Convex Primal Formulation for Convex Hull
 *  Pricing" IEEE Transactions on Power Systems 32(5), 3814 - 3823, 2017
 *
 * for one generator (their Theorem 2, stated without initial state): the
 * commitment polytope given by (1) to (3), with the fixings of the initial
 * state and the windows cut at the instant 0, is integral (in the
 * cumulative sums of \f$ v \f$ and \f$ w \f$ every row is a difference
 * constraint, whose matrix is totally unimodular), and a point
 * \f$ ( u , v , w , p ) \f$ of the relaxation is the convex combination
 * \f$ \sum_k \lambda_k ( u^k , v^k , w^k ) \f$ of integral commitments,
 * with which the powers split as \f$ p^k_t = u^k_t \, p_t / u_t \f$ if
 * \f$ u_t > 0 \f$ and \f$ p^k_t = 0 \f$ otherwise: every remaining row is
 * homogeneous of degree one in \f$ ( p_t , r_t , u_t ) \f$ and involves a
 * single instant, hence each piece is feasible, and the pieces average to
 * the point. Thus, with linear costs the continuous relaxation of the
 * problem of the unit alone, with any prices on its generators, has the
 * value of the problem with integer commitment, which is the case of the
 * subproblem of the unit in a Lagrangian decomposition of the UCBlock; with
 * a quadratic cost the feasible set is still the hull, while the value of
 * the relaxation is not that of the integer problem, unless the cost is
 * written in perspective form. The ramps and the start-up and shut-down
 * limits couple an instant to the next one, or the powers to the
 * indicators, so that the pieces of the split are no longer feasible: the
 * formulation is then tight and compact, in the sense of
 *
 *  G. Morales-Espana, J.M. Latorre, A. Ramos "Tight and Compact MILP
 *  Formulation for the Thermal Unit Commitment Problem" IEEE Transactions
 *  on Power Systems 28(4), 4897 - 4908, 2013
 *
 * but not ideal, and the convex hull of a unit with ramps on several
 * outputs is not known.
 *
 * <b>The unit in the UCBlock.</b> Each generator \f$ g \f$ enters the
 * linking rows of the UCBlock as the generator of any other unit (see
 * \ref ucblock_model): its active power in the node injection rows of its
 * node, its spinning reserves (if any) in the reserve requirements of the
 * zone of its node, its inertia per unit of power \f$ h^p_{t,g} \f$ in the
 * inertia requirement of that zone, and its emission rate (the
 * "PollutantRho" of the UCBlock, which is per generator) in the pollutant
 * budgets of that zone. The commitment, which is one for the whole unit, is
 * exposed through one generator only, the commitment generator
 * \f$ g^u \f$ ("CommitmentGenerator", 0 by default): get_commitment()
 * returns \f$ u \f$ for \f$ g = g^u \f$ and nullptr for every other
 * generator, which the UCBlock then treats as a generator without
 * commitment, i.e., one that has no fixed consumption and gives no inertia
 * per unit of commitment. This is because a linking row of the UCBlock that
 * holds two generators of the unit, say an inertia requirement of a zone
 * that contains the nodes of both, would otherwise have the same Variable
 * twice, which a LinearFunction does not allow; the fixed consumption
 * \f$ P^{au}_t \f$ of the unit when it is off and its inertia per unit of
 * commitment \f$ h^u_t \f$ are therefore data of the commitment generator,
 * and are counted at its node and in the zones of its node. The scale
 * factor \f$ \sigma \f$ of the unit (see UnitBlock::scale()) multiplies
 * every term of the unit, as for any UnitBlock.
 *
 * <b>Modules.</b> A unit may be made of \f$ N \geq 1 \f$ identical modules
 * ("Modules"), committed independently of one another: \f$ u_t \f$,
 * \f$ v_t \f$ and \f$ w_t \f$ are then the numbers of modules on, started
 * up and shut down, integer in \f$ \{ 0 , \ldots , N \} \f$, and every
 * datum (bound, row side, limit, ramp, reserve fraction, cost) is that of
 * one module, as in the clustered formulation of
 *
 *  B.S. Palmintier, M.D. Webster "Heterogeneous Unit Clustering for
 *  Efficient Operational Flexibility Modeling" IEEE Transactions on Power
 *  Systems 29(3), 1089 - 1098, 2014
 *
 * Every row of the unit but the minimum down time is homogeneous of degree
 * one in \f$ ( p , r , u , v , w ) \f$, hence it is the sum of the rows of
 * the modules, each with its share of the powers; the right-hand side of
 * (3) becomes \f$ N - u_t \f$, the number of modules that are off, and
 * \f$ u_{-1} \f$ is the number of modules on before the horizon
 * ("InitialCommitment", all of them by default), the others being off
 * since long enough to start, with the shut-downs fixed to 0 before
 * \f$ t_0 \f$ since the modules on before the horizon cannot shut down
 * then. With these rows the integer \f$ ( u , v , w ) \f$ are exactly the
 * sums of \f$ N \f$ feasible commitments of one module each (start up the
 * modules that have been off for the longest time and shut down those that
 * have been on for the longest time), and with equal shares of the powers
 * the region, the bounds and the reserves of the unit are exactly those of
 * the modules; the ramps and the start-up and shut-down limits of the unit
 * are the sums of those of the modules, a relaxation of the modules ramped
 * one by one, as in the clustered formulations. The cumulative
 * sums of \f$ v \f$ and \f$ w \f$ still make every commitment row a
 * difference constraint, so that what is said below about the convex hull
 * holds for any \f$ N \f$. With \f$ N = 1 \f$ the model is the one of a
 * single unit, row by row. A fleet of modules all on or all off at the
 * same time is instead given by "Scale" (see UnitBlock::scale()), and the
 * two are not combined; a unit with several modules has no fixed
 * consumption, which the UCBlock counts as \f$ P^{au}_t ( 1 - u_t ) \f$.
 *
 * <b>Investment.</b> As in the other units, a design variable \f$ x \f$
 * ("x_conversion") exists if the cost \f$ c^{inv} \f$ of one unit of it
 * ("InvestmentCost") is nonzero, it costs \f$ \sigma c^{inv} x \f$ and its
 * domain is given by \f$ \underline{X} \f$ ("MinCapacityDesign", 0 by
 * default) and \f$ \bar{X} \f$ ("MaxCapacityDesign", \f$ -1 \f$ by
 * default), with the sign convention of IntermittentUnitBlock; the unit
 * then has no fixed consumption. Three cases follow, and a unit on before
 * the horizon is built: \f$ x = 1 \f$ in (a), \f$ x \geq u_{-1} \f$ (the
 * modules on before the horizon) in (b), both written in the bounds of
 * \f$ x \f$, and in (c) the rows below.
 *
 * (a) \f$ \bar{X} = -1 \f$: \f$ x \f$ is binary, the unit is built or not
 * (all its \f$ N \f$ modules together), and the rows \f$ u_t \leq N x \f$
 * ("CommitmentDesign_Const_Conversion") are added, as in
 * ThermalUnitBlock.
 *
 * (b) \f$ \bar{X} = -N < -1 \f$: \f$ x \in \{ \underline{X} , \ldots , N
 * \} \f$ is the number of modules built, each committed by itself as above
 * (the unit has \f$ N \f$ modules, "Modules" being either absent or
 * \f$ N \f$), and the rows are \f$ u_t \leq x \f$ and the minimum down time
 * rows (3) with \f$ x - u_t \f$, the number of modules built and off, in
 * place of \f$ N - u_t \f$. Every other row is homogeneous in
 * \f$ ( p , r , u , v , w ) \f$ and does not involve \f$ x \f$, which
 * enters the powers only through \f$ u_t \leq x \f$; as for \f$ N \f$
 * modules (see above, with \f$ x \f$ in place of \f$ N \f$), the integer
 * points are the sums of \f$ x \f$ feasible schedules of one module, and
 * the model is the one of \f$ N \f$ copies of the unit with a binary
 * design each, in fewer rows and Variable and without their symmetry.
 *
 * (c) \f$ \bar{X} \geq 0 \f$: \f$ x \in [ \underline{X} , \bar{X} ] \f$
 * is continuous (say, the nominal capacity of the input, as an extendable
 * unit of an energy-system model), every bound and side of the operating
 * rows is per unit of design, and the commitment is binary (one module).
 * The product \f$ x u_t \f$ then appears in every row that is
 * proportional to the commitment, and it is replaced by the big-M
 * relaxation on the largest capacity \f$ \bar{X} \f$: a row
 * \f$ e_t \leq P u_t \f$ ((5), (8), (11), (12)) becomes
 * \f[
 *   e_t \leq P x + ( - P )^+ \bar{X} ( 1 - u_t ) \; , \qquad
 *   e_t \leq P^+ \bar{X} u_t \; ,
 * \f]
 * and a row \f$ e_t \geq P u_t \f$ ((4), (8), (11), (12)) becomes
 * \f[
 *   e_t \geq P x - P^+ \bar{X} ( 1 - u_t ) \; , \qquad
 *   e_t \geq - ( - P )^+ \bar{X} u_t \; ,
 * \f]
 * e.g., \f$ p_{t,g} \geq P^{mn}_{t,g} x - P^{mn}_{t,g} \bar{X} ( 1 - u_t )
 * \f$ for a minimum power \f$ P^{mn}_{t,g} \geq 0 \f$ (the second row of
 * each pair is written only for the bounds of the powers, the other rows
 * being implied by them when the unit is off). Each pair is exact for
 * \f$ u_t \in \{ 0 , 1 \} \f$, but its continuous relaxation is weak:
 * with \f$ u_t \f$ fractional the minimum power is lost almost entirely,
 * and the bound of the relaxation can be far below the optimum. The
 * start-up and shut-down limits and the ramps, which multiply the
 * indicators, are written in the same way, as in the extendable units
 * with commitment of the energy-system models: for a limit \f$ L \f$ per
 * unit of design of the power \f$ e_t \f$ of a generator (or of a deployed
 * power, (11)) at a start-up, \f$ e_t \leq L x + M ( 1 - v_t ) \f$ and
 * \f$ e_t \geq L x - M ( 1 - v_t ) \f$ for a lower one, with \f$ w_{t+1}
 * \f$ for a shut-down; for the ramps (9) and (13),
 * \f{align*}{
 *   y_{t,k} - y_{t-1,k} & \leq \Delta^+_{t,k} x + M ( 1 - u_t + v_t ) \; ,
 *   & y_{t-1,k} - y_{t,k} & \leq \Delta^-_{t,k} x + M ( 1 - u_t + v_t )
 *   \; , \\
 *   y_{t,k} - y_{t-1,k} & \leq Y^{su,hi}_{t,k} x + M ( 1 - v_t ) \; ,
 *   & y_{t-1,k} - y_{t,k} & \leq Y^{sd,hi}_{t,k} x + M ( 1 - w_t ) \; ,
 *   \\
 *   y_{t,k} - y_{t-1,k} & \leq - Y^{sd,lo}_{t,k} x + M ( 1 - w_t ) \; ,
 *   & y_{t-1,k} - y_{t,k} & \leq - Y^{su,lo}_{t,k} x + M ( 1 - v_t ) \; ,
 * \f}
 * the ranges \f$ Y \f$ being per unit of design; each \f$ M \f$ is
 * \f$ \bar{X} \f$ times the largest value the row may need, i.e., the
 * largest change of \f$ y_{t,k} \f$ over the bounds of two instants (and
 * 0), or the largest bound per unit of the power, plus the negative part of
 * the coefficient of \f$ x \f$. These rows are exact for an integer
 * commitment and weaker than (6), (7) and (9) in the continuous relaxation.
 * A unit on before the horizon has a capacity that covers its InitialPower,
 * \f$ P^{mx}_{0,g} x \geq p_{-1,g} \f$ for \f$ p_{-1,g} > 0 \f$ and
 * \f$ P^{mn}_{0,g} x \leq p_{-1,g} \f$ for \f$ p_{-1,g} < 0 \f$, and its
 * shut-down limits of the instant 0 with \f$ w_0 \f$ in the same big-M form
 * ("InitialDesign_Const_Conversion"); the minimum power and the operating
 * rows are not imposed on InitialPower, which is the operation before the
 * horizon, with the capacity that it had then. Several modules are
 * refused with a continuous design. For investment studies the
 * recommended forms remain (b), exact and linear,
 * and the investment from outside, through "Scale" changed by an
 * InvestmentBlock (see UnitBlock::scale()), which multiplies the data and
 * involves no product; (c) is meant for translating models that describe
 * an extendable unit with commitment in this way.
 *
 * The continuous relaxation of a unit of \f$ N \f$ modules and that of one
 * module with \f$ \sigma = N \f$ coincide: every row is homogeneous in
 * \f$ ( p , r , u , v , w ) \f$, and \f$ \sigma u_t \f$ ranges in
 * \f$ [ 0 , N ] \f$ as the number of modules on does, the data being those
 * of one module in both. An InvestmentBlock, whose linearization is the one
 * of the relaxation of the units, therefore sees the same function of an
 * integer number of modules whether it changes "Scale" of one module (with
 * integer investment Variable) or the number of modules, and the former is
 * the form for investment studies through it; the two differ only with an
 * integer commitment, where the modules of a scaled unit are all on or
 * all off, and the exact form with modules committed one by one is (b) in
 * the UCBlock.
 *
 * "Scale" \f$ \neq 1 \f$ is refused together with \f$ | \bar{X} | > 1 \f$
 * and with \f$ N > 1 \f$, since they all replicate the unit.
 *
 * <b>Delays.</b> A generator \f$ g \f$ may deliver its power with a delay
 * of \f$ d_g \geq 0 \f$ instants ("Delay"), as the heat of a district
 * heating pipe or the gas of a long pipeline reaches its node some time
 * after the plant has produced it: the operation of the unit at the instant
 * \f$ t \f$ (its commitment and the point of its region) gives the power
 * \f$ q_{t,g} \f$ of the generator, which the node receives at the instant
 * \f$ t + d_g \f$, i.e., the power at the node at the instant \f$ s \f$ is
 * \f$ p_{s,g} = q_{s - d_g , g} \f$. Every row of the unit is written in
 * the instants of operation, hence a delay is a renaming of the Variable
 * and changes neither the form of the rows nor what is said below about
 * their tightness: an output delayed by \f$ d \f$ with respect to its
 * input is a generator with delay \f$ d \f$ and an input with none, and
 * outputs with different delays (say, two pipes of different lengths) have
 * different delays. Two cases differ at the ends of the horizon. If the
 * delay is not cyclic, the node receives at the instants
 * \f$ s < \min\{ d_g , T \} \f$ the powers \f$ \hat{p}_{s,g} \f$ of the
 * operation before the horizon ("InitialDelayedPower", 0 by default), which
 * are data, and the powers \f$ q_{t,g} \f$ of the instants
 * \f$ t \geq T - d_g \f$ reach the node after the horizon, so that they are
 * Variable of the unit that enter no linking row of the UCBlock; if it is
 * cyclic ("CyclicDelay"), the horizon wraps around, \f$ p_{s,g} =
 * q_{( s - d_g ) \bmod T , g} \f$, as for a horizon that repeats itself.
 * A generator with a delay offers no spinning reserve (the reserve would be
 * deployed at an instant other than the one of the operation). The
 * commitment, the fixed consumption and the inertia per unit of
 * commitment are of the instants of operation, and the inertia per unit of
 * power and the emissions of a generator are of the power its node
 * receives.
 *
 * As for ThermalUnitBlock, all the data are given per instant and no length
 * of the time step appears in any row: a ramp is the largest change of a
 * combination of the powers between two consecutive instants, a minimum up
 * or down time is a number of instants, and a cost is the cost of one
 * instant. The reactive power is not represented, and the unit has no
 * reactive power Variable even in a UCBlock with an AC network. */

class ConversionUnitBlock : public UnitBlock
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father block
 /** Constructor of ConversionUnitBlock, taking possibly a pointer of its
  * father Block. */

 explicit ConversionUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of ConversionUnitBlock

 virtual ~ConversionUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the ConversionUnitBlock. Besides the mandatory "type" attribute of any
  * :Block, the group must contain all the data required by the base
  * UnitBlock, as described in the comments to UnitBlock::deserialize(
  * netCDF::NcGroup ), to which we refer for the dimensions "TimeHorizon"
  * (\f$ T \f$), "NumberIntervals" and "ChangeIntervals".
  *
  * A datum "indexed over the instants" is a netCDF::NcDouble variable whose
  * first dimension is either of size 1, and then it has the same value at
  * every instant, or "NumberIntervals" (or "TimeHorizon" if
  * "NumberIntervals" is not given), and then its entry \f$ i \f$ is the
  * value at every instant \f$ t \f$ with ChangeIntervals[ \f$ i - 1 \f$ ]
  * \f$ < t \leq \f$ ChangeIntervals[ \f$ i \f$ ], as in ThermalUnitBlock;
  * a datum indexed over the instants and over the generators has
  * "NumberGenerators" (or a dimension of size 1, for the same value for all
  * the generators) as its second dimension, and similarly for the other
  * second dimensions below. Each datum is expanded to one value per instant
  * (and per generator, row or ramp) when it is read. The group contains:
  *
  * - The dimension "NumberGenerators", the number \f$ G \geq 1 \f$ of
  *   generators; optional with default 1.
  *
  * - The scalar variable "CommitmentGenerator", of type netCDF::NcUint,
  *   the generator \f$ g^u \in \mathcal{G} \f$ through which the commitment
  *   is exposed to the UCBlock (see the class documentation); optional with
  *   default 0.
  *
  * - The variables "MinPower" and "MaxPower", indexed over the instants
  *   and over the generators, the bounds \f$ P^{mn}_{t,g} \leq P^{mx}_{t,g}
  *   \f$ of the active power of the generator \f$ g \f$ when the unit is on
  *   at the instant \f$ t \f$; both may be negative (an input) and they
  *   must be finite. "MaxPower" is mandatory, "MinPower" is optional with
  *   default 0. An outage of the unit at an instant is described by bounds
  *   that leave only \f$ p_{t,g} = 0 \f$, with which the unit may still be
  *   on, paying \f$ c_t \f$ and giving its inertia.
  *
  * - The dimension "NumberOperatingRows", the number \f$ M \geq 0 \f$ of
  *   operating rows; optional with default 0, in which case the region is
  *   the box of the bounds.
  *
  * - The variable "OperatingMatrix", the coefficients \f$ a_{t,m,g} \f$ of
  *   the operating rows; mandatory if \f$ M > 0 \f$. It is either indexed
  *   over "NumberOperatingRows" and "NumberGenerators", and then the
  *   coefficients are the same at every instant, or it has a first
  *   dimension indexed over the instants as the other data (of size 1,
  *   "NumberIntervals" or "TimeHorizon", with the meaning given above)
  *   followed by "NumberOperatingRows" and "NumberGenerators". Whether a
  *   coefficient is zero may change with the instant, since the rows of
  *   different instants are distinct.
  *
  * - The variables "OperatingLHS" and "OperatingRHS", indexed over the
  *   instants and over "NumberOperatingRows", the sides
  *   \f$ \underline{b}_{t,m} \leq \bar{b}_{t,m} \f$ of the operating rows,
  *   \f$ \underline{b}_{t,m} u_t \leq \sum_g a_{t,m,g} p_{t,g} \leq
  *   \bar{b}_{t,m} u_t \f$; each is optional, with default \f$ - \infty \f$,
  *   respectively \f$ + \infty \f$, but each row must have a finite side at
  *   some instant. Whether a side is finite is a property of the row and
  *   of the instant that cannot change after the generation of the
  *   abstract representation.
  *
  * - The dimension "NumberRampRows", the number \f$ K \f$ of ramp
  *   combinations \f$ \mathcal{K} = \{ 0 , \ldots , K - 1 \} \f$, and the
  *   variable "RampMatrix", indexed over "NumberRampRows" and
  *   "NumberGenerators", their coefficients \f$ c_{k,g} \f$: the ramp
  *   \f$ k \f$ limits the change of \f$ \sum_g c_{k,g} p_{t,g} \f$ between
  *   two consecutive instants (e.g., of the fuel burnt by a plant with
  *   several outputs). Both are optional: without them \f$ K = G \f$ and
  *   \f$ c_{k,g} = 1 \f$ if \f$ k = g \f$ and 0 otherwise, i.e., the ramps
  *   are on each generator.
  *
  * - The variables "DeltaRampUp" and "DeltaRampDown", indexed over the
  *   instants and over "NumberRampRows" (over "NumberGenerators" if
  *   "RampMatrix" is not given), the ramp-up limit \f$ \Delta^+_{t,k} \geq
  *   0 \f$ and the ramp-down limit \f$ \Delta^-_{t,k} \geq 0 \f$ of the
  *   combination \f$ k \f$ from \f$ t - 1 \f$ to \f$ t \f$ (from
  *   InitialPower if \f$ t = 0 \f$). Each is optional: without it there are
  *   no ramp rows of that kind.
  *
  * - The variables "StartUpLimit" and "StartUpLowerLimit", indexed over the
  *   instants and over the generators, the largest and the smallest active
  *   power \f$ \bar{S}^{su}_{t,g} \f$ and \f$ \underline{S}^{su}_{t,g} \f$
  *   of the generator at the instant \f$ t \f$ when the unit starts up at
  *   \f$ t \f$, and the variables "ShutDownLimit" and "ShutDownLowerLimit",
  *   the largest and the smallest active power \f$ \bar{S}^{sd}_{t,g} \f$
  *   and \f$ \underline{S}^{sd}_{t,g} \f$ of the generator at the instant
  *   \f$ t - 1 \f$ when the unit shuts down at \f$ t \f$. All are optional,
  *   with default the bounds of the power at the instant they refer to
  *   (\f$ P^{mx}_{t,g} \f$ and \f$ P^{mn}_{t,g} \f$ for the start-up,
  *   \f$ P^{mx}_{t-1,g} \f$ and \f$ P^{mn}_{t-1,g} \f$ for the shut-down,
  *   those of the instant 0 if \f$ t = 0 \f$), i.e., no limit; a lower
  *   limit must not exceed the upper one. Unlike in ThermalUnitBlock, the
  *   default is no limit rather than the minimum power, since the minimum
  *   power of an input is its largest consumption. Whether a generator has
  *   limits (an upper or a lower one) is fixed at the generation of the
  *   abstract representation.
  *
  * - The variables "PrimaryRho" and "SecondaryRho", indexed over the
  *   instants and over the generators, the largest fractions
  *   \f$ \rho^{pr}_{t,g} \geq 0 \f$ and \f$ \rho^{sc}_{t,g} \geq 0 \f$ of
  *   the absolute value of the active power of the generator that it can
  *   offer as primary and secondary spinning reserve. Each is optional; a
  *   generator whose fraction is zero at every instant offers no reserve of
  *   that kind and has no Variable for it. A generator that offers a reserve
  *   must have a power of the same sign at every instant, i.e.,
  *   \f$ P^{mn}_{t,g} \geq 0 \f$ for all \f$ t \f$ (we then write
  *   \f$ \varsigma_g = 1 \f$) or \f$ P^{mx}_{t,g} \leq 0 \f$ for all
  *   \f$ t \f$ (\f$ \varsigma_g = -1 \f$), the fraction being of
  *   \f$ \varsigma_g p_{t,g} \f$.
  *
  * - The variable "ReserveDirection", indexed over "NumberGenerators" twice,
  *   whose entry \f$ D_{h,g} \f$ is the change of the power of the
  *   generator \f$ h \f$ per unit of reserve of the generator \f$ g \f$
  *   deployed, the column of a generator \f$ g \f$ that offers a reserve
  *   having \f$ D_{g,g} = 1 \f$: e.g., for a back-pressure plant whose heat
  *   \f$ h \f$ follows its electricity \f$ e \f$ with ratio \f$ \eta \f$,
  *   \f$ D_{e,e} = 1 \f$ and \f$ D_{h,e} = \eta \f$, so that a deployed
  *   reserve keeps the plant on its back-pressure line. Optional, with
  *   default the identity (the reserve of a generator moves that generator
  *   only, which an equality row of the region then forbids).
  *
  * - The variables "PrimarySpinningReserveCost" and
  *   "SecondarySpinningReserveCost", indexed over the instants and over the
  *   generators, the costs \f$ c^{pr}_{t,g} \f$ and \f$ c^{sc}_{t,g} \f$ of
  *   one unit of reserve; optional with default 0.
  *
  * - The variables "QuadTerm" and "LinearTerm", indexed over the instants
  *   and over the generators, the coefficients \f$ a_{t,g} \geq 0 \f$ and
  *   \f$ b_{t,g} \f$ of the cost \f$ a_{t,g} p_{t,g}^2 + b_{t,g} p_{t,g} \f$
  *   of the generator (the price of a fuel can be either the linear cost of
  *   the fuel generator, or the price at its node); both optional with
  *   default 0.
  *
  * - The variables "ConstTerm", "StartUpCost" and "ShutDownCost", indexed
  *   over the instants, the cost \f$ c_t \f$ of the unit being on at
  *   \f$ t \f$ and the costs \f$ c^{su}_t \f$ and \f$ c^{sd}_t \f$ of a
  *   start-up and of a shut-down at \f$ t \f$; all optional with default 0.
  *
  * - The variable "InitialPower", indexed over "NumberGenerators", the
  *   active powers \f$ p_{-1,g} \f$ of the generators at the instant
  *   \f$ -1 \f$, before the horizon; optional with default 0, and used only
  *   if the unit is on before the horizon (InitUpDownTime \f$ > 0 \f$), in
  *   which case it must lie in the region \f$ \mathcal{P}_0 \f$ of the
  *   instant 0 (an exception is thrown otherwise).
  *
  * - The variable "Delay", of type netCDF::NcUint and indexed over
  *   "NumberGenerators" (or a scalar, for the same delay for all the
  *   generators), the delays \f$ d_g \f$ in instants (see the class
  *   documentation); optional with default 0.
  *
  * - The variable "CyclicDelay", of type netCDF::NcInt and indexed over
  *   "NumberGenerators" (or a scalar), nonzero for the generators whose
  *   delay wraps around the horizon; optional with default 0.
  *
  * - The variable "InitialDelayedPower", indexed over a dimension of size
  *   1 (the same value at every instant) or at least
  *   \f$ H = \max_g \min\{ d_g , T \} \f$ over the generators with a delay
  *   that is not cyclic (say, "NumberDelayedInstants") and over
  *   "NumberGenerators", whose entry \f$ [ s ][ g ] \f$ is the power
  *   \f$ \hat{p}_{s,g} \f$ the generator delivers at the instant
  *   \f$ s < \min\{ d_g , T \} \f$ as the result of the operation before
  *   the horizon (the other entries are ignored); optional with default 0.
  *   It is independent of InitialPower, which is the operating point at
  *   \f$ -1 \f$ for the rows that tie the instant 0 to the past.
  *
  * - The scalar variable "InitUpDownTime", of type netCDF::NcInt, the
  *   number \f$ \tau_0 \f$ of instants the unit has been on
  *   (\f$ \tau_0 > 0 \f$) or off (\f$ \tau_0 \leq 0 \f$, for \f$ - \tau_0
  *   \f$ instants) before the instant 0, with the meaning and the default
  *   of ThermalUnitBlock (\f$ - \tau^- \f$ if InitialPower is 0 for all the
  *   generators and \f$ \tau^+ \f$ otherwise).
  *
  * - The scalar variables "MinUpTime" and "MinDownTime", of type
  *   netCDF::NcUint, the minimum up time \f$ \tau^+ \f$ and the minimum down
  *   time \f$ \tau^- \f$ of the unit, with the meaning, the default and the
  *   clamping of ThermalUnitBlock.
  *
  * - The variable "FixedConsumption", indexed over the instants, the power
  *   \f$ P^{au}_t \geq 0 \f$ that the unit draws at the node of the
  *   commitment generator when it is off; optional with default 0.
  *
  * - The variable "InertiaCommitment", indexed over the instants, the
  *   inertia \f$ h^u_t \f$ that the unit gives when it is on, counted in the
  *   zone of the node of the commitment generator; optional with default 0.
  *
  * - The variable "InertiaPower", indexed over the instants and over the
  *   generators, the inertia \f$ h^p_{t,g} \f$ per unit of active power of
  *   each generator; optional with default 0.
  *
  * - The scalar variable "Scale", the scale factor \f$ \sigma \f$ of the
  *   unit (see UnitBlock::scale()); optional with default 1.
  *
  * - The scalar variable "Modules", of type netCDF::NcUint, the number
  *   \f$ N \geq 1 \f$ of identical modules of the unit, every other datum
  *   being of one module (see the class documentation); optional with
  *   default 1. It is refused together with Scale \f$ \neq 1 \f$ (both
  *   replicate the unit) and with FixedConsumption.
  *
  * - The scalar variables "InvestmentCost", "MinCapacityDesign" and
  *   "MaxCapacityDesign", of type netCDF::NcDouble, the cost
  *   \f$ c^{inv} \f$ of one unit of design and the bounds
  *   \f$ \underline{X} \geq 0 \f$ and \f$ \bar{X} \f$ of the design (see
  *   the class documentation); InvestmentCost is optional with default 0
  *   (no design), the other two are read only if it is given, with
  *   defaults 0 and \f$ -1 \f$ (binary design). A negative
  *   \f$ \bar{X} \f$ must be an integer, \f$ \underline{X} \leq
  *   | \bar{X} | \f$, and \f$ \bar{X} < -1 \f$ fixes the number of
  *   modules; with \f$ \bar{X} \geq 0 \f$ the data are per unit of design
  *   and several modules are refused. A unit with a design has no
  *   FixedConsumption; if it is on before the horizon it is built (see
  *   the class documentation), and with \f$ \bar{X} \geq 0 \f$ its
  *   InitialPower is not checked against the region, which depends on
  *   \f$ x \f$: the design is only required to cover it (see the class
  *   documentation).
  *
  * - The scalar variable "InitialCommitment", of type netCDF::NcUint, the
  *   number \f$ u_{-1} \in \{ 1 , \ldots , N \} \f$ of modules on before
  *   the horizon if \f$ \tau_0 > 0 \f$, InitialPower being their total
  *   power, which must lie in \f$ u_{-1} \f$ times the region of the
  *   instant 0; optional, with default \f$ N \f$, and not used if
  *   \f$ \tau_0 \leq 0 \f$.
  *
  * The data are checked when they are read: an exception is thrown if a
  * bound is not finite, a lower bound or limit exceeds the upper one, an
  * operating row has no finite side, a reserve is offered by a generator
  * whose power changes sign or with a direction whose diagonal entry is not
  * 1, or InitialPower is outside the region of the instant 0 for a unit on
  * before the horizon. That the region of an instant is not empty is not
  * checked (a Solver finds it out). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /// extends UnitBlock::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends UnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the ConversionUnitBlock
 /** Generates the abstract Variable of the ConversionUnitBlock. Let
  * \f$ \tau_0 \f$, \f$ \tau^+ \f$ and \f$ \tau^- \f$ be InitUpDownTime,
  * MinUpTime and MinDownTime (see deserialize()); as in ThermalUnitBlock,
  * the first instant \f$ t_0 \f$ at which the commitment is free is
  * \f[
  *   t_0 = \min\{ T , \max\{ 0 , \tau^+ - \tau_0 \} \}
  *   \;\text{ if } \tau_0 > 0 \; , \qquad
  *   t_0 = \min\{ T , \max\{ 0 , \tau^- + \tau_0 \} \}
  *   \;\text{ if } \tau_0 \leq 0 \; ,
  * \f]
  * and we write \f$ u_{-1} = 1 \f$ if \f$ \tau_0 > 0 \f$ and
  * \f$ u_{-1} = 0 \f$ otherwise, and \f$ p_{-1,g} \f$ for InitialPower if
  * \f$ \tau_0 > 0 \f$ and 0 otherwise. The Variable are
  *
  * - the design \f$ x \f$ ("x_conversion"), if InvestmentCost is nonzero:
  *   binary, integer or continuous as in the class documentation, its
  *   bounds being the BoxConstraint "DesignBound_Conversion";
  *
  * - the binary commitment \f$ u_t \f$, \f$ t \in \mathcal{T} \f$
  *   ("u_conversion"), the number of modules on, integer in
  *   \f$ \{ 0 , \ldots , N \} \f$, if \f$ N > 1 \f$;
  *
  * - the binary start-up and shut-down indicators \f$ v_t \f$
  *   ("v_conversion") and \f$ w_t \f$ ("w_conversion"),
  *   \f$ t \in \mathcal{T} \f$, \f$ v_t = 1 \f$ if the unit is off at
  *   \f$ t - 1 \f$ and on at \f$ t \f$, \f$ w_t = 1 \f$ if it is on at
  *   \f$ t - 1 \f$ and off at \f$ t \f$ (the numbers of modules started up
  *   and shut down, integer, if \f$ N > 1 \f$);
  *
  * - the active powers \f$ p_{t,g} \f$, \f$ t \in \mathcal{T} \f$,
  *   \f$ g \in \mathcal{G} \f$ ("p_conversion", one vector per generator),
  *   continuous, their bounds being rows; for a generator with a delay
  *   \f$ d_g > 0 \f$ that is not cyclic, the vector has the
  *   \f$ T \f$ powers delivered at the instants of the horizon, the first
  *   \f$ \min\{ d_g , T \} \f$ fixed to InitialDelayedPower by their
  *   bounds, followed by the \f$ \min\{ d_g , T \} \f$ powers
  *   \f$ q_{t,g} \f$ of the instants of operation
  *   \f$ t \geq T - d_g \f$, which are delivered after the horizon (see
  *   get_tail_power());
  *
  * - the primary reserves \f$ p^{pr}_{t,g} \geq 0 \f$ ("pr_conversion") of
  *   the generators with a nonzero PrimaryRho, if the UCBlock has a primary
  *   reserve requirement (bit 0 of the reserve_vars set by the UCBlock), and
  *   the secondary reserves \f$ p^{sc}_{t,g} \geq 0 \f$ ("sc_conversion")
  *   under the same conditions with SecondaryRho and bit 1; a reserve with
  *   no Variable is zero.
  *
  * The Variable do not depend on \p stvv, nor on
  * f_BlockConfig->f_static_variables_Configuration, the unit having one
  * formulation only. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the ConversionUnitBlock
 /** Generates the abstract Constraint of the ConversionUnitBlock, in the
  * notation of deserialize() and generate_abstract_variables(); each group
  * of rows is given with the name under which it is added to the Block, and
  * a group whose data are absent is not generated. A term in \f$ w_{t+1}
  * \f$ is absent for \f$ t = T - 1 \f$. Every row is written in the
  * instants of operation: for a generator with a delay \f$ d_g \f$,
  * \f$ p_{t,g} \f$ below is the power \f$ q_{t,g} \f$ of the operation at
  * \f$ t \f$, delivered at \f$ t + d_g \f$ (see the class documentation),
  * and the two coincide if \f$ d_g = 0 \f$.
  *
  * <b>Commitment.</b> For \f$ t \in \mathcal{T} \f$,
  * \f{align*}{
  *   & u_t - u_{t-1} = v_t - w_t \tag{1} \\
  *   & \sum_{ s = \max\{ 0 , t - \tau^+ + 1 \} }^{ t } v_s \leq u_t
  *     \tag{2} \\
  *   & \sum_{ s = \max\{ 0 , t - \tau^- + 1 \} }^{ t } w_s \leq N - u_t
  *     \tag{3}
  * \f}
  * ("Logical_Const_Conversion", "MinUp_Const_Conversion",
  * "MinDown_Const_Conversion"), with the constant \f$ u_{-1} \f$ in (1) at
  * \f$ t = 0 \f$, and the bounds ("Commitment_Bound_Conversion")
  * \f$ u_{-1} \leq u_t \leq N \f$ for \f$ t < t_0 \f$ if \f$ \tau_0 > 0
  * \f$, \f$ u_t = 0 \f$ for \f$ t < t_0 \f$ if \f$ \tau_0 \leq 0 \f$, and
  * \f$ 0 \leq u_t \leq N \f$ otherwise, \f$ N \f$ being the number of
  * modules (1 by default, see the class documentation); with
  * \f$ N > 1 \f$ and \f$ \tau_0 > 0 \f$ also \f$ w_t \leq 0 \f$ for
  * \f$ t < t_0 \f$ ("ShutDown_Initial_Bound_Conversion" for
  * \f$ t \geq 1 \f$, whose bound is \f$ w_t \leq N \f$ otherwise, and
  * "ShutDownZero_Bound_Conversion" below for \f$ t = 0 \f$). These are the
  * rows of the 3bin formulation of ThermalUnitBlock, with the windows of
  * (2) and (3) cut at the instant 0 rather than at \f$ t_0 \f$, the bounds of
  * \f$ u \f$ saying the rest; they force \f$ v_t = w_t = 0 \f$ for
  * \f$ t < t_0 \f$.
  *
  * <b>Bounds of the powers.</b> For \f$ t \in \mathcal{T} \f$ and
  * \f$ g \in \mathcal{G} \f$ ("MinPower_Const_Conversion",
  * "MaxPower_Const_Conversion"),
  * \f{align*}{
  *   & P^{mn}_{t,g} u_t \leq p_{t,g} \; , \tag{4} \\
  *   & p_{t,g} \leq P^{mx}_{t,g} u_t \; , \tag{5}
  * \f}
  * and \f$ N \min\{ 0 , P^{mn}_{t,g} \} \leq p_{t,g} \leq N \max\{ 0 ,
  * P^{mx}_{t,g} \} \f$ as a BoxConstraint ("ActivePower_Bound_Conversion"),
  * which the rows imply but a Solver that only reads the bounds of the
  * Variable needs; the powers delivered at \f$ s < \min\{ d_g , T \} \f$
  * by the operation before the horizon are fixed to
  * \f$ \hat{p}_{s,g} \f$ in the same BoxConstraint. For a generator with an
  * upper start-up or shut-down limit, let
  * \f[
  *   \hat{U}^{su}_{t,g} = \min\{ \bar{S}^{su}_{t,g} , P^{mx}_{t,g} \} \; ,
  *   \qquad
  *   \hat{U}^{sd}_{t,g} = \min\{ \bar{S}^{sd}_{t+1,g} , P^{mx}_{t,g} \}
  *   \;\text{ for } t < T - 1 \; , \qquad
  *   \hat{U}^{sd}_{T-1,g} = P^{mx}_{T-1,g} \; ;
  * \f]
  * then (5) is replaced, if \f$ \tau^+ \geq 2 \f$, by
  * \f[
  *   p_{t,g} \leq P^{mx}_{t,g} u_t - ( P^{mx}_{t,g} - \hat{U}^{su}_{t,g} )
  *   v_t - ( P^{mx}_{t,g} - \hat{U}^{sd}_{t,g} ) w_{t+1} \; , \tag{6}
  * \f]
  * and, if \f$ \tau^+ = 1 \f$, when a run of one instant is possible, by
  * the two rows
  * \f{align*}{
  *   p_{t,g} & \leq P^{mx}_{t,g} u_t - ( P^{mx}_{t,g} - \hat{U}^{sd}_{t,g} )
  *     w_{t+1} - ( \hat{U}^{sd}_{t,g} - \hat{U}^{su}_{t,g} )^+ v_t \; , \\
  *   p_{t,g} & \leq P^{mx}_{t,g} u_t - ( \hat{U}^{su}_{t,g} -
  *     \hat{U}^{sd}_{t,g} )^+ w_{t+1} - ( P^{mx}_{t,g} - \hat{U}^{su}_{t,g}
  *     ) v_t \; , \tag{6'}
  * \f}
  * the first of which is absent at \f$ t = T - 1 \f$, where it coincides with
  * the second; these are the rows (20) and (21) of ThermalUnitBlock, and
  * they cap the power at \f$ \hat{U}^{su}_{t,g} \f$ at a start-up, at
  * \f$ \hat{U}^{sd}_{t,g} \f$ at the last instant before a shut-down and at
  * the smaller of the two for a run of one instant. Symmetrically, for a
  * generator with a lower limit, with \f$ \hat{L}^{su}_{t,g} = \max\{
  * \underline{S}^{su}_{t,g} , P^{mn}_{t,g} \} \f$ and \f$ \hat{L}^{sd}_{t,g}
  * = \max\{ \underline{S}^{sd}_{t+1,g} , P^{mn}_{t,g} \} \f$ (and
  * \f$ \hat{L}^{sd}_{T-1,g} = P^{mn}_{T-1,g} \f$), (4) is replaced by
  * \f[
  *   p_{t,g} \geq P^{mn}_{t,g} u_t + ( \hat{L}^{su}_{t,g} - P^{mn}_{t,g} )
  *   v_t + ( \hat{L}^{sd}_{t,g} - P^{mn}_{t,g} ) w_{t+1} \tag{7}
  * \f]
  * if \f$ \tau^+ \geq 2 \f$, and by the two rows mirroring (6') if
  * \f$ \tau^+ = 1 \f$.
  *
  * <b>Operating region</b> ("Operating_Const_Conversion"). For
  * \f$ t \in \mathcal{T} \f$ and \f$ m \in \mathcal{M} \f$, the row
  * \f[
  *   \sum_{ g \in \mathcal{G} } a_{t,m,g} p_{t,g} - \bar{b}_{t,m} u_t \leq 0
  *   \; , \qquad
  *   \sum_{ g \in \mathcal{G} } a_{t,m,g} p_{t,g} - \underline{b}_{t,m} u_t
  *   \geq 0 \tag{8}
  * \f]
  * for each finite side, one equality row if the two sides coincide.
  *
  * <b>Ramps</b> ("RampUp_Const_Conversion", "RampDown_Const_Conversion",
  * generated if DeltaRampUp, respectively DeltaRampDown, is given). Let
  * \f$ y_{t,k} = \sum_g c_{k,g} p_{t,g} \f$, with \f$ y_{-1,k} = \sum_g
  * c_{k,g} p_{-1,g} \f$, and let \f$ [ Y^{su,lo}_{t,k} , Y^{su,hi}_{t,k} ]
  * \f$ be the range of \f$ y_{t,k} \f$ over the box of the start-up limits
  * at \f$ t \f$ intersected with the bounds, i.e., \f$ Y^{su,hi}_{t,k} =
  * \sum_g \max\{ c_{k,g} \hat{L}^{su}_{t,g} , c_{k,g} \hat{U}^{su}_{t,g}
  * \} \f$ and \f$ Y^{su,lo}_{t,k} \f$ the same with the minimum, and
  * \f$ [ Y^{sd,lo}_{t,k} , Y^{sd,hi}_{t,k} ] \f$ the range of
  * \f$ y_{t-1,k} \f$ over the box of the shut-down limits at \f$ t \f$
  * (that of the default limits if none is given). Then, for
  * \f$ t \in \mathcal{T} \f$ and \f$ k \in \mathcal{K} \f$,
  * \f{align*}{
  *   y_{t,k} - y_{t-1,k} & \leq \Delta^+_{t,k} ( u_t - v_t ) +
  *     Y^{su,hi}_{t,k} v_t - Y^{sd,lo}_{t,k} w_t \; , \\
  *   y_{t-1,k} - y_{t,k} & \leq \Delta^-_{t,k} ( u_t - v_t ) +
  *     Y^{sd,hi}_{t,k} w_t - Y^{su,lo}_{t,k} v_t \; , \tag{9}
  * \f}
  * where \f$ u_t - v_t = u_{t-1} - w_t \f$ is 1 exactly when the unit is on
  * at \f$ t - 1 \f$ and at \f$ t \f$: the rows bound the change of
  * \f$ y \f$ by the ramps then, keep \f$ y_{t,k} \f$ in the start-up range
  * at a start-up, keep \f$ y_{t-1,k} \f$ in the shut-down range at a
  * shut-down, and are void when the unit is off at both instants. With one
  * generator, the default limits and nonnegative bounds, (9) is the ramp row
  * of the T formulation of ThermalUnitBlock. At \f$ t = 0 \f$ the terms in
  * \f$ y_{-1,k} \f$ are constants; a unit on before the horizon whose
  * InitialPower is outside the shut-down box at 0 (\f$ u_{-1} \f$ times
  * it, with several modules) cannot shut down at 0, and the bound
  * \f$ w_0 \leq 0 \f$ ("ShutDownZero_Bound_Conversion") says so also
  * without ramps (the bound is \f$ 0 \leq w_0 \leq N \f$ otherwise).
  *
  * <b>Spinning reserves.</b> Let \f$ \mathcal{R} \f$ be the generators that
  * have reserve Variable, \f$ r_{t,g} = p^{pr}_{t,g} + p^{sc}_{t,g} \f$ for
  * \f$ g \in \mathcal{R} \f$ (a reserve with no Variable being 0), and
  * \f$ q^s_{t,h} = p_{t,h} + s \sum_{ g \in \mathcal{R} } D_{h,g} r_{t,g}
  * \f$ the power of the generator \f$ h \f$ when the reserves are deployed
  * upwards (\f$ s = 1 \f$) or downwards (\f$ s = -1 \f$), all of them in the
  * same direction, as a frequency deviation deploys them. The rows are
  * ("PrimaryRho_Const_Conversion", "SecondaryRho_Const_Conversion")
  * \f[
  *   p^{pr}_{t,g} \leq \rho^{pr}_{t,g} \varsigma_g p_{t,g} \; , \qquad
  *   p^{sc}_{t,g} \leq \rho^{sc}_{t,g} \varsigma_g p_{t,g} \tag{10}
  * \f]
  * for \f$ g \in \mathcal{R} \f$, which give no reserve when the unit is
  * off; then, for \f$ s \in \{ 1 , -1 \} \f$, the deployed powers satisfy
  * the bounds ("Reserve_Bound_Const_Conversion")
  * \f[
  *   P^{mn}_{t,h} u_t \leq q^s_{t,h} \leq P^{mx}_{t,h} u_t \tag{11}
  * \f]
  * and the operating rows ("Reserve_Operating_Const_Conversion")
  * \f[
  *   \underline{b}_{t,m} u_t \leq \sum_h a_{t,m,h} q^s_{t,h} \leq
  *   \bar{b}_{t,m} u_t \; , \tag{12}
  * \f]
  * where only the sides that the deployment can violate are written: the
  * upper side of (11) for \f$ h \f$ if \f$ s D_{h,g} > 0 \f$ for some
  * \f$ g \in \mathcal{R} \f$, the lower one if \f$ s D_{h,g} < 0 \f$ for
  * some, and the same for (12) with the entries \f$ \sum_h a_{t,m,h} D_{h,g}
  * \f$ in place of \f$ D_{h,g} \f$, the other sides being implied by (4),
  * (5) and (8) since \f$ r \geq 0 \f$. An equality operating row is written
  * as \f$ \sum_{ g \in \mathcal{R} } ( \sum_h a_{t,m,h} D_{h,g} ) r_{t,g} = 0
  * \f$, once for both directions, if any of these entries is nonzero. For a
  * generator \f$ h \f$ with start-up or shut-down limits the bounds of (11)
  * are each replaced by two rows, as the rows (30) of ThermalUnitBlock,
  * \f{align*}{
  *   q^s_{t,h} & \leq P^{mx}_{t,h} u_t - ( P^{mx}_{t,h} - \hat{U}^{su}_{t,h}
  *     ) v_t \; , &
  *   q^s_{t,h} & \leq P^{mx}_{t,h} u_t - ( P^{mx}_{t,h} - \hat{U}^{sd}_{t,h}
  *     ) w_{t+1} \; , \\
  *   q^s_{t,h} & \geq P^{mn}_{t,h} u_t + ( \hat{L}^{su}_{t,h} -
  *     P^{mn}_{t,h} ) v_t \; , &
  *   q^s_{t,h} & \geq P^{mn}_{t,h} u_t + ( \hat{L}^{sd}_{t,h} -
  *     P^{mn}_{t,h} ) w_{t+1} \; ,
  * \f}
  * the second of each pair being absent at \f$ t = T - 1 \f$. Finally, with
  * the ramps, the deployed powers are within the ramps from the instant
  * before ("Reserve_Ramp_Const_Conversion"): with \f$ y^s_{t,k} = \sum_h
  * c_{k,h} q^s_{t,h} \f$,
  * \f{align*}{
  *   y^s_{t,k} - y_{t-1,k} & \leq \Delta^+_{t,k} ( u_t - v_t ) +
  *     Y^{su,hi}_{t,k} v_t - Y^{sd,lo}_{t,k} w_t \; , \\
  *   y_{t-1,k} - y^s_{t,k} & \leq \Delta^-_{t,k} ( u_t - v_t ) +
  *     Y^{sd,hi}_{t,k} w_t - Y^{su,lo}_{t,k} v_t \; , \tag{13}
  * \f}
  * the first if \f$ s \sum_h c_{k,h} D_{h,g} > 0 \f$ for some
  * \f$ g \in \mathcal{R} \f$ and the second if it is \f$ < 0 \f$ for some.
  * With one generator, \f$ D = 1 \f$ and the default limits, (10) and (11)
  * are the rows (28)-(30) of ThermalUnitBlock and (13) its rows (31) and
  * (32), which there are relaxed to \f$ P^{mx}_t \f$ at a start-up while
  * here they are relaxed to the start-up range. Of course, the effect of a
  * deployed reserve on the balance of the other nodes (say, the heat that a
  * back-pressure plant produces with its reserve) is not represented, as
  * the deployment of the reserves is not represented in the UCBlock for any
  * unit.
  *
  * <b>Changes after the generation.</b> The rows are written by a single
  * method, from the data as they are, both here and when a setter changes a
  * datum after the generation: the Variable of every row depend only on the
  * data that cannot change then (the operating and ramp matrices, the
  * reserve directions, the finite sides, the presence of the ramps and of
  * the limits, the minimum up time, the generators that offer a reserve),
  * hence a change enters the coefficients and the sides of the rows, and the
  * abstract representation is the one a fresh load of the changed data
  * would generate. If \p stcc (or, if \p stcc is nullptr,
  * f_BlockConfig->f_static_constraints_Configuration) is a
  * SimpleConfiguration< int > with a nonzero value and \f$ N = 1 \f$, the
  * bounds \f$ 0 \leq v_t \leq 1 \f$ and \f$ 0 \leq w_t \leq 1 \f$ are added
  * as explicit ZOConstraint as well ("StartUp_Binary_Bound_Conversion",
  * "ShutDown_Binary_Bound_Conversion"). */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the ConversionUnitBlock
 /** Generates the Objective of the ConversionUnitBlock, a DQuadFunction to
  * be minimized,
  * \f[
  *   \sigma \sum_{ t \in \mathcal{T} } \Bigl( c^{su}_t v_t + c^{sd}_t w_t
  *   + c_t u_t + \sum_{ g \in \mathcal{G} } \bigl( a_{t,g} p_{t,g}^2 +
  *   b_{t,g} p_{t,g} \bigr) + \sum_{ g \in \mathcal{R} } \bigl(
  *   c^{pr}_{t,g} p^{pr}_{t,g} + c^{sc}_{t,g} p^{sc}_{t,g} \bigr) \Bigr)
  *   \; , \tag{14}
  * \f]
  * where every Variable is in the Objective, also those whose cost is 0, so
  * that every cost can be changed after the generation (a reserve with no
  * Variable having no term). The order of the Variable in the DQuadFunction
  * is \f$ v \f$, \f$ w \f$, \f$ u \f$ (each over \f$ \mathcal{T} \f$), then
  * the powers of each generator in turn (those of the instants of operation,
  * \f$ p_{t,g} = q_{t,g} \f$ in the order of \f$ t \f$, the powers
  * delivered from before the horizon having no term), then the primary and
  * the secondary reserves of each generator that has them, in turn, and
  * last the design \f$ x \f$, with the term \f$ \sigma c^{inv} x \f$, if
  * it exists. The parameter \p objc is not used. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*------------- Methods for checking the ConversionUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** Returns true if every Variable of the ConversionUnitBlock is feasible
  * and the violation of every Constraint, relative or absolute, is not
  * greater than the tolerance, which are read from \p fsbc or from
  * f_BlockConfig->f_is_feasible_Configuration as in
  * ThermalUnitBlock::is_feasible() (Block::DefaultFeasTol and the relative
  * violation by default). Only the abstract representation is used, hence
  * \p useabstract is ignored; with no abstract Variable the method returns
  * true. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE UNIT -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ConversionUnitBlock
 * @{ */

 /// returns the number of generators \f$ G \f$
 Index get_number_generators( void ) const override { return( f_G ); }

 /// returns the number of operating rows \f$ M \f$
 Index get_number_operating_rows( void ) const { return( f_M ); }

 /// returns the number of ramp combinations \f$ K \f$ (0 if no ramp)
 Index get_number_ramp_rows( void ) const {
  return( ( v_DeltaRampUp.empty() && v_DeltaRampDown.empty() ) ? 0 : f_K );
  }

 /// returns the commitment generator \f$ g^u \f$
 Index get_commitment_generator( void ) const {
  return( f_commitment_generator );
  }

 /// returns the minimum power \f$ P^{mn}_{t,g} \f$
 double get_min_power( Index t , Index generator = 0 ) const override {
  return( v_MinPower[ t ][ generator ] );
  }

 /// returns the maximum power \f$ P^{mx}_{t,g} \f$
 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower[ t ][ generator ] );
  }

 /// returns the coefficient \f$ a_{t,m,g} \f$ of the operating row m
 double get_operating_coefficient( Index m , Index generator ,
                                   Index t = 0 ) const {
  return( v_OperatingMatrix[ v_OperatingMatrix.shape()[ 0 ] > 1 ? t : 0 ][
                                                          m ][ generator ] );
  }

 /// true if the coefficients of the operating rows change with the instant
 bool has_varying_operating_matrix( void ) const {
  return( v_OperatingMatrix.shape()[ 0 ] > 1 );
  }

 /// returns the left-hand side \f$ \underline{b}_{t,m} \f$ (may be -INF)
 double get_operating_lhs( Index t , Index m ) const;

 /// returns the right-hand side \f$ \bar{b}_{t,m} \f$ (may be +INF)
 double get_operating_rhs( Index t , Index m ) const;

 /// returns the coefficient \f$ c_{k,g} \f$ of the ramp combination k
 double get_ramp_coefficient( Index k , Index generator ) const;

 /// returns the initial power \f$ p_{-1,g} \f$ (as read)
 double get_initial_power( Index generator ) const {
  return( v_InitialPower[ generator ] );
  }

 /// returns the delay \f$ d_g \f$ of the generator, in instants
 Index get_delay( Index generator ) const {
  return( v_Delay.empty() ? 0 : v_Delay[ generator ] );
  }

 /// true if the delay of the generator wraps around the horizon
 bool is_cyclic_delay( Index generator ) const {
  return( ( ! v_CyclicDelay.empty() ) && v_CyclicDelay[ generator ] &&
          get_delay( generator ) );
  }

 /// returns the power \f$ \hat{p}_{s,g} \f$ delivered at s < \f$ d_g \f$
 /** Returns the power that the generator delivers at the instant \p s of the
  * horizon as the result of the operation before it (InitialDelayedPower),
  * 0 if it is not given. */

 double get_initial_delayed_power( Index s , Index generator ) const;

 /// returns the initial up/down time \f$ \tau_0 \f$
 int get_init_up_down_time( void ) const { return( f_InitUpDownTime ); }

 /// returns the minimum up time \f$ \tau^+ \f$
 Index get_min_up_time( void ) const { return( f_MinUpTime ); }

 /// returns the minimum down time \f$ \tau^- \f$
 Index get_min_down_time( void ) const { return( f_MinDownTime ); }

 /// returns the number \f$ N \f$ of identical modules of the unit
 /** Returns "Modules", or \f$ | \bar{X} | \f$ for an integer design (see
  * the class documentation). */

 Index get_number_modules( void ) const {
  return( design_mode() == 'b' ? Index( - f_MaxCapacityDesign )
                               : f_Modules );
  }

 /// true if the unit has a design variable \f$ x \f$
 bool has_design( void ) const { return( f_has_design ); }

 /// the kind of design: 0 none, 'a' binary, 'b' integer, 'c' continuous
 char design_mode( void ) const {
  if( ! has_design() )
   return( 0 );
  if( f_MaxCapacityDesign >= 0 )
   return( 'c' );
  return( f_MaxCapacityDesign < -1 ? 'b' : 'a' );
  }

 /// returns the investment cost \f$ c^{inv} \f$ of one unit of design
 double get_investment_cost( void ) const { return( f_InvestmentCost ); }

 /// returns the lower bound \f$ \underline{X} \f$ of the design
 double get_min_capacity_design( void ) const {
  return( f_MinCapacityDesign );
  }

 /// returns "MaxCapacityDesign" \f$ \bar{X} \f$, signed as it is read
 double get_max_capacity_design( void ) const {
  return( f_MaxCapacityDesign );
  }

 /// returns \f$ u_{-1} \f$, the number of modules on before the horizon
 Index get_initial_commitment( void ) const;

 /// returns the largest multiplier of the data of one module
 /** Returns the largest number by which the data of the unit are
  * multiplied, which the UCBlock uses in the bounds of the injections at
  * the nodes (see UnitBlock::get_design_ub()): \f$ \bar{X} \f$ with a
  * continuous design, the number \f$ N \f$ of modules otherwise. */

 double get_design_ub( void ) const override {
  return( design_mode() == 'c' ? f_MaxCapacityDesign
                               : double( get_number_modules() ) );
  }

 /// returns the linear cost \f$ b_{t,g} \f$
 double get_linear_term( Index t , Index generator ) const {
  return( v_LinearTerm.empty() ? 0 : v_LinearTerm[ t ][ generator ] );
  }

 /// returns the quadratic cost \f$ a_{t,g} \f$
 double get_quad_term( Index t , Index generator ) const {
  return( v_QuadTerm.empty() ? 0 : v_QuadTerm[ t ][ generator ] );
  }

 /// returns the cost \f$ c_t \f$ of being on
 double get_const_term( Index t ) const {
  return( v_ConstTerm.empty() ? 0 : v_ConstTerm[ t ] );
  }

 /// returns the start-up cost \f$ c^{su}_t \f$
 double get_start_up_cost( Index t ) const {
  return( v_StartUpCost.empty() ? 0 : v_StartUpCost[ t ] );
  }

 /// returns the shut-down cost \f$ c^{sd}_t \f$
 double get_shut_down_cost( Index t ) const {
  return( v_ShutDownCost.empty() ? 0 : v_ShutDownCost[ t ] );
  }

 /// returns the fixed consumption, of the commitment generator only
 const double * get_fixed_consumption( Index generator ) const override {
  if( ( generator != f_commitment_generator ) ||
      v_FixedConsumption.empty() )
   return( nullptr );
  return( v_FixedConsumption.data() );
  }

 /// returns the inertia commitment, of the commitment generator only
 const double * get_inertia_commitment( Index generator ) const override {
  if( ( generator != f_commitment_generator ) ||
      v_InertiaCommitment.empty() )
   return( nullptr );
  return( v_InertiaCommitment.data() );
  }

 /// returns the inertia per unit of power of the given generator
 const double * get_inertia_power( Index generator ) const override {
  if( ( generator >= v_InertiaPower.size() ) ||
      v_InertiaPower[ generator ].empty() )
   return( nullptr );
  return( v_InertiaPower[ generator ].data() );
  }

 /// returns the scale factor \f$ \sigma \f$
 double get_scale( void ) const override { return( f_scale ); }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE Variable OF THE UNIT ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the ConversionUnitBlock
 * @{ */

 /// the unit is committed
 bool has_commitment( void ) const override { return( true ); }

 /// the unit has no reactive power
 bool has_reactive_power( void ) const override { return( false ); }

 /// true if some generator offers primary reserve and the UCBlock asks it
 bool has_primary_reserve( void ) const override {
  return( ( reserve_vars & 1u ) && f_has_primary );
  }

 /// true if some generator offers secondary reserve and the UCBlock asks it
 bool has_secondary_reserve( void ) const override {
  return( ( reserve_vars & 2u ) && f_has_secondary );
  }

 /// returns the commitment for the commitment generator, nullptr otherwise
 /** Returns the array of the commitment Variable \f$ u_t \f$ if
  * \p generator is the commitment generator \f$ g^u \f$, and nullptr for
  * every other generator, as well as before the generation of the
  * Variable (see the class documentation). */

 ColVariable * get_commitment( Index generator ) override {
  if( ( generator != f_commitment_generator ) || v_commitment.empty() )
   return( nullptr );
  return( v_commitment.data() );
  }

 /// returns the commitment Variable \f$ u_t \f$, whatever the generator
 ColVariable * get_unit_commitment( void ) {
  return( v_commitment.empty() ? nullptr : v_commitment.data() );
  }

 /// returns the power Variable of the generator at the operation instant t
 /** Returns the Variable \f$ q_{t,g} \f$ of the power of the generator due
  * to the operation of the unit at the instant \p t, i.e., the one
  * delivered at \f$ t + d_g \f$ (see generate_abstract_variables()). */

 ColVariable & get_operating_power( Index generator , Index t ) {
  return( op( generator , t ) );
  }

 /// the number of powers of the generator delivered after the horizon
 Index get_number_tail_powers( Index generator ) const {
  return( ( generator < v_active_power.size() ) &&
          ( v_active_power[ generator ].size() > f_time_horizon ) ?
          v_active_power[ generator ].size() - f_time_horizon : 0 );
  }

 /// the powers of the generator delivered after the horizon, if any
 /** Returns the array of the get_number_tail_powers() Variable of the powers
  * of the operation instants \f$ t \geq T - d_g \f$, in the order of
  * \f$ t \f$, which are delivered after the horizon (nullptr if there are
  * none). */

 ColVariable * get_tail_power( Index generator ) {
  return( get_number_tail_powers( generator ) ?
          v_active_power[ generator ].data() + f_time_horizon : nullptr );
  }

 /// returns the design Variable \f$ x \f$ ("x_conversion")
 /** Returns the design Variable, which is in the abstract representation
  * only if has_design() (its value is NaN otherwise). */

 ColVariable & get_design( void ) { return( design ); }

 /// returns the design Variable \f$ x \f$, const version
 const ColVariable & get_const_design( void ) const { return( design ); }

 /// returns the start-up Variable \f$ v_t \f$
 ColVariable * get_start_up( void ) {
  return( v_start_up.empty() ? nullptr : v_start_up.data() );
  }

 /// returns the shut-down Variable \f$ w_t \f$
 ColVariable * get_shut_down( void ) {
  return( v_shut_down.empty() ? nullptr : v_shut_down.data() );
  }

 /// returns the active power Variable of the given generator
 /** Returns the array of the powers \f$ p_{s,g} \f$ that the generator
  * delivers at its node at the instants \f$ s \in \mathcal{T} \f$ of the
  * horizon (see generate_abstract_variables()). */

 ColVariable * get_active_power( Index generator ) override {
  if( ( generator >= v_active_power.size() ) ||
      v_active_power[ generator ].empty() )
   return( nullptr );
  return( v_active_power[ generator ].data() );
  }

 /// returns the primary reserve Variable of the generator, if any
 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( ( generator >= v_primary_reserve.size() ) ||
      v_primary_reserve[ generator ].empty() )
   return( nullptr );
  return( v_primary_reserve[ generator ].data() );
  }

 /// returns the secondary reserve Variable of the generator, if any
 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( ( generator >= v_secondary_reserve.size() ) ||
      v_secondary_reserve[ generator ].empty() )
   return( nullptr );
  return( v_secondary_reserve[ generator ].data() );
  }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a ConversionUnitBlockSolution of this ConversionUnitBlock
 /** Returns a ConversionUnitBlockSolution, the parts of which are decided
  * as in UnitBlock::get_Solution(), the start-up and shut-down indicators
  * being saved with the commitment. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

 /// returns an empty ConversionUnitBlockSolution
 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UNIT -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the ConversionUnitBlock
 * @{ */

 /// extends UnitBlock::serialize( netCDF::NcGroup )
 /** Writes the ConversionUnitBlock in the format of deserialize(), every
  * datum indexed over the instants being written over "TimeHorizon". */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR CHANGING THE DATA ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the ConversionUnitBlock
 *
 * Each setter changes the values of a datum at the instants of \p subset
 * (unordered if \p ordered is false) or \p rng, \p values pointing to as
 * many values; the data indexed also over the generators, the operating
 * rows or the ramp combinations take the index of the column in
 * \p generator, \p row or \p ramp. Each of these also has a form with no
 * such index, which is the one registered in the methods factory: its
 * \p subset or \p rng then runs over the entries [ g ][ t ] taken
 * generator after generator (row after row, ramp after ramp), i.e., the
 * index of the entry of \f$ t \f$ and \f$ g \f$ is \f$ g T + t \f$, as for
 * the data of the arcs of a HydroUnitBlock. If nothing changes the method
 * does
 * nothing. Otherwise, as for the other units, the physical representation
 * is changed and a ConversionUnitBlockMod is issued according to
 * \p issuePMod, and the abstract representation, if it is there, is changed
 * according to \p issueAMod, by changing the coefficients of the Objective
 * or those and the sides of the rows (see generate_abstract_constraints()),
 * all in one GroupModification. A change that the abstract representation
 * cannot follow, i.e., one that would change which rows or which Variable
 * there are (a finite side becoming infinite, a ramp or a reserve cost that
 * has no row or no Variable), or one that makes the data inconsistent (a
 * minimum above a maximum, a reserve generator changing sign, an
 * InitialPower outside the region), is refused with an exception and
 * changes nothing.
 * @{ */

 /// sets the linear cost \f$ b_{t,g} \f$ of a generator
 void set_linear_term( MF_dbl_it values , Index generator ,
                       Subset && subset , const bool ordered = false ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

 /// sets the linear cost \f$ b_{t,g} \f$ of a generator, range version
 void set_linear_term( MF_dbl_it values , Index generator ,
                       Range rng = INFRange ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

 /// sets the quadratic cost \f$ a_{t,g} \geq 0 \f$ of a generator
 void set_quad_term( MF_dbl_it values , Index generator ,
                     Subset && subset , const bool ordered = false ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

 /// sets the quadratic cost of a generator, range version
 void set_quad_term( MF_dbl_it values , Index generator ,
                     Range rng = INFRange ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

 /// sets the primary reserve cost \f$ c^{pr}_{t,g} \f$ of a generator
 void set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                         Index generator ,
                                         Range rng = INFRange ,
                                         ModParam issuePMod = eNoBlck ,
                                         ModParam issueAMod = eNoBlck );

 /// sets the secondary reserve cost \f$ c^{sc}_{t,g} \f$ of a generator
 void set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                           Index generator ,
                                           Range rng = INFRange ,
                                           ModParam issuePMod = eNoBlck ,
                                           ModParam issueAMod = eNoBlck );

 /// sets the cost \f$ c_t \f$ of being on
 void set_const_term( MF_dbl_it values , Subset && subset ,
                      const bool ordered = false ,
                      ModParam issuePMod = eNoBlck ,
                      ModParam issueAMod = eNoBlck );

 /// sets the cost \f$ c_t \f$ of being on, range version
 void set_const_term( MF_dbl_it values , Range rng = INFRange ,
                      ModParam issuePMod = eNoBlck ,
                      ModParam issueAMod = eNoBlck );

 /// sets the start-up cost \f$ c^{su}_t \f$
 void set_startup_costs( MF_dbl_it values , Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the start-up cost \f$ c^{su}_t \f$, range version
 void set_startup_costs( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the shut-down cost \f$ c^{sd}_t \f$
 void set_shutdown_costs( MF_dbl_it values , Subset && subset ,
                          const bool ordered = false ,
                          ModParam issuePMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

 /// sets the shut-down cost \f$ c^{sd}_t \f$, range version
 void set_shutdown_costs( MF_dbl_it values , Range rng = INFRange ,
                          ModParam issuePMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

 /// sets the maximum power \f$ P^{mx}_{t,g} \f$ of a generator
 void set_maximum_power( MF_dbl_it values , Index generator ,
                         Subset && subset , const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the maximum power of a generator, range version
 void set_maximum_power( MF_dbl_it values , Index generator ,
                         Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the minimum power \f$ P^{mn}_{t,g} \f$ of a generator
 void set_minimum_power( MF_dbl_it values , Index generator ,
                         Subset && subset , const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the minimum power of a generator, range version
 void set_minimum_power( MF_dbl_it values , Index generator ,
                         Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the right-hand side \f$ \bar{b}_{t,m} \f$ of an operating row
 void set_operating_rhs( MF_dbl_it values , Index row ,
                         Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the left-hand side \f$ \underline{b}_{t,m} \f$ of an operating row
 void set_operating_lhs( MF_dbl_it values , Index row ,
                         Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the ramp-up limit \f$ \Delta^+_{t,k} \f$ of a combination
 void set_delta_ramp_up( MF_dbl_it values , Index ramp ,
                         Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the ramp-down limit \f$ \Delta^-_{t,k} \f$ of a combination
 void set_delta_ramp_down( MF_dbl_it values , Index ramp ,
                           Range rng = INFRange ,
                           ModParam issuePMod = eNoBlck ,
                           ModParam issueAMod = eNoBlck );

 /// sets the initial powers \f$ p_{-1,g} \f$ of the generators in rng
 /** Here \p rng is a range of generators, not of instants. */

 void set_initial_power( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the powers \f$ \hat{p}_{s,g} \f$ delivered from before the horizon
 /** Sets the powers that the generator \p generator delivers at the
  * instants \p rng of the horizon as the result of the operation before it
  * (InitialDelayedPower), the instants at or after \f$ \min\{ d_g , T \}
  * \f$ being ignored; they are the bounds of the delivered powers of those
  * instants. */

 void set_initial_delayed_power( MF_dbl_it values , Index generator ,
                                 Range rng = INFRange ,
                                 ModParam issuePMod = eNoBlck ,
                                 ModParam issueAMod = eNoBlck );

 /// sets the initial up/down time \f$ \tau_0 \f$
 /** Sets \f$ \tau_0 \f$ to *values, \p rng being ignored unless it is
  * empty. The change moves the bounds of \f$ u \f$ and \f$ w_0 \f$ and the
  * constants of the rows at the instant 0; the minimum up and down times
  * are not clamped anew, since their clamping only matters for the
  * Variable, which do not depend on it. */

 void set_init_updown_time( MF_int_it values , Range rng = INFRange ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );


 /// the forms of the setters over the entries [ g ][ t ] (see above)

 void set_linear_term( MF_dbl_it values , Subset && subset ,
                       const bool ordered = false ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );
 void set_linear_term( MF_dbl_it values , Range rng = INFRange ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

 void set_quad_term( MF_dbl_it values , Subset && subset ,
                     const bool ordered = false ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );
 void set_quad_term( MF_dbl_it values , Range rng = INFRange ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

 void set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                         Subset && subset ,
                                         const bool ordered = false ,
                                         ModParam issuePMod = eNoBlck ,
                                         ModParam issueAMod = eNoBlck );
 void set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                         Range rng = INFRange ,
                                         ModParam issuePMod = eNoBlck ,
                                         ModParam issueAMod = eNoBlck );

 void set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                           Subset && subset ,
                                           const bool ordered = false ,
                                           ModParam issuePMod = eNoBlck ,
                                           ModParam issueAMod = eNoBlck );
 void set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                           Range rng = INFRange ,
                                           ModParam issuePMod = eNoBlck ,
                                           ModParam issueAMod = eNoBlck );

 void set_maximum_power( MF_dbl_it values , Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );
 void set_maximum_power( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_minimum_power( MF_dbl_it values , Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );
 void set_minimum_power( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_operating_rhs( MF_dbl_it values , Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );
 void set_operating_rhs( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_operating_lhs( MF_dbl_it values , Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );
 void set_operating_lhs( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_delta_ramp_up( MF_dbl_it values , Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );
 void set_delta_ramp_up( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_delta_ramp_down( MF_dbl_it values , Subset && subset ,
                           const bool ordered = false ,
                           ModParam issuePMod = eNoBlck ,
                           ModParam issueAMod = eNoBlck );
 void set_delta_ramp_down( MF_dbl_it values , Range rng = INFRange ,
                           ModParam issuePMod = eNoBlck ,
                           ModParam issueAMod = eNoBlck );

 /// sets the initial powers \f$ p_{-1,g} \f$ of the generators in subset
 void set_initial_power( MF_dbl_it values , Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 /// sets the powers delivered from before the horizon, entries [ g ][ s ]
 void set_initial_delayed_power( MF_dbl_it values , Subset && subset ,
                                 const bool ordered = false ,
                                 ModParam issuePMod = eNoBlck ,
                                 ModParam issueAMod = eNoBlck );
 void set_initial_delayed_power( MF_dbl_it values , Range rng = INFRange ,
                                 ModParam issuePMod = eNoBlck ,
                                 ModParam issueAMod = eNoBlck );

 /// sets \f$ \tau_0 \f$ to *values if \p subset is not empty
 void set_init_updown_time( MF_int_it values , Subset && subset ,
                            const bool ordered = false ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

 /// sets the investment cost \f$ c^{inv} \f$ to *values
 /** Sets \f$ c^{inv} \f$ to *values, \p rng being ignored unless it is
  * empty; refused if the unit has no design variable. A zero cost keeps
  * the design variable, which then costs nothing (a unit with a zero
  * InvestmentCost written and read back has none). */

 void set_investment_cost( MF_dbl_it values , Range rng = INFRange ,
                           ModParam issuePMod = eNoBlck ,
                           ModParam issueAMod = eNoBlck );

 /// sets \f$ c^{inv} \f$ to *values if \p subset is not empty
 void set_investment_cost( MF_dbl_it values , Subset && subset ,
                           const bool ordered = false ,
                           ModParam issuePMod = eNoBlck ,
                           ModParam issueAMod = eNoBlck );

 /// sets the scale factor \f$ \sigma \f$ [see UnitBlock::scale()]
 void scale( MF_dbl_it values , Subset && subset ,
             const bool ordered = false ,
             c_ModParam issuePMod = eNoBlck ,
             c_ModParam issueAMod = eNoBlck ) override;

 using UnitBlock::scale;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 using MAdbl = boost::multi_array< double , 2 >;  ///< a matrix of data

 struct RowCmp;  ///< what update_rows() collects [see the .cpp]

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// reads "OperatingMatrix", static or indexed over the instants
 void read_operating_matrix( const netCDF::NcGroup & group ,
                             const std::string & who );

 /// reads the delays and the powers delivered from before the horizon
 void read_delays( const netCDF::NcGroup & group , const std::string & who );

 /// the number of instants at which g delivers the operation before 0
 Index head_size( Index g ) const;

 /// the position in v_active_power[ g ] of the operation instant t
 Index position_of( Index g , Index t ) const;

 /// the operation instant of the position i >= head_size( g )
 Index instant_of( Index g , Index i ) const;

 /// the power Variable of g at the operation instant t
 ColVariable & op( Index g , Index t ) {
  return( v_active_power[ g ][ position_of( g , t ) ] );
  }

 /// the first instant at which the commitment is free
 Index first_free_instant( void ) const;

 /// checks the consistency of the data, throws if they are not
 void check_data( const std::string & who ) const;

 /// computes the sign \f$ \varsigma_g \f$ of the reserve generators
 void compute_reserve_signs( const std::string & who );

 /// the box of the start-up (su = true) or shut-down limits of g at t
 void limit_box( bool su , Index t , Index g , double & lo ,
                 double & hi ) const;

 /// the range of the ramp combination k over the box of the limits
 void combination_range( bool su , Index t , Index k , double & lo ,
                         double & hi ) const;

 /// whether the generator g has an upper (lower) start-up/shut-down limit
 bool has_upper_limit( Index g ) const;
 bool has_lower_limit( Index g ) const;

 /// writes (or compares, if update_rows() is running) all the rows
 void build_rows( void );

 /// rewrites the rows after a change of the data [see build_rows()]
 void update_rows( ModParam issueAMod );

 /// true if build_rows() is generating rather than comparing
 bool generating_rows( void ) const { return( f_row_cmp == nullptr ); }

 /// writes the next row of a group, or compares it
 void push_row( std::vector< FRowConstraint > & rows ,
                LinearFunction::v_coeff_pair && vars , double lhs ,
                double rhs );

 /// sets the bounds of a BoxConstraint, or compares them
 void put_box( BoxConstraint & box , ColVariable * var , double lhs ,
               double rhs );

 /// registers a group of rows when generating
 void add_rows( std::vector< FRowConstraint > & rows , std::string && name );

 /// rewrites the coefficient of the design in the Objective
 void update_design_cost( ModParam issueAMod );

 /// rewrites the coefficients of the Objective for the given instants
 void update_objective( const Subset & instants , ModParam issueAMod );

 /// issues a ConversionUnitBlockMod
 void issue_mod( int type , Index index , Subset && subset ,
                 ModParam issuePMod );

 /// the common part of the setters of a datum indexed over [ T ][ col ]
 void guts_of_set( MAdbl & data , double dflt , int type , Index col ,
                   MF_dbl_it values , Subset && subset , bool ordered ,
                   ModParam issuePMod , ModParam issueAMod ,
                   bool rows , const std::string & who );

 /// the setters over the entries [ col ][ t ], column after column
 void flat_set( MAdbl & data , double dflt , int type , Index ncol ,
                MF_dbl_it values , Subset && subset , bool ordered ,
                ModParam issuePMod , ModParam issueAMod , bool rows ,
                const std::string & who );

 /// the common part of the setters of a datum indexed over [ T ]
 void guts_of_set( std::vector< double > & data , int type ,
                   MF_dbl_it values , Subset && subset , bool ordered ,
                   ModParam issuePMod , ModParam issueAMod ,
                   const std::string & who );

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS -----------------------------*/
/*--------------------------------------------------------------------------*/

 Index f_G = 1;                       ///< the number of generators
 Index f_M = 0;                       ///< the number of operating rows
 Index f_K = 0;                       ///< the number of ramp combinations
 Index f_commitment_generator = 0;    ///< the commitment generator

 MAdbl v_MinPower;                    ///< [ T ][ G ] minimum power
 MAdbl v_MaxPower;                    ///< [ T ][ G ] maximum power
 boost::multi_array< double , 3 > v_OperatingMatrix;
 ///< [ T or 1 ][ M ][ G ] operating rows, the first index the instant
 MAdbl v_OperatingLHS;                ///< [ T ][ M ], empty = -INF
 MAdbl v_OperatingRHS;                ///< [ T ][ M ], empty = +INF
 MAdbl v_RampMatrix;                  ///< [ K ][ G ], empty = identity
 MAdbl v_DeltaRampUp;                 ///< [ T ][ K ], empty = none
 MAdbl v_DeltaRampDown;               ///< [ T ][ K ], empty = none
 MAdbl v_StartUpLimit;                ///< [ T ][ G ], empty = none
 MAdbl v_StartUpLowerLimit;           ///< [ T ][ G ], empty = none
 MAdbl v_ShutDownLimit;               ///< [ T ][ G ], empty = none
 MAdbl v_ShutDownLowerLimit;          ///< [ T ][ G ], empty = none
 MAdbl v_PrimaryRho;                  ///< [ T ][ G ], empty = 0
 MAdbl v_SecondaryRho;                ///< [ T ][ G ], empty = 0
 MAdbl v_ReserveDirection;            ///< [ G ][ G ], empty = identity
 MAdbl v_PrimaryCost;                 ///< [ T ][ G ], empty = 0
 MAdbl v_SecondaryCost;               ///< [ T ][ G ], empty = 0
 MAdbl v_LinearTerm;                  ///< [ T ][ G ], empty = 0
 MAdbl v_QuadTerm;                    ///< [ T ][ G ], empty = 0

 std::vector< double > v_ConstTerm;          ///< [ T ], empty = 0
 std::vector< double > v_StartUpCost;        ///< [ T ], empty = 0
 std::vector< double > v_ShutDownCost;       ///< [ T ], empty = 0
 std::vector< double > v_FixedConsumption;   ///< [ T ], empty = none
 std::vector< double > v_InertiaCommitment;  ///< [ T ], empty = none

 std::vector< std::vector< double > > v_InertiaPower;
 ///< [ G ][ T ] inertia per unit of power, an empty row being none

 std::vector< double > v_InitialPower;  ///< [ G ] initial powers
 std::vector< Index > v_Delay;          ///< [ G ] delays, empty = 0
 std::vector< int > v_CyclicDelay;      ///< [ G ] cyclic, empty = none
 MAdbl v_InitialDelayedPower;           ///< [ H ][ G ], empty = 0
 int f_InitUpDownTime = 0;              ///< initial up/down time
 Index f_MinUpTime = 1;                 ///< minimum up time
 Index f_MinDownTime = 1;               ///< minimum down time
 double f_scale = 1;                    ///< the scale factor
 Index f_Modules = 1;                   ///< the number of modules
 Index f_InitialCommitment = 0;         ///< modules on before, 0 = all
 bool f_has_design = false;             ///< the unit has a design
 double f_InvestmentCost = 0;           ///< cost of one unit of design
 double f_MinCapacityDesign = 0;        ///< lower bound of the design
 double f_MaxCapacityDesign = -1;       ///< upper bound, signed

 std::vector< int > v_sign;  ///< [ G ] the sign of a reserve generator
 bool f_has_primary = false;    ///< some generator offers primary reserve
 bool f_has_secondary = false;  ///< some generator offers secondary one

 // the Variable
 ColVariable design;                       ///< x, if has_design()
 std::vector< ColVariable > v_commitment;  ///< u[ t ]
 std::vector< ColVariable > v_start_up;    ///< v[ t ]
 std::vector< ColVariable > v_shut_down;   ///< w[ t ]
 std::vector< std::vector< ColVariable > > v_active_power;
 ///< p[ g ][ s ], the T delivered powers then the tail [see op()]
 std::vector< std::vector< ColVariable > > v_primary_reserve;
 ///< pr[ g ][ t ], an empty row for a generator with no primary reserve
 std::vector< std::vector< ColVariable > > v_secondary_reserve;
 ///< sc[ g ][ t ], an empty row for a generator with no secondary reserve

 // the Constraint
 std::vector< FRowConstraint > Logical_Const;
 std::vector< FRowConstraint > MinUp_Const;
 std::vector< FRowConstraint > MinDown_Const;
 std::vector< FRowConstraint > MinPower_Const;
 std::vector< FRowConstraint > MaxPower_Const;
 std::vector< FRowConstraint > Operating_Const;
 std::vector< FRowConstraint > RampUp_Const;
 std::vector< FRowConstraint > RampDown_Const;
 std::vector< FRowConstraint > PrimaryRho_Const;
 std::vector< FRowConstraint > SecondaryRho_Const;
 std::vector< FRowConstraint > ReserveBound_Const;
 std::vector< FRowConstraint > ReserveOperating_Const;
 std::vector< FRowConstraint > ReserveRamp_Const;
 std::vector< FRowConstraint > CommitmentDesign_Const;
 std::vector< FRowConstraint > InitialDesign_Const;

 std::vector< BoxConstraint > Commitment_Bound;     ///< bounds of u[ t ]
 BoxConstraint ShutDownZero_Bound;                  ///< bound of w[ 0 ]
 std::vector< BoxConstraint > ShutDownInitial_Bound;
 ///< bounds of w[ t ], t >= 1, with several modules
 BoxConstraint DesignBound;                         ///< bounds of x
 std::vector< std::vector< BoxConstraint > > ActivePower_Bound;
 ///< [ g ][ t ] bounds of the power
 std::vector< ZOConstraint > StartUp_Binary_Bound;   ///< 0 <= v <= 1
 std::vector< ZOConstraint > ShutDown_Binary_Bound;  ///< 0 <= w <= 1

 FRealObjective objective;  ///< the Objective

 std::vector< std::vector< FRowConstraint > * > f_row_groups;
 ///< the groups of rows, in the order build_rows() writes them

 RowCmp * f_row_cmp = nullptr;  ///< not nullptr while comparing the rows

 bool f_counting = false;  ///< true while the rows are being counted

 std::map< std::vector< FRowConstraint > * , Index > f_row_count;
 ///< the rows of each group: their number, then the next to be written

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;  ///< registration in the factory

 /// registers the setters in the methods factory
 static void static_initialization( void )
 {
  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_linear_term" ,
   & ConversionUnitBlock::set_linear_term );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_linear_term" ,
   & ConversionUnitBlock::set_linear_term );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_quad_term" ,
   & ConversionUnitBlock::set_quad_term );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_quad_term" ,
   & ConversionUnitBlock::set_quad_term );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_primary_spinning_reserve_cost" ,
   & ConversionUnitBlock::set_primary_spinning_reserve_cost );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_primary_spinning_reserve_cost" ,
   & ConversionUnitBlock::set_primary_spinning_reserve_cost );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_secondary_spinning_reserve_cost" ,
   & ConversionUnitBlock::set_secondary_spinning_reserve_cost );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_secondary_spinning_reserve_cost" ,
   & ConversionUnitBlock::set_secondary_spinning_reserve_cost );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_maximum_power" ,
   & ConversionUnitBlock::set_maximum_power );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_maximum_power" ,
   & ConversionUnitBlock::set_maximum_power );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_minimum_power" ,
   & ConversionUnitBlock::set_minimum_power );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_minimum_power" ,
   & ConversionUnitBlock::set_minimum_power );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_operating_rhs" ,
   & ConversionUnitBlock::set_operating_rhs );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_operating_rhs" ,
   & ConversionUnitBlock::set_operating_rhs );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_operating_lhs" ,
   & ConversionUnitBlock::set_operating_lhs );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_operating_lhs" ,
   & ConversionUnitBlock::set_operating_lhs );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_delta_ramp_up" ,
   & ConversionUnitBlock::set_delta_ramp_up );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_delta_ramp_up" ,
   & ConversionUnitBlock::set_delta_ramp_up );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_delta_ramp_down" ,
   & ConversionUnitBlock::set_delta_ramp_down );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_delta_ramp_down" ,
   & ConversionUnitBlock::set_delta_ramp_down );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_const_term" ,
   & ConversionUnitBlock::set_const_term );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_const_term" ,
   & ConversionUnitBlock::set_const_term );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_startup_costs" ,
   & ConversionUnitBlock::set_startup_costs );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_startup_costs" ,
   & ConversionUnitBlock::set_startup_costs );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_shutdown_costs" ,
   & ConversionUnitBlock::set_shutdown_costs );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_shutdown_costs" ,
   & ConversionUnitBlock::set_shutdown_costs );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_initial_power" ,
   & ConversionUnitBlock::set_initial_power );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_initial_power" ,
   & ConversionUnitBlock::set_initial_power );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_investment_cost" ,
   & ConversionUnitBlock::set_investment_cost );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_investment_cost" ,
   & ConversionUnitBlock::set_investment_cost );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::scale" ,
   & ConversionUnitBlock::scale );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::scale" ,
   & ConversionUnitBlock::scale );

  // scale() again under the name of the way of sizing a Block it is: k
  // identical copies of this unit, committed together

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::replicate" ,
   & ConversionUnitBlock::scale );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::replicate" ,
   & ConversionUnitBlock::scale );

  register_method< ConversionUnitBlock , MF_dbl_it , Subset && , bool >(
   "ConversionUnitBlock::set_initial_delayed_power" ,
   & ConversionUnitBlock::set_initial_delayed_power );

  register_method< ConversionUnitBlock , MF_dbl_it , Range >(
   "ConversionUnitBlock::set_initial_delayed_power" ,
   & ConversionUnitBlock::set_initial_delayed_power );

  register_method< ConversionUnitBlock , MF_int_it , Subset && , bool >(
   "ConversionUnitBlock::set_init_updown_time" ,
   & ConversionUnitBlock::set_init_updown_time );

  register_method< ConversionUnitBlock , MF_int_it , Range >(
   "ConversionUnitBlock::set_init_updown_time" ,
   & ConversionUnitBlock::set_init_updown_time );
  }

 };  // end( class( ConversionUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ConversionUnitBlockMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// the physical Modification of a ConversionUnitBlock
/** A ConversionUnitBlockMod says which datum of a ConversionUnitBlock has
 * changed (type()), for which generator, operating row or ramp combination
 * (index(), ignored for the data indexed over the instants only, and for
 * the initial powers, whose generators are in nms()), and at which instants
 * (nms()). */

class ConversionUnitBlockMod : public UnitBlockMod
{
 public:

 /// the types of ConversionUnitBlockMod
 enum CUB_mod_type
 {
  eSetLinT = eUBModLastParam ,  ///< linear cost of a generator
  eSetQuadT ,                   ///< quadratic cost of a generator
  eSetPrCost ,                  ///< primary reserve cost of a generator
  eSetScCost ,                  ///< secondary reserve cost of a generator
  eSetConstT ,                  ///< cost of being on
  eSetSUC ,                     ///< start-up cost
  eSetSDC ,                     ///< shut-down cost
  eSetMaxP ,                    ///< maximum power of a generator
  eSetMinP ,                    ///< minimum power of a generator
  eSetOpRHS ,                   ///< right-hand side of an operating row
  eSetOpLHS ,                   ///< left-hand side of an operating row
  eSetRampUp ,                  ///< ramp-up limit of a combination
  eSetRampDown ,                ///< ramp-down limit of a combination
  eSetInitP ,                   ///< initial powers (nms() = generators)
  eSetInitUD ,                  ///< initial up/down time
  eSetInitDelP ,                ///< powers delivered from before the horizon
  eSetInvCost ,                 ///< investment cost
  eCUBModLastParam    ///< first allowed parameter value for derived classes
  };

 /// what a Modification of the given type changes
 /** What a Modification of the given type changes of the model [see
  * Modification::ModConcern], which is what changes() says and what the
  * methods that issue it ask before issuing it [see Observer::issue_pmod()].
  * When it is not known for sure that a datum only enters the sides of the
  * Constraint, its change is said to be one of the Constraint. */

 static ModConcern changes_of( int type ) {
  switch( type ) {
   case( eSetLinT ):
   case( eSetQuadT ):
   case( eSetPrCost ):
   case( eSetScCost ):
   case( eSetConstT ):
   case( eSetSUC ):
   case( eSetSDC ):
   case( eSetInvCost ):
    return( ePhysObj );
   case( eSetMaxP ):
   case( eSetMinP ):
   case( eSetRampUp ):
   case( eSetRampDown ):
   case( eSetInitP ):
   case( eSetInitDelP ):
    return( ePhysCns );
   case( eSetInitUD ):
    return( ePhysCns | ePhysVars );
   case( eSetOpRHS ):
   case( eSetOpLHS ):
    return( ePhysSides );
   default:
    return( UnitBlockMod::changes_of( type ) );
   }
  }

 /// what the Modification changes, as its type says [see changes_of()]
 [[nodiscard]] ModConcern changes( void ) const override {
  return( changes_of( type() ) );
  }

 /// whether the Modification changes the power range of the unit
 /** Of the types of this class, those of the MaxPower and MinPower [see
  * UnitBlockMod::changes_power_range()]. */

 [[nodiscard]] bool changes_power_range( void ) const override {
  return( ( f_type == eSetMaxP ) || ( f_type == eSetMinP ) ||
          UnitBlockMod::changes_power_range() );
  }

 /// constructor: the ConversionUnitBlock, the type, the index, the instants
 ConversionUnitBlockMod( ConversionUnitBlock * const fblock , const int type ,
                         Block::Index index , Block::Subset && nms )
  : UnitBlockMod( fblock , type ) , f_index( index ) ,
    f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~ConversionUnitBlockMod() override = default;

 /// the generator, operating row or ramp combination of the change
 Block::Index index( void ) const { return( f_index ); }

 /// the instants of the change (the generators for eSetInitP)
 Block::c_Subset & nms( void ) const { return( f_nms ); }

 protected:

 /// prints the ConversionUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "ConversionUnitBlockMod[" << this << "]: type " << f_type
         << ", index " << f_index << " (# " << f_nms.size() << ")"
         << std::endl;
  }

 Block::Index f_index;  ///< the generator, row or combination
 Block::Subset f_nms;   ///< the instants

 };  // end( class( ConversionUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS ConversionUnitBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/

/// a [UnitBlock]Solution of a ConversionUnitBlock
/** The ConversionUnitBlockSolution class derives from UnitBlockSolution and
 * adds to the information stored there (the active powers, the commitment
 * of the commitment generator and the reserves) the start-up and shut-down
 * indicators of the unit, which are saved with the commitment (deriving
 * them from the commitment when the Solution is written back is right for
 * one schedule but not for a convex combination of several, whose
 * start-ups are the average of those of the schedules, not those of their
 * averaged commitment), the powers that the generators with a delay
 * deliver after the horizon, which are saved with the active powers, and
 * the design variable, if any. */

class ConversionUnitBlockSolution : public UnitBlockSolution
{
 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend ConversionUnitBlock;  ///< make ConversionUnitBlock friend

/*------- CONSTRUCTING AND DESTRUCTING ConversionUnitBlockSolution ---------*/

 /// constructor, it has nothing to do
 explicit ConversionUnitBlockSolution( void ) = default;

 /// destructor: it is virtual, and empty
 ~ConversionUnitBlockSolution() override = default;

/*---------- METHODS DESCRIBING THE BEHAVIOR OF THE SOLUTION ---------------*/

 /// reads the Solution, the format of serialize()
 void deserialize( const netCDF::NcGroup & group ) override;

 /// reads the solution of the ConversionUnitBlock \p block
 void read( const Block * block ) override;

 /// writes the solution into the ConversionUnitBlock \p block
 void write( Block * block ) override;

 /// writes the Solution in a netCDF::NcGroup
 /** The format is the one of UnitBlockSolution [see
  * UnitBlockSolution::serialize()], plus, if they are saved, the variables
  * "ConversionStartUp" and "ConversionShutDown", of type netCDF::NcDouble
  * and indexed over "TimeHorizon", the start-up and shut-down indicators of
  * the unit, and "ConversionTailPower", of type netCDF::NcDouble and indexed
  * over "NumberTailPowers", the powers delivered after the horizon by the
  * generators with a delay, generator after generator, and the scalar
  * "ConversionDesign", the value of the design variable, if the unit has
  * one. */

 void serialize( netCDF::NcGroup & group ) const override;

 /// returns a copy of this Solution with all values times \p factor
 ConversionUnitBlockSolution * scale( double factor ) const override;

 /// adds \p multiplier times \p solution to this Solution
 void sum( const Solution * solution , double multiplier ) override;

 /// returns a copy (an empty one if \p empty) of this Solution
 ConversionUnitBlockSolution * clone( bool empty = false ) const override;

 /// the start-up indicators, empty if not saved
 [[nodiscard]] const std::vector< double > & get_start_up( void ) const {
  return( v_start_up );
  }

 /// the shut-down indicators, empty if not saved
 [[nodiscard]] const std::vector< double > & get_shut_down( void ) const {
  return( v_shut_down );
  }

 /// the design, NaN if not saved
 [[nodiscard]] double get_design( void ) const { return( f_design ); }

 /// sets the start-up indicators
 void set_start_up( std::vector< double > && su ) {
  v_start_up = std::move( su );
  }

 /// sets the shut-down indicators
 void set_shut_down( std::vector< double > && sd ) {
  v_shut_down = std::move( sd );
  }

 protected:

 /// prints the ConversionUnitBlockSolution
 void print( std::ostream & output ) const override {
  output << "ConversionUnitBlockSolution [" << this << "]: " << std::endl;
  }

 private:

 std::vector< double > v_start_up;   ///< the start-up indicators, if saved

 std::vector< double > v_shut_down;  ///< the shut-down indicators, if saved

 std::vector< double > v_tail;  ///< the powers delivered after the horizon

 double f_design = std::numeric_limits< double >::quiet_NaN();
 ///< the design, NaN if not saved

 SMSpp_insert_in_factory_h;  ///< registration in the factory

 };  // end( class( ConversionUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __ConversionUnitBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ConversionUnitBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
