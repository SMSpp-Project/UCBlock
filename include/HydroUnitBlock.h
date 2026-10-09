/*--------------------------------------------------------------------------*/
/*------------------------- File HydroUnitBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class HydroUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" hydro unit
 * of a Unit Commitment Problem.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HydroUnitBlock
 #define __HydroUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColVariable.h"

#include "FRowConstraint.h"

#include "DQuadFunction.h"

#include "UnitBlock.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS HydroUnitBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a hydro valley
/** HydroUnitBlock implements the Block concept [see Block.h] for a hydro
 * valley, i.e., a set of reservoirs joined by plants (turbines, pumps and
 * spillways), as a UnitBlock of a UCBlock. Its place in the complete model
 * (the linking constraints that its powers and reserves enter, the
 * conventions on instants, units and signs) is described in
 * \ref ucblock_model. A HydroSystemUnitBlock, in turn, groups several
 * HydroUnitBlock together with the future cost of the water left in their
 * reservoirs. We see the valley as a directed graph whose nodes are the
 * reservoirs \f$ n \in \mathcal{N}^{hy} = \{ 0 , \ldots , R - 1 \} \f$
 * ("NumberReservoirs") and whose arcs \f$ l \in \mathcal{L}^{hy} \f$
 * ("NumberArcs") are the plants, and each arc is one generator of the
 * UnitBlock. Arc \f$ l \f$ goes from its start reservoir \f$ s( l ) \f$
 * ("StartArc") to its end reservoir \f$ e( l ) \f$ ("EndArc"), and
 * \f$ e( l ) = R \f$ means that the water leaves the valley. Several plants
 * between the same two reservoirs are parallel arcs. For each instant \f$ t
 * \in \mathcal{T} = \{ 0 , \ldots , T - 1 \} \f$ the variables are the
 * volumes \f$ v^{hy}_{n,t} \f$ of the reservoirs at the end of the instant,
 * the flows \f$ f_{t,l} \f$ and the active powers \f$ p^{ac}_{t,l} \f$ of the
 * arcs and, when the UCBlock requires them and the data allow them, the
 * primary and secondary reserves \f$ p^{pr}_{t,l} \f$ and
 * \f$ p^{sc}_{t,l} \f$ of the arcs. Their constraints (see
 * generate_abstract_constraints()) are the water balance of each reservoir,
 * with the delays of the arcs, the bounds on the volumes and on the flows,
 * and the ramps of the flows. Then come the relation between the flow and the
 * power of turbines and pumps, and the bounds on the power together with the
 * reserves. The model is linear and continuous, and its cost (see
 * generate_objective()) is linear in the powers.
 *
 * All the data are given per instant, since no length of the time step
 * appears anywhere (see \ref ucbm_conv_time). The volumes, the inflows and
 * the flows share one unit of volume: the flow of an arc is the volume of
 * water that goes through it during one instant, i.e., with instants of
 * \f$ \Delta t \f$ hours, \f$ 3600 \, \Delta t \f$ times a flow rate in
 * m\f$ ^3 \f$/s. Hence, a coefficient of the power curve given in MW per
 * m\f$ ^3 \f$/s has to be divided by \f$ 3600 \, \Delta t \f$. A ramp is the
 * largest change of the flow between two consecutive instants, in the same
 * unit, and a delay is an integer number of instants (a travel time of
 * \f$ d \f$ hours is given as \f$ \lceil d / \Delta t \rceil \f$ instants).
 * HydroUnitBlock does not implement the scale factor of UnitBlock, i.e.,
 * get_scale() is 1.
 *
 * Some features of a valley are represented through the data alone. A
 * reservoir without storage (a channel) has
 * \f$ V^{mn}_{n,t} = V^{mx}_{n,t} = 0 \f$, and therefore its water balance
 * imposes that the water entering it leaves it at the same instant. A target
 * on a volume at some instant (e.g., a level to be reached in the middle or
 * at the end of the horizon) is a time-dependent bound on that volume. A
 * reversible plant is a turbine and a pump between the same two reservoirs. A
 * spillway is an arc with a single piece whose linear term \f$ \rho^{hy} \f$
 * is zero and whose constant term is zero as well (otherwise the arc would
 * produce that power whatever its flow); hence, it releases water without
 * producing power, and its capacity is its maximum flow. Finally, the value
 * of the water left at the end of the horizon is the polyhedral function of a
 * HydroSystemUnitBlock.
 *
 * We do not model the following features. The power of a turbine is a concave
 * function of its flow alone, the same at each instant. Hence, operating
 * points that are discrete (with the staircase or incremental formulations
 * that describe them through integer variables) and a power curve that is not
 * concave, or that changes over time, are not represented. Also, the concave
 * curve is an outer approximation of the true one, and it is exact at an
 * optimum where producing less than the curve allows is never profitable.
 * Neither are the dependence of the power on the head (i.e., on the volumes
 * of the reservoirs), a coefficient of a pump that changes over time, and a
 * reserve that is a datum of the operating point (rather than a variable
 * bounded by (1)-(6) of generate_abstract_constraints()). Since the model has
 * no integer variable, three rules are not imposed either, i.e., (i) that a
 * plant spills only when it turbines at its maximum, (ii) that the turbine
 * and the pump of a reversible plant do not work at the same instant, nor at
 * two consecutive instants when the mode changes, and (iii) that the
 * operating point of a plant does not change twice in three consecutive
 * instants. Pumping and turbining at the same time wastes energy, since a
 * pump usually consumes more than the turbine produces with the same water,
 * and it is therefore unlikely (but not excluded) at an optimum. Finally, the
 * ramps (10) and (11) of generate_abstract_constraints() bound the flow of
 * each arc, while the sum of the flows of a turbine and of its spillway is
 * not bounded, and the water in transit at the borders of the horizon is not
 * accounted for (see (12) there). */

class HydroUnitBlock : public UnitBlock
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
/** Constructor of HydroUnitBlock, taking possibly a pointer of its father
 * Block. */

 explicit HydroUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of HydroUnitBlock

 virtual ~HydroUnitBlock() override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HydroUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals": a variable below
 * that has a dimension "NumberIntervals" may have, in that dimension, size
 * 1 (the value holds at every instant), size "NumberIntervals" (each value
 * holds in one of the intervals of "ChangeIntervals") or, if
 * "NumberIntervals" is not provided, size "TimeHorizon", and it is turned
 * into a matrix with one entry per instant. A variable indexed over
 * "NumberArcs" or "NumberReservoirs" alone, "StartArc" and "EndArc"
 * excepted, may also be a scalar, which is then the value of every arc or
 * reservoir. With the notation of the class
 * (see generate_abstract_constraints() for the constraints), the
 * netCDF::NcGroup must then also contain:
 *
 * - The dimension "NumberReservoirs", the number \f$ R \f$ of the
 *   reservoirs of the valley. The dimension is optional, if it is not
 *   provided then it is taken to be 1: a single reservoir can still have
 *   several plants (see "NumberArcs").
 *
 * - The dimension "NumberArcs", the number of arcs (plants) of the valley,
 *   i.e., the number of generators of the UnitBlock. The dimension is
 *   optional, if it is not provided then it is taken to be 1.
 *
 * - The variable "StartArc", of type netCDF::NcUint and indexed over the
 *   dimension "NumberArcs": StartArc[ l ] is the reservoir \f$ s( l ) \in
 *   \{ 0 , \ldots , R - 1 \} \f$ from which arc \f$ l \f$ starts. A positive
 *   flow along arc \f$ l \f$ (a turbine) takes water from \f$ s( l ) \f$ and
 *   delivers it to \f$ e( l ) \f$, a negative flow (a pump) does the
 *   opposite.
 *
 * - The variable "EndArc", of type netCDF::NcUint and indexed over the
 *   dimension "NumberArcs": EndArc[ l ] is the reservoir \f$ e( l ) \in
 *   \{ 0 , \ldots , R \} \f$ at which arc \f$ l \f$ ends, the value \f$ R \f$
 *   (not a reservoir) meaning that the water flowing along the arc leaves
 *   the valley, as it happens after the most downstream plants, and is no
 *   longer counted. An arc with \f$ s( l ) = e( l ) \f$ (a self-loop) is not
 *   allowed, while several arcs between the same two reservoirs are; often
 *   the same equipment can work both as a turbine and as a pump, and it is
 *   then represented by two parallel arcs with different flow bounds (see
 *   "MinFlow" and "MaxFlow").
 *
 *   "StartArc" and "EndArc" are optional, but they are given together, and
 *   they are required if \f$ R > 1 \f$: without them every arc leaves the
 *   only reservoir and the valley, and the delays (see "UphillFlow" and
 *   "DownhillFlow") play no role.
 *
 * - The variables "MinFlow" and "MaxFlow", of type netCDF::NcDouble and
 *   indexed over the dimensions "NumberIntervals" and "NumberArcs": the
 *   matrices \f$ F^{mn}_{t,l} \f$ and \f$ F^{mx}_{t,l} \f$ of the bounds on
 *   the flow of arc \f$ l \f$ at instant \f$ t \f$, which must satisfy
 *   \f$ F^{mn}_{t,l} \leq F^{mx}_{t,l} \f$. Both are optional, and 0 if not
 *   provided. Their signs give the kind of the arc at instant \f$ t \f$: it
 *   is a turbine if \f$ 0 \leq F^{mn}_{t,l} \f$ and \f$ F^{mx}_{t,l} > 0
 *   \f$, a pump if \f$ F^{mn}_{t,l} < 0 \f$ and \f$ F^{mx}_{t,l} \leq 0 \f$,
 *   and it is idle (\f$ f_{t,l} = 0 \f$) if \f$ F^{mn}_{t,l} =
 *   F^{mx}_{t,l} = 0 \f$. An arc with \f$ F^{mn}_{t,l} < 0 < F^{mx}_{t,l}
 *   \f$ at some instant is rejected, and so is an arc that is a turbine at
 *   some instant and a pump at another (it may be idle at some instants):
 *   the flow-to-power function of a turbine is a concave piecewise-linear
 *   function with possibly many pieces, while that of a pump is linear (see
 *   "NumberPieces", "LinearTerm" and "ConstantTerm"). A plant that can work
 *   in both modes is split, when the data are prepared, into a turbine and a
 *   pump, and nothing then prevents the two from working at the same
 *   instant (see the description of the class).
 *
 * - The variables "MinVolumetric" and "MaxVolumetric", of type
 *   netCDF::NcDouble and indexed over the dimensions "NumberReservoirs" and
 *   "NumberIntervals": the matrices \f$ V^{mn}_{n,t} \f$ and
 *   \f$ V^{mx}_{n,t} \f$ of the bounds on the volume of reservoir \f$ n \f$
 *   at the end of instant \f$ t \f$, which must satisfy \f$ 0 \leq
 *   V^{mn}_{n,t} \leq V^{mx}_{n,t} \f$. Both are optional, and 0 if not
 *   provided, which for "MaxVolumetric" means a reservoir without storage.
 *   The two bounds may coincide at some instants, e.g., to fix the volume at
 *   the end of the horizon or at a given instant, and at all instants, for a
 *   channel in which the water does not stop.
 *
 * - The variable "Inflows", of type netCDF::NcDouble and indexed over the
 *   dimensions "NumberReservoirs" and "TimeHorizon": the matrix
 *   \f$ A_{n,t} \f$ of the volume of water that enters reservoir \f$ n \f$
 *   during instant \f$ t \f$ by natural causes (rain, melting ice, rivers
 *   that are not controlled), net of what leaves it in the same way
 *   (evaporation, withdrawals for other uses), and is therefore available at
 *   the end of instant \f$ t \f$; it is a volume per instant, in the unit of
 *   the volumes, and it may be negative. The variable is optional, and 0 if
 *   not provided.
 *
 * - The variables "MinPower" and "MaxPower", of type netCDF::NcDouble and
 *   indexed over the dimensions "NumberIntervals" and "NumberArcs": the
 *   matrices \f$ P^{mn}_{t,l} \f$ and \f$ P^{mx}_{t,l} \f$ of the bounds on
 *   the active power of arc \f$ l \f$ at instant \f$ t \f$, which must
 *   satisfy \f$ P^{mn}_{t,l} \leq P^{mx}_{t,l} \f$ (the power of a pump
 *   being nonpositive). Both are optional, and 0 if not provided.
 *
 * - The variables "DeltaRampUp" and "DeltaRampDown", of type
 *   netCDF::NcDouble and indexed over the dimensions "NumberIntervals" and
 *   "NumberArcs": the matrices \f$ \Delta^+_{t,l} \f$ and
 *   \f$ \Delta^-_{t,l} \f$ of the largest increase and decrease of the flow
 *   of arc \f$ l \f$ between instants \f$ t - 1 \f$ and \f$ t \f$, in the
 *   unit of the flows (a gradient of \f$ G \f$ m\f$ ^3 \f$/s per hour is
 *   \f$ 3600 \, G \, \Delta t^2 \f$ with instants of \f$ \Delta t \f$
 *   hours). Both are optional: if "DeltaRampUp" (respectively,
 *   "DeltaRampDown") is not provided, there is no ramp-up (respectively,
 *   ramp-down) constraint.
 *
 * - The variables "PrimaryRho" and "SecondaryRho", of type
 *   netCDF::NcDouble and indexed over the dimensions "NumberIntervals" and
 *   "NumberArcs": the matrices \f$ \rho^{pr}_{t,l} \f$ and
 *   \f$ \rho^{sc}_{t,l} \f$ of the largest fraction of the active power of
 *   arc \f$ l \f$ at instant \f$ t \f$ that can be primary and secondary
 *   reserve. Both are optional: if "PrimaryRho" (respectively,
 *   "SecondaryRho") is not provided, the unit gives no primary
 *   (respectively, secondary) reserve and has no such Variable. Only
 *   turbines give reserves: \f$ \rho^{pr}_{t,l} = \rho^{sc}_{t,l} = 0 \f$
 *   whenever arc \f$ l \f$ is a pump at instant \f$ t \f$.
 *
 * - The dimension "TotalNumberPieces" and the variable "NumberPieces", of
 *   type netCDF::NcUint and indexed over the dimension "NumberArcs":
 *   NumberPieces[ l ] is the number \f$ | \mathcal{J}_l | \geq 1 \f$ of the
 *   pieces of the flow-to-power function of arc \f$ l \f$, which is 1 for
 *   a pump, and "TotalNumberPieces" is \f$ \sum_l | \mathcal{J}_l | \f$.
 *   Both are optional: without "NumberPieces" every arc has one piece, and
 *   without "TotalNumberPieces" the sum is computed.
 *
 * - The variables "LinearTerm" and "ConstantTerm", of type netCDF::NcDouble
 *   and indexed over the set \f$ \{ 0 , \ldots ,
 *   \mathrm{TotalNumberPieces} - 1 \} \f$: the coefficients
 *   \f$ \rho^{hy}_j \f$ and \f$ P^{hy}_j \f$ of the pieces
 *   \f$ \rho^{hy}_j f + P^{hy}_j \f$ of the flow-to-power functions, whose
 *   minimum over \f$ j \in \mathcal{J}_l \f$ bounds the power of turbine
 *   \f$ l \f$ (see (7) of generate_abstract_constraints()), while a pump
 *   uses only \f$ \rho^{hy}_j \f$ of its piece. The pieces are numbered arc
 *   after arc, in the order of the arcs:
 *   - index 0 is the first piece of arc 0, index 1 its second piece, and
 *     so on up to index NumberPieces[ 0 ] - 1, its last piece;
 *   - index NumberPieces[ 0 ] is the first piece of arc 1, and so on,
 *
 *   so that the index of a piece is that of its arc when every arc has one
 *   piece. "LinearTerm" is required if some arc is a turbine at some
 *   instant; if it is not provided, a pump has \f$ p^{ac}_{t,l} = f_{t,l}
 *   \f$. "ConstantTerm" is optional, and 0 if not provided. A piece known
 *   as the tangent \f$ P_j + \rho_j ( f - \bar f_j ) \f$ to the curve at the
 *   flow \f$ \bar f_j \f$ is given by \f$ \rho^{hy}_j = \rho_j \f$ and
 *   \f$ P^{hy}_j = P_j - \rho_j \bar f_j \f$, with \f$ \rho_j \f$ and
 *   \f$ \bar f_j \f$ converted to the unit of the flows.
 *
 * - The variable "ActivePowerCost", of type netCDF::NcDouble and indexed
 *   over the dimension "NumberArcs": the cost \f$ b_l \f$ of one unit of
 *   active power of arc \f$ l \f$ during one instant. The variable is
 *   optional, and 0 if not provided.
 *
 * - The variable "InertiaPower", of type netCDF::NcDouble and indexed over the
 *   dimensions "NumberArcs" and "NumberIntervals": the matrix
 *   \f$ h^p_{t,l} \f$ of the coefficients of the active power of arc \f$ l \f$
 *   at instant \f$ t \f$ in the inertia constraints of the UCBlock, i.e.,
 *   \f$ h^p_{t,l} = 1.2 H_{t,l} \f$ for a plant whose inertia constant is
 *   \f$ H_{t,l} \f$ (see \ref ucbm_link_in). The variable is optional, and if
 *   it is not provided the unit gives no inertia.
 *
 * - The variable "InitialFlowRate", of type netCDF::NcDouble and indexed
 *   over the dimension "NumberArcs": the flow \f$ F^0_l \f$ of arc \f$ l \f$
 *   at instant \f$ -1 \f$, which only the ramp constraints of instant 0
 *   use. The variable is optional, and 0 if not provided.
 *
 * - The variable "InitialVolumetric", of type netCDF::NcDouble and indexed
 *   over the dimension "NumberReservoirs": the volume \f$ V^0_n \f$ of
 *   reservoir \f$ n \f$ at instant \f$ -1 \f$, i.e., at the beginning of the
 *   horizon. A negative value means the cyclic closure, i.e., the volume at
 *   instant \f$ -1 \f$ is the variable \f$ v^{hy}_{n,T-1} \f$. The variable
 *   is optional, and 0 if not provided.
 *
 *   Under the cyclic closure the water balances (12) of the reservoir, summed
 *   over the instants, give \f$ \sum_t ( v^{hy}_{n,t} - v^{hy}_{n,t-1} ) =
 *   0 \f$ on the left (with \f$ v^{hy}_{n,-1} = v^{hy}_{n,T-1} \f$), hence
 *   the equality between the total inflow \f$ \sum_t A_{n,t} \f$ and the
 *   total net outflow, i.e., the flows of the arcs that leave the reservoir
 *   minus those of the arcs that reach it, each counted at the instants at
 *   which (12) of generate_abstract_constraints() counts it. Only this
 *   aggregate condition follows: within the horizon the reservoir stores
 *   the inflow of an instant that its arcs cannot carry away, and the arcs
 *   that reach it add water. With the flow bounds it implies
 *   \f[
 *     \sum_{ ( t , l ) \in O_n } F^{mn}_{t,l}
 *       - \sum_{ ( t , l ) \in I_n } F^{mx}_{t,l} \; \leq \;
 *     \sum_{ t \in \mathcal{T} } A_{n,t} \; \leq \;
 *     \sum_{ ( t , l ) \in O_n } F^{mx}_{t,l}
 *       - \sum_{ ( t , l ) \in I_n } F^{mn}_{t,l} \; ,
 *   \f]
 *   where \f$ O_n \f$ (resp. \f$ I_n \f$) are the pairs of an arc that
 *   leaves (resp. reaches) \f$ n \f$ and of an instant at which its flow
 *   is in a balance of \f$ n \f$ (all the arcs and instants for a single
 *   reservoir without "StartArc" and "EndArc"); this necessary condition
 *   is checked when the data are read (see check_data_consistency()), and
 *   a reservoir whose total inflow exceeds what its arcs can carry, e.g.,
 *   because its spillway is too small, makes deserialize() throw.
 *
 * - The variables "UphillFlow" and "DownhillFlow", of type netCDF::NcInt
 *   and netCDF::NcUint and indexed over the dimension "NumberArcs": the
 *   delays \f$ \tau^{up}_l \f$ and \f$ \tau^{dn}_l \geq 0 \f$, in instants,
 *   of arc \f$ l \f$. The flow of arc \f$ l \f$ at instant \f$ s \f$ is
 *   withdrawn from \f$ s( l ) \f$ at instant \f$ s + \tau^{up}_l \f$ and it
 *   reaches \f$ e( l ) \f$ at instant \f$ s + \tau^{dn}_l \f$ (see (12) of
 *   generate_abstract_constraints()). Both are optional, and 0 if not
 *   provided.
 *
 *   The plant can be far enough from its reservoirs for the water to take
 *   one or more instants (especially if these are short, say 5 or 15
 *   minutes) to go from the start reservoir to the plant, and from the plant
 *   to the end reservoir, as depicted below for a turbine:
 *   \verbatim
      s(l) >==================> [ PLANT ] >==================> e(l)
            water leaves s(l)   works at s   water reaches e(l)
            at s + UphillFlow                at s + DownhillFlow
     \endverbatim
 *   The delay \f$ \tau^{dn}_l \f$ is the travel time from the plant to
 *   \f$ e( l ) \f$. When the pipe from \f$ s( l ) \f$ to the plant has to be
 *   filled before the plant can work, the water leaves \f$ s( l ) \f$ before
 *   it reaches the plant, and \f$ \tau^{up}_l \f$ is the opposite of the
 *   travel time, i.e., negative. When, instead, the pipe is full, it works
 *   as a small reservoir of its own: starting the plant creates a
 *   depression that travels up the pipe, and the water starts leaving
 *   \f$ s( l ) \f$ only when the depression reaches it, i.e.,
 *   \f$ \tau^{up}_l \f$ is positive. These are rather crude approximations
 *   of the physical behavior, which are however accurate enough for this
 *   setting.
 *
 * - The variable "ReferenceSchedule", of type netCDF::NcDouble and indexed
 *   over the dimension "NumberIntervals": a total active power
 *   \f$ \hat p_t \f$ of the valley at each instant, from which the objective
 *   penalizes the absolute deviation (see generate_objective()). The
 *   variable is optional, and without it there is no such penalty.
 *
 * - The variables "MinReactivePower" and "MaxReactivePower", of type
 *   netCDF::NcDouble and indexed over the dimensions "NumberIntervals" and
 *   "NumberArcs": the bounds on the reactive power of each arc, which is a
 *   Variable only if the UCBlock asks for it (see
 *   UnitBlock::set_reactive_power()). Both are optional; if one of them is
 *   provided with only zero entries, it is taken as not provided.
 *
 * The consistency of the data is checked as described in
 * check_data_consistency(), which throws an exception if it does not
 * hold. */

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
 /// generate the abstract variables of the HydroUnitBlock
 /** The HydroUnitBlock has the following groups of Variable, all of them
  * boost::multi_array< ColVariable , 2 > but the last one:
  *
  * - the volumes \f$ v^{hy}_{n,t} \geq 0 \f$, indexed [ n ][ t ], group
  *   "v_hydro";
  *
  * - the flows \f$ f_{t,l} \f$, indexed [ l ][ t ], group "f_hydro",
  *   nonnegative if \f$ F^{mn}_{t,l} \geq 0 \f$ and nonpositive if
  *   \f$ F^{mx}_{t,l} \leq 0 \f$ (both, i.e., zero, for an idle arc);
  *
  * - the active powers \f$ p^{ac}_{t,l} \f$, indexed [ l ][ t ], group
  *   "p_hydro", nonnegative if \f$ P^{mn}_{t,l} \geq 0 \f$ and nonpositive
  *   if \f$ P^{mx}_{t,l} \leq 0 \f$;
  *
  * - the reactive powers, indexed [ l ][ t ], group "q_hydro", only if the
  *   UCBlock asks for them (see UnitBlock::set_reactive_power());
  *
  * - the primary reserves \f$ p^{pr}_{t,l} \geq 0 \f$, indexed [ l ][ t ],
  *   group "pr_hydro", only if the UCBlock asks for the primary reserve
  *   (see UnitBlock::set_reserve_vars()) and "PrimaryRho" is provided;
  *
  * - the secondary reserves \f$ p^{sc}_{t,l} \geq 0 \f$, indexed
  *   [ l ][ t ], group "sr_hydro", only if the UCBlock asks for the
  *   secondary reserve and "SecondaryRho" is provided;
  *
  * - the deviations \f$ \delta^{rs}_t \geq 0 \f$ from the reference schedule,
  *   a std::vector< ColVariable > indexed by \f$ t \f$, group
  *   "v_absh_refschd", only if "ReferenceSchedule" is provided. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraints of the HydroUnitBlock
 /** This method generates the static constraints of the HydroUnitBlock,
  * with the notation of the description of the class and of deserialize().
  * Each group of constraints is a boost::multi_array, registered as a static
  * constraint under the name given below; the indices of the arrays start
  * at 0, and an array indexed by instant and arc has
  * get_time_horizon() \f$ \times \f$ get_number_generators() entries. At
  * every instant \f$ t \f$ an arc \f$ l \f$ is a turbine, a pump or idle,
  * according to the signs of \f$ F^{mn}_{t,l} \f$ and \f$ F^{mx}_{t,l} \f$
  * (see deserialize()); the reserve Variable exist only if the UCBlock asks
  * for the reserve and "PrimaryRho" (respectively, "SecondaryRho") is
  * provided, and a term with a reserve that does not exist is absent from
  * the rows below.
  *
  * - Power and reserves (MaxPowerPrimarySecondary_Const and
  *   MinPowerPrimarySecondary_Const, both boost::multi_array<
  *   FRowConstraint , 2 > indexed [ t ][ l ], groups
  *   "MaxPowerPrimarySecondary_HydroUnit" and
  *   "MinPowerPrimarySecondary_HydroUnit"): if the unit has some reserve
  *   Variable, for \f$ t \in \mathcal{T} \f$ and \f$ l \in \mathcal{L}^{hy}
  *   \f$
  *   \f[
  *     p^{ac}_{t,l} + p^{pr}_{t,l} + p^{sc}_{t,l} \leq P^{mx}_{t,l} \; ,
  *     \tag{1}
  *   \f]
  *   \f[
  *     p^{ac}_{t,l} - p^{pr}_{t,l} - p^{sc}_{t,l} \geq P^{mn}_{t,l} \; ,
  *     \tag{2}
  *   \f]
  *   i.e., the reserves are symmetric: the plant must be able to move its
  *   power by their amount both upward and downward. If the unit has no
  *   reserve Variable, (1) and (2) are replaced by the bounds
  *   \f[
  *     P^{mn}_{t,l} \leq p^{ac}_{t,l} \leq P^{mx}_{t,l}
  *     \qquad t \in \mathcal{T} , \; l \in \mathcal{L}^{hy}
  *     \tag{14}
  *   \f]
  *   (ActivePower_Bound_Const, boost::multi_array< BoxConstraint , 2 >
  *   indexed [ t ][ l ], group "ActivePower_HydroUnit"). Besides, the sign
  *   of \f$ p^{ac}_{t,l} \f$ is fixed by the Variable itself, nonnegative if
  *   \f$ P^{mn}_{t,l} \geq 0 \f$ and nonpositive if \f$ P^{mx}_{t,l} \leq 0
  *   \f$, and so is that of \f$ f_{t,l} \f$, nonnegative if
  *   \f$ F^{mn}_{t,l} \geq 0 \f$ and nonpositive if
  *   \f$ F^{mx}_{t,l} \leq 0 \f$, while \f$ p^{pr}_{t,l} \f$ and
  *   \f$ p^{sc}_{t,l} \f$ are nonnegative.
  *
  * - Reserves and power (ActivePowerPrimary_Const and
  *   ActivePowerSecondary_Const, both boost::multi_array< FRowConstraint , 2
  *   > indexed [ t ][ l ], groups "ActivePowerPrimary_HydroUnit" and
  *   "ActivePowerSecondary_HydroUnit", each generated only if the
  *   corresponding reserve Variable exist): if arc \f$ l \f$ is a turbine at
  *   instant \f$ t \f$
  *   \f[
  *     p^{pr}_{t,l} \leq \rho^{pr}_{t,l} \, p^{ac}_{t,l} \; , \tag{3}
  *   \f]
  *   \f[
  *     p^{sc}_{t,l} \leq \rho^{sc}_{t,l} \, p^{ac}_{t,l} \; , \tag{4}
  *   \f]
  *   if it is a pump
  *   \f[
  *     p^{pr}_{t,l} = 0 \; , \tag{5}
  *   \f]
  *   \f[
  *     p^{sc}_{t,l} = 0 \; , \tag{6}
  *   \f]
  *   and if it is idle \f$ p^{pr}_{t,l} = f_{t,l} \f$ and
  *   \f$ p^{sc}_{t,l} = f_{t,l} \f$, the flow being zero.
  *
  * - Flow to power (FlowActivePower_Const, a boost::multi_array<
  *   std::vector< FRowConstraint > , 2 > indexed [ t ][ l ], group
  *   "FlowActivePower_HydroUnit", whose entry holds one row for each piece
  *   of arc \f$ l \f$ if it is a turbine at instant \f$ t \f$, and a single
  *   row otherwise): if arc \f$ l \f$ is a turbine at instant \f$ t \f$
  *   \f[
  *     p^{ac}_{t,l} \leq \rho^{hy}_j f_{t,l} + P^{hy}_j
  *     \qquad j \in \mathcal{J}_l \; , \tag{7}
  *   \f]
  *   i.e., \f$ p^{ac}_{t,l} \leq \min_{ j \in \mathcal{J}_l } \{
  *   \rho^{hy}_j f_{t,l} + P^{hy}_j \} \f$, the minimum of affine functions
  *   being a concave piecewise-linear function of the flow. If the pieces
  *   are the tangents \f$ P_j + \rho_j ( f - \bar f_j ) \f$ to a concave
  *   curve at the flows \f$ \bar f_j \f$ (see deserialize()), the function is
  *   the cutting-plane model of the curve, which it bounds from above; at an
  *   optimum the power of a turbine is usually on the curve, since water
  *   turbined without producing power is wasted, and (7) only allows it to
  *   produce less. A single-piece turbine with \f$ \rho^{hy}_j \neq 0 \f$
  *   and \f$ P^{hy}_j = 0 \f$ has the equality \f$ p^{ac}_{t,l} =
  *   \rho^{hy}_j f_{t,l} \f$ in place of (7) (which lets a Solver substitute
  *   the flow away) if a spillway can take the water that the turbine would
  *   pass without producing, with the same effect on the reservoirs, i.e.,
  *   if (i) some other arc \f$ k \f$ is a spillway (a single piece with
  *   \f$ \rho^{hy} = 0 \f$) with \f$ s( k ) = s( l ) \f$,
  *   \f$ e( k ) = e( l ) \f$, \f$ \tau^{up}_k = \tau^{up}_l \f$ and
  *   \f$ \tau^{dn}_k = \tau^{dn}_l \f$ (any spillway if there are no
  *   "StartArc" and "EndArc"), the first such arc being the spillway of
  *   \f$ l \f$, (ii) \f$ F^{mn}_{t,l} \leq 0 \f$ at every instant (i.e.,
  *   0 where it is a turbine) and neither
  *   "DeltaRampUp" nor "DeltaRampDown" is provided, so that nothing else
  *   than its power bounds the flow of the turbine, and (iii) at every
  *   instant \f$ F^{mx}_{t,k} \f$ is at least the sum of the
  *   \f$ F^{mx}_{t,l} \f$ of the turbines whose spillway is \f$ k \f$
  *   (otherwise none of them has the equality). Even then the equality
  *   restricts the feasible set, since the spillway carries its own flow as
  *   well: a solution of (7) is excluded if at some instant the water that
  *   the turbines and their spillway release exceeds the maximum flow of the
  *   spillway plus the water that the turbines need for their power. If arc
  *   \f$ l \f$ is a pump at instant \f$ t \f$
  *   \f[
  *     p^{ac}_{t,l} = \rho^{hy}_j f_{t,l} \; , \tag{8}
  *   \f]
  *   with \f$ j \f$ the only piece of the arc and \f$ P^{hy}_j \f$ ignored
  *   (\f$ \rho^{hy}_j = 1 \f$ if "LinearTerm" is not provided), so that a
  *   pump with \f$ \rho^{hy}_j > 0 \f$ consumes power,
  *   \f$ p^{ac}_{t,l} \leq 0 \f$, in proportion to the water it lifts,
  *   \f$ - f_{t,l} \geq 0 \f$; if the arc is idle,
  *   \f$ p^{ac}_{t,l} = f_{t,l} \f$, both being zero.
  *
  * - Flow bounds (FlowRateBounds_Const, a boost::multi_array<
  *   BoxConstraint , 2 > indexed [ t ][ l ], group
  *   "FlowRateBounds_HydroUnit"):
  *   \f[
  *     F^{mn}_{t,l} \leq f_{t,l} \leq F^{mx}_{t,l}
  *     \qquad t \in \mathcal{T} , \; l \in \mathcal{L}^{hy} \; . \tag{9}
  *   \f]
  *
  * - Ramps (RampUp_Const and RampDown_Const, both boost::multi_array<
  *   FRowConstraint , 2 > indexed [ t ][ l ], groups "RampUp_HydroUnit" and
  *   "RampDown_HydroUnit", each generated only if "DeltaRampUp"
  *   (respectively, "DeltaRampDown") is provided): for
  *   \f$ t \in \mathcal{T} \f$ and \f$ l \in \mathcal{L}^{hy} \f$
  *   \f[
  *     f_{t,l} - f_{t-1,l} \leq \Delta^+_{t,l} \; , \tag{10}
  *   \f]
  *   \f[
  *     f_{t-1,l} - f_{t,l} \leq \Delta^-_{t,l} \; , \tag{11}
  *   \f]
  *   where \f$ f_{-1,l} = F^0_l \f$, i.e., the rows of instant 0 are the
  *   bounds \f$ F^0_l - \Delta^-_{0,l} \leq f_{0,l} \leq F^0_l +
  *   \Delta^+_{0,l} \f$, which a change of the initial flows moves (see
  *   set_initial_flow_rate()). The ramps bound the flow of each arc, pumps
  *   included.
  *
  * - Water balance (FinalVolumeReservoir_Const, a boost::multi_array<
  *   FRowConstraint , 2 > indexed [ t ][ n ], group
  *   "FinalVolumeReservoir_HydroUnit"): for \f$ t \in \mathcal{T} \f$ and
  *   \f$ n \in \mathcal{N}^{hy} \f$
  *   \f[
  *     v^{hy}_{n,t} = v^{hy}_{n,t-1} + A_{n,t}
  *       + \sum_{ l \in \mathcal{L}^{hy} \, : \, e( l ) = n , \;
  *                t - \tau^{dn}_l \in \mathcal{T} } f_{t - \tau^{dn}_l , l}
  *       - \sum_{ l \in \mathcal{L}^{hy} \, : \, s( l ) = n , \;
  *                t - \tau^{up}_l \in \mathcal{T} } f_{t - \tau^{up}_l , l}
  *     \; , \tag{12}
  *   \f]
  *   where \f$ v^{hy}_{n,-1} = V^0_n \f$ if \f$ V^0_n \geq 0 \f$ and
  *   \f$ v^{hy}_{n,-1} = v^{hy}_{n,T-1} \f$ (the cyclic closure) if
  *   \f$ V^0_n < 0 \f$; the row is written with the Variable on the left
  *   and \f$ A_{n,t} \f$, plus \f$ V^0_n \f$ if \f$ t = 0 \f$ and
  *   \f$ V^0_n \geq 0 \f$, on the right. That is, the flow of arc \f$ l \f$
  *   at instant \f$ s \f$ leaves \f$ s( l ) \f$ at instant
  *   \f$ s + \tau^{up}_l \f$ and reaches \f$ e( l ) \f$ at instant
  *   \f$ s + \tau^{dn}_l \f$, and every term counts the water of one instant;
  *   a pump, whose flow is negative, thus takes water from \f$ e( l ) \f$
  *   and gives it to \f$ s( l ) \f$. Without "StartArc" and "EndArc" there
  *   is a single reservoir, which every arc leaves with no delay:
  *   \f$ v^{hy}_{0,t} = v^{hy}_{0,t-1} + A_{0,t} - \sum_{ l } f_{t,l} \f$.
  *
  *   The flows before instant 0 and after instant \f$ T - 1 \f$ are not
  *   part of the model, and a term of (12) is present only if both the
  *   instant of the flow and the instant at which the water leaves or
  *   reaches the reservoir are in \f$ \mathcal{T} \f$. Hence, the water in
  *   transit at the borders of the horizon is not accounted for: with
  *   \f$ \tau^{dn}_l > 0 \f$ the water released before the horizon is not
  *   delivered in its first \f$ \tau^{dn}_l \f$ instants, and the water
  *   released in its last \f$ \tau^{dn}_l \f$ instants never reaches
  *   \f$ e( l ) \f$, so that the final volume of \f$ e( l ) \f$ (hence the
  *   future value of the water, see HydroSystemUnitBlock) does not count
  *   it; with \f$ \tau^{up}_l > 0 \f$ the flows of the last
  *   \f$ \tau^{up}_l \f$ instants are never withdrawn from \f$ s( l ) \f$,
  *   and with \f$ \tau^{up}_l < 0 \f$ neither are those of the first
  *   \f$ - \tau^{up}_l \f$ instants, the water having left the reservoir
  *   before the horizon. The water of a valley is therefore conserved over
  *   the horizon, whatever the flows, only if all its delays are zero.
  *
  *   Summing (12) over the instants gives the final volume
  *   \f[
  *     v^{hy}_{n,T-1} = V^0_n + \sum_{ t \in \mathcal{T} } \Bigl( A_{n,t}
  *       + \sum_{ l : e( l ) = n , \, t - \tau^{dn}_l \in \mathcal{T} }
  *         f_{t - \tau^{dn}_l , l}
  *       - \sum_{ l : s( l ) = n , \, t - \tau^{up}_l \in \mathcal{T} }
  *         f_{t - \tau^{up}_l , l} \Bigr)
  *   \f]
  *   if \f$ V^0_n \geq 0 \f$, an affine function of the initial volumes, of
  *   the flows and of the inflows, in which \f$ V^0_n \f$ appears only
  *   through the right-hand side of the row of instant 0. This is what
  *   allows a multistage model to use the volumes as its state: the initial
  *   volumes of a stage are written in the right-hand sides of
  *   FinalVolumeReservoir_Const[ 0 ] (see set_initial_volume()), and the
  *   dual values of these rows give the cuts of the future cost of the
  *   previous stage (see HydroSystemUnitBlock).
  *
  * - Volume bounds (VolumetricBounds_Const, a boost::multi_array<
  *   BoxConstraint , 2 > indexed [ n ][ t ], group
  *   "VolumetricBounds_HydroUnit"):
  *   \f[
  *     V^{mn}_{n,t} \leq v^{hy}_{n,t} \leq V^{mx}_{n,t}
  *     \qquad n \in \mathcal{N}^{hy} , \; t \in \mathcal{T} \; , \tag{13}
  *   \f]
  *   besides \f$ v^{hy}_{n,t} \geq 0 \f$, which the Variable states.
  *
  * - Reference schedule (Reference_Schedule_Const, a std::vector<
  *   FRowConstraint > of size \f$ 2 T \f$, group
  *   "Norm1_H_Reference_Schedule", generated only if "ReferenceSchedule" is
  *   provided): for \f$ t \in \mathcal{T} \f$, rows \f$ t \f$ and
  *   \f$ T + t \f$ are
  *   \f[
  *     \sum_{ l \in \mathcal{L}^{hy} } p^{ac}_{t,l} - \delta^{rs}_t \leq
  *     \hat p_t \; , \qquad
  *     - \sum_{ l \in \mathcal{L}^{hy} } p^{ac}_{t,l} - \delta^{rs}_t \leq
  *     - \hat p_t \; , \tag{15}
  *   \f]
  *   so that \f$ \delta^{rs}_t \geq | \sum_l p^{ac}_{t,l} - \hat p_t | \f$,
  *   with \f$ \delta^{rs}_t \geq 0 \f$ the deviation Variable of instant
  *   \f$ t \f$ (see generate_objective()).
  *
  * - Reactive power (ReactivePower_Bound_Const, a boost::multi_array<
  *   BoxConstraint , 2 > indexed [ l ][ t ], group "ReactivePowerBound",
  *   generated only if the reactive power Variable exist and some bound on
  *   them is provided): the bounds "MinReactivePower" and
  *   "MaxReactivePower" on the reactive power of each arc.
  *
  * The inertia that the arcs give, \f$ h^p_{t,l} p^{ac}_{t,l} \f$, is a term
  * of the inertia constraints of the UCBlock, not a constraint of this
  * Block. */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the HydroUnitBlock
 /** Method that generates the objective of the HydroUnitBlock, a
  * FRealObjective with a LinearFunction to be minimized:
  * \f[
  *   \min \; \sum_{ t \in \mathcal{T} } \Bigl( \sum_{ l \in
  *     \mathcal{L}^{hy} } b_l \, p^{ac}_{t,l} + \delta^{rs}_t \Bigr) \; ,
  * \f]
  * where \f$ b_l \f$ is "ActivePowerCost" (the term is absent if it is not
  * provided) and the deviation \f$ \delta^{rs}_t \geq 0 \f$ from the reference
  * schedule \f$ \hat p_t \f$, bounded by (15) of
  * generate_abstract_constraints(), is present only if "ReferenceSchedule"
  * is provided. The value of the energy, of the reserves and of the inertia
  * that the valley gives is not part of this objective: it comes from the
  * linking constraints of the UCBlock, whose dual values are the prices that
  * a Lagrangian relaxation of those constraints puts on \f$ p^{ac} \f$,
  * \f$ p^{pr} \f$ and \f$ p^{sc} \f$ (see \ref ucbm_dual), and the value of
  * the water left at the end of the horizon is the polyhedral function of
  * the enclosing HydroSystemUnitBlock. */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------- Methods for checking the HydroUnitBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the HydroUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this HydroUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this HydroUnitBlock is not
  *      greater than the tolerance.
  *
  * Every Constraint of this HydroUnitBlock is a RowConstraint and its
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
  * - Otherwise, by default, the tolerance is 0 and the relative violation
  *   is considered.
  *
  * This function considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed (with
  * Constraint::compute()).
  *
  * @param useabstract This parameter is ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE HydroUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of hydro units
 * @{ */

 /// returns the number of reservoirs
 Index get_number_reservoirs( void ) const {
  return( f_NumberReservoirs ? f_NumberReservoirs : 1 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of arcs / generators

 Index get_number_generators( void ) const override {
  return( f_NumberArcs ? f_NumberArcs : 1 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start arcs
 /** Returns the vector of the start reservoirs \f$ s( l ) \f$ of the arcs
  * ("StartArc"), which has size get_number_generators(), or is empty if
  * "StartArc" is not provided, in which case there is a single reservoir,
  * which every arc leaves. */

 const std::vector< Index > & get_start_arc( void ) const {
  return( v_StartArc );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of end arcs
 /** Returns the vector of the end reservoirs \f$ e( l ) \f$ of the arcs
  * ("EndArc"), get_number_reservoirs() meaning that the arc leaves the
  * valley, which has size get_number_generators(), or is empty if "EndArc"
  * is not provided, in which case every arc leaves the valley. */

 const std::vector< Index > & get_end_arc( void ) const {
  return( v_EndArc );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia power values of the given generator
 /** The returned value U = get_inertia_power( generator ) contains the
  * coefficients \f$ h^p_{t,l} \f$ of the active power of the given arc
  * (generator) \f$ l \f$ in the inertia constraints of the UCBlock: U[ t ]
  * is that of instant t, for t = 0 , ... , get_time_horizon() - 1. If
  * "InertiaPower" is not provided, the unit gives no inertia and nullptr is
  * returned. */

 const double * get_inertia_power( Index generator ) const override {
  if( v_InertiaPower.empty() )
   return( nullptr );
  return( v_InertiaPower.data() + generator * f_time_horizon );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of minimum volumetric
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the minimum volumetric of the reservoir n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then the minimum volumes are
  *   0;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_reservoirs() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] gives the minimum
  *   volumetric of reservoir n at time instant t. */

 const boost::multi_array< double , 2 > & get_min_volumetric( void ) const {
  return( v_MinVolumetric );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of maximum volumetric
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the maximum volumetric of the reservoir n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then the maximum volumes are
  *   0, i.e., the reservoirs have no storage;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_reservoirs() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] gives the maximum
  *   volumetric of reservoir n at time instants t. */

 const boost::multi_array< double , 2 > & get_max_volumetric( void ) const {
  return( v_MaxVolumetric );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inflows
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the inflows of the reservoir n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then the inflows are 0;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_reservoirs() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] represents the inflows
  *   of reservoir n at time instant t. */

 const boost::multi_array< double , 2 > & get_inflows( void ) const {
  return( v_inflows );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power of \p generator at time \p t

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( v_MinPower.empty() ? 0  :
	  *( v_MinPower.data() + t * f_NumberArcs + generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power of \p generator at time \p t

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower.empty() ? 0 :
	  *( v_MaxPower.data() + t * f_NumberArcs + generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum reactive power of \p generator at time \p t

 double get_min_reactive_power( Index t , Index generator = 0 )
  const override {
  return( v_MinReactivePower.empty() ? 0 :
	  *( v_MinReactivePower.data() + t * f_NumberArcs + generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of \p generator at time \p t

 double get_max_reactive_power( Index t , Index generator = 0 )
  const override {
  return( v_MaxReactivePower.empty() ? 0 :
	  *( v_MaxReactivePower.data() + t * f_NumberArcs + generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum flow of \p generator at time \p t

 double get_min_flow( Index t , Index generator = 0 ) const {
  return( v_MinFlow.empty() ? 0 :
	  *( v_MinFlow.data() + t * f_NumberArcs + generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum flow of \p generator at time \p t

 double get_max_flow( Index t , Index generator = 0 ) const {
  return( v_MaxFlow.empty() ? 0 :
	  *( v_MaxFlow.data() + t * f_NumberArcs + generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of delta ramp up
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the delta ramp up at each time t associated with unit
  * (arc) i. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no ramping constraints;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_generators() columns, and
  *   M[ t , i ] gives the delta ramp up value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_delta_ramp_up( void ) const {
  return( v_DeltaRampUp );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of delta ramp down
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the delta ramp down at each time t associated with unit
  * (arc) i. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no ramping constraints;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_generators() columns, and
  *   M[ t , i ] gives the delta ramp down value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_delta_ramp_down( void ) const {
  return( v_DeltaRampDown );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of primary rho
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the primary rho at each time t associated with unit
  * (arc) i. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no primary reserve
  *   constraints;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_generators() columns, and
  *   M[ t , i ] gives the primary rho value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_primary_rho( void ) const {
  return( v_PrimaryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of secondary rho
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the secondary rho at each time t associated with unit
  * (arc) i. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no secondary reserve
  *   constraints;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_generators() columns, and
  *   M[ t , i ] gives the secondary rho value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_secondary_rho( void ) const {
  return( v_SecondaryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of number pieces
 /** The returned vector V contains the number of pieces of the
  * flow-to-power function of each arc: either V is empty, and every arc has
  * one piece, or V has size get_number_generators() and V[ l ] is the
  * number of pieces of arc l. */

 const std::vector< Index > & get_number_pieces( void ) const {
  return( v_NumberPieces );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of constant term
 /** The returned vector V contains the constant terms \f$ P^{hy}_j \f$ of
  * the pieces, numbered as described in deserialize(): either V is empty,
  * and every constant term is 0, or V has one entry per piece. */

 const std::vector< double > & get_const_term( void ) const {
  return( v_ConstTerm );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of linear term
 /** The returned vector V contains the linear terms \f$ \rho^{hy}_j \f$ of
  * the pieces, numbered as described in deserialize(): either V has one
  * entry per piece, or V is empty, which is allowed only if no arc is ever
  * a turbine, the power of a pump then being equal to its flow. */

 const std::vector< double > & get_linear_term( void ) const {
  return( v_LinearTerm );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of uphill delay
 /** The returned vector V contains the uphill delays \f$ \tau^{up}_l \f$
  * of the arcs (see deserialize()): either V is empty, and every delay is
  * 0, or V has size get_number_generators() and V[ l ] is the delay of arc
  * l. */

 const std::vector< int > & get_uphill_delay( void ) const {
  return( v_UphillDelay );
  }

/*--------------------------------------------------------------------------*/
 /// returns the uphill delay for the given \p arc
 /** This function returns the uphill delay for the given \p arc.
  *
  * @return The uphill delay for the given \p arc. */

 int get_uphill_delay( Index arc ) const {
  return( v_UphillDelay.empty() ? 0 : v_UphillDelay[ arc ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of downhill delay
 /** The returned vector V contains the downhill delays \f$ \tau^{dn}_l \f$
  * of the arcs (see deserialize()): either V is empty, and every delay is
  * 0, or V has size get_number_generators() and V[ l ] is the delay of arc
  * l. */

 const std::vector< Index > & get_downhill_delay( void ) const {
  return( v_DownhillDelay );
  }

/*--------------------------------------------------------------------------*/
 /// returns the downhill delay for the given \p arc
 /** This function returns the downhill delay for the given \p arc.
  *
  * @return The downhill delay for the given \p arc. */

 Index get_downhill_delay( Index arc ) const {
  return( v_DownhillDelay.empty() ? 0 : v_DownhillDelay[ arc ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the initial volumetric at the given \p reservoir
 /** This method returns the initial volumetric at the given \p reservoir.
  *
  * @param reservoir The index of the reservoir whose initial volumetric is
  *                  desired.
  *
  * @return The initial volume \f$ V^0_n \f$ of the given \p reservoir,
  *         negative under the cyclic closure (see deserialize()). If
  *         "InitialVolumetric" is not provided, then 0.0 is returned. */

 double get_initial_volumetric( Index reservoir ) const {
  assert( reservoir < get_number_reservoirs() );
  return( v_InitialVolumetric.empty() ? 0.0 :
          ( ( v_InitialVolumetric.size() == 1 ) ?
            v_InitialVolumetric.front() :
            v_InitialVolumetric[ reservoir ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the initial flow rate at the given \p arc
 /** This method returns the initial flow rate at the given \p arc.
  *
  * @param arc The index of the arc whose initial flow rate is desired.
  *
  * @return The initial flow rate \f$ F^0_l \f$ of the given \p arc, 0.0
  *         if "InitialFlowRate" is not provided. */

 double get_initial_flow_rate( Index arc ) const {
  assert( arc < get_number_generators() );
  return( v_InitialFlowRate.empty() ? 0.0 :
          ( ( v_InitialFlowRate.size() == 1 ) ? v_InitialFlowRate.front() :
            v_InitialFlowRate[ arc ] ) );
  }

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE HydroUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the HydroUnitBlock
 *
 * These methods allow to read the five groups of Variable that any
 * HydroUnitBlock in principle has (although some may not):
 *
 * - the volumetric variables;
 *
 * - the flow rate variables;
 *
 * - the active power variables;
 *
 * - the primary spinning reserve variables;
 *
 * - the secondary spinning reserve variables;
 *
 * All these five groups of variables are (if not empty)
 * boost::multi_array< ColVariable , 2 >. The volumetric variables have the
 * number reservoirs as the first dimension and time horizon as the second
 * dimension; all other variables have number of arcs (generators) as the
 * first dimension and time horizon as the second one.
 * @{ */

 /// returns the array of volume variables of the given \p reservoir
 /** This method returns the array of ColVariable representing the volumes of
  * the given \p reservoir at all time instants t in {0, ..., time_horizon -
  * 1}.
  *
  * @param reservoir The index of a reservoir (a number between 0 and
  *                  get_number_reservoirs() - 1).
  *
  * @return The array of ColVariable representing the volumes of the given \p
  *         reservoir. */

 ColVariable * get_volumetric( Index reservoir ) {
  if( reservoir < get_number_reservoirs() ) {
   const auto offset = reservoir * f_time_horizon;
   if( offset < v_volumetric.num_elements() )
    return( v_volumetric.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const array of volume variables of the given \p reservoir

 const ColVariable * get_const_volumetric( Index reservoir ) const {
  if( reservoir < get_number_reservoirs() ) {
   const auto offset = reservoir * f_time_horizon;
   if( offset < v_volumetric.num_elements() )
    return( v_volumetric.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the volume of the given reservoir at the given time
 /** Returns a pointer to the ColVariable representing the volume of the given
  * \p reservoir at the given \p time.
  *
  * @param reservoir The index of the reservoir whose volume is desired (a
  *                  number between 0 and get_number_reservoirs() - 1).
  *
  * @param time The time at which the volume is desired (a number between 0
  *             and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the volume of the given
  *         \p reservoir at the given \p time. */

 ColVariable * get_volume( Index reservoir , Index time ) {
  if( reservoir < get_number_reservoirs() && time < f_time_horizon ) {
   const auto offset = reservoir * f_time_horizon + time;
   if( offset < v_volumetric.num_elements() )
    return( v_volumetric.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// the storages of a hydro unit are its reservoirs

 Index get_number_storages( void ) const override {
  return( get_number_reservoirs() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the volumes of the given reservoir [see get_storage_level()]

 ColVariable * get_storage_level( Index storage ) override {
  return( get_volumetric( storage ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the array of active power variables of the given \p generator
 /** This method returns the array of ColVariable representing the active
  * power of the given \p generator at all time instants t in {0, ...,
  * time_horizon - 1}.
  *
  * @param generator The index of a generator (a number between 0 and
  *                  get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the active power of the
  *         given \p generator. */

 ColVariable * get_active_power( Index generator ) override {
  if( generator < get_number_generators() ) {
   const auto offset = generator * f_time_horizon;
   if( offset < v_active_power.num_elements() )
    return( v_active_power.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the array of reactive power variables of the given \p generator
 /** This method returns the array of ColVariable representing the reactive
  * power of the given \p generator at all time instants t in {0, ...,
  * time_horizon - 1}.
  *
  * @param generator The index of a generator (a number between 0 and
  *                  get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the reactive power of the
  *         given \p generator. */

 ColVariable * get_reactive_power( Index generator ) override {
  if( generator < get_number_generators() ) {
   const auto offset = generator * f_time_horizon;
   if( offset < v_reactive_power.num_elements() )
    return( v_reactive_power.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the active power of the given generator at the given time
 /** Returns a pointer to the ColVariable representing the active power of the
  * given \p generator at the given \p time.
  *
  * @param generator The index of the generator whose active power is desired
  *                  (a number between 0 and get_number_generators() - 1).
  *
  * @param time The time at which the active power is desired (a number
  *             between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the active power of the
  *         given \p generator at the given \p time. */

 ColVariable * get_active_power( Index generator , Index time ) {
  if( generator < get_number_generators() && time < f_time_horizon ) {
   const auto offset = generator * f_time_horizon + time;
   if( offset < v_active_power.num_elements() )
    return( v_active_power.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the reactive power of the given generator at the given time
 /** Returns a pointer to the ColVariable representing the reactive power of
  * the given \p generator at the given \p time.
  *
  * @param generator The index of the generator whose reactive power is
  *                  desired (a number between 0 and
  *                  get_number_generators() - 1)
  *
  * @param time The time at which the reactive power is desired (a number
  *             between 0 and get_time_horizon() - 1)
  *
  * @return A pointer to the ColVariable representing the reactive power of
  *         the given \p generator at the given \p time */

 ColVariable * get_reactive_power( Index generator , Index time ) {
  if( generator < get_number_generators() && time < f_time_horizon ) {
   const auto offset = generator * f_time_horizon + time;
   if( offset < v_reactive_power.num_elements() )
    return( v_reactive_power.data() + offset );
   }
  return( nullptr );
  } 

/*--------------------------------------------------------------------------*/
 /// returns the array of flow rate variables along the given \p arc
 /** Returns the array of ColVariable representing the flow rate along the
  * given \p arc at all time instants t in { 0, ..., time_horizon - 1 }.
  *
  * @param arc The index of an arc (a number between 0 and
  *            get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the flow rate along the
  *         given \p arc. */

 ColVariable * get_flow_rate( Index arc ) {
  if( arc < get_number_generators() ) {
   const auto offset = arc * f_time_horizon;
   if( offset < v_flow_rate.num_elements() )
    return( v_flow_rate.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const array of flow rate variables along the given \p arc

 const ColVariable * get_const_flow_rate( Index arc ) const {
  if( arc < get_number_generators() ) {
   const auto offset = arc * f_time_horizon;
   if( offset < v_flow_rate.num_elements() )
    return( v_flow_rate.data() + offset );
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the flow rate along the given arc at the given time
 /** Returns a pointer to the ColVariable representing the flow rate along the
  * given \p arc at the given \p time.
  *
  * @param arc The index of the arc whose flow rate is desired (a number
  *            between 0 and get_number_generators() - 1).
  *
  * @param time The time at which the flow rate is desired (a number
  *             between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the flow rate along the
  *         given \p arc at the given \p time. */

 ColVariable * get_flow_rate( Index arc , Index time ) {
  if( arc < get_number_generators() && time < f_time_horizon ) {
   const auto offset = arc * f_time_horizon + time;
   if( offset < v_flow_rate.num_elements() )
    return( v_flow_rate.data() + offset );
    }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// the unit provides the reserve only if it has any to give
 /** The enclosing UCBlock asking for the reserve is not enough, the unit has
  * to have some to give: this is the very condition with which the Variable
  * are generated, said in terms of the data alone. */

 bool has_primary_reserve( void ) const override {
  return( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the unit provides the secondary reserve only if it has any to give

 bool has_secondary_reserve( void ) const override {
  return( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the array of primary spinning reserve of the given \p generator
 /** This method returns the array of ColVariable representing the primary
  * spinning reserve of the given \p generator at all time instants t in {0,
  * ..., time_horizon - 1}.
  *
  * @param generator The index of a generator (a number between 0 and
  *                  get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the primary spinning
  *         reserve of the given \p generator. */

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( generator < get_number_generators() ) {
   const auto offset = generator * f_time_horizon;
   if( offset < v_primary_spinning_reserve.num_elements() )
    return( v_primary_spinning_reserve.data() + offset );
    }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the primary spinning reserve of a generator at the given time
 /** Returns a pointer to the ColVariable representing the primary spinning
  * reserve of the given \p generator at the given \p time.
  *
  * @param generator The index of the generator whose primary spinning reserve
  *                  is desired (a number between 0 and
  *                  get_number_generators() - 1).
  *
  * @param time The time at which the primary spinning reserve is desired (a
  *             number between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the primary spinning
  *         reserve of the given \p generator at the given \p time. */

 ColVariable * get_primary_spinning_reserve( Index generator , Index time ) {
  if( generator < get_number_generators() && time < f_time_horizon ) {
   const auto offset = generator * f_time_horizon + time;
   if( offset < v_primary_spinning_reserve.num_elements() )
    return( v_primary_spinning_reserve.data() + offset );
    }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the array of secondary spinning reserve of the given \p generator
 /** This method returns the array of ColVariable representing the secondary
  * spinning reserve of the given \p generator at all time instants t in {0,
  * ..., time_horizon - 1}.
  *
  * @param generator The index of a generator (a number between 0 and
  *        get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the secondary spinning
  *         reserve of the given \p generator. */

 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( generator < get_number_generators() ) {
   const auto offset = generator * f_time_horizon;
   if( offset < v_secondary_spinning_reserve.num_elements() )
    return( v_secondary_spinning_reserve.data() + offset );
    }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the secondary spinning reserve of a generator at the given time
 /** Returns a pointer to the ColVariable representing the secondary spinning
  * reserve of the given \p generator at the given \p time.
  *
  * @param generator The index of the generator whose secondary spinning
  *        reserve is desired (a number between 0 and get_number_generators()
  *        - 1).
  *
  * @param time The time at which the secondary spinning reserve is desired (a
  *        number between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the secondary spinning
  *         reserve of the given \p generator at the given \p time. */

 ColVariable * get_secondary_spinning_reserve( Index generator , Index time )
 {
  if( generator < get_number_generators() && time < f_time_horizon ) {
   const auto offset = generator * f_time_horizon + time;
   if( offset < v_secondary_spinning_reserve.num_elements() )
    return( v_secondary_spinning_reserve.data() + offset );
    }
  return( nullptr );
  }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution storing the current solution of this HydroUnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this HydroUnitBlock. This
  * is a HydroUnitBlockSolution extending UnitBlockSolution with the
  * specific extra solution information of HydroUnitBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - the first four bits (bit 0 to bit 3) are "taken" by the base
  *   UnitBlock[Solution]
  *
  * - bit 4 (& 16) means "store the volumetric levels"
  *
  * - bit 5 (& 32) means "store the flow rates"
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr, and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 63 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [Hydro]UnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE HydroUnitBlock------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the HydroUnitBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * HydroUnitBlock. See HydroUnitBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE HydroUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the HydroUnitBlock
 * @{ */

 /// loading from a stream is not implemented: it throws

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "HydroUnitBlock::load() not implemented yet" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for changing the data of the HydroUnitBlock
 *  @{ */

 /// set the inflow values
 /** This function sets the inflow values \f$ A_{n,t} \f$ of this
  * HydroUnitBlock, i.e., the right-hand sides of the water balances (12)
  * [see generate_abstract_constraints()], \p values being in the order of
  * the entries [ n ][ t ] (reservoir after reservoir). The right-hand side
  * of instant 0 is \f$ V^0_n + A_{n,0} \f$, and \f$ A_{n,0} \f$ alone for a
  * reservoir under the cyclic closure (\f$ V^0_n < 0 \f$).
  *
  * @param values  Iterator to a vector containing the inflow values.
  * @param subset  If non-empty, the inflow values corresponding to the
  *                indices in \p subset are set to the values pointed by
  *                \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_inflow( MF_dbl_it values ,
                  Subset && subset , bool ordered = false ,
                  c_ModParam issuePMod = eNoBlck ,
                  c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the inflow values
 /** This function sets the inflow values \f$ A_{n,t} \f$ of this
  * HydroUnitBlock, i.e., the right-hand sides of the water balances (12)
  * [see generate_abstract_constraints()], \p values being in the order of
  * the entries [ n ][ t ] (reservoir after reservoir). The right-hand side
  * of instant 0 is \f$ V^0_n + A_{n,0} \f$, and \f$ A_{n,0} \f$ alone for a
  * reservoir under the cyclic closure (\f$ V^0_n < 0 \f$).
  *
  * @param values Iterator to a vector containing the inflow values.
  * @param rng    If non-empty, the inflow values corresponding to the
  *               indices in \p rng are set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_inflow( MF_dbl_it values ,
                  Range rng = Range( 0 , Inf< Index >() ) ,
                  c_ModParam issuePMod = eNoBlck ,
                  c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the inertia power values
 /** This function sets the inertia power values \f$ h^p_{t,l} \f$ of this
  * HydroUnitBlock, \p values being in the order of the entries [ l ][ t ]
  * of "InertiaPower" (arc after arc). No Constraint of this Block uses
  * them: the enclosing UCBlock, when it sees the HydroUnitBlockMod of type
  * eSetInerP, changes the coefficients of the inertia constraints, which
  * however cannot gain a term that they did not have when they were
  * generated, i.e., setting an inertia power for a unit that had none then
  * makes the UCBlock throw std::logic_error.
  *
  * @param values  Iterator to a vector containing the inertia power values.
  * @param subset  If non-empty, the inertia power values corresponding to the
  *                indices in \p subset are set to the values pointed by
  *                \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_inertia_power( MF_dbl_it values ,
                         Subset && subset , bool ordered = false ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the inertia power values
 /** This function sets the inertia power values \f$ h^p_{t,l} \f$ of this
  * HydroUnitBlock, \p values being in the order of the entries [ l ][ t ]
  * of "InertiaPower" (arc after arc). No Constraint of this Block uses
  * them: the enclosing UCBlock, when it sees the HydroUnitBlockMod of type
  * eSetInerP, changes the coefficients of the inertia constraints, which
  * however cannot gain a term that they did not have when they were
  * generated, i.e., setting an inertia power for a unit that had none then
  * makes the UCBlock throw std::logic_error.
  *
  * @param values Iterator to a vector containing the inertia power values.
  * @param rng    If non-empty, the inertia power values corresponding to the
  *               indices in \p rng are set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_inertia_power( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the initial volume
 /** This function sets the initial volumes \f$ V^0_n \f$ of this
  * HydroUnitBlock, i.e., the right-hand sides of the water balances (12) of
  * instant 0 [see generate_abstract_constraints()]. A negative value means
  * the cyclic closure, which shapes those rows when they are generated:
  * once they are, a value whose sign differs from that of the current one
  * makes the method throw std::invalid_argument.
  *
  * @param values  Iterator to a vector containing the initial volume values.
  * @param subset  If non-empty, the initial volume corresponding to the
  *                indices in \p subset is set to the values pointed by
  *                \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_initial_volume( MF_dbl_it values ,
                          Subset && subset , bool ordered = false ,
                          c_ModParam issuePMod = eNoBlck ,
                          c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the initial volume
 /** This function sets the initial volumes \f$ V^0_n \f$ of this
  * HydroUnitBlock, i.e., the right-hand sides of the water balances (12) of
  * instant 0 [see generate_abstract_constraints()]. A negative value means
  * the cyclic closure, which shapes those rows when they are generated:
  * once they are, a value whose sign differs from that of the current one
  * makes the method throw std::invalid_argument.
  *
  * @param values Iterator to a vector containing the initial volume values.
  * @param rng    If non-empty, the initial volume corresponding to the
  *               indices in \p rng is set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_initial_volume( MF_dbl_it values ,
                          Range rng = Range( 0 , Inf< Index >() ) ,
                          c_ModParam issuePMod = eNoBlck ,
                          c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the initial flow rate
 /** This function sets the initial flows \f$ F^0_l \f$ of this
  * HydroUnitBlock, i.e., the sides of the ramp rows (10) and (11) of
  * instant 0 [see generate_abstract_constraints()], which change if they
  * exist. If "InitialFlowRate" is absent, it is first given one entry per
  * arc, all 0.
  *
  * @param values  Iterator to a vector containing the initial flow rate
  *                values.
  * @param subset  If non-empty, the initial flow rate corresponding to the
  *                indices in \p subset is set to the values pointed by
  *                \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_initial_flow_rate( MF_dbl_it values ,
                             Subset && subset , bool ordered = false ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the initial flow rate
 /** This function sets the initial flows \f$ F^0_l \f$ of this
  * HydroUnitBlock, i.e., the sides of the ramp rows (10) and (11) of
  * instant 0 [see generate_abstract_constraints()], which change if they
  * exist. If "InitialFlowRate" is absent, it is first given one entry per
  * arc, all 0.
  *
  * @param values Iterator to a vector containing the initial flow rate values.
  * @param rng    If non-empty, the initial flow rate corresponding to the
  *               indices in \p rng is set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_initial_flow_rate( MF_dbl_it values ,
                             Range rng = Range( 0 , Inf< Index >() ) ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the active power cost values
 /** This function sets the active power cost values of this HydroUnitBlock.
  * If the Objective is generated without "ActivePowerCost", and hence
  * without terms in the active power, the terms of the arcs whose cost is
  * set are added to it.
  *
  * @param values  Iterator to a vector containing the active power cost
  *                values.
  * @param subset  If non-empty, the active power cost values corresponding
  *                to the indices in \p subset are set to the values pointed
  *                by \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_active_power_cost( MF_dbl_it values ,
                             Subset && subset ,
                             bool ordered = false ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the active power cost values
 /** This function sets the active power cost values of this HydroUnitBlock.
  * If the Objective is generated without "ActivePowerCost", and hence
  * without terms in the active power, the terms of the arcs whose cost is
  * set are added to it.
  *
  * @param values Iterator to a vector containing the active power cost
  *               values.
  * @param rng    If non-empty, the active power cost values corresponding
  *               to the indices in \p rng are set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_active_power_cost( MF_dbl_it values ,
                             Range rng = Range( 0 , Inf< Index >() ) ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the same active power cost for every arc
 /** This function sets the active power cost \f$ b_l \f$ of every arc
  * \f$ l \f$ of this HydroUnitBlock to \p value.
  *
  * @param value     The value of the active power cost.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_active_power_cost( double value ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector( get_number_generators() , value );
  set_active_power_cost( vector.cbegin() ,
                         Range( 0 , Inf< Index >() ) ,
                         issuePMod , issueAMod );
 }

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// verify whether the data in this HydroUnitBlock is consistent
 /** This function checks whether the data in this HydroUnitBlock is
  * consistent, with the notation of deserialize(), and throws an exception
  * if it is not. The data are consistent if all the following conditions
  * are met.
  *
  * - "StartArc" and "EndArc" are both provided or both not provided, and
  *   they are provided if there is more than one reservoir; if they are,
  *   \f$ s( l ) < R \f$, \f$ e( l ) \leq R \f$ and \f$ s( l ) \neq e( l ) \f$
  *   for every arc \f$ l \f$.
  *
  * - \f$ P^{mn}_{t,l} \leq P^{mx}_{t,l} \f$ and
  *   \f$ F^{mn}_{t,l} \leq F^{mx}_{t,l} \f$ for every arc \f$ l \f$ and
  *   instant \f$ t \f$, if both bounds are provided.
  *
  * - An arc \f$ l \f$ that is a pump at instant \f$ t \f$ (i.e.,
  *   \f$ F^{mn}_{t,l} < 0 \f$ and \f$ F^{mx}_{t,l} \leq 0 \f$) has one
  *   piece and \f$ \rho^{pr}_{t,l} = \rho^{sc}_{t,l} = 0 \f$.
  *
  * - No arc has \f$ F^{mn}_{t,l} < 0 < F^{mx}_{t,l} \f$ at some instant,
  *   or is a turbine at some instant and a pump at another, and "LinearTerm"
  *   is provided if some arc is a turbine at some instant.
  *
  * - \f$ 0 \leq V^{mn}_{n,t} \leq V^{mx}_{n,t} \f$ for every reservoir
  *   \f$ n \f$ and instant \f$ t \f$, if both bounds are provided.
  *
  * - For every reservoir \f$ n \f$ under the cyclic closure
  *   (\f$ V^0_n < 0 \f$), the total inflow \f$ \sum_t A_{n,t} \f$ (0 without
  *   "Inflows") lies between the smallest and the largest total net outflow
  *   that the flow bounds allow, up to a tolerance of
  *   \f$ 10^{-9} \max \{ 1 , | \sum_t A_{n,t} | \} \f$, as explained in
  *   deserialize(); this is necessary for the balances of the reservoir to
  *   have a solution. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the number of reservoirs (nodes) of the problem
 Index f_NumberReservoirs{};

 /// the number of arcs connecting the reservoirs in the cascading system
 Index f_NumberArcs{};

 /// the total number of pieces
 Index f_TotalNumberPieces{};

 /// the vector of UphillDelay
 std::vector< int > v_UphillDelay;

 /// the vector of DownhillDelay
 std::vector< Index > v_DownhillDelay;

 /// the vector of starting arc
 std::vector< Index > v_StartArc;

 /// the vector of ending arcs
 std::vector< Index > v_EndArc;

 /// the vector of initial volumetric
 std::vector< double > v_InitialVolumetric;

 /// the vector of initial flow rate
 std::vector< double > v_InitialFlowRate;

 /// the vector of NumberPieces
 std::vector< Index > v_NumberPieces;

 /// the vector of LinearTerm
 std::vector< double > v_LinearTerm;

 /// the vector of ConstTerm
 std::vector< double > v_ConstTerm;

 /// the vector of ActivePowerCost
 std::vector< double > v_ActivePowerCost;

 /// the matrix of inertia power of generators
 /** Indexed over the dimensions NumberArcs and TimeHorizon. */
 boost::multi_array< double , 2 > v_InertiaPower;

 /// the matrix of MinVolumetric
 /** Indexed over the dimensions NumberReservoirs and TimeHorizon. */
 boost::multi_array< double , 2 > v_MinVolumetric;

 /// the matrix of MaxVolumetric
 /** Indexed over the dimensions NumberReservoirs and TimeHorizon. */
 boost::multi_array< double , 2 > v_MaxVolumetric;

 /// the matrix of Inflows
 /** Indexed over the dimensions NumberReservoirs and TimeHorizon. */
 boost::multi_array< double , 2 > v_inflows;

 /// the matrix of MinPower
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_MinPower;

 /// the matrix of MaxPower
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_MaxPower;

 /// the matrix of MinReactivePower
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_MinReactivePower;

 /// the matrix of MaxReactivePower
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_MaxReactivePower;

 /// the matrix of MinFlow
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_MinFlow;

 /// the matrix of MaxFlow
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_MaxFlow;

 /// the matrix of DeltaRampUp
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_DeltaRampUp;

 /// the matrix of DeltaRampDown
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_DeltaRampDown;

 /// the matrix of PrimaryRho
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_PrimaryRho;

 /// the matrix of SecondaryRho
 /** Indexed over the dimensions TimeHorizon and NumberArcs. */
 boost::multi_array< double , 2 > v_SecondaryRho;

 /// the reference schedule to deviate minimally from if there
 std::vector< double > v_RefSchedule;

/*-------------------------------- variables -------------------------------*/

 /// the matrix of volumetric variables
 boost::multi_array< ColVariable , 2 > v_volumetric;

 /// the matrix of flow rate variables
 boost::multi_array< ColVariable , 2 > v_flow_rate;

 /// the active power variables
 boost::multi_array< ColVariable , 2 > v_active_power;

 /// the reactive power variables
 boost::multi_array< ColVariable , 2 > v_reactive_power;

 /// the primary spinning reserve variables
 boost::multi_array< ColVariable , 2 > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 boost::multi_array< ColVariable , 2 > v_secondary_spinning_reserve;

  /// the variables for deviation to reference schedule
 std::vector< ColVariable > v_abs_ref_schedule;

/*------------------------------- constraints ------------------------------*/

 /// the reference schedule constraints
 std::vector< FRowConstraint > Reference_Schedule_Const;

 /// maximum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint , 2 > MaxPowerPrimarySecondary_Const;

 /// minimum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint , 2 > MinPowerPrimarySecondary_Const;

 /// power output relation with to primary reserves constraints
 boost::multi_array< FRowConstraint , 2 > ActivePowerPrimary_Const;

 /// power output relation with to secondary reserves constraints
 boost::multi_array< FRowConstraint , 2 > ActivePowerSecondary_Const;

 /// flow to active power function constraints
 /** FlowActivePower_Const[ t ][ l ] has a row for each piece of arc l if it
  * is a turbine at time t, and a single row otherwise. */
 boost::multi_array< std::vector< FRowConstraint > , 2 > FlowActivePower_Const;

 /// not generated, the bounds on the power being ActivePower_Bound_Const
 boost::multi_array< FRowConstraint , 2 > ActivePowerBounds_Const;

 /// active power bounds when the unit produces no reserve
 boost::multi_array< BoxConstraint , 2 > ActivePower_Bound_Const;

 /// the reactive power bound constraints
 boost::multi_array< BoxConstraint , 2 > ReactivePower_Bound_Const;

 /// ramp-up constraints
 boost::multi_array< FRowConstraint , 2 > RampUp_Const;

 /// ramp-down constraints
 boost::multi_array< FRowConstraint , 2 > RampDown_Const;

 /// the water balance constraints (12), indexed [ t ][ n ]
 boost::multi_array< FRowConstraint , 2 > FinalVolumeReservoir_Const;

 /// flow rate bounds constraints
 boost::multi_array< BoxConstraint , 2 > FlowRateBounds_Const;

 /// volumetric bounds constraints
 boost::multi_array< BoxConstraint , 2 > VolumetricBounds_Const;

 /*!! Q <= P
 boost::multi_array< FRowConstraint , 2 > Reactive_2_Active_Const;
 !!*/

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
 /// updates the constraints for the given arcs at time 0
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints associated with the given
  * \p arcs at time 0 (which are the ramp constraints that depend on the
  * initial flow rate).
  *
  * @param arcs The indices of the arcs whose constraints must be updated
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void update_initial_flow_rate_in_cnstrs( Range arcs ,
                                          c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the given arcs at time 0
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints associated with the given
  * \p arcs at time 0 (which are the ramp constraints that depend on the
  * initial flow rate).
  *
  * @param arcs The indices of the arcs whose associated constraints must
  *             be updated
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void update_initial_flow_rate_in_cnstrs( const Block::Subset & arcs ,
                                          c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 static void static_initialization( void )
 {
  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_inflow" , & HydroUnitBlock::set_inflow );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_inflow" , & HydroUnitBlock::set_inflow );

  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_inertia_power" ,
   & HydroUnitBlock::set_inertia_power );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_inertia_power" ,
   & HydroUnitBlock::set_inertia_power );

  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_initial_volume" ,
   & HydroUnitBlock::set_initial_volume );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_initial_volume" ,
   & HydroUnitBlock::set_initial_volume );

  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_initial_flow_rate" ,
   & HydroUnitBlock::set_initial_flow_rate );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_initial_flow_rate" ,
   & HydroUnitBlock::set_initial_flow_rate );

  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_active_power_cost" ,
   & HydroUnitBlock::set_active_power_cost );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_active_power_cost" ,
   & HydroUnitBlock::set_active_power_cost );
 }

/*--------------------------------------------------------------------------*/

 };  // end( class( HydroUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS HydroUnitBlockMod -------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a HydroUnitBlock
class HydroUnitBlockMod : public UnitBlockMod
{
 public:

 /// public enum for the types of HydroUnitBlockMod
 enum HUB_mod_type {
  eSetInf = eUBModLastParam , ///< set inflow values
  eSetInerP ,                 ///< set inertia power values
  eSetInitF ,                 ///< set initial flow rate values
  eSetInitV ,                 ///< set initial volumetric values
  eSetActPCost ,              ///< set active power cost values
  eHUBModLastParam            ///< first allowed parameter for derived classes
  /**< Convenience value to easily allow derived classes to extend the set
   * of types of HydroUnitBlockMod. */
  };

 /// constructor, takes the HydroUnitBlock and the type
 HydroUnitBlockMod( HydroUnitBlock * fblock , int type )
  : UnitBlockMod( fblock , type ) , f_Block( fblock ) {}

 /// destructor, does nothing
 virtual ~HydroUnitBlockMod( void ) override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the HydroUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "HydroUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetInf ):
    output << "set inflow values";
    break;
   case( eSetInerP ):
    output << "set inertia power values";
    break;
   case( eSetInitF ):
    output << "set initial flow rate values";
    break;
   case( eSetInitV ):
    output << "set initial volumetric values";
    break;
   case( eSetActPCost ):
    output << "set active power cost values";
    break;
   default:;
   }
  }

 HydroUnitBlock * f_Block;
 ///< pointer to the Block to which the Modification refers

 };  // end( class( HydroUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS HydroUnitBlockRngdMod -----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from HydroUnitBlockMod for "ranged" modifications
class HydroUnitBlockRngdMod : public HydroUnitBlockMod
{

 public:

 /// constructor: takes the HydroUnitBlock, the type, and the range
 HydroUnitBlockRngdMod( HydroUnitBlock * fblock , int type ,
                        Block::Range rng )
  : HydroUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~HydroUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the HydroUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  HydroUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

 };  // end( class( HydroUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS HydroUnitBlockSbstMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from HydroUnitBlockMod for "subset" modifications
class HydroUnitBlockSbstMod : public HydroUnitBlockMod
{
 public:

 /// constructor: takes the HydroUnitBlock, the type, and the subset
 HydroUnitBlockSbstMod( HydroUnitBlock * fblock , int type ,
                        Block::Subset && nms )
  : HydroUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~HydroUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the HydroUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  HydroUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

 Block::Subset f_nms;  ///< the subset

 };  // end( class( HydroUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS HydroUnitBlockSolution -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [UnitBlock]Solution of a HydroUnitBlock
/** The HydroUnitBlockSolution class derives from UnitBlockSolution and
 * adds the "standard" information stored in there (active power, possibly
 * commitment and primary/secondary reserve) the other information that is
 * typical of the HydroUnitBlock, i.e.,
 *
 * - [possibly] the volumetric level of each reservoir at each time instant
 *
 * - [possibly] the volumetric flow on each arc at each time instant */

class HydroUnitBlockSolution : public UnitBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend HydroUnitBlock;  ///< make HydroUnitBlock friend

/*--------- CONSTRUCTING AND DESTRUCTING HydroUnitBlockSolution ----------*/

 /// constructor, it has nothing to do
 explicit HydroUnitBlockSolution( void )
  : UnitBlockSolution() , f_reservoirs( 0 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~HydroUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*----- METHODS DESCRIBING THE BEHAVIOR OF A HydroUnitBlockSolution -----*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a HydroUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a HydroUnitBlockSolution into a netCDF::NcGroup. The format
  * is the one of UnitBlockSolution [cf. UnitBlockSolution::serialize()],
  * plus:
  *
  * - The dimension "NumberReservoirs" containing the number of reservoirs
  *   (nodes) in the HydroUnitBlock. The dimension is optional, if it is not
  *   provided then it is taken to be == 1. Note, however, that a single
  *   reservoir can still have multiple hydro generating units (arc). Yet,
  *   there is no need of a dimension for the number of arcs since "number
  *   of arcs == f_number_generators" and this is already a field of
  *   UnitBlockSolution.
  *
  * - The variable "VolumetricLevel", of type netCDF::NcDouble and indexed
  *   over both the dimension "TimeHorizon" (that is mandatory for the base
  *   UnitBlockSolution) and "NumberReservoirs"; VolumetricLevel[ i ][ t ]
  *   is the optimal value of the volumetric level of the reservoir i at
  *   time t. The variable is optional.
  *
  * - The variable "VolumetricFlow", of type netCDF::NcDouble and indexed
  *   over both the dimension "TimeHorizon" (that is mandatory for the base
  *   UnitBlockSolution) and "NumberGenerators" (that is optional for the
  *   base UnitBlockSolution and provides the number of generators, which
  *   in hydro units correspond to arcs): VolumetricFlow[ i ][ t ] is the
  *   amount of water flowing down arc i at time t. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 HydroUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 HydroUnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "HydroUnitBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 Index f_reservoirs;   ///< the number of reservoirs

 boost::multi_array< double , 2 > v_volume;
 ///< v_volume[ i ][ t ] = volume of reservoir i at time t

 boost::multi_array< double , 2 > v_flow;
 ///< v_flow[ i ][ t ] = flow down arc i at time t

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( HydroUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __HydroUnitBlock */

/*--------------------------------------------------------------------------*/
/*---------------------- End File HydroUnitBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
