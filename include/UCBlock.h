/*--------------------------------------------------------------------------*/
/*--------------------------- File UCBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class UCBlock, which implements the Block concept [see
 * Block.h] for the Unit Commitment (UC) problem in electrical power
 * production. This is typically a short-term (across, for instance, one
 * week or one-day time horizon) deterministic problem regarding finding
 * an optimal production schedule of electrical generators satisfying a
 * (large) set of technical constraints. The model as a whole is described
 * in \ref ucblock_model.
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Kostas Tavlaridis-Gyparakis,
 *                      Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __UCBlock
 #define __UCBlock  ///< self-identification: endif at the end of the file

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "UnitBlock.h"

#include "NetworkBlock.h"

#include "FRowConstraint.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS UCBlock --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Unit Commitment problem
/** UCBlock implements the Block concept [see Block.h] for the Unit Commitment
 * (UC) problem in electrical power production. This is typically a short-term
 * (across, for instance, one week or one-day time horizon) deterministic
 * problem regarding finding an optimal schedule of the production of
 * electrical generators satisfying a (large) set of technical constraints.
 *
 * The model is quite flexible, since the class manages son Block of type
 * UnitBlock and NetworkBlock, and therefore one can use different types of
 * units and network constraints. Also, UCBlock handles a quite large variety
 * of constraints, regarding not only active power but also primary and
 * secondary reserve and inertia.
 *
 * In particular, UCBlock handles the following main elements:
 *
 * - A time horizon, i.e., a discrete set of (typically, equally-spaced) time
 *   instants at which decisions are made (like, the 24 hours in a day).
 *
 * - A set of electricity generating units, represented by derived classes
 *   of the base class UnitBlock.
 *
 * - A set of NetworkBlock that altogether cover (or may be created to cover)
 *   the entire time horizon which represents the constraints on the
 *   electricity demand satisfaction and the technical constraints on the
 *   network. These can be basically "empty" if the capacity of the
 *   transmission network is such as to never really impact generation
 *   decisions (a "bus").
 *
 * - Constraints linking the production decisions at the units and ensuring:
 *
 *   - balance between production of active power and injection in the
 *     transmission network, at each node and for each time instant;
 *
 *   - possibly, primary and secondary reserve constraints for each "zone"
 *     (appropriately defined subset of the nodes of the transmission
 *     network) and for each time instant;
 *
 *   - possibly, constraints about inertia  for each "zone" (appropriately
 *     defined subset of the nodes of the transmission network) and for
 *     each time instant;
 *
 *   - possibly, constraints on the maximum (and minimum) emission of
 *     different kinds of pollutant, for each "zone" (appropriately defined
 *     subset of the nodes of the transmission network) and across the
 *     whole time horizon.
 *
 * Note that UCBlock has no Variable, and its Objective is a constant (see
 * generate_objective()); hence, the cost is that of its sub-Blocks, and its
 * only constraints are the linking ones (see
 * generate_abstract_constraints()), which is the structure that a
 * LagrangianDualSolver exploits. The capacity of the assets enters the model
 * through the scale of the UnitBlock (see UnitBlock::scale()), which
 * represents a number of identical copies of a unit, and through the capacity
 * multipliers of IntermittentUnitBlock, BatteryUnitBlock and DCNetworkBlock.
 * Since either can be changed after the Block is built, one can evaluate the
 * optimal cost, together with a subgradient, as a function of the capacities;
 * alternatively, the same capacities can be made variables of the model
 * through the design variables of the units and of DesignNetworkBlock. In the
 * page \ref ucblock_model we describe the complete model, the conventions it
 * follows (instants, units of measure, signs), its Lagrangian dual and the
 * meaning of the dual values of the linking constraints, the representation
 * of the capacities and what is not modeled.
 */

class UCBlock : public Block
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

 /// constructor of UCBlock, taking possibly a pointer of its father Block

 explicit UCBlock( Block * father = nullptr )
  : Block( father ) , f_time_horizon( 0 ) , f_has_reactive( false ) ,
    f_number_networks( 0 ) , f_number_units( 0 ) ,
    f_number_elc_generators( 0 ) , f_total_number_pollutant_zones( 0 ) ,
    f_number_primary_zones( 0 ) , f_number_secondary_zones( 0 ) ,
    f_number_inertia_zones( 0 ) , f_number_pollutants( 0 ) ,
    f_number_storages( 0 ) ,
    f_NetworkData( nullptr ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of UCBlock

 virtual ~UCBlock() override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the UCBlock. Besides the mandatory "type" attribute of any :Block, the
  * group should contain the following:
  *
  * - The dimension "TimeHorizon" containing the number of time steps in the
  *   problem.
  *
  * - The dimension "NumberUnits" containing the number of units (UnitBlock)
  *   in the problem.
  *
  * - The dimension "NumberNetworks" containing the number of networks
  *   (NetworkBlock) in the problem. This dimension is optional. If it is not
  *   provided, it is taken to be the same as "TimeHorizon".
  *
  * - The dimension "NumberNodes" containing the number of nodes in the
  *   problem. This dimension is optional. If it is not provided, it is taken
  *   to be 1, i.e., the network is a single node (a bus), in which case the
  *   UCBlock has no NetworkBlock (see "NetworkBlock_i" below).
  *
  * - The dimension "NumberElectricalGenerators" containing the overall number
  *   of electrical generators for each UnitBlock of the problem. This
  *   dimension is optional. If it is not provided, it is taken to be the sum
  *   of the electrical generators of each UnitBlock, so, if each of these
  *   will have just one electrical generator, i.e., no HydroUnitBlock or
  *   HydroSystemUnitBlock will be provided, it is taken to be the same as
  *   "NumberUnits". If it is provided, it must be that sum, otherwise
  *   std::invalid_argument is thrown.
  *
  * - The groups "UnitBlock_0", "UnitBlock_1", ..., "UnitBlock_n" with
  *   n == NumberUnits - 1, containing each one UnitBlock corresponding to
  *   one or more electricity generating units (electrical generators). It is
  *   an error if the corresponding groups are not there. Each UnitBlock can
  *   have more than one electrical generator (cf.
  *   UnitBlock::get_number_generators()); a value that is useful in the
  *   following (cf. "GeneratorNode") is the total number of those, which we
  *   will refer to as "NumberElectricalGenerators". This is computed by just
  *   calling get_number_generators() on each of the UnitBlock and summing all
  *   the results. Clearly, NumberElectricalGenerators >= NumberUnits. Indeed,
  *   most of the UnitBlock can be expected to have just one electrical
  *   generator. If this happens for all the units, then
  *   NumberElectricalGenerators == NumberUnits. If, instead, some UnitBlock
  *   (like cascades of hydro generators or combined cycle plants) actually
  *   has more than one electrical generator, then NumberElectricalGenerators
  *   > NumberUnits (each UnitBlock must have at least one). It is then useful
  *   (cf. "GeneratorNode") to be able to assign a unique index g = 0, 1, ...,
  *   NumberElectricalGenerators - 1 to each of the electrical generators in
  *   the UCBlock. When NumberElectricalGenerators == NumberUnits, the index
  *   is the same as i = 0, 1, ..., NumberUnits - 1 (there is a one-to-one
  *   correspondence between UnitBlock and electrical generators). When,
  *   instead, NumberElectricalGenerators > NumberUnits, a mapping must be
  *   defined. The mapping is the obvious one: UnitBlock has an ordering
  *   i = 0, 1, ..., NumberUnits - 1 (cf. the groups "UnitBlock_0",
  *   "UnitBlock_1", ... above), and the electrical generators into each
  *   UnitBlock also have some natural ordering (corresponding to the columns
  *   of the matrices of variables, cf. e.g., UnitBlock::get_commitment()).
  *   Thus, in general, with k = UnitBlock_0->get_number_generators(), the
  *   mapping is:
  *
  *   - electrical generator 0 is the first generator of UnitBlock_0;
  *
  *   - electrical generator 1 is the second generator of UnitBlock_0;
  *
  *   - electrical generator k - 1 is the k-th (last) generator of
  *     UnitBlock_0;
  *
  *   - electrical generator k is the first generator of UnitBlock_1;
  *
  *   - electrical generator k + 1 is the second generator of UnitBlock_1;
  *
  *   and so on, which of course boils down to "g = i" when each UnitBlock
  *   has exactly one electrical generator.
  *
  * - Optionally, the dimensions and variables necessary to deserialize a
  *   NetworkData object that describes the transmission network; see
  *   NetworkBlock::NetworkData::deserialize() for details. If that is not
  *   provided (basically, "NumberNodes" is not provided, or it is == 1),
  *   then the transmission network is taken to have only one node (a bus).
  *   If the data of a NetworkData is specified, the NetworkData is passed
  *   to each of the NetworkBlock (see below) of the UCBlock, if any. However,
  *   if a NetworkBlock also has a NetworkData specified in its own group,
  *   then the NetworkData inside the NetworkBlock overrules that inside the
  *   UCBlock, which is ignored by that NetworkBlock.
  *
  * - The variable "ActivePowerDemand", of type netCDF::NcDouble and indexed
  *   both over the dimensions "NumberNodes" and "TimeHorizon". This variable
  *   is optional if:
  *
  *   - a NetworkBlock is defined for all necessary intervals so that they
  *     collectively cover the entire "TimeHorizon" (see "NetworkBlock_i"
  *     below), and
  *
  *   - each defined NetworkBlock has the "ActiveDemand" variable specified
  *     in the corresponding group.
  *
  *   Otherwise, it is mandatory. When it is defined, the entry
  *   ActivePowerDemand[ n , t ] is assumed to contain the active power demand
  *   \f$ D^{ac}_{t,n} \f$ of node n of the transmission network at the
  *   given time instant t (a power, see \ref ucbm_conv), where the first
  *   dimension "NumberNodes" can be read via get_number_nodes(). The demand
  *   of a node is meant to include the reference consumption profile of any
  *   load-shifting flexibility connected to it, whose deviation from that
  *   profile is the active power of the BatteryUnitBlock that represents it.
  *   When "ActivePowerDemand" is defined, and also the "ActiveDemand"
  *   variable is defined in the group of some NetworkBlock, then
  *   "ActiveDemand" overrules the value in the corresponding row of
  *   "ActivePowerDemand", which is ignored.
  *
  * - The variable "ReactivePowerDemand", of type netCDF::NcDouble and
  *   indexed as "ActivePowerDemand", the reactive power demand
  *   \f$ D^{rc}_{t,n} \f$, with the same rules as "ActivePowerDemand" and
  *   "ReactiveDemand" of the NetworkBlock; it is only used if the network
  *   handles reactive power (see NetworkBlock::handles_reactive()), or, in
  *   the bus case, if it is given, in which case the units are asked for
  *   their reactive power and the reactive balance (1') of
  *   generate_abstract_constraints() is written.
  *
  * - The groups "NetworkBlock_0", "NetworkBlock_1", ..., "NetworkBlock_T"
  *   with T = NumberNetworks - 1 (or a subset of these, see below), each
  *   containing the network constraints for some subset of time intervals in
  *   the horizon. These can be handled in two different ways:
  *
  *   - If "NumberNodes" is == 1 (e.g., it is not provided), then the network
  *     is a "bus", i.e., there are no constraints on how the data flows ("no
  *     network") and all the units can contribute to satisfy the energy
  *     demand, which is then a unique value for each time instant t. In this
  *     case one would assume that the (total) demand for each time instant is
  *     provided by "ActivePowerDemand", and this would be enough. However,
  *     there is an extra mechanism: if some "NetworkBlock_i" is still
  *     provided (see below for details on how they are treated), then they
  *     are read. This is done in order to get the corresponding matrix of
  *     active demands (see NetworkBlock::get_active_demand()): if it is
  *     present (which it may not), then the corresponding values are used as
  *     active power demand for these intervals, superseding the values found
  *     in "ActivePowerDemand". It is therefore possible that then
  *     "ActivePowerDemand" is not there at all, but this requires that
  *     "NetworkBlock_i" are given that cover the whole time horizon (this
  *     is otherwise not necessary, see below) so that active power demand
  *     data is specified for each time instant; if this fails to happen,
  *     exception is thrown. The NetworkBlock read are then deleted, and the
  *     UCBlock has none (see get_network_blocks()); their constant terms, if
  *     any, are kept in the Objective of the UCBlock (see
  *     "NetworkConstantTerms" and generate_objective()).
  *
  *   - If, instead, "NumberNodes" > 1, then there are constraints (and,
  *     possibly, variables) limiting how energy flows between producing
  *     units and demand. Then, the sub-group "NetworkBlock_0",
  *     "NetworkBlock_1", ..., "NetworkBlock_T" with T = NumberNetworks - 1
  *     are attempted to be read. They are assumed to be arranged in
  *     chronological order according to their index. Each NetworkBlock
  *     specifies the number of consecutive time intervals it spans, see
  *     NetworkBlock::get_number_intervals(). Thus, "NetworkBlock_0" covers
  *     the intervals 0, ..., t - 1 with
  *     t = "NetworkBlock_0"->get_number_intervals(), "NetworkBlock_1" covers
  *     the intervals t, ..., t + w - 1 with
  *     w = "NetworkBlock_1"->get_number_intervals(), and so on. The last
  *     "NetworkBlock_T" has to cover the final intervals up to
  *     "TimeHorizon"; otherwise, an exception is raised. However, it is
  *     possible that some (or even, in principle, all) "NetworkBlock_i" is
  *     not specified, which is useful when (as it often happens, since the
  *     network may easily not change in the short time horizon typical of
  *     the UC problem) some (or all) of them are "equal". If any
  *     "NetworkBlock_i" is missing, a NetworkBlock is automatically built by
  *     using the "global" NetworkData in UCBlock. The corresponding value of
  *     NetworkBlock::get_number_intervals() is used to identify the subset of
  *     time instants it covers (starting from the initial time t identified as
  *     previously specified), and the "ActivePowerDemand" data for these
  *     intervals t is used to set its demand. Hence, an exception is raised
  *     if "ActivePowerDemand" or "NetworkData" are missing. In this way it
  *     is possible, e.g., to just specify "ActivePowerDemand" and one
  *     "NetworkData" and have UCBlock to automatically construct all the
  *     necessary NetworkBlock (which are then identical copies of one another
  *     save possibly for the active power demand). A typical case is that in
  *     which "NumberNetworks" is equal to "TimeHorizon", i.e., each
  *     "NetworkBlock_t" covers exactly one interval. In this case it is clear
  *     which time instant each "NetworkBlock_t" covers (i.e., "t"), and
  *     obviously it must be NetworkBlock::get_number_intervals() == 1 (the
  *     default).
  *
  * - The variable "GeneratorNode", of type netCDF::NcUint and indexed over
  *   the set { 0 , ... , NumberElectricalGenerators - 1 }; GeneratorNode[ g ]
  *   tells to which node of the transmission network the specified electrical
  *   generator g belongs. Note that this means that different electrical
  *   generators in the same UnitBlock can belong to different nodes of the
  *   transmission network. This is justified, e.g., by hydro cascade units
  *   where different turbines can be rather far apart geographically, but
  *   still linked by (long) stretches of rivers. If NumberElectricalGenerators
  *   == NumberUnits (all UnitBlock have exactly one electrical generator),
  *   then this variable is indexed over NumberUnits. If NumberNodes == 1
  *   (say, it is not provided at all), then this variable need not be
  *   defined, since it is not used.
  *
  * - The dimension "NumberPrimaryZones" tells how many "primary spinning
  *   reserve zones" are there in the problem. The dimension is optional, if it
  *   is not provided then it is taken to be 0, which means that no primary
  *   reserve constraints are present in the problem.
  *
  * - The variable "PrimaryZones", of type netCDF::NcUint and indexed over
  *   the dimension "NumberNodes" (a scalar meaning the same value for all
  *   the nodes). The entry PrimaryZones[ n ] tells to which primary zone the
  *   node n belongs: if PrimaryZones[ n ] >= NumberPrimaryZones, this means
  *   that node n does not belong to any primary zone, and hence the
  *   corresponding electrical generators are not involved in the primary
  *   reserve constraints; this holds whatever the number of zones, one
  *   included. If NumberPrimaryZones == 0 (say, it is not provided at all)
  *   then this variable need not be defined, since it is not used. If
  *   NumberPrimaryZones == 1 and this variable is not defined, then all the
  *   nodes belong to the only primary zone; if NumberPrimaryZones > 1 the
  *   variable is mandatory, and an exception is thrown if it is not there.
  *
  * - The variable "PrimaryDemand", of type netCDF::NcDouble and indexed both
  *   over the dimensions "NumberPrimaryZones" and "TimeHorizon": entry
  *   PrimaryDemand[ z , t ] is assumed to contain the primary reserve
  *   requirement \f$ D^{pr}_{t,z} \f$ of the primary reserve zone z at the
  *   time instant t. If
  *   NumberPrimaryZones == 0 (say, it is not provided at all), then this
  *   variable need not be defined, since it is not loaded.
  *
  * - The dimension "NumberSecondaryZones" tells how many "secondary spinning
  *   reserve zones" are there in the problem. The dimension is optional, if it
  *   is not provided then it is taken to be 0, which means that no secondary
  *   reserve constraints are present in the problem.
  *
  * - The variable "SecondaryZones", of type netCDF::NcUint and indexed over
  *   the dimension "NumberNodes", with the same meaning and rules as
  *   "PrimaryZones" for the secondary zones: SecondaryZones[ n ] >=
  *   NumberSecondaryZones means that node n belongs to no secondary zone,
  *   the variable is not needed with at most one zone (all the nodes being
  *   then in it), and it is mandatory with more than one.
  *
  * - The variable "SecondaryDemand", of type netCDF::NcDouble and indexed
  *   both over the dimensions "NumberSecondaryZones" and "TimeHorizon": entry
  *   SecondaryDemand[ z , t ] is assumed to contain the secondary reserve
  *   requirement \f$ D^{sc}_{t,z} \f$ of the secondary reserve zone z at
  *   the time instant t.
  *   If NumberSecondaryZones == 0 (say, it is not provided at all), then
  *   this variable need not be defined, since it is not loaded.
  *
  * - The dimension "NumberInertiaZones" tells how many "inertia constraints
  *   zones" are there in the problem. The dimension is optional, if it is not
  *   provided then it is taken to be 0, which means that no inertia
  *   constraints are present in the problem.
  *
  * - The variable "InertiaZones", of type netCDF::NcUint and indexed over
  *   the dimension "NumberNodes", with the same meaning and rules as
  *   "PrimaryZones" for the inertia zones: InertiaZones[ n ] >=
  *   NumberInertiaZones means that node n belongs to no inertia zone, the
  *   variable is not needed with at most one zone (all the nodes being then
  *   in it), and it is mandatory with more than one.
  *
  * - The variable "InertiaDemand", of type netCDF::NcDouble and indexed both
  *   over the dimensions "NumberInertiaZones" and "TimeHorizon": entry
  *   InertiaDemand[ z , t ] is assumed to contain the inertia requirement
  *   \f$ D^{in}_{t,z} \f$ of the inertia zone z at the time instant t.
  *   If NumberInertiaZones == 0 (say, it is not provided at all), then this
  *   variable need not be defined, since it is not loaded.
  *
  * - The dimension "NumberPollutants" containing the number of pollutants in
  *   the problem. The dimension is optional, if it is not provided then it is
  *   taken to be 0, which means that no pollutant constraints are present in
  *   the problem.
  *
  * - The variable "NumberPollutantZones" of type netCDF::NcUint and indexed
  *   over the dimension "NumberPollutants": the entry
  *   NumberPollutantZones[ p ] is assumed to contain the number of pollutant
  *   zones associated with pollutant p. If NumberPollutants == 0 (say, it is
  *   not provided) then this variable need not be defined, since it is not
  *   loaded. If NumberPollutants > 0 and this variable is not defined, then
  *   each pollutant has exactly one zone. The total number of pollutant zones
  *   is useful (cf. PollutantBudget); it will be referred to as
  *   "TotalNumberPollutantZones", and it is computed simply as
  *   TotalNumberPollutantZones = NumberPollutantZones[ 0 ] + ... +
  *   NumberPollutantZones[ NumberPollutants - 1 ]. A dimension with this name
  *   may be present (typically, to index PollutantBudget), in which case it
  *   is an error if its size is not that sum.
  *
  * - The variable "PollutantZones", of type netCDF::NcUint and indexed over
  *   the dimensions "NumberPollutants" and "NumberNodes": the entry
  *   PollutantZones[ p , n ] tells to which pollutant zone associated with
  *   pollutant p the node n belongs. If PollutantZones[ p , n ] >=
  *   NumberPollutantZones[ p ], this means that node n does not belong to any
  *   pollutant zone, and hence the corresponding units are not involved in the
  *   pollutant budget constraints associated with pollutant p. The first
  *   dimension can have size 1, in which case the same zones are used for
  *   all pollutants. If NumberPollutants == 0 (say, it is not provided) then
  *   this variable need not be defined, since it is not loaded. If each
  *   pollutant has exactly one zone and this variable is not defined, then
  *   all the nodes belong to the unique zone of each pollutant.
  *
  * - The variable "PollutantBudget", of type netCDF::NcDouble and indexed
  *   over the set { 0 , ... , TotalNumberPollutantZones - 1 }: the entry
  *   PollutantBudget[ n ] for n == 0, ..., TotalNumberPollutantZones - 1 is
  *   assumed to contain the pollutant budget (across all the time horizon)
  *   for the pair (pollutant zone, pollutant) corresponding to n. In other
  *   words, since the number of pollutant zones of each pollutant may differ,
  *   it is useful (to avoid having to store PollutantBudget as an irregular
  *   matrix) to be able to assign a unique index
  *   n = 0, 1, ..., TotalNumberPollutantZones - 1 to each pollutant budget
  *   of each pollutant zone. A mapping must be defined between each entry n
  *   and the pair (pollutant zone, pollutant). The mapping is the obvious
  *   one: UCBlock has a set of pollutants p = 0, 1, ..., NumberPollutants -
  *   1, and each pollutant p may have several pollutant zones (see comments
  *   of variable "NumberPollutantZones" above). Thus, in general, with
  *   k = NumberPollutantZones[ 0 ], the mapping is:
  *
  *   - n = 0 corresponds to zone 0 of pollutant 0;
  *
  *   - n = 1 corresponds to zone 1 of pollutant 0;
  *
  *   - n = k - 1 corresponds to zone k - 1 (the last) of pollutant 0;
  *
  *   - n = k corresponds to zone 0 of pollutant 1;
  *
  *   - n = k + 1 corresponds to zone 1 of pollutant 1;
  *
  *   and so on. The budget of zone z of pollutant p is denoted
  *   \f$ O_{z,p} \f$. If NumberPollutants == 0 (say, it is not provided)
  *   then this variable need not be defined, since it is not loaded;
  *   otherwise it is mandatory.
  *
  * - The variable "PollutantRho", of type netCDF::NcDouble and indexed over
  *   three dimensions which are "TimeHorizon" and "NumberPollutants" and
  *   the set { 0, ..., NumberElectricalGenerators - 1 } (see comments above).
  *   The first dimension can have size either 1 or "TimeHorizon". In the
  *   former case the entry PollutantRho[ 0 , p , g ] is assumed to contain
  *   the conversion factor of pollutant p due to the electrical generator g
  *   which is equal for all time instants t. Otherwise, the first dimension
  *   has full size "TimeHorizon" and the entry PollutantRho[ t , p , g ]
  *   gives the conversion factor of pollutant p due to the electrical
  *   generator g for time t. The conversion factor \f$ \rho_{t,p,g} \f$
  *   multiplies the active power of the generator, hence it also accounts
  *   for the duration of the time instant (it is, say, in tonnes per MW over
  *   one time instant, see \ref ucbm_conv). If
  *   NumberPollutants == 0 (it is not provided) then this variable need not
  *   be defined, since it's not loaded; otherwise it is mandatory.
  *
  * - The variable "PollutantMinBudget", of type netCDF::NcDouble and indexed
  *   over the set { 0 , ... , TotalNumberPollutantZones - 1 } as
  *   PollutantBudget: the entry PollutantMinBudget[ n ] is the lower bound
  *   \f$ O^{mn}_{z,p} \f$ on
  *   the same emission that PollutantBudget[ n ] bounds from above (hence a
  *   limit that the emission must reach, or with PollutantMinBudget[ n ] ==
  *   PollutantBudget[ n ] that it must match). The variable is optional, if
  *   it is not provided then there is no lower bound, i.e., all the entries
  *   are -INF.
  *
  * - The dimension "NumberStorages", the number of storages of all the units
  *   (see UnitBlock::get_number_storages()), those of unit 0 first, then
  *   those of unit 1, and so on. The dimension is optional, since the number
  *   is known from the units, and if it is provided then it is an error if
  *   it is a different one.
  *
  * - The variable "PollutantStorageRho", of type netCDF::NcDouble and indexed
  *   over three dimensions which are "TimeHorizon", "NumberPollutants" and
  *   "NumberStorages", with the first that can have size either 1 or
  *   "TimeHorizon" as in PollutantRho. The entry PollutantStorageRho[ t , p ,
  *   k ] is the factor \f$ \rho^{v}_{t,p,k} \f$ of the level of storage k
  *   at the end of time t in the
  *   pollutant budget constraints of pollutant p (with t == 0 for all t if
  *   the first dimension has size 1), which is how the change of the level
  *   of a storage that holds a pollutant (or a fuel that emits one) is
  *   accounted for, typically with a nonzero factor at the last time instant
  *   only. The storages of a unit belong to the node of the first electrical
  *   generator of the unit. The variable is optional, if it is not provided
  *   then all the factors are 0; if NumberPollutants == 0 it is not loaded.
  *
  * - The variable "NetworkConstantTerms", of type netCDF::NcDouble and
  *   indexed over the dimension "NumberNetworks"; the entry
  *   NetworkConstantTerms[ k ] tells the constant term \f$ c^0_k \f$,
  *   i.e., typically the fixed costs, of the NetworkBlock k. The variable is
  *   optional, and the "ConstantTerm" of a "NetworkBlock_k" group, if any,
  *   prevails over NetworkConstantTerms[ k ], which applies to the
  *   NetworkBlock whose group does not give one (and to those that the
  *   UCBlock builds); without either, the constant term is 0. In the bus case
  *   the constant terms, read in the same way, are kept in the Objective of
  *   the UCBlock (see generate_objective()), since the NetworkBlock are not
  *   built, and serialize() writes them as "NetworkConstantTerms"; with more
  *   than one node they are in the NetworkBlock, whose groups serialize()
  *   writes with them.
  *
  * - The variable "NetworkBlockClassname", of type netCDF::NcString specifies
  *   the classname of the specific NetworkBlock to be instantiated, if no one
  *   is explicitly given in input; the default value is "DCNetworkBlock".
  *
  * - The variable "NetworkDataClassname", of type netCDF::NcString specifies
  *   the classname of the specific NetworkData to be instantiated, if no one
  *   is explicitly given in input; the default value is "DCNetworkData". */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /// extends Block::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends UnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generates the static constraints of the UCBlock
 /** This method generates the abstract constraints of the sub-Blocks and
  * then those of the UCBlock, i.e., the constraints that link the
  * sub-Blocks with each other; the model they belong to, with the notation
  * used here, is described in \ref ucblock_model.
  *
  * Consider a network defined by a set of nodes \f$ \mathcal{N} \f$ and a
  * set of lines \f$ \mathcal{L} \f$ connecting the nodes. Three families of
  * zones are given, each zone being a set of nodes and the zones of a family
  * being pairwise disjoint: \f$ \mathcal{Z}^{pr} \f$ for the primary
  * spinning reserve, \f$ \mathcal{Z}^{sc} \f$ for the secondary spinning
  * reserve and \f$ \mathcal{Z}^{in} \f$ for inertia, each zone \f$ z \f$
  * carrying one requirement per instant. The families may coincide, and
  * they need not cover \f$ \mathcal{N} \f$: a node may belong to no zone of
  * a family (see "PrimaryZones" in deserialize()), in which case its
  * generators do not enter the corresponding constraints, and when every
  * node belongs to some zone the family is a partition of
  * \f$ \mathcal{N} \f$. In the same way, for each pollutant
  * \f$ p \in \mathcal{P} \f$ a family \f$ \mathcal{Z}^{p} \f$ of disjoint
  * zones carries the emission budgets. Active power needs no zones, its
  * demand being given node by node.
  *
  * The electrical system contains a set of "units" (e.g., power plants,
  * load flexibilities or storage devices), the UnitBlock
  * \f$ i \in \mathcal{I} \f$. Each unit contains one (or more) electrical
  * generators, the set of all electrical generators being
  * \f$ \mathcal{G} \f$, and \f$ \mathcal{G}_n \f$ the set of those
  * connected to node \f$ n \in \mathcal{N} \f$. Let \f$ \sigma_g \f$ be the
  * scale of the UnitBlock of generator \f$ g \f$ (see UnitBlock::scale()),
  * which multiplies every term of the generator below, since the Variable
  * of a scaled unit are those of one copy of it.
  *
  * The UCBlock does not have any Variable of its own, but it defines
  * Constraint on the Variable of the UnitBlock and of the NetworkBlock, as
  * follows:
  *
  * - \f$ p^{ac}_{t,g} \f$: the active power of generator
  *   \f$ g \in \mathcal{G} \f$ at instant \f$ t \in \mathcal{T} \f$,
  *   obtained from UnitBlock::get_active_power();
  *
  * - \f$ p^{rc}_{t,g} \f$: the reactive power of generator \f$ g \f$ at
  *   instant \f$ t \f$, obtained from UnitBlock::get_reactive_power() (may
  *   not be defined);
  *
  * - \f$ S_{t,n} \f$: the active power node injection of node
  *   \f$ n \in \mathcal{N} \f$ at instant \f$ t \f$, obtained from
  *   NetworkBlock::get_node_injection();
  *
  * - \f$ R_{t,n} \f$: the reactive power node injection of node \f$ n \f$
  *   at instant \f$ t \f$, obtained from
  *   NetworkBlock::get_reactive_node_injection() (may not be defined);
  *
  * - \f$ p^{pr}_{t,g} \f$: the primary spinning reserve of generator
  *   \f$ g \f$ at instant \f$ t \f$, obtained from
  *   UnitBlock::get_primary_spinning_reserve() (may not be defined);
  *
  * - \f$ p^{sc}_{t,g} \f$: the secondary spinning reserve of generator
  *   \f$ g \f$ at instant \f$ t \f$, obtained from
  *   UnitBlock::get_secondary_spinning_reserve() (may not be defined);
  *
  * - \f$ u_{t,g} \in \{ 0 , 1 \} \f$: the commitment of generator \f$ g \f$
  *   at instant \f$ t \f$, obtained from UnitBlock::get_commitment() (may
  *   not be defined, in which case the generator is always on);
  *
  * - \f$ v^{sl}_{t,k} \f$: the level of storage \f$ k \f$ at the end of
  *   instant \f$ t \f$, obtained from UnitBlock::get_storage_level() (may
  *   not be defined).
  *
  * The linking constraints of the unit commitment problem are the following;
  * each is written with all the Variable on the left-hand side, which fixes
  * the sign of its dual value (see get_Solution()).
  *
  * - Active power node injection constraints. Let \f$ P^{au}_{t,g} \geq 0
  *   \f$ be the power that generator \f$ g \f$ draws when it is off (see
  *   UnitBlock::get_fixed_consumption()), counted only if the generator has
  *   commitment variables. If get_number_nodes() > 1, for each
  *   \f$ t \in \mathcal{T} \f$ and \f$ n = 0 , \ldots , N - 1 \f$
  *   \f[
  *     \sum_{ g \in \mathcal{G}_n } \sigma_g \bigl( p^{ac}_{t,g}
  *       - P^{au}_{t,g} ( 1 - u_{t,g} ) \bigr) - S_{t,n} = 0 \; , \tag{1}
  *   \f]
  *   written as \f$ \sum_g \sigma_g ( p^{ac}_{t,g} + P^{au}_{t,g} u_{t,g} )
  *   - S_{t,n} = \sum_g \sigma_g P^{au}_{t,g} \f$; \f$ S_{t,n} \f$ is the
  *   node injection of the NetworkBlock that covers instant \f$ t \f$, whose
  *   node balances, with the demand \f$ D^{ac}_{t,n} \f$, complete the
  *   balance of active power (see \ref ucbm_link_ac). If the network has a
  *   single node there is no NetworkBlock, and the row is
  *   \f[
  *     \sum_{ g \in \mathcal{G} } \sigma_g \bigl( p^{ac}_{t,g}
  *       - P^{au}_{t,g} ( 1 - u_{t,g} ) \bigr) = D^{ac}_{t} \; , \tag{1b}
  *   \f]
  *   \f$ D^{ac}_{t} \f$ being "ActivePowerDemand". The power of a pump or of
  *   a charging battery, being negative, enters as a consumption. The rows
  *   are stored in a boost::multi_array< FRowConstraint , 2 > indexed by
  *   [ t ][ n ] (see get_node_injection_constraints()).
  *
  * - Reactive power node injection constraints. If the network has more
  *   than one node and its NetworkBlock handle reactive power (see
  *   NetworkBlock::handles_reactive(); either all of them do or none does,
  *   and only the first one is asked), then for each
  *   \f$ t \in \mathcal{T} \f$ and \f$ n = 0 , \ldots , N - 1 \f$
  *   \f[
  *     \sum_{ g \in \mathcal{G}_n } \sigma_g \, p^{rc}_{t,g} - R_{t,n} = 0
  *     \; , \tag{1'}
  *   \f]
  *   while in the bus case, if "ReactivePowerDemand" is given, the single
  *   row of each instant is \f$ \sum_{ g \in \mathcal{G} } \sigma_g \,
  *   p^{rc}_{t,g} = D^{rc}_t \f$. The fixed consumption is an active power,
  *   and it does not appear in these rows. They are stored as those of (1)
  *   (see get_reactive_node_injection_constraints()).
  *
  * - Primary spinning reserve constraints: if get_number_primary_zones()
  *   > 0, for \f$ t \in \mathcal{T} \f$ and \f$ z \in \mathcal{Z}^{pr} \f$
  *   (the zones being numbered \f$ z = 0 , \ldots , \f$
  *   get_number_primary_zones() \f$ - 1 \f$)
  *   \f[
  *     \sum_{ n \in z } \sum_{ g \in \mathcal{G}_n } \sigma_g \,
  *       p^{pr}_{t,g} \geq D^{pr}_{t,z} \; , \tag{2}
  *   \f]
  *   where a generator whose UnitBlock has no primary reserve variables
  *   does not appear; the rows are stored in a
  *   boost::multi_array< FRowConstraint , 2 > indexed by [ t ][ z ] (see
  *   get_primary_demand_constraints()).
  *
  * - Secondary spinning reserve constraints: if
  *   get_number_secondary_zones() > 0, for \f$ t \in \mathcal{T} \f$ and
  *   \f$ z \in \mathcal{Z}^{sc} \f$
  *   \f[
  *     \sum_{ n \in z } \sum_{ g \in \mathcal{G}_n } \sigma_g \,
  *       p^{sc}_{t,g} \geq D^{sc}_{t,z} \; , \tag{3}
  *   \f]
  *   stored as (2) (see get_secondary_demand_constraints()).
  *
  * - Inertia constraints. The inertia that generator \f$ g \f$ provides at
  *   instant \f$ t \f$ is the affine function \f$ h^u_{t,g} u_{t,g} +
  *   h^p_{t,g} p^{ac}_{t,g} \f$, with \f$ h^u_{t,g} \f$ given by
  *   UnitBlock::get_inertia_commitment() and \f$ h^p_{t,g} \f$ by
  *   UnitBlock::get_inertia_power(), each term being absent if the
  *   corresponding datum or Variable is. For a synchronous machine of
  *   inertia constant \f$ H_t \f$ and rated power \f$ \hat P^{mx}_t \f$
  *   the contribution is usually \f$ 1.2 H_t \hat P^{mx}_t \f$ whenever it
  *   is on, i.e., \f$ h^u_{t,g} = 1.2 H_t \hat P^{mx}_t \f$ and
  *   \f$ h^p_{t,g} = 0 \f$,
  *   while for a turbine or an intermittent unit it is taken proportional
  *   to the power produced, i.e., \f$ h^u_{t,g} = 0 \f$ and
  *   \f$ h^p_{t,g} = 1.2 H_t \f$; the factor 1.2 converts active into
  *   apparent power under a constant phase angle and belongs to the data,
  *   which the UCBlock uses as they are. A BatteryUnitBlock provides no
  *   inertia. If get_number_inertia_zones() > 0, for
  *   \f$ t \in \mathcal{T} \f$ and \f$ z \in \mathcal{Z}^{in} \f$
  *   \f[
  *     \sum_{ n \in z } \sum_{ g \in \mathcal{G}_n } \sigma_g \bigl(
  *       h^u_{t,g} u_{t,g} + h^p_{t,g} p^{ac}_{t,g} \bigr)
  *       \geq D^{in}_{t,z} \; , \tag{4}
  *   \f]
  *   stored as (2) (see get_inertia_demand_constraints()).
  *
  * - Pollutant budget constraints. For each pollutant
  *   \f$ p \in \mathcal{P} \f$ and each zone \f$ z \in \mathcal{Z}^{p} \f$
  *   of it, the emission of pollutant \f$ p \f$ in \f$ z \f$, summed over
  *   the whole time horizon, is bounded from above by the budget
  *   \f$ O_{z,p} \f$ and from below by \f$ O^{mn}_{z,p} \f$ (\f$ - \infty
  *   \f$ if there is no lower bound, and equal to \f$ O_{z,p} \f$ for a
  *   limit the emission has to match):
  *   \f[
  *     O^{mn}_{z,p} \leq \sum_{ t \in \mathcal{T} } \sum_{ n \in z }
  *       \Bigl( \sum_{ g \in \mathcal{G}_n } \sigma_g \, \rho_{t,p,g}
  *       \, p^{ac}_{t,g} + \sum_{ k \in \mathcal{K}_n } \sigma_k \,
  *       \rho^{v}_{t,p,k} \, v^{sl}_{t,k} \Bigr) \leq O_{z,p} \; . \tag{5}
  *   \f]
  *   The conversion factor \f$ \rho_{t,p,g} \f$ (see get_pollutant_rho())
  *   is the emission of pollutant \f$ p \f$ per unit of active power of
  *   generator \f$ g \f$ over instant \f$ t \f$, i.e., the emission per unit
  *   of energy times the duration \f$ \Delta t \f$ of the instant (and
  *   divided by the efficiency of the generator, if the emission is that of
  *   the fuel it burns). The second sum is on the storages
  *   \f$ \mathcal{K}_n \f$ of the units at node \f$ n \f$ (see
  *   UnitBlock::get_number_storages(); the storages of a unit are at the
  *   node of its first electrical generator), \f$ \sigma_k \f$ being the
  *   scale of the unit of storage \f$ k \f$ and \f$ \rho^{v}_{t,p,k} \f$ its
  *   factor (see get_pollutant_storage_rho()). This is how a limit accounts
  *   for the change of the level of a storage over the horizon, as that of
  *   a store of CO2 or of a fuel that emits it, typically with a nonzero
  *   factor at the last time instant only; the contribution of its initial
  *   level, which is a constant, is then to be subtracted from both bounds
  *   by whoever writes the data. The nodes of a zone need not be electric:
  *   the generators on a heat node (a boiler, the heat of a
  *   ConversionUnitBlock) have their own factors, so that a budget over a
  *   zone that holds the heat nodes counts the emissions due to heat. The
  *   rows are stored in a std::vector< std::vector< FRowConstraint > > C,
  *   where C[ p ] has get_number_pollutant_zones()[ p ] entries (the
  *   pollutants may have a
  *   different number of zones, and no row is wasted on the zones one of
  *   them does not have), C[ p ][ z ] being (5) (see
  *   get_pollutant_constraints()). A row with
  *   \f$ - \infty < O^{mn}_{z,p} < O_{z,p} < \infty \f$ cannot be relaxed
  *   by a LagrangianDualSolver.
  *
  * When the scale of a unit changes after the constraints are generated,
  * the coefficients of its terms, and the right-hand sides of (1) and (1b),
  * are changed accordingly, those of (1') included (see
  * add_Modification()). */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// generates the objective of the UCBlock
 /** This method generates the Objective of all the sub-Blocks, and then
  * that of the UCBlock itself, which is a LinearFunction with no Variable
  * whose constant term \f$ c^{U} \f$ is the sum of the constant terms of the
  * NetworkBlock that are not built because the network is a bus (see
  * "NetworkConstantTerms" in deserialize()), and zero otherwise. As for any
  * Block, the objective of the problem is the sum of the Objective of the
  * UCBlock and of those of all its sub-Blocks, hence the problem minimizes
  * \f[
  *   \sum_{ i \in \mathcal{I} } f_i + \sum_{ k } \Bigl( c^0_k
  *     + \sum_{ l \in \mathcal{L} } c^{net}_l V_{k,l} \Bigr) + c^{U} \; ,
  * \f]
  * where \f$ f_i \f$ is the cost of UnitBlock \f$ i \f$, already
  * multiplied by its scale (see UnitBlock::scale()) and, for a
  * HydroSystemUnitBlock, including the future cost of the water left in its
  * reservoirs (see HydroSystemUnitBlock), and the second sum runs over the
  * NetworkBlock \f$ k \f$, \f$ c^0_k \f$ being its constant term and
  * \f$ c^{net}_l V_{k,l} \f$ the optional cost of the flow on line
  * \f$ l \f$ (see DCNetworkBlock). All the costs are per instant, no length
  * of the time step appearing in the model (see \ref ucbm_conv). The
  * Configuration \p objc is not used. */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Checking the UCBlock and handling its Solution
 *  @{ */

 /// checks whether the current solution is feasible for the UCBlock
 /** The solution is feasible if every sub-Block is, and the violation of
  * each Constraint of the UCBlock is not greater than the tolerance. The
  * tolerance and the type of violation are taken from \p fsbc if it is a
  * SimpleConfiguration< double > (tolerance, relative violation) or a
  * SimpleConfiguration< std::pair< double , int > > (tolerance, relative
  * violation if the second is nonzero), otherwise from
  * f_BlockConfig->f_is_feasible_Configuration in the same way, otherwise
  * they are 1e-6 and the relative violation. Each sub-Block is checked with
  * the same tolerance and type of violation, unless its BlockConfig has its
  * own f_is_feasible_Configuration, which is then the one used. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


/*--------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/
/*
 Block * get_R3_Block( Configuration *r3bc , Block * base , Block * father )
  override;

 void map_back_solution( Block *R3B , Configuration *r3bc ,
				      Configuration *solc ) override;

 void map_forward_solution( Block *R3B , Configuration *r3bc ,
				         Configuration *solc ) override;

 bool map_forward_Modification( Block * R3B , c_p_Mod mod ,
                                Configuration * r3bc ,
				ModParam issuePMod , ModParam issueAMod )
  override;

 bool map_back_Modification( Block *R3B , c_p_Mod mod , Configuration *r3bc ,
			     ModParam issuePMod , ModParam issueAMod )
  override;
*/

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
 /// returns a UCBlockSolution with the current solution of this UCBlock
 /** Returns a UCBlockSolution representing the current solution status of
  * this UCBlock. What kind of solution is saved depends on the integer value
  * ws, obtained as follows:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then
  *   ws == solc->f_value:
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, ws is its f_value
  *
  * - otherwise ws is 7 (only save the UnitBlock(s) and NetworkBlock Solution,
  *   with the latter in compressed format)
  *
  * The encoding of ws is bit-wise:
  *
  * - bit 0 (& 1): save the Solution of all the UnitBlock;
  *
  * - bit 1 (& 2): save the Solution of all the NetworkBlock; this bit (and,
  *   therefore, the next one) is ignored if the network only has one node,
  *   since then there is no NetworkBlock to be saved;
  *
  * - bit 2 (& 4): save the NetworkBlock in compressed format;
  *
  * - bit 3 (& 8): save the dual values of the node injection constraints
  *   (1);
  *
  * - bit 4 (& 16): save the dual values of the primary reserve constraints
  *   (2);
  *
  * - bit 5 (& 32): save the dual values of the secondary reserve
  *   constraints (3);
  *
  * - bit 6 (& 64): save the dual values of the inertia constraints (4);
  *
  * - bit 7 (& 128): save the dual values of the pollutant budget
  *   constraints (5).
  *
  * The dual values are those of the RowConstraint (see RowConstraint), as
  * the Solver attached to the UCBlock has written them: with the rows
  * written as in generate_abstract_constraints(), the dual value of a row
  * is the coefficient of its left-hand side in the Lagrangian, i.e., minus
  * the derivative of the optimal value with respect to the finite side that
  * is active. Hence the dual value \f$ y^{ac}_{t,n} \f$ of (1) is free,
  * \f$ - y^{ac}_{t,n} \f$ being the marginal cost of the demand of node n at
  * instant t (the locational marginal price), those of (2)-(4) are
  * nonpositive, their opposites being the marginal costs of the
  * requirements, and that of (5) is nonnegative when the budget is active
  * and nonpositive when the lower bound is; all of them are values per
  * instant, but that of (5), which is per unit of emission over the
  * horizon (see \ref ucbm_dual_sign).
  *
  * Note that UCBlock may not contain some or all of the required solution,
  * if the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception, unless emptys = true, in which case the
  * UCBlockSolution object is only prepped for getting a solution, but it is
  * not really getting one yet.
  *
  * Note that, although the method clearly returns a UCBlockSolution,
  * formally the return type is Solution *. This is because it is not
  * possible to forward declare UCBlockSolution as a derived class from
  * Solution, nor to define UCBlockSolution before UCBlock because the former
  * uses some type information declared in the latter. */

 Solution * get_Solution( Configuration *solc = nullptr ,
			  bool emptys = true ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the related blocks to the UC problem.
 * @{ */

 /// returns the sense of the Objective of this UCBlock
 /** This function returns the sense of the Objective of this UCBlock, which
  * is minimization (Objective::eMin). */

 int get_objective_sense( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the time horizon of the problem

 Index get_time_horizon( void ) const { return( f_time_horizon ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of UnitBlock

 Index get_number_units( void ) const { return( f_number_units ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of NetworkBlocks

 Index get_number_networks( void ) const { return( f_number_networks ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of nodes

 Index get_number_nodes( void ) const {
  return( f_NetworkData ? f_NetworkData->get_number_nodes() : 1 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of primary zones of the problem

 Index get_number_primary_zones( void ) const {
  return( f_number_primary_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of secondary zones of the problem

 Index get_number_secondary_zones( void ) const {
  return( f_number_secondary_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of inertia zones of the problem

 Index get_number_inertia_zones( void ) const {
  return( f_number_inertia_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of pollutants of the problem

 Index get_number_pollutants( void ) const { return( f_number_pollutants ); }

 /// returns true if the given node belongs to the given primary zone
 /** This function returns true if and only if the node identified by \p
  * node_id belongs to the primary zone identified by \p zone_id.
  *
  * @param node_id The ID of a node.
  *
  * @param zone_id The ID of a primary zone.
  *
  * @return True if and only if the given node belongs to the given primary
  *         zone. */

 bool node_belongs_to_primary_zone( Index node_id , Index zone_id ) const {
  // without the vector there is at most one zone, with all the nodes
  // [deserialize() refuses the vector absent with more than one zone]
  return( v_primary_zones.empty() ||
          ( v_primary_zones[ node_id ] == zone_id ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if the given node belongs to the given secondary zone
 /** This function returns true if and only if the node identified by \p
  * node_id belongs to the secondary zone identified by \p zone_id.
  *
  * @param node_id The ID of a node.
  *
  * @param zone_id The ID of a secondary zone.
  *
  * @return True if and only if the given node belongs to the given secondary
  *         zone. */

 bool node_belongs_to_secondary_zone( Index node_id , Index zone_id ) const {
  // without the vector there is at most one zone, with all the nodes
  // [deserialize() refuses the vector absent with more than one zone]
  return( v_secondary_zones.empty() ||
          ( v_secondary_zones[ node_id ] == zone_id ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if the given node belongs to the given inertia zone
 /** This function returns true if and only if the node identified by \p
  * node_id belongs to the inertia zone identified by \p zone_id.
  *
  * @param node_id The ID of a node.
  * @param zone_id The ID of an inertia zone.
  *
  * @return True if and only if the given node belongs to the given inertia
  *         zone. */

 bool node_belongs_to_inertia_zone( Index node_id , Index zone_id ) const {
  // without the vector there is at most one zone, with all the nodes
  // [deserialize() refuses the vector absent with more than one zone]
  return( v_inertia_zones.empty() ||
          ( v_inertia_zones[ node_id ] == zone_id ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if the given \p generator belongs to the given node

 bool generator_belongs_to_node( Index generator , Index node_id )
  const {
  if( ( get_number_nodes() > 1 ) &&
      ( node_id != v_generator_node[ generator ] ) )
   return( false );
  return( true );
  }

/*--------------------------------------------------------------------------*/
 /// returns the NetworkData object
 /** Note that no NetworkData may be defined (see comments to deserialize()),
  * which means that the transmission network is a "bus"; in this case, this
  * method will return nullptr. */

 NetworkBlock::NetworkData * get_NetworkData( void ) const {
  return( f_NetworkData );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of (pointers to) NetworkBlock elements.
 /** This vector contains one entry for each NetworkBlock defined or created
  * for the UCBlock. Its size matches get_number_networks(), and each
  * NetworkBlock may cover one or more time intervals within the
  * "TimeHorizon". In scenarios where the network is a transmission network,
  * i.e., a DCNetworkBlock, with a single bus (NumberNodes == 1), this vector
  * will be empty. */

 const std::vector< NetworkBlock * > & get_network_blocks( void ) const {
  return( v_network_blocks );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of active power demand
 /** This method returns a two-dimensional boost::multi_array<> M such that,
  * if it is not empty, M[ n , t ] gives the active power demand of node n at
  * the time instant t. If it is empty, the active power demand can be found
  * in each NetworkBlock of this UCBlock. */

 const boost::multi_array< double , 2 > &
 get_active_power_demand( void ) const {
  return( v_active_power_demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary zones
 /** Returns the vector V of the primary zones of the nodes, as read from
  * "PrimaryZones" (see deserialize()):
  *
  * - if V is empty, then either there are no primary zones, or there is
  *   only one and all the nodes belong to it;
  *
  * - otherwise, V has size get_number_nodes(), and V[ n ] tells to which
  *   primary zone node n belongs; if V[ n ] >= get_number_primary_zones(),
  *   then node n does not belong to any primary zone, and hence the
  *   corresponding electrical generators are not involved in the primary
  *   constraints, whatever the number of zones (one included). */

 const std::vector< Index > & get_primary_zone( void ) const {
  return( v_primary_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary zones
 /** Returns the vector V of the secondary zones of the nodes, as read from
  * "SecondaryZones" (see deserialize()):
  *
  * - if V is empty, then either there are no secondary zones, or there is
  *   only one and all the nodes belong to it;
  *
  * - otherwise, V has size get_number_nodes(), and V[ n ] tells to which
  *   secondary zone node n belongs; if V[ n ] >= get_number_secondary_zones(),
  *   then node n does not belong to any secondary zone, and hence the
  *   corresponding electrical generators are not involved in the secondary
  *   constraints, whatever the number of zones (one included). */

 const std::vector< Index > & get_secondary_zone( void ) const {
  return( v_secondary_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inertia zones
 /** Returns the vector V of the inertia zones of the nodes, as read from
  * "InertiaZones" (see deserialize()):
  *
  * - if V is empty, then either there are no inertia zones, or there is
  *   only one and all the nodes belong to it;
  *
  * - otherwise, V has size get_number_nodes(), and V[ n ] tells to which
  *   inertia zone node n belongs; if V[ n ] >= get_number_inertia_zones(),
  *   then node n does not belong to any inertia zone, and hence the
  *   corresponding electrical generators are not involved in the inertia
  *   constraints, whatever the number of zones (one included). */

 const std::vector< Index > & get_inertia_zone( void ) const {
  return( v_inertia_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of primary demand
 /** Returns the two-dimensional boost::multi_array<> M such that M[ z ][ t ]
  * is the primary reserve requirement of zone z at instant t (see
  * "PrimaryDemand" in deserialize()); M is empty() if there are no
  * primary zones, and has get_number_primary_zones() rows of
  * get_time_horizon() entries otherwise. */

 const boost::multi_array< double , 2 > & get_primary_demand( void ) const {
  return( v_primary_demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of secondary demand
 /** Returns the two-dimensional boost::multi_array<> M such that M[ z ][ t ]
  * is the secondary reserve requirement of zone z at instant t (see
  * "SecondaryDemand" in deserialize()); M is empty() if there are no
  * secondary zones, and has get_number_secondary_zones() rows of
  * get_time_horizon() entries otherwise. */

 const boost::multi_array< double , 2 > & get_secondary_demand( void ) const {
  return( v_secondary_demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inertia demand
 /** Returns the two-dimensional boost::multi_array<> M such that M[ z ][ t ]
  * is the inertia requirement of zone z at instant t (see
  * "InertiaDemand" in deserialize()); M is empty() if there are no
  * inertia zones, and has get_number_inertia_zones() rows of
  * get_time_horizon() entries otherwise. */

 const boost::multi_array< double , 2 > & get_inertia_demand( void ) const {
  return( v_inertia_demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of pollutant zones associated with each pollutant
 /** This method returns a vector containing the number of pollutant zones for
  * each pollutant. The i-th entry of this vector is the number of pollutant
  * zones associated with pollutant i. */

 const std::vector< Index > & get_number_pollutant_zones( void ) const {
  return( v_number_pollutant_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the total number of pollutant zones
 /** This method returns the sum, over all pollutants, of the number of
  * pollutant zones associated with each pollutant (cf.
  * get_number_pollutant_zones()). */

 Index get_total_number_pollutant_zones( void ) const {
  return( f_total_number_pollutant_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of pollutant zones
 /** The method returns a two-dimensional boost::multi_array<> M such that
  * M[ p , n ] tells to which pollutant zone associated with pollutant p the
  * node n belongs. There are two possible cases:
  *
  * - if the boost::multi_array<> M is empty() then either there are no
  *   pollutants, or each pollutant has exactly one zone and all the nodes
  *   belong to it;
  *
  * - otherwise, M has get_number_pollutants() rows and get_number_nodes()
  *   columns, and M[ p , n ] tells to which pollutant zone associated with
  *   pollutant p the node n belongs; if M[ p , n ] >=
  *   get_number_pollutant_zones()[ p ], node n belongs to no zone of
  *   pollutant p. */

 const boost::multi_array< Index , 2 > & get_pollutant_zone( void ) const {
  return( v_pollutant_zones );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of pollutant budget
 /** The method returns the vector V of the pollutant budgets (across all the
  * time horizon) of the zones of all the pollutants, one after the other as
  * in the netCDF variable PollutantBudget (see deserialize()): the budget of
  * zone z of pollutant p is
  *
  *    V[ get_number_pollutant_zones()[ 0 ] + ... +
  *       get_number_pollutant_zones()[ p - 1 ] + z ] .
  *
  * V is empty() if there are no pollutants, otherwise its size is
  * get_total_number_pollutant_zones(). */

 const std::vector< double > & get_pollutant_budget( void ) const {
  return( v_pollutant_budget );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of pollutant rho
 /** The method returns a three-dimensional boost::multi_array<> M such that
  * M[ t , p , g ] gives the production of pollutant p from electrical
  * generator g at time t. This three-dimensional boost::multi_array<> M
  * considers two possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no pollutant zones are
  *   defined, and there are no pollutant budget constraints;
  *
  * - otherwise, two possible cases may happen to the first dimension of the
  *   M[ t , p , g ];
  *
  *   - if the first dimension of the boost::multi_array<> M has size one,
  *     then each element of the matrix M [ 0 , p , g ] gives the conversion
  *     factor of pollutant p due to the electrical generator g;
  *
  *   - if the first dimension of the boost::multi_array<> M has full size
  *     then, each element of the matrix M[ t , p , g ] gives the conversion
  *     factor of pollutant p due to the electrical generator g for time
  *     instant t. */

 const boost::multi_array< double , 3 > & get_pollutant_rho( void ) const {
  return( v_pollutant_rho );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the conversion factor of pollutant p for generator g at time t
 /** This method returns the conversion factor of pollutant \p p due to the
  * electrical generator \p g at time \p t, whether or not the factors
  * depend on time (cf. get_pollutant_rho()). It must not be called if
  * get_number_pollutants() == 0. */

 double get_pollutant_rho( Index t , Index p , Index g ) const {
  return( v_pollutant_rho[ v_pollutant_rho.shape()[ 0 ] > 1 ? t : 0 ]
                         [ p ][ g ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of the lower bounds of the pollutant budget
 /** The method returns the vector of the lower bounds on the emissions of the
  * zones of all the pollutants, ordered as get_pollutant_budget(); an entry
  * is -INF if there is no lower bound. */

 const std::vector< double > & get_pollutant_min_budget( void ) const {
  return( v_pollutant_min_budget );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of storages of all the units

 Index get_number_storages( void ) const { return( f_number_storages ); }

/*--------------------------------------------------------------------------*/
 /// returns the factors of the storage levels in the pollutant budget
 /** The method returns the three-dimensional boost::multi_array<> M such that
  * M[ t , p , k ] is the factor of the level of storage k at the end of time
  * t in the pollutant budget constraints of pollutant p, with the first
  * dimension of size 1 if the factors do not depend on time; M is empty() if
  * no storage level is in those constraints. */

 const boost::multi_array< double , 3 > & get_pollutant_storage_rho( void )
  const {
  return( v_pollutant_storage_rho );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the factor of the level of storage k for pollutant p at time t
 /** It must not be called if get_pollutant_storage_rho() is empty(). */

 double get_pollutant_storage_rho( Index t , Index p , Index k ) const {
  return( v_pollutant_storage_rho[ v_pollutant_storage_rho.shape()[ 0 ] > 1 ?
                                   t : 0 ][ p ][ k ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the u-th UnitBlock

 UnitBlock * get_unit_block( Index u ) const {
  if( u >= f_number_units )
   throw( std::invalid_argument( "invalid unit index" ) );
  return( static_cast< UnitBlock * >( v_Block[ u ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the n-th NetworkBlock

 NetworkBlock * get_network_block( Index n ) const {
  if( v_network_blocks.empty() )
   return( nullptr );
  if( n >= f_number_networks )
   throw( std::invalid_argument( "invalid network index" ) );

  return( v_network_blocks[ n ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of generator node
 /** This method returns the vector V that indicates to which node of the
  * transmission network each electrical generator belongs, as read from
  * "GeneratorNode" (see deserialize()): V[ g ] is the node of electrical
  * generator g, for g = 0, ..., NumberElectricalGenerators - 1. V may be
  * empty in the bus case, where all the generators are at the only node. */

 const std::vector< Index > & get_generator_node( void ) const {
  return( v_generator_node );
  }

/*--------------------------------------------------------------------------*/
 /// returns the node injection constraints
 /** This method returns the boost multi_array C containing the node
  * injection constraints (1) (see generate_abstract_constraints()).
  * C[ t ][ n ] is the node injection constraint associated with time t and
  * node n, with n = 0 only in the bus case, where it is (1b). The dual value
  * of each row (RowConstraint::get_dual()) follows the convention described
  * in get_Solution(): minus the locational marginal price, per instant. */

 boost::multi_array< FRowConstraint , 2 > &
  get_node_injection_constraints( void ) {
  return( v_node_injection_Const );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the reactive power node injection constraints
 /** This method returns the boost multi_array C containing the reactive
  * power node injection constraints (1') (see
  * generate_abstract_constraints()), indexed as those of
  * get_node_injection_constraints(); C is empty if there is no reactive
  * power in the model. */

 boost::multi_array< FRowConstraint , 2 > &
  get_reactive_node_injection_constraints( void ) {
  return( v_reactive_node_injection_Const );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the (const) node injection constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the node injection constraints. C[ t ][ n ] is the node
  * injection constraint associated with time t and node n. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_const_node_injection_constraints( void ) const {
  return( v_node_injection_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the primary demand constraints
 /** This method returns the boost multi_array C containing the primary
  * reserve constraints (2) (see generate_abstract_constraints()).
  * C[ t ][ z ] is the primary reserve constraint associated with time t and
  * primary zone z; its dual value is nonpositive, its opposite being the
  * marginal cost of the requirement (see get_Solution()). */

 boost::multi_array< FRowConstraint , 2 > &
 get_primary_demand_constraints( void ) {
  return( v_PrimaryDemand_Const );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the (const) primary demand constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the primary demand constraints. C[ t ][ z ] is the primary
  * demand constraint associated with time t and primary zone z. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_const_primary_demand_constraints( void ) const {
  return( v_PrimaryDemand_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the secondary demand constraints
 /** This method returns the boost multi_array C containing the secondary
  * reserve constraints (3) (see generate_abstract_constraints()).
  * C[ t ][ z ] is the secondary reserve constraint associated with time t
  * and secondary zone z; its dual value is nonpositive, its opposite being
  * the marginal cost of the requirement (see get_Solution()). */

 boost::multi_array< FRowConstraint , 2 > &
 get_secondary_demand_constraints( void ) {
  return( v_SecondaryDemand_Const );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the (const) secondary demand constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the secondary demand constraints. C[ t ][ z ] is the secondary
  * demand constraint associated with time t and secondary zone z. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_const_secondary_demand_constraints( void ) const {
  return( v_SecondaryDemand_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia demand constraints
 /** This method returns the boost multi_array C containing the inertia
  * constraints (4) (see generate_abstract_constraints()). C[ t ][ z ] is
  * the inertia constraint associated with time t and inertia zone z; its
  * dual value is nonpositive, its opposite being the marginal cost of the
  * requirement (see get_Solution()). */

 boost::multi_array< FRowConstraint , 2 > &
 get_inertia_demand_constraints( void ) {
  return( v_InertiaDemand_Const );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the (const) inertia demand constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the inertia demand constraints. C[ t ][ z ] is the inertia
  * demand constraint associated with time t and inertia zone z. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_const_inertia_demand_constraints( void ) const {
  return( v_InertiaDemand_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the pollutant budget constraints
 /** This method returns the vector C of the pollutant budget constraints
  * (5) (see generate_abstract_constraints()): C[ p ] has
  * get_number_pollutant_zones()[ p ] entries, C[ p ][ z ] being the
  * constraint of zone z of pollutant p, whose budget is
  * get_pollutant_budget()[ get_number_pollutant_zones()[ 0 ] + ... +
  * get_number_pollutant_zones()[ p - 1 ] + z ]. Its dual value is
  * nonnegative when the budget is active and nonpositive when the lower
  * bound is (see get_Solution()). */

 std::vector< std::vector< FRowConstraint > > &
  get_pollutant_constraints( void ) {
  return( v_PollutantBudget_Const );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the (const) pollutant budget constraints
 /** This method returns (a const reference to) the pollutant budget
  * constraints, indexed as in get_pollutant_constraints(). */

 const std::vector< std::vector< FRowConstraint > > &
 get_const_pollutant_constraints( void ) const {
  return( v_PollutantBudget_Const );
  }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UCBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing and saving the UCBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * UCBlock. See UCBlock::deserialize( netCDF::NcGroup ) for the format's
  * details of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UCBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UCBlock
 *  @{ */

 /// loading a UCBlock from a stream is not implemented: it throws

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "UCBlock::load() not implemented yet" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for changing the data of the UCBlock
 *  @{ */

 /// method for handling Modification
 /** This method has to intercept any "abstract Modification" that modifies
  * the "abstract representation" of the UCBlock, and "translate" them into
  * both changes of the actual data structures and corresponding "physical
  * Modification". These Modification(s) are those for which
  * Modification::concerns_Block() is true. This method handles the
  * UnitBlockMod of type UnitBlockMod::eScale, also inside a
  * GroupModification: the coefficients of the terms of the scaled units in
  * the constraints (1)-(5) and (1') of generate_abstract_constraints(), and
  * the right-hand sides of (1), are changed to the new scale. It also handles
  * the change of the inertia power of a HydroUnitBlock (also inside a
  * HydroSystemUnitBlock), whose coefficients in (4) follow it, a term that
  * the rows do not have being refused with std::logic_error; and the
  * changes of the scale, of the kappa and of the maximum power of the units
  * make it give the NetworkBlock the bounds on the node injections again
  * (see \ref ucbm_link_bounds and set_node_injection_bounds()). The rows
  * are rewritten by walking them as they are generated, node by node and
  * unit by unit, so that a unit with generators at several nodes of a zone
  * is handled. */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// update the active power demand
 /** Updates the active power demand for the entries specified by \p subset.
  *
  * The active power demand is conceptually stored as a matrix
  * ActivePowerDemand[ n ][ t ], where n is the node index and t is the time
  * instant. The iterator \p values must therefore provide the new demand
  * values in the same order as the entries identified by \p subset, where
  * each element of \p subset is interpreted as the flattened index
  *
  *    index = n * get_time_horizon() + t .
  *
  * Hence, for a generic element \c index in \p subset, the corresponding
  * demand value read from \p values is assigned to
  *
  *    ActivePowerDemand[ index / get_time_horizon() ]
  *                     [ index % get_time_horizon() ] .
  *
  * Therefore, the underlying flattened layout is node-major, i.e., all time
  * instants of node 0 come first, then all time instants of node 1, and so
  * on:
  *
  *    [ (0,0), (0,1), ... , (0,T-1), (1,0), ... , (N-1,T-1) ] .
  *
  * If \p ordered is true, \p subset is assumed to be already sorted in
  * increasing order; otherwise it may be reordered internally before issuing
  * the corresponding Modification. If NetworkBlock-s are present, the method
  * forwards each updated entry to the corresponding NetworkBlock covering the
  * associated time interval. */

 void set_active_power_demand( MF_dbl_it values , Subset && subset = { 0 } ,
                               bool ordered = false ,
                               ModParam issuePMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the active power demand
 /** Updates the active power demand for the flattened range \p rng.
  *
  * The active power demand is conceptually stored as a matrix
  * ActivePowerDemand[ n ][ t ], where n is the node index and t is the time
  * instant. The iterator \p values must provide one value for each flattened
  * entry in the half-open range [ rng.first , rng.second ), in increasing
  * index order, with each flattened index interpreted as
  *
  *    index = n * get_time_horizon() + t .
  *
  * Hence, the k-th value read from \p values is assigned to the entry
  * corresponding to the k-th flattened index in the range, namely
  *
  *    ActivePowerDemand[ index / get_time_horizon() ]
  *                     [ index % get_time_horizon() ] .
  *
  * Therefore, the underlying flattened layout is node-major, i.e., all time
  * instants of node 0 come first, then all time instants of node 1, and so
  * on:
  *
  *    [ (0,0), (0,1), ... , (0,T-1), (1,0), ... , (N-1,T-1) ] .
  *
  * If NetworkBlock-s are present, each updated entry is forwarded to the
  * corresponding NetworkBlock covering the associated time interval. */

 void set_active_power_demand( MF_dbl_it values , Range rng = Range( 0 , 1 ) ,
                               ModParam issuePMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the pollutant budget
 /** Updates the pollutant budget for the (flattened) indices in \p subset.
  *
  * The budget is indexed as get_pollutant_budget(), i.e., as the netCDF
  * variable "PollutantBudget" (see deserialize()): all the zones of pollutant
  * 0 come first, then all the zones of pollutant 1, and so on, i.e., index k
  * corresponds to zone z of pollutant p with
  *
  *    k = NumberPollutantZones[ 0 ] + ... + NumberPollutantZones[ p - 1 ]
  *        + z .
  *
  * The iterator \p values must provide one value for each element of
  * \p subset, in the same order. If \p ordered is true, \p subset is assumed
  * to be already sorted in increasing order; otherwise it may be reordered
  * internally before issuing the corresponding Modification. */

 void set_pollutant_budget( MF_dbl_it values , Subset && subset ,
                            bool ordered = false ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the pollutant budget
 /** Updates the pollutant budget for the (flattened) indices in the
  * half-open range [ rng.first , rng.second ), with the same flattening as
  * set_pollutant_budget( subset ). The iterator \p values must provide one
  * value for each index in the range, in increasing index order. */

 void set_pollutant_budget( MF_dbl_it values ,
                            Range rng = Range( 0 , Inf< Index >() ) ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the lower bound of the pollutant budget
 /** As set_pollutant_budget( subset ), but for the lower bounds of the
  * emissions [see get_pollutant_min_budget()]. */

 void set_pollutant_min_budget( MF_dbl_it values , Subset && subset ,
                                bool ordered = false ,
                                ModParam issuePMod = eNoBlck ,
                                ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the lower bound of the pollutant budget
 /** As set_pollutant_budget( range ), but for the lower bounds of the
  * emissions [see get_pollutant_min_budget()]. */

 void set_pollutant_min_budget( MF_dbl_it values ,
                                Range rng = Range( 0 , Inf< Index >() ) ,
                                ModParam issuePMod = eNoBlck ,
                                ModParam issueAMod = eNoBlck );

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED TYPES OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 using MAdouble = boost::multi_array< double , 2 >;  ///< matrix of data
 using MAdouble_ext = MAdouble::extent_gen;  ///< extents of a MAdouble

 using MAFRC = boost::multi_array< FRowConstraint , 2 >;  ///< matrix of rows
 using MAFRC_ext = MAdouble::extent_gen;  ///< extents of a MAFRC

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// states that the Variable of the UCBlock have been generated
 void set_variables_generated( void ) { AR |= HasVar; }

 /// states that the Constraint(s) of the UCBlock have been generated
 void set_constraints_generated( void ) { AR |= HasCst; }

 /// states that the Objective of the UCBlock has been generated
 void set_objective_generated( void ) { AR |= HasObj; }

 /// indicates whether the Variable(s) of the UCBlock have been generated
 bool variables_generated( void ) const { return( AR & HasVar ); }

 /// indicates whether the Constraint(s) of the UCBlock has been generated
 bool constraints_generated( void ) const { return( AR & HasCst ); }

 /// indicates whether the Objective of the UCBlock has been generated
 bool objective_generated( void ) const { return( AR & HasObj ); }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the NetworkData object
 NetworkBlock::NetworkData * f_NetworkData;

 /// the classname of the NetworkBlock that the UCBlock builds
 /** E.g., "DCNetworkBlock" or "ECNetworkBlock" (see
  * "NetworkBlockClassname" in deserialize()). */
 std::string network_block_classname;

 /// the classname of the NetworkData that the UCBlock builds
 std::string network_data_classname;

 /// the time horizon of the problem
 Index f_time_horizon;

 /// the number of the networks of the problem
 Index f_number_networks;

 /// the number of units of the problem
 Index f_number_units;

 /// the number of electrical generators of the problem
 Index f_number_elc_generators;

 /// true if reactive power constraints have to be managed
 bool f_has_reactive;
 
 /// the total number of pollutant zones of the problem
 Index f_total_number_pollutant_zones;

 /// the number of primary zones
 Index f_number_primary_zones;

 /// the number of secondary zones
 Index f_number_secondary_zones;

 /// the number of inertia zones
 Index f_number_inertia_zones;

 /// the number of pollutants
 Index f_number_pollutants;

 /// the constant terms of the NetworkBlock in the bus case
 /** In the bus case, the constant terms of the NetworkBlock, which are not
  * built, indexed over "NumberNetworks"; their sum is the constant of the
  * Objective of the UCBlock. Empty with more than one node, where the
  * constant terms are in the NetworkBlock, or if there is none. */
 std::vector< double > v_network_constant_terms;

 /// the number of pollutant zones of each pollutant
 std::vector< Index > v_number_pollutant_zones;

 /// the matrix of pollutant zones
 /** Indexed over the dimensions NumberPollutants and NumberNodes. */
 boost::multi_array< Index , 2 > v_pollutant_zones;

 /// vector of pointers to the NetworkBlock.
 /** This vector has size get_number_networks(). */
 std::vector< NetworkBlock * > v_network_blocks;

 /// the matrix of ActivePowerDemand
 /** Indexed over the dimensions "NumberNodes" and "TimeHorizon"; empty with
  * more than one node, where the demand is in the NetworkBlock. */
 boost::multi_array< double , 2 > v_active_power_demand;

 /// the matrix of ReactivePowerDemand, indexed as v_active_power_demand
 boost::multi_array< double , 2 > v_reactive_power_demand;

 /// the vector of PrimaryZones
 std::vector< Index > v_primary_zones;

 /// the matrix of PrimaryDemand
 /** Indexed over the dimensions "PrimaryZones" and "TimeHorizon". */
 boost::multi_array< double , 2 > v_primary_demand;

 /// the vector of SecondaryZones
 std::vector< Index > v_secondary_zones;

 /// the matrix of SecondaryDemand
 /** Indexed over the dimensions "SecondaryZones" and "TimeHorizon". */
 boost::multi_array< double , 2 > v_secondary_demand;

 /// the vector of InertiaZones
 std::vector< Index > v_inertia_zones;

 /// the matrix of InertiaDemand
 /** Indexed over the dimensions "InertiaZones" and "TimeHorizon". */
 boost::multi_array< double , 2 > v_inertia_demand;

 /// the PollutantBudget
 /** The budgets of the zones of all the pollutants, one after the other. */
 std::vector< double > v_pollutant_budget;

 /// the matrix of PollutantRho
 /** Indexed over "TimeHorizon" (or a singleton), "NumberPollutants", and
  * "NumberElectricalGenerators". */
 boost::multi_array< double , 3 > v_pollutant_rho;

 /// the PollutantMinBudget, ordered as v_pollutant_budget
 std::vector< double > v_pollutant_min_budget;

 /// the number of storages of all the units
 Index f_number_storages;

 /// the matrix of PollutantStorageRho
 /** Indexed over "TimeHorizon" (or a singleton), "NumberPollutants", and
  * "NumberStorages"; empty if no storage level is in the constraints. */
 boost::multi_array< double , 3 > v_pollutant_storage_rho;

 /// v_generator_node[ g ] tells to which node generator g belongs
 std::vector< Index > v_generator_node;

/*-------------------------------- variables -------------------------------*/

/*------------------------------- constraints ------------------------------*/

 /// node injection constraints for each time and node
 boost::multi_array< FRowConstraint , 2 > v_node_injection_Const;

 /// reactive node injection constraints for each time and node
 boost::multi_array< FRowConstraint , 2 > v_reactive_node_injection_Const;

 /// primary demand constraints for each time and primary zone
 boost::multi_array< FRowConstraint , 2 > v_PrimaryDemand_Const;

 /// secondary demand constraints for each time and secondary zone
 boost::multi_array< FRowConstraint , 2 > v_SecondaryDemand_Const;

 /// inertia demand constraints for each time and inertia zone
 boost::multi_array< FRowConstraint , 2 > v_InertiaDemand_Const;

 /// pollutant budget constraints for each pollutant and pollutant zone
 /** v_PollutantBudget_Const[ p ][ z ] is the constraint of zone z of
  * pollutant p. */
 std::vector< std::vector< FRowConstraint > > v_PollutantBudget_Const;

 FRealObjective objective;  ///< the objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 unsigned char AR{};  ///< bit-wise coded: what abstract is there

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variables have been constructed

 static constexpr unsigned char HasCst = 2;
 ///< second bit of AR == 1 if the Constraints have been constructed

 static constexpr unsigned char HasObj = 4;
 ///< third bit of AR == 1 if the Objective has been constructed

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  ///< registration in the factory

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// deserialize the sub-blocks of UCBlock that have the given prefix name

 void deserialize_sub_blocks( const netCDF::NcGroup & group ,
                              const std::string & prefix ,
                              Index num_sub_blocks );

/*--------------------------------------------------------------------------*/
 /// deserialize the Network Blocks of UCBlock

 void deserialize_network_blocks( const netCDF::NcGroup & group );

/*--------------------------------------------------------------------------*/
 /// generate the node injection constraints (1)

 void generate_node_injection_constraints( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the reactive node injection constraints (1')

 void generate_reactive_node_injection_constraints( void );

/*--------------------------------------------------------------------------*/
 /// generate the primary demand constraints

 void generate_primary_demand_constraints( void );

/*--------------------------------------------------------------------------*/
 /// generate the secondary demand constraints

 void generate_secondary_demand_constraints( void );

/*--------------------------------------------------------------------------*/
 /// generate the inertia demand constraints

 void generate_inertia_demand_constraints( void );

/*--------------------------------------------------------------------------*/
 /// generate the pollutant budget constraints

 void generate_pollutant_budget_constraints( void );

/*--------------------------------------------------------------------------*/
 /// gives each NetworkBlock the bounds of its node injections
 /** The bounds of the injection at a node are the sums, over the generators
  * at the node, of get_scale() * get_kappa() * get_design_ub() times the
  * minimum (or the opposite of the fixed consumption, if smaller) and the
  * maximum power; the same for the reactive power, if any. Called by
  * deserialize() and by add_Modification() whenever one of these data of a
  * unit changes; the NetworkBlock rewrite the rows of the bounds, if they
  * have them. */

 void set_node_injection_bounds( void );

/*--------------------------------------------------------------------------*/
 /// updates the node injection constraints
 /** This function updates the node injection constraints considering that the
  * scale factors of the given units may have been modified. The vector \p
  * modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered.
  *
  * @param issueMod How the Modification of the rows are issued.
  *
  * @param reactive If true, the reactive node injection constraints are
  *        updated instead of the active ones. */

 void update_node_injection_constraints(
			       const std::vector< Index > & modified_units ,
			       ModParam issueMod = eNoBlck ,
			       bool reactive = false );

/*--------------------------------------------------------------------------*/
 /// updates the primary demand constraints
 /** This function updates the primary demand constraints considering that the
  * scale factors of the given units may have been modified. The vector \p
  * modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered.
  *
  * @param issueMod How the Modification of the rows are issued. */

 void update_primary_demand_constraints(
			       const std::vector< Index > & modified_units ,
			       ModParam issueMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the secondary demand constraints
 /** This function updates the secondary demand constraints considering that
  * the scale factors of the given units may have been modified. The vector
  * \p modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered.
  *
  * @param issueMod How the Modification of the rows are issued. */

 void update_secondary_demand_constraints(
			       const std::vector< Index > & modified_units ,
			       ModParam issueMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the inertia demand constraints
 /** This function updates the inertia demand constraints considering that
  * the scale factors of the given units may have been modified. The vector
  * \p modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered.
  *
  * @param issueMod How the Modification of the rows are issued. */

 void update_inertia_demand_constraints(
			       const std::vector< Index > & modified_units ,
			       ModParam issueMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the pollutant budget constraints
 /** This function updates the pollutant budget constraints considering that
  * the scale factors of the given units may have been modified. The vector
  * \p modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered.
  *
  * @param issueMod How the Modification of the rows are issued. */

 void update_pollutant_budget_constraints(
			       const std::vector< Index > & modified_units ,
			       ModParam issueMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the pollutant budget of the (flattened) index k
 /** Sets the budget of the pollutant zone with (flattened) index \p k [see
  * set_pollutant_budget()], its lower bound if \p lower, to \p budget, and
  * changes the corresponding side of the constraint, if any, according to
  * \p issueAMod. Returns true if the budget has changed. */

 bool set_pollutant_budget_k( Index k , double budget , bool lower ,
                              ModParam issuePMod , ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// sets the pollutant budget, or its lower bound, of the given indices
 /** Does what set_pollutant_budget( subset ) does, or
  * set_pollutant_min_budget( subset ) if \p lower is true. */

 void set_pollutant_bounds( MF_dbl_it values , Subset && subset ,
                            bool ordered , bool lower ,
                            c_ModParam issuePMod , c_ModParam issueAMod );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the pollutant budget, or its lower bound, of the given range
 /** Does what set_pollutant_budget( range ) does, or
  * set_pollutant_min_budget( range ) if \p lower is true. */

 void set_pollutant_bounds( MF_dbl_it values , Range rng , bool lower ,
                            c_ModParam issuePMod , c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// calls visit on each term of the pollutant budget constraints of p
 /** Calls visit( unit_id , unit_block , zone_id , var , factor ) for each
  * term of the constraints of pollutant \p p, unit by unit and in the order
  * in which they are in the constraints: the active power of each generator
  * of the unit at each time, then the level of each of its storages at each
  * time; \p factor is the one of the data, i.e., not yet multiplied by the
  * scale of the unit, and the terms with a zero factor are not visited. */

 template< class F >
 void for_each_pollutant_term( Index p , F && visit );

/*--------------------------------------------------------------------------*/
 /// updates a node injection constraint for the given demand
 /** This function updates the node injection constraint at the given \p time
  * for the node whose index is \p node_index considering the given \p demand.
  *
  * @param time A time between 0 and get_time_horizon() - 1.
  *
  * @param node_index The index of a node.
  *
  * @param demand The demand at the given node at the given time. */

 void update_node_injection_constraints( Index time , Index node_index ,
                                         double demand );

/*--------------------------------------------------------------------------*/
 /// returns the primary zone to which the given generator belongs

 Index get_primary_zone( Index elc_generator ) const {
  // without the vector all the nodes are in zone 0 [see
  // node_belongs_to_primary_zone()]; the zone returned is
  // >= get_number_primary_zones() if the node is in no zone
  if( v_primary_zones.empty() )
   return( 0 );

  // Node to which the given electrical generator belongs
  Index node = 0;
  if( get_number_nodes() > 1 )
   node = v_generator_node[ elc_generator ];

  return( v_primary_zones[ node ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the secondary zone to which the given generator belongs

 Index get_secondary_zone( Index elc_generator ) const {
  // without the vector all the nodes are in zone 0 [see
  // node_belongs_to_secondary_zone()]; the zone returned is
  // >= get_number_secondary_zones() if the node is in no zone
  if( v_secondary_zones.empty() )
   return( 0 );

  // Node to which the given electrical generator belongs
  Index node = 0;
  if( get_number_nodes() > 1 )
   node = v_generator_node[ elc_generator ];

  return( v_secondary_zones[ node ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia zone to which the given generator belongs

 Index get_inertia_zone( Index elc_generator ) const {
  // without the vector all the nodes are in zone 0 [see
  // node_belongs_to_inertia_zone()]; the zone returned is
  // >= get_number_inertia_zones() if the node is in no zone
  if( v_inertia_zones.empty() )
   return( 0 );

  // Node to which the given electrical generator belongs
  Index node = 0;
  if( get_number_nodes() > 1 )
   node = v_generator_node[ elc_generator ];

  return( v_inertia_zones[ node ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the zone of pollutant p to which the given node belongs

 Index get_pollutant_zone_of_node( Index p , Index node ) const {
  if( v_pollutant_zones.empty() )
   return( 0 );  // the unique zone of each pollutant
  return( v_pollutant_zones[ p ][ node ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the zone of pollutant p to which the given generator belongs

 Index get_pollutant_zone( Index p , Index elc_generator ) const {
  // Node to which the given electrical generator belongs
  Index node = 0;
  if( get_number_nodes() > 1 )
   node = v_generator_node[ elc_generator ];

  return( get_pollutant_zone_of_node( p , node ) );
  }

/*--------------------------------------------------------------------------*/

 /// registers the methods of the UCBlock in the method factory

 static void static_initialization( void )
 {
  register_method< UCBlock , MF_dbl_it , Subset && , bool >(
   "UCBlock::set_active_power_demand" , & UCBlock::set_active_power_demand );

  register_method< UCBlock , MF_dbl_it , Range >(
   "UCBlock::set_active_power_demand" , & UCBlock::set_active_power_demand );

  register_method< UCBlock , MF_dbl_it , Subset && , bool >(
   "UCBlock::set_pollutant_budget" , & UCBlock::set_pollutant_budget );

  register_method< UCBlock , MF_dbl_it , Range >(
   "UCBlock::set_pollutant_budget" , & UCBlock::set_pollutant_budget );

  register_method< UCBlock , MF_dbl_it , Subset && , bool >(
   "UCBlock::set_pollutant_min_budget" , & UCBlock::set_pollutant_min_budget );

  register_method< UCBlock , MF_dbl_it , Range >(
   "UCBlock::set_pollutant_min_budget" , & UCBlock::set_pollutant_min_budget );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( UCBlock ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS UCBlockMod ----------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a UCBlock

 class UCBlockMod : public Modification
{
 public:

 /// public enum for the types of UCBlockMod
 enum UCB_mod_type
 {
  eSetActD = 0 ,  ///< set active power demand
  eSetPolB ,      ///< set pollutant budget
  eSetPolMinB     ///< set the lower bound of the pollutant budget
  };

 /// constructor, takes the UCBlock and the type
 UCBlockMod( UCBlock * const fblock , const int type )
  : f_Block( fblock ) , f_type( type ) {}

 /// destructor, does nothing
 virtual ~UCBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 /// accessor to the type of modification
 int type( void ) const { return( f_type ); }

 protected:

 /// prints the UCBlockMod
 void print( std::ostream & output ) const override {
  output << "UCBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetPolB ):
    output << "Set pollutant budget";
    break;
   case( eSetPolMinB ):
    output << "Set pollutant minimum budget";
    break;
   default:
    output << "Set active power demand";
   }
  }

 UCBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

 int f_type;  ///< type of modification

 };  // end( class( UCBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS UCBlockRngdMod ---------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from UCBlockMod for "ranged" modifications

class UCBlockRngdMod : public UCBlockMod
{
 public:

 /// constructor: takes the UCBlock, the type, and the range
 UCBlockRngdMod( UCBlock * const fblock , const int type , Block::Range rng )
  : UCBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~UCBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the UCBlockRngdMod
 void print( std::ostream & output ) const override {
  UCBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

 };  // end( class( UCBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS UCBlockSbstMod -------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from UCBlockMod for "subset" modifications

class UCBlockSbstMod : public UCBlockMod
{
 public:

 /// constructor: takes the UCBlock, the type, and the subset
 UCBlockSbstMod( UCBlock * const fblock , const int type ,
                 Block::Subset && nms )
  : UCBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~UCBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the UCBlockSbstMod
 void print( std::ostream & output ) const override {
  UCBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

 Block::Subset f_nms;  ///< the subset

 };  // end( class( UCBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS UCBlockSolution ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a Solution of a UCBlock
/** The UCBlockSolution class, derived from Solution, represents a solution
 * of a UCBlock, i.e.:
 *
 * - [optionally] the [:UnitBlock]Solution of all the units;
 *
 * - [optionally] the [:NetworkBlock]Solution of all the networks;
 *
 * - [optionally] the dual values of the node injection constraints (1);
 *
 * - [optionally] the dual values of the primary reserve constraints (2);
 *
 * - [optionally] the dual values of the secondary reserve constraints (3);
 *
 * - [optionally] the dual values of the inertia constraints (4);
 *
 * - [optionally] the dual values of the pollutant budget constraints (5),
 *
 * the constraints being numbered as in
 * UCBlock::generate_abstract_constraints(), and the dual values following
 * the convention described in UCBlock::get_Solution(). */

class UCBlockSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 using Index = Block::Index;  ///< "import" Index

/*------------------------------- FRIENDS ----------------------------------*/

 friend UCBlock;  ///< make UCBlock friend

/*-------------- CONSTRUCTING AND DESTRUCTING UCBlockSolution --------------*/

 explicit UCBlockSolution( void ) : f_time_horizon( 0 ) ,
  f_number_nodes( 0 ) , f_number_primary_zones( 0 ) ,
  f_number_secondary_zones( 0 ) , f_number_inertia_zones( 0 ) ,
  f_total_number_pollutant_zones( 0 ) , f_compressed_network( true ) {}
 /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~UCBlockSolution() {
  for( auto ubs : v_unit_Solution )
   delete ubs;
  for( auto nbs : v_network_Solution )
   delete nbs;
  }

/*----------- METHODS DESCRIBING THE BEHAVIOR OF A UCBlockSolution ---------*/

 /// reads the solution of the UCBlock \p block into this UCBlockSolution
 void read( const Block * block ) override final;

 /// writes this UCBlockSolution into the UCBlock \p block
 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a UCBlockSolution into a netCDF::NcGroup
 /** Serialize a UCBlockSolution into a netCDF::NcGroup, with the following
  * format:
  *
  * - The dimension "TimeHorizon" containing the number of time steps in the
  *   problem. It is mandatory.
  *
  * - The dimension "NumberUnits" containing the number of units (UnitBlock)
  *   in the problem; the dimension is optional, if it is missing then no
  *   unit Solution (see "UnitBlock_i" below) is present.
  *
  * - The groups "UnitBlock_0", "UnitBlock_1", ..., "UnitBlock_n" with n ==
  *   NumberUnits - 1, containing each the UnitBlockSolution corresponding
  *   to that electrical generator. If NumberUnits is present, it is an error
  *   if the corresponding groups are not there.
  *
  * - The dimension "NumberNetworks" containing the number of networks
  *   (NetworkBlock) in the problem; the dimension is optional, if it is
  *   missing then no network Solution (see "NetworkBlock*" below) is
  *   present.
  *
  * - There are two possible versions of network information:
  *
  *   - The "standard" one, i.e., groups "NetworkBlock_0", "NetworkBlock_1",
  *     ..., "NetworkBlock_T" with T = NumberNetworks - 1, with
  *     "NetworkBlock_t" containing each the NetworkBlockSolution
  *     corresponding to that network constraints for some specific subset
  *     of time instants t. This is the most flexible case, as it allows to
  *     have different types of :NetworkBlock for each t. However, since T
  *     may be large, this may create performance problems since netCDF does
  *     not like to have many groups.
  *
  *   - The single group "NetworkBlock" that contains all the data of all the
  *     NetworkBlockSolution corresponding to all the time instants t in
  *     0, ..., TimeHorizon - 1; see NetworkBlock::serialize( group & , int )
  *     for a description of the format. This has much better performances,
  *     but requires that all :NetworkBlock are actually of the same type
  *     for each time instant t.
  *
  *   If NumberNetworks is present, it is an error if one of the two
  *   representations is not there.
  *
  * - The dimension "NumberNodes" containing the number of nodes in the
  *   networks, and therefore the number of active power demand constraints
  *   for each time instants. The dimension is optional, if it is missing
  *   then no dual solution for the active power demand constraints is
  *   present.
  *
  * - The variable "ActivePowerDuals", of type netCDF::NcDouble and indexed
  *   over the dimensions "TimeHorizon" and "NumberNodes", in this order.
  *   This variable is only required to be present if "NumberNodes" is
  *   present, otherwise it is optional (since it is ignored).
  *   ActivePowerDuals[ t , n ] is the dual value \f$ y^{ac}_{t,n} \f$ of
  *   the node injection constraint (1) of node n at the time instant t (of
  *   (1b) if there is one node), i.e., minus the locational marginal price
  *   per instant (see UCBlock::get_Solution()).
  *
  * - The dimension "NumberPrimaryZones" tells how many "primary spinning
  *   reserve zones" are there in the problem. The dimension is optional, if
  *   it is not provided then it is taken to be 0, which means that no dual
  *   solution for the primary reserve constraints is present.
  *
  * - The variable "PrimaryDuals", of type netCDF::NcDouble and indexed over
  *   the dimensions "TimeHorizon" and "NumberPrimaryZones", in this order.
  *   This variable is only required to be present if "NumberPrimaryZones"
  *   is present, otherwise it is optional (since it is ignored). Entry
  *   PrimaryDuals[ t , z ] is the dual value \f$ y^{pr}_{t,z} \leq 0 \f$
  *   of the primary reserve constraint (2) of zone z at the time instant t.
  *
  * - The dimension "NumberSecondaryZones" tells how many "secondary spinning
  *   reserve zones" are there in the problem. The dimension is optional, if
  *   it is not provided then it is taken to be 0, which means that no dual
  *   solution for the secondary reserve constraints is present.
  *
  * - The variable "SecondaryDuals", of type netCDF::NcDouble and indexed
  *   over the dimensions "TimeHorizon" and "NumberSecondaryZones", in this
  *   order. This variable is only required to be present if
  *   "NumberSecondaryZones" is present, otherwise it is optional (since it
  *   is ignored). Entry SecondaryDuals[ t , z ] is the dual value
  *   \f$ y^{sc}_{t,z} \leq 0 \f$ of the secondary reserve constraint (3)
  *   of zone z at the time instant t.
  *
  * - The dimension "NumberInertiaZones" tells how many "inertia constraints
  *   zones" are there in the problem. The dimension is optional, if it is not
  *   provided then it is taken to be 0, which means that no dual solution for
  *   the inertia constraints is present.
  *
  * - The variable "InertiaDuals", of type netCDF::NcDouble and indexed over
  *   the dimensions "TimeHorizon" and "NumberInertiaZones", in this order.
  *   This variable is only required to be present if "NumberInertiaZones"
  *   is present, otherwise it is optional (since it is ignored). Entry
  *   InertiaDuals[ t , z ] is the dual value \f$ y^{in}_{t,z} \leq 0 \f$
  *   of the inertia constraint (4) of zone z at the time instant t.
  *
  * - The dimension "TotalNumberPollutantZones" tells how many pollutant
  *   budget constraints are there in the problem (cf.
  *   UCBlock::get_total_number_pollutant_zones()). The dimension is
  *   optional, if it is not provided then it is taken to be 0, which means
  *   that no dual solution for the pollutant budget constraints is present.
  *
  * - The variable "PollutantDuals", of type netCDF::NcDouble and indexed
  *   over the dimension "TotalNumberPollutantZones". This variable is only
  *   required to be present if "TotalNumberPollutantZones" is present,
  *   otherwise it is optional (since it is ignored). The entries are the
  *   dual values \f$ y^{p}_{z} \f$ of the pollutant budget constraints (5),
  *   ordered as in the PollutantBudget of the UCBlock: all the zones of
  *   pollutant 0, then all the zones of pollutant 1, and so on. */

 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// returns a copy of this UCBlockSolution with all values times \p factor
 UCBlockSolution * scale( double factor ) const override final;

 /// adds \p multiplier times \p solution to this UCBlockSolution
 void sum( const Solution * solution , double multiplier ) override final;

 /// returns a copy (an empty one if \p empty) of this UCBlockSolution
 UCBlockSolution * clone( bool empty = false ) const override final;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 /// prints the UCBlockSolution
 void print( std::ostream &output ) const override final {
  output << "UCBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 Index f_time_horizon;            ///< the time horizon
 Index f_number_nodes;            ///< the number of nodes
 Index f_number_primary_zones;    ///< the number of primary zones
 Index f_number_secondary_zones;  ///< the number of secondary zones
 Index f_number_inertia_zones;    ///< the number of inertia zones
 Index f_total_number_pollutant_zones;
 ///< the total number of pollutant zones

 bool f_compressed_network;
 ///< true if using the "compressed" format for NetworkBlock

 std::vector< UnitBlockSolution * > v_unit_Solution;
 ///< the Solution for each UnitBlock

 std::vector< NetworkBlockSolution * > v_network_Solution;
 ///< the Solution for each NetworkBlock

 boost::multi_array< double , 2 > v_demand_duals;
 ///< the dual variables for the node injection constraints
 /**< v_demand_duals[ t ][ n ] is the dual variable of the injection
  * constraint for node n at time t. */

 boost::multi_array< double , 2 > v_primary_duals;
 ///< the dual variables for the primary demand constraints
 /**< v_primary_duals[ t ][ z ] is the dual value of the primary reserve
  * constraint of zone z at time t. */

 boost::multi_array< double , 2 > v_secondary_duals;
 ///< the dual variables for the secondary demand constraints
 /**< v_secondary_duals[ t ][ z ] is the dual value of the secondary reserve
  * constraint of zone z at time t. */

 boost::multi_array< double , 2 > v_inertia_duals;
 ///< the dual variables for the inertia demand constraints
 /**< v_inertia_duals[ t ][ z ] is the dual value of the inertia constraint
  * of zone z at time t. */

 std::vector< double > v_pollutant_duals;
 ///< the dual variables for the pollutant budget constraints
 /**< The dual variables of the pollutant budget constraints of all the zones
  * of pollutant 0, then of all those of pollutant 1, and so on. */

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  ///< registration in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( UCBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __UCBlock */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File UCBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
