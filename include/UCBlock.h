/*--------------------------------------------------------------------------*/
/*--------------------------- File UCBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class UCBlock, which implements the Block concept [see
 * Block.h] for the Unit Commitment (UC) problem in electrical power
 * production. This is typically a short-term (across for instance one week or
 * one day time horizon) *deterministic* problem regarding finding an optimal
 * schedule of the production of electrical generators satisfying a (large)
 * set of technical constraints.
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
 #define __UCBlock    /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "UnitBlock.h"

#include "NetworkBlock.h"

#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

// TODO commented away until HeatBlock are properly managed
//class HeatBlock;     // forward declaration of HeatBlock

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS UCBlock --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Unit Commitment problem
/** The class UCBlock, implements the Block concept [see Block.h] for the
 * Unit Commitment (UC) problem in electrical power production. This is
 * typically a short-term (across for instance one week or one day time
 * horizon) *deterministic* problem regarding finding an optimal schedule
 * of the production of electrical generators satisfying a (large) set of
 * technical constraints.
 *
 * The model is quite flexible due to the fact that different types of units
 * and network constraints can be used by means of the fact that the class
 * manages son Block of type UnitBlock and NetworkBlock. Also UCBlock handles
 * a reasonably large variety of constraints, regarding not only active power
 * but also primary and secondary reserve and inertia. Admittedly, some
 * choices in UCBlock (like HeatBlock, pollution constraints, ...) are quite
 * specific of the UC of the plan4res project; however, all the "nonstandard"
 * aspects of UC can be switched away from the model (by simply not providing
 * the data describing them).
 *
 * The main elements that UCBlock handles are:
 *
 * - The time horizon of the problem, i.e., a discrete set of (typically,
 *   equally-spaced) time instants at which decisions are made (like, the
 *   24 hours in a day).
 *
 * - A set of electricity generating units, represented by derived classes
 *   of the base class UnitBlock.
 *
 * - A set of NetworkBlock, one for each time instant in the time horizon,
 *   which represent the constraints on the electricity demand satisfaction
 *   and the technical constraints on the transmission network. These can be
 *   basically "empty" if the capacity of the transmission network is such
 *   as to never really impact generation decisions (a "bus").
 *
 * - An optional set of HeatBlock, each representing the satisfaction of
 *   some specific "type of heat" on a close geographical area by
 *   heat-generating units possibly coupled with a heat storage. The link
 *   with the rest of the UC model lies in the fact that some of the
 *   heat-generating units in a HeatBlock may also be electricity generating
 *   ones (i.e., a UnitBlock); actually, the same UnitBlock can generate
 *   heat of "different types", and therefore appear as a heat-generating
 *   units in more than one HeatBlock.
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
 *   - possibly, constraints maximum pollutants emission for different kinds
 *     of pollutant, each "zone" (appropriately defined subset of the nodes
 *     of the transmission network) and for each time instant;
 *
 *   - possibly, constraints linking the electricity production of some
 *     UnitBlock with the heat production of some unit in a HeatBlock,
 *     for the appropriate units and for each time instant. */

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
  : Block( father ) , f_NetworkData( nullptr ) {}

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
 * - The dimension "NumberUnits" containing the number of units (UnitBlock) in
 *   the problem;
 *
 * - The dimension "NumberElectricalGenerators" containing the overall number
 *   of electrical generators for each UnitBlock of the problem;
 *
 * - The dimension "NumberNetworks" containing the number of networks
 *   (NetworkBlock) in the problem;
 *
 * - The groups "UnitBlock_0", "UnitBlock_1", ..., "UnitBlock_n" with n ==
 *   NumberUnits - 1, containing each one UnitBlock corresponding to one or
 *   more electricity generating unit (electrical generator). It is an error
 *   if the corresponding groups are not there. Each UnitBlock can have more
 *   than one electrical generator (cf. UnitBlock::get_number_generators());
 *   a value that is useful in the following (cf. "GeneratorNode") is the
 *   total number of those, which we will refer to as
 *   "NumberElectricalGenerators". This is computed by just calling
 *   get_number_generators() on each of the UnitBlock and summing all the
 *   results. Clearly, NumberElectricalGenerators >= NumberUnits; indeed,
 *   most of the UnitBlock can be expected to have just one electrical
 *   generator. If this happens for all the units then
 *   NumberElectricalGenerators == NumberUnits. If, instead, some UnitBlock
 *   (like cascades of hydro generators or combined cycle plants) actually
 *   has more than one electrical generator, then NumberElectricalGenerators >
 *   NumberUnits (each UnitBlock must have at least one). It is then useful
 *   (cf. "GeneratorNode") to be able to assign a unique index g = 0, 1, ...,
 *   NumberElectricalGenerators - 1 to each of the electrical generators in
 *   the UCBlock. When NumberElectricalGenerators == NumberUnits, the index
 *   is the same as i = 0, 1, ..., NumberUnits - 1 (there is a one-to-one
 *   correspondence between UnitBlock and electrical generators). When,
 *   instead, NumberElectricalGenerators > NumberUnits, a mapping must be
 *   defined. The mapping is the obvious one: UnitBlock have an ordering
 *   i = 0, 1, ..., NumberUnits - 1 (cf. the groups "UnitBlock_0",
 *   "UnitBlock_1", ... above), and the electrical generators into each
 *   UnitBlock also have some natural ordering (corresponding to the columns
 *   of the matrices of variables, cf. e.g., UnitBlock::get_commitment()).
 *   Thus, in general the mapping is:
 *     electrical generator 0 = first generator of UnitBlock_0
 *     electrical generator 1 = second generator of UnitBlock_0
 *     ...
 *     electrical generator k - 1 = k-th generator of UnitBlock_0
 *                          k     = UnitBlock_0->get_number_generators()
 *     electrical generator k     = first generator of UnitBlock_1
 *     electrical generator k + 1 = second generator of UnitBlock_1
 *     ...
 *   which of course boils down to "g = i" when each UnitBlock has exactly one
 *   electrical generator.
 *
 * - The dimension "NumberHeatBlocks" containing the number of heat blocks in
 *   the problem. The dimension is optional: if it is not provided then it is
 *   taken to be 0, which means that there is no heat block in the problem.
 *
 * - The groups "HeatBlock_0", "HeatBlock_1", ..., "HeatBlock_n" with
 *   n == NumberHeatBlocks - 1, containing each one a HeatBlock. When
 *   NumberHeatBlocks == 0, these groups need not be there since they are not
 *   read. If, instead, NumberHeatBlocks > 0, it is an error if the
 *   corresponding groups are not there. Since each HeatBlock can have more
 *   than one heat generator (cf. HeatBlock::get_number_heat_generators()),
 *   a value that is useful in the following (cf. "HeatSet") is the total
 *   number of those. We will refer to such number as "NumberHeatGenerators",
 *   which is computed by just calling get_number_heat_generators() on each
 *   of the HeatBlock and summing all the results. Clearly,
 *   NumberHeatGenerators >= NumberHeatBlock. Some of the HeatBlock may have
 *   just one heat generator; if this happens for all the heat blocks (but
 *   this is not likely), then NumberHeatGenerators ==  NumberHeatBlocks. It
 *   is then useful (cf. "HeatNode") to be able to assign a unique index
 *   h = 0, 1, ..., NumberHeatGenerators - 1 to each of the heat generators in
 *   the UCBlock. When NumberHeatGenerators == NumberHeatBlocks the index is
 *   the same as b = 0, 1, ..., NumberHeatBlock - 1 (there is a one-to-one
 *   correspondence between HeatBlock and heat generators, but this is not
 *   likely to happen). When, instead, NumberHeatGenerators > HeatBlock, a
 *   mapping must be defined. The mapping is the obvious one: HeatBlock have
 *   an ordering b = 0, 1, ..., NumberHeatBlocks - 1 (cf. the groups
 *   "HeatBlock_0", "HeatBlock_1", ... above), and the heat generators into
 *   each HeatBlock also have some natural ordering (corresponding to the
 *   columns of the matrices of variables, cf. e.g., HeatBlock::get_heat()).
 *   Thus, in general the mapping is:
 *     heat generator 0 = first generator of HeatBlock_0
 *     heat generator 1 = second generator of HeatBlock_0
 *     ...
 *     heat generator k = k-th generator of HeatBlock_0
 *                    k = HeatBlock_0->get_number_heat_units()
 *     heat generator k + 1 = first generator of HeatBlock_1
 *     heat generator k + 2 = second generator of HeatBlock_1
 *     ...
 *   which of course boils down to "h = b" when each HeatBlock has exactly
 *   one heat generator (but this is not assumed to happen).
 *
 * - Optionally, the dimensions and variables necessary to deserialize a
 *   NetworkData object that describes the transmission network; see
 *   NetworkBlock::NetworkData::deserialize() for details. If that is not
 *   provided (basically, "NumberNodes" is not provided or it is == 1), then
 *   the transmission network is taken to have only one node (a bus). If the
 *   data of a NetworkData is specified, the NetworkData is passed to each
 *   of the NetworkBlock (see below) of the UCBlock, if any. However, if
 *   a NetworkBlock also has a NetworkData specified in its own group, then
 *   the NetworkData inside the NetworkBlock overrules that inside the
 *   UCBlock, which is ignored by that NetworkBlock.
 *
 * - The variable "ActivePowerDemand", of type netCDF::NcDouble and indexed
 *   both over the dimensions "NumberNodes" and "TimeHorizon". This variable
 *   is optional if:
 *
 *   - a NetworkBlock is defined for each time instant (see below), and
 *
 *   - each of the defined NetworkBlock has the "ActiveDemand" variable
 *     specified in the corresponding group.
 *
 *   Otherwise it is mandatory. When it is defined, the entry
 *   ActivePowerDemand[ n , t ] is  assumed to contain the active power
 *   demand of node n of the transmission network at the given time
 *   instant t, where the first dimension "NumberNodes" can be read via
 *   NetworkBlock::NetworkData::get_number_nodes() from either the
 *   NetworkData object in UCBlock, or these in the NetworkBlock. When
 *   "ActivePowerDemand" is defined, and also the "ActiveDemand" variable
 *   is defined in the group of some NetworkBlock, then "ActiveDemand"
 *   overrules the value in the corresponding row of "ActivePowerDemand",
 *   which is ignored.
 *
 * - The groups "NetworkBlock_0", "NetworkBlock_1", ..., "NetworkBlock_T"
 *   with T = NumberNetworks - 1, with "NetworkBlock_t" containing the
 *   constraints on the transmission or community network at time t. The
 *   NetworkBlocks are optional, but if any of them are missing, then:
 *
 *   - "ActivePowerDemand" (see above) is mandatory in UCBlock;
 *
 *   - also the NetworkData (see above) is mandatory in UCBlock, unless
 *     the transmission network is a bus (that is, "NumberNodes" is not
 *     provided or it is == 1).
 *
 *   In particular, when a NetworkBlock is not defined for a given time
 *   instant t then a DCNetworkBlock is automatically constructed for that
 *   time instant, it is provided with the NetworkData object (which must be
 *   present in UCBlock) and the row ActivePowerDemand[ ... , t ] contains
 *   the active demand of each node at time instant t.
 *
 * - The variable "GeneratorNode", of type netCDF::NcUint and indexed over
 *   the set { 0 , ... , NumberElectricalGenerators - 1 }; GeneratorNode[ g ]
 *   tells to which node of the transmission network, the specified electrical
 *   generator g belongs. Note that this means that different electrical
 *   generators in the same UnitBlock can belong to different nodes of the
 *   transmission network. This is justified e.g., by hydro cascade units where
 *   different turbines can be rather far apart geographically, but still
 *   linked by (long) stretches of rivers. If NumberElectricalGenerators ==
 *   NumberUnits (all UnitBlock have exactly one electrical generator), then
 *   this variable is indexed over NumberUnits. If NumberNodes == 1 (say, it
 *   is not provided at all), then this variable need not be defined, since it
 *   is not loaded.
 *
 * - The variable "HeatNode", of type netCDF::NcUint and indexed over the
 *   dimension "NumberHeatBlocks"; the entry HeatNode[ h ] tells to which node
 *   of the transmission network *all* the heat generators that are also
 *   electrical generators in the HeatBlock h belong. Note that this implies
 *   that, unlike electrical generators in a UnitBlock, heat generators in a
 *   HeatBlock are always "geographically near to each other" so that they are
 *   necessarily attached to the same node of the transmission network (which
 *   is a reasonable assumption since heat, unlike say water, usually cannot
 *   travel much far). If NumberHeatBlocks == 0 (say, it is not provided at
 *   all), then this variable need not be defined, since it is not loaded.
 *   This information is actually only used to determine in which Pollutant
 *   Zone a HeatBlock is located, in order to add the corresponding heat
 *   generators (that are not also electricity generators) to the pollution
 *   constraints. This means that also if NumberPollutants == 0 (say, it is
 *   not provided at all) this variable is useless and therefore need not be
 *   defined, since it is not loaded. Finally, notice that for heat generators
 *   into a HeatBlock that also are electrical generators, this variable
 *   provides again an information that is already known, i.e., to which node
 *   they belong to (cf. "GeneratorNode" above). Of course *the two
 *   information must agree*, otherwise the input file is ill-defined and
 *   exception is thrown.
 *
 * - The variable "HeatSet", of type netCDF::NcUint and indexed over the set
 *   { 0 , ... , NumberHeatGenerators - 1 }. If HeatSet[ h ] = k <
 *   NumberElectricalGenerators, then k is the unique name of the electrical
 *   generator corresponding to the heat generator h; otherwise the heat
 *   generator h is not also an electrical generator. If NumberHeatBlocks == 0
 *   (there is no HeatBlock) then this variable need not be defined, since it
 *   is not loaded. It is possible that different heat generators correspond
 *   to the same electrical generator (that is, HeatSet[ h1 ] == HeatSet[ h2 ]
 *   for some h1 != h2), but not vice-versa: a heat generator either
 *   corresponds to a specific electrical generator, or to no electrical
 *   generator (it is an heat-only generator).
 *
 * - The variable "PowerHeatRho", of type netCDF::NcDouble and indexed over
 *   the dimension "NumberUnits": entry PowerHeatRho[ i ] is assumed to
 *   contain the electrical-power-to-heat ratio for *all generators* within
 *   unit i (most likely, a single generator). This is only used in the
 *   constraints linking electrical units to heat units, hence if
 *   NumberHeatBlocks == 0 (there are no HeatBlock) then this variable need
 *   not be defined, since it is not loaded.
 *
 * - The dimension "NumberPrimaryZones" tells how many "primary spinning
 *   reserve zones" are there in the problem. The dimension is optional, if it
 *   is not provided then it is taken to be 0, which means that no primary
 *   reserve constraints are present in the problem.
 *
 * - The variable "PrimaryZones", of type netCDF::NcUint and indexed over
 *   the dimension "NumberNodes". The entry PrimaryZones[ i ] tells to which
 *   primary zone the node i belongs: if PrimaryZones[ i ] >=
 *   NumberPrimaryZones, this means that node i does not belong to any primary
 *   zone, and hence the corresponding electrical generators are not involved
 *   into the primary reserve constraints. If NumberPrimaryZones == 0 (say, it
 *   is not provided at all) then this variable need not be defined, since it
 *   is not loaded. If NumberPrimaryZones == 1 and this variable is not
 *   defined, then there is only one primary zone and all the nodes belong to
 *   it.
 *
 * - The variable "PrimaryDemand", of type netCDF::NcDouble and indexed both
 *   over the dimensions "NumberPrimaryZones" and "TimeHorizon": entry
 *   PrimaryDemand[ i , t ] is assumed to contain the primary reserves
 *   requirement which are specified on the primary reserve zone i in the time
 *   t. If NumberPrimaryZones == 0 (say, it is not provided at all), then this
 *   variable need not be defined, since it is not loaded.
 *
 * - The dimension "NumberSecondaryZones" tells how many "secondary spinning
 *   reserve zones" are there in the problem. The dimension is optional, if it
 *   is not provided then it is taken to be 0, which means that no secondary
 *   reserve constraints are present in the problem.
 *
 * - The variable "SecondaryZones", of type netCDF::NcUint and indexed over
 *   the dimension "NumberNodes"; the entry SecondaryZones[ i ] tells to which
 *   secondary zone the node i belongs. If SecondaryZones[ i ] >=
 *   NumberSecondaryZones, this means that node i does not belong to any
 *   secondary zone, and hence the corresponding units are not involved into
 *   the secondary reserve constraints. If NumberSecondaryZones == 0 (say, it
 *   is not provided at all) then this variable need not be defined, since it
 *   is not loaded. If NumberSecondaryZones == 1 and this variable is not
 *   defined, then there is only one secondary zone and all the nodes belong
 *   to it.
 *
 * - The variable "SecondaryDemand", of type netCDF::NcDouble and indexed
 *   both over the dimensions "NumberSecondaryZones" and "TimeHorizon": entry
 *   SecondaryDemand[ i , t ] is assumed to contain the secondary reserve
 *   requirement which are specified on the secondary reserve zone i in the
 *   time t. If NumberSecondaryZones == 0 (say, it is not provided at all),
 *   then this variable need not be defined, since it is not loaded.
 *
 * - The dimension "NumberInertiaZones" tells how many "inertia constraints
 *   zones" are there in the problem. The dimension is optional, if it is not
 *   provided then it is taken to be 0, which means that no inertia
 *   constraints are present in the problem.
 *
 * - The variable "InertiaZones", of type netCDF::NcUint and indexed over
 *   the dimension "NumberNodes"; the entry InertiaZones[ n ] tells to which
 *   inertia zone the node n belongs. If InertiaZones[ n ] >=
 *   NumberInertiaZones, this means that node n does not belong to any inertia
 *   zone, and hence the corresponding units are not involved into the inertia
 *   reserve constraints. If NumberInertiaZones == 0 (say, it is not provided
 *   at all) then this variable need not be defined, since it is not
 *   loaded. If NumberInertiaZones == 1 and this variable is not defined, then
 *   there is only one inertia zone and all the nodes belong to it.
 *
 * - The variable "InertiaDemand", of type netCDF::NcDouble and indexed both
 *   over the dimensions "NumberInertiaZones" and "TimeHorizon": entry
 *   InertiaDemand[ i , t ] is assumed to contain the inertia reserves
 *   requirement which are specified on the inertia constraints zone i in the
 *   time t. If NumberInertiaZones == 0 (say, it is not provided at all), then
 *   this variable need not be defined, since it is not loaded.
 *
 * - The dimension "NumberPollutants" containing the number of pollutants in
 *   the problem. The dimension is optional, if it is not provided then it is
 *   taken to be 0, which means that no pollutants constraints are present in
 *   the problem.
 *
 * - The variable "NumberPollutantZones" of type netCDF::NcUint and indexed
 *   over the dimension "NumberPollutants": the entry NumberPollutantZones[ p
 *   ] is assumed to contain the number of pollutant zones associated with
 *   pollutant p. If NumberPollutants == 0 (say, it is not provided) then this
 *   variable need not be defined, since it is not loaded. The total number of
 *   pollutant zones is useful (cf. PollutantBudget); it will be referred to
 *   as "TotalNumberPollutantZones", and it is computed simply as
 *   TotalNumberPollutantZones = NumberPollutantZone[ 0 ] + ... +
 *   NumberPollutantZone[ NumberPollutants - 1 ].
 *
 * - The variable "PollutantZones", of type netCDF::NcUint and indexed over
 *   the dimensions "NumberPollutants" and "NumberNodes": the entry
 *   PollutantZones[ p , n ] tells to which pollutant zone associated with
 *   pollutant p the node n belongs. If PollutantZones[ p , n ] >=
 *   NumberPollutantZones[ p ], this means that node n does not belong to any
 *   pollutant zone, and hence the corresponding units are not involved into
 *   the pollutant budget constraints associated with pollutant p. If
 *   NumberPollutants == 0 (say, it is not provided) then this variable need
 *   not be defined, since it is not loaded.
 *
 * - The variable "PollutantBudget", of type netCDF::NcDouble and indexed
 *   over the set { 0 , ... , TotalNumberPollutantZones - 1 }: the entry
 *   PollutantBudget[ n ] for n == 0 , ..., TotalNumberPollutantZones - 1 is
 *   assumed to contain the pollutant budget (across all the time horizon) for
 *   the pair (zone of the pollutant , pollutant) corresponding to n. In
 *   another word, since the number of pollutant zones of each pollutant may
 *   not be equal with each other it is useful (to avoid having to store
 *   PollutantBudget as an irregular matrix) to be able to assign a unique
 *   index n = 0, 1, ..., TotalNumberPollutantZones - 1 to each pollutant
 *   budget of each pollutant zone. A mapping must be defined between each
 *   entry n and the pair (zone of the pollutant , pollutant). The mapping is
 *   the obvious one: UCBlock has a set of pollutants p = 0, 1, ...,
 *   NumberPollutants - 1, and each pollutant p may have several pollutant
 *   zones (see comments of variable "NumberPollutantZones" above). Thus, in
 *   general the mapping is: n = 0 corresponds to the zone 0 of pollutant 0 n
 *   = 1 corresponds to the zone 1 of pollutant 0 ... n =
 *   NumberPollutantZone[ 0 ] - 1 corresponds to the zone NumberPollutantZone[
 *   0 ] - 1 of pollutant 0 n = NumberPollutantZone[ 0 ] corresponds to the
 *   zone 0 of pollutant 1 n = NumberPollutantZone[ 0 ] + 1 corresponds to the
 *   zone 1 of pollutant 1 ... . If NumberPollutants == 0 (say, it is not
 *   provided) then this variable need not be defined, since it is not loaded.
 *
 * - The variable "PollutantRho", of type netCDF::NcDouble and indexed over
 *   three dimensions which are "TimeHorizon" and "NumberPollutants" and the
 *   set { 0 , ... , NumberElectricalGenerators - 1 } (see comments above). The
 *   first dimension can have size either 1 or TimeHorizon. In the former case
 *   the entry PollutantRho[ 0 , p , g ] is assumed to contain the conversion
 *   factor of pollutant p due to the electrical generator g which is equal
 *   for all time instants t. Otherwise, the first dimension has full size
 *   TimeHorizon and the entry PollutantRho[ t , p , g ] gives the conversion
 *   factor of pollutant p due to the electrical generator g for time t. If
 *   NumberPollutants == 0 (it is not provided) then this variable need not be
 *   defined, since it's not loaded.
 *
 * - The variable "StartNetworkIntervals", of type netCDF::NcUint and indexed
 *   over the dimension "NumberNetworks"; the entry StartNetworkIntervals[ n ]
 *   tells at which time horizon w.r.t. the given dimension "TimeHorizon"
 *   starts the NetworkBlock n, and the difference StartNetworkIntervals[ n +
 *   1 ] - StartNetworkIntervals[ n ] indicates the number sub time horizons
 *   spanned by the NetworkBlock n. If it is not given in input, each
 *   NetworkBlock is supposed to span one single time horizon, so
 *   StartNetworkIntervals[ n ] is set equal to n.
 *
 * - The variable "NetworkConstantTerms", of type netCDF::NcDouble and
 *   indexed over the dimension "NumberNetworks"; the entry
 *   NetworkConstantTerms[ n ] tells the constant term, i.e., typically the
 *   fixed costs, of the NetworkBlock n.
 *
 * - The variable "NetworkBlockClassname", of type netCDF::NcString specify
 *   the classname of the specific NetworkBlock to be instantiate if no one is
 *   explicitly given in input. For backward compatibility reasons w.r.t. the
 *   netCDF input data files already given, the default value is
 *   "DCNetworkBlock".
 *
 * - The variable "NetworkDataClassname", of type netCDF::NcString specify
 *   the classname of the specific NetworkData to be instantiate. For backward
 *   compatibility reasons w.r.t. the netCDF input data files already given,
 *   the default value is "DCNetworkData". */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generates the static constraint of the UCBlock
 /** This method generates the abstract constraints of the UCBlock.
  *
  * Consider a network defined by a set of nodes \f$ \mathcal{N} \f$ and a set
  * of arcs connecting the nodes \f$ \mathcal{L} \f$. There are moreover given
  * three partitions of the set of nodes which may or may not be identical:
  *
  * (i). \f$ \mathcal{B}^{pr}(\mathcal{N}) \f$ partitions \f$ \mathcal{N} \f$
  * in several zones (sets of nodes) each one being associated with one
  * specific primary spinning reserve requirement;
  *
  * (ii). \f$ \mathcal{B}^{sc}(\mathcal{N}) \f$ partitions \f$ \mathcal{N}
  * \f$ in several zones each one being associated with one specific secondary
  * spinning reserve requirement;
  *
  * (iii). \f$ \mathcal{B}^{in}(\mathcal{N}) \f$ partitions \f$ \mathcal{N}
  * \f$ in several zones each one being associated with one specific inertia
  * requirement;
  *
  * Optionally we are given a partition of the nodes \f$ B^{p}(\mathcal{N})
  * \f$ corresponding to zones which are associated with an emissions
  * constraint on the specific pollutant \f$ p \f$ in the set of pollutants
  * \f$ \mathcal{P} \f$.
  *
  * The electrical system contains a set of “units” (e.g., power plants, load
  * flexibilities or storage devices) indexed by \f$ i \in \mathcal{I} \f$.
  * The set \f$ \mathcal{I}_n \f$ will indicate units connected to node \f$ n
  * \in \mathcal{N} \f$. Each unit contains one (or more) electrical generator
  * where the set of all electrical generators is defined by
  * \f$ g \in \mathcal{G} \f$ and the set \f$ \mathcal{G}_n \f$ will indicate
  * electrical generators connected to node \f$ n \in \mathcal{N} \f$.
  *
  * In the UCBlock there is not defined any variable but the decision
  * variables here are present by using method get_variable() from UnitBlock,
  * HeatBlock and NetworkBlock as follow:
  *
  * - \f$ p^{ac}_{t,g} \f$ : the active power variable for each time period
  *   \f$ t \in \mathcal{T} \f$ and each electrical generator
  *   \f$ g \in \mathcal{G} \f$ is called from UnitBlock by get_active_power()
  *   method;
  *
  * - \f$ S_{t,n} \f$ : the node injection variable for each time period
  *   \f$ t \in \mathcal{T} \f$ and each node \f$ n \in \mathcal{N} \f$ is
  *   called from NetworkBlock by get_node_injection() method;
  *
  * - \f$ p^{pr}_{t,g} \f$ : the primary spinning reserves variable for each
  *   time period \f$ t \in \mathcal{T} \f$ and each electrical generator
  *   \f$ g \in \mathcal{G} \f$ is called from UnitBlock by
  *   get_primary_spinning_reserve() method;
  *
  * - \f$ p^{sc}_{t,g} \f$ : the secondary spinning reserves variable for each
  *   time period \f$ t \in \mathcal{T} \f$ and each electrical generator
  *   \f$ g \in \mathcal{G} \f$ is called from UnitBlock by
  *   get_secondary_spinning_reserve() method;
  *
  * - \f$ u_{t,g}  \in \{ 0 , 1 \} \f$ : the commitment state at time period
  *   \f$ t \in \mathcal{T} \f$ and each electrical generator
  *   \f$ g \in \mathcal{G} \f$ is called from UnitBlock by get_commitment()
  *   method;
  *
  * - \f$ p^{he}_{t,i} \f$ : the heat variable for each time period
  *   \f$ t \in \mathcal{T} \f$ and each heat unit
  *   \f$ i \in \mathcal{I}(h)= \mathcal{I}^{ho}(h) \cup \mathcal{I}^{ec}(h)\f$
  *   is called from HeatBlock by get_heat() method, where
  *   \f$\mathcal{I}^{ho}(h)\f$  and \f$\mathcal{I}^{ec}(h)\f$ are heat-only-
  *   producing unit and electricity-producing one, respectively;
  *
  * The global constraints of unit commitment problem, on the time horizon
  * \f$ \mathcal{T} \f$ write as follow:
  *
  * - Node injection Constraints:
  *   In the unit commitment problem, \f$ P^{au}_{t , g} \f$ denotes the fixed
  *   consumption of the power plant when it is off, and \f$ S_{t,n} \f$ is the
  *   node injection variable for each time period \f$ t \in \mathcal{T} \f$
  *   and each node \f$ n \in \mathcal{N} \f$. Therefore, if the NetworkBlock::
  *   get_number_nodes() > 0, a boost::multi_array< FRowConstraint , 2 >; with
  *   two dimensions which are get_time_horizon() and
  *   NetworkBlock::get_number_nodes() entries, where the entry
  *   t = 0, ..., f_time_horizon - 1 and the entry
  *   n = 1, ..., get_number_nodes() being the node injection constraints at
  *   time t and node n as follow:
  *
  *   \f[
  *    \sum_{ g \in \mathcal{G}_n } (p^{ac}_{t,g} + P^{au}_{t , g}(1 - u_{t,g}))
  *     = S_{t,n}     \quad t \in \mathcal{T} \quad n \in \mathcal{N} \quad (1)
  *   \f]
  *
  * - Primary Demand Constraints:
  *   In the unit commitment problem, the primary demand
  *   \f$ D^{pr}_{\mathcal{B} , t} \f$ which are specified on the primary
  *   reserve zones \f$ \mathcal{B} \in \mathcal{B}^{pr}(\mathcal{N}) \f$ will
  *   be satisfied. Therefore, if the f_number_primary_zones > 0,
  *   a boost::multi_array< FRowConstraint , 2 >; with two dimensions which are
  *   f_time_horizon, and f_number_primary_zones entries, where the entry
  *   t = 0, ..., f_time_horizon - 1 and the entry
  *   \f$ \mathcal{B}\f$ = 0, ..., f_number_primary_zones - 1 being the primary
  *   demand constraints at time t and primary zone \f$ \mathcal{B}\f$ as below;
  *
  *   \f[
  *    \sum_{n \in \mathcal{B}}\sum_{ g \in \mathcal{G}_n } p^{pr}_{t,g} \geq
  *     D^{pr}_{\mathcal{B} , t} \quad t \in \mathcal{T}
  *      \quad \mathcal{B} \in \mathcal{B}^{pr}(\mathcal{N})          \quad (2)
  *   \f]
  *
  * - Secondary Demand Constraints:
  *   In the unit commitment problem, the secondary demand
  *   \f$ D^{sc}_{\mathcal{B} , t} \f$ which are specified on the secondary
  *   reserve zones \f$ \mathcal{B} \in \mathcal{B}^{sc}(\mathcal{B}) \f$ will
  *   be satisfied. So, if the f_number_secondary_zones > 0,
  *   a boost::multi_array< FRowConstraint , 2 >; with two dimensions which are
  *   f_time_horizon, and f_number_secondary_zones entries, which the entry
  *   t = 0, ..., f_time_horizon - 1 and the entry
  *   \f$ \mathcal{B}\f$ = 0, ..., f_number_secondary_zones - 1 being the
  *   secondary demand constraints at time t and secondary zone
  *   \f$ \mathcal{B}\f$ as follow;
  *
  *   \f[
  *    \sum_{n \in \mathcal{B}}\sum_{ g \in \mathcal{G}_n } p^{sc}_{t,g} \geq
  *       D^{sc}_{\mathcal{B} , t} \quad t \in \mathcal{T}
  *       \quad \mathcal{B} \in \mathcal{B}^{sc}(\mathcal{N})         \quad (3)
  *   \f]
  *
  * - Inertia Demand Constraints:
  *   In the unit commitment problem, the inertia demand
  *   \f$ D^{in}_{\mathcal{B} , t}\f$ which are specified on the inertia zones
  *   \f$ \mathcal{B} \in \mathcal{B}^{in}(\mathcal{N}) \f$ with defined
  *   parameters \f$ \alpha_{t , g} \f$ and \f$ \beta_{t , g} \f$ will be
  *   satisfied. Therefore, if the f_number_inertia_zones > 0,
  *   a boost::multi_array< FRowConstraint , 2 >; with two dimensions which are
  *   f_time_horizon, and f_number_inertia_zones entries, where the entry
  *   t = 0, ..., f_time_horizon - 1 and the entry
  *   \f$ \mathcal{B}\f$ = 0, ..., f_number_inertia_zones - 1 being the inertia
  *   demand constraints at time t and inertia zone \f$ \mathcal{B}\f$ as below;
  *
  *   \f[
  *    \sum_{n \in \mathcal{B}}\sum_{ g \in \mathcal{G}_n } (\alpha_{t , g}
  *     u_{t,g} + \beta_{t , g} p^{ac}_{t,g}) \geq D^{in}_{\mathcal{B} , t}
  *        \quad t \in \mathcal{T}
  *        \quad \mathcal{B} \in \mathcal{B}^{in}(\mathcal{N})        \quad (4)
  *   \f]
  *
  * - Pollutant Budget Constraints:
  *   In the unit commitment problem, the pollutant budget
  *   \f$ \mathcal{O}_{\mathcal{B},p} \f$ which is specified on each pollutant
  *   \f$ p \in \mathcal{P} \f$ in each pollutant zone
  *   \f$ \mathcal{B} \in \mathcal{B}^{p}(\mathcal{N}) \f$  with two parameters
  *   \f$ \rho_{t , p , g} \f$ and \f$ \gamma_{t , p , h} \f$ where considered
  *   as pollutant ratio and pollutant heat ratio respectively. So, if the
  *   f_number_pollutants > 0, a std::vector< std::vector< FRowConstraint > >;
  *   with two dimensions which are f_number_pollutants and
  *   v_number_pollutant_zones entries, that the entry
  *   p = 0, ..., f_number_pollutants - 1 and the entry
  *   \f$ \mathcal{B}\f$ = 0, ..., v_number_pollutant_zones - 1 being the
  *   pollutant budget constraints at pollutant \f$ \mathcal{B}\f$ and pollutant
  *   zones b as below;
  *
  *   \f[
  *    \sum_{n \in \mathcal{B}}\sum_{ t \in \mathcal{T} }( \sum_{ g \in
  *    \mathcal{G}_n } \rho_{t , p , g} p^{ac}_{t,g} + \sum_{h \in \mathcal{H}_n}
  *    \sum_{ j \in \mathcal{I}^{ho}(h)} \gamma_{t , p , h} p^{h,he}_{t,j} )
  *    \leq \mathcal{O}_{\mathcal{B},p}  \quad \mathcal{B} \in
  *    \mathcal{B}^{p}(\mathcal{N}) \quad p \in \mathcal{P}           \quad (5)
  *   \f]
  *
  *   where \f$ \mathcal{H} \f$ is the set of Heat Blocks.
  *
  * - Heat Constraints:
  *   In the unit commitment problem, the Heat Constraints link the UCBlock
  *   variables with the HeatBlock, where for each heat block
  *   \f$ h \in \mathcal{H} \f$ and each electrical-power-to-heat ratio \f$
  *   \varrho_{g} \f$ of each electrical generator \f$ g \in \mathcal{G} \f$.
  *   Therefore, if the f_number_heat_blocks > 0, a
  *   boost::multi_array< FRowConstraint , 2 > with two dimensions which are
  *   f_time_horizon and the number of electrical generators that belong to
  *   some HeatBlock; the constraint at position ( t, g ) being the heat
  *   constraints at time t and heat generator M[ g ], where M maps the
  *   constraint into an electricity generator that belongs to some HeatBlock.
  *   The Heat Constraints are defined as below:
  *
  *   \f[
  *    \sum_{h \in \mathcal{H} , j \in \mathcal{G}^{ec}(h): e^h(j)=g}
  *     p^{h , he}_{t , j}  \leq \varrho_g p^{ac}_{t,g} \quad g \in \mathcal{G}
  *                                 \quad t \in \mathcal{T}           \quad (6)
  *   \f]
  *
  *   where \f$ j \in \mathcal{G}^{ec}(h) \f$ is an electricity generator in a
  *   heat block \f$ h \in \mathcal{H} \f$. For \f$ j \in
  *   \mathcal{G}^{ec}(h) \f$, there is the need to know which electrical
  *   generator \f$ j \f$ is representing. Thus, we need a mapping
  *   \f$ e^h : \mathcal{G}^{ec}(h) \to \mathcal{G} \f$, where \f$ \mathcal{G}
  *   \f$ is the set of electricity generators (standard electrical generators
  *   in UC parlance). */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the related blocks to the UC problem.
 * @{ */

 /// returns the sense of the Objective of this UCBlock
 /** This function returns the sense of the Objective of this UCBlock, which
  * is currently defined to be minimization (Objective::eMin). */

 int get_objective_sense( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the time horizon of the problem

 Index get_time_horizon( void ) const { return( f_time_horizon ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of UnitBlock

 Index get_number_units( void ) const { return( f_number_units ); }

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
 /** Since there always is a NetworkBlock for each time instant t, this vector
  * should have size of get_time_horizon() where each element of the vector
  * gives the network block at time instant t. */

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
 /** The method returned a std::vector< Index > V and each element of V tells
  * to which primary zone node n belongs. There are three possible cases:
  *
  * - if V is empty, then no primary zones are defined, and there are no
  *   primary reserve constraints;
  *
  * - if V has only one element, then the transmission network is a bus and
  *   that unique node belongs to the primary zone;
  *
  * - otherwise, V must have size NetworkBlock::get_number_nodes() and V[ n ]
  *   tells to which primary zone the node n belongs; if
  *   V[ n ] >= get_number_primary_zones(), then node n does not belong to
  *   any primary zone, and hence the corresponding electrical generators are
  *   not involved into the primary reserve constraints. */

 const std::vector< Index > & get_primary_zone( void ) const {
  return( v_primary_zones );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary zones
 /** The method returned a std::vector< Index > V and each element of V tells
  * to which secondary zone node n belongs. There are three possible cases:
  *
  * - if V is empty, then no secondary zones are defined, and there are no
  *   primary reserve constraints;
  *
  * - if V has only one element, then the transmission network is a bus and
  *   that unique node belongs to the secondary zone;
  *
  * - otherwise, V must have size NetworkBlock::get_number_nodes() and V[ n ]
  *   tells to which secondary zone the node n belongs; if
  *   V[ n ] >= get_number_secondary_zones(), then node n does not belong to
  *   any secondary zone, and hence the corresponding electrical generators
  *   are not involved into the secondary reserve constraints. */

 const std::vector< Index > & get_secondary_zone( void ) const {
  return( v_secondary_zones );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inertia zones
 /** The method returned a std::vector< Index > V and each element of V tells
  * to which inertia zone node n belongs. There are three possible cases:
  *
  * - if V is empty, then no inertia zones are defined, and there are no
  *   inertia reserve constraints;
  *
  * - if V has only one element, then the transmission network is a bus and
  *   that unique node belongs to the inertia zone;
  *
  * - otherwise, V must have size NetworkBlock::get_number_nodes() and V[ n ]
  *   tells to which inertia zone the node n belongs; if
  *   V[ n ] >= get_number_inertia_zones(), then node n does not belong to
  *   any inertia zone, and hence the corresponding electrical generators
  *   are not involved into the inertia reserve constraints. */

 const std::vector< Index > & get_inertia_zone( void ) const {
  return( v_inertia_zones );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of primary demand
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the primary demand of the primary zone n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no primary zones are
  *   defined, and there are no primary reserve constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_time_horizon() and it
  *   means there exist just one primary zone in the problem where the
  *   node(s) belongs to that primary zone. Each element of this vector gives
  *   the primary demand of the unique primary zone at time instant t;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_primary_zones() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] gives the primary
  *   demand of primary zone n at time instant t. */

 const boost::multi_array< double , 2 > & get_primary_demand( void ) const {
  return( v_primary_demand );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of secondary demand
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the secondary demand of the secondary zone n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no secondary zones are
  *   defined, and there are no secondary reserve constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_time_horizon() and it
  *   means there exist just one secondary zone in the problem where the
  *   node(s) belongs to that secondary zone. Each element of this vector
  *   gives the secondary demand of the unique secondary zone at time instant
  *   t;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_secondary_zones() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] gives the secondary
  *   demand of secondary zone n at time instant t. */

 const boost::multi_array< double , 2 > & get_secondary_demand( void ) const {
  return( v_secondary_demand );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inertia demand
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the inertia demand of the inertia zone n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no inertia zones are
  *   defined, and there are no inertia reserve constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size of f_time_horizon and it
  *   means there exist just one inertia zone in the problem where the
  *   node(s) belongs to that inertia zone. Each element of this vector gives
  *   the inertia demand of the unique inertia zone at time instant t;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_inertia_zones() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] gives the inertia
  *   demand of inertia zone n at time instant t. */

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
 /// returns the matrix of pollutant zones
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ p , n ] tells to which pollutant zone associated with pollutant p the
  * node n belongs. There are four possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no pollutant zones are
  *   defined, and there are no pollutant budget constraints;
  *
  * - if the boost::multi_array<> M has only one row, it is a vector
  *   with size NetworkBlock::get_number_nodes(). In this case only one
  *   pollutant zone exists in the problem, and the nodes may belong (or not)
  *   to that pollutant zone. Therefore, each n_th element of the vector
  *   tells if the node n belongs to the unique pollutant zone or not;
  *
  * - if the boost::multi_array<> M has only one element, the transmission
  *   network is a bus with one pollutant zone;
  *
  * - otherwise, the two-dimensional boost::multi_array<> M must have
  *   get_number_pollutants() rows and NetworkBlock::get_number_nodes()
  *   columns, and each element of matrix M[ p , n ] tells to which pollutant
  *   zone associated with pollutant p the node n belongs. */

 const boost::multi_array< Index , 2 > & get_pollutant_zone( void ) const {
  return( v_pollutant_zones );
 }

/*--------------------------------------------------------------------------*/
 /// returns the two-dimensional vector of pollutant budget
 /** The method returned a std::vector< std::vector< double > > V such that
  * V[ p ] contains the pollutant budget (across all the time horizon) for
  * all the pollutant zone of the pollutant p. There are two possible cases:
  *
  * - If V is empty() then no pollutant zones are defined, and there are no
  *   pollutant budget constraints (in this case, get_number_pollutants()
  *   must return zero).
  *
  * - Otherwise, V.size() == get_number_pollutants(). For each pollutant
  *   p = 0, ..., get_number_pollutants() - 1, V[ p ].size() is the number
  *   of different pollutant areas for p, and V[ p ][ z ] is the pollutant
  *   budget (across all the time horizon) for zone z of pollutant p. Note
  *   that, therefore, get_pollutant_zone()[ p ][ i ] is either a number <
  *   V[ p ].size(), which means that node i belongs to one particular
  *   pollutant zone, or get_pollutant_zone()[ p ][ i ] >= V[ p ].size(),
  *   which means that node i does not belong to any pollutant zone for
  *   pollutant p. */

 const std::vector< std::vector< double >> &
 get_pollutant_budget( void ) const {
  return( v_pollutant_budget );
 }

/*--------------------------------------------------------------------------*/
 // TODO commented away until HeatBlock are properly managed

 // /// returns the vector of HeatSet
 // /** The method returns a std::vector< Index > V such that each element
 //  * implies which heat generator is also an electrical generator. There are
 //  * three possible cases:
 //  *
 //  * - if V is empty, then no HeatSet is defined, and there is no heat linking
 //  *   constraints;
 //  *
 //  * - if V has only one element, then there is just one heat generator in
 //  *   heat blocks which is also an electrical generator;
 //  *
 //  * - otherwise, V.size() == NumberHeatGenerators (see the comments to
 //  *   deserialize()), and each element of V[ h ] tells which heat generator
 //  *   is also an electrical generator. */

 // const std::vector< Index > & get_heat_set( void ) const {
 //  return( v_heat_set );
 // }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of pollutant rho
 /** The method returned a three-dimensional boost::multi_array<> M such that
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

/*--------------------------------------------------------------------------*/
 // TODO commented away until HeatBlock are properly managed

 // /// returns the matrix of pollutant heat rho
 // /** The method returned a three-dimensional boost::multi_array<> M such that
 //  * M[ t , p , i ] gives the conversion factor of the given pollutant p due to
 //  * the generation of every heat-only unit i in the given heat block h at the
 //  * given time t. This three-dimensional boost::multi_array<> M considers two
 //  * possible cases:
 //  *
 //  * - if the boost::multi_array<> M is empty() then two possible cases are:
 //  *
 //  *   - there is no pollutant zone, hence there are not defined any pollutant
 //  *     budget constraints;
 //  *
 //  *   - there is no HeatBlock, hence in pollutant budget constraints there
 //  *     is not heat-rho-linking part;
 //  *
 //  * - otherwise, two possible cases may happen to the first dimension of the
 //  *   M[ t , p , i ];
 //  *
 //  *   - if the first dimension of the boost::multi_array<> M has size one,
 //  *     then each element of the matrix M [ 0 , p , i ] gives the conversion
 //  *     factor of pollutant p due to the of every heat-only unit i in the
 //  *     given heat block h;
 //  *
 //  *   - if the first dimension of the boost::multi_array<> M has full size
 //  *     then, each element of the matrix M[ t , p , i ] gives the conversion
 //  *     factor of pollutant p due to the generation of every heat-only unit i
 //  *     in the given heat block h for time instant t. */

 // const boost::multi_array< double, 3 > & get_pollutant_heat_rho( void )
 // const {
 //  return( v_pollutant_heat_rho );
 // }

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
 // TODO commented away until HeatBlock are properly managed

 // /// returns the vector of (pointers to) HeatBlock elements.
 // /** The vector of heat blocks in the problem. There are three possible cases:
 //  *
 //  * - if the vector is empty, then the there is no heat block;
 //  *
 //  * - if the vector only has one element, then there is just one heat block in
 //  *   the problem;
 //  *
 //  * - otherwise the vector must have the size of the number of heat blocks,
 //  *   and the h-th entry gives the corresponding heat block h. */

 // const std::vector< HeatBlock * > & get_heat_block( void ) const {
 //  return( v_heat_blocks );
 // }

/*--------------------------------------------------------------------------*/
 /// returns the vector of generator node
 /** This method returns a vector V that indicates to which node of the
  * transmission network each electrical generator belongs. There are two
  * possible cases:
  *
  * - if V has only one element, then the transmission network is a bus and
  *   all the units (electrical generators) belong to that unique node;
  *
  * - otherwise, the size of the vector V must be the number of units and
  *   V[ g ] indicates to which node the electrical generator g belongs. */

 const std::vector< Index > & get_generator_node( void ) const {
  return( v_generator_node );
 }

/*--------------------------------------------------------------------------*/
 // TODO commented away until HeatBlock are properly managed

 // /// returns the vector of electrical-power-to-heat ratio
 // /** The method returned a std::vector< double > V and each element of V
 //  * contains the electrical-power-to-heat ratio for each unit i of heat block
 //  * h. There are three possible cases:
 //  *
 //  * - if V is empty, then the there is no heat block;
 //  *
 //  * - if V only has one element, then the power heat rho is always equal to
 //  *   the value of that element;
 //  *
 //  * - otherwise the vector must have the size of the number of units, and the
 //  *   V[ i ] gives the power heat rho for each heat block h. */

 // const std::vector< double > & get_power_heat_rho( void ) const {
 //  return( v_power_heat_rho );
 // }

/*--------------------------------------------------------------------------*/
 /// returns the node injection constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the node injection constraints. C[ t ][ n ] is the node
  * injection constraint associated with time t and node n. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_node_injection_constraints( void ) const {
  return( v_node_injection_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the primary demand constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the primary demand constraints. C[ t ][ z ] is the primary
  * demand constraint associated with time t and primary zone z. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_primary_demand_constraints( void ) const {
  return( v_PrimaryDemand_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the secondary demand constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the secondary demand constraints. C[ t ][ z ] is the secondary
  * demand constraint associated with time t and secondary zone z. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_secondary_demand_constraints( void ) {
  return( v_SecondaryDemand_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the inertia demand constraints
 /** This method returns (a const reference to) the boost multi_array C
  * containing the inertia demand constraints. C[ t ][ z ] is the inertia
  * demand constraint associated with time t and inertia zone z. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_inertia_demand_constraints( void ) const {
  return( v_InertiaDemand_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum pollutant emission constraints
 /** This method returns (a const reference to) the vector C containing the
  * maximum pollutant emission constraints. C[ p ][ z ] is the maximum
  * pollutant emission constraint associated with pollutant p and pollutant
  * zone z. */

 const std::vector< std::vector< FRowConstraint > > &
 get_pollutant_constraints( void ) const {
  return( v_PollutantBudget_Const );
 }

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR OBTAINING INFORMATION ABOUT THE UCBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for obtaining information about the UCBlock
 * @{ */

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
  if( ( f_number_primary_zones > 1 ) &&
      ( zone_id != v_primary_zones[ node_id ] ) )
   return( false );
  return( true );
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
  if( ( f_number_secondary_zones > 1 ) &&
      ( zone_id != v_secondary_zones[ node_id ] ) )
   return( false );
  return( true );
 }

/*--------------------------------------------------------------------------*/
 /// returns true if the given node belongs to the given inertia zone
 /** This function returns true if and only if the node identified by \p
  * node_id belongs to the inertia zone identified by \p zone_id.
  *
  * @param node_id The ID of a node.
  *
  * @param zone_id The ID of a inertia zone.
  *
  * @return True if and only if the given node belongs to the given inertia
  *         zone. */

 bool node_belongs_to_inertia_zone( Index node_id , Index zone_id ) const {
  if( ( f_number_inertia_zones > 1 ) &&
      ( zone_id != v_inertia_zones[ node_id ] ) )
   return( false );
  return( true );
 }

/*--------------------------------------------------------------------------*/
 /// returns true if the given electrical generator belongs to the given node

 bool generator_belongs_to_node( Index elc_generator , Index node_id ) const {
  if( ( get_number_nodes() > 1 ) &&
      ( node_id != v_generator_node[ elc_generator ] ) )
   return( false );
  return( true );
 }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UCBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing and saving the UCBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * UCBlock. See UCBlock::deserialize( netCDF::NcGroup ) for details of the
  * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UCBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UCBlock
    @{ */

 /**
  * @brief It loads a UCBlock from a input standard stream.
  * @warning This method is not implemented yet.
  * @param input an input stream
  */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "UCBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 /// method for handling Modification
 /** Method for handling Modification.
  *
  * This method has to intercept any "abstract Modification" that modifies the
  * "abstract representation" of the UCBlock, and "translate" them into both
  * changes of the actual data structures and corresponding "physical
  * Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Currently, this method only
  * handles UnitBlockMod whose modification is associated with the changing of
  * the scale factor of a UnitBlock. */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// update the active power demand

 void set_active_power_demand( MF_dbl_it values ,
                               Subset && subset = { 0 } ,
                               const bool ordered = false ,
                               ModParam issuePMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the active power demand

 void set_active_power_demand( MF_dbl_it values ,
                               Range rng = Range( 0 , 1 ) ,
                               ModParam issuePMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// states that the Variable of the UCBlock have been generated
 void set_variables_generated( void ) { AR |= HasVar; }

 /// states that the Constraint of the UCBlock have been generated
 void set_constraints_generated( void ) { AR |= HasCst; }

 /// states that the Objective of the UCBlock has been generated
 void set_objective_generated( void ) { AR |= HasObj; }

 /// indicates whether the Variable of the UCBlock have been generated
 bool variables_generated( void ) const { return( AR & HasVar ); }

 /// indicates whether the Constraint of the UCBlock have been generated
 bool constraints_generated( void ) const { return( AR & HasCst ); }

 /// indicates whether the Objective of the UCBlock has been generated
 bool objective_generated( void ) const { return( AR & HasObj ); }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the time horizon of the problem
 Index f_time_horizon{};

 /// the specific classname of the networks that need to be instantiated,
 /// e.g., `DCNetworkBlock`, `ECNetworkBlock`, ecc.
 // Used in case of no `NetworkBlock`s were given since there is just one
 // node, i.e., the network is a bus, or all the nodes share the same data in
 // `NetworkData`.
 std::string network_block_classname;
 std::string network_data_classname;

 /// the number of the networks of the problem
 Index f_number_networks{};

 /// the number of units of the problem
 Index f_number_units{};

 /// the number of electrical generators of the problem
 Index f_number_elc_generators{};

 /* TODO commented away until HeatBlock are properly managed
 /// the number of heat generators of the problem
 Index f_number_heat_generators;
 */

 /// the total number of pollutant zones of the problem
 Index f_total_number_pollutant_zones{};

 /// the NetworkData object
 NetworkBlock::NetworkData * f_NetworkData;

 /* TODO commented away until HeatBlock are properly managed
 /// the number of heat block
 Index f_number_heat_blocks;
 */

 /// the number of nodes in primary zones of the network
 Index f_number_primary_zones{};

 /// the number of nodes in secondary zones of the network
 Index f_number_secondary_zones{};

 /// the number of nodes in inertia zones of the network
 Index f_number_inertia_zones{};

 /// the number of pollutants
 Index f_number_pollutants{};

 /// the starting index of each NetworkBlock
 std::vector< Index > v_start_network_intervals;

 /// the constant terms of each NetworkBlock
 std::vector< double > v_network_constant_terms;

 /* TODO commented away until HeatBlock are properly managed
 /// the set of HeatBlock
 std::vector< HeatBlock * > v_heat_blocks;
 */

 /// the number of pollutant zones of each pollutant
 std::vector< Index > v_number_pollutant_zones;

 /// the matrix of pollutant zones
 /** Indexed over the dimensions NumberPollutants and NumberNodes. */
 boost::multi_array< Index , 2 > v_pollutant_zones;

 /// vector of pointers to the NetworkBlock.
 /** This vector has size f_time_horizon. So the NetworkBlock at
  * position t in this vector refers to the network at the t-th time step.
  */
 std::vector< NetworkBlock * > v_network_blocks;

 /// the matrix of ActivePowerDemand
 /** Indexed over the dimensions NumberNodes and TimeHorizon. */
 boost::multi_array< double , 2 > v_active_power_demand;

 /// the vector of PrimaryZones
 std::vector< Index > v_primary_zones;

 /// the matrix of PrimaryDemand
 /** Indexed over the dimensions PrimaryZones and TimeHorizon. */
 boost::multi_array< double , 2 > v_primary_demand;

 /// the vector of SecondaryZones
 std::vector< Index > v_secondary_zones;

 /// the matrix of SecondaryDemand
 /** Indexed over the dimensions SecondaryZones and TimeHorizon. */
 boost::multi_array< double , 2 > v_secondary_demand;

 /// the vector InertiaZones
 std::vector< Index > v_inertia_zones;

 /// the matrix of InertiaDemand
 /** Indexed over the dimensions InertiaZones and TimeHorizon. */
 boost::multi_array< double , 2 > v_inertia_demand;

 /// the vector of PollutantBudget
 /** Indexed over the pair of each NumberPollutantZone and NumberPollutant*/
 std::vector< std::vector< double >> v_pollutant_budget;

 /// the PollutantRho matrix
 /** Indexed over TimeHorizon, NumberPollutants, and NumberElcGenerators. */
 boost::multi_array< double , 3 > v_pollutant_rho;

 /* TODO commented away until HeatBlock are properly managed
 /// the PollutantHeatRho matrix
 /// Indexed over TimeHorizon, NumberPollutants, and NumberHeatBlocks.
 boost::multi_array< double, 3 > v_pollutant_heat_rho;
 */

 /// v_generator_node[ g ] tells to which node generator g belongs
 std::vector< Index > v_generator_node;

 /* TODO commented away until HeatBlock are properly managed
 /// v_heat_node[ h ] tells to which node the heat block h belongs
 std::vector< Index > v_heat_node;
 */

 /* TODO commented away until HeatBlock are properly managed
 /// the HeatSet vector, indexed over NumberHeatGenerators
 std::vector< Index > v_heat_set;
 */

 /* TODO commented away until HeatBlock are properly managed
 /// vector of heat rho
 std::vector< double > v_power_heat_rho;
 */

/*-------------------------------- variables -------------------------------*/



/*------------------------------- constraints ------------------------------*/

 /// node injection constraints for each time and node
 boost::multi_array< FRowConstraint , 2 > v_node_injection_Const;

 /// primary demand constraints for each time and primary zone
 boost::multi_array< FRowConstraint , 2 > v_PrimaryDemand_Const;

 /// secondary demand constraints for each time and secondary zone
 boost::multi_array< FRowConstraint , 2 > v_SecondaryDemand_Const;

 /// inertia demand constraints for each time and inertia zone
 boost::multi_array< FRowConstraint , 2 > v_InertiaDemand_Const;

 /* TODO commented away until HeatBlock are properly managed
 /// heat constraints for each time and index unit
 boost::multi_array< FRowConstraint , 2 > v_power_Heat_Rho_Const;
 */

 /// pollutant demand constraints for each pollutant and pollutant zone
 std::vector< std::vector< FRowConstraint > > v_PollutantBudget_Const;

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

 boost::multi_array< Range , 2 > primary_var_index;
 ///< indices of active Variable in the primary demand constraints
 /**< The active Variables of each LinearFunction defining a primary demand
  * constraint are grouped by UnitBlocks. That is, all active Variables of a
  * given UnitBlock have consecutive indices in the LinearFunction that
  * defines each constraint. The element at position (i, z) in the multi-array
  * will store the range of indices of active Variable (that belong to the
  * i-th UnitBlock) in the constraint associated with zone "z". If no Variable
  * of the i-th UnitBlock is active in the constraint associated with zone
  * "z", then the element at position (i, z) is
  * ( Inf< Index >() , Inf< Index >() ).
  *
  * Notice that these indices do not depend on the time instant. This is
  * because we assume that, if a generator has primary spinning reserve for
  * some time instant, then it has primary spinning reserve for all time
  * instants. */

 boost::multi_array< Range , 2 > secondary_var_index;
 ///< indices of active Variable in the secondary demand constraints
 /**< The active Variables of each LinearFunction defining a secondary demand
  * constraint are grouped by UnitBlocks. That is, all active Variables of a
  * given UnitBlock have consecutive indices in the LinearFunction that
  * defines each constraint. The element at position (i, z) in the multi-array
  * will store the range of indices of active Variable (that belong to the
  * i-th UnitBlock) in the constraint associated with zone "z". If no Variable
  * of the i-th UnitBlock is active in the constraint associated with zone
  * "z", then the element at position (i, z) is
  * ( Inf< Index >() , Inf< Index >() ).
  *
  * Notice that these indices do not depend on the time instant. This is
  * because we assume that, if a generator has secondary spinning reserve for
  * some time instant, then it has secondary spinning reserve for all time
  * instants. */

 boost::multi_array< Index , 2 > inertia_var_index;
 ///< indices of active Variable in the inertia demand constraints
 /**< The active Variables of each LinearFunction defining a inertia demand
  * constraint are grouped by UnitBlocks. That is, all active Variables of a
  * given UnitBlock have consecutive indices in the LinearFunction that
  * defines each constraint. The element at position (i, z) in the multi-array
  * will store the smallest index of an active Variable (that belong to the
  * i-th UnitBlock) in the constraint associated with zone "z". If no Variable
  * of the i-th UnitBlock is active in the constraint associated with zone
  * "z", then the element at position (i, z) is Inf< Index >().
  *
  * Notice that these indices do not depend on the time instant. This is
  * because we assume that, if a generator has commitment variable, inertia
  * commitment, inertia power, or active power variable for some time instant,
  * then it has the same thing for all time instants. */

 SMSpp_insert_in_factory_h;

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
 /// generate the node injection constraints

 void generate_node_injection_constraints( void );

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
 /// generate the heat constraints

 void generate_heat_constraints( void );

/*--------------------------------------------------------------------------*/
 /// updates the node injection constraints
 /** This function updates the node injection constraints considering that the
  * scale factors of the given units may have been modified. The vector \p
  * modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered. */

 void update_node_injection_constraints
  ( const std::vector< Index > & modified_units );

/*--------------------------------------------------------------------------*/
 /// updates the primary demand constraints
 /** This function updates the primary demand constraints considering that the
  * scale factors of the given units may have been modified. The vector \p
  * modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered. */

 void update_primary_demand_constraints
  ( const std::vector< Index > & modified_units );

/*--------------------------------------------------------------------------*/
 /// updates the secondary demand constraints
 /** This function updates the secondary demand constraints considering that the
  * scale factors of the given units may have been modified. The vector \p
  * modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered. */

 void update_secondary_demand_constraints
  ( const std::vector< Index > & modified_units );

/*--------------------------------------------------------------------------*/
 /// updates the inertia demand constraints
 /** This function updates the inertia demand constraints considering that the
  * scale factors of the given units may have been modified. The vector \p
  * modified_units is assumed to be ordered.
  *
  * @param modified_units The indices of the UnitBlocks that may have been
  *        modified. This vector is assumed to be ordered. */

 void update_inertia_demand_constraints
  ( const std::vector< Index > & modified_units );

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
 /// returns the number of nodes

 Index get_number_nodes( void ) const {
  return( f_NetworkData ? f_NetworkData->get_number_nodes() : 1 );
 }

/*--------------------------------------------------------------------------*/
 /// returns the primary zone to which the given electrical generator belongs

 Index get_primary_zone( Index elc_generator ) const {
  if( f_number_primary_zones == 0 )
   return( 0 );

  // Node to which the given electrical generator belongs
  Index node = 0;
  if( get_number_nodes() > 1 )
   node = v_generator_node[ elc_generator ];

  return( v_primary_zones[ node ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the secondary zone to which the given electrical generator belongs

 Index get_secondary_zone( Index elc_generator ) const {
  if( f_number_secondary_zones == 0 )
   return( 0 );

  // Node to which the given electrical generator belongs
  Index node = 0;
  if( get_number_nodes() > 1 )
   node = v_generator_node[ elc_generator ];

  return( v_secondary_zones[ node ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the inertia zone to which the given electrical generator belongs

 Index get_inertia_zone( Index elc_generator ) const {
  if( f_number_inertia_zones == 0 )
   return( 0 );

  // Node to which the given electrical generator belongs
  Index node = 0;
  if( get_number_nodes() > 1 )
   node = v_generator_node[ elc_generator ];

  return( v_inertia_zones[ node ] );
 }

/*--------------------------------------------------------------------------*/

 static void static_initialization( void ) {
  register_method< UCBlock , MF_dbl_it , Subset && , bool >(
   "UCBlock::set_active_power_demand" , &UCBlock::set_active_power_demand );

  register_method< UCBlock , MF_dbl_it , Range >(
   "UCBlock::set_active_power_demand" , &UCBlock::set_active_power_demand );
 }

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
  eSetActD = 0    ///< set active power demand
 };

 /// constructor, takes the UCBlock and the type
 UCBlockMod( UCBlock * const fblock , const int type )
  : f_Block( fblock ) , f_type( type ) {}

 /// destructor, does nothing
 virtual ~UCBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 /// accessor to the type of modification
 int type( void ) { return( f_type ); }

 protected:

 /// prints the UCBlockMod
 void print( std::ostream & output ) const override {
  output << "UCBlockMod[" << this << "]: ";
  switch( f_type ) {
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
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* UCBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File UCBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
