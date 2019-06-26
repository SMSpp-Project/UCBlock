/*--------------------------------------------------------------------------*/
/*--------------------------- File UCBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class UCBlock, which derives from the Block, in order
 * to define a base class for the unit commitment problem that can be
 * considered as a father block of any other possible blocks(such as
 * UnitBlocks and NetworkBlocks) which are attached to it. It has very basic
 * information that can characterize almost any different different variation
 * of the Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 23 - 06 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
#ifndef __UCBlock
#define __UCBlock /* self-identification: #endif at the end of the file */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <boost/multi_array.hpp>
#include <vector>
#include "Block.h"
#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

// forward declaration
class HeatBlock;
class NetworkBlock;
class UnitBlock;

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
 * manages son Block of type UnitBlock and NetworkBlock. Also, UCBlock handles
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
 *   = balance between production of active power and injection in the
 *     transmission network, at each node and for each time instant;
 *
 *   = possibly, primary and secondary reserve constraints for each "zone"
 *     (appropriately defined subset of the nodes of the transmission
 *     network) and for each time instant;
 *
 *   = possibly, constraints about inertia  for each "zone" (appropriately
 *     defined subset of the nodes of the transmission network) and for
 *     each time instant;
 *
 *   = possibly, constraints maximum pollutants emission for different kinds
 *     of pollutant, each "zone" (appropriately defined subset of the nodes
 *     of the transmission network) and for each time instant;
 *
 *   = possibly, constraints linking the electricity production of some
 *     UnitBlock with the heat production of some unit in a HeatBlock,
 *     for the appropriate units and for each time instant.
 */

class UCBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * UCBlock defines a main public type:
 *
 * - Index, the type of indices;
 @{ */

 typedef std::size_t Index;

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of UCBlock, taking possibly a pointer of its father Block
 explicit UCBlock( Block *father = nullptr ) : Block( father ) {}

/*--------------------------------------------------------------------------*/

 ~UCBlock() override = default;   ///< destructor of UCBlock: it is virtual, and empty

/**@} ----------------------------------------------------------------------*/
/*---------------------------- sub-CLASS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup the NetworkData sub_Class in UCBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NetworkData -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// The NetworkData class
/** The NetworkData class is a nested sub_class into the UCBlock. It is
 * characterized by a set of nodes (NetworkNodes) and a set of lines (or arcs)
 * and all the data needed to construct the network (instead of the demand
 * of each node).
 * */
 class NetworkData {

/*--------------------------------------------------------------------------*/
/*----------------- PUBLIC PART OF THE NetworkData CLASS -------------------*/
/*--------------------------------------------------------------------------*/

 public:

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of NetworkData, does nothing
 NetworkData() = default;

 /// destructor of NetworkData: it is virtual, and empty
 virtual ~NetworkData() = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// deserialize a NetworkData out of a netCDF::NcGroup
/** Deserialize a NetworkData out of a netCDF::NcGroup, which should contain
 * the following:
 *
 * - The dimension "NumberNodes" containing the number of nodes in the
 *   problem; this dimension is optional, if it is not provided then it is
 *   taken to be == 1.
 *
 * If NumberNodes == 1 (equivalently, it is not provided), the network is a
 * "bus" formed of only one node, and therefore all the subsequent information
 * need not to be present since it is not loaded. If NumberNodes > 1, then all
 * the subsequent information is mandatory:
 *
 * - The dimension "NumberLines" containing the number of lines in the
 *   transmission network.
 *
 * - The variable "StartLine", of type int and indexed over the dimension
 *   "NumberNodes"; the i-th entry of the variable is the starting point of
 *   the line (a number in 0, ..., NumberNodes - 1). Note that lines are not
 *   oriented, but the flow of energy is; that is, a positive flow along
 *   line i means that energy is being taken away from StartLine[ i ] and
 *   delivered to EndLine[ i ] (see next), a negative flow means vice-versa.
 *
 * - The variable "EndLine", of type int and indexed over the dimension
 *   "NumberNodes"; the i-th entry of the variable is the ending point of the
 *   line (a number in 0, ..., NumberNodes - 1; lines are not oriented, but
 *   see above). StartLine[ i ] == EndLine[ i ] (a self-loop) is not allowed,
 *   but multiple lines between the same pair of nodes are.
 *
 * - The variable "MinPowerFlow", of type double and indexed over the
 *   dimension "NumberLines"; the i-th entry of the variable is assumed to
 *   contain the minimum power flow on line i (note that this is typically
 *   a negative number as lines are bi-directional, see above).
 *
 * - The variable "MaxPowerFlow", of type double and indexed over the
 *   dimension "NumberLines"; the i-th entry of the variable is assumed to
 *   contain the maximum power flow at line i (a non-negative number).
 *
 * - The variable "Susceptance", of type double and indexed over the dimension
 *   "NumberLines"; the i-th entry of this variable is assumed to contain the
 *   susceptance of line i.
 *
 * //TODO: Is Susceptance[ i ] assumed ro be > 0? If so let's say it.
 *    //I HAVENT FOUND IT INTO THE DOCUMENTATION
 *
 * //TODO: In UCBlock::NetworkData::deserialize(), NumberLines need not be
 *       read if NumberNodes == 1 (or not present). Also, we have to make
 *       the basic checks on data:
 *       - nodes starting and ending names are in 0 ... NumNodes - 1
 *       - self loops are not allowed
 *       - min capacity <= 0 <= max capacity
 *       - susceptance > 0 (if it is)
 */

 virtual void deserialize( netCDF::NcGroup & group );

/*--------------------------------------------------------------------------*/

 virtual void load( std::istream &input )  {
     throw( std::logic_error( "NetworkData::load() not implemented yet" ) );
 };

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE NetworkData -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NetworkData
 *  @{ */

/// serialize a NetworkData out of a netCDF::NcGroup
/** Serialize a NetworkData out of a netCDF::NcGroup to the specific format of a
 * NetworkData. See NetworkBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 */
 virtual void serialize( netCDF::NcGroup & group ) const;

/*--------------------------------------------------------------------------*/
/*---------------- PUBLIC FIELDS OF THE NetworkData CLASS ------------------*/
/*--------------------------------------------------------------------------*/

 /// returns the number of nodes of the problem
 Index get_number_nodes( void ) const { return f_number_nodes; }

 /// returns the number of nodes of the problem
 Index get_number_lines( void ) const { return f_number_lines; }

 // TODO: I like it better??

 const std::vector< int > & get_start_line( void ) const {
  return( v_start_line );
  }


 /// returns the start line where the given node belongs to
  inline Index get_start_line( Index node ) const {
    if( !v_start_line.empty() )
      return v_start_line[ node ];
        return 0;
 }

 /// returns the end line where the given node belongs to
  inline Index get_end_line( Index node ) const {
   if( !v_end_line.empty() )
     return v_end_line[ node ];
   return 0;
 }

 /// returns the minimum power flow for the given line l
  inline double get_min_power_flow( Index line ) const {
    if( !v_min_power_flow.empty() )
      return v_min_power_flow[ line ];
    return 0;
 }

 /// returns the maximum power flow for the given line l
 inline double get_max_power_flow( Index line ) const {
   if( !v_max_power_flow.empty() )
     return v_max_power_flow[ line ];
   return 0;
 }

 /// returns the Susceptance for the given line l
 inline double get_susceptance( Index line ) const {
   if( !v_susceptance.empty() )
     return v_susceptance[ line ];
   return 0;
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------- METHODS FOR MODIFYING THE NetworkData -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkData
 *  @{ */

 /// number of nodes of the network
 Index f_number_nodes;

/// number of lines of the network
 Index f_number_lines;

/// set starting lines
 std::vector< int > v_start_line;

/// set ending lines
 std::vector< int > v_end_line;

/// vector to store the susceptance of each line of the network
 std::vector< double > v_susceptance;

/// vector to store the minimum power flow at each line
 std::vector< double > v_min_power_flow;

/// vector to store the maximum power flow at each line
 std::vector< double > v_max_power_flow;


 };   // end( class( NetworkData) )

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 *  the UCBlock. Besides the mandatory "type" attribute of any :Block, the
 *  group should contain the following:
 *
 * - The dimension "TimeHorizon" containing the number of time steps in the
 *   problem.
 *
 * - The dimension "NumberUnits" containing the number of electricity
 *   generating units (UnitBlock) in the problem;
 *
 * - The dimension "NumberHeatBlocks" containing the number of heat blocks in
 *   the problem. The dimension is optional: if it is not provided then it is
 *   taken to be 0, which means that there is no heat block in the problem.
 *
 * - The groups "UnitBlock_0", "UnitBlock_1", ... , "UnitBlock_n" with
 *   n == NumberUnits - 1, containing each one UnitBlock corresponding
 *   to one electricity generating unit.
 *
 * - The groups "HeatBlock_0", "HeatBlock_1", ... , "HeatBlock_n" with
 *   n == NumberHeatBlocks - 1, containing each one a HeatBlock. When
 *   NumberHeatBlocks == 0, these groups need not be there since they are
 *   not read.
 *
 * - Possibly, the dimensions and variables necessary to a NetworkData object,
 *   that describe the transmission network; see NetworkData::deserialize()
 *   for details. All that is optional, if it is not provided (basically,
 *   "NumberNodes" is not provided or it is == 1) then the transmission
 *   network is taken to have only one node (a bus).
 *
 * - The groups "NetworkBlock_0", "NetworkBlock_1", ... , "NetworkBlock_t"
 *   with t = TimeHorizon - 1, containing each the constraints on the
 *   transmission network at time t.
 *
 * - The variable "UnitNode", of type int and indexed over the dimension
 *   "NumberUnits"; the entry UnitNode[ i ] tells to which node of the
 *   transmission network unit i belongs. If NumberNodes == 1 (say, it is
 *   not provided at all), then this variable need not be defined, since it
 *   is not loaded.
 *
 * - The variable "HeatNode", of type int and indexed over the dimension
 *   "NumberHeatBlocks"; the entry HeatNode[ h ] tells to which node of the
 *   transmission network all the heat-generating units that are also
 *   electricity-generating ones in HeatBlock h belong. If NumberHeatBlocks
 *   == 0 (say, it is not provided at all), then this variable need not be
 *   defined, since it is not loaded. This information is actually only used
 *   to determine in which Pollutant Zone a HeatBlock is located, in order to
 *   add the corresponding heat-generating units (that are not also
 *   electricity-generating ones) to the pollution constraints. This means
 *   that also if NumberPollutants == 0 (say, it is not provided at all)
 *   this variable is useless and therefore need not be defined, since it is
 *   not loaded. Finally, notice that for units into a HeatBlock that also
 *   are electrical units, this variable provides another time an information
 *   that is already known, i.e., to which node they belong to. Of course *the
 *   two information must agree*, otherwise the input file is ill-defined and
 *   exception is thrown.
 *
 * - The variable "HeatSet", of type int and indexed both over the dimensions
 *   "NumberUnits" and "NumberHeatBlocks". If HeatSet[ i , h ] = k, with
 *   k < number of heat units in HeatBlock h, then electrical unit i is
 *   represented into HeatBlock h as the heat unit k. If, instead,
 *   HeatSet[ i , h ] = k, with k >= number of heat units in HeatBlock h, then
 *   none of the heat units in HeatBlock h represents the electrical unit i.
 *   If NumberHeatBlocks == 0 (say, there is no HeatBlock) then this variable
 *   need not be defined, since it is not loaded.
 *
 * - The variable "PowerHeatRho", of type double and indexed over the
 *   dimension "NumberUnits": entry PowerHeatRho[ i ] is assumed to contain
 *   the electrical-power-to-heat ratio for unit i.
 *
 * - The dimension "NumberPrimaryZones" tells how many "primary spinning
 *   reserve zones" are there in the problem. The dimension is optional, if it
 *   is not provided then it is taken to be 0, which means that no primary
 *   reserve constraints are present in the problem.
 *
 * - The variable "PrimaryZones", of type int and indexed over the dimension
 *   "NumberNodes". The entry PrimaryZones[ i ] tells to which primary zone
 *   the node i belongs: if PrimaryZones[ i ] >= NumberPrimaryZones, this
 *   means that node i does not belong to any primary zone, and hence the
 *   corresponding units are not involved into the primary reserve
 *   constraints. If NumberPrimaryZones == 0  (say, it is not provided at all)
 *   then this variable need not be defined, since it is not loaded. If
 *   NumberPrimaryZones == 1 and this variable is not defined, then there is
 *   only one primary zone and all the nodes belong to it.
 *
 * - The variable "PrimaryDemand", of type double and indexed both over the
 *   dimensions "PrimaryZones" and "TimeHorizon": entry PrimaryDemand[ i , t ]
 *   is assumed to contain the primary reserves requirement which are
 *   specified on the primary reserve zone i in the time t. If
 *   NumberPrimaryZones == 0 (say, it is not provided at all), then this
 *   variable need not be defined, since it is not loaded.
 *
 * - The dimension "NumberSecondaryZones" tells how many "secondary spinning
 *   reserve zones" are there in the problem. The dimension is optional, if it
 *   is not provided then it is taken to be 0, which means that no secondary
 *   reserve constraints are present in the problem.
 *
 * - The variable "SecondaryZones", of type int and indexed over the dimension
 *   "NumberNodes"; the entry SecondaryZones[ i ] tells to which secondary
 *   zone the node i belongs. If SecondaryZones[ i ] >= NumberSecondaryZones,
 *   this means that node i does not belong to any secondary zone, and hence
 *   the corresponding units are not involved into the secondary reserve
 *   constraints. If NumberSecondaryZones == 0 (say, it is not provided at
 *   all) then this variable need not be defined, since it is not loaded. If
 *   NumberSecondaryZones == 1 and this variable is not defined, then there is
 *   only one secondary zone and all the nodes belong to it.
 *
 * - The variable "SecondaryDemand", of type double and indexed both over the
 *   dimensions "SecondaryZones" and "TimeHorizon": entry
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
 * - The variable "InertiaZones", of type int and indexed over the dimension
 *   "NumberNodes"; the entry InertiaZones[ n ] tells to which inertia zone
 *   the node n belongs. If InertiaZones[ n ] >= NumberInertiaZones, this
 *   means that node n does not belong to any inertia zone, and hence the
 *   corresponding units are not involved into the inertia reserve
 *   constraints. If NumberInertiaZones == 0 (say, it is not provided at all)
 *   then this variable need not be defined, since it is not loaded. If
 *   NumberInertiaZones == 1 and this variable is not defined, then there is
 *   only one inertia zone and all the nodes belong to it.
 *
 * - The variable "InertiaDemand", of type double and indexed both over the
 *   dimensions "InertiaZones" and "TimeHorizon": entry InertiaDemand[ i , t ]
 *   is assumed to contain the inertia reserves requirement which are
 *   specified on the inertia constraints zone i in the time t. If
 *   NumberInertiaZones == 0 (say, it is not provided at all), then this
 *   variable need not be defined, since it is not loaded.
 *
 * - The dimension "NumberPollutants" containing the number of pollutants in
 *   the problem. The dimension is optional, if it is not provided then it is
 *   taken to be 0, which means that no pollutants  constraints are present in
 *   the problem.
 *
 * - The variable "NumberPollutantZones" of type int indexed over the
 *   dimension "NumberPollutants": the entry NumberPollutantZones[ p ] is
 *   assumed to contain the number of pollutant zones associated with
 *   pollutant p. If NumberPollutants == 0 (say, it is not provided) then
 *   this variable need not be defined, since it is not loaded.
 *
 * - The variable "PollutantZones", of type int and indexed over the
 *   dimensions "NumberPollutants" and "NumberNodes": the entry
 *   PollutantZones[ p , n ] tells to which pollutant zone associated
 *   with pollutant p the node n belongs. If PollutantZones[ p , n ] >=
 *   NumberPollutantZones[ p ], this means that node n does not belong to
 *   any pollutant zone, and hence the corresponding units are not involved
 *   into the pollutant demand constraints associated with pollutant p. If
 *   NumberPollutants == 0 (say, it is not provided) then this variable
 *   need not be defined, since it is not loaded.
 *
 * - The variable "PollutantBudget", of type double and indexed over the
 *   dimension "NumberPollutants": the entry PollutantBudget[ i ] is assumed
 *   to contain the total limit (across all the time horizon) of pollutant i.
 *   If NumberPollutants == 0 (say, it is not provided) then this variable
 *   need not be defined, since it is not loaded.
 *
 * - The variable "PollutantRho", of type double and indexed over three
 *   dimensions. The first dimension can have either size 1 or TimeHorizon.
 *   The second and third dimensions have sizes "NumberPollutants" and
 *   "NumberUnits", respectively. The entry PollutantRho[ t , p , i ] is
 *   assumed to contain the conversion factor of pollutant p due to the
 *   generation of unit i for "each" time instant t (when the first dimension
 *   has full size TimeHorizon) or for "all" time instants (when the first
 *   dimension has size 1). If NumberPollutants == 0 (say, it is not
 *   provided) then this variable need not be defined, since it's not loaded.
 *
 * - The variable "PollutantHeatRho", of type double and indexed over three
 *   dimensions. The first dimension can have size 1 or TimeHorizon. The
 *   second and third dimensions have sizes "NumberPollutants" and
 *   "NumberHeatBlocks", respectively. The entry PollutantRho[ t , p , h ] is
 *   assumed to contain the conversion factor of pollutant p due to the
 *   generation of every heat-only unit in HeatBlock h for "each" time instant
 *   t (when the first dimension has full size TimeHorizon) or for "all" time
 *   instants (when the first dimension has size 1). If NumberPollutants == 0
 *   (say, it is not provided) then this variable does not need be defined,
 *   since it is not loaded. If NumberHeatBlocks == 0 (say, there is no
 *   heat-only unit) then this variable need not be defined, since it is not
 *   loaded. */
 void deserialize( netCDF::NcGroup & group ) override;



/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the UCBlock
 /** Method that generates the abstract constraint of the UCBlock. These are:
  *
  * - if f_number_nodes > 0, a boost::multi_array<FRowConstraint, 2>; with two
  *   dimensions which are f_time_horizon, and f_number_nodes entries, the
  *   entry t = 0, ..., f_time_horizon - 1 and the entry
  *   n = 1, ..., f_number_nodes being the node injection constraints at time
  *   t and node n;
  *
  * - if f_number_primary_zones > 0, a boost::multi_array<FRowConstraint, 2>;
  *   with two dimensions which are f_time_horizon, and f_number_primary_zones
  *   entries, the entry t = 0, ..., f_time_horizon - 1 and the entry
  *   z = 1, ..., f_number_primary_zones being the primary demand constraints
  *   at time t and primary zone z;
  *
  * - if f_number_secondary_zones > 0, a boost::multi_array<FRowConstraint,2>;
  *   with two dimensions which are f_time_horizon, and
  *   f_number_secondary_zones entries, the entry
  *   t = 0, ..., f_time_horizon - 1 and the entry
  *   z = 1, ..., f_number_secondary_zones being the secondary demand
  *   constraints at time t and secondary zone z;
  *
  * - if f_number_inertia_zones > 0, a boost::multi_array<FRowConstraint, 2>;
  *   with two dimensions which are f_time_horizon, and f_number_inertia_zones
  *   entries, the entry t = 0, ..., f_time_horizon - 1 and the entry
  *   z = 1, ..., f_number_inertia_zones being the inertia demand constraints
  *   at time t and inertia zone z;
  *
  * - if f_number_pollutants > 0, a std::vector<std::vector<FRowConstraint>>;
  *   with two dimensions which are f_number_pollutants, and
  *   v_number_pollutant_zones entries, the entry
  *   p = 1, ..., f_number_pollutants and the entry
  *   z = 1, ..., v_number_pollutant_zones being the pollutant budget
  *   constraints at pollutant p and pollutant zones z;
  *
  * - if f_number_heat_blocks > 0, a boost::multi_array<FRowConstraint, 2>
  *   with two dimensions which are f_time_horizon and the number of
  *   UnitBlocks that produce electricity and belong to some
  *   HeatBlock; the constraint at position ( t, i ) being the heat
  *   constraints at time t and unit M[ i ], where M maps the
  *   constraint into an electricity-producing unit that belongs to
  *   some HeatBlock.
  *
 *  Consider a network defined by a set of nodes \f$ \mathcal{N} \f$ and a set
 *  of arcs connecting the nodes \f$ \mathcal{L} \f$. There are moreover given
 *  three partitions of the set of nodes which may or may not be identical:
 *
 *  (i). \f$ \mathcal{B}^{pr}(\mathcal{N}) \f$ partitions \f$ \mathcal{N} \f$
 *  in several zones (sets of nodes) each one being associated with one
 *  specific primary spinning reserve requirement;
 *
 *  (ii). \f$ \mathcal{B}^{sc}(\mathcal{N}) \f$ partitions \f$ \mathcal{N}
 *  \f$ in several zones each one being associated with one specific secondary
 *  spinning reserve requirement;
 *
 *  (iii). \f$ \mathcal{B}^{in}(\mathcal{N}) \f$ partitions \f$ \mathcal{N}
 *  \f$ in several zones each one being associated with one specific inertia
 *  requirement;
 *
 *  Optionally we are given a partition of the nodes \f$ B^{p}(\mathcal{N})
 *  \f$ corresponding to zones which are associated with an emissions
 *  constraint on the specific pollutant \f$ p \f$ in the set of pollutants
 *  \f$ \mathcal{P} \f$.
 *
 *  The electrical system contains a set of “units” (e.g., power plants, load
 *  flexibilities or storage devices) indexed by \f$ i \in \mathcal{I} \f$.
 *  The set \f$ \mathcal{I}_n \f$ will indicate units connected to node \f$ n
 *  \in \mathcal{N} \f$.
 *
 *  The decision variables are introduced as:
 *
 * - \f$ p^{ac}_{t,i} \f$ : the active power variable for each time period
 *   \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I} \f$;
 *
 * - \f$ S_{t,n} \f$ : the node injection variable for each time period
 *   \f$ t \in \mathcal{T} \f$ and each node \f$ n \in \mathcal{N} \f$;
 *
 * - \f$ p^{pr}_{t,i} \f$ : the primary spinning reserves variable for each
 *   time period \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I}
 *   \f$;
 *
 * - \f$ p^{sc}_{t,i} \f$ : the secondary spinning reserves variable for each
 *   time period \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I}
 *   \f$;
 *
 * - \f$ u_{t,i}  \in \{ 0 , 1 \} \f$ : the commitment state at time period
 *   \f$ t \in \mathcal{T} \f$ for each unit \f$ i \f$;
 *
 * - \f$ p^{he}_{t,i} \f$ : the heat variable for each time period
 *   \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I} \f$;
 *
 *  The global constraints of unit commitment problem, on the time horizon
 *  \f$ \mathcal{T} \f$ write as follow:
 *
 * - Node injection Constraints:
 *   In the unit commitment problem, \f$ P^{au}_{t , i} \f$ denotes the fixed
 *   consumption of the power plant when it is off, and \f$ S_{t,n} \f$ is the
 *   node injection variable for each time period \f$ t \in \mathcal{T} \f$
 *   and each node \f$ n \in \mathcal{N} \f$. The node injection constraints
 *   will be satisfied as follow:
 *
 * \f[
 *  \sum_{ i \in \mathcal{I}_n } (p^{ac}_{t,i} + P^{au}_{t , i}(1 - u_{t,i}))
 *     = S_{t,n} \quad t \in \mathcal{T} \quad n \in \mathcal{N} \quad     (1)
 * \f]
 *
 * - Primary Demand Constraints:
 *   In the unit commitment problem, the primary demand
 *   \f$ D^{pr}_{\mathcal{B} , t} \f$ which are specified on the primary
 *   reserve zones \f$ \mathcal{B} \in \mathcal{B}^{pr}(\mathcal{N}) \f$ will
 *   be satisfied as follow:
 *
 * \f[
 *  \sum_{n \in \mathcal{B}}\sum_{ i \in I_n } p^{pr}_{t,i} \geq
 *   D^{pr}_{\mathcal{B} , t} \quad t \in \mathcal{T}
 *      \quad \mathcal{B} \in \mathcal{B}^{pr}(\mathcal{N}) \quad          (2)
 * \f]
 *
 * - Secondary Demand Constraints:
 *   In the unit commitment problem, the secondary demand
 *   \f$ D^{sc}_{\mathcal{B} , t} \f$ which are specified on the secondary
 *   reserve zones \f$ \mathcal{B} \in \mathcal{B}^{sc}(\mathcal{B}) \f$ will
 *   be satisfied as follow:
 *
 * \f[
 *  \sum_{n \in \mathcal{B}}\sum_{ i \in \mathcal{I}_n } p^{sc}_{t,i} \geq
 *       D^{sc}_{\mathcal{B} , t} \quad t \in \mathcal{T}
 *       \quad \mathcal{B} \in \mathcal{B}^{sc}(\mathcal{N}) \quad         (3)
 * \f]
 *
 * - Inertia Demand Constraints:
 *   In the unit commitment problem, the inertia demand
 *   \f$ D^{in}_{\mathcal{B} , t}\f$ which are specified on the inertia zones
 *   \f$ \mathcal{B} \in \mathcal{B}^{in}(\mathcal{N}) \f$ with defined
 *   parameters \f$ \alpha_{t , i} \f$ and \f$ \beta_{t , i} \f$ will be
 *   satisfied as follow:
 *
 * \f[
 *  \sum_{n \in \mathcal{B}}\sum_{ i \in \mathcal{I}_n } (\alpha_{t , i}
 *  u_{t,i} + \beta_{t , i} p^{ac}_{t,i}) \geq D^{in}_{\mathcal{B} , t}
 *        \quad t \in \mathcal{T}
 *        \quad \mathcal{B} \in \mathcal{B}^{in}(\mathcal{N}) \quad        (4)
 * \f]
 *
 * - Pollutant Budget Constraints:
 *   In the unit commitment problem, the pollutant budget \f$ \mathcal{O}_p
 *   \f$ which is specified on each pollutant \f$ p \in \mathcal{P} \f$ in
 *   each pollutant zone \f$ \mathcal{B} \in \mathcal{B}^{p}(\mathcal{N}) \f$
 *   and each pollutant heat zone \f$ \mathcal{B'} \in \mathcal{B}^{p}
 *   (\mathcal{H}) \f$ with two parameters \f$ \rho_{t , p , i} \f$ and \f$
 *   \rho'_{t , p , i} \f$ where considered as pollutant ratio and
 *   pollutant heat ratio respectively; is defined  as follow:
 *
 * \f[
 *
 *  \sum_{n \in \mathcal{B}}\sum_{ t \in \mathcal{T} }( \sum_{ i \in
 *  \mathcal{I}_n } \rho_{t , p , i} p^{ac}_{t,i} + \sum_{h \in \mathcal{H}_n}
 *  \sum_{ j \in \mathcal{I}^{ho}(h)} \rho'_{t , p , h} p^{h,he}_{t,j} )
 *  \leq \mathcal{O}_p  \quad \mathcal{B} \in \mathcal{B}^{p}(\mathcal{N})
 *  \quad p \in \mathcal{P} \quad                                          (5)
 * \f]
 *
 *   where \f$ \mathcal{H} \f$ is the set of Heat Blocks.
 *
 * - Heat Constraints:
 *   In the unit commitment problem, the Heat Constraints link the UCBlock
 *   variables with the HeatBlock, where for each heat block
 *   \f$ h \in \mathcal{H} \f$ and each electrical-power-to-heat ratio \f$
 *   \varrho_{i} \f$ of each heat producing unit \f$ i \f$; the Heat
 *   Constraints are defined as below:
 *
 * \f[
 *  \sum_{h \in \mathcal{H} , j \in \mathcal{I}^{ec}(h): e^h(j)=i}
 *   p^{h , he}_{t , j}  \leq \varrho_i p^{ac}_{t,i} \quad i \in \mathcal{I}
 *                                 \quad t \in \mathcal{T} \quad           (6)
 * \f]
 *   where \f$ j \in \mathcal{I}^{ec}(h) \f$ is an electricity producing unit
 *   in heat block \f$ h \in \mathcal{H} \f$. For \f$ j \in
 *   \mathcal{I}^{ec}(h) \f$, there is the need to know which electrical unit
 *   \f$ j \f$ is representing. Thus, we need a mapping
 *   \f$ e^h : \mathcal{I}^{ec}(h) \to \mathcal{I} \f$, where \f$ \mathcal{I}
 *   \f$ is the set of electricity producing units (standard units in UC
 *   parlance).
  */

 void generate_abstract_constraints( Configuration *stcc ) override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
 *  @{ */

 /// returns the time horizon of the problem
 Index get_time_horizon() const { return f_time_horizon; }

 /// returns the NetworkData object
 // TODO: a bit more comments regarding the fact that it may return nullptr
 NetworkData * get_NetworkData() const { return f_NetworkData; }

 /// returns the vector of (pointers to) NetworkBlocks
 const std::vector<NetworkBlock *> & get_network_blocks() const {
   return v_network_blocks;
 }

 // TODO: again consider returning const references to the vector
 //       but I'm open to discussing this
 /// returns the primary demand of the given zone at the given time
 inline double get_primary_demand( Index zone, Index time ) const {
   return v_primary_demand[ zone * f_time_horizon + time ];
 }

 /// returns the secondary demand of the given zone at the given time
 inline double get_secondary_demand( Index zone, Index time ) const {
   return v_secondary_demand[ zone * f_time_horizon + time ];
 }

 /// returns the inertia demand of the given zone at the given time
 inline double get_inertia_demand( Index zone, Index time ) const {
   return v_inertia_demand[ zone * f_time_horizon + time ];
 }

 /** returns the conversion factor of the given pollutant due to the
  * generation of the given unit at the given time */
 inline double get_pollutant_rho( Index time, Index pollutant, Index unit )
   const {
   auto index = time * f_number_pollutants * f_number_units +
     pollutant * f_number_units + unit;
   return v_pollutant_rho[ index ];
 }

 /** returns the conversion factor of the given pollutant due to the
  *  generation of every heat-only unit in the given heat block at the
  *  given time. */
 inline double get_pollutant_heat_rho( Index time, Index pollutant,
                                       Index heat_block )
   const {
   auto index = time * f_number_pollutants * f_number_heat_blocks +
     pollutant * f_number_heat_blocks + heat_block;
   return v_pollutant_heat_rho[ index ];
 }

 /** returns the pollutant zone associated with the given pollutant
  * the given node belongs to */
 inline Index get_pollutant_zone( Index pollutant, Index node ) const {
   return v_pollutant_zones[ pollutant * f_NetworkData->get_number_nodes()
                             + node ];
 }

 /** returns the heat unit that represents the given unit in the given
  * HeatBlock */
 inline Index get_heat_unit( Index unit, Index heat_block ) const {
   return v_heat_set[ unit * f_number_heat_blocks + heat_block ];
 }

 /// returns the i-th UnitBlock
 UnitBlock * get_unit_block( Index i ) const;

 /// returns the t-th NetworkBlock
 NetworkBlock * get_network_block( Index t ) const;

 /// returns the i-th HeatBlock
 HeatBlock * get_heat_block( Index i ) const;

 /// returns the node where the given unit belongs to
 inline Index get_unit_node( Index unit ) const {
   if( f_NetworkData->get_number_nodes() > 1 )
     return v_unit_node[ unit ];
   return 0;
 }

 /// returns the primary zone where the given node belongs to
 inline Index get_primary_zone( Index node ) const {
   if( !v_primary_zones.empty() )
     return v_primary_zones[ node ];
   return 0;
 }

 /// returns the secondary zone where the given node belongs to
 inline Index get_secondary_zone( Index node ) const {
   if( !v_secondary_zones.empty() )
     return v_secondary_zones[ node ];
   return 0;
 }


 /// returns the inertia zone where the given node belongs to
 inline Index get_inertia_zone( Index node ) const {
   if( !v_inertia_zones.empty() )
     return v_inertia_zones[ node ];
   return 0;
 }

 /// returns the electrical-power-to-heat ratio of the given unit
 inline double get_power_heat_rho( Index unit ) const {
   return v_power_heat_rho[ unit ];
 }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UCBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the UCBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * UCBlock. See UCBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UCBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UCBlock
    @{ */

 virtual void load( std::istream &input ) override {
   throw( std::logic_error( "UCBlock::load() not implemented yet" ) );
 };

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
protected:

 /// The time horizon of the problem
 Index f_time_horizon{};

 /// The number of units of the problem
 Index f_number_units{};

 /// the NetworkData object
 NetworkData * f_NetworkData{};

 /// The number of heat block
 Index f_number_heat_blocks{};

 /// The number of nodes in primary zones of the network
 Index f_number_primary_zones{};

 /// The number of nodes in secondary zones of the network
 Index f_number_secondary_zones{};

 /// The number of nodes in inertia zones of the network
 Index f_number_inertia_zones{};

 /// The number of pollutants
 Index f_number_pollutants{};

 /// The set of UnitBlocks
 std::vector<UnitBlock *> v_unit_blocks;

 /// The set of HeatBlocks
 std::vector<HeatBlock *> v_heat_blocks;

 /// The number of pollutant zones of each pollutant
 std::vector<Index> v_number_pollutant_zones;

 /** The matrix PollutantZones indexed over the dimensions
  *  NumberPollutants and NumberNodes */
 std::vector<Index> v_pollutant_zones;

 /** Vector of pointers to the NetworkBlocks. This vector either is
  * empty or has size f_time_horizon. If it is empty, it means there
  * is no network. If it has positive size, then the NetworkBlock at
  * position i in this vector refers to the network at the i-th time
  * step. */
 std::vector<NetworkBlock *> v_network_blocks;

 /// the vector of PrimaryZones
 std::vector<Index> v_primary_zones;

 /** the matrix of PrimaryDemand indexed over the dimensions
  * PrimaryZones and TimeHorizon */
 std::vector<double> v_primary_demand;

 /// the vector of SecondaryZones
 std::vector<Index> v_secondary_zones;

 /** the matrix of SecondaryDemand indexed over the dimensions
  * SecondaryZones and TimeHorizon */
 std::vector<double> v_secondary_demand;

 /// the vector InertiaZones
 std::vector<Index> v_inertia_zones;

 /** the matrix of InertiaDemand indexed over the dimensions
  * InertiaZones and TimeHorizon */
 std::vector< double > v_inertia_demand;

 /// the vector of PollutantDemand
 std::vector< double > v_pollutant_budget;

 /** the PollutantRho matrix index over the dimensions
  * TimeHorizon, NumberPollutants, and NumberUnits */
 std::vector< double > v_pollutant_rho;

 /** the PollutantHeatRho matrix index over the dimensions
  * TimeHorizon, NumberPollutants, and NumberHeatBlocks */
 std::vector< double > v_pollutant_heat_rho;

 /// The entry v_unit_node[ i ] tells to which node unit i belongs
 std::vector<Index> v_unit_node;

 /// The entry v_heat_node[ h ] tells to which node the HeatBlock h belongs
 std::vector<Index> v_heat_node;

 /** the HeatSet matrix index over the dimensions
  * NumberUnits, and NumberHeatBlocks */
 std::vector< Index > v_heat_set;

 /// Vector of heat rho
 std::vector< double > v_power_heat_rho;

 /// Node injection constraints for each time and node
 boost::multi_array<FRowConstraint, 2> v_node_injection_constraints;

 /// Primary demand constraints for each time and primary zone
 boost::multi_array<FRowConstraint, 2> v_PrimaryDemand_Const;

 /// Secondary demand constraints for each time and secondary zone
 boost::multi_array<FRowConstraint, 2>  v_SecondaryDemand_Const;

 /// Inertia demand constraints for each time and inertia zone
 boost::multi_array<FRowConstraint, 2>  v_InertiaDemand_Const;

 /// heat constraints for each time and index unit
 boost::multi_array<FRowConstraint, 2>  v_power_Heat_Rho_Const;

 /// Pollutant demand constraints for each pollutant and pollutant zone
 std::vector<std::vector<FRowConstraint>> v_PollutantBudget_Const;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Deserialize the sub-blocks of UCBlock
 void deserialize_sub_blocks( const netCDF::NcGroup & group );

 /// Deserialize the sub-blocks of UCBlock that have the given prefix name
 void deserialize_sub_blocks( const netCDF::NcGroup & group,
                              const std::string& sub_group_name_prefix,
                              int num_sub_blocks );


  };   // end( class( UCBlock ) )

} /* namespace SMSpp_di_unipi_it */

#endif /* UCBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File UCBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
