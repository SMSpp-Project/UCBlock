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
 * \version 0.11
 *
 * \date 01 - 07 - 2019
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
 #define __UCBlock  /* self-identification: #endif at the end of the file */

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

class HeatBlock;     // forward declaration of HeatBlock
 
class NetworkBlock;  // forward declaration of NetworkBlock

class UnitBlock;     // forward declaration of UnitBlock

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
 * UCBlock defines two main public types:
 *
 * - Index, the type of indices;
 *
 * - NetworkData, a small auxiliary class to bunch together the basic
 *   data (topology and electrical characteristics) of the transmission
 *   network.
 *  @{ */

 typedef std::size_t Index;

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS UCBlock::NetworkData -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// auxiliary class to hold basic data about the transmission network
/** The NetworkData class is a nested sub_class into the UCBlock, which
 * only serves to have a quick way to load all the basic data (topology and
 * electrical characteristics) that describes the transmission network.
 * The rationale is that while often the network does not change during the
 * (short) time horizon of UC, it makes sense to allow for this to happen.
 * This means that individual NetworkBlock objects may in principle have
 * different NetworkData, but most often they can share the same. By bunching
 * all the information together we make it easy for this sharing to happen. */

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
 *   Note that node names here go from 0 to NNodes.getSize() - 1;
 *
 * - The variable "EndLine", of type int and indexed over the dimension
 *   "NumberNodes"; the i-th entry of the variable is the ending point of the
 *   line (a number in 0, ..., NumberNodes - 1; lines are not oriented, but
 *   see above). StartLine[ i ] == EndLine[ i ] (a self-loop) is not allowed,
 *   but multiple lines between the same pair of nodes are. Note that node
 *   names here go from 0 to NNodes.getSize() - 1;
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
 *   susceptance of line i. Note that this is strictly a positive value.
 *
 * //TODO: In UCBlock::NetworkData::deserialize(), NumberLines need not be
 *       read if NumberNodes == 1 (or not present). Also, we have to make
 *       the basic checks on data:
 *       - self loops are not allowed
 *       - min capacity <= 0 <= max capacity
 *       - susceptance > 0 (if it is)
 */
 virtual void deserialize( netCDF::NcGroup & group );

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE NetworkData ------------*/
/*--------------------------------------------------------------------------*/
 /// returns the number of nodes of the network
 /** Method for returning the number of nodes of the network.
  * */
 Index get_number_nodes( void ) const { return f_number_nodes; }

/*--------------------------------------------------------------------------*/
 /// returns the number of lines of the network
 /** Method for returning the number of lines of the network.
  * */
 Index get_number_lines( void ) const { return f_number_lines; }

/*--------------------------------------------------------------------------*/
 /// returns the start node of the given line
 /** Method for returning the vector of starting point of each line. There are
  *  two possible cases:
  *
  *  - if f_number_nodes == 1, this vector have empty size which means there is
  *    no line at network (bus network).
  *
  *  - if f_number_nodes > 1, this vector have size of f_number_nodes and each
  *    element of the vectors gives starting point of each line in the network.
  * */
 Index get_start_line( Index node ) const {
  return( v_start_line.empty() ? 0 : v_start_line[ node ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the vector of start nodes
 /** Method for returning the vector of starting point of each line.*/
 const std::vector< int > & get_start_line( void ) const {
  return( v_start_line );
  }

/*--------------------------------------------------------------------------*/
 /// returns the end node of the given line
 /** Method for returning the vector of ending point of each line. There are
  *  two possible cases:
  *
  *  - if f_number_nodes == 1, this vector have empty size which means there is
  *    no line at network (bus network).
  *
  *  - if f_number_nodes > 1, this vector have size of f_number_nodes and each
  *    element of the vectors gives ending point of each line in the network.
  * */
 Index get_end_line( Index node ) const {
  return( v_end_line.empty() ? 0 : v_end_line[ node ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns vector of end nodes
 /** Method for returning the vector of ending point of each line.*/
 const std::vector< int > & get_end_line( void ) const {
  return( v_end_line );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power flow for the given line l
 /** Method for returning the vector of minimum power flow of each line. There
  *  are two possible cases:
  *
  *  - if f_number_lines == 0, this vector have empty size which means there is
  *    no line at network (bus network).
  *
  *  - if f_number_lines >= 1, this vector have size of f_number_lines and each
  *    element of the vectors gives minimum power flow of each line in the network.
  * */
 double get_min_power_flow( Index line ) const {
  return( v_min_power_flow.empty() ? 0 : v_min_power_flow[ line ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns vector of the minimum power flow
 /** Method for returning the vector of minimum power flow of each line.*/
 const std::vector< double > & get_min_power_flow( void ) const {
  return( v_min_power_flow );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power flow for the given line l
 /** Method for returning the vector of maximum power flow of each line. There
  *  are two possible cases:
  *
  *  - if f_number_lines == 0, this vector have empty size which means there is
  *    no line at network (bus network).
  *
  *  - if f_number_lines >= 1, this vector have size of f_number_lines and each
  *    element of the vectors gives maximum power flow of each line in the network.
  * */
 double get_max_power_flow( Index line ) const {
  return( v_max_power_flow.empty() ? 0 : v_max_power_flow[ line ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns vector of the maximum power flow
 /** Method for returning the vector of maximum power flow of each line.*/
 const std::vector< double > & get_max_power_flow( void ) const {
  return( v_max_power_flow );
  }

/*--------------------------------------------------------------------------*/
 /// returns the Susceptance for the given line l
 /** Method for returning the vector of Susceptance of each line. There
  *  are two possible cases:
  *
  *  - if f_number_lines == 0, this vector have empty size which means there is
  *    no line at network (bus network).
  *
  *  - if f_number_lines >= 1, this vector have size of f_number_lines and each
  *    element of the vectors gives the Susceptance value for each line in the
  *    network.
  * */
 double get_susceptance( Index line ) const {
  return( v_susceptance.empty() ? 0 : v_susceptance[ line ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns vector of the Susceptances
 /** Method for returning the vector of Susceptances for each line.*/
 const std::vector< double > & get_susceptance( void ) const {
  return( v_susceptance );
  }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE NetworkData -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NetworkData
 *  @{ */

/// serialize a NetworkData out of a netCDF::NcGroup
/** Serialize a NetworkData out of a netCDF::NcGroup to the specific format of
 * a NetworkData. See NetworkBlock::deserialize( netCDF::NcGroup ) for details
 * of the format of the created netCDF group. */

 virtual void serialize( netCDF::NcGroup & group ) const;

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------- PROTECTED FIELDS OF THE NetworkData --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkData
 *  @{ */

 Index f_number_nodes;    ///< number of nodes of the network

 Index f_number_lines;    ///< number of lines of the network

 std::vector< int > v_start_line;  ///< vector of starting nodes

 std::vector< int > v_end_line;    ///< vector of ending nodes

 /// vector to store the susceptance of each line of the network
 std::vector< double > v_susceptance;

 /// vector to store the minimum power flow at each line
 std::vector< double > v_min_power_flow;

 /// vector to store the maximum power flow at each line
 std::vector< double > v_max_power_flow;

/*--------------------------------------------------------------------------*/

 };   // end( class( NetworkData ) )

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of UCBlock, taking possibly a pointer of its father Block

 explicit UCBlock( Block *father = nullptr ) : Block( father ) {}

/*--------------------------------------------------------------------------*/
  /// destructor of UCBlock: it is virtual, and empty

 ~UCBlock() override = default;
                                  
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
 *   If NumberHeatBlocks == 0 (there is no HeatBlock) then this variable
 *   need not be defined, since it is not loaded.
 *
 * - The variable "PowerHeatRho", of type double and indexed over the
 *   dimension "NumberUnits": entry PowerHeatRho[ i ] is assumed to contain
 *   the electrical-power-to-heat ratio for unit i. This is only used in
 *   the constraints linking electrical units to heat units, hence if
 *   NumberHeatBlocks == 0 (there is no HeatBlock) then this variable need
 *   not be defined, since it is not loaded.
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
/** The formulation of the problem and the methods that generate the abstract
 * constraint of the UCBlock.
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
 *  In the UCBlock there is not defined any variable but the decision
 *  variables here are present by using method get_variable() from UnitBlock,
 *  HeatBlock and NetworkBlock as follow:
 *
 * - \f$ p^{ac}_{t,i} \f$ : the active power variable for each time period
 *   \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I} \f$ is
 *   called from UnitBlock by get_active_power() method;
 *
 * - \f$ S_{t,n} \f$ : the node injection variable for each time period
 *   \f$ t \in \mathcal{T} \f$ and each node \f$ n \in \mathcal{N} \f$ is
 *   called from NetworkBlock by get_node_injection() method;
 *
 * - \f$ p^{pr}_{t,i} \f$ : the primary spinning reserves variable for each
 *   time period \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I}
 *   \f$ is called from UnitBlock by get_primary_spinning_reserve() method;
 *
 * - \f$ p^{sc}_{t,i} \f$ : the secondary spinning reserves variable for each
 *   time period \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I}
 *   \f$ is called from UnitBlock by get_secondary_spinning_reserve() method;
 *
 * - \f$ u_{t,i}  \in \{ 0 , 1 \} \f$ : the commitment state at time period
 *   \f$ t \in \mathcal{T} \f$ for each unit \f$ i \f$ is called from
 *   UnitBlock by get_commitment() method;
 *
 * - \f$ p^{he}_{t,i} \f$ : the heat variable for each time period
 *   \f$ t \in \mathcal{T} \f$ and each unit \f$ i \in \mathcal{I} \f$ is
 *   called from HeatBlock by get_heat() method;
 *
 *  The global constraints of unit commitment problem, on the time horizon
 *  \f$ \mathcal{T} \f$ write as follow:
 *
 * - Node injection Constraints:
 *   In the unit commitment problem, \f$ P^{au}_{t , i} \f$ denotes the fixed
 *   consumption of the power plant when it is off, and \f$ S_{t,n} \f$ is the
 *   node injection variable for each time period \f$ t \in \mathcal{T} \f$
 *   and each node \f$ n \in \mathcal{N} \f$. Therefore, if the
 *   f_number_nodes > 0, a boost::multi_array<FRowConstraint, 2>; with two
 *   dimensions which are f_time_horizon and f_number_nodes entries, where the
 *   entry t = 0, ..., f_time_horizon - 1 and the entry
 *   n = 1, ..., f_number_nodes being the node injection constraints at time
 *   t and node n as follow:
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
 *   be satisfied. Therefore, if the f_number_primary_zones > 0,
 *   a boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_primary_zones entries, where the entry
 *   t = 0, ..., f_time_horizon - 1 and the entry
 *   z = 1, ..., f_number_primary_zones being the primary demand constraints
 *   at time t and primary zone z as below;
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
 *   be satisfied. So, if the f_number_secondary_zones > 0,
 *   a boost::multi_array<FRowConstraint,2>; with two dimensions which are
 *   f_time_horizon, and f_number_secondary_zones entries, which the entry
 *   t = 0, ..., f_time_horizon - 1 and the entry
 *   z = 1, ..., f_number_secondary_zones being the secondary demand
 *   constraints at time t and secondary zone z as follow;
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
 *   satisfied. Therefore, if the f_number_inertia_zones > 0,
 *   a boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_inertia_zones entries, where the entry
 *   t = 0, ..., f_time_horizon - 1 and the entry
 *   z = 1, ..., f_number_inertia_zones being the inertia demand constraints
 *   at time t and inertia zone z as below;
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
 *   pollutant heat ratio respectively. So, if the f_number_pollutants > 0,
 *   a std::vector<std::vector<FRowConstraint>>; with two dimensions which are
 *   f_number_pollutants, and v_number_pollutant_zones entries, that the entry
 *   p = 1, ..., f_number_pollutants and the entry
 *   z = 1, ..., v_number_pollutant_zones being the pollutant budget
 *   constraints at pollutant p and pollutant zones z as below;
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
 *   \varrho_{i} \f$ of each heat producing unit \f$ i \f$. Therefore, if the
 *   f_number_heat_blocks > 0, a boost::multi_array<FRowConstraint, 2> with
 *   two dimensions which are f_time_horizon and the number of UnitBlocks that
 *   produce electricity and belong to some HeatBlock; the constraint at
 *   position ( t, i ) being the heat constraints at time t and unit M[ i ],
 *   where M maps the constraint into an electricity-producing unit that
 *   belongs to some HeatBlock. The Heat Constraints are defined as below:
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
 *
 * These methods allow to read data that must be common to (in principle) all
 * the related blocks to the UC problem.
 * @{ */

 /// returns the time horizon of the problem
 /** Method for returning the time horizon of the problem */
 Index get_time_horizon() const { return f_time_horizon; }

/*--------------------------------------------------------------------------*/
 /// returns the NetworkData object
 /** Method for returning the NetworkData object. Note that no NetworkData
  * may be defined (see comments to deserialize()), which means that the
  * transmission network is a "bus"; in this case, this method will return
  * nullptr. */

 NetworkData * get_NetworkData() const { return f_NetworkData; }

/*--------------------------------------------------------------------------*/
 /// returns the vector of (pointers to) NetworkBlocks
 /** Method for returning the vector of network blocks. This vector may have
  *  empty size which means there is no network block or may have size of
  *  f_time_horizon where in this case each element of the vector gives the
  *  network block at time instant t.
  * */
 const std::vector<NetworkBlock *> & get_network_blocks() const {
  return v_network_blocks;
  }

/*--------------------------------------------------------------------------*/
 /// method for returning the matrix of primary demand
 /** Method for returning the matrix of primary demand of the given primary
  *  zone at the time instant t. This matrix may have f_number_primary_zones
  *  rows and must have f_time_horizon columns, otherwise it's empty and there
  *  is no primary demand. Therefore, there are two possible cases:
  *
  * - if f_number_primary_zones  == 0 then the number of rows of this matrix
  *   is zero, in this case the matrix is empty since there is no primary
  *   zone.
  *
  * - if f_number_primary_zones != 0 , in this case each row must have size
  *   f_time_horizon and each element of that row gives the primary demand at
  *   time instant t.*/
 const boost::multi_array< double, 2 >& get_primary_demand() const {
  return( v_primary_demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of secondary demand
 /** Method for returning the matrix of secondary demand of the given
  *  secondary zone at the time instant t. This matrix may have
  *  f_number_secondary_zones rows and must have f_time_horizon columns,
  *  otherwise it's empty and there is no secondary demand. Therefore, there
  *  are two possible cases:
  *
  * - if f_number_secondary_zones  == 0 then the number of rows of this matrix
  *   is zero, in this case the matrix is empty since there is no secondary
  *   zone.
  *
  * - if f_number_secondary_zones != 0 , in this case each row must have size
  *   f_time_horizon and each element of that row gives the secondary demand
  *   at time instant t.*/
 const boost::multi_array< double, 2 > & get_secondary_demand() const {
  return( v_secondary_demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inertia demand of the given zone at the given time
 /** Method for returning the matrix of inertia demand of the given inertia
  *  zone at the time instant t. This matrix may have  f_number_inertia_zones
  *  rows and must have f_time_horizon columns, otherwise it's empty and there
  *  is no inertia demand. Therefore, there are two possible cases:
  *
  * - if f_number_inertia_zones  == 0 then the number of rows of this matrix
  *   is zero, in this case the matrix is empty since there is no inertia
  *   zone.
  *
  * - if f_number_inertia_zones != 0 , in this case each row must have size
  *   f_time_horizon and each element of that row gives the inertia demand
  *   at time instant t.*/
 const boost::multi_array< double, 2 > & get_inertia_demand() const {
  return( v_inertia_demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of conversion factor of the given pollutant
 /** Method for returning the matrix of conversion factor of the given
  *  pollutant due to the generation of the given unit at the given time t.
  *  This matrix considers two possible cases:
  *
  * - if f_time_horizon has size 1, then depending on the size of
  *   f_number_pollutants and f_number_units there exist these possibilities:
  *
  *   - if size of f_number_pollutants == 0 or f_number_units == 0 it means
  *     this matrix is empty.
  *   - if size of f_number_pollutants != 0 and or f_number_units != 0 each
  *     element of this matrix gives the conversion factor of pollutant due to
  *     the generation of unit i which for all time instants would be the same
  *     value.
  *
  * - if f_time_horizon has full size, then depending on the size of
  *   f_number_pollutants and f_number_units may happen these cases:
  *
  *   - if size of f_number_pollutants == 0 or f_number_units == 0 it means
  *     this matrix is empty.
  *   - if size of f_number_pollutants != 0 and or f_number_units != 0 each
  *     element of this matrix gives the conversion factor of pollutant due to
  *     the generation of unit i for each time instants t. */
 const boost::multi_array< double, 3 > & get_pollutant_rho() const {
  return(v_pollutant_rho);
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of conversion factor of the given pollutant due to the
 /// generation of every heat-only unit
 /** Method for returning the matrix of conversion factor of the given
  *  pollutant due to the pollutant due to the generation of every heat-only
  *  unit i in the given heat block h at the given time t.
  *  This matrix considers two possible cases:
  *
  * - if f_time_horizon has size 1, then depending on the size of
  *   f_number_pollutants and f_number_heat_blocks there exist these
  *   possibilities:
  *
  *   - if size of f_number_pollutants == 0 or f_number_heat_blocks == 0 it
  *     means this matrix is empty.
  *   - if size of f_number_pollutants != 0 and or f_number_heat_blocks != 0
  *     each element of this matrix gives the conversion factor of pollutant
  *     due to the generation of every heat-only unit in heat block h which
  *     for all time instants would be the same value.
  *
  * - if f_time_horizon has full size, then depending on the size of
  *   f_number_pollutants and f_number_heat_blocks may happen these cases:
  *
  *   - if size of f_number_pollutants == 0 or f_number_heat_blocks == 0 it
  *     means this matrix is empty.
  *   - if size of f_number_pollutants != 0 and or f_number_units != 0 each
  *     element of this matrix gives the conversion factor of pollutant due to
  *     the generation of every heat-only unit in HeatBlock h for each time
  *     instants t. */

 const boost::multi_array< double, 3 > & get_pollutant_heat_rho( ) const {
  return v_pollutant_heat_rho ;
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of pollutant zone
 /** Method for returning the matrix of pollutant zone associated with the
  *  given pollutant the given node belongs to. This matrix indexed over the
  *  dimensions NumberPollutants and NumberNodes. and there are these possible
  *  cases:
  *
  * - if each dimension of the matrix is not present then the matrix is empty
  *   and there is not any pollutant zone.
  *
  * - if the matrix has non-zero f_number_pollutants rows then:
  *
  *   - if f_number_nodes == 1, then the matrix is a vector with the
  *     f_number_pollutants rows and one column which implies that all the
  *     pollutants belong to just one node.
  *
  *   - if f_number_nodes > 1, then the matrix has f_number_pollutants rows
  *     and f_number_nodes columns where each element of the matrix shows
  *     which pollutant p belong to which node.
  *
  * - if the f_number_pollutants == 0, this matrix is empty and there is not
  *    any pollutant zone. */
 const boost::multi_array< Index, 2 > & get_pollutant_zone() const {
   return v_pollutant_zones;
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of heat set
 /** Method for returning the matrix of heat set that represents the given
  *  unit in the given HeatBlock. This matrix indexed over the dimensions
  *  NumberUnits and NumberHeatBlocks. and there are these possible cases:
  *
  * - if each dimension of the matrix is not present then the matrix is empty
  *   and there is not any heat set.
  *
  * - if the matrix has non-zero f_number_units rows then:
  *
  *   - if f_number_heat_blocks == 0, then this matrix is empty and there is
  *     not any heat set.
  *
  *   - if non-zero f_number_heat_blocks is present, then the matrix has
  *     f_number_units rows and f_number_heat_blocks columns, and each element
  *     of the matrix tells to which heat block h, electrical unit i belongs.
  *
  * - if the f_number_units == 0, this matrix is empty and there is not any heat set.
  *  */
 const boost::multi_array< Index, 2 > & get_heat_set() const {
   return v_heat_set;
 }

/*--------------------------------------------------------------------------*/
 /// returns the i-th UnitBlock
 /** Method for returning the unit block i. */
 UnitBlock * get_unit_block( Index i ) const;

/*--------------------------------------------------------------------------*/
 /// returns the t-th NetworkBlock
 /** Method for returning the network block at time instant t. */
 NetworkBlock * get_network_block( Index t ) const;

/*--------------------------------------------------------------------------*/
 /// returns the i-th HeatBlock
 /** Method for returning the heat block h. */
 HeatBlock * get_heat_block( Index h ) const;

/*--------------------------------------------------------------------------*/
 /// returns the node where the given unit belongs to
 /** Method for returning the node where the given unit belongs to.
  *
  *  - if f_number_units is not present or its present but f_number_units == 0
  *    this vector is empty.
  *
  *  - if f_number_units >= 1, then this vector have size f_number_units and
  *    there are two possible cases:
  *
  *    - if the NetworkData is not defined (see comments to deserialize()) or
  *      it's defined and f_number_nodes == 1 then transmission network is a
  *      "bus", then all the elements of this vector belong to that node.
  *
  *    - if the NetworkData is defined (see comments to deserialize()) and
  *      f_number_nodes > 1 then each element of this vector gives which unit
  *      belongs to which node.
  * */
 Index get_unit_node( Index unit ) const {
   if( f_NetworkData ? f_NetworkData->get_number_nodes() : 1 )
     return v_unit_node[ unit ];
   return 0;
 }

/*--------------------------------------------------------------------------*/
 /// returns the primary zone where the given node belongs to
 /** Method for returning the vector of primary zones. This vector may have
  *  empty size which means there is no primary zone or may have size of
  *  f_number_nodes where in this case each element of the vector implies the
  *  primary zone at node n.
  * */
 Index get_primary_zone( Index node ) const {
  return( v_primary_zones.empty() ? 0 : v_primary_zones[ node ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the secondary zone where the given node belongs to
 /** Method for returning the vector of secondary zones. This vector may have
  *  empty size which means there is no secondary zone or may have size of
  *  f_number_nodes where in this case each element of the vector implies the
  *  secondary zone at node n.
  * */
 Index get_secondary_zone( Index node ) const {
  return( v_secondary_zones.empty() ? 0 : v_secondary_zones[ node ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia zone where the given node belongs to
 /** Method for returning the vector of inertia zones. This vector may have
  *  empty size which means there is no inertia zone or may have size of
  *  f_number_nodes where in this case each element of the vector implies the
  *  inertia zone at node n.
  * */
 Index get_inertia_zone( Index node ) const {
  return( v_inertia_zones.empty() ? 0 : v_inertia_zones[ node ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the electrical-power-to-heat ratio of the given unit
 /** Method for returning the vector of electrical power to heat. This vector
  *  may have empty size which means there is no inertia zone or may have size
  *  of f_number_nodes where in this case each element of this vector implies
  *  the inertia zone at node n.
  * */
 double get_power_heat_rho( Index unit ) const {
  return v_power_heat_rho[ unit ];
  }
/*--------------------------------------------------------------------------*/
 /// returns the vector of electrical-power-to-heat ratio
 /** Method for returning the vector electrical-power-to-heat ratio for each
  *  unit i of heat block h. This vector may have  empty size which means
  *  there no power heat rho or may have size of f_number_units where in this
  *  case each element of this vector gives the power heat rho for each heat
  *  block h. Note that if the f_number_heat_blocks is not present or it's
  *  equal to zero, this vector is empty.
  * */
 const std::vector< double > & get_power_heat_rho( void ) const {
  return v_power_heat_rho ;
  }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UCBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the UCBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 *  UCBlock. See UCBlock::deserialize( netCDF::NcGroup ) for details of the
 *  format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UCBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UCBlock
    @{ */

 virtual void load( std::istream &input ) override {
  throw( std::logic_error( "UCBlock::load() not implemented yet" ) );
  }

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 /// The time horizon of the problem
 Index f_time_horizon;

 /// The number of units of the problem
 Index f_number_units;

 /// the NetworkData object
 NetworkData * f_NetworkData;

 /// The number of heat block
 Index f_number_heat_blocks;

 /// The number of nodes in primary zones of the network
 Index f_number_primary_zones;

 /// The number of nodes in secondary zones of the network
 Index f_number_secondary_zones;

 /// The number of nodes in inertia zones of the network
 Index f_number_inertia_zones;

 /// The number of pollutants
 Index f_number_pollutants;

 /// The set of UnitBlocks
 std::vector<UnitBlock *> v_unit_blocks;

 /// The set of HeatBlocks
 std::vector<HeatBlock *> v_heat_blocks;

 /// The number of pollutant zones of each pollutant
 std::vector<Index> v_number_pollutant_zones;

 /** The matrix PollutantZones indexed over the dimensions
  *  NumberPollutants and NumberNodes */
 boost::multi_array< Index, 2 > v_pollutant_zones;

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
 boost::multi_array< double, 2 > v_primary_demand;

 /// the vector of SecondaryZones
 std::vector<Index> v_secondary_zones;

 /** the matrix of SecondaryDemand indexed over the dimensions
  * SecondaryZones and TimeHorizon */
 boost::multi_array< double, 2 > v_secondary_demand;

 /// the vector InertiaZones
 std::vector<Index> v_inertia_zones;

 /** the matrix of InertiaDemand indexed over the dimensions
  * InertiaZones and TimeHorizon */
 boost::multi_array< double, 2 >v_inertia_demand;

 /// the vector of PollutantDemand
 std::vector< double > v_pollutant_budget;

 /** the PollutantRho matrix index over the dimensions
  * TimeHorizon, NumberPollutants, and NumberUnits */
 boost::multi_array< double, 3 > v_pollutant_rho;

 /** the PollutantHeatRho matrix index over the dimensions
  * TimeHorizon, NumberPollutants, and NumberHeatBlocks */
 boost::multi_array< double, 3 > v_pollutant_heat_rho;

 /// The entry v_unit_node[ i ] tells to which node unit i belongs
 std::vector<Index> v_unit_node;

 /// The entry v_heat_node[ h ] tells to which node the HeatBlock h belongs
 std::vector<Index> v_heat_node;

 /** the HeatSet matrix index over the dimensions
  * NumberUnits, and NumberHeatBlocks */
 boost::multi_array< Index, 2 > v_heat_set;

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

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

  };   // end( class( UCBlock ) )

/*--------------------------------------------------------------------------*/

} /* namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* UCBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File UCBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
