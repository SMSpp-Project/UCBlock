/*--------------------------------------------------------------------------*/
/*--------------------------- File UCBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class UCBlock, which implements the
 * base class Block, in order to define a very basic Unit Commitment
 * Block, that will be able to fit and be used as a base for almost
 * any different variation of the Unit Commitment Problem. As a
 * result of this, the UCBlock class is characterized by the following
 * ingredients:
 *
 * - A time horizon of the optimization problem.
 *
 * - A set of units (UnitBlock; that may be referring to any different
 *   kind of unit, such as thermal, hydro, etc).
 *
 * - A set of networks (NetworkBlock; that can potentially refer to
 *   demand satisfaction, transmission line capacities, etc.). This
 *   set is either empty, which means that there is no network in the
 *   model, or has size equals the time horizon.
 *
 * - A multi_array of node injection constraints.
 *
 * - Some basic public methods to read the data and initialize the
 *   optimisation problem.
 *
 * \version 0.11
 *
 * \date 21 - 05 - 2019
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
#ifndef __UCBlock
#define __UCBlock
/* self-identification: #endif at the end of the file */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <boost/multi_array.hpp>
#include <vector>
#include "Block.h"
#include "NetworkBlock.h"
#include "UnitBlock.h"
#include "HeatBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class FRowConstraint; ///< forward declaration of FRowConstraint
class NetworkBlock;   ///< forward declaration of NetworkBlock
class UnitBlock;      ///< forward declaration of UnitBlock

class UCBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */

 typedef unsigned int Index;

/*@} -----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of UCBlock, taking possibly a pointer of its father Block
 UCBlock( Block *father = nullptr ) : Block( father ) {}

/*--------------------------------------------------------------------------*/

 virtual ~UCBlock() {};   ///< destructor of UCBlock: it is virtual, and empty

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/*--------------------------------------------------------------------------*/
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific
 * format of the UCBlock. Besides the mandatory "type" attribute of
 * any :Block, the group should contain the following:
 *
 * - the dimension "TimeHorizon" containing the number of time steps in
 *   the problem;
 *
 * - the dimension "NumberUnits" containing the number of units in
 *   the problem;
 *
 * - the dimension "NumberHeatBlocks" containing the number of
 *   heat blocks in in the problem; the dimension is
 *   optional: if it is not provided then it is taken to be 0, which
 *   means that there is no heat block in the problem;
 *
 * - the groups "UnitBlock_0", "UnitBlock_1", ... , "UnitBlock_n" with
 *   n == NumberUnits - 1, containing each one UnitBlock corresponding
 *   to one unit;
 *
 * - the groups "HeatBlock_0", "HeatBlock_1", ... , "HeatBlock_n" with
 *   n == NumberHeatBlocks - 1, containing each one HeatBlock corresponding to
 *   one energy cell; when NumberHeatBlocks == 0, there is no heat constraint
 *   anywhere in the problem;
 *
 * - the dimension "NumberNodes" containing the number of nodes in
 *   the problem; the dimension is optional, if it is not provided then it is
 *   taken to be 1, which means that all the Unit belong to the same node
 *   (the network is a bus);
 *
 * NO, THIS WE DON'T NEED
 * - the dimension "NumberHeatOnlyUnits" containing the number of
 *   heat-only generation units in the problem; the dimension is
 *   optional: if it is not provided then it is taken to be 0, which
 *   means that there is no heat-only generation unit;
 *
 * - the groups "NetworkBlock_0", "NetworkBlock_1", ... , "NetworkBlock_t"
 *   with t = TimeHorizon - 1, containing each the state of the interconnect
 *   network at time t;
 *
 * - the variable "Node", of type int and indexed over the dimension
 *   "NumberUnits"; the entry Node[ i ] tells to which node unit i belongs;
 *   if NumberNodes == 1  (say, it is not provided at all), then this
 *   variable need not be defined, since it is not loaded;
 *
 * NO, THIS WE DON'T NEED
 * - the variable "HeatOnlyUnits", of type int and indexed over the
 *   dimension "NumberHeatOnlyUnits"; it contains the indices of the
 *   heat-only generation units; if NumberHeatOnlyUnits == 0, then
 *   this variable need not be defined, since there is no heat-only
 *   generation unit and, therefore, this variable is not loaded;
 *
 * WE NEED A MATRIX THAT, FOR EACH HEAT UNIT OF EACH HEAT BLOCK, TELLS
 * IF THIS IS A HEAT-ONLY UNIT (INDEX >= NumberUnits), OR IF THIS IS
 * ALSO AN ELECTRICITY-PRODUCING UNIT (INDEX < NumberUnits)
 *
 * NO, THIS WE DON'T NEED
 * - the variable "UnitEnergyCell", of type int and indexed over the dimension
 *   "NumberHeatBlocks"; it contains the indices of the heat blocks; if
 *   NumberHeatBlocks == 0, then this variable need not be defined, since
 *   there is heat block and, therefore, this variable is not loaded;
 *
 * - the dimension "NumberPrimaryZones" is associated with one specific primary
 *   spinning reserve in the problem. The dimension is optional, if it is not
 *   provided then it is taken to be 0, which means that no primary reserve
 *   constraints are presented in the problem;
 *
 * - the variable "PrimaryZones", of type int and indexed over the dimension
 *   "NumberNodes"; the entry PrimaryZones[ i ] tells to which primary zone
 *   the node i belongs; if PrimaryZones[ i ] >= NumberPrimaryZones, this
 *   means that node i does not belong to any primary zone, and hence the
 *   corresponding units are not involved into the primary reserve
 *   constraints; if NumberPrimaryZones == 0  (say, it is not provided
 *   at all) then this variable need not be defined, since it is not loaded;
 *   if NumberPrimaryZones == 1 and this variable is not defined, then there
 *   is only one primary zone and all the nodes belong to it;
 *
 * - the variable "PrimaryDemand", of type double and indexed both
 *   over the dimensions "PrimaryZones" and "TimeHorizon": entry
 *   PrimaryDemand[ i , t ] is assumed to contain the primary reserves
 *   requirement which are specified on the primary reserves zone i in
 *   the time t; if NumberPrimaryZones == 0 (say, it is not provided
 *   at all), then this variable need not be defined, since it is not
 *   loaded;
 *
 * - the dimension "NumberSecondaryZones" is associated with one
 *   specific secondary spinning reserve in the problem. The dimension
 *   is optional, if it is not provided then it is taken to be 0,
 *   which means that no secondary reserve constraints are presented in
 *   the problem;
 *
 * - the variable "SecondaryZones", of type int and indexed over the
 *   dimension "NumberNodes"; the entry SecondaryZones[ i ] tells to
 *   which secondary zone the node i belongs; if SecondaryZones[ i ]
 *   >= NumberSecondaryZones, this means that node i does not belong
 *   to any secondary zone, and hence the corresponding units are not
 *   involved into the secondary reserve constraints; if
 *   NumberSecondaryZones == 0 (say, it is not provided at all) then
 *   this variable need not be defined, since it is not loaded; if
 *   NumberSecondaryZones == 1 and this variable is not defined, then
 *   there is only one secondary zone and all the nodes belong to it;
 *
 * - the variable "SecondaryDemand", of type double and indexed both
 *   over the dimensions "SecondaryZones" and "TimeHorizon": entry
 *   SecondaryDemand[ i , t ] is assumed to contain the secondary
 *   reserves requirement which are specified on the secondary
 *   reserves zone i in the time t; if NumberSecondaryZones == 0 (say,
 *   it is not provided at all), then this variable need not be
 *   defined, since it is not loaded;
 *
 * - the dimension "NumberInertiaZones" is associated with one
 *   specific inertia zone in the problem. The dimension is optional,
 *   if it is not provided then it is taken to be 0;
 *
 * - the variable "InertiaZones", of type int and indexed over the
 *   dimension "NumberNodes"; the entry InertiaZones[ i ] tells to
 *   which inertia zone the node i belongs; if InertiaZones[ i ] >=
 *   NumberInertiaZones, this means that node i does not belong to any
 *   inertia zone, and hence the corresponding units are not involved
 *   into the inertia reserve constraints; if NumberInertiaZones == 0
 *   (say, it is not provided at all) then this variable need not be
 *   defined, since it is not loaded; if NumberInertiaZones == 1 and
 *   this variable is not defined, then there is only one inertia zone
 *   and all the nodes belong to it;
 *
 * - the variable "InertiaDemand", of type double and indexed both
 *   over the dimensions "InertiaZones" and "TimeHorizon": entry
 *   InertiaDemand[ i , t ] is assumed to contain the inertia reserves
 *   requirement which are specified on the inertia reserves zone i in
 *   the time t; if NumberInertiaZones == 0 (say, it is not provided
 *   at all), then this variable need not be defined, since it is not
 *   loaded;
 *
 * - the dimension "NumberPollutants" containing the number of
 *   pollutants in the problem. The dimension is optional, if it is
 *   not provided then it is taken to be 0;
 *
 * - the variable "NumberPollutantZones" of type int indexed over the
 *   dimension "NumberPollutants"; the i-th entry of the variable is
 *   assumed to contain the number of pollutant zones associated with
 *   pollutant i. If NumberPollutants == 0 (say, there is no
 *   pollutant) then this variable need not be defined, since it is
 *   not loaded.
 *
 * - the variable "PollutantZones", of type int and indexed over the
 *   dimensions "NumberPollutants" and "NumberNodes"; the entry
 *   PollutantZones[ i , j ] tells to which pollutant zone associated
 *   with pollutant i the node j belongs; if PollutantZones[ i , j ]
 *   >= NumberPollutantZones[ i ], this means that node j does not
 *   belong to any pollutant zone, and hence the corresponding units
 *   are not involved into the pollutant demand constraints associated
 *   with pollutant i; if NumberPollutants == 0 (say, there is no
 *   pollutant) then this variable need not be defined, since it is
 *   not loaded;
 *
 * EITHER WE NEED A NEW VARIABLE
 * - the variable "PollutantHeatZones", of type int and indexed over the
 *   dimensions "NumberPollutants" and "NumberHeatBlocks"; the entry
 *   PollutantHeatZones[ i , j ] tells to which pollutant zone associated
 *   with pollutant i the HeatBlock j belongs; if PollutantHeatZones[ i , j ]
 *   >= NumberPollutantZones[ i ], this means that HeatBlock j does not
 *   belong to any pollutant zone, and hence the corresponding units
 *   are not involved into the pollutant demand constraints associated
 *   with pollutant i; if NumberPollutants == 0 (say, there is no
 *   pollutant) then this variable need not be defined, since it is
 *   not loaded;
 * OR WE NEED TO CONCATENATE PollutantZones and PollutantHeatZones, IF THIS
 * IS POSSIBLE IN netCDF
 *
 * DIFFERENT NAME
 * - the variable "PollutantBudget", of type double and indexed over
 *   the dimension "NumberPollutants"; the i-th entry of the variable
 *   is assumed to contain the limit of pollutant i; if
 *   NumberPollutants == 0 (say, there is no pollutant) then this
 *   variable need not be defined, since it is not loaded;
 *
 * - the variable "PollutantRho", of type double and indexed over the
 *   dimensions "TimeHorizon", "NumberPollutants", and "NumberUnits";
 *   the entry PollutantRho[ t , i , j ] is assumed to contain the
 *   conversion factor of pollutant i due to the generation of unit j
 *   at time t. If NumberPollutants == 0 (say, there is no pollutant)
 *   then this variable need not be defined, since it is not loaded;
 *
 * HERE WE NEED THE SAME INFORMATION FOR HEAT BLOCKS: FOR EACH HEAT-ONLY
 * UNIT IN EACH HEAT BLOCK (IF IT APPEARS IN A POLUTANT), WE NEED THE
 * COEFFICIENT IN THE POLLUTANT CONSTRAINT
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
   override;

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
 *  @{ */

 /// Method that initiazes the instance and passes all the needed data
 void instance( std::istream& inStream );

 /// returns the time horizon of the problem
 Index get_time_horizon( void ) const { return f_time_horizon; }

 /// returns the vector of (pointers to) NetworkBlocks
 const std::vector<NetworkBlock *> & get_network_blocks( void ) const {
   return v_network_blocks;
 }

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

 /** returns the pollutant zone associated with the given pollutant
  * the given node belongs to */
 inline Index get_pollutant_zone( Index pollutant, Index node ) const {
   return v_pollutant_zones[ pollutant * f_number_nodes + node ];
 }

 /// returns the i-th UnitBlock
 inline UnitBlock * get_unit_block( Index i ) const {
   return static_cast<UnitBlock *>( v_Block[ i ] );
 }

 /// returns the t-th NetworkBlock
 inline NetworkBlock * get_network_block( Index t ) const {
   return static_cast<NetworkBlock *>( v_Block[ f_number_units + t ] );
 }

 /// returns the node where the given unit belongs to
 inline Index get_node( Index unit ) const {
   if( f_number_nodes > 1 )
     return v_node[ unit ];
   return 0;
 }

 /// returns the primary zone where the given node belongs to
 inline Index get_primary_zone( Index node ) const {
   if( v_primary_zones.size() > 0 )
     return v_primary_zones[ node ];
   return 0;
 }

 /// returns the secondary zone where the given node belongs to
 inline Index get_secondary_zone( Index node ) const {
   if( v_secondary_zones.size() > 0 )
     return v_secondary_zones[ node ];
   return 0;
 }

/// returns the inertia zone where the given node belongs to
inline Index get_inertia_zone( Index node ) const {
    if( v_inertia_zones.size() > 0 )
      return v_inertia_zones[ node ];
    return 0;
}

/// returns the i-th HeatBlock
inline HeatBlock * get_heat_block( Index i ) const {
    return static_cast<HeatBlock *>( v_Block[ i ] );
}

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UCBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the UCBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * UCBlock. See UCBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 *
 * */

 virtual void serialize( netCDF::NcGroup & group ) const override final;

/*@} -----------------------------------------------------------------------*/
/*------------------ METHODS FOR MODIFYING THE UCBlock ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the UCBlock
 *  @{ */

 /// sets the time horizon of the problem
 void set_time_horizon( int t ) { f_time_horizon = t; }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

 /// The time horizon of the problem
 Index f_time_horizon;

 /// The number of units of the problem
 Index f_number_units;

 /// The number of nodes in the network
 Index f_number_nodes;

 /// The number of heat-only generation units
 Index f_number_heat_only_units;

/// The number of heat block
Index f_number_heat_block;

 /// The entry v_node[ i ] tells to which node unit i belongs
 std::vector<Index> v_node;

 /// Contains the indices of the heat-only generation units
 std::vector<Index> v_heat_only_units;

/// Contains the indices of the energy cell
std::vector<Index> v_energy_cell;

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
 std::vector<double> v_inertia_demand;

 /// the vector of PollutantDemand
 std::vector<double> v_pollutant_demand;

 /** the PollutantRho matrix index over the dimensions
  * TimeHorizon, NumberPollutants, and NumberUnits */
 std::vector<double> v_pollutant_rho;

 /// Node injection constraints for each time and node
 boost::multi_array<FRowConstraint, 2> v_node_injection_constraints;

 /// Primary demand constraints for each time and primary zone
 boost::multi_array<FRowConstraint, 2> v_PrimaryDemand_Const;

 /// Secondary demand constraints for each time and secondary zone
 boost::multi_array<FRowConstraint, 2>  v_SecondaryDemand_Const;

 /// Inertia demand constraints for each time and inertia zone
 boost::multi_array<FRowConstraint, 2>  v_InertiaDemand_Const;

/// heat constraints for each time and index unit
boost::multi_array<FRowConstraint, 2>  v_Heat_Const;

 /// Pollutant demand constraints for each pollutant and pollutant zone
 std::vector<std::vector<FRowConstraint>> v_PollutantDemand_Const;

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
                              const std::string sub_group_name_prefix,
                              const int num_sub_blocks );

};   // end( class( UCBlock ) )

} /* namespace SMSpp_di_unipi_it */

#endif /* UCBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File UCBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
