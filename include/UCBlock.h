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
 * - A multi_array of primary demand constraints.
 *
 * - A multi_array of secondary demand constraints.
 *
 * - A multi_array of inertia demand constraints.
 *
 * - A multi_array of pollutant budget constraints.
 *
 * - A multi_array of heat constraints.
 *
 * - Some basic public methods to read the data and initialize the
 *   optimisation problem.
 *
 * \version 0.11
 *
 * \date 03 - 06 - 2019
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

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS UCBlock --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// implementation of the Block concept for the unit commitment problem
/** The UCBlock class implements the Block concept [see Block.h] for the EDF
 *  Unit Commitment Problem.
 *
 *  The unit commitment corresponds to the short-term problem (with for
 *  instance one week or one day time horizon or even longer when the SSV is
 *  not considered and feasible regarding computational limitations) arising
 *  at each stage of the mid-term problem. This consists of a deterministic
 *  optimization problem involving a large number of power plants, load
 *  flexibilities and storage devices intended to jointly satisfy a
 *  deterministic demand at each node of the network. The goal of these
 *  problems is to find an optimal (or near optimal) schedule satisfying the
 *  set of technical constraints.
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
 *   \rho^{'}_{t , p , i} \f$ where considered as pollutant ratio and
 *   pollutant heat ratio respectively; is defined  as follow:
 *
 * \f[
 *  \sum_{n \in \mathcal{B}}\sum_{ t \in \mathcal{T} }\sum_{ i \in
 *  \mathcal{I}_n } \rho_{t , p , i} p^{ac}_{t,i} + \sum_{n \in \mathcal{B'}}
 *  \sum_{ t \in \mathcal{T} }\sum_{ i \in \mathcal{I}(h)
 *  } \rho^{'}_{t , p , h} p^{he}_{t,i}  \leq \mathcal{O}_p  \quad
 *  \mathcal{B} \in \mathcal{B}^{p}(\mathcal{N})  \quad
 *  \mathcal{B'} \in \mathcal{B}^{p}(\mathcal{H})
 *          \quad h \in \mathcal{H} \quad p \in \mathcal{P} \quad          (5)
 * \f]
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
 * - the dimension "NumberHeatBlocks" containing the number of heat blocks in
 *   in the problem; the dimension is optional: if it is not provided then it
 *   is taken to be 0, which means that there is no heat block in the problem;
 *
 * - the groups "UnitBlock_0", "UnitBlock_1", ... , "UnitBlock_n" with
 *   n == NumberUnits - 1, containing each one UnitBlock corresponding
 *   to one unit;
 *
 * - the groups "HeatBlock_0", "HeatBlock_1", ... , "HeatBlock_n" with
 *   n == NumberHeatBlocks - 1, containing each one a HeatBlock corresponding
 *   to one energy cell; when NumberHeatBlocks == 0, there is no heat
 *   constraint anywhere in the problem;
 *
 * - the dimension "NumberNodes" containing the number of nodes in
 *   the problem; the dimension is optional, if it is not provided then it is
 *   taken to be 1, which means that all the Unit belong to the same node
 *   (the network is a bus);
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
 * - the variable "HeatSet", of type int and indexed both over the dimensions
 *   "NumberUnits" and "NumberHeatBlocks"; if HeatSet[ i , h ] = k, with
 *   k < number of heat units in HeatBlock h, then electrical unit i is
 *   represented into HeatBlock h and the heat unit k; if, instead,
 *   HeatSet[ i , h ] = k, with k >= number of heat units in HeatBlock h,
 *   then none of the heat units in HeatBlock h represents the electrical unit
 *   i; if NumberHeatBlocks == 0 (say, there is no HeatBlock) then
 *   this variable need not be defined, since it is not loaded;
 *
 * - the variable "PowerHeatRho", of type double and indexed over the
 *   dimensions "NumberUnits": entry PowerHeatRho[ i ] is assumed to
 *   contain the electrical-power-to-heat ratio for each unit i;
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
 *   dimension "NumberNodes"; the entry InertiaZones[ n ] tells to
 *   which inertia zone the node n belongs; if InertiaZones[ n ] >=
 *   NumberInertiaZones, this means that node n does not belong to any
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
 *   dimension "NumberPollutants"; the p-th entry of the variable is
 *   assumed to contain the number of pollutant zones associated with
 *   pollutant p; if NumberPollutants == 0 (say, there is no
 *   pollutant) then this variable need not be defined, since it is
 *   not loaded.
 *
 * - the variable "PollutantZones", of type int and indexed over the
 *   dimensions "NumberPollutants" and "NumberNodes"; the entry
 *   PollutantZones[ p , n ] tells to which pollutant zone associated
 *   with pollutant p the node n belongs; if PollutantZones[ p , n ]
 *   >= NumberPollutantZones[ p ], this means that node n does not
 *   belong to any pollutant zone, and hence the corresponding units
 *   are not involved into the pollutant demand constraints associated
 *   with pollutant p; if NumberPollutants == 0 (say, there is no
 *   pollutant) then this variable need not be defined, since it is
 *   not loaded;
 *
 * - the variable "PollutantHeatZones", of type int and indexed over the
 *   dimensions "NumberPollutants" and "NumberHeatBlocks"; the entry
 *   PollutantHeatZones[ p , h ] tells to which pollutant zone associated
 *   with pollutant p, the HeatBlock h belongs; if PollutantHeatZones[ p , h ]
 *   >= NumberPollutantZones[ p ], this means that HeatBlock h does not
 *   belong to any pollutant zone, and hence the corresponding units
 *   are not involved into the pollutant budget constraints associated
 *   with pollutant p; if NumberPollutants == 0 (say, there is no
 *   pollutant) then this variable need not be defined, since it is
 *   not loaded; if NumberHeatBlocks == 0 (say, there is no heat-only unit)
 *   then this variable need not be defined, since it is not loaded;
 *
 * - the variable "PollutantBudget", of type double and indexed over
 *   the dimension "NumberPollutants"; the i-th entry of the variable
 *   is assumed to contain the limit of pollutant i; if
 *   NumberPollutants == 0 (say, there is no pollutant) then this
 *   variable need not be defined, since it is not loaded;
 *
 * - the variable "PollutantRho", of type double and indexed over the
 *   dimensions "TimeHorizon", "NumberPollutants", and "NumberUnits";
 *   the entry PollutantRho[ t , p , i ] is assumed to contain the
 *   conversion factor of pollutant p due to the generation of unit i
 *   at time t; if NumberPollutants == 0 (say, there is no pollutant)
 *   then this variable need not be defined, since it is not loaded;
 *
 * - the variable "PollutantHeatRho", of type double and indexed over the
 *   dimensions "TimeHorizon", "NumberPollutants", and "NumberHeatBlocks"; the
 *   entry PollutantRho[ t , p , h ] is assumed to contain the conversion
 *   factor of pollutant p due to the generation of every heat-only unit in
 *   HeatBlock h at time t. If NumberPollutants == 0 (say, there is no
 *   pollutant) then this variable does not need be defined, since it is not
 *   loaded; if NumberHeatBlocks == 0 (say, there is no heat-only unit) then
 *   this variable need not be defined, since it is not loaded;
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

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
  * - if f_number_heat_block > 0, a boost::multi_array<FRowConstraint, 2>;
  *   with two dimensions which are f_time_horizon, and f_number_units
  *   entries, the entry t = 0, ..., f_time_horizon - 1 and the entry
  *   i = 1, ..., f_number_units - 1 being the heat constraints
  *   at time t and unit i;  */

 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
   override;

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
 *  @{ */

 /// Method that initializes the instance and passes all the needed data
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

/** returns the conversion factor of the given pollutant due to the generation
 * of the given heat-only unit of each HeatBlock at the given time */
 inline double get_pollutant_heat_rho( Index time, Index pollutant, Index unit )
   const {
   auto index = time * f_number_pollutants * f_number_heat_block +
          pollutant * f_number_heat_block + unit;
   return v_pollutant_heat_rho[ index ];
 }

 /** returns the pollutant zone associated with the given pollutant
  * the given node belongs to */
 inline Index get_pollutant_zone( Index pollutant, Index node ) const {
   return v_pollutant_zones[ pollutant * f_number_nodes + node ];
 }

 /** returns the HeatSet associated with the given unit
  * the given HeatBlocks belongs to */
 inline Index get_heat_set( Index unit, Index block ) const {
     return v_heat_set[ f_number_units * f_number_heat_block + block ];
 }

 /** returns the pollutant heat zone associated with the given pollutant
  * the given node belongs to */
 inline Index get_pollutant_heat_zone( Index pollutant, Index block ) const {
     return v_pollutant_heat_zones[ pollutant * f_number_heat_block + block ];
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

/// returns the electrical-power-to-heat ratio of the given unit i
inline double get_power_heat_rho( Index unit ) const {
 return v_power_heat_rho[ unit * f_number_units];
}

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UCBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the UCBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * UCBlock. See UCBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

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

 /// The number of heat block
 Index f_number_heat_block;

 /// The entry v_node[ i ] tells to which node unit i belongs
 std::vector<Index> v_node;

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

 /** The matrix PollutantHeatZones indexed over the dimensions
  *  NumberPollutants and NumberHeatBlocks */
  std::vector<Index> v_pollutant_heat_zones;

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
                              const std::string sub_group_name_prefix,
                              const int num_sub_blocks );

};   // end( class( UCBlock ) )

} /* namespace SMSpp_di_unipi_it */

#endif /* UCBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File UCBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
