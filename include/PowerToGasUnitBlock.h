/*--------------------------------------------------------------------------*/
/*----------------- File PowerToGasUnitBlock.h -----------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class PowerToGasUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" Power-to-Gas
 * unit at unit commitment Problem.
 *
 * \version 0.11
 *
 * \date 21 - 08 - 2019
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
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PowerToGasUnitBlock
#define __PowerToGasUnitBlock
                     /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "UnitBlock.h"
#include "FRowConstraint.h"
#include "FRealObjective.h"
#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// Namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS PowerToGasUnitBlock ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the power to gas unit problem
/** The PowerToGasUnitBlock class implements the Block concept [see Block.h]
 * for a "reasonably standard" PowerToGas unit of a unit commitment Problem.
 * That is, the class is designed in order to give mathematical formulation to
 * describe the operation of large set of PowerToGasUnitBlock. Besides the
 * coupling of electricity and heat, the gas sector is another way of coupling
 * energy sectors. To do so, power-to-gas units use the surplus of electricity
 * at times of high renewable generation to disassemble water into hydrogen
 * and oxygen. Storing hydrogen or using the gas grid as an alternative energy
 * infrastructure provides further flexibility to the overall energy system.
 * To ensure the feasibility of extra gas inflow to the gas grid, the final
 * power-to-gas schedules (re-dispatched schedules after the transmission grid
 * operation model) are validated by means of the gas network model.
 * The technical and physical constraints are mainly divided in four different
 * categories:
 * - maximum and minimum level of power output constraints;
 *
 * - ramp-up and ramp-down constraints;
 *
 * - active power relation with storing energy levels constraints;
 *
 * - storage level constraints;
 */
class PowerToGasUnitBlock : public UnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// constructor, takes the father and the time horizon
/** Constructor of PowerToGasUnitBlock, taking possibly a pointer of its
 * father Block.
 */

 explicit PowerToGasUnitBlock( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/

/// destructor of PowerToGasUnitBlock

 ~PowerToGasUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the PowerToGasUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group must contain all the data required by the base UnitBlock,
 * as described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:

 * - The variable "MinStorage", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MinS[ t ] that, for each time instant t, contains the minimum
 *   storage level of the unit for the corresponding time step. If
 *   "MinStorage" has length 1 then MinS[ t ] contains the same value for all
 *   t. Otherwise, MinStorage[ i ] is the fixed value of MinS[ t ] for all t
 *   in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must be
 *   always that MinS[ t ] >= 0, and MinS[ t ] <= MaxS[ t ] for all t. If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxStorage", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MaxS[ t ] that, for each time instant t, contains the maximum
 *   storage level of the unit for the corresponding time step. If
 *   "MaxStorage" has length 1 then MaxS[ t ] contains the same value for all
 *   t. Otherwise, MaxStorage[ i ] is the fixed value of MaxS[ t ] for all t
 *   in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must be
 *   always that MaxS[ t ] >= 0, and MinS[ t ] <= MaxS[ t ] for all t. If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MinPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MinP[ t ] that, for each time instant t, contains the minimum
 *   active power output value of the unit for the corresponding time step.
 *   If "MinPower" has length 1 then MinP[ t ] contains the same value for all
 *   t. Otherwise, MinPower[ i ] is the fixed value of MinP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must be
 *   MinP[ t ] <= MaxP[ t ] for all t. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MaxP[ t ] that, for each time instant t, contains the maximum
 *   active power output value of the unit for the corresponding time step.
 *   If "MaxPower" has length 1 then MaxP[ t ] contains the same value for all
 *   t. Otherwise, MaxPower[ i ] is the fixed value of MaxP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must be
 *   MaxP[ t ] >= MinP[ t ] for all t. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The scalar variable "InitialPower", of type double and not indexed over
 *   any dimension. This variable indicates the amount of the power that the
 *   unit was producing at time instant -1, i.e., before the start of the time
 *   horizon; this is necessary to compute the ramp-up and ramp-down
 *   constraints. This variable is optional; if "DeltaRampUp" and
 *   "DeltaRampDown" are not present, "InitialPower" should not be read and it
 *   means there are no ramping constraints. If "DeltaRampUp" and
 *   "DeltaRampDown" are present but "InitialPower" is not provided, its
 *   initial value is taken to be 0.
 *
 * - The variable "DeltaRampUp", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector DP[ t ] that, for each time instant t, contains the ramp-up value
 *   of the unit for the corresponding time step, i.e., the maximum possible
 *   increase of active power production w.r.t. the power that had been
 *   produced in time instant t - 1, if any. This variable is optional; if it
 *   is not provided then it is assumed that DP[ t ] == MaxP[ t ], i.e., the
 *   unit can ramp up by an arbitrary amount, i.e., there are no ramp-up
 *   constraints. If "DeltaRampUp" has length 1 then DP[ t ] contains the same
 *   value for all t. Otherwise, DeltaRampUp[ i ] is the fixed value of DP[ t ]
 *   for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
 *   0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "DeltaRampDown", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector DM[ t ] that, for each time instant t, contains the ramp-down
 *   value of the unit for the corresponding time step, i.e., the maximum
 *   possible decrease of active power production w.r.t. the power that had
 *   been produced in time instant t - 1, if any. This variable is optional;
 *   if it is not provided then it is assumed that DP[ t ] == MaxP[ t ], i.e.,
 *   the unit can ramp down an arbitrary amount, i.e., there are no
 *   ramp-down constraints. If "DeltaRampDown" has length 1 then DM[ t ]
 *   contains the same value for all t. Otherwise, DeltaRampDown[ i ] is the
 *   fixed value of DM[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "StorageLevelRho", of type double and to be either of size
 *   1 or indexed over the dimension "NumberIntervals". This is meant to
 *   represent the vector SLR[ t ] that, for each time instant t, contains the
 *   inefficiency of storing energy level of the unit for the corresponding
 *   time step. This variable is optional; if it is not provided then it is
 *   assumed that this unit may not be capable of having any storing energy
 *   levels, which correspond to SLR[ t ] == 0 for all t. If "StorageLevelRho"
 *   has length 1 then SLR[ t ] contains the same value for all t. Otherwise,
 *   StorageLevelRho[ i ] is the fixed value of SLR[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ] with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The scalar variable "GasPrice", of type double and not indexed over any
 *   dimension. This variable indicates the gas price which is fixed external
 *   input parameter.
 * */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the PowerToGasUnitBlock
/** The PowerToGasUnitBlock class use get_variable() method to access to
 *  each "group" of active power variable that may create in UnitBlock class.
 *  Moreover, PowerToGasUnitBlock is defined one more group of variable as
 *  follow:
 *
 * - the storage level variable;
 *
 *  This group of variable may have size f_time_horizon or empty size, and
 *  it's optional, and it is also possible to restrict it of is generated with
 *  the parameter stvv. If stvv is not nullptr and it is a
 *  SimpleConfiguration<int>, or if
 *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 *  SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 *  of the optional variables should be created. If the Configuration is not
 *  available, the default value is taken to be 0.*/
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the PowerToGasUnitBlock
/** Method that generates the static constraint of the PowerToGasUnitBlock.
 * The power-to-gas units are optimized as a sub-model to the EUC with an
 * aggregated flexibility. The objective is to maximize the “profit” generated
 * by producing and selling gas given the time horizon which indicated with
 * \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$. The main
 * power-to-gas unit constraints are define as:
 *
* - the maximum and minimum level of power output constraints. This is a
 *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
 *   the entry t = 0, ...,f_time_horizon - 1 being the maximum and minimum
 *   power output bounds at time t. these ensure the maximum or minimum amount
 *   of energy that unit can produce(or use) when it is on(or off).
 *
 *   \f[
 *    p^{ac}_{t} \in [ P^{mn}_{t} , P^{mx}_{t}]
 *                              \quad t \in \mathcal{T}          \quad (1)
 *   \f]
 *   where \f$ P^{mx}_{t} \f$ and \f$ P^{mn}_{t} \f$ are the maximum and
 *   minimum power output parameters for each time t of the time horizon
 *   \f$ \mathcal{T} \f$ respectively.
 *
 * - the ramp-up and ramp-down constraints. Each of them are a
 *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
 *   the entry t = 0, ...,f_time_horizon - 1 being the ramp up and ramp down
 *   constraints which are presented as:
 *
 *   \f[
 *    p^{ac}_{t} - p^{ac}_{t-1} \leq \Delta^{up}_{t}
 *         \quad t \in \mathcal{T}                               \quad (2)
 *   \f]
 *
 *   \f[
 *    p^{ac}_{t} - p^{ac}_{t-1} \geq - \Delta^{dn}_{t}
 *         \quad t \in \mathcal{T}                               \quad (3)
 *   \f]
 *   where \f$ \Delta^{up}_{t} \f$ and \f$ \Delta^{dn}_{t} \f$ are the delta
 *   ramp-up and delta ramp down threshold for each time t of the time horizon
 *   \f$ \mathcal{T} \f$ respectively.
 *
 * - active power relation with storing energy levels constraints.  That is a
 *   std::vector<FRowConstraint>; with the dimension of f_time_horizon,
 *   where the entry t = 0,...,f_time_horizon - 1 being the storage level
 *   relation with active power at time t.
 *
 *   \f[
 *    v^{ptg}_{t} = v^{ptg}_{t-1} + \rho^{ptg}_{t}p^{ac}
 *                  \quad t \in \mathcal{T}                      \quad (4)
 *   \f]
 * - the storage level bounds constraints. it gives the storage levels upper
 *   bound and lower bound at each time instant t.
 *
 *   \f[
 *    v^{ptg}_{t} \in [ V^{ptg,mn}_{t} , V^{ptg,mx}_{t}]
 *                              \quad t \in \mathcal{T}          \quad (5)
 *   \f]
 *
 *
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the objective of the PowerToGasUnitBlock
/** Method that generates the objective of the PowerToGasUnitBlock.
 * - Objective function: the objective function of the PowerToGasUnitBlock
 *   is given as follow:
 *
 *   \f[
 *     \max ( \sum_{ t \in  [0 , \mathcal{T}]  }
 *     ( C^{gas} p^{ac}_t )
 *   \f]
 *
 *   where \f$ C^{gas} \f$, is the gas price.
 *
*/
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE PowerToGasUnitBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the PowerToGasUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of PowerToGasUnitBlock units
 * @{ */

 /// Returns the gas price value
 double get_gas_price() const { return f_gas_price; }

 /// Returns the initial power value
 double get_initial_power() const { return f_initial_power; }
/*--------------------------------------------------------------------------*/
/// returns the vector of minimum storage
/** The method returned a std::vector< double > V and each element of V
 * contains the minimum storage at time t. There are three possible cases:
 *
 * - if the vector is empty, then the minimum storage of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the minimum storage of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the minimum storage value at time t. */

 const std::vector< double > & get_minimum_storage() const {
  return( v_minimum_storage );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of maximum storage
/** The method returned a std::vector< double > V and each element of V
 * contains the maximum storage at time t. There are three possible cases:
 *
 * - if the vector is empty, then the maximum storage of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the maximum storage of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the maximum storage value at time t. */

 const std::vector< double > & get_maximum_storage() const {
  return( v_maximum_storage );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of minimum power
/** The method returned a std::vector< double > V and each element of V
 * contains the minimum power at time t. There are three possible cases:
 *
 * - if the vector is empty, then the minimum power of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the minimum power of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the minimum power value at time t. */

 const std::vector< double > & get_minimum_power() const {
  return( v_minimum_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of maximum power
/** The method returned a std::vector< double > V and each element of V
 * contains the maximum power at time t. There are three possible cases:
 *
 * - if the vector is empty, then the maximum power of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the maximum power of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the maximum power value at time t. */

 const std::vector< double > & get_maximum_power() const {
  return( v_maximum_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of delta ramp up
/** The method returned a std::vector< double > V and each element of V
 * contains the delta ramp up at time t. There are three possible cases:
 *
 * - if the vector is empty, then the delta ramp up of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the delta ramp up of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the delta ramp up value at time t. */

 const std::vector< double > & get_delta_ramp_up() const {
  return( v_delta_ramp_up );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of delta ramp down
/** The method returned a std::vector< double > V and each element of V
 * contains the delta ramp down at time t. There are three possible cases:
 *
 * - if the vector is empty, then the delta ramp down of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the delta ramp down of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the delta ramp down value at time t. */

 const std::vector< double > & get_delta_ramp_down() const {
  return( v_delta_ramp_down );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of storage level rho
/** The method returned a std::vector< double > V and each element of V
 * contains the storage level rho at time t. There are three possible cases:
 *
 * - if the vector is empty, then the storage level rho of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the storage level rho
 *   of the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the storage level rho value at time t. */

 const std::vector< double > & get_storage_level_rho() const {
  return( v_storage_level_rho );
 }

/**@} ----------------------------------------------------------------------*/
/*------ METHODS FOR READING THE Variable OF THE PowerToGasUnitBlock -------*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the PowerToGasUnitBlock
 *
 * These methods allow to read the two groups of Variable that any
 * PowerToGasUnitBlock in principle has (although some may not):
 *
 * - the storage level variables
 *
 * This group of variable is (if not empty) std::vector< ColVariable > with
 * the dimension time horizon.
 * @{ */

/// returns the vector of storage level variables
/** The returned std::vector< ColVariable >, say V, contains the
 * storage level variables and is indexed over the dimension time horizon.
 * There are two possible cases:
 *
 * - if V is empty(), then these variables are not defined;
 *
 * - otherwise, V must have size of get_time_horizon() and V[ t ] is the
 *   storage level variable for time step t.*/

 const std::vector< ColVariable > & get_storage_level() const {
  return v_storage_level;
 }

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE PowerToGasUnitBlock----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the PowerToGasUnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * PowerToGasUnitBlock. See
 * PowerToGasUnitBlock::deserialize( netCDF::NcGroup ) for details of the
 * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR INITIALIZING THE PowerToGasUnitBlock -----------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the PowerToGasUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "PowerToGasUnitBlock::load() not "
                            "implemented yet") );
 };

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/
 /// The vector of minimum storage
 std::vector< double > v_minimum_storage;

 /// The vector of maximum storage
 std::vector< double > v_maximum_storage;

 /// the vector of MinPower
 std::vector< double >  v_minimum_power;

 /// the vector of MaxPower
 std::vector< double >  v_maximum_power;

 /// the vector of RampUp
 std::vector< double >  v_delta_ramp_up;

 /// the vector of RampDown
 std::vector< double >  v_delta_ramp_down;

 /// the vector of StorageLevelRho
 std::vector< double >  v_storage_level_rho;

 /// the gas price value
 double f_gas_price;

 /// the InitialPower value
 double f_initial_power;

/*-----------------------------variables------------------------------------*/
 /// the vector of storage level variables
 std::vector< ColVariable > v_storage_level;

/*----------------------------constraints-----------------------------------*/
/// the active power bound constraints
 std::vector< FRowConstraint > active_power_bound_Constraints;

/// the ramp up constraints
 std::vector< FRowConstraint > ramp_up_Constraints;

/// the ramp down constraints
 std::vector< FRowConstraint > ramp_down_Constraints;

/// the storage and active power relation
 std::vector< FRowConstraint > storage_active_power_Constraints;

 /// the storage level bound constraints
 std::vector< FRowConstraint > storage_level_bound_Constraints;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

};  // end( class( PowerToGasUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* PowerToGasUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File PowerToGasUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
