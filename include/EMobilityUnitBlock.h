/*--------------------------------------------------------------------------*/
/*---------------------- File EMobilityUnitBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class EMobilityUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define a "reasonably standard"
 * E-mobility unit at Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 02 - 08 - 2019
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

#ifndef __EMobilityUnitBlock
#define __EMobilityUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/
#include "ColVariable.h"
#include "FRowConstraint.h"
#include "OneVarConstraint.h"
#include "FRealObjective.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// Namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS EMobilityUnitBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the E-mobility unit problem
/** The EMobilityUnitBlock class implements the Block concept [see Block.h]
 * for a "reasonably standard" E-mobility unit of a Unit Commitment Problem.
 * That is, the class is designed in order to give mathematical formulation to
 * describe the operation of large set of E-mobility. Since the transport
 * sector is moving towards electrification, electric mobility will have a
 * rising impact on the electricity system. First, electricity demand is
 * growing due to a higher amount of electric vehicles that need to be
 * charged. On the other hand, vehicles are used only a small amount of time
 * while being charged over a much longer timespan (e.g. at night). This
 * allows to shift the charging process in time and provide this flexibility
 * to the overall energy system by means of an additional generator
 * (vehicle-to-grid) or an additional load (power-to-vehicle). Within plan4res
 * electric vehicles can be considered as:
 * - A static electric demand that does not provide flexibility to the overall
 *   system;
 * - A flexible electric demand that provides flexibility to the overall
 *   system while accounting for storage level constraints.
 *
 * The first case described above (static electric demand) is trivial and has
 * no need for a mathematical description, since the electric demand needed by
 * electric vehicles is added to the “normal” electric demand. In the second
 * case the E-mobility model is a sub-model of the European unit commitment
 * model and acts similar to a BatteryStorageUnitBlock. The technical and
 * physical constraints are mainly divided in several different categories as:
 * - the active power limit constraints;
 * - the ramping constraints;
 * - the demand constraints;
 * - the storage level bounds constraints. */
class EMobilityUnitBlock : public UnitBlock {

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
/** Constructor of EMobilityUnitBlock, taking possibly a pointer of its
 * father Block.
 */

 explicit EMobilityUnitBlock( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/

/// destructor of EMobilityUnitBlock, it is empty

 ~EMobilityUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the EMobilityUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group must contain all the data required by the base UnitBlock,
 * as described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
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
 * - The variable "DeltaRampUp", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
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
 * - The variable "EMobilityRho", of type double and to be either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector EMR[ t ] that, for each time instant t, contains the possible
 *   fraction of storage level that can be used as charging/discharging of
 *   the unit for the corresponding time step. This variable is optional; if
 *   it is not provided then it is assumed that this unit may not be capable
 *   of charging/discharging, which correspond to EMR[ t ] == 0 for all t. If
 *   "EMobilityRho" has length 1 then EMR[ t ] contains the same value for all
 *   t. Otherwise, EMobilityRho[ i ] is the fixed value of EMR[ t ] for all t
 *   in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ] with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1
 *   or "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "Demand", of type double and indexed over the dimension
 *   "TimeHorizon": entry Demand[ t ] is assumed to contain the energy needed
 *   to discharge of a battery in the time t.
 *
 * - The scalar variable "InitialStorage", of type double and not indexed over
 *   any dimension. This variable indicates the amount of the storage level
 *   that the unit was producing at time instant -1, i.e., before the start of
 *   the time horizon; this is necessary to compute the storage level
 *   connection with active power of the aggregated vehicles and fulfill the
 *   driving profiles constraints.
 *
 * */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the EMobilityUnitBlock
/** The EMobilityUnitBlock class use get_variable() method to access to
 *  the "group" of  active power variable that may create in UnitBlock class.
 *  Moreover, EMobilityUnitBlock is defined one more group of variable as
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
 *  available, the default value is taken to be 0. */
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the EMobilityUnitBlock
/** Method that generates the static constraint of the EMobilityUnitBlock. The
 * operations of the E-mobility units are described on a discrete time horizon
 * as dictated by the UnitBlock interface. In this description we indicate it
 * with \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$. The main
 * E-mobility unit constraints are define as:
 * - the charging/discharging limit constraints;
 *
 *   \f[
 *    p^{ac}_{t} \in [ P^{mn}_{t} , P^{mx}_{t}]
 *                              \quad t \in \mathcal{T}          \quad (1)
 *   \f]
 *
 * - the rammping constraints;
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
 *
 * - the fulfill the driving profiles constraints;
 *
 *   \f[
 *    v^{emob}_{t} = v^{emob}_{t-1} - \rho^{emob}_{t}p^{ac}+_{t} -
 *    d^{emob,dch}_t            \quad t \in \mathcal{T}          \quad (4)
 *   \f]
 * - the storage level bounds constraints.
 *
 *   \f[
 *    v^{emob}_{t} \in [ V^{emob,mn}_{t} , V^{emob,mx}_{t}]
 *                              \quad t \in \mathcal{T}          \quad (5)
 *   \f]
 *
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the objective of the EMobilityUnitBlock
/** Method that generates the objective of the EMobilityUnitBlock.
 * //TODO I should put the objective function here
 *
*/
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE EMobilityUnitBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the EMobilityUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of E-mobility units
 * @{  * */

 /// Returns the initial storage level value
 double get_initial_storage_level() const { return f_initial_storage_level; }
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
/// returns the vector of e-mobility rho
/** The method returned a std::vector< double > V and each element of V
 * contains the e-mobility rho at time t. There are three possible cases:
 *
 * - if the vector is empty, then the e-mobility rho of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the e-mobility rho of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the e-mobility rho value at time t. */

 const std::vector< double > & get_emobility_rho() const {
  return( v_emobility_rho);
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of demand
/** The method returned a std::vector< double > V and each element of V
 * contains the demand at time t. There are two possible cases:
 *
 * - if the vector is empty, then the demand of the unit is 0;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the demand value at time t. */

 const std::vector< double > & get_demand() const {
  return( v_demand);
 }
/**@} ----------------------------------------------------------------------*/
/*------- METHODS FOR READING THE Variable OF THE EMobilityUnitBlock -------*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the EMobilityUnitBlock
 *
 * These methods allow to read the a group of Variable that any
 * EMobilityUnitBlock in principle has (although may not):
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
/*---------------- METHODS FOR SAVING THE EMobilityUnitBlock----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the EMobilityUnitBlock
 *  @{ */
/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * EMobilityUnitBlock. See EMobilityUnitBlock::deserialize( netCDF::NcGroup )
 * for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR INITIALIZING THE EMobilityUnitBlock -------------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the EMobilityUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "EMobilityUnitBlock::load() not "
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

 /// the vector of e-mobility rho
 std::vector< double >  v_emobility_rho;

 /// the vector of demand
 std::vector< double >  v_demand;

 /// the initial storage level value
 double f_initial_storage_level;

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

/// the demand constraints
 std::vector< FRowConstraint > demand_Constraints;

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

};  // end( class( EMobilityUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* EMobilityUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File EMobilityUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
