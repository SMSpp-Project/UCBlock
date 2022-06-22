/*--------------------------------------------------------------------------*/
/*------------------------- File BatteryUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BatteryUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" Battery
 * storage, E-mobility, Centralized demand response, Distributed load
 * management, Distributed storage, and Power to gas units in a single class
 * at Unit commitment problem.
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

#ifndef __BatteryUnitBlock
 #define __BatteryUnitBlock
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

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS BatteryUnitBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the BatteryUnit problem
/** The BatteryUnitBlock class implements the Block concept [see Block.h] for
 * a large class of units that allow direct storage of electrical energy. This
 * can be the case of actual physical batteries, either "large" (battery
 * storage) or "small" (e-mobility, distributed storage), of methods that use
 * some intermediate energy vector with limited local storage/production
 * (power-to-gas units), as well as of "logical" mechanisms that allow to
 * temporally shift production/consumption in a limited way, thereby acting
 * like an energy storage (centralized demand response, distributed load
 * management). BatteryUnitBlock provides a quite general concept of battery
 * that covers different units which mostly fit the same mathematical
 * equations pattern. For instance, a BatteryUnitBlock may or may not have a
 * fixed demand (e-mobility has, other units have not) and it may or may not
 * provide primary and secondary reserve (battery storage may do, but other
 * units don't).
 *
 * Battery storage provide an additional flexibility to the system by shifting
 * a surplus of electric energy (e.g., due to high renewable feeding) to times
 * with high demand or lower renewable generation. The distributed battery
 * storage can be aggregated in the energy cells or directly placed in a
 * single node of the network. We will therefore not stress this dependency in
 * the subsequent equations. We emphasize that potential contribution of
 * batteries to inertia is still a subject of active research and should be
 * considered as optional. Besides, since the transport sector is moving
 * towards electrification, electric mobility will have a rising impact on the
 * electricity system. First, electricity demand is growing due to a higher
 * amount of electric vehicles that need to be charged. On the other hand,
 * vehicles are used only a small amount of time while being charged over a
 * much longer timespan (e.g., at night). This allows to shift the charging
 * process in time and provide this flexibility to the overall energy system
 * by means of an additional generator (vehicle-to-grid) or an additional load
 * (power-to-vehicle). Two main differences between battery storages unit and
 * other existing units in this class are:
 *
 * - battery storages unit can do primary and secondary reserve, while
 *   other units cannot;
 *
 * - some of the units may have a fixed demand that battery storages unit has
 *   not.
 *
 * Moreover, as the considered storage cycle is small w.r.t. the EUC time
 * horizon, distributed storage is not considered as seasonal storage. Hence,
 * the associated mathematical description follows the same equations as the
 * one provided for battery storages unit. The specificity of distributed
 * storage only relies on the fact that it is connected to a distribution grid
 * node.
 *
 * To model the BatteryUnitBlock systems several technical parameters have to
 * be considered. These are divided into the battery storage level parameters,
 * the ramping parameters, the active power bound parameters, and a flexible
 * electric demand that provides flexibility to the overall system while
 * accounting for storage level constraints. The technical and physical
 * constraints are mainly divided in several different categories as:
 *
 * - the maximum and minimum power output constraints according to primary and
 *   secondary spinning reserves (if any);
 *
 * - the ramp-up and ramp-down constraints;
 *
 * - the active power relation with storing and extracting energy levels
 *   constraints (if any);
 *
 * - the intake upper bound (if any);
 *
 * - the storage level constraints;
 *
 * - the binary variable relation with intake and outtake level constraints
 *   (if any);
 *
 * - the primary reserve upper bound (if any);
 *
 * - the secondary reserve upper bound (if any);
 *
 * - the demand constraints (if any). */

class BatteryUnitBlock : public UnitBlock
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 enum battery_type
 {
  ///< when ExtractingBatteryRho >= 1 and StoringBatteryRho <=1
  ASSUME_POSITIVE_PRICES ,
  ///< when ExtractingBatteryRho and StoringBatteryRho not defined(both == 1)
  NO_Binary_Variables_Constraints ,
  ///< otherwise binary variables with related constraints are needed
  Binary_Variables_Constraints
 };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father and the time horizon
 /** Constructor of BatteryUnitBlock, taking possibly a pointer of its father
  * Block. */

 explicit BatteryUnitBlock( Block * f_block = nullptr , Index t = 0 )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of BatteryUnitBlock

 virtual ~BatteryUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the BatteryUnitBlock. Besides the mandatory "type" attribute of any
  * :Block, the group must contain all the data required by the base
  * UnitBlock, as described in the comments to UnitBlock::deserialize(
  * netCDF::NcGroup ). In particular, we refer to that description for the
  * crucial dimensions "TimeHorizon", "NumberIntervals" and
  * "ChangeIntervals". The netCDF::NcGroup must then also contain:
  *
  * - The variable "MinStorage", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MinS[ t ] that,
  *   for each time instant t, contains the minimum storage level of the unit
  *   for the corresponding time step. If "MinStorage" has length 1 then MinS[
  *   t ] contains the same value for all t. Otherwise, MinStorage[ i ] is the
  *   fixed value of MinS[ t ] for all t in the interval [ ChangeIntervals[ i
  *   - 1 ] , ChangeIntervals[ i ] ], with the assumption that
  *   ChangeIntervals[ - 1 ] = 0. Note that it must always be 0 <= MinS[ t ] <
  *   MaxS[ t ] for all t. If NumberIntervals <= 1 or NumberIntervals >=
  *   TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "MaxStorage", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MaxS[ t ] that,
  *   for each time instant t, contains the maximum storage level of the unit
  *   for the corresponding time step. If "MaxStorage" has length 1 then MaxS[
  *   t ] contains the same value for all t. Otherwise, MaxStorage[ i ] is the
  *   fixed value of MaxS[ t ] for all t in the interval [ ChangeIntervals[ i
  *   - 1 ] , ChangeIntervals[ i ] ], with the assumption that
  *   ChangeIntervals[ - 1 ] = 0. Note that it must always be [0 <=] MinS[ t ]
  *   < MaxS[ t ] for all t. If NumberIntervals <= 1 or NumberIntervals >=
  *   TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "MinPower", of type double and either of size 1 or indexed
  *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
  *   provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MinP[ t ] that,
  *   for each time instant t, contains the minimum active power output value
  *   of the unit for the corresponding time step.  If "MinPower" has length 1
  *   then MinP[ t ] contains the same value for all t. Otherwise, MinPower[ i
  *   ] is the fixed value of MinP[ t ] for all t in the interval [
  *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. Note that it must be MinP[ t ] <= 0 for
  *   all t. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
  *   the mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded.
  *
  * - The variable "MaxPower", of type double and either of size 1 or indexed
  *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
  *   provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MaxP[ t ] that,
  *   for each time instant t, contains the maximum active power output value
  *   of the unit for the corresponding time step.  If "MaxPower" has length 1
  *   then MaxP[ t ] contains the same value for all t. Otherwise, MaxPower[ i
  *   ] is the fixed value of MaxP[ t ] for all t in the interval [
  *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. Note that it must be MinP[ t ] < MaxP[
  *   t ] for all t. If NumberIntervals <= 1 or NumberIntervals >=
  *   TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The scalar variable "InitialPower", of type double and not indexed over
  *   any dimension. This variable indicates the amount of the power that the
  *   unit was producing at time instant -1, i.e., before the start of the
  *   time horizon; this is necessary to compute the ramp-up and ramp-down
  *   constraints. This variable is optional; if "DeltaRampUp" and
  *   "DeltaRampDown" are not present, "InitialPower" should not be read,
  *   since there are no ramping constraints. If "DeltaRampUp" and
  *   "DeltaRampDown" are present but "InitialPower" is not provided, its
  *   value is taken to be 0.
  *
  * - The variable "MaxPrimaryPower", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MaxPP[ t ] that,
  *   for each time instant t, contains the maximum active power that can be
  *   used as primary reserve of the unit for the corresponding time step. If
  *   "MaxPrimaryPower" has length 1 then MaxPP[ t ] contains the same value
  *   for all t. Otherwise, MaxPrimaryPower[ i ] is the fixed value of MaxPP[
  *   t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
  *   = 0. This variable is optional, if is not provided then MaxPP[ t ] == 0
  *   for all t. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon,
  *   then the mapping clearly does not require "ChangeIntervals", which in
  *   fact is not loaded.
  *
  * - The variable "MaxSecondaryPower", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
  *   provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MaxSP[ t ] that,
  *   for each time instant t, contains the maximum active power that can be
  *   used as secondary reserve of the unit for the corresponding time step. If
  *   "MaxSecondaryPower" has length 1 then MaxSP[ t ] contains the same value
  *   for all t. Otherwise, MaxSecondaryPower[ i ] is the fixed value of MaxSP[
  *   t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
  *   0. This variable is optional, if is not provided then MaxSP[ t ] == 0 for
  *   all t. Note that MaxPP[ t ] == 0 implies MaxSP[ t ] == 0 (that is, if
  *   MaxPrimaryPower is not defined then neither should MaxSecondaryPower). If
  *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
  *   clearly does not require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "DeltaRampUp", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector DP[ t ] that, for
  *   each time instant t, contains the ramp-up value of the unit for the
  *   corresponding time step, i.e., the maximum possible increase of active
  *   power production w.r.t. the power that had been produced in time instant
  *   t - 1, if any. If "DeltaRampUp" has length 1 then DP[ t ] contains the
  *   same value for all t. Otherwise, DeltaRampUp[ i ] is the fixed value of
  *   DP[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
  *   = 0. This variable is optional; if it is not provided then it is assumed
  *   that DP[ t ] == MaxP[ t ], i.e., the unit can ramp up by an arbitrary
  *   amount, i.e., there are no ramp-up constraints. If NumberIntervals <= 1
  *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "DeltaRampDown", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector DM[ t ] that, for
  *   each time instant t, contains the ramp-down value of the unit for the
  *   corresponding time step, i.e., the maximum possible decrease of active
  *   power production w.r.t. the power that had been produced in time instant
  *   t - 1, if any. If "DeltaRampDown" has length 1 then DM[ t ] contains the
  *   same value for all t. Otherwise, DeltaRampDown[ i ] is the fixed value
  *   of DM[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
  *   = 0. This variable is optional; if it is not provided then it is assumed
  *   that DM[ t ] == MaxP[ t ], i.e., the unit can ramp down an arbitrary
  *   amount, i.e., there are no ramp-down constraints. If NumberIntervals <=
  *   1 or NumberIntervals >= TimeHorizon, then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "StoringBatteryRho", of type double and to be either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector SBR[
  *   t ] that, for each time instant t, contains the inefficiency of storing
  *   energy of the unit for the corresponding time step. This variable is
  *   optional; if it is not provided then it is assumed that SBR[ t ] == 1
  *   for all t, i.e., no (significant) energy is spent just for storing it in
  *   the battery (this simplifies the model somewhat, see below). If
  *   "StoringBatteryRho" has length 1 then SBR[ t ] contains the same value
  *   for all t. Otherwise, StoringBatteryRho[ i ] is the fixed value of SBR[
  *   t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ] with the assumption that ChangeIntervals[ - 1 ] =
  *   0. Note that it must be always such that SBR[ t ] <= 1 for all t (as
  *   SBR[ t ] is the amount of energy actually going in the battery for each
  *   1 unit of input energy). If NumberIntervals <= 1 or NumberIntervals >=
  *   TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "ExtractingBatterRho", of type double and to be either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector EBR[
  *   t ] that, for each time instant t, contains the inefficiency of
  *   extracting energy of the unit for the corresponding time step. This
  *   variable is optional; if it is not provided, then it is assumed that
  *   EBR[ t ] == 1 for all t, i.e., no (significant) energy is spent just for
  *   extracting it from the battery (this simplifies the model somewhat, see
  *   below). If "ExtractingBatterRho" has length 1 then EBR[ t ] contains the
  *   same value for all t. Otherwise, ExtractingBatterRho[ i ] is the fixed
  *   value of EBR[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ] with the assumption that ChangeIntervals[ - 1 ] =
  *   0. Note that it must be always such that EBR[ t ] >= 1 [>= SBR[ t ]] for
  *   all t (as EBR[ t ] is the amount of energy that is taken away from the
  *   battery to obtain 1 unit of output energy). If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded.
  *
  * Note: the special case in which EBR[ t ] == SBR[ t ] == 1 for all t, i.e.,
  * no energy is spent for storing it in / retrieving it from the battery,
  * leads to significantly simpler mathematical models. In particular, one
  * single variable can be used to represent both storing and retrieving,
  * rather than requiring two separate ones (unless primary and secondary
  * reserve are allowed and/or the cost is defined, since this also requires
  * using two), and the binary variables need not to be defined. For details,
  * see the comments to generate_abstract_variables() and
  * generate_abstract_constraints().
  *
  * - The scalar variable "InitialStorage", of type double and not indexed
  *   over any dimension. This variable indicates the amount of the storage
  *   level that the unit was producing at time instant -1, i.e., before the
  *   start of the time horizon; this is necessary to compute the storage
  *   level connection with intake and outtake constraints.
  *
  * - The variable "Cost", of type double and either of size 1 or indexed over
  *   the dimension "NumberIntervals" (if "NumberIntervals" is not provided,
  *   then this variable can also be indexed over "TimeHorizon"). This is
  *   meant to represent the vector C[ t ] that, for each time instant t,
  *   contains the monetary cost of storing one unit or energy into, or
  *   extracting it from, the battery (the cost is the same in both cases) at
  *   the corresponding time step.  This variable is optional; if it is not
  *   provided then it's taken to be zero. If "Cost" has length 1 then C[ t ]
  *   contains the same value for all t. Otherwise, Cost[ i ] is the fixed
  *   value of C[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
  *   = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "Demand", of type double and indexed over the dimension
  *   "TimeHorizon": the entry Demand[ t ] is assumed to contain the amount of
  *   energy that must be discharged from the battery and "sent away for some
  *   other purpose" (say, driving your e-car) at time t. This variable is
  *   optional; if it isn't defined, then Demand[ t ] == 0. Otherwise the
  *   Demand[ t ] contains the demand value for each time instant t.
  *
  * - The scalar variable "Kappa", of type netCDF::NcDouble(). This variable
  *   contains the factor that multiplies the minimum and maximum active
  *   power, maximum primary and secondary reserve, and the minimum and
  *   maximum storage levels, at each time instant t. This variable is
  *   optional, if it is not provided it is taken to be Kappa == 1. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the BatteryUnitBlock
 /** The BatteryUnitBlock class use get_variable() method to access to each
  *  "group" of variable that may create in UnitBlock class which are:
  *
  *  - the primary spinning reserve variables;
  *
  *  - the secondary spinning reserve variables;
  *
  *  - the active power variables; it can be positive or negative, if it is
  *    positive the unit is giving energy to the system, if it is negative it
  *    is taking energy away and adding to the storage. Since the storing and
  *    extracting amount of active power are not always equal, to deal with
  *    this issue, the usual trick of splitting the active power variable by
  *    two new non-negative variables which are called intake and outtake
  *    levels for each time t (see equation (5)) is used. If
  *    "StoringBatteryRho" == "ExtractingBatterRho" == 1, we don not need to
  *    split the active power and the constraint (5-7 and 10-11) will be
  *    replaced by (8)).
  *
  *  All of those variables are optional except the active power variables in
  *  the sense that the model may just not have them and whenever a group of
  *  above variables is created, its size will be the time horizon. Moreover,
  *  BatteryUnitBlock defines four more groups of variables as follows:
  *
  *  - the storage level variables;
  *
  *  - the intake and outtake levels variable; they are needed to split the
  *    active power variable (if it's needed);
  *
  *  - the binary variables; when "StoringBatteryRho" == "ExtractingBatterRho"
  *    == 1, then this binary variable and all constraints which are depended
  *    on this variable is not required to be define.
  *
  *  These three groups of variables may have size f_time_horizon or empty
  *  size.  All of these variables are optional, and it is also possible to
  *  restrict which of the subsets are generated with the parameter stvv. If
  *  stvv is not nullptr and it is a SimpleConfiguration<int>, or if
  *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is
  *  a SimpleConfiguration<int>, then the f_value (an int) indicates whether
  *  each of the optional variables should be created. If the Configuration is
  *  not available, the default value is taken to be 0. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the BatteryUnitBlock
 /** Method that generates the static constraint of the BatteryUnitBlock. The
  * operations of the battery storage unit are described on a discrete time
  * horizon as dictated by the UnitBlock interface. In this description we
  * indicate it with \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\}
  * \f$. The main constraints of this unit are define as:
  *
  * - maximum and minimum power output constraints according to primary and
  *   secondary spinning reserves are presented in (1)-(2). Each of them is a
  *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
  *   the entry t = 0, ...,f_time_horizon - 1 being the maximum and minimum
  *   power output value according to the primary and the secondary spinning
  *   reserves at time t. these ensure the maximum(or minimum) amount of
  *   energy that unit can produce(or use) when it is on(or off).
  *
  *   \f[
  *      p^{ac}_{t} + p^{pr}_{t} + p^{sc}_{t} \leq P^{mx}_{t}
  *          \quad t \in \mathcal{T}                              \quad (1)
  *   \f]
  *
  *   \f[
  *     P^{mn}_{t} \leq p^{ac}_{t} - p^{pr}_{t} - p^{sc}_{t}
  *         \quad t \in \mathcal{T}                               \quad (2)
  *   \f]
  *
  *   where \f$ P^{mx}_{t} \f$ and \f$ P^{mn}_{t} \f$ are the maximum and
  *   minimum power output parameters for each time t of the time horizon \f$
  *   \mathcal{T} \f$ respectively.
  *
  * - ramp-up and ramp-down constraints are presented in (3)-(4). Each of them
  *   is a std::vector<FRowConstraint>; with the dimension of f_time_horizon,
  *   where the entry t = 0, ...,f_time_horizon - 1 being the ramp up and ramp
  *   down constraints which are presented as:
  *
  *   \f[
  *    p^{ac}_{t} - p^{ac}_{t-1} \leq \Delta^{up}_{t}
  *         \quad t \in \mathcal{T}                               \quad (3)
  *   \f]
  *
  *   \f[
  *    p^{ac}_{t} - p^{ac}_{t-1} \geq - \Delta^{dn}_{t}
  *         \quad t \in \mathcal{T}                               \quad (4)
  *   \f]
  *
  *   where \f$ \Delta^{up}_{t} \f$ and \f$ \Delta^{dn}_{t} \f$ are the delta
  *   ramp-up and delta ramp down threshold for each time t of the time
  *   horizon \f$ \mathcal{T} \f$ respectively.
  *
  * - active power relation with intake and outtake levels constraints are
  *   presented in (5). Each of them is a std::vector<FRowConstraint>; with
  *   the dimension of f_time_horizon, where the entry t =
  *   0,...,f_time_horizon - 1 being the active power relation with intake and
  *   outtake levels at time t.  These ensure the active power at each time
  *   should be equal to the intake and outtake difference. The equation (6)
  *   also indicates the upper bound of intake level at each time instant t.
  *
  *   \f[
  *    p^{ac}_{t} = p^+_t - p^-_{t}
  *         \quad t \in \mathcal{T}                               \quad (5)
  *   \f]
  *
  *   \f[
  *     p^+_t \leq  P^{mx}_{t}
  *         \quad t \in \mathcal{T}                               \quad (6)
  *   \f]
  *
  * - storage level relation with intake and outtake levels (if any)
  *   constraints in Battery unit are presented in (7). That is a
  *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
  *   the entry t = 0,...,f_time_horizon - 1 being the storage level relation
  *   with intake and outtake levels at time t.
  *
  *   \f[
  *    v^{ba}_{t} = v^{ba}_{t-1} - \rho^+_{t}p^+_{t} +
  *    \rho^-_{t}p^-_{t}-d^{ba}_t      \quad t \in \mathcal{T}    \quad (7)
  *   \f]
  *
  *   Note that if the equation (7) will change as below which is a
  *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
  *   the entry t = 0,...,f_time_horizon - 1 being the storage level relation
  *   with battery demand (if any) at time t.
  *
  *   \f[
  *    v^{ba}_{t} = v^{ba}_{t-1} - p^{ac}_{t} - d^{ba}_t
  *               \quad t \in \mathcal{T}          \quad (8)
  *   \f]
  *
  *   The equation (9) gives the storage levels upper bound and lower bound at
  *   each time instant t.
  *
  *   \f[
  *    v^{ba}_{t} \in [ V^{mn}_{t} , V^{mx}_{t}]
  *                              \quad t \in \mathcal{T}          \quad (9)
  *   \f]
  *
  *   where \f$ \rho^+_{t} \f$ and \f$ \rho^-_{t} \f$ are the
  *   ExtractingBatteryRho and StoringBatteryRho, and \f$ V^{mn}_t\f$ and \f$
  *   V^{mx}_t\f$ are the minimum and maximum storage level for each time t of
  *   the time horizon \f$ \mathcal{T} \f$ respectively.
  *
  * - binary variable relation with storing and extracting energy level (if
  *   any) constraints are presented in (10-11). Each of them is a
  *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
  *   the entry t = 0,...,f_time_horizon - 1 being the binary variable
  *   relation with storing and extracting energy levels at time t.
  *
  *   \f[
  *    p^+_{t} \leq u^+_t P^{mx}_{t}
  *                              \quad t \in \mathcal{T}          \quad (10)
  *   \f]
  *
  *   \f[
  *    p^-_{t} \leq -(1 - u^+_t) P^{mn}_{t}
  *                              \quad t \in \mathcal{T}         \quad (11)
  *   \f]
  *
  *   Note that when \f$ \rho^+_t = \rho^-_t = 1 \f$, the binary variable \f$
  *   u^+_t \f$ is not required and neither are the last two constraints
  *   (10-11).
  *
  * - primary and secondary reserve upper bounds (if any) are presented by the
  *   equations (12-13). Each of them is a std::vector<FRowConstraint>; with
  *   the dimension of f_time_horizon, where the entry
  *   t = 0,...,f_time_horizon - 1 being the primary and secondary reserve
  *   upper bounds at time t.
  *
  *   \f[
  *     p^{pr}_{t} \leq P^{mx, pr}_{t}
  *                              \quad t \in \mathcal{T}          \quad (12)
  *   \f]
  *
  *   \f[
  *     p^{sc}_{t} \leq P^{mx, sc}_{t}
  *                              \quad t \in \mathcal{T}          \quad (13)
  *   \f]
  */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the BatteryUnitBlock
 /** Method that generates the objective of the BatteryUnitBlock.
  * - Objective function: the objective function of the BatteryUnitBlock
  *   is given as follow:
  *
  *   \f[
  *     \min ( \sum_{ t \in  [0 , \mathcal{T}]  }
  *     ( C_t (p^+_t +  p^-_t) )
  *   \f]
  *
  *   where \f$ C_t \f$, is a certain proportion cost function. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE BatteryUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the BatteryUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of battery storage units
 * @{ */

 /// returns the initial storage value
 double get_initial_storage( void ) const { return( f_initial_storage ); }

 /// returns the initial power value
 double get_initial_power( void ) const { return( f_initial_power ); }

 /// returns the operation and maintenance costs
 double get_oem_cost( void ) const { return( f_oem_cost ); }

 /// returns the battery investment cost, i.e., the capital expenditure cost
 double get_batt_investment_cost( void ) const {
  return( f_batt_investment_cost );
 }

 /// returns the converter investment cost, i.e., the capital expenditure cost
 double get_converter_investment_cost( void ) const {
  return( f_conv_investment_cost );
 }

 /// returns the replacement cost of the battery at the end of the lifetime
 double get_batt_replacement_cost( void ) const {
  return( f_batt_replacement_cost );
 }

 /// returns the replacement cost of the battery at the end of the lifetime
 double get_converter_replacement_cost( void ) const {
  return( f_conv_replacement_cost );
 }

 /// returns the residual value of the battery at the end of the lifetime
 double get_batt_residual_value( void ) const {
  return( f_batt_residual_value );
 }

 /// returns the residual value of the converter at the end of the lifetime
 double get_converter_residual_value( void ) const {
  return( f_conv_residual_value );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of minimum storage
 /** This method returns a vector V containing the minimum storage at each
  * time instant. There are three possible cases:
  *
  * - if the vector is empty, then the minimum storage of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the minimum storage of
  *   the unit for all time instants;
  *
  * - otherwise, the vector must have size get_time_horizon() and each V[ t ]
  *   represents the minimum storage value at time t.
  *
  * @return The vector containing the minimum storage. */

 const std::vector< double > & get_minimum_storage( void ) const {
  return( v_minimum_storage );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum storage at the given time instant
 /** This method returns the minimum storage at the given time instant \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The minimum storage at the given time instant \p t. */

 double get_minimum_storage( Index t ) const {
  if( v_minimum_storage.empty() )
   return( 0 );
  if( v_minimum_storage.size() == 1 )
   return( v_minimum_storage.front() );
  assert( t < v_minimum_storage.size() );
  return( v_minimum_storage[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum storage
 /** This method returns a vector V containing the maximum storage at all time
  * instants. There are three possible cases:
  *
  * - if the vector is empty, then the maximum storage of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the maximum storage of
  *   the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each V[ t ]
  *   represents the maximum storage value at time t. */

 const std::vector< double > & get_maximum_storage( void ) const {
  return( v_maximum_storage );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum storage at the given time instant
 /** This method returns the maximum storage at the given time instant \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The maximum storage at the given time instant \p t. */

 double get_maximum_storage( Index t ) const {
  if( v_maximum_storage.empty() )
   return( 0 );
  if( v_maximum_storage.size() == 1 )
   return( v_maximum_storage.front() );
  assert( t < v_maximum_storage.size() );
  return( v_maximum_storage[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of minimum power
 /** This method returns a vector V containing the minimum power at all time
  * instants. There are three possible cases:
  *
  * - if the vector is empty, then the minimum power of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the minimum power of
  *   the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each V[ t
  *   ] represents the minimum power value at time t. */

 const std::vector< double > & get_minimum_power( void ) const {
  return( v_minimum_power );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power at the given time instant
 /** This method returns the minimum power at the given time \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The minimum power at the given time instant \p t. */

 double get_minimum_power( Index t ) const {
  if( v_minimum_power.empty() )
   return( 0 );
  if( v_minimum_power.size() == 1 )
   return( v_minimum_power.front() );
  assert( t < v_minimum_power.size() );
  return( v_minimum_power[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum power
 /** This method returns a vector V containing the maximum power at time all
  * time instants. There are three possible cases:
  *
  * - if the vector is empty, then the maximum power of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the maximum power of
  *   the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each V[ t
  *   ] represents the maximum power value at time t. */

 const std::vector< double > & get_maximum_power( void ) const {
  return( v_maximum_power );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power at the given time instant
 /** This method returns the maximum power at the given time \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The maximum power at the given time instant \p t. */

 double get_maximum_power( Index t ) const {
  if( v_maximum_power.empty() )
   return( 0 );
  if( v_maximum_power.size() == 1 )
   return( v_maximum_power.front() );
  assert( t < v_maximum_power.size() );
  return( v_maximum_power[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum primary reserve power
 /** This method returns a vector containing the maximum active power that can
  * be used as primary reserve. There are three possible cases:
  *
  * - if this vector is empty, then the maximum primary power of the
  *   unit is 0;
  *
  * - if this vector has only one element, then the maximum primary power is
  *   equal to that value at all time instants;
  *
  * - otherwise, the vector must have size get_time_horizon() and its t-th
  *   element represents the maximum primary power at time t. */

 const std::vector< double > & get_maximum_primary_power( void ) const {
  return( v_maximum_primary_rho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum primary reserve power at the given time instant
 /** This method returns the maximum active power that can be used as primary
  * reserve at the given time \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The maximum primary reserve power at the given time instant \p
  * t. */

 double get_maximum_primary_power( Index t ) const {
  if( v_maximum_primary_rho.empty() )
   return( 0 );
  if( v_maximum_primary_rho.size() == 1 )
   return( v_maximum_primary_rho.front() );
  assert( t < v_maximum_primary_rho.size() );
  return( v_maximum_primary_rho[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum secondary reserve power
 /** This method returns a vector containing the maximum active power that can
  * be used as secondary reserve. There are three possible cases:
  *
  * - if this vector is empty, then the maximum secondary power of the
  *   unit is 0;
  *
  * - if this vector has only one element, then the maximum secondary power is
  *   equal to that value at all time instants;
  *
  * - otherwise, the vector must have size get_time_horizon() and its t-th
  *   element represents the maximum secondary power at time t. */

 const std::vector< double > & get_maximum_secondary_power( void ) const {
  return( v_maximum_secondary_rho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum secondary power at the given time instant
 /** This method returns the maximum active power that can be used as
  * secondary reserve at the given time \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The maximum secondary reserve power at the given time instant \p
  * t. */

 double get_maximum_secondary_power( Index t ) const {
  if( v_maximum_secondary_rho.empty() )
   return( 0 );
  if( v_maximum_secondary_rho.size() == 1 )
   return( v_maximum_secondary_rho.front() );
  assert( t < v_maximum_secondary_rho.size() );
  return( v_maximum_secondary_rho[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of delta ramp up
 /** This method returns a V containing the delta ramp up at time all time
  * instants. There are three possible cases:
  *
  * - if the vector is empty, then the delta ramp up of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the delta ramp up of
  *   the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the delta ramp up value at time t. */

 const std::vector< double > & get_delta_ramp_up( void ) const {
  return( v_delta_ramp_up );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of delta ramp down
 /** This method returns a vector V containing the delta ramp down at time all
  * time instants. There are three possible cases:
  *
  * - if the vector is empty, then the delta ramp down of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the delta ramp down of
  *   the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the delta ramp down value at time t. */

 const std::vector< double > & get_delta_ramp_down( void ) const {
  return( v_delta_ramp_down );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inefficiency of storing energy
 /** This method returns a vector V containing the storing battery rho
  * (inefficiency of storing energy) at all time instants. There are three
  * possible cases:
  *
  * - if the vector is empty, then the storing battery of the unit is 1;
  *
  * - if the vector has only one element, then V[ 0 ] is the storing battery
  *   of the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the storing battery rho value at time t. */

 const std::vector< double > & get_storing_battery_rho( void ) const {
  return( v_storing_battery_rho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inefficiency of extracting energy of the unit
 /** This method returns a vector V containing the extracting battery rho
  * (inefficiency of extracting energy of the unit) at all time
  * instants. There are three possible cases:
  *
  * - if the vector is empty, then the extracting battery rho of the unit is 1;
  *
  * - if the vector has only one element, then V[ 0 ] is the extracting battery
  *   rho of the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the extracting battery rho value at time t. */

 const std::vector< double > & get_extracting_battery_rho( void ) const {
  return( v_extracting_battery_rho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of storage and extraction cost of energy
 /** This method returns a V containing the monetary cost of storing one unit
  * or energy into, or extracting it from, the battery (the cost is the same
  * in both cases) at all time instants. There are three possible cases:
  *
  * - if the vector is empty, then the cost of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the cost of the unit
  *   for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the cost of the unit at time t. */

 const std::vector< double > & get_cost( void ) const {
  return( v_cost );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of E-mobility demand
 /** This method returns a vector V containing the demand at all time
  * instants. There are two possible cases:
  *
  * - if the vector is empty, then the demand of the unit is 0;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the demand value at time t. */

 const std::vector< double > & get_demand( void ) const {
  return( v_demand );
 }

/*--------------------------------------------------------------------------*/
 /// returns the scale factor of this BatteryUnitBlock
 double get_scale( void ) const override { return( f_scale ); }

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE Variable OF THE BatteryUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the BatteryUnitBlock
 *
 * These methods allow to read the two groups of Variable that any
 * BatteryUnitBlock in principle has (although some may not):
 *
 * - the storage level variables
 *
 * - the intake and outtake variables
 *
 * All these two groups of variables are (if not empty)
 * std::vector< ColVariable > with the dimension time horizon.
 * @{ */

 /// returns the vector of storage level variables
 /** This method returns a vector V containing the storage level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the storage
  *   level variable for time step t. */

 const std::vector< ColVariable > & get_storage_level( void ) const {
  return( v_storage_level );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of intake level variables
 /** This method returns a vector V containing the intake level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size of get_time_horizon() and V[ t ] is the
  *   intake level variable for time step t. */

 const std::vector< ColVariable > & get_intake_level( void ) const {
  return( v_intake_level );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of outtake level variables
 /** This method returns a vector V containing the outtake level variables.
  * There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size of get_time_horizon() and V[ t ] is the
  *   outtake level variable for time step t.*/

 const std::vector< ColVariable > & get_outtake_level( void ) const {
  return( v_outtake_level );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of battery binary variables
 /** This method returns a vector V containing the battery binary variables.
  * There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size of get_time_horizon() and V[ t ] is the
  *   battery binary variable for time step t. */

 const std::vector< ColVariable > & get_battery_binary( void ) const {
  return( v_battery_binary );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active power variables

 ColVariable * get_active_power( Index generator ) override {
  if( v_active_power.empty() )
   return( nullptr );
  return( &( v_active_power.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary spinning reserve variables

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( v_primary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_primary_spinning_reserve.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary spinning reserve variables

 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( v_secondary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_secondary_spinning_reserve.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of battery_design variables, or nullptr if not defined

 ColVariable * get_battery_design( void ) {
  if( v_battery_design.empty() )
   return( nullptr );
  return( &( v_battery_design.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of converter_design variables, or nullptr if not defined

 ColVariable * get_converter_design( void ) {
  if( v_converter_design.empty() )
   return( nullptr );
  return( &( v_converter_design.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the type of this battery unit
 /** This method returns the type of this battery unit. */

 battery_type get_battery_type( void ) const {

  if( std::all_of( v_storing_battery_rho.cbegin() ,
                   v_storing_battery_rho.cend() ,
                   []( double s ) { return( s <= 1.0 ); } ) &&
      std::all_of( v_extracting_battery_rho.cbegin() ,
                   v_extracting_battery_rho.cend() ,
                   []( double s ) { return( s >= 1.0 ); } ) )
   return( ASSUME_POSITIVE_PRICES );

  if( ( ! v_storing_battery_rho.empty() ) &&
      ( ! v_extracting_battery_rho.empty() ) )
   return( NO_Binary_Variables_Constraints );

  return( Binary_Variables_Constraints );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power output constraints
 /** This method returns the minimum power output constraints.
  *
  * @return The minimum power output constraints. */

 const std::vector< FRowConstraint > & get_min_power_constraints( void ) const {
  return( active_power_lower_bound_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power output constraints
 /** This method returns the maximum power output constraints.
  *
  * @return The maximum power output constraints. */

 const std::vector< FRowConstraint > & get_max_power_constraints( void ) const {
  return( active_power_upper_bound_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraints
 /** This method returns the intake upper bound constraints.
  *
  * @return The intake upper bound constraints. */

 const std::vector< BoxConstraint > & get_max_intake_constraints( void ) const {
  return( intake_upper_bound_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraints with binary variables
 /** This method returns the intake upper bound constraints with binary
  * variables.
  *
  * @return The intake upper bound constraints with binary variables. */

 const std::vector< FRowConstraint > &
 get_max_intake_binary_const( void ) const {
  return( intake_binary_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraints with binary variables
 /** This method returns the outtake upper bound constraints with binary
  * variables.
  *
  * @return The outtake upper bound constraints with binary variables. */

 const std::vector< FRowConstraint > &
 get_max_outtake_binary_const( void ) const {
  return( outtake_binary_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the storage level bound constraints
 /** This method returns the storage level bound constraints.
  *
  * @return The storage level bound constraints. */

 const std::vector< BoxConstraint > &
 get_storage_level_bound_constraints( void ) const {
  return( storage_level_bounds_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the primary reserve bound constraints
 /** This method returns the primary reserve bound constraints.
  *
  * @return The primary reserve bound constraints. */

 const std::vector< BoxConstraint > & get_primary_reserve_bounds( void ) const {
  return( primary_upper_bound_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the secondary reserve bound constraints
 /** This method returns the secondary reserve bound constraints.
  *
  * @return The secondary reserve bound constraints. */

 const std::vector< BoxConstraint > &
 get_secondary_reserve_bounds( void ) const {
  return( secondary_upper_bound_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake/outtake binary variables
 /** This method returns the intake/outtake binary variables.
  *
  * @return The intake/outtake binary variables. */

 const std::vector< ColVariable > &
 get_intake_outtake_binary_variables( void ) const {
  return( v_battery_binary );
  }

/** @} ---------------------------------------------------------------------*/
/*---------------- METHODS FOR SAVING THE BatteryUnitBlock------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the BatteryUnitBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * BatteryUnitBlock. See BatteryUnitBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE BatteryUnitBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the BatteryUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "BatteryUnitBlock::load not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_initial_storage( std::vector< double >::const_iterator it ,
                           Subset && subset , bool ordered = false ,
                           c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck );

 void set_initial_storage( std::vector< double >::const_iterator it ,
                           Range rng = Range( 0 , Inf< Index >() ) ,
                           c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given \p subset contains the 0 index, this function sets the
  * initial power. If the given \p subset does not contain the index 0, this
  * function does nothing. Since \p subset can have multiple zeros, only the
  * last one is considered, which means that the value for the initial power
  * will be that in the vector pointed by \p it associated with this last
  * zero.
  */

 void set_initial_power( std::vector< double >::const_iterator it ,
                         Subset && subset ,
                         const bool ordered = false ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given Range \p rng contains 0, this function sets the initial
  * power. In this case, if the first element of \p rng is 0, the initial
  * power will be set to the value pointed by the given iterator. In general,
  * the initial power will be the one found at position -rng.first in the
  * vector pointed by \p it if this Range contains the 0 index. If the given
  * Range \p rng does not contain the 0 index, this function does nothing.
  */

 void set_initial_power( std::vector< double >::const_iterator it ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the scale factor of this BatteryUnitBlock
 /** This method sets the scale factor of this BatteryUnitBlock.
  *
  * @param values An iterator to a vector containing the scale factor.
  *
  * @param subset If non-empty, the scale factor is set to the value
  *        pointed by \p values. If empty, no operation is performed.
  *
  * @param ordered This parameter is ignored.
  *
  * @param issuePMod Controls how physical Modification are issued.
  *
  * @param issueAMod Controls how abstract Modification are issued. */

 void scale( std::vector< double >::const_iterator values ,
             Subset && subset ,
             const bool ordered = false ,
             c_ModParam issuePMod = eNoBlck ,
             c_ModParam issueAMod = eNoBlck ) override;

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets the kappa constant, which multiplies the minimum and
  * maximum active power, maximum primary and secondary reserve, and the
  * minimum and maximum storage levels in the constraints of this
  * BatteryUnitBlock.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param subset If non-empty, the kappa constant is set to the value pointed
  *        by \p values. If empty, no operation is performed.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_kappa( std::vector< double >::const_iterator values ,
                 Subset && subset , const bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets the kappa constant, which multiplies the minimum and
  * maximum active power, maximum primary and secondary reserve, and the
  * minimum and maximum storage levels in the constraints of this
  * BatteryUnitBlock.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param rng If non-empty, the kappa constant is set to the value pointed by
  *        \p values. If empty, no operation is performed.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_kappa( std::vector< double >::const_iterator values ,
                 Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets the kappa constant, which multiplies the minimum and
  * maximum active power, maximum primary and secondary reserve, and the
  * minimum and maximum storage levels in the constraints of this
  * BatteryUnitBlock.
  *
  * @param value The value of the kappa constant.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_kappa( double value , c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector = { value };
  set_kappa( vector.cbegin() , Range( 0 , Inf< Index >() ) ,
             issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

 // For the Range version, use the default implementation defined in UnitBlock
 using UnitBlock::scale;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 static void static_initialization( void ) {

  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args(),
   *
   * register_method< BatteryUnitBlock >( "BatteryUnitBlock::set_initial_storage",
   *                                      &BatteryUnitBlock::set_initial_storage,
   *                                      MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< BatteryUnitBlock, MF_dbl_it , Subset && , bool >(
   "BatteryUnitBlock::set_initial_storage" ,
   &BatteryUnitBlock::set_initial_storage );

  register_method< BatteryUnitBlock , MF_dbl_it , Range >(
   "BatteryUnitBlock::set_initial_storage" ,
   &BatteryUnitBlock::set_initial_storage );

  register_method< BatteryUnitBlock , MF_dbl_it , Subset && , bool >(
   "BatteryUnitBlock::set_kappa" ,
   &BatteryUnitBlock::set_kappa );

  register_method< BatteryUnitBlock , MF_dbl_it , Range >(
   "BatteryUnitBlock::set_kappa" ,
   &BatteryUnitBlock::set_kappa );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the vector of minimum storage
 std::vector< double > v_minimum_storage;

 /// the vector of maximum storage
 std::vector< double > v_maximum_storage;

 /// the vector of MinPower
 std::vector< double > v_minimum_power;

 /// the vector of MaxPower
 std::vector< double > v_maximum_power;

 /// the vector of MaxPrimaryRho
 std::vector< double > v_maximum_primary_rho;

 /// the vector of MaxSecondaryRho
 std::vector< double > v_maximum_secondary_rho;

 /// the vector of RampUp
 std::vector< double > v_delta_ramp_up;

 /// the vector of RampDown
 std::vector< double > v_delta_ramp_down;

 /// the vector of storing battery rho
 std::vector< double > v_storing_battery_rho;

 /// the vector of extracting battery rho
 std::vector< double > v_extracting_battery_rho;

 /// the vector of Cost
 std::vector< double > v_cost;

 /// the vector of demand
 std::vector< double > v_demand;

 /// the InitialStorage value
 double f_initial_storage{};

 /// the InitialPower value
 double f_initial_power{};

 /// the kappa value
 double f_kappa = 1;

 /// the scale factor of this BatteryUnitBlock
 double f_scale = 1;

 /// the operation and maintenance costs
 double f_oem_cost{};

 /// the battery investment cost, i.e., the capital expenditure cost
 double f_batt_investment_cost{};

 /// the converter investment cost, i.e., the capital expenditure cost
 double f_conv_investment_cost{};

 /// the replacement cost of the battery at the end of the lifetime
 double f_batt_replacement_cost{};

 /// the replacement cost of the converter at the end of the lifetime
 double f_conv_replacement_cost{};

 /// the residual value of the battery at the end of the lifetime
 double f_batt_residual_value{};

 /// the residual value of the converter at the end of the lifetime
 double f_conv_residual_value{};

/*-------------------------------- variables -------------------------------*/

 /// the vector of storage level variables
 std::vector< ColVariable > v_storage_level;

 /// the vector of intake level variables
 std::vector< ColVariable > v_intake_level;

 /// the vector of outtake level variables
 std::vector< ColVariable > v_outtake_level;

 /// the vector of binary variables
 std::vector< ColVariable > v_battery_binary;

 /// the active power variables
 std::vector< ColVariable > v_active_power;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

 /// the battery design binary variables
 std::vector< ColVariable > v_battery_design;

 /// the converter design binary variables
 std::vector< ColVariable > v_converter_design;

/*------------------------------- constraints ------------------------------*/

 /// the active power upper bound constraints
 std::vector< FRowConstraint > active_power_upper_bound_Const;

 /// the active power lower bound constraints
 std::vector< FRowConstraint > active_power_lower_bound_Const;

 /// the ramp up constraints
 std::vector< FRowConstraint > ramp_up_Const;

 /// the ramp down constraints
 std::vector< FRowConstraint > ramp_down_Const;

 /// the active power, intake and outtake relation constraints
 std::vector< FRowConstraint > power_intake_outtake_Const;

 /// the storage , intake and outtake level relation constraints
 std::vector< FRowConstraint > storage_intake_outtake_Const;

 /// the intake and binary variable relation constraints
 std::vector< FRowConstraint > intake_binary_Const;

 /// the outtake and binary variable relation constraints
 std::vector< FRowConstraint > outtake_binary_Const;

 /// the demand constraints
 std::vector< FRowConstraint > demand_Const;


 /// the storage level bounds constraints
 std::vector< BoxConstraint > storage_level_bounds_Const;

 /// the intake upper bound constraints
 std::vector< BoxConstraint > intake_upper_bound_Const;

 /// primary upper bound constraints
 std::vector< BoxConstraint > primary_upper_bound_Const;

 /// secondary upper bound constraints
 std::vector< BoxConstraint > secondary_upper_bound_Const;


 /// the vector of binary variables
 std::vector< ZOConstraint > battery_binary_bound_Const;


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

 /// updates the constraints for the current initial storage
 /** This function updates both sides of the demand constraint at time 0
  * (which is the constraint that depends on the initial storage). */

 void update_initial_storage_in_constraints( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current initial power
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints at time 0 (which are the
  * constraints that depend on the initial power). */

 void update_initial_power_in_constraints( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current kappa
 /** This function updates the constraints to take into account the current
  * value of the kappa constant. */

 void update_kappa_in_constraints( ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the coefficients of the Objective
 /** This method updates the coefficients of the Objective.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective( c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this BatteryUnitBlock is consistent
 /** This function checks whether the data in this BatteryUnitBlock is
  * consistent. The data is consistent if all of the following conditions are
  * met.
  *
  * - The maximum power is greater than or equal to the minimum power.
  *
  * - The maximum storage level is greater than or equal to the minimum
  *   storage level.
  *
  * - The minimum storage level is nonnegative.
  *
  * - The inefficiency of storing energy is less than or equal to 1.
  *
  * - The inefficiency of extracting energy is greater than or equal to 1.
  *
  * - The inefficiency of extracting energy is greater than or equal to the
  *   inefficiency of storing energy.
  *
  * - The demand is nonnegative.
  *
  * - The maximum active power that can be used as primary and secondary
  *   reserves are nonnegative.
  *
  * - The initial storage is nonnegative.
  *
  * If any of the above conditions are not met, an exception is thrown. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS BatteryUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a BatteryUnitBlock
 class BatteryUnitBlockMod : public UnitBlockMod
 {

  public:

  /// public enum for the types of BatteryUnitBlockMod
  enum BUB_mod_type
  {
   eSetInitS = eUBModLastParam , ///< set initial storage values
   eSetInitP ,                   ///< set initial power values
   eSetKappa ,                   ///< set the kappa constant
  };

  /// constructor, takes the BatteryUnitBlock and the type
  BatteryUnitBlockMod( BatteryUnitBlock * const fblock , const int type )
   : UnitBlockMod( fblock , type ) {}

  /// destructor, does nothing
  virtual ~BatteryUnitBlockMod() override = default;

  /// returns the Block to which the Modification refers
  Block * get_Block( void ) const override { return( f_Block ); }

  protected:

  /// prints the BatteryUnitBlockMod
  void print( std::ostream & output ) const override {
   output << "BatteryUnitBlockMod[" << this << "]: ";
   switch( f_type ) {
    case( eSetInitS ):
     output << "set initial storage values ";
     break;
    case( eSetInitP ):
     output << "set initial power values ";
     break;
    case( eSetKappa ):
     output << "set kappa ";
     break;
   }
  }

  BatteryUnitBlock * f_Block{};
  ///< pointer to the Block to which the Modification refers

 };  // end( class( BatteryUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS BatteryUnitBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from BatteryUnitBlockMod for "ranged" modifications
 class BatteryUnitBlockRngdMod : public BatteryUnitBlockMod
 {

  public:

  /// constructor: takes the BatteryUnitBlock, the type, and the range
  BatteryUnitBlockRngdMod( BatteryUnitBlock * const fblock ,
                           const int type ,
                           Block::Range rng )
   : BatteryUnitBlockMod( fblock , type ) , f_rng( rng ) {}

  /// destructor, does nothing
  virtual ~BatteryUnitBlockRngdMod() override = default;

  /// accessor to the range
  Block::c_Range & rng( void ) { return( f_rng ); }

  protected:

  /// prints the BatteryUnitBlockRngdMod
  void print( std::ostream & output ) const override {
   BatteryUnitBlockMod::print( output );
   output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

  Block::Range f_rng; ///< the range

 };  // end( class( BatteryUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS BatteryUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived from BatteryUnitBlockMod for "subset" modifications
 class BatteryUnitBlockSbstMod : public BatteryUnitBlockMod
 {

  public:

  /// constructor: takes the BatteryUnitBlock, the type, and the subset
  BatteryUnitBlockSbstMod( BatteryUnitBlock * const fblock ,
                           const int type ,
                           Block::Subset && nms )
   : BatteryUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

  /// destructor, does nothing
  virtual ~BatteryUnitBlockSbstMod() override = default;

  /// accessor to the subset
  Block::c_Subset & nms( void ) { return( f_nms ); }

  protected:

  /// prints the BatteryUnitBlockSbstMod
  void print( std::ostream & output ) const override {
   BatteryUnitBlockMod::print( output );
   output << "(# " << f_nms.size() << ")" << std::endl;
  }

  Block::Subset f_nms; ///< the subset

 };  // end( class( BatteryUnitBlockSbstMod ) )

};  // end( class( BatteryUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* BatteryUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File BatteryUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
