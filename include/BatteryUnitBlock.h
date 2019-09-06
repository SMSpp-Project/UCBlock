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
 * \version 0.11
 *
 * \date 02 - 09 - 2019
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

/// Namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS BatteryUnitBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
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
 * management). BatteryUnitBlock provides a quite general concept of "battery"
 * that covers different units which mostly fit the same mathematical
 * equations pattern. For instance, a BatteryUnitBlock may or may not have a
 * fixed demand (e-mobility has, other units have not) and it may or may not
 * provide primary and secondary reserve (battery storage may do, but other
 * units don't).
 *
 * Battery storage provide an additional flexibility to the system by shifting
 * a surplus of electric energy (e.g. due to high renewable feeding) to times
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
 * much longer timespan (e.g. at night). This allows to shift the charging
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
 *   secondary spinning reserves(if any);
 *
 * - the ramp-up and ramp-down constraints;
 *
 * - the active power relation with storing and extracting energy levels
 *   constraints(if any);
 *
 * - the intake upper bound(if any);
 *
 * - the storage level constraints;
 *
 * - the binary variable relation with intake and outtake level constraints(if
 *   any);
 *
 * - the primary reserve upper bound(if any);
 *
 * - the secondary reserve upper bound(if any);
 *
 * - the demand constraints(if any).*/

class BatteryUnitBlock : public UnitBlock {

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
/** Constructor of BatteryUnitBlock, taking possibly a pointer of its father
 * Block.
 */

 explicit BatteryUnitBlock( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/

/// destructor of BatteryUnitBlock

 ~BatteryUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the BatteryUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ). In
 * particular, we refer to that description for the crucial dimensions
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
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must
 *   always be 0 <= MinS[ t ] < MaxS[ t ] for all t. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxStorage", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MaxS[ t ] that, for each time instant t, contains the maximum
 *   storage level of the unit for the corresponding time step. If
 *   "MaxStorage" has length 1 then MaxS[ t ] contains the same value for all
 *   t. Otherwise, MaxStorage[ i ] is the fixed value of MaxS[ t ] for all t
 *   in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must always
 *   be [0 <=] MinS[ t ] < MaxS[ t ] for all t. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MinPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MinP[ t ] that, for each time instant t, contains the minimum
 *   active power output value of the unit for the corresponding time step.
 *   If "MinPower" has length 1 then MinP[ t ] contains the same value for all
 *   t. Otherwise, MinPower[ i ] is the fixed value of MinP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must
 *   be MinP[ t ] <= 0 for all t. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MaxP[ t ] that, for each time instant t, contains the maximum
 *   active power output value of the unit for the corresponding time step.
 *   If "MaxPower" has length 1 then MaxP[ t ] contains the same value for all
 *   t. Otherwise, MaxPower[ i ] is the fixed value of MaxP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. Note that it must
 *   be MinP[ t ] < MaxP[ t ] for all t. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The scalar variable "InitialPower", of type double and not indexed over
 *   any dimension. This variable indicates the amount of the power that the
 *   unit was producing at time instant -1, i.e., before the start of the time
 *   horizon; this is necessary to compute the ramp-up and ramp-down
 *   constraints. This variable is optional; if "DeltaRampUp" and
 *   "DeltaRampDown" are not present, "InitialPower" should not be read,
 *   since there are no ramping constraints. If "DeltaRampUp" and
 *   "DeltaRampDown" are present but "InitialPower" is not provided, its
 *   value is taken to be 0.
 *
 * - The variable "MaxPrimaryPower", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector MaxPP[ t ] that, for each time instant t, contains the maximum
 *   active power that can be used as primary reserve of the unit for the
 *   corresponding time step. If "MaxPrimaryPower" has length 1 then
 *   MaxPP[ t ] contains the same value for all t. Otherwise,
 *   MaxPrimaryPower[ i ] is the fixed value of MaxPP[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional,
 *   if is not provided then MaxPP[ t ] == 0 for all t. If NumberIntervals
 *   <= 1 or NumberIntervals >= TimeHorizon, then the mapping clearly does
 *   not require "ChangeIntervals", which in fact is not loaded. 
 *
 * - The variable "MaxSecondaryPower", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector MaxSP[ t ] that, for each time instant t, contains the maximum
 *   active power that can be used as secondary reserve of the unit for the
 *   corresponding time step. If "MaxSecondaryPower" has length 1 then
 *   MaxSP[ t ] contains the same value for all t. Otherwise,
 *   MaxSecondaryPower[ i ] is the fixed value of MaxSP[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional,
 *   if is not provided then MaxSP[ t ] == 0 for all t. Note that
 *   MaxPP[ t ] == 0 implies MaxSP[ t ] == 0 (that is, if MaxPrimaryPower is
 *   not defined then neither should MaxSecondaryPower). If NumberIntervals
 *   <= 1 or NumberIntervals >= TimeHorizon, then the mapping clearly does
 *   not require "ChangeIntervals", which in fact is not loaded. 
 *
 * - The variable "DeltaRampUp", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector DP[ t ] that, for each time instant t, contains the ramp-up
 *   value of the unit for the corresponding time step, i.e., the maximum
 *   possible increase of active power production w.r.t. the power that had
 *   been produced in time instant t - 1, if any. If "DeltaRampUp" has length
 *   1 then DP[ t ] contains the same value for all t. Otherwise,
 *   DeltaRampUp[ i ] is the fixed value of DP[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. This variable is optional; if it is not
 *   provided then it is assumed that DP[ t ] == MaxP[ t ], i.e., the unit
 *   can ramp up by an arbitrary amount, i.e., there are no ramp-up
 *   constraints. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon,
 *   then the mapping clearly does not require "ChangeIntervals", which in
 *   fact is not loaded.
 *
 * - The variable "DeltaRampDown", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector DM[ t ] that, for each time instant t, contains the ramp-down
 *   value of the unit for the corresponding time step, i.e., the maximum
 *   possible decrease of active power production w.r.t. the power that had
 *   been produced in time instant t - 1, if any. If "DeltaRampDown" has
 *   length 1 then DM[ t ] contains the same value for all t. Otherwise,
 *   DeltaRampDown[ i ] is the fixed value of DM[ t ] for all t in the 
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional;
 *   if it is not provided then it is assumed that DM[ t ] == MaxP[ t ], i.e.,
 *   the unit can ramp down an arbitrary amount, i.e., there are no
 *   ramp-down constraints. If NumberIntervals <= 1 or NumberIntervals >=
 *   TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "StoringBatteryRho", of type double and to be either of size
 *   1 or indexed over the dimension "NumberIntervals". This is meant to
 *   represent the vector SBR[ t ] that, for each time instant t, contains the
 *   inefficiency of storing energy of the unit for the corresponding time
 *   step. This variable is optional; if it is not provided then it is
 *   assumed that SBR[ t ] == 1 for all t, i.e., no (significant) energy is
 *   spent just for storing it in the battery (this simplifies the model
 *   somewhat, see below). If "StoringBatteryRho" has length 1 then
 *   SBR[ t ] contains the same value for all t. Otherwise,
 *   StoringBatteryRho[ i ] is the fixed value of SBR[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ] with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. Note that it must be always
 *   such that 1 <= SBR[ t ] for all t (as SBR[ t ] is the amount of energy
 *   that has to be used to store 1 unit of energy in the battery). If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "ExtractingBatterRho", of type double and to be either of
 *   size 1 or indexed over the dimension "NumberIntervals". This is meant to
 *   represent the vector EBR[ t ] that, for each time instant t, contains the
 *   inefficiency of extracting energy of the unit for the corresponding time
 *   step. This variable is optional; if it is not provided, then it is
 *   assumed that EBR[ t ] == 1 for all t, i.e., no (significant) energy is
 *   spent just for extracting it from the battery (this simplifies the model
 *   somewhat, see below). If "ExtractingBatterRho" has length 1 then
 *   EBR[ t ] contains the same value for all t. Otherwise,
 *   ExtractingBatterRho[ i ] is the fixed value of EBR[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ] with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. Note that it must be always
 *   such that EBR[ t ] <= 1 [<= SBR[ t ]] for all t (as EBR[ t ] is the
 *   amount of energy that is obtained when 1 unit of energy is removed from
 *   the battery). If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon,
 *   then the mapping clearly does not require "ChangeIntervals", which in
 *   fact is not loaded.
 *
 * Note: the special case in which EBR[ t ] == SBR[ t ] == 1 for all t, i.e.,
 * no energy is spent for storing it in / retrieving it from the battery,
 * leads to significantly simpler mathematical models. In particular, one
 * single variable can be used to represent both storing and retrieving,
 * rather than requiring two separate ones (unless primary and secondary
 * reserve are allowed and/or the cost is defined, since this also requires
 * using two), and the binary variables need not to be defined. For details,
 * see the comments to generate_abstract_variable() and
 * generate_abstract_constraints().
 *
 * - The scalar variable "InitialStorage", of type double and not indexed over
 *   any dimension. This variable indicates the amount of the storage level
 *   that the unit was producing at time instant -1, i.e., before the start of
 *   the time horizon; this is necessary to compute the storage level
 *   connection with intake and outtake constraints.
 *
 * - The variable "Cost", of type double and either of size 1 or indexed over
 *   the dimension "NumberIntervals". This is meant to represent the vector
 *   C[ t ] that, for each time instant t, contains the monetary cost of
 *   storing one unit or energy into, or extracting it from, the battery
 *   (the cost is the same in both cases) at the corresponding time step.
 *   This variable is optional; if it is not provided then it's taken to be
 *   zero. If "Cost" has length 1 then C[ t ] contains the same value for
 *   all t. Otherwise, Cost[ i ] is the fixed value of C[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "Demand", of type double and indexed over the dimension
 *   "TimeHorizon": the entry Demand[ t ] is assumed to contain the amount of
 *   energy that must be discharged from the battery and "sent away for some
 *   other purpose" (say, driving your e-car) at time t. This variable is
 *   optional; if it isn't defined, then Demand[ t ] == 0. Otherwise the
 *   Demand[ t ] contains the demand value for each time instant t.
 * */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the BatteryUnitBlock
/** The BatteryUnitBlock class use get_variable() method to access to
 *  each "group" of variable that may create in UnitBlock class which are:
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
 *    active power variable(if it's needed);
 *
 *  - the binary variables; when "StoringBatteryRho" == "ExtractingBatterRho"
 *    == 1, then this binary variable and all constraints which are depended
 *    on this variable is not required to be define.
 *
 *  These three groups of variables may have size f_time_horizon or empty size.
 *  All of these variables are optional, and it is also possible to restrict
 *  which of the subsets are generated with the parameter stvv. If stvv is not
 *  nullptr and it is a SimpleConfiguration<int>, or if
 *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 *  SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 *  of the optional variables should be created. If the Configuration is not
 *  available, the default value is taken to be 0.
 * */
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the BatteryUnitBlock
/** Method that generates the static constraint of the BatteryUnitBlock. The
 * operations of the battery storage unit are described on a discrete time
 * horizon as dictated by the UnitBlock interface. In this description we
 * indicate it with \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$. The
 * main constraints of this unit are define as:
 *
 * - maximum and minimum power output constraints according to primary and
 *   secondary spinning reserves are presented in (1)-(2). Each of them is a
 *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
 *   the entry t = 0, ...,f_time_horizon - 1 being the maximum and minimum
 *   power output value according to the primary and the secondary spinning
 *   reserves at time t. these ensure the maximum(or minimum) amount of energy
 *   that unit can produce(or use) when it is on(or off).
 *   \f[
 *      p^{ac}_{t} + p^{pr}_{t} + p^{sc}_{t} \leq P^{mx}_{t}
 *          \quad t \in \mathcal{T}                              \quad (1)
 *   \f]
 *
 *   \f[
 *     P^{mn}_{t} \leq p^{ac}_{t} - p^{pr}_{t} - p^{sc}_{t}
 *         \quad t \in \mathcal{T}                               \quad (2)
 *   \f]
 *   where \f$ P^{mx}_{t} \f$ and \f$ P^{mn}_{t} \f$ are the maximum and
 *   minimum power output parameters for each time t of the time horizon
 *   \f$ \mathcal{T} \f$ respectively.
 *
 * - ramp-up and ramp-down constraints are presented in (3)-(4). Each of them
 *   is a std::vector<FRowConstraint>; with the dimension of f_time_horizon,
 *   where the entry t = 0, ...,f_time_horizon - 1 being the ramp up and ramp
 *   down constraints which are presented as:
 *   \f[
 *    p^{ac}_{t} - p^{ac}_{t-1} \leq \Delta^{up}_{t}
 *         \quad t \in \mathcal{T}                               \quad (3)
 *   \f]
 *
 *   \f[
 *    p^{ac}_{t} - p^{ac}_{t-1} \geq - \Delta^{dn}_{t}
 *         \quad t \in \mathcal{T}                               \quad (4)
 *   \f]
 *   where \f$ \Delta^{up}_{t} \f$ and \f$ \Delta^{dn}_{t} \f$ are the delta
 *   ramp-up and delta ramp down threshold for each time t of the time horizon
 *   \f$ \mathcal{T} \f$ respectively.
 *
 * - active power relation with intake and outtake levels constraints are
 *   presented in (5). Each of them is a std::vector<FRowConstraint>; with the
 *   dimension of f_time_horizon, where the entry t = 0,...,f_time_horizon - 1
 *   being the active power relation with intake and outtake levels at time t.
 *   These ensure the active power at each time should be equal to the intake
 *   and outtake difference. The equation (6) also indicates the upper bound
 *   of intake level at each time instant t.
 *   \f[
 *    p^{ac}_{t} = p^+_t - p^-_{t}
 *         \quad t \in \mathcal{T}                               \quad (5)
 *   \f]
 *   \f[
 *     p^+_t \leq  P^{mx}_{t}
 *         \quad t \in \mathcal{T}                               \quad (6)
 *   \f]
 * - storage level relation with intake and outtake levels(if any) constraints
 *   in Battery unit are presented in (7). That is a
 *   std::vector<FRowConstraint>; with the dimension of f_time_horizon,
 *   where the entry t = 0,...,f_time_horizon - 1 being the storage level
 *   relation with intake and outtake levels at time t.
 *   \f[
 *    v^{ba}_{t} = v^{ba}_{t-1} - \rho^+_{t}p^+_{t} +
 *    \rho^-_{t}p^-_{t}-d^{ba}_t      \quad t \in \mathcal{T}    \quad (7)
 *   \f]
 *
 *   Note that if the equation (7) will change as
 *   below which is a std::vector<FRowConstraint>; with the dimension of
 *   f_time_horizon, where the entry t = 0,...,f_time_horizon - 1 being the
 *   storage level relation with battery demand(if any) at time t.
 *   \f[
 *    v^{ba}_{t} = v^{ba}_{t-1} - p^{ac}_{t} - d^{ba}_t
 *               \quad t \in \mathcal{T}          \quad (8)
 *   \f]
 *
 *   The equation (9) gives the storage levels upper bound and
 *   lower bound at each time instant t.
 *
 *   \f[
 *    v^{ba}_{t} \in [ V^{mn}_{t} , V^{mx}_{t}]
 *                              \quad t \in \mathcal{T}          \quad (9)
 *   \f]
 *   where \f$ \rho^+_{t} \f$ and \f$ \rho^-_{t} \f$ are the intake and
 *   outtake rho and \f$ V^{mn}_t\f$ and \f$ V^{mx}_t\f$ are the minimum and
 *   maximum storage level for each time t of the time horizon
 *   \f$ \mathcal{T} \f$ respectively.
 *
 * - binary variable relation with storing and extracting energy level(if any)
 *   constraints are presented in (10-11). Each of them is a
 *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
 *   the entry t = 0,...,f_time_horizon - 1 being the binary variable relation
 *   with storing and extracting energy levels at time t.
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
 *   Note that when \f$ \rho^+_t = \rho^-_t = 1 \f$, the binary variable
 *   \f$ u^+_t \f$ is not required and neither are the last two constraints
 *   (10-11).
 *
 * - primary and secondary reserve upper bounds(if any) are presented by the
 *   equations (12-13). Each of them is a std::vector<FRowConstraint>; with
 *   the dimension of f_time_horizon, where the entry
 *   t = 0,...,f_time_horizon - 1 being the primary and secondary reserve
 *   upper bounds at time t.
 *
 *   \f[
 *     p^{pr}_{t} \leq P^{mx, pr}_{t}
 *                              \quad t \in \mathcal{T}          \quad (12)
 *   \f]
 *   \f[
 *     p^{sc}_{t} \leq P^{mx, sc}_{t}
 *                              \quad t \in \mathcal{T}          \quad (13)
 *   \f]
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
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
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE BatteryUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the BatteryUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of battery storage units
 * @{ */

 /// Returns the initial storage value
 double get_initial_storage() const { return f_initial_storage; }

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
/// returns the vector of storing battery rho
/** The method returned a std::vector< double > V and each element of V
 * contains the storing battery at time t. There are three possible cases:
 *
 * - if the vector is empty, then the storing battery of the unit is 1;
 *
 * - if the vector has only one element, then V[ 0 ] is the storing battery of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the storing battery value at time t. */

 const std::vector< double > & get_storing_battery_rho() const {
  return( v_storing_battery_rho );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of extracting battery rho
/** The method returned a std::vector< double > V and each element of V
 * contains the extracting battery rho at time t. There are three possible
 * cases:
 *
 * - if the vector is empty, then the extracting battery rho of the unit is 1;
 *
 * - if the vector has only one element, then V[ 0 ] is the extracting battery
 *   rho of the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the extracting battery rho value at time t.
 *   */

 const std::vector< double > & get_extracting_battery_rho() const {
  return( v_extracting_battery_rho );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of cost
/** The method returned a std::vector< double > V and each element of V
 * contains the cost of the unit at time t. There are three possible cases:
 *
 * - if the vector is empty, then the cost of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the cost of the unit
 *   for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the cost value of the unit at time t. */

 const std::vector< double > & get_cost() const {
  return( v_cost );
 }

/*--------------------------------------------------------------------------*/
/// returns the vector of E-mobility demand
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
/*--------------------------------------------------------------------------*/
/// returns the vector of intake level variables
/** The returned std::vector< ColVariable >, say V, contains the
 * intake level variables and is indexed over the dimension time horizon.
 * There are two possible cases:
 *
 * - if V is empty(), then these variables are not defined;
 *
 * - otherwise, V must have size of get_time_horizon() and V[ t ] is the
 *   intake level variable for time step t.*/

 const std::vector< ColVariable > & get_intake_level() const {
  return v_intake_level;
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of outtake level variables
/** The returned std::vector< ColVariable >, say V, contains the
 * outtake level variables and is indexed over the dimension time horizon.
 * There are two possible cases:
 *
 * - if V is empty(), then these variables are not defined;
 *
 * - otherwise, V must have size of get_time_horizon() and V[ t ] is the
 *   outtake level variable for time step t.*/

 const std::vector< ColVariable > & get_outtake_level() const {
  return v_outtake_level;
 }
/**@} ----------------------------------------------------------------------*/
/*---------------- METHODS FOR SAVING THE BatteryUnitBlock------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the BatteryUnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * BatteryUnitBlock. See BatteryUnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE BatteryUnitBlock ------------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the BatteryUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "BatteryUnitBlock::load() not "
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

 /// the vector of intake rho
 std::vector< double >  v_storing_battery_rho;

 /// the vector of outtake rho
 std::vector< double >  v_extracting_battery_rho;

 /// the vector of Cost
 std::vector< double >  v_cost;

 /// the InitialStorage value
 double f_initial_storage;

 /// the InitialPower value
 double f_initial_power;

 /// the vector of demand
 std::vector< double >  v_demand;
/*-----------------------------variables------------------------------------*/
 /// the vector of storage level variables
 std::vector< ColVariable > v_storage_level;

 /// the vector of intake level variables
 std::vector< ColVariable > v_intake_level;

 /// the vector of outtake level variables
 std::vector< ColVariable > v_outtake_level;

 /// the vector of analogous variables
 std::vector< ColVariable > v_analogous;
/*----------------------------constraints-----------------------------------*/
/// the active power upper bound constraints
 std::vector< FRowConstraint > active_power_upper_bound_Constraints;

/// the active power lower bound constraints
 std::vector< FRowConstraint > active_power_lower_bound_Constraints;

/// the ramp up constraints
 std::vector< FRowConstraint > ramp_up_Constraints;

/// the ramp down constraints
 std::vector< FRowConstraint > ramp_down_Constraints;

/// the active power, intake and outtake relation constraints
 std::vector< FRowConstraint > power_intake_outtake_Constraints;

/// the intake upper bound constraints
 std::vector< FRowConstraint > intake_upper_bound_Constraints;

/// the storage , intake and outtake level relation constraints
 std::vector< FRowConstraint > storage_intake_outtake_Constraints;

/// the storage level bounds constraints
 std::vector< FRowConstraint > storage_level_bounds_Constraints;

/// the intake and binary variable relation constraints
 std::vector< FRowConstraint > intake_binary_Constraints;

/// the outtake and binary variable relation constraints
 std::vector< FRowConstraint > outtake_binary_Constraints;

/// the demand constraints
 std::vector< FRowConstraint > demand_Constraints;

/// primary upper bound constraints
 std::vector< FRowConstraint > primary_upperbound_Constraints;

/// secondary upper bound constraints
 std::vector< FRowConstraint > secondary_upperbound_Constraints;
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
