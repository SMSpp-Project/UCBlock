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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli
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
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father block
 /** Constructor of BatteryUnitBlock, taking possibly a pointer of its father
  * Block. */

 explicit BatteryUnitBlock( Block * f_block = nullptr )
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
  * - The variable "MinStorage", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MinS[ t ] that,
  *   for each time instant t, contains the minimum storage level of the unit
  *   for the corresponding time step. If "MinStorage" has length 1 then
  *   MinS[ t ] contains the same value for all t. Otherwise, MinStorage[ i ]
  *   is the fixed value of MinS[ t ] for all t in the interval
  *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. Note that it must always be
  *   0 <= MinS[ t ] < MaxS[ t ] for all t. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "MaxStorage", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MaxS[ t ] that,
  *   for each time instant t, contains the maximum storage level of the unit
  *   for the corresponding time step. If "MaxStorage" has length 1 then
  *   MaxS[ t ] contains the same value for all t. Otherwise, MaxStorage[ i ]
  *   is the fixed value of MaxS[ t ] for all t in the interval
  *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. Note that it must always be
  *   [0 <=] MinS[ t ] < MaxS[ t ] for all t. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "MinPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MinP[ t ] that,
  *   for each time instant t, contains the minimum active power output value
  *   of the unit for the corresponding time step. If "MinPower" has length 1
  *   then MinP[ t ] contains the same value for all t. Otherwise,
  *   MinPower[ i ] is the fixed value of MinP[ t ] for all t in the interval
  *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. Note that it must be MinP[ t ] <= 0 for
  *   all t. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
  *   the mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded. If "MinPower" is not provided, then MinP[ t ] is taken to
  *   be equal to -MaxP[ t ] for all t.
  *
  * - The variable "MaxPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MaxP[ t ] that,
  *   for each time instant t, contains the maximum active power output value
  *   of the unit for the corresponding time step. If "MaxPower" has length 1
  *   then MaxP[ t ] contains the same value for all t. Otherwise,
  *   MaxPower[ i ] is the fixed value of MaxP[ t ] for all t in the interval
  *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. Note that it must be
  *   MinP[ t ] < MaxP[ t ] for all t. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The scalar variable "InitialPower", of type netCDF::NcDouble and not
  *   indexed over any dimension. This variable indicates the amount of the
  *   power that the unit was producing at time instant -1, i.e., before the
  *   start of the time horizon; this is necessary to compute the ramp-up and
  *   ramp-down constraints. This variable is optional; if "DeltaRampUp" and
  *   "DeltaRampDown" are not present, "InitialPower" should not be read,
  *   since there are no ramping constraints. If "DeltaRampUp" and
  *   "DeltaRampDown" are present but "InitialPower" is not provided, its
  *   value is taken to be 0.
  *
  * - The variable "MaxPrimaryPower", of type netCDF::NcDouble and either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   MaxPP[ t ] that, for each time instant t, contains the maximum active
  *   power that can be used as primary reserve of the unit for the
  *   corresponding time step. If "MaxPrimaryPower" has length 1 then
  *   MaxPP[ t ] contains the same value for all t. Otherwise,
  *   MaxPrimaryPower[ i ] is the fixed value of MaxPP[ t ] for all t in the
  *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional,
  *   if is not provided then MaxPP[ t ] == 0 for all t. If NumberIntervals
  *   <= 1 or NumberIntervals >= TimeHorizon, then the mapping clearly does
  *   not require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "MaxSecondaryPower", of type netCDF::NcDouble and either
  *   of size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   MaxSP[ t ] that, for each time instant t, contains the maximum active
  *   power that can be used as secondary reserve of the unit for the
  *   corresponding time step. If "MaxSecondaryPower" has length 1 then
  *   MaxSP[ t ] contains the same value for all t. Otherwise,
  *   MaxSecondaryPower[ i ] is the fixed value of MaxSP[ t ] for all t in
  *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
  *   the assumption that ChangeIntervals[ - 1 ] = 0. This variable is
  *   optional, if is not provided then MaxSP[ t ] == 0 for all t. Note that
  *   MaxPP[ t ] == 0 implies MaxSP[ t ] == 0 (that is, if MaxPrimaryPower is
  *   not defined then neither should MaxSecondaryPower). If
  *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
  *   clearly does not require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "DeltaRampUp", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
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
  * - The variable "DeltaRampDown", of type netCDF::NcDouble and either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   DM[ t ] that, for each time instant t, contains the ramp-down value of
  *   the unit for the corresponding time step, i.e., the maximum possible
  *   decrease of active power production w.r.t. the power that had been
  *   produced in time instant t - 1, if any. If "DeltaRampDown" has length 1
  *   then DM[ t ] contains the same value for all t. Otherwise,
  *   DeltaRampDown[ i ] is the fixed value of DM[ t ] for all t in the
  *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional;
  *   if it is not provided then it is assumed that DM[ t ] == MaxP[ t ],
  *   i.e., the unit can ramp down an arbitrary amount, i.e., there are no
  *   ramp-down constraints. If NumberIntervals <= 1 or NumberIntervals >=
  *   TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "StoringBatteryRho", of type netCDF::NcDouble and to be
  *   either of size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   SBR[ t ] that, for each time instant t, contains the inefficiency of
  *   storing energy of the unit for the corresponding time step. This
  *   variable is optional; if it is not provided then it is assumed that
  *   SBR[ t ] == 1 for all t, i.e., no (significant) energy is spent just
  *   for storing it in the battery (this simplifies the model somewhat, see
  *   below). If "StoringBatteryRho" has length 1 then SBR[ t ] contains the
  *   same value for all t. Otherwise, StoringBatteryRho[ i ] is the fixed
  *   value of SBR[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ] with the assumption that ChangeIntervals[ - 1 ] =
  *   0. Note that it must be always such that SBR[ t ] <= 1 for all t (as
  *   SBR[ t ] is the amount of energy actually going in the battery for each
  *   1 unit of input energy). If NumberIntervals <= 1 or NumberIntervals >=
  *   TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "ExtractingBatterRho", of type netCDF::NcDouble and to be
  *   either of size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   EBR[ t ] that, for each time instant t, contains the inefficiency of
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
  * - The scalar variable "InitialStorage", of type netCDF::NcDouble and not
  *   indexed over any dimension. This variable indicates the amount of the
  *   storage level that the unit was producing at time instant -1, i.e.,
  *   before the start of the time horizon; this is necessary to compute the
  *   storage level connection with intake and outtake constraints. If
  *   InitialStorage <= 0 then the cyclical notation is used, so the
  *   constraint v_storage_level[ 0 ] = v_storage_level[ t - 1 ] is added
  *   to handle the unknown storage level of the battery at time zero; but
  *   since negative values for this data does not make sense in domain
  *   since the batteries cannot have a negative storage level obviously, it is
  *   used as if it was 0.
  *
  * - The variable "Cost", of type netCDF::NcDouble and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector C[ t ] that, for
  *   each time instant t, contains the monetary cost of storing one unit or
  *   energy into, or extracting it from, the battery (the cost is the same
  *   in both cases) at the corresponding time step. This variable is
  *   optional; if it is not provided then it's taken to be zero. If "Cost"
  *   has length 1 then C[ t ] contains the same value for all t. Otherwise,
  *   Cost[ i ] is the fixed value of C[ t ] for all t in the interval
  *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "Demand", of type netCDF::NcDouble and indexed over the
  *   dimension "TimeHorizon": the entry Demand[ t ] is assumed to contain
  *   the amount of energy that must be discharged from the battery and "sent
  *   away for some other purpose" (say, driving your e-car) at time t. This
  *   variable is optional; if it isn't defined, then Demand[ t ] == 0.
  *   Otherwise the Demand[ t ] contains the demand value for each time
  *   instant t.
  *
  * - The scalar variable "Kappa", of type netCDF::NcDouble. This variable
  *   contains the factor that multiplies the minimum and maximum active
  *   power, maximum primary and secondary reserve, and the minimum and
  *   maximum storage levels, at each time instant t. This variable is
  *   optional, if it is not provided it is taken to be Kappa == 1. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the BatteryUnitBlock
 /** This function generates the static variables of The BatteryUnitBlock,
  * which are:
  *
  * - The primary spinning reserve variables.
  *
  * - The secondary spinning reserve variables.
  *
  * - The active power variables, which can be positive or negative. If it is
  *   positive, the unit is giving energy to the system. If it is negative,
  *   it is taking energy away and adding to the storage. Since the storing
  *   and extracting amount of active power are not always equal, to deal
  *   with this issue, the usual trick of splitting the active power variable
  *   in two new non-negative variables which are called intake and outtake
  *   levels for each time t (see equation (5)) is used. If
  *   "StoringBatteryRho" == "ExtractingBatterRho" == 1, we do not need to
  *   split the active power and the constraints ((5)-(7) and (10)-(11)) will
  *   be replaced by (8)).
  *
  * - The storage level variables.
  *
  * - The intake and outtake levels variable. They are needed to split the
  *   active power variable (if it is needed).
  *
  * - The binary variables. When "StoringBatteryRho" == "ExtractingBatterRho"
  *   == 1, then this binary variable and all constraints which are depended
  *   on this variable are not required to be define.
  *
  * Each of these groups of variables either has size #f_time_horizon or is
  * empty (in case the variables have not been generated).
  *
  * The primary and secondary spinning reserve and the binary variables are
  * optional:
  *
  * - The primary spinning reserve variables are generated only if they were
  *   instructed to be (see set_reserve_vars()) and "MaxPrimaryPower" is not
  *   zero.
  *
  * - The secondary spinning reserve variables are generated only if they
  *   were instructed to be (see set_reserve_vars()) and "MaxSecondaryPower"
  *   is not zero.
  *
  * - The binary variables are generated only if negative prices may occur
  *   (which can be informed via a Configuration; see below) and there exists
  *   t such that StoringBatteryRho[ t ] < 1 and ExtractingBatterRho[ t ] > 1.
  *
  * Notice that despite the BatteryUnitBlock being a single logical unit, in
  * fact a battery is made up of the battery itself responsible for the
  * energy storage, and the converter responsible for the intake and outtake
  * of the energy from the battery. For this reason, in the design scenario
  * of the UC problem, i.e., if an investment cost is given for both the
  * battery and the converter, additional two binary variables are needed in
  * order to let the model infer how much capacity to install of either.
  *
  * The parameter \p stvv and the Configuration for this function presented in
  * the BlockConfig (namely, #f_BlockConfig->f_static_variables_Configuration)
  * can be used to indicate whether negative prices may occur. The parameter
  * \p stvv has priority over the BlockConfig in the sense that the
  * Configuration in the BlockConfig is only considered if no valid
  * Configuration has been provided in \p stvv. By default, it is assumed that
  * negative prices do not occur and, therefore, the binary variables are not
  * generated. If the Configuration is a SimpleConfiguration< int >, then a
  * nonzero value stored in this Configuration indicates that negative prices
  * may occur. The value zero indicates that negative prices do not occur. */

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
  *   std::vector< FRowConstraint >; with the dimension of f_time_horizon, where
  *   the entry t = 0, ..., f_time_horizon - 1 being the maximum and minimum
  *   power output value according to the primary and the secondary spinning
  *   reserves at time t. these ensure the maximum (or minimum) amount of
  *   energy that unit can produce (or use) when it is on (or off).
  *
  *   \f[
  *      p^{ac}_t + p^{pr}_t + p^{sc}_t \leq \kappa P^{mx,b}_t
  *                                           \quad t \in \mathcal{T} \quad (1)
  *   \f]
  *
  *   \f[
  *     \kappa P^{mn,b}_t \leq p^{ac}_t - p^{pr}_t - p^{sc}_t
  *                                           \quad t \in \mathcal{T} \quad (2)
  *   \f]
  *
  *   where \f$ P^{mx,b}_t \f$ and \f$ P^{mn,b}_t \f$ are the maximum and
  *   minimum power output parameters for each time t of the time horizon \f$
  *   \mathcal{T} \f$ respectively. In the design scenario of the UC problem,
  *   they become respectively:
  *
  *   \f[
  *      p^{ac}_t + p^{pr}_t + p^{sc}_t \leq x_b ( \kappa P^{mx,b}_t )
  *                                           \quad t \in \mathcal{T} \quad (1)
  *   \f]
  *
  *   \f[
  *     x_b ( \kappa P^{mn,b}_t ) \leq p^{ac}_t - p^{pr}_t - p^{sc}_t
  *                                           \quad t \in \mathcal{T} \quad (2)
  *   \f]
  *
  *   where \f$ x_b \f$ is the design variable of the battery.
  *
  * - ramp-up and ramp-down constraints are presented in (3)-(4). Each of them
  *   is a std::vector< FRowConstraint >; with the dimension of f_time_horizon,
  *   where the entry t = 0, ..., f_time_horizon - 1 being the ramp up and ramp
  *   down constraints which are presented as:
  *
  *   \f[
  *    p^{ac}_t - p^{ac}_{t-1} \leq \Delta^{up}_t
  *                                           \quad t \in \mathcal{T} \quad (3)
  *   \f]
  *
  *   \f[
  *    p^{ac}_t - p^{ac}_{t-1} \geq - \Delta^{dn}_t
  *                                           \quad t \in \mathcal{T} \quad (4)
  *   \f]
  *
  *   where \f$ \Delta^{up}_t \f$ and \f$ \Delta^{dn}_t \f$ are the delta
  *   ramp-up and delta ramp down threshold for each time t of the time
  *   horizon \f$ \mathcal{T} \f$ respectively.
  *
  * - active power relation with intake and outtake levels constraints are
  *   presented in (5). Each of them is a std::vector< FRowConstraint >; with
  *   the dimension of f_time_horizon, where the entry t =
  *   0, ..., f_time_horizon - 1 being the active power relation with intake and
  *   outtake levels at time t. These ensure the active power at each time
  *   should be equal to the intake and outtake difference, i.e.:
  *
  *   \f[
  *    p^{ac}_t = p^+_t - p^-_t
  *                                           \quad t \in \mathcal{T} \quad (5)
  *   \f]
  *
  *   The equations (6.1-6.2) also indicates the upper bound of intake and
  *   outtake level at each time instant t:
  *
  *   \f[
  *    p^+_t \leq \kappa C^- P^{mn,b}_t
  *                                         \quad t \in \mathcal{T} \quad (6.1)
  *   \f]
  *
  *   \f[
  *     p^-_t \leq \kappa C^+ P^{mx,b}_t
  *                                         \quad t \in \mathcal{T} \quad (6.2)
  *   \f]
  *
  *   where \f$ C^+ \f$ and \f$ C^- \f$ are the C-rate of the battery in
  *   charge and discharge respectively and which, in the design scenario of
  *   the UC problem, become:
  *
  *   \f[
  *     p^+_t \leq \kappa x_b ( C^- P^{mn,b}_t )
  *                                         \quad t \in \mathcal{T} \quad (6.3)
  *   \f]
  *
  *   \f[
  *     p^-_t \leq \kappa x_b ( C^+ P^{mx,b}_t )
  *                                         \quad t \in \mathcal{T} \quad (6.4)
  *   \f]
  *
  *   \f[
  *     p^+_t + p^-_t \leq x_c ( \kappa P^{mx,c}_t )
  *                                         \quad t \in \mathcal{T} \quad (6.5)
  *   \f]
  *
  *   where \f$ P^{mx,c}_t \f$ is the maximum power of the converter, and
  *   \f$ x_b \f$ and \f$ x_c \f$ are the battery and converter design variable
  *   respectively.
  *
  * - storage level relation with intake and outtake levels (if any)
  *   constraints in Battery unit are presented in (7). That is a
  *   std::vector< FRowConstraint >; with the dimension of f_time_horizon, where
  *   the entry t = 0, ..., f_time_horizon - 1 being the storage level relation
  *   with intake and outtake levels at time t.
  *
  *   \f[
  *    v^{ba}_t = v^{ba}_{t-1} - \rho^+_tp^+_t +
  *    \rho^-_tp^-_t-d^{ba}_t
  *                                           \quad t \in \mathcal{T} \quad (7)
  *   \f]
  *
  *   Note that if the equation (7) will change as below which is a
  *   std::vector< FRowConstraint >; with the dimension of f_time_horizon, where
  *   the entry t = 0, ..., f_time_horizon - 1 being the storage level relation
  *   with battery demand (if any) at time t.
  *
  *   \f[
  *    v^{ba}_t = v^{ba}_{t-1} - p^{ac}_t - d^{ba}_t
  *                                           \quad t \in \mathcal{T} \quad (8)
  *   \f]
  *
  *   The equation (9a) gives the storage levels upper bound and lower bound at
  *   each time instant t.
  *
  *   \f[
  *    v^{ba}_t \in [ \kappa V^{mn}_t , \kappa V^{mx}_t]
  *                                          \quad t \in \mathcal{T} \quad (9a)
  *   \f]
  *
  *   which, in the design scenario of the UC problem, become:
  *
  *   \f[
  *    x_b ( \kappa V^{mn}_t ) \leq v^{ba}_t \leq x_b ( \kappa V^{mx}_t )
  *                                          \quad t \in \mathcal{T} \quad (9b)
  *   \f]
  *
  *   where \f$ \rho^+_t \f$ and \f$ \rho^-_t \f$ are the
  *   ExtractingBatteryRho and StoringBatteryRho, \f$ V^{mn}_t \f$ and
  *   \f$ V^{mx}_t \f$ are the minimum and maximum storage level for each time
  *   t of the time horizon \f$ \mathcal{T} \f$ respectively; and \f$ x_b \f$
  *   is the design variable of the battery.
  *
  * - binary variable relation with storing and extracting energy level (if
  *   any) constraints are presented in (10-11). Each of them is a
  *   std::vector< FRowConstraint >; with the dimension of f_time_horizon, where
  *   the entry t = 0, ..., f_time_horizon - 1 being the binary variable
  *   relation with storing and extracting energy levels at time t.
  *
  *   \f[
  *    p^+_t \leq u^+_t P^{mx,b}_t
  *                                          \quad t \in \mathcal{T} \quad (10)
  *   \f]
  *
  *   \f[
  *    p^-_t \leq -(1 - u^+_t) P^{mn,b}_t
  *                                          \quad t \in \mathcal{T} \quad (11)
  *   \f]
  *
  *   Note that when \f$ \rho^+_t = \rho^-_t = 1 \f$, the binary variable \f$
  *   u^+_t \f$ is not required and neither are the last two constraints
  *   (10-11).
  *
  * - primary and secondary reserve upper bounds (if any) are presented by the
  *   equations (12-13). Each of them is a std::vector< FRowConstraint >; with
  *   the dimension of f_time_horizon, where the entry
  *   t = 0, ..., f_time_horizon - 1 being the primary and secondary reserve
  *   upper bounds at time t.
  *
  *   \f[
  *     p^{pr}_t \leq P^{mx, pr}_t
  *                                          \quad t \in \mathcal{T} \quad (12)
  *   \f]
  *
  *   \f[
  *     p^{sc}_t \leq P^{mx, sc}_t
  *                                          \quad t \in \mathcal{T} \quad (13)
  *   \f]
  */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the BatteryUnitBlock
 /** Method that generates the objective of the BatteryUnitBlock.
  *
  * - Objective function: the objective function of the BatteryUnitBlock
  *   is given as follow:
  *
  *   \f[
  *     \min ( I_b x_b + I_c x_c +
  *     \sum_{ t \in \mathcal{T} } C_t (p^+_t +  p^-_t) )
  *   \f]
  *
  *   where \f$ I_b \f$ and \f$ I_c \f$ are the the investment costs of the
  *   battery and the converter respectively, \f$ x_b \f$ and \f$ x_c \f$ are
  *   the battery and converter design variable respectively, and \f$ C_t \f$
  *   is a certain proportion cost function. */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------- Methods for checking the BatteryUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the BatteryUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this BatteryUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this BatteryUnitBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this BatteryUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or #f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  * - If \p fsbc is not a nullptr and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr and it is a
  *   SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *   fsbc->f_value.first and the type of violation is determined by
  *   fsbc->f_value.second (any nonzero number for relative violation and
  *   zero for absolute violation);
  *
  * - Otherwise, if both #f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of the
  *   parameters are obtained analogously as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation
  *   is considered.
  *
  * This function currently considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * currently ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed (Constraint::compute()).
  *
  * @param useabstract This parameter is currently ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE BatteryUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the BatteryUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of battery storage units
 * @{ */

 /// returns the initial storage value
 double get_initial_storage( void ) const { return( f_InitialStorage ); }

 /// returns the initial power value
 double get_initial_power( void ) const { return( f_InitialPower ); }

 /// returns the maximum C-rate of the battery in charge
 double get_max_C_rate_charge( void ) const { return( f_MaxCRateCharge ); }

 /// returns the maximum C-rate of the battery in discharge
 double get_max_C_rate_discharge( void ) const { return( f_MaxCRateDischarge ); }

 /// returns the battery investment cost
 double get_batt_investment_cost( void ) const {
  return( f_BattInvestmentCost );
 }

 /// returns the converter investment cost
 double get_conv_investment_cost( void ) const {
  return( f_ConvInvestmentCost );
 }

 /// returns the maximum battery installable capacity by the user
 double get_batt_max_capacity( void ) const { return( f_BattMaxCapacity ); }

 /// returns the maximum converter installable capacity by the user
 double get_conv_max_capacity( void ) const { return( f_ConvMaxCapacity ); }

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

 const std::vector< double > & get_min_storage( void ) const {
  return( v_MinStorage );
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

 const std::vector< double > & get_max_storage( void ) const {
  return( v_MaxStorage );
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
  * - otherwise, the vector V must have size get_time_horizon() and each
  * V[ t ] represents the minimum power value at time t. */

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( v_MinPower[ t ] );
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
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the maximum power value at time t. */

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower[ t ] );
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

 const std::vector< double > & get_max_primary_power( void ) const {
  return( v_MaxPrimaryPower );
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

 const std::vector< double > & get_max_secondary_power( void ) const {
  return( v_MaxSecondaryPower );
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
  return( v_DeltaRampUp );
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
  return( v_DeltaRampDown );
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
  return( v_StoringBatteryRho );
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
  return( v_ExtractingBatteryRho );
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
  return( v_Cost );
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
  return( v_Demand );
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

 /// returns the kappa factor
 /** This function returns the kappa factor, which multiplies the minimum and
  * maximum active power, maximum primary and secondary reserve, and the
  * minimum and maximum storage levels.
  *
  * @return The kappa factor. */

 double get_kappa( void ) const { return( f_kappa ); }

/*--------------------------------------------------------------------------*/
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
  *   outtake level variable for time step t. */

 const std::vector< ColVariable > & get_outtake_level( void ) const {
  return( v_outtake_level );
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
 /// returns the battery design variable

 ColVariable & get_batt_design( void ) {
  return( batt_design );
 }

/*--------------------------------------------------------------------------*/
 /// returns the converter design variable

 ColVariable & get_conv_design( void ) {
  return( conv_design );
 }

/*--------------------------------------------------------------------------*/
 /// returns the intake/outtake binary variables

 const std::vector< ColVariable > &
 get_intake_outtake_binary_variables( void ) const {
  return( v_battery_binary );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power output constraints

 const FRowConstraint * get_min_power_constraints( void ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 0 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const.data()[ 0 ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power output constraint associated with time t

 const FRowConstraint * get_min_power_constraint( Index t ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 0 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const[ 0 ][ t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power output constraints

 const FRowConstraint * get_max_power_constraints( void ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 1 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const.data()[ 1 ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power output constraint associated with time t

 const FRowConstraint * get_max_power_constraint( Index t ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 1 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const[ 1 ][ t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraints with binary variables

 const FRowConstraint * get_max_intake_binary_constraints( void ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 0 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const.data()[ 0 ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraint with binary variables for time t

 const FRowConstraint * get_max_intake_binary_constraint( Index t ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 0 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const[ 0 ][ t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraints with binary variables

 const FRowConstraint * get_max_outtake_binary_constraints( void ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 1 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const.data()[ 1 ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraint with binary variables for time t

 const FRowConstraint * get_max_outtake_binary_constraints( Index t ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 1 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const[ 1 ][ t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the storage level bound constraints

 const std::vector< BoxConstraint > & get_storage_level_bounds( void ) const {
  return( storage_level_bounds_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraints

 const LB0Constraint * get_max_intake_bounds( void ) const {
  if( intake_outtake_bounds_Const.empty() ||
      intake_outtake_bounds_Const[ 0 ].empty() )
   return( nullptr );
  return( &( intake_outtake_bounds_Const.data()[ 0 ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraint associated with time t

 const LB0Constraint * get_max_intake_bound( Index t ) const {
  if( intake_outtake_bounds_Const.empty() ||
      intake_outtake_bounds_Const[ 0 ].empty() )
    return( nullptr );
  return( & ( intake_outtake_bounds_Const[ 0 ][ t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraints

 const LB0Constraint * get_max_outtake_bounds( void ) const {
  if( intake_outtake_bounds_Const.empty() ||
      intake_outtake_bounds_Const[ 1 ].empty() )
   return( nullptr );
  return( &( intake_outtake_bounds_Const.data()[ 1 ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraint associated with time t

 const LB0Constraint * get_max_outtake_bound( Index t ) const {
  if( intake_outtake_bounds_Const.empty() ||
      intake_outtake_bounds_Const[ 1 ].empty() )
    return( nullptr );
  return( & ( intake_outtake_bounds_Const[ 1 ][ t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the primary reserve bound constraints

 const std::vector< LB0Constraint > & get_primary_reserve_bounds( void ) const {
  return( primary_upper_bound_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the secondary reserve bound constraints

 const std::vector< LB0Constraint > & get_secondary_reserve_bounds( void ) const {
  return( secondary_upper_bound_Const );
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
/** @name Methods for changing the data of the BatteryUnitBlock
 *  @{ */

 void set_initial_storage( MF_dbl_it it ,
                           Subset && subset ,
                           const bool ordered = false ,
                           c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck );

 void set_initial_storage( MF_dbl_it it ,
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

 void set_initial_power( MF_dbl_it it ,
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

 void set_initial_power( MF_dbl_it it ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the scale factor of this BatteryUnitBlock
 /** This method sets the scale factor of this BatteryUnitBlock.
  *
  * @param values An iterator to a vector containing the scale factor.
  *
  * @param subset If non-empty, the scale factor is set to the value pointed
  *               by \p values. If empty, no operation is performed.
  *
  * @param ordered This parameter is ignored.
  *
  * @param issuePMod Controls how physical Modification are issued.
  *
  * @param issueAMod Controls how abstract Modification are issued. */

 void scale( MF_dbl_it values ,
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

 void set_kappa( MF_dbl_it values ,
                 Subset && subset ,
                 const bool ordered = false ,
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

 void set_kappa( MF_dbl_it values ,
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

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the vector of minimum storage
 std::vector< double > v_MinStorage;

 /// the vector of maximum storage
 std::vector< double > v_MaxStorage;

 /// the vector of MinPower
 std::vector< double > v_MinPower;

 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the vector of ConverterMaxPower
 std::vector< double > v_ConvMaxPower;

 /// the vector of MaxPrimaryPower
 std::vector< double > v_MaxPrimaryPower;

 /// the vector of MaxSecondaryPower
 std::vector< double > v_MaxSecondaryPower;

 /// the vector of RampUp
 std::vector< double > v_DeltaRampUp;

 /// the vector of RampDown
 std::vector< double > v_DeltaRampDown;

 /// the vector of StoringBatteryRho
 std::vector< double > v_StoringBatteryRho;

 /// the vector of ExtractingBatteryRho
 std::vector< double > v_ExtractingBatteryRho;

 /// the vector of Cost
 std::vector< double > v_Cost;

 /// the vector of demand
 std::vector< double > v_Demand;


 /// the battery investment cost
 double f_BattInvestmentCost{};

 /// the converter investment cost
 double f_ConvInvestmentCost{};

 /// the maximum battery installable capacity by the user
 double f_BattMaxCapacity{};

 /// the maximum converter installable capacity by the user
 double f_ConvMaxCapacity{};

 /// the InitialStorage value
 double f_InitialStorage{};

 /// the InitialPower value
 double f_InitialPower{};

 /// the MaxCRateCharge value
 double f_MaxCRateCharge = 1;

 /// the MaxCRateDischarge value
 double f_MaxCRateDischarge = 1;

 /// the kappa value
 double f_kappa = 1;

 /// the scale factor
 double f_scale = 1;

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

 /// the battery design variable
 ColVariable batt_design;

 /// the converter design variable
 ColVariable conv_design;

/*------------------------------- constraints ------------------------------*/

 /// the active power bounds constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_Const;

 /// the active power bounds design constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_design_Const;

 /// the intake outtake upper bounds design constraints
 boost::multi_array< FRowConstraint , 2 > intake_outtake_upper_bounds_design_Const;

 /// the storage level bounds design constraints
 boost::multi_array< FRowConstraint , 2 > storage_level_bounds_design_Const;

 /// the intake and outtake binary variable relation constraints
 boost::multi_array< FRowConstraint , 2 > intake_outtake_binary_Const;

 /// the active power, intake and outtake relation constraints
 std::vector< FRowConstraint > power_intake_outtake_Const;

 /// the ramp up constraints
 std::vector< FRowConstraint > ramp_up_Const;

 /// the ramp down constraints
 std::vector< FRowConstraint > ramp_down_Const;

 /// the demand constraints
 std::vector< FRowConstraint > demand_Const;


 /// the storage level bound constraints
 std::vector< BoxConstraint > storage_level_bounds_Const;


 /// the intake and outtake bounds constraints
 boost::multi_array< LB0Constraint , 2 > intake_outtake_bounds_Const;

 /// primary upper bound constraints
 std::vector< LB0Constraint > primary_upper_bound_Const;

 /// secondary upper bound constraints
 std::vector< LB0Constraint > secondary_upper_bound_Const;


 /// the vector of binary bound constraints
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

 void update_initial_storage_in_cnstrs( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current initial power
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints at time 0 (which are the
  * constraints that depend on the initial power). */

 void update_initial_power_in_cnstrs( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current kappa
 /** This function updates the constraints to take into account the current
  * value of the kappa constant. */

 void update_kappa_in_cnstrs( ModParam issueAMod = eNoBlck );

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
  * If any of the above conditions are not met, an exception is thrown. */

 void check_data_consistency( void ) const;

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

};  // end( class( BatteryUnitBlock ) )

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

 Block::Range f_rng;  ///< the range

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

 Block::Subset f_nms;  ///< the subset

};  // end( class( BatteryUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* BatteryUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File BatteryUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
