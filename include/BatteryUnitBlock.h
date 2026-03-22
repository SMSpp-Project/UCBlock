/*--------------------------------------------------------------------------*/
/*------------------------- File BatteryUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BatteryUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" Battery
 * storage, E-mobility, Centralized demand response, Distributed load
 * management, Distributed storage, and Power-to-gas units in a single class
 * in the Unit Commitment problem.
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
 * equation pattern. For instance, a BatteryUnitBlock may or may not have a
 * fixed demand (e-mobility has, other units have not) and it may or may not
 * provide primary and secondary reserve (battery storage may do, but other
 * units don't).
 *
 * Battery storage provides additional flexibility to the system by shifting
 * a surplus of electric energy (e.g., due to high renewable feeding) to times
 * with high demand or lower renewable generation. Distributed battery storage
 * can be aggregated in the energy cells or directly placed in a single node
 * of the network. We will therefore not stress this dependency in the
 * subsequent equations. We emphasize that the potential contribution of
 * batteries to inertia is still a subject of active research and should be
 * considered optional. Besides, since the transport sector is moving towards
 * electrification, electric mobility will have a rising impact on the
 * electricity system. First, electricity demand is growing due to a higher
 * amount of electric vehicles that need to be charged. On the other hand,
 * vehicles are used only a small amount of time while being charged over a
 * much longer timespan (e.g., at night). This allows shifting the charging
 * process in time and providing flexibility to the overall energy system by
 * means of an additional generator (vehicle-to-grid) or an additional load
 * (power-to-vehicle). Two main differences between battery storage units and
 * other existing units in this class are:
 *
 * - battery storage units can do primary and secondary reserve, while other
 *   units cannot;
 *
 * - some of the units may have a fixed demand that battery storage units do
 *   not.
 *
 * Moreover, as the considered storage cycle is small w.r.t. the EUC time
 * horizon, distributed storage is not considered as seasonal storage. Hence,
 * the associated mathematical description follows the same equations as the
 * one provided for battery storage units. The specificity of distributed
 * storage only relies on the fact that it is connected to a distribution grid
 * node.
 *
 * To model BatteryUnitBlock systems several technical parameters have to be
 * considered. These are divided into the battery storage level parameters,
 * the ramping parameters, the active power bound parameters, and a flexible
 * electric demand that provides flexibility to the overall system while
 * accounting for storage level constraints. The technical and physical
 * constraints are mainly divided in several different categories:
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

 explicit BatteryUnitBlock( Block * f_block = nullptr ) :
  UnitBlock( f_block ), f_BattInvestmentCost( 0 ), f_ConvInvestmentCost( 0 ),
  f_BattMinCapacityDesign( 0 ), f_BattMaxCapacityDesign( 1 ),
  f_ConvMinCapacityDesign( 0 ), f_ConvMaxCapacityDesign( 1 ),
  f_InitialStorage( 0 ), f_InitialPower( 0 ), f_MaxCRateCharge( 1 ),
  f_MaxCRateDischarge( 1 ), f_kappa( 1 ), f_scale( 1 ) {}

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
  * - The scalar variable "BatteryInvestmentCost", of type netCDF::NcDouble
  *   and not indexed over any dimension. When provided and nonzero, a
  *   battery design variable \( x_b \) is created and contributes
  *   \( I_b x_b \) to the objective.
  *
  * - The scalar variable "ConverterInvestmentCost", of type netCDF::NcDouble
  *   and not indexed over any dimension. When provided and nonzero, a 
  *   converter design variable \( x_c \) is created and contributes
  *   \( I_c x_c \) to the objective.
  *
  * - The scalar variable "BatteryMaxCapacityDesign", of type netCDF::NcDouble
  *   and not indexed over any dimension. This limits the design variable
  *   \f$ x_b \f$:
  *   - if \f$ \mathrm{BatteryMaxCapacityDesign} < 0 \f$ then \f$ x_b \f$ is
  *     binary;
  *   - otherwise \f$ x_b \f$ is a nonnegative continuous variable with
  *     \f$ 0 \le x \le \mathrm{BatteryMaxCapacityDesign} \f$.
  *   If not provided, the default is 1.
  *
  * - The scalar variable "BatteryMinCapacityDesign", of type
  *   netCDF::NcDouble and not indexed over any dimension. This sets the
  *   lower bound of the battery design variable \( x_b \) in design mode.
  *   If not provided, the default is 0. Its meaning depends on
  *   "BatteryMaxCapacityDesign":
  *
  *   - if \( \mathrm{BatteryMaxCapacityDesign} < 0 \) (binary design), then
  *     \( x_b \in \{ 0 , 1 \} \) and
  *     \( \mathrm{BatteryMinCapacityDesign} > 0 \) implies \( x_b = 1 \);
  *
  *   - otherwise (continuous design), \( x_b \) is nonnegative continuous
  *     with
  *     \$[
  *       \mathrm{BatteryMinCapacityDesign} \le x_b \le
  *       \mathrm{BatteryMaxCapacityDesign}
  *     \$]
  *
  * - The scalar variable "ConverterMaxCapacityDesign", of type 
  *   netCDF::NcDouble and not indexed over any dimension. This limits the
  *   design variable \f$ x_c \f$:
  *
  *   - if \f$ \mathrm{ConverterMaxCapacityDesign} < 0 \f$ then \f$ x_c \f$
  *     is binary;
  *
  *   - otherwise \f$ x_c \f$ is a nonnegative continuous variable with
  *     \f$ 0 \le x \le \mathrm{ConverterMaxCapacityDesign} \f$.
  *
  *   If not provided, the default is 1.
  *
  * - The scalar variable "ConverterMinCapacityDesign", of type
  *   netCDF::NcDouble and not indexed over any dimension. This sets the
  *   lower bound of the converter design variable \( x_c \) in design mode.
  *   If not provided, the default is 0. Its meaning depends on
  *   "ConverterMaxCapacityDesign":
  *
  *   - if \( \mathrm{ConverterMaxCapacityDesign} < 0 \) (binary design),
  *     then \( x_c \in \{ 0 , 1 \} \) and
  *     \( \mathrm{ConverterMinCapacityDesign} > 0 \) implies \( x_c = 1 \);
  *
  *   - otherwise (continuous design), \( x_c \) is nonnegative continuous
  *     with
  *     \$[
  *       \mathrm{ConverterMinCapacityDesign} \le x_c \le
  *       \mathrm{ConverterMaxCapacityDesign}
  *     \$]
  *
  * - The scalar variable "BatteryMaxCapacity", of type netCDF::NcDouble and
  *   not indexed over any dimension. This is the maximum installable battery
  *   capacity chosen by the user.
  *
  * - The scalar variable "ConverterMaxCapacity", of type netCDF::NcDouble
  *   and not indexed over any dimension. This is the maximum installable
  *   converter capacity chosen by the user.
  *
  * - The variable "MinStorage", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector \f$ \mathrm{MinS}
  *   [ t ] \f$ that, for each time instant \f$ t \f$, contains the minimum
  *   storage level of the unit for the corresponding time step. If
  *   "MinStorage" has length 1 then \f$ \mathrm{MinS}[ t ] \f$ contains the
  *   same value for all \f$ t \f$. Otherwise, \f$ \mathrm{MinStorage}[ i ] \f$
  *   is the fixed value of \f$ \mathrm{MinS}[ t ] \f$ for all \f$ t \f$ in
  *   the interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. Note that it must always be
  *   \f$ 0 \le \mathrm{MinS}[ t ] < \mathrm{MaxS}[ t ] \f$ for all \f$ t \f$.
  *   If \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "MaxStorage", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector \f$ \mathrm{MaxS}
  *   [ t ] \f$ that, for each time instant \f$ t \f$, contains the maximum
  *   storage level of the unit for the corresponding time step. If
  *   "MaxStorage" has length 1 then \f$ \mathrm{MaxS}[ t ] \f$ contains the
  *   same value for all \f$ t \f$. Otherwise, \f$ \mathrm{MaxStorage}[ i ] \f$
  *   is the fixed value of \f$ \mathrm{MaxS}[ t ] \f$ for all \f$ t \f$ in
  *   the interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. Note that it must always be
  *   \f$ \mathrm{MinS}[ t ] < \mathrm{MaxS}[ t ] \f$ for all \f$ t \f$. If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "MinPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector \f$ \mathrm{MinP}
  *   [ t ] \f$ that, for each time instant \f$ t \f$, contains the minimum
  *   active power output value of the unit for the corresponding time step.
  *   If "MinPower" has length 1 then \f$ \mathrm{MinP}[ t ] \f$ contains the
  *   same value for all \f$ t \f$. Otherwise, \f$ \mathrm{MinPower}[ i ] \f$
  *   is the fixed value of \f$ \mathrm{MinP}[ t ] \f$ for all \f$ t \f$ in
  *   the interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. Note that it must be
  *   \f$ \mathrm{MinP}[ t ] \le 0 \f$ for all \f$ t \f$. If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded. If "MinPower" is not provided, then
  *   \f$ \mathrm{MinP}[ t ] = - \mathrm{MaxP}[ t ] \f$ for all \f$ t \f$.
  *
  * - The variable "MaxPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector \f$ \mathrm{MaxP}
  *   [ t ] \f$ that, for each time instant \f$ t \f$, contains the maximum
  *   active power output value of the unit for the corresponding time step.
  *   If "MaxPower" has length 1 then \f$ \mathrm{MaxP}[ t ] \f$ contains the
  *   same value for all \f$ t \f$. Otherwise, \f$ \mathrm{MaxPower}[ i ] \f$
  *   is the fixed value of \f$ \mathrm{MaxP}[ t ] \f$ for all \f$ t \f$ in
  *   the interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. Note that it must be
  *   \f$ \mathrm{MinP}[ t ] < \mathrm{MaxP}[ t ] \f$ for all \f$ t \f$. If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "ConverterMaxPower", of type netCDF::NcDouble and either
  *   of size 1 or indexed over "NumberIntervals" (if not provided, it can be
  *   indexed over "TimeHorizon"). It represents the vector
  *   \f$ P^{mx,c}_t \f$, i.e., the converter maximum power at time \f$ t \f$.
  *   If length is 1 the same value applies to all \f$ t \f$. Otherwise,
  *   \f$ \mathrm{ConverterMaxPower}[ i ] \f$ fixes \f$ P^{mx,c}_t \f$ for
  *   \f$ t \in [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$. This variable is optional; if
  *   it is not provided then it is taken to be \f$ P^{mx,c}_t = P^{mx,b}_t \f$
  *   (i.e., equal to "MaxPower" at each time).
  *
  * - The scalar variable "InitialPower", of type netCDF::NcDouble and not
  *   indexed over any dimension. This variable indicates the amount of power
  *   that the unit was producing at time instant -1, i.e., before the start
  *   of the time horizon; this is necessary to compute the ramp-up and
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
  *   \f$ \mathrm{MaxPP}[ t ] \f$ that, for each time instant \f$ t \f$,
  *   contains the maximum active power that can be used as primary reserve
  *   for the corresponding time step. If "MaxPrimaryPower" has length 1 then
  *   \f$ \mathrm{MaxPP}[ t ] \f$ contains the same value for all \f$ t \f$.
  *   Otherwise, \f$ \mathrm{MaxPrimaryPower}[ i ] \f$ is the fixed value of
  *   \f$ \mathrm{MaxPP}[ t ] \f$ for all \f$ t \f$ in the interval
  *   \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *         \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. This variable is optional;
  *   if it is not provided then \f$ \mathrm{MaxPP}[ t ] = 0 \f$ for all
  *   \f$ t \f$. If \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "MaxSecondaryPower", of type netCDF::NcDouble and either
  *   of size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   \f$ \mathrm{MaxSP}[ t ] \f$ that, for each time instant \f$ t \f$,
  *   contains the maximum active power that can be used as secondary reserve
  *   for the corresponding time step. If "MaxSecondaryPower" has length 1
  *   then \f$ \mathrm{MaxSP}[ t ] \f$ contains the same value for all
  *   \f$ t \f$. Otherwise, \f$ \mathrm{MaxSecondaryPower}[ i ] \f$ is the
  *   fixed value of \f$ \mathrm{MaxSP}[ t ] \f$ for all \f$ t \f$ in the
  *   interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. This variable is optional;
  *   if it is not provided then \f$ \mathrm{MaxSP}[ t ] = 0 \f$ for all
  *   \f$ t \f$. Note that \f$ \mathrm{MaxPP}[ t ] = 0 \f$ implies
  *   \f$ \mathrm{MaxSP}[ t ] = 0 \f$ (that is, if MaxPrimaryPower is not
  *   defined then neither should MaxSecondaryPower). If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "DeltaRampUp", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector \f$ \mathrm{DP}
  *   [ t ] \f$ that, for each time instant \f$ t \f$, contains the ramp-up
  *   value of the unit for the corresponding time step, i.e., the maximum
  *   possible increase of active power production w.r.t. the power that had
  *   been produced in time instant \f$ t - 1 \f$, if any. If "DeltaRampUp"
  *   has length 1 then \f$ \mathrm{DP}[ t ] \f$ contains the same value for
  *   all \f$ t \f$. Otherwise, \f$ \mathrm{DeltaRampUp}[ i ] \f$ is the fixed
  *   value of \f$ \mathrm{DP}[ t ] \f$ for all \f$ t \f$ in the interval
  *   \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *        \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that 
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. This variable is optional;
  *   if it is not provided then it is assumed that
  *   \f$ \mathrm{DP}[ t ] = \mathrm{MaxP}[ t ] \f$, i.e., the unit can ramp
  *   up by an arbitrary amount (no ramp-up constraints). If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "DeltaRampDown", of type netCDF::NcDouble and either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   \f$ \mathrm{DM}[ t ] \f$ that, for each time instant \f$ t \f$, contains
  *   the ramp-down value of the unit for the corresponding time step, i.e.,
  *   the maximum possible decrease of active power production w.r.t. the
  *   power that had been produced in time instant \f$ t - 1 \f$, if any. If
  *   "DeltaRampDown" has length 1 then \f$ \mathrm{DM}[ t ] \f$ contains the
  *   same value for all \f$ t \f$. Otherwise, \f$ \mathrm{DeltaRampDown}[ i ]
  *   \f$ is the fixed value of \f$ \mathrm{DM}[ t ] \f$ for all \f$ t \f$ in
  *   the interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. This variable is optional;
  *   if it is not provided then it is assumed that
  *   \f$ \mathrm{DM}[ t ] = \mathrm{MaxP}[ t ] \f$, i.e., the unit can ramp
  *   down an arbitrary amount (no ramp-down constraints). If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "StoringBatteryRho", of type netCDF::NcDouble and either
  *   of size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   \f$ \mathrm{SBR}[ t ] \f$ that, for each time instant \f$ t \f$,
  *   contains the inefficiency of storing energy of the unit for the
  *   corresponding time step. This variable is optional; if it is not
  *   provided then \f$ \mathrm{SBR}[ t ] = 1 \f$ for all \f$ t \f$, i.e., no
  *   (significant) energy is spent just for storing it in the battery (this
  *   simplifies the model somewhat, see below). If "StoringBatteryRho" has
  *   length 1 then \f$ \mathrm{SBR}[ t ] \f$ contains the same value for all
  *   \f$ t \f$. Otherwise, \f$ \mathrm{StoringBatteryRho}[ i ] \f$ is the
  *   fixed value of \f$ \mathrm{SBR}[ t ] \f$ for all \f$ t \f$ in the
  *   interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$ with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. Note that it must always
  *   be \f$ \mathrm{SBR}[ t ] \le 1 \f$ for all \f$ t \f$ (as \f$ \mathrm{SBR}
  *   [ t ] \f$ is the amount of energy actually going in the battery for each
  *   1 unit of input energy). If \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "ExtractingBatteryRho", of type netCDF::NcDouble and either
  *   of size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   \f$ \mathrm{EBR}[ t ] \f$ that, for each time instant \f$ t \f$,
  *   contains the inefficiency of extracting energy of the unit for the
  *   corresponding time step. This variable is optional; if it is not
  *   provided, then \f$ \mathrm{EBR}[ t ] = 1 \f$ for all \f$ t \f$, i.e., no
  *   (significant) energy is spent just for extracting it from the battery
  *   (this simplifies the model somewhat, see below). If
  *   "ExtractingBatteryRho" has length 1 then \f$ \mathrm{EBR}[ t ] \f$
  *   contains the same value for all \f$ t \f$. Otherwise,
  *   \f$ \mathrm{ExtractingBatteryRho}[ i ] \f$ is the fixed value of
  *   \f$ \mathrm{EBR}[ t ] \f$ for all \f$ t \f$ in the interval
  *   \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] , \mathrm{ChangeIntervals}
  *   [ i ] ] \f$ with the assumption that \f$ \mathrm{ChangeIntervals}[ - 1 ]
  *   = 0 \f$. Note that it must always be \f$ \mathrm{EBR}[ t ] \ge 1 \f$
  *   [\f$ \ge \mathrm{SBR}[ t ] \f$] for all \f$ t \f$ (as \f$ \mathrm{EBR}
  *   [ t ] \f$ is the amount of energy that is taken away from the battery to
  *   obtain 1 unit of output energy). If \f$ \mathrm{NumberIntervals} \le 1
  *   \f$ or \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then
  *   the mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded.
  *
  * Note: the special case in which \f$ \mathrm{EBR}[ t ] = \mathrm{SBR}[ t ]
  * = 1 \f$ for all \f$ t \f$, i.e., no energy is spent for storing it in /
  * retrieving it from the battery, leads to significantly simpler
  * mathematical models. In particular, one single variable can be used to
  * represent both storing and retrieving, rather than requiring two separate
  * ones (unless primary and secondary reserve are allowed and/or the cost is
  * defined, since this also requires using two), and the binary variables
  * need not be defined. For details, see the comments to
  * generate_abstract_variables() and generate_abstract_constraints().
  *
  * - The scalar variable "InitialStorage", of type netCDF::NcDouble and not
  *   indexed over any dimension. This variable indicates the amount of
  *   storage level that the unit had at time instant -1, i.e., before the
  *   start of the time horizon; this is necessary to compute the storage
  *   level connection with intake and outtake constraints. If
  *   \f$ \mathrm{InitialStorage} < 0 \f$ then the cyclical notation is used,
  *   so the constraint \f$ v_{\mathrm{storage\_level}}[ 0 ] =
  *   v_{\mathrm{storage\_level}}[ \mathrm{t} - 1 ] \f$ is added to handle the
  *   unknown storage level of the battery at time zero; but since negative
  *   values for this datum do not make sense (batteries cannot have a
  *   negative storage level), it is used as if it were 0.
  *
  * - The variable "Cost", of type netCDF::NcDouble and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector \f$ C[ t ] \f$
  *   that, for each time instant \f$ t \f$, contains the monetary cost of
  *   storing one unit of energy into, or extracting it from, the battery (the
  *   cost is the same in both cases) at the corresponding time step. This
  *   variable is optional; if it is not provided then it's taken to be zero.
  *   If "Cost" has length 1 then \f$ C[ t ] \f$ contains the same value for
  *   all \f$ t \f$. Otherwise, \f$ \mathrm{Cost}[ i ] \f$ is the fixed value
  *   of \f$ C[ t ] \f$ for all \f$ t \f$ in the interval
  *   \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] , \mathrm{ChangeIntervals}
  *   [ i ] ] \f$, with the assumption that \f$ \mathrm{ChangeIntervals}[ - 1 ]
  *   = 0 \f$. If \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is not
  *   loaded.
  *
  * - The variable "Demand", of type netCDF::NcDouble and indexed over the
  *   dimension "TimeHorizon": the entry \f$ \mathrm{Demand}[ t ] \f$ is
  *   assumed to contain the amount of energy that must be discharged from the
  *   battery and "sent away for some other purpose" (say, driving your e-car)
  *   at time \f$ t \f$. This variable is optional; if it isn't defined, then
  *   \f$ \mathrm{Demand}[ t ] = 0 \f$. Otherwise, \f$ \mathrm{Demand}[ t ] \f$
  *   contains the demand value for each time instant \f$ t \f$.
  *
  * - The scalar variable "Kappa", of type netCDF::NcDouble. This variable
  *   contains the factor that multiplies the minimum and maximum active
  *   power, maximum primary and secondary reserve, and the minimum and
  *   maximum storage levels, at each time instant \f$ t \f$. This variable is
  *   optional; if it is not provided it is taken to be \f$ \kappa = 1 \f$. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends UnitBlock::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends UnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the BatteryUnitBlock
 /** This function generates the static variables of the BatteryUnitBlock,
  * which are:
  *
  * - The primary spinning reserve variables.
  *
  * - The secondary spinning reserve variables.
  *
  * - The active power variables, which can be positive or negative. If it is
  *   positive, the unit is giving energy to the system. If it is negative, it
  *   is taking energy away and adding to the storage. Since the storing and
  *   extracting amount of active power are not always equal, to deal with
  *   this the usual trick of splitting the active power variable in two new
  *   non-negative variables called intake and outtake levels for each time
  *   \f$ t \f$ (see equation (5)) is used. If
  *   "StoringBatteryRho" == "ExtractingBatteryRho" == 1, we do not need to
  *   split the active power and constraints ((5)–(7) and (10)–(11)) are
  *   replaced by (8).
  *
  * - The storage level variables.
  *
  * - The intake and outtake level variables. They are needed to split the
  *   active power variable (if needed).
  *
  * - The binary variables. When "StoringBatteryRho" == "ExtractingBatteryRho"
  *   == 1, then this binary variable and all constraints that depend on it
  *   are not required.
  *
  * Each of these groups of variables either has size #f_time_horizon or is
  * empty (if the variables have not been generated).
  *
  * The primary and secondary spinning reserve and the binary variables are
  * optional:
  *
  * - The primary spinning reserve variables are generated only if they were
  *   instructed to be (see set_reserve_vars()) and "MaxPrimaryPower" is not
  *   zero.
  *
  * - The secondary spinning reserve variables are generated only if they were
  *   instructed to be (see set_reserve_vars()) and "MaxSecondaryPower" is not
  *   zero.
  *
  * - The binary variables are generated only if negative prices may occur
  *   (which can be informed via a Configuration; see below) and there exists
  *   \f$ t \f$ such that \f$ \mathrm{SBR}[ t ] < 1 \f$ and
  *   \f$ \mathrm{EBR}[ t ] > 1 \f$.
  *
  * Notice that despite the BatteryUnitBlock being a single logical unit, in
  * fact a battery is made up of the battery itself responsible for the energy
  * storage, and the converter responsible for the intake and outtake of the
  * energy from the battery. For this reason, in the design scenario of the UC
  * problem, i.e., if an investment cost is given for both the battery and the
  * converter, two additional design variables are needed in order to let the
  * model infer how much capacity to install of either. Denote them by
  * \f$ x_b \f$ (battery) and \f$ x_c \f$ (converter). Their **domains are
  * controlled** by the scalar parameters
  * \f$ \mathrm{BatteryMaxCapacityDesign} \f$ and
  * \f$ \mathrm{ConverterMaxCapacityDesign} \f$, respectively, as follows:
  *
  * - if \f$ \mathrm{BatteryMaxCapacityDesign} < 0 \f$ then \f$ x_b \f$ is
  *   **binary** (\f$ x_b \in \{0,1\} \f$); otherwise \f$ x_b \f$ is a
  *   **nonnegative continuous** variable with
  *   \f$ 0 \le x_b \le \mathrm{BatteryMaxCapacityDesign} \f$;
  *
  * - if \f$ \mathrm{ConverterMaxCapacityDesign} < 0 \f$ then \f$ x_c \f$ is
  *   **binary** (\f$ x_c \in \{0,1\} \f$); otherwise \f$ x_c \f$ is a
  *   **nonnegative continuous** variable with
  *   \f$ 0 \le x_c \le \mathrm{ConverterMaxCapacityDesign} \f$.
  *
  * This mirrors the behavior used in IntermittentUnitBlock so that the
  * "design mode" is consistent across unit types.
  *
  * In the design scenario (i.e., when a nonzero investment cost is provided),
  * additional design variables are created:
  * - \( x_b \) (battery), if "BatteryInvestmentCost" ≠ 0;
  * - \( x_c \) (converter), if "ConverterInvestmentCost" ≠ 0.
  *
  * Their domains are controlled by the corresponding *MaxCapacityDesign* and
  * *MinCapacityDesign* parameters, analogously to IntermittentUnitBlock:
  *
  * - if \( \mathrm{BatteryMaxCapacityDesign} < 0 \) then \( x_b \in \{0,1\} \);
  *   if moreover \( \mathrm{BatteryMinCapacityDesign} > 0 \) then \( x_b = 1 \).
  *   Otherwise, \( \mathrm{BatteryMinCapacityDesign} \le x_b \le \mathrm{BatteryMaxCapacityDesign} \).
  *
  * - if \( \mathrm{ConverterMaxCapacityDesign} < 0 \) then \( x_c \in \{0,1\} \);
  *   if moreover \( \mathrm{ConverterMinCapacityDesign} > 0 \) then \( x_c = 1 \).
  *   Otherwise, \( \mathrm{ConverterMinCapacityDesign} \le x_c \le \mathrm{ConverterMaxCapacityDesign} \).
  *
  * The parameter \p stvv and the Configuration for this function presented in
  * the BlockConfig (namely,
  * #f_BlockConfig->f_static_variables_Configuration) can be used to indicate
  * whether negative prices may occur. The parameter \p stvv has priority over
  * the BlockConfig in the sense that the Configuration in the BlockConfig is
  * only considered if no valid Configuration has been provided in \p stvv. By
  * default, it is assumed that negative prices do not occur and, therefore,
  * the binary variables are not generated. If the Configuration is a
  * SimpleConfiguration< int >, then a nonzero value stored in this
  * Configuration indicates that negative prices may occur. The value zero
  * indicates that negative prices do not occur. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the BatteryUnitBlock
 /** Method that generates the static constraint of the BatteryUnitBlock. The
  * operations of the battery storage unit are described on a discrete time
  * horizon as dictated by the UnitBlock interface. In this description we
  * indicate it with
  * \f$ \mathcal{T} = \{ 0 , \ldots , | \mathcal{T} | - 1 \} \f$. The main
  * constraints of this unit are defined as:
  *
  * - maximum and minimum power output constraints according to primary and
  *   secondary spinning reserves, as presented in (1)–(2). Each of them is a
  *   std::vector< FRowConstraint > with dimension f_time_horizon; the entry
  *   \f$ t = 0 , \ldots , \mathrm{f\_time\_horizon} - 1 \f$ is the maximum or
  *   minimum power output value according to the primary and the secondary
  *   spinning reserves at time \f$ t \f$. These ensure the maximum (or
  *   minimum) amount of energy that the unit can produce (or use).
  *
  *   \f[
  *     p^{ac}_t + p^{pr}_t + p^{sc}_t \leq \kappa P^{mx,b}_t
  *                                         \quad t \in \mathcal{T} \quad (1)
  *   \f]
  *
  *   \f[
  *     \kappa P^{mn,b}_t \leq p^{ac}_t - p^{pr}_t - p^{sc}_t
  *                                         \quad t \in \mathcal{T} \quad (2)
  *   \f]
  *
  *   where \f$ P^{mx,b}_t \f$ and \f$ P^{mn,b}_t \f$ are the maximum and
  *   minimum power output parameters for each time \f$ t \f$ of the time
  *   horizon \f$ \mathcal{T} \f$, respectively. In the design scenario of the
  *   UC problem, they become respectively:
  *
  *   \f[
  *     p^{ac}_t + p^{pr}_t + p^{sc}_t \leq x_b \, ( \kappa P^{mx,b}_t )
  *                                         \quad t \in \mathcal{T} \quad (1)
  *   \f]
  *
  *   \f[
  *     x_b \, ( \kappa P^{mn,b}_t ) \leq p^{ac}_t - p^{pr}_t - p^{sc}_t
  *                                         \quad t \in \mathcal{T} \quad (2)
  *   \f]
  *
  *   where \f$ x_b \f$ is the design variable of the battery.
  *
  * - ramp-up and ramp-down constraints are presented in (3)–(4). Each of
  *   them is a std::vector< FRowConstraint > with dimension f_time_horizon;
  *   the entry \f$ t = 0 , \ldots , \mathrm{f\_time\_horizon} - 1 \f$ gives
  *   the ramp-up and ramp-down constraints:
  *
  *   \f[
  *     p^{ac}_t - p^{ac}_{t-1} \leq \Delta^{up}_t
  *                                         \quad t \in \mathcal{T} \quad (3)
  *   \f]
  *
  *   \f[
  *     p^{ac}_t - p^{ac}_{t-1} \geq - \Delta^{dn}_t
  *                                         \quad t \in \mathcal{T} \quad (4)
  *   \f]
  *
  *   where \f$ \Delta^{up}_t \f$ and \f$ \Delta^{dn}_t \f$ are the ramp-up
  *   and ramp-down thresholds for each time \f$ t \f$.
  *
  * - active power relation with intake and outtake levels constraints are
  *   presented in (5). Each of them is a std::vector< FRowConstraint > with
  *   dimension f_time_horizon; the entry \f$ t = 0 , \ldots ,
  *   \mathrm{f\_time\_horizon} - 1 \f$ enforces the relation:
  *
  *   \f[
  *     p^{ac}_t = p^+_t - p^-_t
  *                                         \quad t \in \mathcal{T} \quad (5)
  *   \f]
  *
  *   The equations (6.1)–(6.2) indicate the upper bounds of intake and
  *   outtake levels at each time instant \f$ t \f$:
  *
  *   \f[
  *     p^+_t \leq \kappa C^+ P^{mx,b}_t
  *                                       \quad t \in \mathcal{T} \quad (6.1)
  *   \f]
  *
  *   \f[
  *     p^-_t \leq \kappa C^- ( -P^{mn,b}_t )
  *                                       \quad t \in \mathcal{T} \quad (6.2)
  *   \f]
  *
  *   where \f$ C^+ \f$ and \f$ C^- \f$ are the C-rates of the battery in
  *   charge and discharge, respectively; in the design scenario they become:
  *
  *   \f[
  *     p^+_t \leq \kappa x_b \, ( C^+ P^{mx,b}_t )
  *                                       \quad t \in \mathcal{T} \quad (6.3)
  *   \f]
  *
  *   \f[
  *     p^-_t \leq \kappa x_b \, ( C^- ( -P^{mn,b}_t ) )
  *                                       \quad t \in \mathcal{T} \quad (6.4)
  *   \f]
  *
  *   \f[
  *     p^+_t + p^-_t \leq x_c \, ( \kappa P^{mx,c}_t )
  *                                       \quad t \in \mathcal{T} \quad (6.5)
  *   \f]
  *
  *   where \f$ P^{mx,c}_t \f$ is the maximum power of the converter, and
  *   \f$ x_b \f$ and \f$ x_c \f$ are the battery and converter design
  *   variables, respectively.
  *
  * - storage level relation with intake and outtake levels (if any) is
  *   presented in (7), a std::vector< FRowConstraint > with dimension
  *   f_time_horizon:
  *
  *   \f[
  *     v^{ba}_t = v^{ba}_{t-1} + ρ^+_t · p^+_t − ρ^-_t · p^-_t − d^{ba}_t
  *                                         \quad t \in \mathcal{T} \quad (7)
  *   \f]
  *
  *   Note that (7) changes as below in the special case (see note above),
  *   giving (8):
  *
  *   \f[
  *     v^{ba}_t = v^{ba}_{t-1} + p^{ac}_t − d^{ba}_t
  *                                         \quad t \in \mathcal{T} \quad (8)
  *   \f]
  *
  *   The bounds on storage levels at each time instant \f$ t \f$ are:
  *
  *   \f[
  *     v^{ba}_t \in [ \kappa V^{mn}_t , \kappa V^{mx}_t ]
  *                                       \quad t \in \mathcal{T} \quad (9a)
  *   \f]
  *
  *   which, in the design scenario, become:
  *
  *   \f[
  *     x_b \, ( \kappa V^{mn}_t ) \leq v^{ba}_t \leq
  *     x_b \, ( \kappa V^{mx}_t ) \quad t \in \mathcal{T} \quad (9b)
  *   \f]
  *
  *   where \f$ \rho^+_t \f$ and \f$ \rho^-_t \f$ are
  *   StoringBatteryRho and ExtractingBatteryRho, \f$ V^{mn}_t \f$ and
  *   \f$ V^{mx}_t \f$ are the minimum and maximum storage levels, and
  *   \f$ x_b \f$ is the battery design variable.
  *
  * - binary variable relation with storing and extracting energy levels (if
  *   any) are presented in (10)–(11), each a std::vector< FRowConstraint >
  *   with dimension f_time_horizon:
  *
  *   \f[
  *     p^+_t \leq u^+_t P^{mx,b}_t
  *                                       \quad t \in \mathcal{T} \quad (10)
  *   \f]
  *
  *   \f[
  *     p^-_t \leq - ( 1 - u^+_t ) P^{mn,b}_t
  *                                       \quad t \in \mathcal{T} \quad (11)
  *   \f]
  *
  *   When \f$ \rho^+_t = \rho^-_t = 1 \f$, the binary variable \f$ u^+_t \f$
  *   is not required and neither are (10)–(11).
  *
  * - primary and secondary reserve upper bounds (if any), presented in (12)
  *   and (13), each a std::vector< FRowConstraint > with dimension
  *   f_time_horizon:
  *
  *   \f[
  *     p^{pr}_t \leq P^{mx,pr}_t
  *                                       \quad t \in \mathcal{T} \quad (12)
  *   \f]
  *
  *   \f[
  *     p^{sc}_t \leq P^{mx,sc}_t
  *                                       \quad t \in \mathcal{T} \quad (13)
  *   \f]
  *
  * - design bounds for the installation variables:
  *
  *   Continuous design case:
  *   \[
  *     \mathrm{BatteryMinCapacityDesign} \le x_b \le \mathrm{BatteryMaxCapacityDesign} \qquad (14)
  *   \]
  *   \[
  *     \mathrm{ConverterMinCapacityDesign} \le x_c \le \mathrm{ConverterMaxCapacityDesign} \qquad (15)
  *   \]
  *
  *   Binary design case (when the corresponding MaxCapacityDesign < 0):
  *   \[
  *     x_b \in \{0,1\} \quad (\text{and } \mathrm{BatteryMinCapacityDesign} > 0 \Rightarrow x_b = 1)
  *   \]
  *   \[
  *     x_c \in \{0,1\} \quad (\text{and } \mathrm{ConverterMinCapacityDesign} > 0 \Rightarrow x_c = 1)
  *   \]
  */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the BatteryUnitBlock
 /** Method that generates the objective of the BatteryUnitBlock. The
  * objective can include:
  *
  * - a linear term on intake and outtake power, if the vector "Cost" is
  *   provided:
  *   \f[
  *     \min \ \sum_{t \in \mathcal{T}} C[t] \cdot (p^+_t + p^-_t)
  *   \f]
  *   (coefficients are also scaled by the Block scale factor, if any);
  *
  * - an investment term in design mode, if "BatteryInvestmentCost" and/or
  *   "ConverterInvestmentCost" are nonzero:
  *   \f$ + \ I_b \cdot x_b \;+\; I_c \cdot x_c \f$.
  *
  * Hence, in the design scenario the full objective is:
  * \f[
  *   \min \ \sum_{t \in \mathcal{T}} C[t] \cdot (p^+_t + p^-_t)
  *          \;+\; I_b \cdot x_b \;+\; I_c \cdot x_c \; .
  * \f]
  * If "Cost" is not provided, \f$ C[t] = 0 \f$ and only the investment terms
  * remain in design mode. If both investment costs are zero or not provided,
  * the objective reduces to the operational cost term only.
  */

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
  * - If \p fsbc is not nullptr and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr and it is a pointer to a
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
  * - Otherwise, by default, the tolerance is 0 and the relative violation is
  *   considered.
  *
  * This function currently considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * currently ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variables only. Notice also that, before checking if the solution
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
 * kinds of battery storage units
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
 double get_batt_max_capacity_design( void ) const {
  return ( f_BattMaxCapacityDesign );
  }

 /// returns the maximum converter installable capacity by the user
 double get_conv_max_capacity_design( void ) const {
  return ( f_ConvMaxCapacityDesign );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of minimum storage
 /** This method returns a vector V containing the minimum storage at each
  * time instant. There are three possible cases:
  *
  * - if the vector is empty, then the minimum storage of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the minimum storage
  *   of the unit for all time instants;
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
  * - if the vector has only one element, then V[ 0 ] is the maximum storage
  *   of the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each V[ t ]
  *   represents the maximum storage value at time t. */

 const std::vector< double > & get_max_storage( void ) const {
  return( v_MaxStorage );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power at a given time
 /** This method returns the minimum active power of the unit at time t.
  *
  * There are three possible cases:
  *
  * - if the internal vector is empty, the minimum power is 0 for all t;
  *
  * - if the internal vector has only one element, that value is the
  *   minimum power for all time instants;
  *
  * - otherwise, the internal vector has size get_time_horizon() and
  *   the entry at position t represents the minimum power at time t.
  */

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( v_MinPower[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power at a given time
 /** This method returns the maximum active power of the unit at time t.
  *
  * There are three possible cases:
  *
  * - if the internal vector is empty, the maximum power is 0 for all t;
  *
  * - if the internal vector has only one element, that value is the
  *   maximum power for all time instants;
  *
  * - otherwise, the internal vector has size get_time_horizon() and
  *   the entry at position t represents the maximum power at time t.
  */

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum converter power at a given time
 /** This method returns the maximum active power of the converter at time t.
  *
  * There are three possible cases:
  *
  * - if the internal vector is empty, the maximum converter power is 0 for all t;
  *
  * - if the internal vector has only one element, that value is the
  *   maximum converter power for all time instants;
  *
  * - otherwise, the internal vector has size get_time_horizon() and
  *   the entry at position t represents the maximum converter power at
  *   time t. */

 const std::vector< double > & get_converter_max_power( void ) const {
  return( v_ConvMaxPower );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum reactive power of the \p generator at time \p t

 double get_min_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MinReactivePower.size() > t ) ? v_MinReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of the \p generator at time \p t

 double get_max_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MaxReactivePower.size() > t ) ? v_MaxReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum primary reserve power
 /** This method returns a vector containing the maximum active power that can
  * be used as primary reserve. There are three possible cases:
  *
  * - if this vector is empty, then the maximum primary power of the unit is
  *   0;
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
  * - if this vector is empty, then the maximum secondary power of the unit is
  *   0;
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
 /** This method returns a vector V containing the delta ramp up at all time
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
 /** This method returns a vector V containing the delta ramp down at all time
  * instants. There are three possible cases:
  *
  * - if the vector is empty, then the delta ramp down of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the delta ramp down
  *   of the unit for all time instants;
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
  * - if the vector is empty, then the storing battery rho of the unit is 1;
  *
  * - if the vector has only one element, then V[ 0 ] is the storing battery
  *   rho of the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the storing battery rho value at time t. */

 const std::vector< double > & get_storing_battery_rho( void ) const {
  return( v_StoringBatteryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inefficiency of extracting energy of the unit
 /** This method returns a vector V containing the extracting battery rho
  * (inefficiency of extracting energy of the unit) at all time instants.
  * There are three possible cases:
  *
  * - if the vector is empty, then the extracting battery rho of the unit is
  *   1;
  *
  * - if the vector has only one element, then V[ 0 ] is the extracting
  *   battery rho of the unit for all time instants;
  *
  * - otherwise, the vector V must have size get_time_horizon() and each
  *   V[ t ] represents the extracting battery rho value at time t. */

 const std::vector< double > & get_extracting_battery_rho( void ) const {
  return( v_ExtractingBatteryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of storage and extraction cost of energy
 /** This method returns a vector V containing the monetary cost of storing
  * one unit of energy into, or extracting it from, the battery (the cost is
  * the same in both cases) at all time instants. There are three possible
  * cases:
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

 std::vector< ColVariable > & get_storage_level( void ) {
  return( v_storage_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const vector of storage level variables
 /** This method returns a const vector V containing the storage level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the storage
  *   level variable for time step t. */

 const std::vector< ColVariable > & get_const_storage_level( void ) const {
  return( v_storage_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of intake level variables
 /** This method returns a vector V containing the intake level variables.
  * There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the intake
  *   level variable for time step t. */

 std::vector< ColVariable > & get_intake_level( void ) {
  return( v_intake_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const vector of intake level variables
 /** This method returns a const vector V containing the intake level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the intake
  *   level variable for time step t. */

 const std::vector< ColVariable > & get_const_intake_level( void ) const {
  return( v_intake_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of outtake level variables
 /** This method returns a vector V containing the outtake level variables.
  * There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the outtake
  *   level variable for time step t. */

 std::vector< ColVariable > & get_outtake_level( void ) {
  return( v_outtake_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const vector of outtake level variables
 /** This method returns a const vector V containing the outtake level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the outtake
  *   level variable for time step t. */

 const std::vector< ColVariable > & get_const_outtake_level( void ) const {
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
 /// returns the vector of reactive power variables

 ColVariable * get_reactive_power( Index generator ) override {
  if( v_reactive_power.empty() )
   return( nullptr );
  return( &( v_reactive_power.front() ) );
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

 ColVariable & get_batt_design( void ) { return( batt_design ); }

/*--------------------------------------------------------------------------*/
 /// returns the const battery design variable

 const ColVariable & get_const_batt_design( void ) const {
  return( batt_design );
  }

/*--------------------------------------------------------------------------*/
 /// returns the converter design variable

 ColVariable & get_conv_design( void ) { return( conv_design ); }

/*--------------------------------------------------------------------------*/
 /// returns the const converter design variable

 const ColVariable & get_const_conv_design( void ) const {
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
 /// returns the outtake upper bound constraint + binary variables for time t
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
  return( &( intake_outtake_bounds_Const[ 0 ][ t ] ) );
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
  return( &( intake_outtake_bounds_Const[ 1 ][ t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the primary reserve bound constraints
 const std::vector< LB0Constraint > & get_primary_reserve_bounds( void )
  const {
  return( primary_upper_bound_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the secondary reserve bound constraints
 const std::vector< LB0Constraint > & get_secondary_reserve_bounds( void )
  const {
  return( secondary_upper_bound_Const );
 }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution storing the current solution of this BatteryUnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this BatteryUnitBlock. This
  * is a BatteryUnitBlockSolution extending UnitBlockSolution with the
  * specific extra solution information of BatteryUnitBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - the first four bits (bit 0 to bit 3) are "taken" by the base
  *   UnitBlock[Solution]
  *
  * - bit 4 (& 16) means "store the storage levels"
  *
  * - bit 5 (& 32) means "store the intakes and outtakes"
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 63 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [Battery]UnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

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
  * zero. */

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
  * Range \p rng does not contain the 0 index, this function does nothing. */

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
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

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
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

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
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

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
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_kappa( double value , c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector = { value };
  set_kappa( vector.cbegin() , Range( 0 , Inf< Index >() ) ,
             issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

 // For the Range version, use the default implementation defined in UnitBlock
 using UnitBlock::scale;

/** @} ---------------------------------------------------------------------*/
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

 /// the vector of MinReactivePower
 std::vector< double > v_MinReactivePower;

 /// the vector of MaxReactivePower
 std::vector< double > v_MaxReactivePower;

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
 double f_BattInvestmentCost;

 /// the converter investment cost
 double f_ConvInvestmentCost;

 /** the minimum battery capacity design allowed (lower bound on x_b in
  * design mode); default 0.
  * If BatteryMaxCapacityDesign < 0 (binary), BatteryMinCapacityDesign > 0 
  * forces x_b = 1 */
 double f_BattMinCapacityDesign;

 /// the maximum battery capacity design allowed
 /** If < 0, x_b is binary; if > 0, x_b is continuous with bounds
  * [BatteryMinCapacityDesign, BatteryMaxCapacityDesign]. */
 double f_BattMaxCapacityDesign;

 /** the minimum converter capacity design allowed (lower bound on x_c in
  * design mode); default 0.
  * If ConverterMaxCapacityDesign < 0 (binary),
  * ConverterMinCapacityDesign > 0 forces x_c = 1 */
 double f_ConvMinCapacityDesign;

 /** the maximum converter capacity design allowed
  * If < 0, x_c is binary; if > 0, x_c is continuous with bounds
  * [ConverterMinCapacityDesign, ConverterMaxCapacityDesign] */
 double f_ConvMaxCapacityDesign;

 /// the InitialStorage value
 double f_InitialStorage;

 /// the InitialPower value
 double f_InitialPower;

 /// the MaxCRateCharge value
 double f_MaxCRateCharge;

 /// the MaxCRateDischarge value
 double f_MaxCRateDischarge;

 /// the kappa value
 double f_kappa;

 /// the scale factor
 double f_scale;

 /// the reference Schedule to deviate minimally from if there
 std::vector< double > v_RefSchedule;

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

 /// the reactive power variables
 std::vector< ColVariable > v_reactive_power;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

 /// the battery design variable
 ColVariable batt_design;

 /// the converter design variable
 ColVariable conv_design;

 /// the variables for deviation to reference schedule
 std::vector< ColVariable > v_abs_ref_schedule;

/*------------------------------- constraints ------------------------------*/

/// the reference schedule constraints
 std::vector< FRowConstraint > Reference_Schedule_Const;

 /// the active power bounds constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_Const;

 /// the battery design bound constraint
 BoxConstraint batt_design_bound_Const;

 /// the converter design bound constraint
 BoxConstraint conv_design_bound_Const;

 /// the active power bounds design constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_design_Const;

 /// the intake outtake upper bounds design constraints
 boost::multi_array< FRowConstraint , 2 >
                                   intake_outtake_upper_bounds_design_Const;

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

 /// the reactive power bound constraints
 std::vector< BoxConstraint > ReactivePower_Bound_Const;

 /*!! Q <= P
 std::vector< FRowConstraint > Reactive_2_Active_Const;
 */

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
  * @param issueAMod Controls how abstract Modifications are issued. */

 void update_objective( c_ModParam issueAMod ) const;

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this BatteryUnitBlock is consistent
 /** This function checks whether the data in this BatteryUnitBlock is
  * consistent. The data is consistent if all the following conditions are met.
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
  *   reserves is nonnegative.
  *
  * - Design bounds consistency (battery):
  *   - \( \mathrm{BatteryMinCapacityDesign} \ge 0 \);
  *   - if \( \mathrm{BatteryMaxCapacityDesign} > 0 \), then
  *     \( \mathrm{BatteryMinCapacityDesign} \le \mathrm{BatteryMaxCapacityDesign} \);
  *   - if \( |\mathrm{BatteryMaxCapacityDesign}| = 1 \), then
  *     \( \mathrm{BatteryMinCapacityDesign} \le 1 \);
  *   - if \( \mathrm{BatteryMaxCapacityDesign} < 0 \) (binary), then
  *     \( \mathrm{BatteryMinCapacityDesign} \le 1 \)
  *     (note: \( \mathrm{BatteryMinCapacityDesign} > 0 \Rightarrow x_b = 1 \)).
  *
  * - Design bounds consistency (converter):
  *   - \( \mathrm{ConverterMinCapacityDesign} \ge 0 \);
  *   - if \( \mathrm{ConverterMaxCapacityDesign} > 0 \), then
  *     \( \mathrm{ConverterMinCapacityDesign} \le \mathrm{ConverterMaxCapacityDesign} \);
  *   - if \( |\mathrm{ConverterMaxCapacityDesign}| = 1 \), then
  *     \( \mathrm{ConverterMinCapacityDesign} \le 1 \);
  *   - if \( \mathrm{ConverterMaxCapacityDesign} < 0 \) (binary), then
  *     \( \mathrm{ConverterMinCapacityDesign} \le 1 \)
  *     (note: \( \mathrm{ConverterMinCapacityDesign} > 0 \Rightarrow x_c = 1 \)).
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

  register_method< BatteryUnitBlock , MF_dbl_it , Subset && , bool >(
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
   default: ;
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
                          const Block::Range & rng )
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
/*-------------------- CLASS BatteryUnitBlockSolution ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [UnitBlock]Solution of a BatteryUnitBlock
/** The BatteryUnitBlockSolution class derives from UnitBlockSolution and
 * adds to the "standard" information stored in there (active power, possibly
 * commitment and primary/secondary reserve) the other information that is
 * typical of the BatteryUnitBlock, i.e.,
 *
 * - [possibly] the storage level of the battery at each time instant
 *
 * - [possibly] the intake/outtake in the battery at each time instant;
 *   since the battery is supposed to never be charged and discharged at
 *   the same time instant, the value is positive if the battery is being
 *   charged (intake) and negative if it is being discharged (outtake)
 *
 * - if defined, the value of the Battery Design Variable
 *
 * - if defined, the value of the Converter Design Variable */

class BatteryUnitBlockSolution : public UnitBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*----------------------------- CONSTANTS ----------------------------------*/

 static constexpr double dNaN = std::numeric_limits< double >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

/*------------------------------- FRIENDS ----------------------------------*/

 friend BatteryUnitBlock;  ///< make BatteryUnitBlock friend

/*--------- CONSTRUCTING AND DESTRUCTING BatteryUnitBlockSolution ----------*/

 /// constructor, it has nothing to do
 explicit BatteryUnitBlockSolution( void ) :
  f_b_design( dNaN ) , f_c_design( dNaN ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~BatteryUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*----- METHODS DESCRIBING THE BEHAVIOR OF A BatteryUnitBlockSolution -----*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a BatteryUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a BatteryUnitBlockSolution into a netCDF::NcGroup. The format
  * is the one of UnitBlockSolution [cf. UnitBlockSolution::serialize()],
  * plus:
  *
  * - The variable "StorageLevel", of type netCDF::NcDouble and indexed over
  *   the dimension "TimeHorizon"; StorageLevel[ t ] is the optimal value of
  *   the storage level of the battery at time \f$ t \f$. The variable is
  *   optional.
  *
  * - The variable "InOutTake", of type netCDF::NcDouble and indexed over
  *   the dimension "TimeHorizon"; InOutTake[ t ] is the amount of energy
  *   being charged in the battery at time \f$ t \f$ (negative if it is
  *   discharged). The variable is optional.
  *
  * Note that, unlike those of the base class, these variables do not need
  * to be indexed over the dimension "NumberGenerators" since
  * BatteryUnitBlock always has exactly one generator.
  *
  * - The scalar variable "BatteryDesign", of type netCDF::NcDouble, that
  *   represent the value of the dimensioning variable of the battery;
  *   the variable is optional in that the battery may not have any
  *   dimensioning variable.
  *
  * - The scalar variable "ConverterDesign", of type netCDF::NcDouble, that
  *   represent the value of the dimensioning variable of the converter;
  *   the variable is optional in that the converter may not have any
  *   dimensioning variable. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 BatteryUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 BatteryUnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "BatteryUnitBlockSolution [" << this << "]: " << std::endl;
 }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 std::vector< double > v_storage;
 ///< v_storage[ t ] = value of stored energy at time t

 std::vector< double > v_intake;  ///< v_intake[ t ] = intake at time t

 double f_b_design;    ///< the value of the battery dimensioning variable

 double f_c_design;    ///< the value of the converter dimensioning variable

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( BatteryUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __BatteryUnitBlock */

/*--------------------------------------------------------------------------*/
/*---------------------- End File BatteryUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
