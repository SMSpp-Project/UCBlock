/*--------------------------------------------------------------------------*/
/*------------------------- File ThermalUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class ThermalUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" thermal unit
 * of a Unit Commitment Problem. A ThermalUnitBlock corresponds to a single
 * electrical generator.
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
 * \author Tiziano Bacci \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli, Tiziano Bacci
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ThermalUnitBlock
 #define __ThermalUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

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
/*----------------------- CLASS ThermalUnitBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a thermal unit
/** The ThermalUnitBlock class derives from UnitBlock and implements a
 * "reasonably standard" thermal unit of a Unit Commitment Problem. That is,
 * the class is designed in order to give mathematical formulation to describe
 * the operation of large set of conventional power plants (such as nuclear,
 * hard coal, gas turbine, gas, combined cycle, oil, ...) which are directly
 * connected to the transmission grid. The technical and physical constraints
 * are mainly divided in four different categories:
 *
 * - minimum up and down time constraints;
 *
 * - ramp-up/down rate constraints;
 *
 * - maximum and minimum power output constraints;
 *
 * - active power relation with primary and secondary spinning reserves. */

class ThermalUnitBlock : public UnitBlock
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
 /** Constructor of ThermalUnitBlock, taking possibly a pointer of its father
  * Block. */

 explicit ThermalUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of ThermalUnitBlock

 virtual ~ThermalUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the ThermalUnitBlock. Besides the mandatory "type" attribute of any :Block,
  * the group must contain all the data required by the base UnitBlock, as
  * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
  * In particular, we refer to that description for the crucial dimensions
  * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
  * must then also contain:
  *
  * - The variable "MinPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MnP[ t ] that, for
  *   each time instant t, contains the nominal minimum active power output
  *   value of the unit for the corresponding time step. If "MinPower" has
  *   length 1 then MnP[ t ] contains the same value for all t. Otherwise,
  *   MinPower[ i ] is the fixed value of MnP[ t ] for all t in the interval [
  *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. Note that it must be MnP[ t ] >= 0 for
  *   all t. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
  *   the mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded.
  *
  * - The variable "MaxPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MxP[ t ] that, for
  *   each time instant t, contains the nominal maximum active power output
  *   value of the unit for the corresponding time step. If "MaxPower" has
  *   length 1 then MxP[ t ] contains the same value for all t. Otherwise,
  *   MaxPower[ i ] is the fixed value of MxP[ t ] for all t in the interval [
  *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. Note that it must be MxP[ t ] >= MnP[ t
  *   ] >= 0 for all t. If NumberIntervals <= 1 or NumberIntervals >=
  *   TimeHorizon, then the mapping clearly does not require "ChangeIntervals",
  *   which in fact is not loaded.
  *
  * - The variable "Availability", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector Av[ t ] that, for
  *   each time instant t, contains the availability of the unit for the
  *   corresponding time step. If "Availability" has length 1 then Av[ t ] is
  *   equal to the single given value in "Availability" for all t. Otherwise,
  *   Availability[ i ] is the fixed value of Av[ t ] for all t in the interval
  *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  *   The availability of the unit is given by a number between 0 and 1. Let t
  *   be a time instant in {0, ..., TimeHorizon - 1}. The operational
  *   (effective) maximum active power output of the unit at time t is given by
  *   Av[ t ] * MxP[ t ] (see the variable "MaxPower" for the definition of
  *   MxP). The operational minimum active power output of the unit at time t
  *   is zero if Av[ t ] == 0 and it is MnP[ t ] if Av[ t ] > 0 (see the
  *   variable "MinPower" for the definition of MnP).
  *
  *   This variable is optional. If it is not provided, then the unit is fully
  *   operational at all time instants, i.e., we assume that Av[ t ] = 1 for
  *   all t in {0, ..., TimeHorizon - 1}.
  *
  * - The variable "DeltaRampUp", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector DP[ t ] that, for
  *   each time instant t, contains the ramp-up value of the unit for the
  *   corresponding time step, i.e., the maximum possible increase of active
  *   power production w.r.t. the power that had been produced in time instant
  *   t - 1, if any. This variable is optional; if it is not provided then it
  *   is assumed that DP[ t ] == MxP[ t ], i.e., the unit can ramp up by an
  *   arbitrary amount, i.e., there are no ramp-up constraints. If
  *   "DeltaRampUp" has length 1 then DP[ t ] contains the same value for all
  *   t. Otherwise, DeltaRampUp[ i ] is the fixed value of DP[ t ] for all t in
  *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
  *   the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1
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
  *   produced in time instant t - 1, if any. This variable is optional; if
  *   it is not provided then it is assumed that DP[ t ] == MxP[ t ], i.e.,
  *   the unit can ramp down an arbitrary amount, i.e., there are no
  *   ramp-down constraints. If "DeltaRampDown" has length 1 then DM[ t ]
  *   contains the same value for all t. Otherwise, DeltaRampDown[ i ] is the
  *   fixed value of DM[ t ] for all t in the interval [ ChangeIntervals[ i -
  *   1 ] , ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[
  *   - 1 ] = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon,
  *   then the mapping clearly does not require "ChangeIntervals", which in
  *   fact is not loaded.
  *
  * - The variable "PrimaryRho", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector PR[ t ] that, for
  *   each time instant t, contains the maximum possible fraction of active
  *   power that can be used as primary reserve value of the unit for the
  *   corresponding time step. This variable is optional; if it is not provided
  *   then it is assumed that this unit may not be capable of producing any
  *   primary reserve, which correspond to PR[ t ] == 0 for all t. If
  *   "PrimaryRho" has length 1 then PR[ t ] contains the same value for all
  *   t. Otherwise, PrimaryRho[ i ] is the fixed value of PR[ t ] for all t in
  *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ] with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
  *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "SecondaryRho", of type netCDF::NcDouble and either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   SR[ t ] that, for each time instant t, contains the maximum possible
  *   fraction of active power that can be used as secondary reserve value of
  *   the unit for the corresponding time step. This variable is optional; if
  *   it is not provided then it is assumed that this unit may not be capable
  *   of producing any secondary reserve, which correspond to SR[ t ] == 0
  *   for all t. If "SecondaryRho" has length 1 then SR[ t ] contains the
  *   same value for all t. Otherwise, SecondaryRho[ i ] is the fixed value
  *   of SR[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ] with the assumption that ChangeIntervals[ - 1 ]
  *   = 0. If "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon"
  *   then the mapping clearly does not require "ChangeIntervals", which in
  *   fact is not loaded.
  *
  * - The variable "QuadTerm", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector A[ t ] that, for
  *   each time instant t, contains the quadratic term of power cost function
  *   of the unit for the corresponding time step. This variable is optional;
  *   if it is not provided then it is assumed that A[ t ] == 0, i.e., the cost
  *   of the unit is linear in the produced power. If "QuadTerm" has length 1
  *   then A[ t ] contains the same value for all t. Otherwise, QuadTerm[ i ]
  *   is the fixed value of A[ t ] for all t in the interval [ ChangeIntervals[
  *   i - 1 ] , ChangeIntervals[ i ] ], with the assumption that
  *   ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
  *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded.
  *
  * - The variable "StartUpCost", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector SC[ t ] that, for
  *   each time instant t, contains the start-up cost value of the unit for the
  *   corresponding time step. This variable is optional; if it is not provided
  *   then it is assumed that SC[ t ] == 0, i.e., this unit may not have any
  *   start-up cost. If "StartUpCost" has length 1 then SC[ t ] contains the
  *   same value for all t. Otherwise, StartUpCost[ i ] is the fixed value of
  *   SC[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
  *   0. If "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then
  *   the mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded.
  *
  * - The variable "LinearTerm", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector B[ t ] that, for
  *   each time instant t, contains the linear term of power cost function of
  *   the unit for the corresponding time step. This variable is optional; if
  *   it is not provided then it is assumed that B[ t ] == 0, i.e., the cost of
  *   the unit has no linear dependence on the produced power (say, only the
  *   quadratic one). If "LinearTerm" has length 1 then A[ t ] contains the
  *   same value for all t. Otherwise, LinearTerm[ i ] is the fixed value of
  *   B[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
  *   0. If "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then
  *   the mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded.
  *
  * - The variable "ConstTerm", of type netCDF::NcDouble and to be either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   C[ t ] that, for each time instant t, contains the constant term of
  *   power cost function of the unit for the corresponding time step. This
  *   variable is optional; if it is not provided then it is assumed that C[
  *   t ] == 0, i.e., the cost of the unit has no fixed term, only those
  *   depending (linearly or quadratically) on the produced power. If
  *   "ConstTerm" has length 1 then C[ t ] contains the same value for all t.
  *   Otherwise, ConstTerm[ i ] is the fixed value of C[ t ] for all t in the
  *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1
  *   or "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded.
  *
  * - The scalar variable "InitialPower", of type netCDF::NcDouble and not
  *   indexed over any dimension. This variable indicates the amount of the
  *   power that the unit was producing at time instant -1, i.e., before the
  *   start of the time horizon; this is necessary to compute the ramp-up and
  *   ramp-down constraints. This variable is optional. If it is not
  *   provided, then it is taken to be 0. Clearly, it must be that
  *   MaxPower >= InitialPower >= MinPower if the unit was "on" at time
  *   instant - 1, and it would be ignored if the unit was "off" at time
  *   instant - 1. The on/off status of the unit is also encoded by the scalar
  *   variable InitUpDownTime: in particular, InitUpDownTime > 0 then the
  *   unit was on at time instant -1, and therefore InitialPower >= MinPower
  *   must hold, while if InitUpDownTime <= 0 then the unit was off at time
  *   instant - 1, and therefore InitialPower is ignored. In fact, if
  *   InitUpDownTime <= 0 then this variable need not be defined since it is
  *   not loaded.
  *
  * - The scalar variable "InitUpDownTime", of type netCDF::NcInt and not
  *   indexed over any dimension and indicates the initial time to generating
  *   the unit. If InitUpDownTime > 0, this means that the unit has been on
  *   for InitUpDownTime timestamps prior to timestamp 0 (the beginning of
  *   the horizon). If, instead, InitUpDownTime <= 0, this means that the unit
  *   has been off for - InitUpDownTime timestamps prior to timestamp 0; note
  *   that InitUpDownTime == 0 means that the unit has been just shut-down at
  *   the end of time instant -1, i.e., the beginning of time instant 0. This
  *   variable is optional. If it is not provided, then it is taken to be
  *   -MinDownTime if InitialPower == 0, and MinUpTime if InitialPower > 0.
  *
  * - The positive scalar variable "MinUpTime", of type netCDF::NcUint and not
  *   indexed over any dimension, which indicates the minimum allowed up time
  *   in this unit. This variable is optional, if it is not provided it is
  *   taken to be MinUpTime == 1. Since MinUpTime == 0 is actually possible,
  *   which means that the unit can shut-down in the very same timestamp in
  *   which it starts up, we also know that starting up a unit only to
  *   power it down again immediately is never a good idea, so we force
  *   up-time periods to be at least of length 1, even if it is given equals
  *   to 0.
  *
  * - The positive scalar variable "MinDownTime", of type netCDF::NcUint and
  *   not indexed over any dimension, which indicates the minimum allowed down
  *   time in this unit. This variable is optional, if it is not provided it is
  *   taken to be MinDownTime == 1. Since MinDownTime == 0 is actually
  *   possible, which means that the unit can start-up in the very same
  *   timestamp in which it shuts down, we also know that shutting down a
  *   unit only to power it up again immediately is never a good idea, so
  *   we force down-time periods to be at least of length 1, even if it is
  *   given equals to 0.
  *
  * - The variable "FixedConsumption", of type netCDF::NcDouble and either
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon"), or having size 1. This is meant to represent the vector
  *   FC[ t ] which, for each time instant t, contains the fixed consumption
  *   of the power plant if it is OFF at time t. The variable is optional; if
  *   it is not defined, FC[ t ] == 0 for all time instants. If it has size
  *   1, then FC[ t ] == FixedConsumption[ 0 ] for all t, regardless to what
  *   "NumberIntervals" says. Otherwise, FixedConsumption[ i ] is the fixed
  *   value of FC[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ],
  *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
  *   = 0.
  *
  * - The variable "InertiaCommitment", of type netCDF::NcDouble and either
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
  *   not provided, then this variable can also be indexed over
  *   "TimeHorizon") or has size 1. This is meant to represent the vector
  *   IC[ t ] which, for each time instant t, contains the contribution that
  *   the unit can give to the inertia constraint for the sole fact that is is
  *   on (basically, the constant to be multiplied to the commitment
  *   variable) at time t. The variable is optional; if it is not defined,
  *   IC[ t ] == 0 for all time instants. If it has size 1, then IC[ t ] ==
  *   InertiaCommitment[ 0 ] for all t, regardless to what "NumberIntervals"
  *   says. Otherwise, InertiaCommitment[ i ] is the fixed value of IC[ t ]
  *   for all t in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[
  *   i ] ], with the assumption that ChangeIntervals[ - 1 ] = 0.
  *
  * - The variable "StartUpLimit", // TODO
  *
  * - The variable "ShutDownLimit", // TODO
  *
  */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the ThermalUnitBlock
 /** Method that generates the abstract Variable of the ThermalUnitBlock,
  * meanwhile deciding which of the different formulations of the problem is
  * produced as the "abstract representation" of the block.
  * The different possible formulations are represented by a single int value
  * "wf" that is obtained as follows:
  *
  * - if either \p stvv is not nullptr and it is a SimpleConfiguration< int >,
  *   or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_static_variables_Configuration is not nullptr,
  *   and it is a SimpleConfiguration< int >, then wf is the f_value of the
  *   SimpleConfiguration< int >
  *
  * - otherwise, wf is 0, i.e., the default is the T formulation.
  *
  * The list of supported formulations is:
  *
  * - wf & 3 == 0 is the "three binaries" (3bin) formulation. This
  * formulation of the ThermalUnitBlock class has six different variables:
  *
  * - the binary commitment variables which takes the value of 1 if unit is ON
  *   at time instant t and 0 otherwise;
  *
  * - the primary spinning reserve variables;
  *
  * - the secondary spinning reserve variables;
  *
  * - the active power variables.
  *
  * All of those variables are optional except the active power variables in
  * the sense that the model may just not have them and whenever a group of
  * above variables is created, its size will be the time horizon. Moreover,
  * ThermalUnitBlock is defined more groups of variables as follow:
  *
  * - the binary variable start_up status of the unit which takes the value of
  *   1 if the unit starts up at time instant t and 0 otherwise;
  *
  * - the binary variable shut_down status of the unit which takes the value
  *   of 1 if the unit shuts down at time instant t and 0 otherwise;
  *
  * These two groups of variables have size f_time_horizon - init_t, and
  * provide the unit commitment problem with a tight 3-binary MIP formulation.
  * Since these two variables may have shorter size (when init_t > 0), the
  * commitment variable needs to be fixed to 0 or 1 for the first init_t time
  * steps 0, ..., init_t - 1 (see initial time step concept in the
  * generate_abstract_constraints()).
  *
  * - wf & 3 == 1: the T formulation.
  *
  * // TODO add here details about T formulation
  *
  * - wf & 3 == 2: the "dynamic programming" inspired formulation (DP).
  *
  * // TODO add here details about DP formulation
  *
  * - wf & 3 == 3: the p_t formulation.
  *
  *  // TODO add here details about p_t formulation
  *
  * - wf & 3 == 4: the "start-up" formulation (SU).
  *
  * // TODO add here details about SU formulation
  *
  * - wf & 3 == 5: the "shut-down" formulation (SD).
  *
  * // TODO add here details about SD formulation
  *
  */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

// TODO following the code-flow / order, add in the method blow, for each
//  constraint, i.e., start-up / shut-down cnstrs, min/max power cnstrs, etc.,
//  the details about all the other formulations, and all the new other cnstrs
//  added, i.e., p/c cnstrs.
/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the ThermalUnitBlock
 /** This method generates the abstract constraints of the ThermalUnitBlock.
  *
  * The operations of the thermal generating unit are described on a discrete
  * time horizon as dictated by the UnitBlock interface. In this description we
  * indicate it with \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$.
  * Considering three parameters: InitUpDownTime \f$ \tau_0 \f$ which can be a
  * positive or negative (or 0) integer number and tells for how many time
  * steps before time step 0 the unit was ON (when \f$ \tau_0 > 0 \f$) or OFF
  * (when \f$ \tau_0 < 0 \f$), MaxUpTime \f$ \tau_+ \f$ which is a positive
  * integer number and indicates for how many time steps after time step 0 the
  * unit can remain ON, and MinDownTime \f$ \tau_- \f$ that is also a positive
  * integer number and indicates for how many time steps after time step 0 the
  * unit can remain off. Therefore, the starting-up (shutting-down) of
  * generating of each unit depends on these three parameters. Then, the first
  * time instant may not be always equal to zero. For this matter the concept
  * of first time instant which is called "init_t" is defined as below:
  *
  * - If \f$ \tau_0 > 0 \f$, this means that the unit has been on for
  *   \f$ \tau_0 \f$ timestamps prior to timestamp 0 (the beginning of the
  *   time horizon):
  *
  *   - if \f$ \tau_0 \geq \tau_+ \f$ then init_t = 0 ;
  *
  *   - otherwise init_t = \f$ \tau_+ \f$ - \f$ \tau_0 \f$ ;
  *
  * - If, instead, \f$ \tau_0 < 0 \f$, this means that the unit has
  *   been off for \f$ - \tau_0 \f$ timestamps prior to timestamp 0:
  *
  *   - if \f$ - \tau_0  \geq \tau_- \f$  then init_t = 0;
  *
  *   - otherwise init_t = \f$ \tau_- \f$ + \f$ \tau_0 \f$;
  *
  * - If the \f$ \tau_0 == 0 \f$ means that the unit has been just shut-down at
  *   the end of time instant -1, i.e., the beginning of time instant 0 and
  *   init_t == 0.
  *
  * - Note: for the entry a = 0, ..., init_t - 1 when commitment variables fix
  *   to 0 then active power variables fix to 0, but when commitment variables
  *   fix to 1 the bound constraints, a std::vector< LB0Constraint > with
  *   exactly init_t entries, the entry a = 0, ..., init_t - 1 being the bound
  *   constraints of the ColVariable corresponding to the active power
  *   production of the unit (put init_t := \f$ t_0 \f$).
  *
  * Then the main thermal unit constraints with three 3 binary variables
  * \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ are are presented as following:
  *
  * - Min Up/Down-time Constraints: a thermal unit may have minimum up and down
  *   time constraints and one possible representation of the constraints could
  *   be as below:
  *
  *   \f[
  *     u_t - u_{t-1} = v_t - w_t
  *          \quad t \in \{ t_0 , ...,\mathcal{T}- 1 \}               \quad (1)
  *   \f]
  *
  *   \f[
  *    \sum_{ s \in [ t - \tau_+  , t ] } v_s \leq
  *           u_t \quad t \in \{ \tau_+ + t_0, ..., \mathcal{T} - 1\}
  *                                                                   \quad (2)
  *   \f]
  *
  *   \f[
  *    \sum_{ s \in [ t - \tau_- , t ] } w_s \leq
  *         1 - u_t \quad t \in \{ \tau_- + t_0, ...,\mathcal{T} - 1\}
  *                                                                   \quad (3)
  *   \f]
  *
  * - since the variables commitment \f$ u_t \f$ have full size of time horizon
  *   and other two remain variables which are start-up \f$ v_t \f$ and shut-down
  *   \f$ w_t \f$ have size (f_time_horizon - init_t), the equalities (1) show
  *   a std::vector< FRowConstraint > with exactly (f_time_horizon - init_t)
  *   entries, the entry a = init_t, ...,(f_time_horizon - 1) being the startup
  *   and shut-down connecting constraints at time t. According to the concept
  *   of init_t when \f$ \tau_0 < 0 \f$  and \f$ -\tau_0 < \tau_- \f$ the
  *   commitment variables \f$ u_t \f$  are fixed to zero for init_t time
  *   steps (starting from zero till init_t - 1). When \f$ \tau_0 > 0 \f$  and
  *   \f$ \tau_0 < \tau_+ \f$ the commitment variables \f$ u_t \f$  are fixed
  *   to one for init_t time steps (starting from zero till init_t - 1). Since
  *   \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ are binary variables, we
  *   can ensure (for all periods \f$ t \in \{ t_0, ..., \mathcal{T} -1 \} \f$)
  *   that \f$ v_t = 1 \f$ if and only if \f$ u_t = 1 \f$ and
  *   \f$ u_{t-1} = 0 \f$.
  *   It also obvious that \f$ w_t = 1 \f$ if and only if \f$ u_t = 0 \f$ and
  *   \f$ u_{t-1} = 1 \f$. These conditions are satisfied by equality (1).
  *
  *   Considering the above description about the size of each existing binary
  *   variable, the inequalities (2) show a std::vector< FRowConstraint > with
  *   exactly (f_time_horizon - init_t - f_MinUpTime) entries, the entry
  *   a = init_t + f_MinUpTime, ...,(f_time_horizon - 1) being the startup
  *   constraints at time t. when unit in time t is OFF (\f$ u_t = 0 \f$), it
  *   could not have been turned on in the last \f$ \tau_+ \f$ periods
  *   (including period t) because of the minimum up constraints. But this is
  *   exactly what the turn on inequalities (2) for time period t say. On the
  *   other hand, when unit in time t is ON (\f$ u_t = 1 \f$), it could have
  *   been turned on at most once in the last \f$ \tau_+ + \tau_- \f$ periods
  *   (including t).
  *
  *   Similarly for turn off inequalities (3), where it is a
  *   std::vector< FRowConstraint > with exactly
  *   (f_time_horizon - init_t - f_MinDownTime) entries, the entry
  *   a = init_t + f_MinDownTime, ...,(f_time_horizon - 1) being the shut-down
  *   constraints at time t. when unit in the time t is OFF (\f$ u_t = 0
  *   \f$), it could have been turned off at most once in the last
  *   \f$ \tau_+ + \tau_-\f$ periods (including t). On the other hand, when
  *   unit in time t is ON (\f$ u_t = 1 \f$), it could not have been turned
  *   off in the last \f$ \tau_- \f$ periods (including period t).
  *
  * - Ramp Up/Down-time Constraints:
  *
  *   Another set of constraints where each thermal unit may have are ramping
  *   constraints. The ramp-up constraints is a std::vector< FRowConstraint >
  *   with exactly f_time_horizon entries, which are
  *   a = 0, ..., (f_time_horizon - 1). The one possible implementation in
  *   terms of the three binary variables for ramp-up constraints is:
  *
  *   \f[
  *     p_{t+1}^{ac} - p_t^{ac} \leq ( - \Delta^+_t)  v_{t+1}
  *        + (\underline{p}_t + \Delta^+_t) u_{t+1} - \underline{p}_t u_t
  *                       \quad t \in \{ t_0, ..., \mathcal{T} - 1 \} \quad (4)
  *   \f]
  *
  *   where \f$ \Delta^+_t \f$ and \f$ \Delta^-_t \f$ are the constants
  *   defining ramp-up threshold and \f$ \underline{p}_t  \f$ and
  *   \f$ \bar{p}_t \f$  are the defining minimum and maximum output
  *   respectively. Let \f$ p_t^{ac} \f$ be the active power variable in time
  *   period t in all time horizon \f$ \mathcal{T} \f$.
  *
  *   According to above definition about the size of variables and since the
  *   variables commitment \f$ u_t \f$ have full size of time horizon and
  *   start-up \f$ v_t \f$ variables have size (f_time_horizon - init_t).
  *   Analyzing the left hand side of the ramp-up constraint (4), in any
  *   integral feasible solution we can see that
  *   \f$ p_{t+1}^{ac} - p_t^{ac} \f$ can be bounded from above based on the
  *   values of \f$ u_{t+1}\f$, \f$ u_t \f$ and \f$ v_{t+1}\f$. Then for each
  *   (0, ..., f_time_horizon - 1) entries of this
  *   std::vector< FRowConstraint >, there are two possible cases for t from
  *   0 until init_t - 1:
  *
  *   - when \f$ u_t = 0 \f$, and \f$ u_{t+1} = 0 \f$ then
  *     \f$ p_{t+1}^{ac} - p_t^{ac}  \leq 0 \f$.
  *
  *   - when \f$ u_t = 1 \f$, and \f$ u_{t+1} = 1 \f$ then
  *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \Delta^+_t \f$.
  *
  *   and four possible cases for each t from init_t until
  *   \f$ \mathcal{T} - 1 \f$:
  *
  *   - when \f$ u_t = 0 \f$, \f$ u_{t+1} = 0 \f$ and \f$ v_{t+1} = 0 \f$ then
  *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq 0 \f$.
  *
  *   - when \f$ u_t = 0 \f$, \f$ u_{t+1} = 1 \f$ and \f$ v_{t+1} = 1 \f$ then
  *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \underline{p}_t \f$.
  *
  *   - when \f$ u_t = 1 \f$, \f$ u_{t+1} = 0 \f$ and \f$ v_{t+1} = 0 \f$ then
  *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq - \underline{p}_t \f$.
  *
  *   - when \f$ u_t = 1 \f$, \f$ u_{t+1} = 1 \f$ and \f$ v_{t+1} = 0 \f$ then
  *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \Delta^+_t \f$.
  *
  *   Using the symmetry between ramp up and ramp down constraints, we can
  *   derive the ramp-down analogues of the ramp-up inequality as below:
  *
  *   \f[
  *     p_t^{ac} - p_{t+1}^{ac} \leq ( - \Delta^-_t) w_{t+1}
  *       + (\underline{p}_t + \Delta^-_t)  u_t - \underline{p}_t u_{t+1}
  *                        \quad t \in \{t_0, ..., \mathcal{T} - 1 \} \quad (5)
  *   \f]
  *
  *   The sam analyzing the left hand side of the ramp-down constraint (5), in
  *   any integral feasible solution we can see that
  *   \f$ p_t^{ac} - p_{t+1}^{ac} \f$ can be bounded from above based on the
  *   values of \f$ u_{t+1}\f$, \f$ u_t \f$ and \f$ w_{t+1}\f$. Then for each
  *   (0, ..., f_time_horizon - 1) entries of this std::vector< FRowConstraint >
  *   there are two possible cases for t from 0 until init_t - 1:
  *
  *   - when \f$ u_t = 0 \f$, and \f$ u_{t+1} = 0 \f$ then
  *     \f$ p_t^{ac} - p_{t+1}^{ac}  \leq 0 \f$.
  *
  *   - when \f$ u_t = 1 \f$, and \f$ u_{t+1} = 1 \f$ then
  *     \f$ p_t^{ac} - p_{t+}^{ac} \leq \Delta^-_t \f$.
  *
  *   and four possible cases for each t from init_t until
  *   \f$ \mathcal{T} - 1 \f$:
  *
  *   - when \f$ u_t = 0 \f$, \f$ u_{t+1} = 0 \f$ and \f$ w_{t+1} = 0 \f$ then
  *     \f$ p_t^{ac} - p_{t+1}^{ac}  \leq 0 \f$.
  *
  *   - when \f$ u_t = 0 \f$, \f$ u_{t+1} = 1 \f$ and \f$ w_{t+1} = 0 \f$ then
  *     \f$ p_t^{ac} - p_{t+1}^{ac} \leq \underline{p}_t \f$.
  *
  *   - when \f$ u_t = 1 \f$, \f$ u_{t+1} = 0 \f$ and \f$ w_{t+1} = 1 \f$ then
  *     \f$ p_t^{ac} - p_{t+1}^{ac} \leq - \underline{p}_t \f$.
  *
  *   - when \f$ u_t = 1 \f$, \f$ u_{t+1} = 1 \f$ and \f$ w_{t+1} = 0 \f$ then
  *     \f$ p_t^{ac} - p_{t+1}^{ac} \leq \Delta^-_t \f$.
  *
  * - Power output Constraints:
  *   Since commitment variable \f$ u_t \f$ is fixed to one or zero
  *   for "init_t" time steps (look above comments), because of power output
  *   constraint (look constraint (6)) when for the (0, ..., init_t - 1) time
  *   steps, commitment variable \f$ u_t \f$ is fixed to zero we must fix
  *   \f$ p_t^{ac} \f$, \f$ p_t^{pr} \f$, and \f$ p_t^{sc}\f$ to zero for the same
  *   time steps.
  *
  *   Maximum and minimum power output constraints according to active power,
  *   primary and secondary spinning reserves variables are presented in
  *   inequalities (6) and (7) respectively. Each of them is a
  *   std::vector< FRowConstraint > with exactly f_time_horizon entries
  *   (0, ..., (f_time_horizon) - 1) and ensures the maximum (or minimum) amount
  *   of energy that unit can produce (or use) when it is on (or off).
  *
  *   \f[
  *      p_t^{ac} + p_t^{pr} + p_t^{sc} \leq \bar{p}_t u_t            \quad (6)
  *   \f]
  *
  *   \f[
  *     \underline{p}_t u_t \leq p_t^{ac} - p_t^{pr} - p_t^{sc}       \quad (7)
  *   \f]
  *
  *   The same as inequalities (6)-(7), the inequalities (8)-(9) ensure that
  *   maximum amount of primary and secondary spinning reserve in the problem
  *   respectively. Each of them is a std::vector< FRowConstraint > with exactly
  *   f_time_horizon entries (0, ..., (f_time_horizon) - 1) as below:
  *
  *   \f[
  *     p_t^{pr} \leq \rho^{pr}_t p_t^{ac}                            \quad (8)
  *   \f]
  *
  *   \f[
  *     p_t^{sc} \leq \rho^{sc}_t p_t^{ac}                            \quad (9)
  *   \f]
  *
  *   There are two more power out put tighter formulations which make the
  *   maximum power output being a function of three binary variables
  *   \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ as below. More specifically in
  *   the case  \f$ \tau_+ \geq 2 \f$, the
  *   following constraint is introduced, which is valid for
  *   \f$ t \in \{t_0 + 2, ..., \mathcal{T} - 1\}  \f$:
  *
  *   \f[
  *     p_t^{ac} \leq \bar{p}_t  u_t  - ( \bar{p}_t - \underline{p}_t ) v_t
  *                   - ( \bar{p}_t - \underline{p}_t ) w_{t+1}
  *            \quad t \in \{t_0 + 2, ..., \mathcal{T} - 1\}         \quad (10)
  *   \f]
  *
  *   and in the case  \f$ \tau_+ = 1 \f$:
  *
  *   \f[
  *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) w_{t+1}
  *            \quad t \in \{t_0 + 2, ..., \mathcal{T} - 1\}         \quad (11)
  *   \f]
  *
  *   \f[
  *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) v_t
  *             \quad t \in \{t_0 + 2, ...,  \mathcal{T} - 1\}       \quad (12)
  *   \f]
  *
  *   These inequalities give the active power output generation limits when
  *   unit is ON or OFF. More precisely, the unit generation limits taking into
  *   account its maximum \f$ \bar{p}_t \f$ and minimum \f$ \underline{p}_t \f$
  *   production, as well as its start-up and shut-down capabilities (here
  *   both of them are assumed be equal with minimum production \f$
  *   \underline{p}_t \f$) in each time step t. Be aware that (10) may be
  *   infeasible in the event that the unit is online for just one period.
  *   That is, \f$ v_t = w_{t+1} = 1 \f$ and the right side of the (10) can
  *   be negative. Consequently, (10) is only valid when \f$ \tau_+ \geq 2
  *   \f$. Therefore, the correct formulation for units with \f$ \tau_+ = 1
  *   \f$ is given by (11) and (12). */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the dynamic constraint of the ThermalUnitBlock
 /**
  * TODO - following the code-flow / order, add in the method blow the cnstrs
  * TODO - details for each formulation.
  */

 void generate_dynamic_constraints( Configuration * dycc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the ThermalUnitBlock
 /** Method that generates the objective of the ThermalUnitBlock. The objective
  * function of the ThermalUnitBlock representing the total power production
  * cost to be minimized has the form:
  *
  * \f[
  *   \min ( \sum_{ t \in  [t_0 , \mathcal{|T|} - 1] } s_t v_t +
  *   \sum_{ t \in \mathcal{T}  } (a_t p_t^2 + b_t p_t + c_t u_t) )
  * \f]
  *
  * where \f$ v_t \f$ indicates that the unit is starting up at time
  * \f$ t \f$, \f$ u_t \f$ indicates that the unit is committed at time
  * \f$ t \f$, \f$ p_t \f$ is the active power produced at time \f$ t \f$,
  * \f$ \sum_{ t \in [t_0, \mathcal{|T|} - 1] } s_t v_t \f$ is the start-up
  * cost of the unit, which we assume to be time-independent
  *
  * Note: time-independent means here start-up cost is "independent from how
  *       long the unit has been off", and it is not meaning "always should
  *       be equal at each time instant"
  *
  * and \f$ a_t \f$, \f$ b_t \f$, and \f$ c_t \f$ are, respectively, the
  * quadratic, linear, and constant terms of the power cost function of the
  * unit at time period \f$ t \in \mathcal{T} \f$.
  *
  * If the primary and/or the secondary spinning reserve variables have been
  * generated, it is also possible to consider them in linear form in the
  * objective function. This can be instructed by using either the parameter
  * \p objc or f_BlockConfig->f_objective_Configuration. If \p objc is not
  * nullptr and it is a SimpleConfiguration< int >, or if
  * f_BlockConfig->f_objective_Configuration is not nullptr and it is a
  * SimpleConfiguration< int >, then the f_value of this SimpleConfiguration
  * (an int) indicates whether the primary and/or the secondary reserve
  * variables should be included in the objective function. If the
  * Configuration is not available, the default value is taken to be 0. If
  * the first bit of this int value is 1, then the primary spinning reserve
  * variables are added to the objective function, i.e., the following term
  * is added to the objective function described above:
  *
  * \f[
  *    \sum_{ t \in \mathcal{T} } c_t^{pr} p_t^{pr}
  * \f]
  *
  * If the second bit of this int value is 1, then the secondary spinning
  * reserve variables are added to the objective function, i.e., the following
  * term is added to the objective function described above:
  *
  * \f[
  *    \sum_{ t \in \mathcal{T} } c_t^{sc} p_t^{sc}
  * \f]
  *
  * If the primary and/or secondary spinning reserve variables are included in
  * the objective function, their coefficients can be set by the
  * set_primary_spinning_reserve_cost() and
  * set_secondary_spinning_reserve_cost() methods, respectively.
  * In the the design scenario of the UC problem, an additional cost is
  * added to the objective, i.e.:
  *
  * \f[
  *    I x
  * \f]
  *
  * where \f$ I \f$ is the investment cost and \f$ x \f$ is the design binary
  * variable. */

 void generate_objective( Configuration * objc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// ignore reserve netCDF variables when a ThermalUnitBlock is deserialized
 /** This function instructs the ThermalUnitBlock to ignore the reserve netCDF
  * variables, namely "PrimaryRho" and "SecondaryRho", when it is
  * deserialized. */

 static void ignore_reserve() {
  f_ignore_netcdf_vars |= 1;
 }

/**@} ----------------------------------------------------------------------*/
/*--------------- Methods for checking the ThermalUnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the ThermalUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this ThermalUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this ThermalUnitBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this ThermalUnitBlock is a RowConstraint and its
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

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE ThermalUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ThermalUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of electrical generation units, i.e.:
 *
 * - fixed consumption when the unit is off;
 *
 * - the contribution to the inertia depending on the commitment status.
 * @{ */

 /// returns the initial power value
 double get_initial_power( void ) const { return( f_InitialPower ); }

 /// returns the init up and down time value
 int get_init_up_down_time( void ) const { return( f_InitUpDownTime ); }

 /// returns the minimum allowed up time value
 Index get_min_up_time( void ) const { return( f_MinUpTime ); }

 /// returns the minimum allowed down time value
 Index get_min_down_time( void ) const { return( f_MinDownTime ); }

 /// returns the investment cost
 double get_investment_cost( void ) const { return( f_InvestmentCost ); }

 /// returns the installable capacity by the user
 double get_capacity( void ) const { return( f_Capacity ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nominal minimum active power output
 /** This method returns (a const reference to) the vector containing the
  * nominal minimum active power output of the unit for all time steps. When
  * the unit is available, get_min_power()[ t ] gives the minimum active power
  * output of the unit at time t, for each t in {0, ..., get_time_horizon() -
  * 1}. */

 const std::vector< double > & get_min_power( void ) const {
  return( v_MinPower );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power of the given generator at the given time

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( v_MinPower[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the operational minimum active power output at the given time
 /** This method returns the operational minimum active power output of the
  * unit at the given time. See get_availability() for the definition of
  * operational minimum power.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The operational minimum active power output of the unit at the
  *         given time.
  */
 double get_operational_min_power( Index t ) const {
  assert( t < get_time_horizon() );
  return( compute_operational_min_power( v_MinPower[ t ] ,
                                         get_availability( t ) ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nominal maximum active power output
 /** This method returns (a const reference to) the vector containing the
  * nominal maximum active power output of the unit for all time steps. When
  * the unit is fully available, get_min_power()[ t ] gives the maximum active
  * power output of the unit at time t, for each t in {0, ...,
  * get_time_horizon() - 1}. See get_availability() to understand the
  * difference between nominal and operational maximum active power. */

 const std::vector< double > & get_max_power( void ) const {
  return( v_MaxPower );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power of the given generator at the given time

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the operational maximum active power output at the given time
 /** This method returns the operational maximum active power output of the
  * unit at the given time. See get_availability() for the definition of
  * operational maximum power.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The operational maximum active power output of the unit at the
  *         given time.
  */
 double get_operational_max_power( Index t ) const {
  assert( t < get_time_horizon() );
  return( compute_operational_max_power( v_MaxPower[ t ] ,
                                         get_availability( t ) ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the availability of the unit at all time instants
 /** This method returns (a const reference to) the vector containing the
  * availability of the unit at all time instants. For each t in {0, ...,
  * get_time_horizon() - 1}, get_availability()[ t ] is the availability of
  * the unit at time t, which is a number between 0 and 1. When the
  * availability of the unit is zero, the unit is not under operation (for
  * instance, due to an outage or maintenance). When the availability of the
  * unit is 1, it is fully available and operating at maximum capacity.
  *
  * The availability of the unit determines its operational minimum and
  * maximum active power output, i.e., the effective bounds on the active
  * power output under which the unit operates. For each t in {0, ...,
  * get_time_horizon() - 1}, let MinPower[ t ] and MaxPower[ t ] be the
  * nominal minimum and maximum active power output of the unit (given by
  * get_min_power() and get_max_power(), respectively) and let AvMinPower[ t ]
  * and AvMaxPower[ t ] be the operational minimum and maximum active power
  * output of the unit, which depends on its availability. Then,
  *
  *   AvMaxPower[ t ] = Availability[ t ] * MaxPower[ t ]
  *
  * and
  *
  *   AvMinPower[ t ] = MinPower[ t ] if Availability[ t ] > 0, and
  *
  *   AvMinPower[ t ] = 0 if Availability[ t ] = 0,
  *
  * where Availability[ t ] denotes the availability of the unit at time t. */

 const std::vector< double > & get_availability( void ) const {
  return( v_Availability );
 }

/*--------------------------------------------------------------------------*/
 /// returns the availability of the unit at the given time

 double get_availability( Index t ) const {
  if( v_Availability.empty() )
   return( 1.0 );
  assert( t < get_time_horizon() );
  return( v_Availability[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary rho
 /** The returned vector contains the primary rho at each time.
  * The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_primary_rho( void ) const {
  return( v_PrimaryRho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary rho
 /** The returned vector contains the secondary rho at each time.
  * The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_secondary_rho( void ) const {
  return( v_SecondaryRho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary spinning reserve costs
 /** This function returns the vector of primary spinning reserve costs. If it
  * is empty, then the costs are all zero. Otherwise, it has size
  * get_time_horizon() and its t-th element is the linear cost of the primary
  * spinning reserve variable at time t.
  *
  * @return The vector containing the primary spinning reserve costs. */

 const std::vector< double > & get_primary_spinning_reserve_cost( void )
 const { return( v_PrimarySpinningReserveCost ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary spinning reserve costs
 /** This function returns the vector of secondary spinning reserve costs. If
  * it is empty, then the costs are all zero. Otherwise, it has size
  * get_time_horizon() and its t-th element is the linear cost of the
  * secondary spinning reserve variable at time t.
  *
  * @return The vector containing the secondary spinning reserve costs. */

 const std::vector< double > & get_secondary_spinning_reserve_cost( void )
 const { return( v_SecondarySpinningReserveCost ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of delta ramp-up
 /** The returned vector contains the delta ramp-up at each time.
  * The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_delta_ramp_up( void ) const {
  return( v_DeltaRampUp );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of delta ramp-down
 /** The returned vector contains the delta ramp-up at each time.
  * The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_delta_ramp_down( void ) const {
  return( v_DeltaRampDown );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of quadratic term
 /** The returned vector contains to quadratic term at time t. There are three
  * possible cases:
  *
  * - if the vector is empty, then the quadratic term of the unit is 0;
  *
  * - if the vector has only one element, then the quadratic term of the unit
  *   for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the amount of quadratic term at time t. */

 const std::vector< double > & get_quad_term( void ) const {
  return( v_QuadTerm );
 }

/*--------------------------------------------------------------------------*/
/// returns the coefficient of the quadratic term of the power cost function
/** This function returns the coefficient of the quadratic term of the
 * quadratic function that represents the cost of the power produced by the
 * unit at the given time instant.
 *
 * @param t A time instant between 0 and get_time_horizon() - 1.
 *
 * @return The coefficient of the quadratic term of the quadratic function
 *         that represents the cost of the power produced by the unit at the
 *         given time instant. */

 double get_quad_term( Index t ) const {
  if( v_QuadTerm.empty() )
   return( 0 );
  if( v_QuadTerm.size() == 1 )
   return( v_QuadTerm.front() );
  assert( v_QuadTerm.size() == f_time_horizon );
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_quad_term: Invalid "
                            "time index: " + std::to_string( t ) ) );
  return( v_QuadTerm[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of linear term
 /** The returned vector contains to linear term at time t. There are three
  * possible cases:
  *
  * - if the vector is empty, then the linear term of the unit is 0;
  *
  * - if the vector has only one element, then the linear term of the unit for
  *   all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the amount of linear term at time t. */

 const std::vector< double > & get_linear_term( void ) const {
  return( v_LinearTerm );
 }

/*--------------------------------------------------------------------------*/
/// returns the coefficient of the linear term of the power cost function
/** This function returns the coefficient of the linear term of the quadratic
 * function that represents the cost of the power produced by the unit at the
 * given time instant.
 *
 * @param t A time instant between 0 and get_time_horizon() - 1.
 *
 * @return The coefficient of the linear term of the quadratic function that
 *         represents the cost of the power produced by the unit at the given
 *         time instant. */

 double get_linear_term( Index t ) const {
  if( v_LinearTerm.empty() )
   return( 0 );
  if( v_LinearTerm.size() == 1 )
   return( v_LinearTerm.front() );
  assert( v_LinearTerm.size() == f_time_horizon );
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_linear_term: Invalid "
                            "time index: " + std::to_string( t ) ) );
  return( v_LinearTerm[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of constant term
 /** The returned vector contains to constant term at time t. There are three
  * possible cases:
  *
  * - if the vector is empty, then the constant term of the unit is 0;
  *
  * - if the vector has only one element, then the constant term of the unit
  *   for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the amount of constant term at time t. */

 const std::vector< double > & get_const_term( void ) const {
  return( v_ConstTerm );
 }

/*--------------------------------------------------------------------------*/
/// returns the constant term of the power cost function
/** This function returns the constant term of the function that represents
 * the cost of the power produced by the unit at the given time instant. This
 * is the fixed cost incurred when the unit is committed at time instant \p t.
 *
 * @param t A time instant between 0 and get_time_horizon() - 1.
 *
 * @return The fixed cost when the unit is committed at the given time
 *         instant. */

 double get_const_term( Index t ) const {
  if( v_ConstTerm.empty() )
   return( 0 );
  if( v_ConstTerm.size() == 1 )
   return( v_ConstTerm.front() );
  assert( v_ConstTerm.size() == f_time_horizon );
  if( t >= f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::get_const_term: Invalid "
                            "time index: " + std::to_string( t ) ) );
  return( v_ConstTerm[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start-up cost
 /** The returned vector contains to start-up cost at time t. There are three
  * possible cases:
  *
  * - if the vector is empty, then the start-up cost of the unit is 0;
  *
  * - if the vector has only one element, then the start-up cost of the unit
  *   for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the start-up cost value at time t. */

 const std::vector< double > & get_start_up_cost( void ) const {
  return( v_StartUpCost );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of fixed consumption
 /** The returned value U = get_fixed_consumption() contains the contribution
  * to fixed consumption (basically, the constants to be multiplied by the
  * commitment variables returned by get_commitment()) of all the generators
  * at all time instants. There are three possible cases:
  *
  * - if the vector is empty, then the fixed consumption is always 0 and this
  *   function returns nullptr;
  *
  * - if the vector only has one element, then the fixed consumption for the
  *   fixed consumption of the unit for all t;
  *
  * - otherwise, the vector must have size get_time_horizon(), and each element
  *   of vector represents the fixed consumption at time t. */

 const double * get_fixed_consumption( Index generator ) const override {
  if( v_FixedConsumption.empty() )
   return( nullptr );
  return( &( v_FixedConsumption.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inertia commitment
 /** The returned value U = get_inertia_commitment() contains the contribution
  * to inertia (basically, the constants to be multiplied by the commitment
  * variables returned by get_commitment()) of all the generators at all time
  * instants. There are three possible cases:
  *
  * - if the vector is empty, then the inertia commitment is always 0 and this
  *   functions returns nullptr;
  *
  * - if the vector only has one element, then the inertia commitment for the
  *   fixed consumption of the unit for all t;
  *
  * - otherwise, the vector must have size get_time_horizon(), and each element
  *   of vector represents the inertia commitment at time t. */

 const double * get_inertia_commitment( Index generator ) const override {
  if( v_InertiaCommitment.empty() )
   return( nullptr );
  return( &( v_InertiaCommitment.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the start-up limit
 const std::vector< double > & get_start_up_limit( void ) const {
  return( v_StartUpLimit );
 }

/*--------------------------------------------------------------------------*/
 /// returns the shut-down limit
 const std::vector< double > & get_shut_down_limit( void ) const {
  return( v_ShutDownLimit );
 }

/*--------------------------------------------------------------------------*/
 /// returns the scale factor
 double get_scale( void ) const override { return( f_scale ); }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Variable OF THE ThermalUnitBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the ThermalUnitBlock
 *
 * These methods allow to read the each group of Variable that any
 * ThermalUnitBlock in principle has (although some may not):
 *
 * - commitment variables;
 *
 * - active_power variables;
 *
 * - primary_spinning_reserve variables;
 *
 * - secondary_spinning_reserve variables;
 *
 * - start_up variables;
 *
 * - shut_down variables.
 * @{ */

 /// returns the vector of commitment variables
 ColVariable * get_commitment( Index generator ) override {
  if( v_commitment.empty() )
   return( nullptr );
  return( &( v_commitment.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active_power variables
 ColVariable * get_active_power( Index generator ) override {
  if( v_active_power.empty() )
   return( nullptr );
  return( &( v_active_power.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary_spinning_reserve variables
 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( v_primary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_primary_spinning_reserve.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary_spinning_reserve variables
 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( v_secondary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_secondary_spinning_reserve.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start_up variables, or nullptr if not defined
 ColVariable * get_start_up( void ) {
  if( v_start_up.empty() )
   return( nullptr );
  return( &( v_start_up.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of shut_down variables, or nullptr if not defined
 ColVariable * get_shut_down( void ) {
  if( v_shut_down.empty() )
   return( nullptr );
  return( &( v_shut_down.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the shut_down variable for time t, or nullptr if not defined
 ColVariable * get_shut_down( Index t ) {
  if( v_shut_down.empty() || ( t < init_t ) )
   return( nullptr );
  return( &( v_shut_down[ t - init_t ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the design binary variable

 ColVariable & get_design( void ) {
  return( design );
 }

/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE ThermalUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the ThermalUnitBlock
 * @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * ThermalUnitBlock. See ThermalUnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR INITIALIZING THE ThermalUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the ThermalUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "ThermalUnitBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
 /** Method for handling Modification.
  *
  * This method has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the ThermalUnitBlock, and
  * "translate" them into both changes of the actual data structures and
  * corresponding "physical Modification". These Modification are those
  * for which Modification::concerns_Block() is true.
  *
  *     THE IMPLEMENTATION OF THIS METHOD IS BOTH PARTIAL AND HORRIBLE,
  *     ONE SINGLE ABSTRACT MODIFICATION CAN GIVE RISE TO MANY MANY MANY
  *     PHYSICAL ONES, IT SHOULD BE COMPLETELY OVERHAULED!!!
  */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// update the availability of the unit
 /** This method updates the availability of the unit. The \p subset parameter
  * contains a list of time instants and \p values contains the availability
  * of the unit at those time instants. The availability of the unit at time
  * subset[ i ] is given by std::next( values , i ) for each i in {0, ...,
  * subset.size() - 1}.
  *
  * Let AvMinPower[ t ] and AvMaxPower[ t ] denote the operational minimum and
  * maximum active power of the unit at time t. Then, the following condition
  * must be satisfied:
  *
  *   AvMinPower[ t ] <= AvMaxPower[ t ]
  *
  * for each t in {0, ..., get_time_horizon() - 1} (see get_availability() for
  * the definition of operational maximum and minimum active power). If the
  * given availability in \p values is such that this condition does not hold,
  * an exception is thrown. */

 void set_availability( MF_dbl_it values ,
                        Subset && subset ,
                        const bool ordered = false ,
                        ModParam issuePMod = eNoBlck ,
                        ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update the availability of the unit
 /** This method updates the availability of the unit. The \p rng parameter
  * contains a range of time instants and \p values contains the availability
  * of the unit at those time instants. The availability of the unit at time
  * rng.first + i is given by std::next( values , i ) for each i in {0, ...,
  * ( std::min( rng.second, get_time_horizon() ) - rng.first - 1 )}.
  *
  * Let AvMinPower[ t ] and AvMaxPower[ t ] denote the operational minimum and
  * maximum active power of the unit at time t. Then, the following condition
  * must be satisfied:
  *
  *   AvMinPower[ t ] <= AvMaxPower[ t ]
  *
  * for each t in {0, ..., get_time_horizon() - 1} (see get_availability() for
  * the definition of operational maximum and minimum active power). If the
  * given availability in \p values is such that this condition does not hold,
  * an exception is thrown. */

 void set_availability( MF_dbl_it values , Range rng = INFRange ,
                        ModParam issuePMod = eNoBlck ,
                        ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_maximum_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_maximum_power( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_startup_costs( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_startup_costs( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_const_term( MF_dbl_it values ,
                      Subset && subset ,
                      const bool ordered = false ,
                      ModParam issuePMod = eNoBlck ,
                      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_const_term( MF_dbl_it values , Range rng = INFRange ,
                      ModParam issuePMod = eNoBlck ,
                      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_linear_term( MF_dbl_it values ,
                       Subset && subset ,
                       const bool ordered = false ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_linear_term( MF_dbl_it values , Range rng = INFRange ,
                       ModParam issuePMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_quad_term( MF_dbl_it values ,
                     Subset && subset ,
                     const bool ordered = false ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_quad_term( MF_dbl_it values , Range rng = INFRange ,
                     ModParam issuePMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                         Subset && subset ,
                                         const bool ordered = false ,
                                         ModParam issuePMod = eNoBlck ,
                                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                         Range rng = INFRange ,
                                         ModParam issuePMod = eNoBlck ,
                                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                           Subset && subset ,
                                           const bool ordered = false ,
                                           ModParam issuePMod = eNoBlck ,
                                           ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                           Range rng = INFRange ,
                                           ModParam issuePMod = eNoBlck ,
                                           ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given \p subset contains the 0 index, this function sets the
  * initial power. If the given \p subset does not contain the index 0, this
  * function does nothing. Since \p subset can have multiple zeros, only the
  * last one is considered, which means that the value for the initial power
  * will be that in the vector pointed by \p it associated with this last
  * zero. */

 void set_initial_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given Range \p rng contains 0, this function sets the initial
  * power. In this case, if the first element of \p rng is 0, the initial
  * power will be set to the value pointed by the given iterator. In general,
  * the initial power will be the one found at position -rng.first in the
  * vector pointed by \p it if this Range contains the 0 index. If the given
  * Range \p rng does not contain the 0 index, this function does nothing. */

 void set_initial_power( MF_dbl_it values , Range rng = INFRange ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_init_updown_time( MF_int_it values ,
                            Subset && subset ,
                            const bool ordered = false ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_init_updown_time( MF_int_it values , Range rng = INFRange ,
                            ModParam issuePMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the scale factor
 /** This method sets the scale factor.
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

 // For the Range version, use the default implementation defined in UnitBlock
 using UnitBlock::scale;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// updates the abstract representation dependent on the availability
 /** This method updates any part of the abstract representation that may
  * depend on the availability of the unit at the given time \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 virtual void update_availability_dependents( Index t , c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current initial power
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints at time 0 (which are the
  * constraints that depend on the initial power). */

 virtual void update_initial_power_in_cnstrs( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// returns true if and only if the given availability is consistent
 /** This method checks whether the given \p availability is consistent at
  * time \p t. An availability is consistent at a given time instant if the
  * resulting operational minimum active power output is less than or equal to
  * the resulting operational maximum active power at that time.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @param availability A number between 0 and 1.
  *
  * @return True if and only if the given \p availability is consistent at
  *         time \p t. */

 bool availability_is_consistent( Index t , double availability ) const {
  assert( t < get_time_horizon() );
  const auto min_power = compute_operational_min_power( v_MinPower[ t ] ,
                                                        availability );
  const auto max_power = compute_operational_max_power( v_MaxPower[ t ] ,
                                                        availability );
  return( min_power <= max_power );
 }

/*--------------------------------------------------------------------------*/
 /// returns the operational minimum power
 /** This method computes the operational minimum power for the given nominal
  * minimum power and availability.
  *
  * @param min_power The nominal minimum power.
  *
  * @param availability A number between 0 and 1.
  *
  * @return The operational minimum power. */

 double compute_operational_min_power( double nominal_min_power ,
                                       double availability ) const {
  return( availability > 0.0 ? nominal_min_power : 0.0 );
 }

/*--------------------------------------------------------------------------*/
 /// returns the operational maximum power
 /** This method computes the operational maximum power for the given nominal
  * maximum power and availability.
  *
  * @param min_power The nominal maximum power.
  *
  * @param availability A number between 0 and 1.
  *
  * @return The operational maximum power. */

 double compute_operational_max_power( double nominal_max_power ,
                                       double availability ) const {
  return( nominal_max_power * availability );
 }

/*--------------------------------------------------------------------------*/
 /// updates the terms of the Objective associated with the start-up cost
 /** This method updates the terms of the Objective that are associated with
  * the start-up cost.
  *
  * @param subset A set of time instants at which the start-up costs must be
  *        updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_start_up( const Subset & subset , c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the terms of the Objective associated with the active power cost
 /** This method updates the terms of the Objective that are associated with
  * the active power cost.
  *
  * @param subset A set of time instants at which the active power costs must
  *        be updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_active_power( const Subset & subset ,
                                     c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the terms of the Objective associated with the fixed cost
 /** This method updates the terms of the Objective that are associated with
  * the fixed cost.
  *
  * @param subset A set of time instants at which the fixed costs must be
  *        updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective_commitment( const Subset & subset ,
                                   c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the coefficients of the Objective
 /** This method updates the coefficients of the Objective.
  *
  * @param subset A set of time instants at which the coefficients must be
  *        updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective( const Subset & subset , c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the coefficients of the Objective
 /** This method updates the coefficients of the Objective.
  *
  * @param subset A set of time instants at which the coefficients must be
  *        updated.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective( Range rng , c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this ThermalUnitBlock is consistent
 /** This function checks whether the data in this ThermalUnitBlock is
  * consistent. The data is consistent if all of the following conditions are
  * met.
  *
  * - The minimum power is not greater than the maximum power.
  *
  * - The availability is between 0 and 1.
  *
  * - The delta ramp-up and ramp-down are nonnegative.
  *
  * - The quadratic term of the objective function is nonnegative.
  *
  * If any of the above conditions are not met, an exception is thrown. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );

 void handle_objective_change( FunctionMod * mod , ChnlName chnl );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the vector of MinPower
 std::vector< double > v_MinPower;

 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the vector of Availability
 std::vector< double > v_Availability;

 /// the vector of PrimaryRho
 std::vector< double > v_PrimaryRho;

 /// the vector of SecondaryRho
 std::vector< double > v_SecondaryRho;

 /// the vector of RampUp
 std::vector< double > v_DeltaRampUp;

 /// the vector of RampDown
 std::vector< double > v_DeltaRampDown;

 /// the vector of QuadTerm
 std::vector< double > v_QuadTerm;

 /// the vector of LinearTerm
 std::vector< double > v_LinearTerm;

 /// the vector of ConstTerm
 std::vector< double > v_ConstTerm;

 /// the vector of StartUpCost
 std::vector< double > v_StartUpCost;

 /// the vector of primary spinning reserve linear costs
 std::vector< double > v_PrimarySpinningReserveCost;

 /// the vector of secondary spinning reserve linear costs
 std::vector< double > v_SecondarySpinningReserveCost;

 /// the vector of fixed consumption of generator
 std::vector< double > v_FixedConsumption;

 /// the vector of inertia commitment of generator
 std::vector< double > v_InertiaCommitment;

 /// the vector of start-up limits
 std::vector< double > v_StartUpLimit;

 /// the vector of shut-down limits
 std::vector< double > v_ShutDownLimit;


 /// the vector of index to map the active power variables of the 3bin
 /// formulation with those of the DP formulation  // TODO check
 std::vector< std::pair< Index , std::pair< Index , Index > > > v_P_h_k;

 /// the vector of index to map the perspective cuts variables of the 3bin
 /// formulation with those of the DP formulation  // TODO check
 std::vector< std::pair< Index , std::pair< Index , Index > > > v_Z_h_k;

 /// the vector of index to map the active power variables of the 3bin
 /// formulation with those of the SU formulation  // TODO check
 std::vector< std::pair< Index , Index > > v_P_h;

 /// the vector of index to map the perspective cuts variables of the 3bin
 /// formulation with those of the SU formulation  // TODO check
 std::vector< std::pair< Index , Index > > v_Z_h;

 /// the vector of index to map the active power variables of the 3bin
 /// formulation with those of the SD formulation  // TODO check
 std::vector< std::pair< Index , Index > > v_P_k;

 /// the vector of index to map the perspective cuts variables of the 3bin
 /// formulation with those of the SD formulation  // TODO check
 std::vector< std::pair< Index , Index > > v_Z_k;

 /// TODO short description
 std::vector< Index > v_nodes_plus;

 /// TODO short description
 std::vector< Index > v_nodes_minus;

 /// the vector of index to map the commitment variables of the 3bin
 /// formulation with those of the SD, DP and SU formulations  // TODO check
 std::vector< std::pair< Index , Index > > v_Y_plus;
 std::vector< std::pair< Index , Index > > v_Y_minus;


 /// the investment cost
 double f_InvestmentCost{};

 /// the installable capacity by the user
 double f_Capacity{};

 /// the InitialPower value
 double f_InitialPower{};

 /// the InitUpDownTime value
 int f_InitUpDownTime{};

 /// the MinUpTime value
 Index f_MinUpTime = 1;

 /// the MinDownTime value
 Index f_MinDownTime = 1;

 /// variable denoting the time-steps unit is subjected to initial conditions
 Index init_t{};

 /// the scale factor
 double f_scale = 1;

 /// this variable indicates which netCDF variables must be ignored
 static bool f_ignore_netcdf_vars;

/*-------------------------------- variables -------------------------------*/

 /// the design binary variable
 ColVariable design;

 /// the start-up binary variables
 std::vector< ColVariable > v_start_up;

 /// the shut-down binary variables
 std::vector< ColVariable > v_shut_down;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;


 /// the commitment binary variables for 3bin, T and pt formulations
 std::vector< ColVariable > v_commitment;

 /// the y^+ commitment binary variables for DP, SU and SD formulations
 std::vector< ColVariable > v_commitment_plus;

 /// the y^- commitment binary variables for DP, SU and SD formulations
 std::vector< ColVariable > v_commitment_minus;


 /// the active power variables for 3bin, T and pt formulations
 std::vector< ColVariable > v_active_power;

 /// the active power variables for DP model
 std::vector< ColVariable > v_active_power_h_k;

 /// the active power variables for SU model
 std::vector< ColVariable > v_active_power_h;

 /// the active power variables for SD model
 std::vector< ColVariable > v_active_power_k;


 /// the perspective cuts variables for 3bin, T and pt formulations
 std::vector< ColVariable > v_cut;

 /// the perspective cuts variables for DP model
 std::vector< ColVariable > v_cut_h_k;

 /// the perspective cuts variables for SU model
 std::vector< ColVariable > v_cut_h;

 /// the perspective cuts variables for SD model
 std::vector< ColVariable > v_cut_k;

/*------------------------------- constraints ------------------------------*/

 /// the commitment design constraints
 std::vector< FRowConstraint > CommitmentDesign_Const;

 /// the connection min up and down time constraints
 std::vector< FRowConstraint > StartUp_ShutDown_Variables_Const;

 /// the turn on min up and down time constraints
 std::vector< FRowConstraint > StartUp_Const;

 /// the shut-down min up and down time constraints
 std::vector< FRowConstraint > ShutDown_Const;

 /// the RampUp time constraints
 std::vector< FRowConstraint > RampUp_Const;

 /// the RampDown time constraints
 std::vector< FRowConstraint > RampDown_Const;

 /// the active power upper bound constraints
 std::vector< FRowConstraint > MinPower_Const;

 /// the active power lower bound constraints
 std::vector< FRowConstraint > MaxPower_Const;

 /// the PrimaryRho fraction constraints
 std::vector< FRowConstraint > PrimaryRho_Const;

 /// the SecondaryRho fraction constraints
 std::vector< FRowConstraint > SecondaryRho_Const;


 /// the constraints connecting power variables of 3bin, T and
 /// pt formulations with those of DP, SU and SD formulations
 std::vector< FRowConstraint > Eq_ActivePower_Const;

 /// the constraints connecting commitment variables of 3bin and T
 /// formulations with those of pt, DP, SU and SD formulations
 std::vector< FRowConstraint > Eq_Commitment_Const;

 /// the constraints connecting start-up variables of 3bin and T
 /// formulations with those of pt, DP, SU and SD formulations
 std::vector< FRowConstraint > Eq_StartUp_Const;

 /// the constraints connecting shut-down variables of 3bin and T
 /// formulations with those of pt, DP, SU and SD formulations
 std::vector< FRowConstraint > Eq_ShutDown_Const;

 /// the network constraints of the pt, DP, SU and SD formulations
 std::vector< FRowConstraint > Network_Const;


 /// the initial perspective cuts constraints
 std::vector< FRowConstraint > Init_PC_Const;

 /// the constraints connecting perspective cuts variables of 3bin, T
 /// and pt formulations with those of DP, SU and SD formulations
 std::vector< FRowConstraint > Eq_PC_Const;


 /// the perspective dynamic cuts constraints
 std::list< FRowConstraint > PC_cuts;


 /// the commitment bound constraints
 std::vector< ZOConstraint > Commitment_bound_Const;

 /// the start-up binary bound constraints
 std::vector< ZOConstraint > StartUp_Binary_bound_Const;

 /// the shut-down binary bound constraints
 std::vector< ZOConstraint > ShutDown_Binary_bound_Const;


 /// the commitment fixed to one BoxConstraints
 std::vector< BoxConstraint > Commitment_fixed_to_One_Const;


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

 static void static_initialization( void ) {

  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args(),
   *
   * register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_availability",
   *                                      &ThermalUnitBlock::set_availability,
   *                                      MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_availability" ,
   &ThermalUnitBlock::set_availability );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_availability" ,
   &ThermalUnitBlock::set_availability );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_maximum_power" ,
   &ThermalUnitBlock::set_maximum_power );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_maximum_power" ,
   &ThermalUnitBlock::set_maximum_power );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::set_initial_power" ,
   &ThermalUnitBlock::set_initial_power );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::set_initial_power" ,
   &ThermalUnitBlock::set_initial_power );

  register_method< ThermalUnitBlock , MF_int_it , Subset && , bool >(
   "ThermalUnitBlock::set_init_updown_time" ,
   &ThermalUnitBlock::set_init_updown_time );

  register_method< ThermalUnitBlock , MF_int_it , Range >(
   "ThermalUnitBlock::set_init_updown_time" ,
   &ThermalUnitBlock::set_init_updown_time );

  register_method< ThermalUnitBlock , MF_dbl_it , Subset && , bool >(
   "ThermalUnitBlock::scale" , &ThermalUnitBlock::scale );

  register_method< ThermalUnitBlock , MF_dbl_it , Range >(
   "ThermalUnitBlock::scale" , &ThermalUnitBlock::scale );
 }

};  // end( class( ThermalUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a ThermalUnitBlock
class ThermalUnitBlockMod : public UnitBlockMod
{

 public:

 /// public enum for the types of ThermalUnitBlockMod
 enum TUBB_mod_type
 {
  eSetMaxP = eUBModLastParam , ///< set max power values
  eSetInitP ,                  ///< set initial power values
  eSetInitUD ,                 ///< set initial up/down times
  eSetAv ,                     ///< set availability
  eSetSUC ,                    ///< set start-up costs
  eSetLinT ,                   ///< set linear term
  eSetQuadT ,                  ///< set quad term
  eSetConstT ,                 ///< set constant term
  eSetPrSpResCost ,            ///< set primary spinning reserve (linear) costs
  eSetSecSpResCost ,
  ///< set secondary spinning reserve (linear) costs
  eTUBBModLastParam       ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of ThermalUnitBlockMod. */
 };

 /// constructor, takes the ThermalUnitBlock and the type
 ThermalUnitBlockMod( ThermalUnitBlock * const fblock , const int type )
  : UnitBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the ThermalUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "ThermalUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetMaxP ):
    output << "set max power values";
    break;
   case( eSetInitP ):
    output << "set initial power values";
    break;
   case( eSetInitUD ):
    output << "Set initial up/down times";
    break;
   case( eSetAv ):
    output << "Set availability";
    break;
   case( eSetSUC ):
    output << "Set start-up costs";
    break;
   case( eSetLinT ):
    output << "Set linear term";
    break;
   case( eSetQuadT ):
    output << "Set quad term";
    break;
   case( eSetConstT ):
    output << "Set constant term";
    break;
   default:;
  }
 }

};  // end( class( ThermalUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS ThermalUnitBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ThermalUnitBlockMod for "ranged" modifications
class ThermalUnitBlockRngdMod : public ThermalUnitBlockMod
{

 public:

 /// constructor: takes the ThermalUnitBlock, the type, and the range
 ThermalUnitBlockRngdMod( ThermalUnitBlock * const fblock ,
                          const int type ,
                          Block::Range rng )
  : ThermalUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the ThermalUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  ThermalUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng;  ///< the range

};  // end( class( ThermalUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ThermalUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ThermalUnitBlockMod for "subset" modifications
class ThermalUnitBlockSbstMod : public ThermalUnitBlockMod
{

 public:

 /// constructor: takes the ThermalUnitBlock, the type, and the subset
 ThermalUnitBlockSbstMod( ThermalUnitBlock * const fblock ,
                          const int type ,
                          Block::Subset && nms )
  : ThermalUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the ThermalUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  ThermalUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

};  // end( class( ThermalUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ThermalUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
