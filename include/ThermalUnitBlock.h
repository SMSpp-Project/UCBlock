/*--------------------------------------------------------------------------*/
/*------------------------- File ThermalUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class ThermalUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" thermal unit
 * of a Unit Commitment Problem. A ThermalUnitBlock corresponds to a single
 * electrical generator.
 *
 * \version 0.11
 *
 * \date 03 - 10 - 2020
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, and Rafael
 * Durbano Lobato
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
#include "DQuadFunction.h"
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

class ThermalUnitBlock : public UnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor, takes the father and the time horizon
 /** Constructor of ThermalUnitBlock, taking possibly a pointer of its
  * father Block and the time horizon. */

 explicit ThermalUnitBlock( Block * f_block = nullptr, Index t = 0 ) :
  UnitBlock( f_block ) {}
/*--------------------------------------------------------------------------*/
 /// destructor of ThermalUnitBlock

 virtual ~ThermalUnitBlock() override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the ThermalUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The variable "MinPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
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
 * - The variable "MaxPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector MxP[ t ] that, for
 *   each time instant t, contains the nominal maximum active power output
 *   value of the unit for the corresponding time step.  If "MaxPower" has
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
 * - The variable "DeltaRampUp", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
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
 * - The variable "DeltaRampDown", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector DM[ t ] that, for
 *   each time instant t, contains the ramp-down value of the unit for the
 *   corresponding time step, i.e., the maximum possible decrease of active
 *   power production w.r.t. the power that had been produced in time instant
 *   t - 1, if any. This variable is optional; if it is not provided then it
 *   is assumed that DP[ t ] == MxP[ t ], i.e., the unit can ramp down an
 *   arbitrary amount, i.e., there are no ramp-down constraints. If
 *   "DeltaRampDown" has length 1 then DM[ t ] contains the same value for all
 *   t. Otherwise, DeltaRampDown[ i ] is the fixed value of DM[ t ] for all t
 *   in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "PrimaryRho", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
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
 * - The variable "SecondaryRho", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector SR[ t ] that, for
 *   each time instant t, contains the maximum possible fraction of active
 *   power that can be used as secondary reserve value of the unit for the
 *   corresponding time step. This variable is optional; if it is not provided
 *   then it is assumed that this unit may not be capable of producing any
 *   secondary reserve, which correspond to SR[ t ] == 0 for all t. If
 *   "SecondaryRho" has length 1 then SR[ t ] contains the same value for all
 *   t. Otherwise, SecondaryRho[ i ] is the fixed value of SR[ t ] for all t
 *   in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ] with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1
 *   or "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "QuadTerm", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector A[ t ] that, for
 *   each time instant t, contains the quadratic term of power cost function
 *   of the unit for the corresponding time step.  This variable is optional;
 *   if it is not provided then it is assumed that A[ t ] == 0, i.e., the cost
 *   of the unit is linear in the produced power.  If "QuadTerm" has length 1
 *   then A[ t ] contains the same value for all t.  Otherwise, QuadTerm[ i ]
 *   is the fixed value of A[ t ] for all t in the interval [ ChangeIntervals[
 *   i - 1 ] , ChangeIntervals[ i ] ], with the assumption that
 *   ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "StartUpCost", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector SC[ t ] that, for
 *   each time instant t, contains the start up cost value of the unit for the
 *   corresponding time step. This variable is optional; if it is not provided
 *   then it is assumed that SC[ t ] == 0, i.e., this unit may not have any
 *   start up cost. If "StartUpCost" has length 1 then SC[ t ] contains the
 *   same value for all t. Otherwise, StartUpCost[ i ] is the fixed value of
 *   SC[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
 *   0. If "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then
 *   the mapping clearly does not require "ChangeIntervals", which in fact is
 *   not loaded.
 *
 * - The variable "LinearTerm", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector B[ t ] that, for
 *   each time instant t, contains the linear term of power cost function of
 *   the unit for the corresponding time step.  This variable is optional; if
 *   it is not provided then it is assumed that B[ t ] == 0, i.e., the cost of
 *   the unit has no linear dependence on the produced power (say, only the
 *   quadratic one). If "LinearTerm" has length 1 then A[ t ] contains the
 *   same value for all t. Otherwise, LinearTerm[ i ] is the fixed value of B[
 *   t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
 *   0. If "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then
 *   the mapping clearly does not require "ChangeIntervals", which in fact is
 *   not loaded.
 *
 * - The variable "ConstTerm", of type double and to be either of size 1 or
 *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector C[ t ] that, for
 *   each time instant t, contains the constant term of power cost function of
 *   the unit for the corresponding time step.  This variable is optional; if
 *   it is not provided then it is assumed that C[ t ] == 0, i.e., the cost of
 *   the unit has no fixed term, only those depending (linearly or
 *   quadratically) on the produced power. If "ConstTerm" has length 1 then C[
 *   t ] contains the same value for all t.  Otherwise, ConstTerm[ i ] is the
 *   fixed value of C[ t ] for all t in the interval [ ChangeIntervals[ i - 1
 *   ] , ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1
 *   ] = 0. If "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon"
 *   then the mapping clearly does not require "ChangeIntervals", which in
 *   fact is not loaded.
 *
 * - The scalar variable "InitialPower", of type double and not indexed over
 *   any dimension. This variable indicates the amount of the power that the
 *   unit was producing at time instant -1, i.e., before the start of the
 *   time horizon; this is necessary to compute the ramp-up and ramp-down
 *   constraints. Clearly, it must be that MaxPower >= InitialPower >=
 *   MinPower if the unit was "on" at time instant -1, and it would be ignored
 *   if the unit was "off" at time instant -1. The on/off status of the unit
 *   is also encoded by the scalar variable InitUpDownTime: in particular,
 *   InitUpDownTime > 0 then the unit was on at time instant -1, and therefore
 *   InitialPower >= MinPower must hold, while if InitUpDownTime <= 0 then the
 *   unit was off at time instant -1, and therefore InitialPower is ignored.
 *   In fact, if InitUpDownTime <= 0 then this variable need not be defined
 *   since it is not loaded.
 *
 * - The scalar variable "InitUpDownTime", of type netCDF::NcInt and not
 *   indexed over any dimension and indicates the initial time to generating
 *   the unit.  If InitUpDownTime > 0, this means that the unit has been on
 *   for InitUpDownTime time stamps prior to time stamp 0 (the beginning of
 *   the horizon). If, instead, InitUpDownTime <= 0, this means that the unit
 *   has been off for - InitUpDownTime time stamps prior to time stamp 0; note
 *   that InitUpDownTime == 0 means that the unit has been just shut down at
 *   the end of time instant -1, i.e., the beginning of time instant 0.
 *
 * - The positive scalar variable "MinUpTime", of type netCDF::NcUint and not
 *   indexed over any dimension, which indicates the minimum allowed up time
 *   in this unit. This variable is optional, if it is not provided it is
 *   taken to be MinUpTime == 0, which mean that the unit can shut down in the
 *   very same time stamp in which it starts up.
 *
 * - The positive scalar variable "MinDownTime", of type netCDF::NcUint and
 *   not indexed over any dimension, which indicates the minimum allowed down
 *   time in this unit.This variable is optional, if it is not provided it is
 *   taken to be MinDownTime == 0, which mean that the unit can start up in
 *   the very same time stamp in which it shuts down.
 *
 * - The variable "FixedConsumption", of type double and either indexed over
 *   the dimension "NumberIntervals" (if "NumberIntervals" is not provided,
 *   then this variable can also be indexed over "TimeHorizon"), or having
 *   size 1. This is meant to represent the vector FC[ t ] which, for each
 *   time instant t, contains the fixed consumption of the power plant if it
 *   is OFF at time t. The variable is optional; if it is not defined, FC[ t ]
 *   == 0 for all time instants. If it has size 1, then FC[ t ] ==
 *   FixedConsumption[ 0 ] for all t, regardless to what "NumberIntervals"
 *   says. Otherwise, FixedConsumption[ i ] is the fixed value of FC[ t ] for
 *   all t in the interval [ ChangeIntervals[ i - 1 ], ChangeIntervals[ i ] ],
 *   with the assumption that ChangeIntervals[ - 1 ] = 0.
 *
 * - The variable "InertiaCommitment", of type double and either indexed over
 *   the dimension "NumberIntervals" (if "NumberIntervals" is not provided,
 *   then this variable can also be indexed over "TimeHorizon") or has size 1.
 *   This is meant to represent the vector IC[ t ] which, for each time
 *   instant t, contains the contribution that the unit can give to the
 *   inertia constraint for the sole fact that is is on (basically, the
 *   constant to be multiplied to the commitment variable) at time t. The
 *   variable is optional; if it is not defined, IC[ t ] == 0 for all time
 *   instants. If it has size 1, then IC[ t ] == InertiaCommitment[ 0 ] for
 *   all t, regardless to what "NumberIntervals" says. Otherwise,
 *   InertiaCommitment[ i ] is the fixed value of IC[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the ThermalUnitBlock
/** The ThermalUnitBlock class has six different variables:
 *
 *  - the binary commitment variables which takes the value of 1 if unit is ON
 *    at time instant t and 0 otherwise;
 *
 *  - the primary spinning reserve variables;
 *
 *  - the secondary spinning reserve variables;
 *
 *  - the active power variables.
 *
 *  All of those variables are optional except the active power variables in
 *  the sense that the model may just not have them and whenever a group of
 *  above variables is created, its size will be the time horizon. Moreover,
 *  ThermalUnitBlock is defined more groups of variables as follow:
 *
 *  - the binary variable start_up status of the unit which takes the value of
 *    1 if the unit starts up at time instant t and 0 otherwise;
 *
 *  - the binary variable shut_down status of the unit which takes the value
 *    of 1 if the unit shuts down at time instant t and 0 otherwise;
 *
 *  These two groups of variables have size f_time_horizon - init_t, and
 *  provide the unit commitment problem with a tight 3-binary MIP formulation.
 *  Since these two variables may have shorter size (when init_t > 0), the
 *  commitment variable needs to be fixed to 0 or 1 for the first init_t time
 *  steps 0, ..., init_t - 1 (see initial time step concept in the
 *  generate_abstract_constraints()).
 *
 *  All of these variables are optional, and it is also possible to restrict
 *  which of the subsets are generated with the parameter stvv. If stvv is not
 *  nullptr and it is a SimpleConfiguration<int>, or if
 *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 *  SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 *  of the optional variables should be created. If the Configuration is not
 *  available, the default value is taken to be 0.
 *
 * Note that there may be other formulations (like the DP one), which will
 * possibly be implemented in the future. */

 void generate_abstract_variables( Configuration *stvv ) override;

/*--------------------------------------------------------------------------*/
/// Generate the static constraint of the ThermalUnitBlock
/** This method generates the abstract constraints of the ThermalUnitBlock.
 *
 * The operations of the thermal generating unit are described on a discrete
 * time horizon as dictated by the UnitBlock interface. In this description
 * we indicate it with \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$.
 * Considering three parameters: InitUpDownTime \f$ \tau_0 \f$ which can be a
 * positive or negative(or 0 )integer number and tells for how many time steps
 * before time step 0 the unit was ON(when\f$ \tau_0 > 0 \f$) or OFF (when
 * \f$ \tau_0 < 0 \f$), MaxUpTime \f$ \tau_+ \f$ which is a positive
 * integer number and indicates for how many time steps after time step
 * 0, the unit can remain ON, and MinDownTime \f$ \tau_- \f$ that is also a
 * positive integer number and indicates for how many time steps after time
 * step 0, the unit can remain off). Therefor, the starting-up (shutting-down)
 * of generating of each unit depends on these three parameters. Then, the
 * firs time instant may not be always equal to zero. For this matter the
 * concept of first time instant which is called "init_t" is defined as below:
 *
 * - If \f$ \tau_0 > 0 \f$, this means that the unit has been on for
 *   \f$ \tau_0 \f$ time stamps prior to time stamp 0 (the beginning of the
 *   time horizon):
 *
 *   - if  \f$ \tau_0 \geq \tau_+ \f$ then init_t = 0 ;
 *
 *   - otherwise init_t = \f$ \tau_+ \f$ - \f$ \tau_0 \f$ ;
 *
 * - If, instead, \f$ \tau_0 < 0 \f$, this means that the unit has
 *   been off for \f$ - \tau_0 \f$ time stamps prior to time stamp 0:
 *
 *   - if \f$ - \tau_0  \geq \tau_- \f$  then init_t = 0;
 *
 *   - otherwise init_t = \f$ \tau_- \f$ + \f$ \tau_0 \f$;
 *
 * - If the \f$ \tau_0 == 0\f$ means that the unit has been just shutdown at
 *   the end of time instant -1, i.e., the beginning of time instant 0 and
 *   init_t == 0.
 *
 * - Note: for the entry a = 0, ..., init_t - 1 when commitment variables fix
 *   to 0 then active power variables fix to 0, but when commitment variables
 *   fix to 1 the bound constraints, a std::vector<LB0Constraint> with exactly
 *   init_t entries, the entry a = 0, ..., init_t - 1 being the bound
 *   constraints of the ColVariable corresponding to the active power
 *   production of the unit(put init_t := \f$ t_0 \f$).
 *
 * Then the main thermal unit constraints with three 3 binary variables
 * \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$are are presented as following:
 *
 * - Min Up/Down-time Constraints: a thermal unit may have minimum up and down
 *   time constraints and one possible representation of the constraints could
 *   be as below:
 *
 *   \f[
 *     u_t - u_{t-1} = v_t - w_t
 *          \quad t \in \{ t_0 , ...,\mathcal{T}- 1 \}             \quad (1)
 *   \f]
 *   \f[
 *    \sum_{ s \in [ t - \tau_+  , t ] } v_s \leq
 *           u_t \quad t \in \{ \tau_+ + t_0, ..., \mathcal{T} - 1\}
 *                                                                   \quad (2)
 *   \f]
 *   \f[
 *    \sum_{ s \in [ t - \tau_- , t ] } w_s \leq
 *         1 - u_t \quad t \in \{ \tau_- + t_0, ...,\mathcal{T} - 1\}
 *                                                                   \quad (3)
 *   \f]
 *
 * - since the variables commitment \f$ u_t \f$ have full size of time horizon
 *   and other two remain variables which are startup \f$ v_t \f$ and shutdown
 *   \f$ w_t \f$ have size (f_time_horizon - init_t), the equalities (1) show
 *   a std::vector<FRowConstraint> with exactly (f_time_horizon - init_t)
 *   entries, the entry a = init_t, ...,(f_time_horizon - 1) being the startup
 *   and shutdown connecting constraints at time t. According to the concept
 *   of init_t when \f$ \tau_0 < 0 \f$  and \f$ -\tau_0 < \tau_- \f$ the
 *   commitment variables \f$ u_t \f$  are fixed to zero for init_t time
 *   steps(starting from zero till init_t - 1). When \f$ \tau_0 > 0 \f$  and
 *   \f$ \tau_0 < \tau_+ \f$ the commitment variables \f$ u_t \f$  are fixed
 *   to one for init_t time steps(starting from zero till init_t - 1). Since
 *   \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ are binary variables, we
 *   can ensure (for all periods \f$ t \in \{ t_0, ..., \mathcal{T} -1 \} \f$)
 *   that \f$ v_t = 1\f$ if and only if \f$ u_t = 1\f$ and \f$ u_{t-1} = 0\f$.
 *   It also obvious that \f$ w_t = 1 \f$ if and only if \f$ u_t = 0 \f$ and
 *   \f$ u_{t-1} = 1 \f$. These conditions are satisfied by equality (1).
 *
 *   Considering the above description about the size of each existing binary
 *   variable, the inequalities (2) show a std::vector<FRowConstraint> with
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
 *   std::vector<FRowConstraint> with exactly
 *   (f_time_horizon - init_t - f_MinDownTime) entries, the entry
 *   a = init_t + f_MinDownTime, ...,(f_time_horizon - 1) being the shutdown
 *   constraints at time t. when unit in the time t is OFF (\f$ u_t = 0
 *   \f$), it could have been turned off at most once in the last
 *   \f$ \tau_+ + \tau_-\f$ periods (including t). On the other hand, when
 *   unit in time t is ON (\f$ u_t = 1 \f$), it could not have been turned
 *   off in the last \f$ \tau_- \f$ periods (including period t).
 *
 * - Ramp Up/Down-time Constraints:
 *
 *   Another set of constraints where each thermal unit may have are ramping
 *   constraints. The ramp-up constraints is a std::vector<FRowConstraint>
 *   with exactly f_time_horizon entries, which are
 *   a = 0, ..., (f_time_horizon - 1). The one possible implementation in
 *   terms of the three binary variables for ramp-up constraints is:
 *   \f[
 *     p_{t+1}^{ac} - p_t^{ac} \leq ( - \Delta^+_t)  v_{t+1}
 *        + (\underline{p}_t + \Delta^+_t) u_{t+1} - \underline{p}_t u_t
 *            \quad t \in \{ t_0, ..., \mathcal{T} - 1 \} \quad (4)
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
 *   startup \f$ v_t \f$ variables have size (f_time_horizon - init_t).
 *   Analyzing the left hand side of the ramp-up constraint(4), in any
 *   integral feasible solution we can see that
 *   \f$ p_{t+1}^{ac} - p_t^{ac} \f$ can be bounded from above based on the
 *   values of \f$ u_{t+1}\f$, \f$ u_{t}\f$ and \f$ v_{t+1}\f$. Then for each
 *   (0, ..., f_time_horizon - 1) entries of this std::vector<FRowConstraint>,
 *   there are two possible cases for t from 0 until init_t - 1:
 *
 *   - when \f$ u_{t} = 0\f$, and \f$ u_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac}  \leq 0 \f$.
 *
 *   - when \f$ u_{t} = 1\f$, and \f$ u_{t+1} = 1 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \Delta^+_t \f$.
 *
 *   and four possible cases for each t from init_t until
 *   \f$ \mathcal{T} - 1\f$:
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 0 \f$ and \f$ v_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq 0 \f$.
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 1 \f$ and \f$ v_{t+1} = 1 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 0 \f$ and \f$ v_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq - \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 1 \f$ and \f$ v_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \Delta^+_t \f$.
 *
 *   Using the symmetry between ramp up and ramp down constraints, we can
 *   derive the ramp-down analogues of the ramp-up inequality as below:
 *   \f[
 *     p_t^{ac} - p_{t+1}^{ac} \leq ( - \Delta^-_t) w_{t+1}
 *       + (\underline{p}_t + \Delta^-_t)  u_t - \underline{p}_t u_{t+1}
 *            \quad t \in \{t_0, ..., \mathcal{T} - 1 \}             \quad (5)
 *   \f]
 *
 *   The sam analyzing the left hand side of the ramp-down constraint(5), in
 *   any integral feasible solution we can see that
 *   \f$ p_t^{ac} - p_{t+1}^{ac} \f$ can be bounded from above based on the
 *   values of \f$ u_{t+1}\f$, \f$ u_{t}\f$ and \f$ w_{t+1}\f$. Then for each
 *   (0, ..., f_time_horizon - 1) entries of this std::vector<FRowConstraint>
 *   there are two possible cases for t from 0 until init_t - 1:
 *
 *   - when \f$ u_{t} = 0\f$, and \f$ u_{t+1} = 0 \f$ then
 *     \f$ p_{t}^{ac} - p_{t+1}^{ac}  \leq 0 \f$.
 *
 *   - when \f$ u_{t} = 1\f$, and \f$ u_{t+1} = 1 \f$ then
 *     \f$ p_{t}^{ac} - p_{t+}^{ac} \leq \Delta^-_t \f$.
 *
 *   and four possible cases for each t from init_t until
 *   \f$ \mathcal{T} - 1\f$:
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 0 \f$ and \f$ w_{t+1} = 0 \f$ then
 *     \f$ p_t^{ac} - p_{t+1}^{ac}  \leq 0 \f$.
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 1 \f$ and \f$ w_{t+1} = 0 \f$ then
 *     \f$ p_t^{ac} - p_{t+1}^{ac} \leq \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 0 \f$ and \f$ w_{t+1} = 1 \f$ then
 *     \f$ p_t^{ac} - p_{t+1}^{ac} \leq - \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 1 \f$ and \f$ w_{t+1} = 0 \f$ then
 *     \f$ p_t^{ac} - p_{t+1}^{ac} \leq \Delta^-_t \f$.
 *
 * - Power output Constraints:
 *   Since commitment variable \f$ u_{t} \f$ is fixed to one or zero
 *   for "init_t" time steps(look above comments), because of power output
 *   constraint (look constraint (6)) when for the (0, ..., init_t - 1) time
 *   steps, commitment variable \f$ u_{t} \f$ is fixed to zero we must fix
 *   \f$ p_t^{ac} \f$, \f$ p_t^{pr} \f$, and \f$ p_t^{sc}\f$ to zero for the same
 *   time steps.
 *
 *   Maximum and minimum power output constraints according to active power,
 *   primary and secondary spinning reserves variables are presented in
 *   inequalities (6) and (7) respectively. Each of them is a
 *   std::vector<FRowConstraint> with exactly f_time_horizon entries
 *   (0, ..., (f_time_horizon) - 1) and ensures the maximum(or minimum) amount
 *   of energy that unit can produce(or use) when it is on(or off).
 *   \f[
 *
 *      p_t^{ac} + p_t^{pr} + p_t^{sc} \leq \bar{p}_t u_t          \quad (6)
 *
 *   \f]
 *
 *   \f[
 *
 *     \underline{p}_t u_t \leq p_t^{ac} - p_t^{pr} - p_t^{sc}   \quad (7)
 *
 *   \f]
 *
 *   The same as inequalities(6)-(7), the inequalities(8)-(9) ensure that
 *   maximum amount of primary and secondary spinning reserve in the problem
 *   respectively. Each of them is a std::vector<FRowConstraint> with exactly
 *   f_time_horizon entries (0, ..., (f_time_horizon) - 1) as below:
 *
 *   \f[
 *     p_t^{pr} \leq \rho^{pr}_t p_t^{ac}                       \quad (8)
 *   \f]
 *   \f[
 *     p_t^{sc} \leq \rho^{sc}_t p_t^{ac}                        \quad (9)
 *   \f]
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
 *            \quad t \in \{t_0 + 2, ..., \mathcal{T} - 1\} \quad  (10)
 *   \f]
 *
 *   and in the case  \f$ \tau_+ = 1 \f$:
 *
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) w_{t+1}
 *            \quad t \in \{t_0 + 2, ..., \mathcal{T} - 1\} \quad  (11)
 *   \f]
 *
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) v_t
 *             \quad t \in \{t_0 + 2, ...,  \mathcal{T} - 1\} \quad  (12)
 *   \f]
 *
 *   These inequalities give the active power output generation limits when
 *   unit is ON or OFF. More precisely, the unit generation limits taking into
 *   account its maximum \f$ \bar{p}_t \f$ and minimum \f$ \underline{p}_t \f$
 *   production, as well as its startup and shutdown capabilities(here both of
 *   them are assumed be equal with minimum production \f$ \underline{p}_t\f$)
 *   in each time step t. Be aware that (10) may be infeasible in the event
 *   that the unit is online for just one period. That is,
 *   \f$ v_t = w_{t+1} = 1 \f$ and the right side of the (10) can be negative.
 *   Consequently, (10) is only valid when \f$ \tau_+ \geq 2 \f$. Therefore,
 *   the correct formulation for units with \f$ \tau_+ = 1 \f$ is given by
 *   (11) and (12).
 *
 * Note that there may be other formulations (like the DP one), which will
 * possibly be implemented in the future. */

 void generate_abstract_constraints( Configuration *stcc ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the ThermalUnitBlock
/** Method that generates the objective of the ThermalUnitBlock.
 *
 * - Objective function: the objective function of the ThermalUnitBlock
 *   representing the total power production cost to be minimized has the
 *   form:
 *
 *   \f[
 *     \min ( \sum_{ t \in  [t_0 , \mathcal{T}]  } s_t v_t +
 *            \sum_{ t \in \mathcal{T}  } (a_t p_t^2 + b_t p_t + c_t u_t) )
 *   \f]
 *
 *   where \f$ \sum_{ t \in \mathcal{T} } s_t  v_t \f$ is the
 *   start-up cost of the unit, which we assume to be time-independent
 *
 *   Note: time-independent means here start-up cost is "independent from how
 *         long the unit has been off", and it is not meaning "always should
 *         be equal at each time instant"
 *
 *   and \f$ a_t \f$, \f$ b_t \f$, and \f$ c_t \f$ are, respectively, the
 *   quadratic, linear, and constant terms of the power cost function of the
 *   unit at time period \f$ t \in \mathcal{T} \f$. */

 void generate_objective( Configuration *objc ) override;

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
 *
 * @{ */

 /// returns the initial power value
 double get_initial_power() const { return f_initial_power; }

 /// returns the init up and down time value
 Index get_init_up_down_time() const { return f_InitUpDownTime; }

 /// returns the minimum allowed up time value
 Index get_min_up_time() const { return f_MinUpTime; }

 /// returns the minimum allowed down time value
 Index get_min_down_time() const { return f_MinDownTime; }

/*--------------------------------------------------------------------------*/

 /// returns the vector of nominal minimum active power output
 /** This method returns (a const reference to) the vector containing the
  * nominal minimum active power output of the unit for all time steps. When
  * the unit is available, get_min_power()[ t ] gives the minimum active power
  * output of the unit at time t, for each t in {0, ..., get_time_hotizon() -
  * 1}.  */
 const std::vector< double > & get_min_power() const {
  return v_MinPower;
 }

/*--------------------------------------------------------------------------*/

 /// returns the operational minimum active power output at the given time \t
 /** This method returns the operational minimum active power output of the
  * unit at the given time \t. See get_availability() for the definition of
  * operational minimum power.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The operational minimum active power output of the unit at the
  *         given time \t.
  */
 double get_operational_min_power( Index t ) const {
  assert( t < get_time_horizon() );
  return compute_operational_min_power( v_MinPower[ t ] ,
                                        get_availability( t ) );
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of nominal maximum active power output
 /** This method returns (a const reference to) the vector containing the
  * nominal maximum active power output of the unit for all time steps. When
  * the unit is fully available, get_min_power()[ t ] gives the maximum active
  * power output of the unit at time t, for each t in {0, ...,
  * get_time_hotizon() - 1}. See get_availability() to understand the
  * difference between nominal and operational maximum active power.  */
 const std::vector< double > & get_max_power() const {
  return v_MaxPower;
 }

/*--------------------------------------------------------------------------*/

 /// returns the operational maximum active power output at the given time \t
 /** This method returns the operational maximum active power output of the
  * unit at the given time \t. See get_availability() for the definition of
  * operational maximum power.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The operational maximum active power output of the unit at the
  *         given time \t.
  */
 double get_operational_max_power( Index t ) const {
  assert( t < get_time_horizon() );
  return compute_operational_max_power( v_MaxPower[ t ] ,
                                        get_availability( t ) );
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

 const std::vector< double > & get_availability() const {
  return v_Availability;
 }

/*--------------------------------------------------------------------------*/

 /// returns the availability of the unit at the given time \t
 double get_availability( Index t ) const {
  if( v_Availability.empty() )
   return 1.0;
  assert( t < get_time_horizon() );
  return v_Availability[ t ];
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of primary rho
 /** The returned vector contains the primary rho at each time.
  * The size of the vector is always get_time_horizon().
  */
 const std::vector< double > & get_primary_rho() const {
  return v_PrimaryRho;
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of secondary rho
 /** The returned vector contains the secondary rho at each time.
  * The size of the vector is always get_time_horizon().
  */
 const std::vector< double > & get_secondary_rho() const {
  return v_SecondaryRho;
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of delta ramp-up
 /** The returned vector contains the delta ramp-up at each time.
  * The size of the vector is always get_time_horizon().
  */
 const std::vector< double > & get_delta_ramp_up() const {
  return v_DeltaRampUp;
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of delta ramp-down
 /** The returned vector contains the delta ramp-up at each time.
  * The size of the vector is always get_time_horizon().
  */
 const std::vector< double > & get_delta_ramp_down() const {
  return v_DeltaRampDown;
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
 *   of vector represents the amount of quadratic term at time t.  */

 const std::vector< double > & get_quad_term() const { return( v_QuadTerm ); }

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
 *   of vector represents the amount of linear term at time t.  */

 const std::vector< double > & get_linear_term() const {
  return( v_LinearTerm );
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
 *   of vector represents the amount of constant term at time t.  */

 const std::vector< double > & get_const_term() const {
  return( v_ConstTerm );
  }

/*--------------------------------------------------------------------------*/
/// returns the vector of startup cost
/** The returned vector contains to startup cost at time t.  There are three
 * possible cases:
 *
 * - if the vector is empty, then the startup cost of the unit is 0;
 *
 * - if the vector has only one element, then the startup cost of the unit
 *   for all time horizon;
 *
 * - otherwise, the vector must have size get_time_horizon() and each element
 *   of vector represents the start up cost value at time t.  */

 const std::vector< double > & get_start_up_cost() const {
  return( v_StartUpCost );
  }

/*--------------------------------------------------------------------------*/
/// returns the vector of fixed consumption
/** The returned value U = get_fixed_consumption() contains the contribution
 *  to fixed consumption (basically, the constants to be multiplied by the
 *  commitment variables returned by get_commitment()) of all the generators
 *  at all time instants. There are three possible cases:
 *
 * - if the vector is empty, then the fixed consumption is always 0;
 *
 * - if the vector only has one element, then the fixed consumption for the
 *   fixed consumption of the unit for all t
 *
 * - otherwise, the vector must have size get_time_horizon(), and each element
 *   of vector represents the fixed consumption at time t. */

 double * get_fixed_consumption( Index generator )
  override {
  return & ( v_fixed_consumption.front() );
  }

/*--------------------------------------------------------------------------*/
/// returns the vector of inertia commitment
/** The returned value U = get_inertia_commitment() contains the contribution
 *  to inertia (basically, the constants to be multiplied by the commitment
 *  variables returned by get_commitment()) of all the generators at all time
 *  instants. There are three possible cases:
 *
 * - if the vector is empty, then the inertia commitment is always 0;
 *
 * - if the vector only has one element, then the inertia commitment for the
 *   fixed consumption of the unit for all t
 *
 * - otherwise, the vector must have size get_time_horizon(), and each element
 *   of vector represents the inertia commitment at time t. */

 double * get_inertia_commitment( Index generator )
  override {
  return & ( v_inertia_commitment.front() );
  }

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
 * - shut_down variables;
 *
 * @{ */
 /// returns the vector of commitment variables
 ColVariable * get_commitment( Index generator  ) override {
  return &( v_commitment.front() );
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of active_power variables
 ColVariable * get_active_power( Index generator )
 override {
  return &( v_active_power.front() );
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of primary_spinning_reserve variables
 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  return &( v_primary_spinning_reserve.front() );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary_spinning_reserve variables
 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  return &( v_secondary_spinning_reserve.front() );
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of start_up variables
 const std::vector< ColVariable > & get_start_up() const {
  return v_start_up;
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of shut_down variables
 const std::vector< ColVariable > & get_shut_down() const {
  return v_shut_down;
  }

/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE ThermalUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the ThermalUnitBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * ThermalUnitBlock. See ThermalUnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*-------------- METHODS FOR INITIALIZING THE ThermalUnitBlock -------------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the ThermalUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "ThermalUnitBlock::load() not implemented yet") );
  }

/**@} ----------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 // update the availability of the unit
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
  * an exception is thrown.  */

 void set_availability( std::vector< double >::const_iterator values,
                        Subset && subset,
                        const bool ordered = false,
                        c_ModParam issuePMod = eNoBlck,
                        c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 // update the availability of the unit
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
  * an exception is thrown.  */

 void set_availability( std::vector< double >::const_iterator values,
                        Range rng = Range( 0, Inf< Index >() ),
                        c_ModParam issuePMod = eNoBlck,
                        c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_maximum_power( std::vector< double >::const_iterator values,
                         Subset && subset,
                         const bool ordered = false,
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck );

 void set_maximum_power( std::vector< double >::const_iterator values,
                         Range rng = Range( 0, Inf< Index >() ),
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck );

 void set_initial_power( std::vector< double >::const_iterator values,
                         Subset && subset,
                         const bool ordered = false,
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck );

 void set_initial_power( std::vector< double >::const_iterator values,
                         Range rng = Range( 0, Inf< Index >() ),
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck );

 void set_init_updown_time( std::vector< int >::const_iterator values,
                            Subset && subset,
                            const bool ordered = false,
                            c_ModParam issuePMod = eNoBlck,
                            c_ModParam issueAMod = eNoBlck );

 void set_init_updown_time( std::vector< int >::const_iterator values,
                            Range rng = Range( 0, Inf< Index >() ),
                            c_ModParam issuePMod = eNoBlck,
                            c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/

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

 /// the InitialPower value
 double f_initial_power{};

 /// the vector of fixed consumption of generator
 std::vector< double > v_fixed_consumption;

 /// the vector of inertia commitment of generator
 std::vector< double > v_inertia_commitment;

 /// the MinUpTime value
 Index f_MinUpTime{};

 /// the MinDownTime value
 Index f_MinDownTime{};

 /// the InitUpDownTime value
 int f_InitUpDownTime{};

 /// variable denoting the time-steps unit is subjected to initial conditions
 Index init_t{};

/*-----------------------------variables------------------------------------*/
 /// the start up binary variables
 std::vector< ColVariable > v_start_up;

 /// the shut down binary variables
 std::vector< ColVariable > v_shut_down;

 /// the commitment variables
 std::vector< ColVariable > v_commitment;

 /// the active power variables
 std::vector< ColVariable > v_active_power;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;
/*----------------------------constraints-----------------------------------*/

 /// the connection power out put constraints
 std::vector< FRowConstraint > Power_StartUp_ShutDown_Variables_Constraints;

 /// the Start Up power out put constraints
 std::vector< FRowConstraint > Power_StartUp_Variable_Constraints;

 /// the Shut Down power out put constraints
 std::vector< FRowConstraint > Power_ShutDown_Variable_Constraints;

 /// the connection min up and down time constraints
 std::vector< FRowConstraint > StartUp_ShutDown_Variables_Constraints;

 /// the TURN ON min up and down time constraints
 std::vector< FRowConstraint > StartUp_Constraints;

 /// the SHUT DOWN min up and down time constraints
 std::vector< FRowConstraint > ShutDown_Constraints;

 /// the RampUp time constraints
 std::vector< FRowConstraint > RampUp_Constraints;

 /// the RampDown time constraints
 std::vector< FRowConstraint > RampDown_Constraints;

 /// the PrimaryRho fraction constraints
 std::vector< FRowConstraint > PrimaryRho_Constraints;

 /// the SecondaryRho fraction constraints
 std::vector< FRowConstraint > SecondaryRho_Constraints;

 std::vector< FRowConstraint > MinPower_Constraints;

 std::vector< FRowConstraint > MaxPower_Constraints;

 /// the commitment bound constraints
 std::vector< ZOConstraint > Commitment_bound_Constraints;

 /// the commitment fixed to one BoxConstraints
 std::vector< BoxConstraint > Commitment_fixed_to_One_Constraints;

 /// the startup binary bound constraints
 std::vector< ZOConstraint > StartUp_Binary_bound_Constraints;

 /// the shout down binary bound constraints
 std::vector< ZOConstraint > ShoutDown_Binary_bound_Constraints;

 /// the objective function
 FRealObjective objective;

 static void static_initialization() {
  /*!!
   * Not all C++ compilers enjoy the template wizardry behing the three-args
   * version of register_method<> with the compact MS_*_*::args(), so we just
   * use the slightly less compact one with the explicit argument and be done
   * with it. !!*/
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_availability",
  //                                      &ThermalUnitBlock::set_availability,
  //                                      MS_dbl_sbst::args() );
  //
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_availability",
  //                                      &ThermalUnitBlock::set_availability,
  //                                      MS_dbl_rngd::args() );
  //
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_maximum_power",
  //                                      &ThermalUnitBlock::set_maximum_power,
  //                                      MS_dbl_sbst::args() );
  //
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_maximum_power",
  //                                      &ThermalUnitBlock::set_maximum_power,
  //                                      MS_dbl_rngd::args() );
  //
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_initial_power",
  //                                      &ThermalUnitBlock::set_initial_power,
  //                                      MS_dbl_sbst::args() );
  //
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_initial_power",
  //                                      &ThermalUnitBlock::set_initial_power,
  //                                      MS_dbl_rngd::args() );
  //
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_init_updown_time",
  //                                      &ThermalUnitBlock::set_init_updown_time,
  //                                      MS_int_sbst::args() );
  //
  // register_method< ThermalUnitBlock >( "ThermalUnitBlock::set_init_updown_time",
  //                                      &ThermalUnitBlock::set_init_updown_time,
  //                                      MS_int_rngd::args() );
  register_method< ThermalUnitBlock, MF_dbl_it, Subset &&, const bool >(
   "ThermalUnitBlock::set_availability",
   &ThermalUnitBlock::set_availability );

  register_method< ThermalUnitBlock, MF_dbl_it, Range >(
   "ThermalUnitBlock::set_availability",
   &ThermalUnitBlock::set_availability );

  register_method< ThermalUnitBlock, MF_dbl_it, Subset &&, const bool >(
   "ThermalUnitBlock::set_maximum_power",
   &ThermalUnitBlock::set_maximum_power );

  register_method< ThermalUnitBlock, MF_dbl_it, Range >(
   "ThermalUnitBlock::set_maximum_power",
   &ThermalUnitBlock::set_maximum_power );

  register_method< ThermalUnitBlock, MF_dbl_it, Subset &&, const bool >(
   "ThermalUnitBlock::set_initial_power",
   &ThermalUnitBlock::set_initial_power );

  register_method< ThermalUnitBlock, MF_dbl_it, Range >(
   "ThermalUnitBlock::set_initial_power",
   &ThermalUnitBlock::set_initial_power );

  register_method< ThermalUnitBlock, MF_int_it, Subset &&, const bool >(
   "ThermalUnitBlock::set_init_updown_time",
   &ThermalUnitBlock::set_init_updown_time );

  register_method< ThermalUnitBlock, MF_int_it, Range >(
   "ThermalUnitBlock::set_init_updown_time",
   &ThermalUnitBlock::set_init_updown_time );
 }

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

 /// Resizes a vector to time_horizon by using change_intervals
 template< typename T > void decompress_vector( std::vector< T > & v );

/*--------------------------------------------------------------------------*/

 /// updates the abstract representation dependent on the availability
 /** This method updates any part of the abstract representation that may
  * depend on the availability of the unit at the given time \p t.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @param issueAMod controls how abstract Modification are issued. */
 void update_availability_dependents( Index t , c_ModParam issueAMod );

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
  * @return true if and only if the given \p availability is consistent at
  *         time \p t. */
 bool availability_is_consistent( Index t , double availability ) const {
  assert( t < get_time_horizon() );
  const auto min_power = compute_operational_min_power
   ( v_MinPower[ t ] , availability );
  const auto max_power = compute_operational_max_power
   ( v_MaxPower[ t ] , availability );
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
  if( availability > 0.0 )
   return nominal_min_power;
  return 0.0;
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
  return nominal_max_power * availability;
 }

};  // end( class( ThermalUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// Derived class from Modification for modifications to a ThermalUnitBlock
class ThermalUnitBlockMod : public Modification {

 public:

 /// Public enum for the types of ThermalUnitBlockMod
 enum TUBB_mod_type {
  eSetMaxP = 0 ,   ///< Set max power values
  eSetInitP    ,   ///< Set initial power values
  eSetInitUD   ,   ///< Set initial up/down times
  eSetAv           ///< Set availability
 };

 /// Constructor, takes the ThermalUnitBlock and the type
 ThermalUnitBlockMod( ThermalUnitBlock * const fblock, const int type )
  : f_Block( fblock ), f_type( type ) {}

 ///< Destructor, does nothing
 virtual ~ThermalUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block() const override { return ( f_Block ); }

 /// Accessor to the type of modification
 int type() { return ( f_type ); }

 protected:

 /// prints the ThermalUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "ThermalUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case ( eSetMaxP ):
    output << "set max power values ";
    break;
   case ( eSetInitP ):
    output << "set initial power values ";
    break;
   default:
    output << "set initial up/down times ";
  }
 }

 ThermalUnitBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

 int f_type; ///< type of modification
}; // end( class( ThermalUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS ThermalUnitBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from ThermalUnitBlockMod for "ranged" modifications
class ThermalUnitBlockRngdMod : public ThermalUnitBlockMod {

 public:

 /// constructor: takes the ThermalUnitBlock, the type, and the range
 ThermalUnitBlockRngdMod( ThermalUnitBlock * const fblock,
                          const int type,
                          Block::Range rng )
  : ThermalUnitBlockMod( fblock, type ), f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng() { return( f_rng ); }

 protected:

 /// prints the ThermalUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  ThermalUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng; ///< the range
};  // end( class( ThermalUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ThermalUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ThermalUnitBlockMod for "subset" modifications
class ThermalUnitBlockSbstMod : public ThermalUnitBlockMod {

 public:

 /// constructor: takes the ThermalUnitBlock, the type, and the subset
 ThermalUnitBlockSbstMod( ThermalUnitBlock * const fblock,
                          const int type,
                          Block::Subset && nms )
  : ThermalUnitBlockMod( fblock, type ), f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~ThermalUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms() { return( f_nms ); }

 protected:

 /// prints the ThermalUnitBlockSbstMod
 void print( std::ostream &output ) const override {
  ThermalUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms; ///< the subset

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
