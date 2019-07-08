/*--------------------------------------------------------------------------*/
/*------------------------- File ThermalUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class ThermalUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" thermal unit
 * of a Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 03 - 07 - 2019
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
#include "MultiUnitBlock.h"

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

/// Implementation of the Block concept for the thermal unit problem
/** The ThermalUnitBlock class implements the Block concept [see Block.h]
 * for a "reasonably standard" thermal unit of a Unit Commitment Problem.
 * That is, the class is designed in order to give mathematical formulation to
 * describe the operation of large set of conventional power plants (such as
 * nuclear, hard coal, gas turbine, gas, combined cycle, oil, ...)  which
 * are directly connected to the transmission grid. The technical and physical
 * constraints are mainly divided in four different categories:
 *
 * - minimum up and down time constraints;
 * - ramp-up/down rate constraints;
 * - maximum and minimum power output constraints;
 * - active power relation with primary and secondary spinning reserves.
 * */

class ThermalUnitBlock : public MultiUnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * ThermalUnitBlock defines the following main public types:
 *
 * - Index, the type of parameters indices;
 *
 * @{ */

/*--------------------------------------------------------------------------*/

 typedef std::size_t Index;  ///< index of parameters

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// Constructor, takes the father and the time horizon
/** Constructor of ThermalUnitBlock, taking possibly a pointer of its
 * father Block.
 *
 * //TODO: if the constructor of MultiUnitBlock takes the time horizon, why this
 *       one does not?
 */

 explicit ThermalUnitBlock( Block * f_block = nullptr ): MultiUnitBlock( f_block ) { }

/*--------------------------------------------------------------------------*/

 /// Destructor of ThermalUnitBlock

 ~ThermalUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the ThermalUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base MultiUnitBlock, as
 * described in the comments to MultiUnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The variable "MinPower", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector MnP[ t ] which,
 *   for each time instant t, contains the minimum power output value of the
 *   unit for the corresponding time steps; it must be that MnP[ t ] >= 0 for
 *   all t. MinPower[ i ] is the fixed value of MnP[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1
 *   or "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxPower", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector MxP[ t ] which,
 *   for each time instant t, contains the maximum power output value of the
 *   unit for the corresponding time steps; it must be that MxP[ t ] >=
 *   MnP[ t ] >= 0 for all t. MaxPower[ i ] is the fixed value of MxP[ t ]
 *   for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0. If "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon"
 *   then the mapping clearly does not require "ChangeIntervals", which in
 *   fact is not loaded.
 *
 * - The variable "DeltaRampUp", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector DP[ t ] which,
 *   for each time instant t, contains the maximum possible increase of power
 *   production w.r.t. the power that had been produced in time instant t - 1,
 *   if any. This variable is optional; if it is not provided then it is
 *   assumed that DP[ t ] == MxP[ t ], i.e., the unit can ramp up by an
 *   arbitrary amount, i.e., there are no ramp-up constraints.
 *   DeltaRampUp[ i ] is the fixed value of DP[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "DeltaRampDown", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector DM[ t ] which,
 *   for each time instant t, contains the maximum possible decrease of power
 *   production w.r.t. the power that had been produced in time instant t - 1,
 *   if any. This variable is optional; if it is not provided then it is
 *   assumed that DM[ t ] == MxP[ t ], i.e., the unit can ramp down by an
 *   arbitrary amount, i.e., there are no ramp-down constraints.
 *   DeltaRampDown[ i ] is the fixed value of DM[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "PrimaryRho", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector PR[ t ] which,
 *   for each time instant t, contains the maximum possible fraction of
 *   active power that can be used as primary reserve.
 *
 *   //TODO: is this variable optional? Is is possible that a unit may not be
 *         capable of producing any primary reserve, which correspond to
 *         PR[ t ] == 0 for all t?
 *
 *   PrimaryRho[ i ] is the fixed value of PR[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "SecondaryRho", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector SR[ t ] which,
 *   for each time instant t, contains the maximum possible fraction of
 *   active power that can be used as secondary reserve.
 *
 *   //TODO: is this variable optional? Is is possible that a unit may not be
 *         capable of producing any primary reserve, which correspond to
 *         SR[ t ] == 0 for all t? Is there a logic relationship with
 *         PrimaryRho, like PR[ i ] == 0 ==> SR[ i ] == 0??
 *
 *   SecondaryRho[ i ] is the fixed value of SR[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "QuadTerm", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector A[ t ] which,
 *   for each time instant t, contains the quadratic term of power cost
 *   function of the unit for the corresponding time steps. This variable is
 *   optional; if it is not provided then it is assumed that A[ t ] == 0,
 *   i.e., the cost of the unit is linear in the produced power. QuadTerm[ i ]
 *   is the fixed value of A[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - the variable "StartUpCost", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector SC[ t ] which,
 *   for each time instant t, contains the start up cost value of the
 *   unit for the corresponding time steps; //todo look better
 *
 * - The variable "LinearTerm", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector B[ t ] which,
 *   for each time instant t, contains the linear term of power cost
 *   function of the unit for the corresponding time steps. This variable is
 *   optional; if it is not provided then it is assumed that B[ t ] == 0,
 *   i.e., the cost of the unit has no linear dependence on the produced power
 *   (say, only the quadratic one). LinearTerm[ i ] is the fixed value of
 *   B[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that
 *   ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "ConstTerm", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector C[ t ] which,
 *   for each time instant t, contains the constant term of power cost
 *   function of the unit for the corresponding time steps. This variable is
 *   optional; if it is not provided then it is assumed that C[ t ] == 0,
 *   i.e., the cost of the unit has no fixed term, only those depening
 *   (linearly or quadratically) on the produced power. ConstTerm[ i ] is
 *   the fixed value of C[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The scalar variable "InitialPower", of type double and not indexed over
 *   any dimension. If InitUpDownTime > 0, it means that the unit was on at
 *   time instant -1 (prior to the beginning of the horizon). It indicates the
 *   amount of the power that the unit was producing at time instant -1; if
 *   InitUpDownTime <= 0 then this variable need not be defined since it is
 *   not loaded, if the variable is provided then it must be that MaxPower
 *   >= its value >= MinPower;
 *
 * - The scalar variable "InitUpDownTime", of type UInt64 and not indexed over
 *   any dimension and indicates the initial time to generating the unit.
 *   If InitUpDownTime > 0, this means that the unit has been on for
 *   InitUpDownTime time stamps prior to time stamp 0 (the beginning of the
 *   horizon). If, instead, InitUpDownTime <= 0, this means that the unit
 *   has been off for - InitUpDownTime time stamps prior to time stamp 0;
 *   note that InitUpDownTime == 0 means that the unit has been just shut
 *   down at the end of time instant -1, i.e., the beginning of time
 *   instant 0;
 *
 * - The positive scalar variable "MinUpTime", of type UInt64 and not indexed
 *   over any dimension, which indicates the minimum allowed down time in this
 *   unit. This variable is optional, if it is not provided it is taken to be
 *   MinUpTime == 0, which mean that the unit can shut down in the very
 *   same time stamp in which it starts up.
 *
 * - The positive scalar variable "MinDownTime", of type UInt64 and not
 *   indexed over any dimension, which indicates the minimum allowed up time
 *   in this unit.This variable is optional, if it is not provided it is taken
 *   to be MinDownTime == 0, which mean that the unit can start up in the very
 *   same time stamp in which it starts up. */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// Generate the abstract variables of the ThermalUnitBlock
/** The ThermalUnitBlock class use get_variable() method to access to each
 *  "group" of variable that may create in MultiUnitBlock class which are:
 *
 *  - the binary commitment variables which takes the value of 1 if unit is ON
 *    at time instant t and 0 otherwise;
 *
 *  - the primary spinning reserve variables;
 *
 *  - the secondary spinning reserve variables;
 *
 *  - the active power variables;
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
 *  All of these variables are optional,and it is also possible to restrict
 *  which of the subsets are generated with the parameter stvv. If stvv is not
 *  nullptr and it is a SimpleConfiguration<int>, or if
 *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 *  SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 *  of the optional variables should be created. If the Configuration is not
 *  available, the default value is taken to be 0.
 * */

 void generate_abstract_variables( Configuration *stvv ) override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Variable OF THE ThermalUnitBlock -------*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the ThermalUnitBlock
 *
 * These methods allow to read the each group of Variable that any
 * ThermalUnitBlock in principle has (although some may not):
 *
 * - start_up variables;
 *
 * - shut_down variables;
 *
 * @{ */

 /// Returns the vector of start_up variables
 const std::vector< ColVariable > & get_start_up() const {
  return v_start_up;
 }
/*--------------------------------------------------------------------------*/
 /// Returns the vector of shut_down variables
 const std::vector< ColVariable > & get_shut_down() const {
  return v_shut_down;
 }

/*--------------------------------------------------------------------------*/
/// Generate the static constraint of the ThermalUnitBlock
/** This method generates the abstract constraints of the ThermalUnitBlock.
 *
 * The operations of the thermal generating unit are described on a discrete
 * time horizon as dictated by the MultiUnitBlock interface. In this description
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
 *   production of the unit.
 *
 * Then the main thermal unit constraints with three 3 binary variables
 * \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$are are presented as following:
 *
 * - Min Up/Down-time Constraints: a thermal unit may have minimum up and down
 *   time constraints and one possible representation of the constraint in
 *   terms of the
 *   defined above is:
 *
 *   \f[
 *     u_t - u_{t-1} = v_t - w_t
 *          \quad t \in \{ 1, ...,\mathcal{T} \}                     \quad (1)
 *   \f]
 *   \f[
 *    \sum_{ s \in ( t - \tau_+ + 1 , t ) } v_s \leq
 *           u_t \quad t \in \{ \tau_+ + 1, ..., \mathcal{T}\}
 *                                                                   \quad (2)
 *   \f]
 *   \f[
 *    \sum_{ s \in ( t - \tau_- + 1 , t ) } w_s \leq
 *         1 - u_t \quad t \in \{ \tau_- + 1, ...,\mathcal{T} \}
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
 *   can ensure (for all periods \f$ t \in \{ 1, ..., \mathcal{T} \} \f$) that
 *   \f$ v_t = 1 \f$ if and only if \f$ u_t = 1 \f$ and \f$ u_{t-1} = 0 \f$.
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
 *            \quad t \in \{ 1, ..., \mathcal{T} - 1 \} \quad (4)
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
 *   (0, ..., f_time_horizon - 1) entries of this std::vector<FRowConstraint>
 *   there are four possible cases:
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 0 \f$ and \f$ v_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac}  \leq 0 \f$.
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
 *            \quad t \in \{1, ..., \mathcal{T} - 1 \}             \quad (5)
 *   \f]
 *
 *   The sam analyzing the left hand side of the ramp-down constraint(5), in
 *   any integral feasible solution we can see that
 *   \f$ p_t^{ac} - p_{t+1}^{ac} \f$ can be bounded from above based on the
 *   values of \f$ u_{t+1}\f$, \f$ u_{t}\f$ and \f$ w_{t+1}\f$. Then for each
 *   (0, ..., f_time_horizon - 1) entries of this std::vector<FRowConstraint>
 *   there are four possible cases:
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 0 \f$ and \f$ w_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac}  \leq 0 \f$.
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 1 \f$ and \f$ w_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 0 \f$ and \f$ w_{t+1} = 1 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq - \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 1 \f$ and \f$ w_{t+1} = 0 \f$ then
 *     \f$ p_{t+1}^{ac} - p_t^{ac} \leq \Delta^-_t \f$.
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
 *     p_t^{pr} \leq \rho{pr}_t p_t^{ac}                       \quad (8)
 *   \f]
 *   \f[
 *     p_t^{sc} \leq \rho{sc}_t p_t^{ac}                        \quad (9)
 *   \f]
 *   There are two more power out put tighter formulations which make the
 *   maximum power output being a function of three binary variables
 *   \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ as below. More specifically in
 *   the case  \f$ \tau_+ \geq 2 \f$, the
 *   following constraint is introduced, which is valid for
 *   \f$ t \in \{2, ..., \mathcal{T} - 1\}  \f$:
 *
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t  u_t  - ( \bar{p}_t - \underline{p}_t ) v_t
 *                   - ( \bar{p}_t - \underline{p}_t ) w_{t+1}
 *            \quad t \in \{2, ..., \mathcal{T} - 1\} \quad  (10)
 *   \f]
 *
 *   and in the case  \f$ \tau_+ = 1 \f$:
 *
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) w_{t+1}
 *            \quad t \in \{2, ..., \mathcal{T} - 1\} \quad  (11)
 *   \f]
 *
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) v_t
 *             \quad t \in \{2, ...,  \mathcal{T} - 1\} \quad  (12)
 *   \f]
 *
 *   for the \f$ t \in \{0, ...,  \mathcal{T} - 1 \}  \f$ following
 *   inequalities ensure the relation between active output and primary and
 *   //TODO: explain what these do.
 *
 * - Objective function: the objective function of the ThermalUnitBlock
 *   representing the total power production cost to be minimized has the
 *   form:
 *
 *   \f[
 *     \min ( \sum_{ t \in  \mathcal{T}  } s_t v_t +
 *            \sum_{ t \mathcal{T}  } (a_t p_t^2 + b_t p_t + c_t u_t) )
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
 *   unit at time period \f$ t \in \mathcal{T} \f$.
 *
 * //TODO: we could allow to only generate a subset of those via stcc.
 *       Maybe we don't want to.
 *
 * //TODO: one day we will do the DP formulation, and we will possibly have
 *       different groups of constraints.
 *
 * These are the:
 *
 *
 * - Power output Constraints, a std::vector<FRowConstraint> with
 *   exactly ((f_time_horizon) - (init_t)) entries, the entry a = init_t, ...,
 *   (f_time_horizon) - 1 being the power output constraints at time t;
 *
 * - Relation between Power output, Primary and Secondary Spinning reserves
 *   Constraints, a std::vector<FRowConstraint> with exactly f_time_horizon
 *   entries, the entry a = 0, ..., (f_time_horizon) - 1 being the relation
 *   between Power output Primary and Secondary Spinning reserves constraints
 *   at time t;
 *
 * - Relation between Power output and Primary Spinning reserves Constraints,
 *   a std::vector<FRowConstraint> with exactly f_time_horizon entries, the
 *   entry a = 0, ..., (f_time_horizon) - 1 being the relation between Power
 *   output and Primary Spinning reserves constraints at time t;
 *
 * - Relation between Power output and Secondary Spinning reserves Constraints
 *   a std::vector<FRowConstraint> with exactly f_time_horizon entries, the
 *   entry a = 0, ..., (f_time_horizon) - 1 being the relation between Power
 *   output and Secondary Spinning reserves constraints at time t;
 */

 void generate_abstract_constraints( Configuration *stcc ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the ThermalUnitBlock
/** Method that generates the objective of the ThermalUnitBlock. */

 void generate_objective( Configuration *objc ) override;

/*@} -----------------------------------------------------------------------*/
/*----------- Methods for reading the data of the ThermalUnitBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the ThermalUnitBlock
 *  @{  */

/*@} -----------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE ThermalUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the ThermalUnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
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
 };

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/

 /// the vector of change interval
 std::vector< Index >  v_change_interval;

 /// the vector of MinPower
 std::vector< double >  v_MinPower;

 /// the vector of MaxPower
 std::vector< double >  v_MaxPower;

 /// the vector of PrimaryRho
 std::vector< double >  v_PrimaryRho;

 /// the vector of SecondaryRho
 std::vector< double >  v_SecondaryRho;

 /// the vector of RampUp
 std::vector< double >  v_DeltaRampUp;

 /// the vector of RampDown
 std::vector< double >  v_DeltaRampDown;

 /// the vector of QuadTerm
 std::vector< double >  v_QuadTerm;

 /// the vector of LinearTerm
 std::vector< double >  v_LinearTerm;

 /// the vector of ConstTerm
 std::vector< double >  v_ConstTerm;

 /// the vector of StartUpCost
 std::vector< double > v_StartUpCost;

 /// the InitialPower value
 double f_initial_power{};

 double f_initial_min_power{}; //TODO remove

 double f_initial_delta_ramp_up{}; //TODO remove

 double f_initial_delta_ramp_down{}; //TODO remove

 /// the MinUpTime value
 Index f_MinUpTime;

 /// the MinDownTime value
 Index f_MinDownTime;

 /// the InitUpDownTime value
 Index f_InitUpDownTime;

 /// variable denoting the time-steps unit is subjected to initial conditions
 Index init_t;

/*-----------------------------variables------------------------------------*/
 /* Each of the following vectors of Variable may either have size
  * f_time_horizon - init_t, meaning that there is not any Variable for each
  * defining time step (init_t , ..., f_time_horizon-1), or be empty, in which
  * case the variables simply do not exist. */

 /// the start up binary variables
 std::vector< ColVariable > v_start_up;

 /// the shut down binary variables
 std::vector< ColVariable > v_shut_down;

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

 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
private:

  SMSpp_insert_in_factory_h;

}; // end( class( ThermalUnitBlock ) )

/*@}  end( class( ThermalUnitBlock ) ) -------------------------------------*/
/*--------------------------------------------------------------------------*/
} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ThermalUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ThermalUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
