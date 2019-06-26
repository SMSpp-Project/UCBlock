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
 * \date 25 - 06 - 2019
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

/// implementation of the Block concept for the thermal unit problem
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
 *
 * The operations of the thermal generating unit are described on a discrete
 * time horizon as dictated by the UnitBlock interface; in this description
 * we indicate it with \f$ \mathcal{T} \f$. For simplicity of notation it is
 * assumed that time steps are homogeneous with size \f$ \delta t\f$ in hours.
 * The first time instant is called "init_t" and depending on the presented
 * parameters InitUpDownTime(\f$ \tau_0 \f$), MinUpTime(\f$ \tau_+ \f$) and
 * MinDownTime(\f$ \tau_- \f$) is defined as below:
 *
 * - If \f$ \tau_0 > 0 \f$, this means that the unit has been on for
 *   \f$ \tau_0 \f$ time stamps prior to time stamp 0 (the beginning of the
 *   time horizon):
 *
 *   init_t = ( \f$ \tau_0 \f$ >= \f$ \tau_+ \f$ ? 0 :
 *              \f$ \tau_+ \f$ - \f$ \tau_0 \f$ );
 *
 * - If, instead, \f$ \tau_0 < 0 \f$, this means that the unit has
 *   been off for \f$ - \tau_0 \f$ time stamps prior to time stamp 0:
 *
 *   init_t = ( - \f$ \tau_0 \f$ >= \f$ \tau_- \f$ ? 0 :
 *              \f$ \tau_- \f$ + \f$ \tau_0 \f$ );
 *
 * - Note that \f$ \tau_0 == 0\f$ means that the unit has been just shut down
 *   at the end of time instant -1, i.e., the beginning of time instant 0.
 *
 * The value of \f$ \tau_0 \f$ impacts the minimum up and down time
 * constraints or ram-up and down constraints which are discussed in details
 * below.
 *
 * A possible MIP formulation of the problem uses three sets of binary
 * variables for each time instant \f$ t \in \mathcal{T} \f$:
 *
 * - \f$ u_t = 1\f$ if the unit is on at time period t;
 *   Note that the dimension of variable \f$ u_t \f$ is equal the dimension of
 *   time horizon \f$ \mathcal{T} \f$, and since other two variables below
 *   have the dimensions (\f$ \mathcal{T} \f$ - init_t) which is shorter than
 *   \f$ \mathcal{T} \f$, we need to fix this variable to 1 or 0 for the first
 *   init_t times(from 0 till init_t - 1), so:
 *
 *   - If \f$ \tau_0 \f$ < 0 and - \f$ \tau_0 \f$ < \f$ \tau_- \f$;
 *     then \f$ u_t\f$ must be fixed to 0 from 0 till init_t - 1 time steps
 *
 *   - If \f$ \tau_0 \f$ > 0 and  \f$ \tau_0 \f$ < \f$ \tau_+ \f$;
 *     then \f$ u_t\f$ must be fixed to 1 from 0 till init_t - 1 time steps
 *
 *
 * - \f$ v_t = 1 \f$ if the unit has started in time period t, i.e.,
 *   \f$ u_t = 1 \f$ but \f$ u_{t-1} = 0 \f$;
 *
 * - \f$  w_t = 1 \f$ if the unit shuts down in time period t, i.e.,
 *        \f$ u_{t-1}= 1 \f$ but \f$ u_t = 0 \f$.
 *
 * The main thermal unit constraints are categorized as following:
 *
 * - Min Up/Down-time Constraints: a thermal unit has a minimum up time
 *   \f$ \tau_+ \f$ and a minimum down time \f$ \tau_- \f$ value. It means
 *   that if thermal unit is started up in time \f$t\f$, then it must
 *   remain ON for the next \f$ \tau_+ - 1 \f$ time periods, it means if
 *   \f$ \tau_0 \geq 1 \f$, and that \f$ \tau_+ = 1 \f$ according to the
 *   definition of the first time step concept init_t there is no constraint.
 *
 *   One possible representation of the constraint in terms of the 3 binary
 *   variables \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ defined above is
 *   \f[
 *     u_t - u_{t-1} = v_t - w_t
 *          \quad t \in \{ 2, ...,\mathcal{T} \}                     \quad (1)
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
 *   When unit in time t is OFF (\f$ u_t = 0 \f$), it could not have been
 *   turned on in the last \f$ \tau_+ \f$ periods (including period t)
 *   because of the minimum up constraints. But this is exactly what the turn on
 *   inequality (2) for time period t\ says. On the other hand, when unit in
 *   time t is ON (\f$ u_t = 1 \f$), it could have been turned on at most
 *   once in the last \f$ \tau_+ + \tau_- \f$ periods (including t). Similarly
 *   for turn off inequality (3), when unit in the time t is OFF (\f$ u_t = 0
 *   \f$), it could have been turned off at most once in the last
 *   \f$ \tau_+ + \tau_-\f$ periods (including t). On the other hand, when
 *   unit in time t is ON (\f$ u_t = 1 \f$), it could not have been turned
 *   off in the last \f$ \tau_- \f$ periods (including period t). Since
 *   \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ are binary variables, we
 *   can ensure (for all periods \f$ t \in \{ 2, ..., \mathcal{T} \} \f$) that
 *   \f$ v_t = 1 \f$ if and only if \f$ u_t = 1 \f$ and \f$ u_{t-1} = 0 \f$.
 *   It also obvious that \f$ w_t = 1 \f$ if and only if \f$ u_t = 0 \f$ and
 *   \f$ u_{t-1} = 1 \f$. These conditions are satisfied by equality (1).
 *
 *   //TODO: discuss what happens when t <  \tau_+ + 1 or t <  \tau_- + 1
 *   does above explanation is enough?
 *
 * - Ramp Up/Down-time Constraints:
 *
 *   TODO: first discuss what the constraints should logically achieve
 *         (p_{t+1} \leq p_t + \Delta^+_t ...), then introduce

 *   Another set of constraints where each thermal unit may has are ramp
 *   constraints. Here the two-period ramp up inequality is defined separately
 *   and the following constraints are proposed and shown to be valid for
 *   \f$ t= \{ 1, ..,\mathcal{T} - 1\}\f$ where
 *   \f$ \Delta^+_t \f$ and \f$ \Delta^-_t \f$ are the constants defining
 *   ramp-up and ramp-down threshold and
 *   \f$ \underline{p}_t  \f$ and \f$ \bar{p}_t \f$  are the defining
 *   minimum and maximum output respectively. Let \f$ p_t^{ac} \f$ be the
 *   active power variable in time period t in all time horizon
 *   \f$ \mathcal{T} \f$.
 *
 *   Note that since commitment variable \f$ u_{t} \f$ is fixed to one or zero
 *   for "init_t" time steps(look above comments), and because of power output
 *   constraint (look constraint (9))we also fixed \f$ p_t^{ac} \f$ to zero
 *   for "init_t" time steps.
 *
 *   Analyzing the left hand side of the ramping constraint, in any integral
 *   feasible solution we can see that \f$ p_{t+1}^{ac} - p_t^{ac} \f$ can be
 *   bounded from above based on the values of \f$ u_{t+1}\f$, \f$ u_{t}\f$
 *   and \f$ v_{t+1}\f$. It may illustrate in four different ways:
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 0 \f$ and \f$ v_{t+1} = 0 \f$ then
 *     upper bound on LHS \f$ p_{t+1}^{ac} - p_t^{ac} = 0 \f$.
 *
 *   - when \f$ u_{t} = 0\f$, \f$ u_{t+1} = 1 \f$ and \f$ v_{t+1} = 1 \f$ then
 *     upper bound on LHS \f$ p_{t+1}^{ac} - p_t^{ac} =  \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 0 \f$ and \f$ v_{t+1} = 0 \f$ then
 *     upper bound on LHS \f$ p_{t+1}^{ac} - p_t^{ac} = - \underline{p}_t \f$.
 *
 *   - when \f$ u_{t} = 1\f$, \f$ u_{t+1} = 1 \f$ and \f$ v_{t+1} = 0 \f$ then
 *     upper bound on LHS \f$ p_{t+1}^{ac} - p_t^{ac} = \Delta^+_t \f$.
 *
 *   Considering the same logic for the ramp-down inequalities, nne possible
 *   implementation in terms of the three binary variables:
 *   \f[
 *     p_{t+1}^{ac} - p_t^{ac} \leq ( - \Delta^+_t)  v_{t+1}
 *        + (\underline{p}_t + \Delta^+_t) u_{t+1} - \underline{p}_t u_t
 *            \quad t \in \{ 1, ..., \mathcal{T} - 1 \} \quad (4)
 *   \f]
 *   Using the symmetry between ramp up and ramp down constraints, we can
 *   derive the ramp-down analogues of the ramp-up inequality as below:
 *   \f[
 *     p_t^{ac} - p_{t+1}^{ac} \leq ( - \Delta^-_t) w_{t+1}
 *       + (\underline{p}_t + \Delta^-_t)  u_t - \underline{p}_t u_{t+1}
 *            \quad t \in \{1, ..., \mathcal{T} - 1 \}             \quad (5)
 *   \f]
 *
 * - Power output Constraints:
 *
 *   TODO: again it is not very clear what is the logical condition that
 *         the power output constraint should satisfy. Please discuss
 *
 *   There are several types of power out put inequalities which are
 *   considered below. More specially in case  \f$ 2 \leq \tau_+ \f$, the
 *   following constraint is introduced, which is valid for
 *   \f$ t \in \{2, ..., \mathcal{T} - 1\}  \f$:
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t  u_t  - ( \bar{p}_t - \underline{p}_t ) v_t
 *                   - ( \bar{p}_t - \underline{p}_t ) w_{t+1}
 *            \quad t \in \{2, ..., \mathcal{T} - 1\} \quad  (6)
 *   \f]
 *   and in the case  \f$ \tau_+ = 1 \f$:
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) w_{t+1}
 *            \quad t \in \{2, ..., \mathcal{T} - 1\} \quad  (7)
 *   \f]
 *   \f[
 *     p_t^{ac} \leq \bar{p}_t u_t - ( \bar{p}_t - \underline{p}_t ) v_t
 *             \quad t \in \{2, ...,  \mathcal{T} - 1\} \quad  (8)
 *   \f]
 *   for the \f$ t \in \{1, ...,  \mathcal{T} \}  \f$ following
 *   inequalities ensure the relation between active output and primary and
 *   secondary spinning reserves:
 *   \f[
 *     p_t^{ac} + p_t^{pr} + p_t^{sc} \leq \bar{p}_t u_t          \quad (9)
 *   \f]
 *   \f[
 *     \underline{p}_t u_t \leq p_t^{ac} - p_t^{pr} - p_t^{sc}   \quad (10)
 *   \f]
 *
 *   TODO: explain what these do.
 *
 * - Objective function: the objective function of the ThermalUnitBlock
 *   representing the total power production cost to be minimized has the
 *   form:
 *   \f[
 *     \min ( \sum_{ t \in  \mathcal{T}  } s_t v_t +
 *            \sum_{ t \mathcal{T}  } (a_t p_t^2 + b_t p_t + c_t u_t) )
 *   \f]
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
 */

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
 * father Block.
 *
 * //TODO: if the constructor of UnitBlock takes the time horizon, why this
 *       one does not?
 */

 ThermalUnitBlock( Block * f_block = nullptr ): UnitBlock( f_block ) { }

/*--------------------------------------------------------------------------*/

 /// destructor of ThermalUnitBlock

 virtual ~ThermalUnitBlock() { };

/*@}------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// loads the ThermalUnitBlock instance from memory
/** Loads the ThermalUnitBlock instance from memory.
 * Like load( std::istream & ), if there is any Solver attached to this
 * ThermalUnitBlock then a NBModification (the "nuclear option") is issued.
 * */

 virtual void load( std::istream &input ) override {
  throw( std::logic_error( "ThermalUnitBlock::load() not implemented yet" ) );
  };

/*--------------------------------------------------------------------------*/
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the ThermalUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
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
 * - The scalar variable "InitialDeltaRampUp", of type double and not indexed
 *   over any dimension; it indicates the delta ramp up value at time instant
 *   zero(the initial condition);
 *
 *   //TODO: I don't agree with the name, and anyway the comment in unclear.
 *         What we need is InitialPower, i.e., the amount of power that the
 *         unit was producing at the time instant before 1. And you already
 *         have it below.
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
 * //TODO: I don't agree with this, ActiveP0 should work both for the ramp-up
 *       and for the ramp-down constraints
 *
 * - the scalar variable "InitialDeltaRampDown", of type double and not
 *   indexed over any dimension; it indicates the delta ramp down value at
 *   time instant zero(the initial condition);
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
 * //TODO: I don't agree, start-up cost can be time-dependent in the sense
 *       of being s_t, although not (for us) in the sense that it depends
 *       on how much the unit has been off beforw restarting
 *
 * - the variable "StartUpCost", of type double and indexed over the dimension
 *   "NumberIntervals". This is meant to represent the vector SC[ t ] which,
 *   for each time instant t, contains the start up cost value of the
 *   unit for the corresponding time steps;
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
 * - The scalar variable "InitialPower", of type double and not indexed over
 *   any dimension. If InitUpDownTime > 0, it means that the unit was on at
 *   time instant -1 (prior to the beginning of the horizon),
 *
 *   and then InitialPower indicates the amount of the power
 *   that the unit was producing at time instant -1; if InitUpDownTime
 *   <= 0 then this variable need not be defined since it is not
 *   loaded, if the variable is provided then it must be that MaxPower
 *   >= its value >= MinPower;
 *
 * //TODO: I don't understand, what't this for??
 *
 * - The scalar variable "InitialMinPower", of type double and not indexed
 *   over any dimension; it indicates the minimum power at time instant zero
 *   (the initial condition);
 *
 * - The scalar variable "MinUpTime", of type UInt64 and not indexed over
 *   any dimension, which indicates the minimum allowed down time in this
 *   unit. This variable is optional, if it is not provided it is taken to be
 *   MinUpTime == 0, which mean that the unit can shut down in the very
 *   same time stamp in which it starts up.
 *
 * - The scalar variable "MinDownTime", of type UInt64 and not indexed over
 *   any dimension, which indicates the minimum allowed up time in this unit.
 *   This variable is optional, if it is not provided it is taken to be
 *   MinDownTime == 0, which mean that the unit can start up in the very
 *   same time stamp in which it starts up. */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the ThermalUnit
/** Method that generates the abstract variables of the ThermalUnitBlock.
 *
 * //TODO: this comment is not very clear, please rewrite. In particular,
 *       you give access to these variables with start_up() and shut_down(),
 *       right? Why don't you mention these?
 *
 * //TODO: in UnitBlock we allow not to generate some of the variables with
 *       the stvv, why don't we here? Even if we don't, let's comment it.
 *
 * //TODO: one day we will do the DP formulation, and we will possibly have
 *       different groups of variables.
 *
 * These are as std::vector< ColVariable >  with exactly :
 *
 *  i). f_time_horizon entries, the entry a = 0, ...,
 *      (f_time_horizon) - 1 corresponding commitment and active power
 *       variables.
 *
 * ii). ((f_time_horizon) - (init_t)) entries, the entry a = init_t, ...,
 *      (f_time_horizon) - 1 corresponding start_up and
 *      shut_down variables.
 *
 *  Note1: commitment variable are fixed to 0 or 1 for the entry a = 0, ...,
 *      init_t - 1.
 *
 *
 * //TODO: this remark is about constraints, so it should go in
 *       generate_abstract_constraints()
 *
 *  Note2: for the entry a = 0, ..., init_t - 1 when commitment variables fix
 *  to 0 then active power variables fix to 0, but when commitment variables
 *  fix to 1 the bound constraints, a std::vector<LB0Constraint> with exactly
 *  init_t entries, the entry a = 0, ..., init_t - 1 being the bound
 *  constraints of the ColVariable corresponding to the active power
 *  production of the unit.
 *
 * */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
   override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the ThermalUnit
/** Method that generates the static constraint of the ThermalUnitBlock.
 *
 * //TODO: we could allow to only generate a subset of those via stcc.
 *       Maybe we don't want to.
 *
 * //TODO: one day we will do the DP formulation, and we will possibly have
 *       different groups of constraints.
 *
 * These are the:
 *
 * - Min Up/Down-time Constraints, a std::vector<FRowConstraint> with
 *   exactly ((f_time_horizon) - (init_t)) entries, the entry a = init_t, ...,
 *   (f_time_horizon) - 1 being the minimum up/down constraints at time t;
 *
 * - Ramp Up/Down-time Constraints, a std::vector<FRowConstraint> with
 *   exactly ((f_time_horizon) - (init_t)) entries, the entry a = init_t, ...,
 *   (f_time_horizon) - 1 being the ramp up/down time constraints at time t;
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

 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
   override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the ThermalUnitBlock
/** Method that generates the objective of the ThermalUnitBlock. */

 virtual void generate_objective( Configuration *objc = nullptr ) override;

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

 virtual void serialize( netCDF::NcGroup & group ) const override;

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
* @{ */

/*@} -----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE ThermalUnitBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ThermalUnitBlock
    @{ */

 /** returns the start up variable associated with time t such that
  * init_t <= t < time_horizon. */
 inline ColVariable & start_up( Index t ) {
   return v_start_up[ t - init_t ];
 }

 /** returns the shut down variable associated with time t such that
  * init_t <= t < time_horizon. */
 inline ColVariable & shut_down( Index t ) {
   return v_shut_down[ t - init_t ];
 }

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/

 /// the vector of change interval
 std::vector< int >  v_change_interval;

 /// the vector of MinPower
 std::vector< double >  v_MinPower;

 /// the vector of MaxPower
 std::vector< double >  v_MaxPower;

 /// the vector of PrimaryRho
 std::vector< double >  v_PrimaryRho;

 /// the vector of SecondaryRho
 std::vector< double >  v_SecondaryRho;

 /// the vector of InitPower
 std::vector< double >  v_StartUpLim;

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
 double f_initial_power;

 double f_initial_min_power;

 double f_initial_delta_ramp_up;

 double f_initial_delta_ramp_down;

 /// the MinUpTime value
 int f_MinUpTime;

 /// the MinDownTime value
 int f_MinDownTime;

 /// the InitUpDownTime value
 int f_InitUpDownTime;

 /// variable denoting the time-steps unit is subjected to initial conditions
 Index init_t;

/*-----------------------------variables------------------------------------*/

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
}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ThermalUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ThermalUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
