/*--------------------------------------------------------------------------*/
/*------------------------- File ThermalUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class ThermalUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define the thermal unit of EDF
 * unit commitment Problem.
 *
 * A ThermalUnitBlock class is designed in order to give mathematical
 * formulation to describe the operation of large conventional power
 * plants (such as nuclear, hard coal, gas turbine, gas, combined
 * cycle, oil, ...)  which directly connected to the transmission
 * grid. The technical and physical constraints are mainly divided in
 * four different categories as bellow:
 *
 * - Maximum and minimum power output constraints;
 * - Ramp-up/down rate constraints;
 * - Minimum up and down time constraints;
 * - Start-up cost constraints.
 *
 * Based on the above description the class has been constructed having the
 * following elements:
 *
 * - A virtual public method is used in order to initialize and read
 *   the data of any possible derived ThermalUnitBlock class.
 *
 * - A set of protected methods, one for the initialization of each different
 *   type of constraints for the mathematical formulations with the methods
 *   also passing the constraints to the vector of static constraints of the
 *   Block.
 *
 * - A number of different vectors of FRowConstraint Objects that are used
 *   to store all the information of the different sets of constraints of
 *   the thermal unit.
 *
 * - An object of DQuadFunction and one of FRealFunction that are used to
 *   store the information of the quadratic or linear objective function of
 *   the unit.
 *
 * - Several different double variables that are used in order to store all
 *   the different values needed to describe the above mentioned constraints
 *   and costs.
 *
 * \version 0.11
 *
 * \date 09 - 06 - 2019
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

#include "ColVariable.h"
#include "FRowConstraint.h"
#include "OneVarConstraint.h"
#include "FRealObjective.h"
#include "DQuadFunction.h"
#include "UnitBlock.h"
#include "Block.h"

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
 * for the EDF Unit Commitment Problem.
 *
 * Consider a thermal generating unit and a set of time horizon \f$ T \f$
 * which is usually divided in a set of discrete time steps \f$ t \in T \f$.
 * Three binary decision variables are introduced as:
 * - \f$ u_t  \in \{ 0 , 1 \} \f$ : the commitment state of thermal unit at
 *   time period \f$ t \in T \f$.
 * - \f$  v_t \in \{ 0 , 1 \}  \f$: the start up of thermal unit at time
 *   period \f$ t \in T \f$.
 * - \f$  w_t \in \{ 0 , 1 \} \f$: the shut down of thermal unit at time
 *   period \f$ t \in T \f$. The main thermal unit constraints are categorized
 *   as following:
 *
 * - Min Up/Down-time Constraints:
 *   In the unit commitment problem, a thermal unit may has a
 *   minimum up time \f$ \tau_+ \f$ and a minimum down time
 *   \f$ \tau_- \f$ value. They refer to the minimum allowed up and down
 *   time of a generating thermal unit. It means, if thermal unit is committed
 *   in time \f$t\f$, then it must remain ON for the next \f$ \tau_+ - 1\f$
 *   time periods (and the same when shut down). By considering 3-binary
 *   variables \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ as defined above,
 *   the following constraints will be introduced as minimum up and down time
 *   constraints:
 *
 * \f[
 *  \sum_{ s \in ( t - \tau_+ + 1 , t ) } v_s \leq
 *               u_t \quad t \in \{ \tau_+ + 1, ..., T \}  \quad  (1)
 * \f]
 * \f[
 *  \sum_{ s \in ( t - \tau_- + 1 , t ) } w_s \leq
 *               1 - u_t \quad t \in \{ \tau_- + 1, ..., T \} \quad (2)
 * \f]
 * \f[
 *  u_t - u_{t-1} = v_t - w_t
 *                             \quad t \in \{ 2, ..., T \}  \quad (3)
 * \f]
 *
 * When unit in time \f$ t \f$ is OFF (\f$u_t = 0\f$), it could not have been
 * turned on in the last \f$ \tau_+ \f$ periods (including period \f$t\f$)
 * because of the minimum up constraints. But this is exactly what the turn on
 * inequality (1) for time period \f$t\f$ says. On the other hand, when unit in
 * time \f$t\f$ is ON (\f$u_t = 1\f$), it could have been turned on at most
 * once in the last \f$\tau_+ + \tau_- \f$ periods (including \f$t\f$).
 * Similarly for turn off inequality (2), when unit in the time \f$t\f$ is OFF
 * (\f$u_t = 0\f$), it could have been turned off at most once in the last
 * \f$ \tau_+ + \tau_-\f$ periods (including \f$t\f$). On the other hand, when
 * unit in time \f$t\f$ is ON (\f$u_t = 1\f$), it could not have been turned
 * off in the last \f$ \tau_-\f$ periods (including period \f$t\f$).
 * Since \f$ u_t \f$, \f$ v_t \f$, and \f$ w_t \f$ are binary variables, we
 * can ensure (for all periods \f$t \in \{ 2, ..., T \} \f$) that
 * \f$v_t = 1\f$ if and only if \f$u_t = 1\f$ and \f$u_{t-1} = 0\f$. It also
 * obvious \f$w_t = 1\f$ if and only if \f$u_t = 0\f$ and \f$u_{t-1} = 1\f$.
 * These conditions are satisfied by equality (3).
 *
 *
 *
 * - Ramp Up/Down-time Constraints:
 *   Another set of constraints where each thermal unit may has are ramp
 *   constraints. Here the two-period ramp up inequality is defined separately
 *   and the following constraints are proposed and shown to be valid for
 *   \f$ t= \{ 1, ..,T-1\}\f$ where \f$ \Delta^+_t \f$ and \f$ \Delta^-_t \f$
 *   are the constants defining ramp-up and ramp-down threshold and
 *   \f$ \underline{p}_t  \f$ and \f$ \bar{p}_t \f$  are the defining
 *   minimum and maximum output respectively:
 *
 *
 * \f[
 *   p_{t+1}^{ac} - p_t ^{ac} \leq
 *  ( - \Delta^+_t)  v_{t+1}
 *  + (\underline{p}_t- \Delta^{+}_t)  u_{t+1} -
 *  \underline{p}_t  u_t \quad t \in \{1, ..., T-1 \} \quad    (4)
 * \f]
 *
 *   Using the symmetry between ramping up and ramping down constraints, we
 *   can derive the ramp-down analogues of the ramp-up inequality as below:
 * \f[
 *   p_t^{ac} - p_{t+1} ^{ac} \leq
 *  ( - \Delta^{-}_t)  w_{t+1}
 *  + (\underline{p}_t - \Delta^{-}_t)  u_t -
 *  \underline{p}_t  u_{t+1} \quad t \in \{1, ..., T-1 \}   \quad   (5)
 * \f]
 *
 *
 * - Power output Constraints:
 *   There are several types of power out put inequalities which are
 *   considered below. More specially in case  \f$ 2 \leq \tau_+ \f$, the
 *   following constraint is introduced, which is valid for
 *   \f$ t \in \{2, ..., T-1\}  \f$:
 * \f[
 *   p_t^{ac}  \leq \bar{p}_t  u_t  - ( \bar{p}_t -
 *   \underline{p}^0_t  ) v_t - ( \bar{p}_t -
 *   \bar{p}^0_t )  w_{t+1}  \quad t \in \{2, ..., T-1\}\quad  (6)
 * \f]
 *
 *  and in the case  \f$ \tau_+ = 1 \f$:
 *
 *  \f[
 *   p_t^{ac} \leq \bar{p}_t  u_t  - ( \bar{p}_t -
 *   \bar{p}^0_t ) w_{t+1} - max( \bar{p}^0_t -
 *   \underline{p}^0_t , 0 ) v_t   \quad t \in \{2, ..., T-1\} \quad  (7)
 * \f]
 * \f[
 *   p_t^{ac} \leq \bar{p}_t  u_t  - ( \bar{p}_t -
 *   \underline{p}^0_t ) v_t - max( \underline{p}^0_t-
 *   \bar{p}^0_t  , 0 ) w_{t+1}   \quad t \in \{2, ..., T-1\} \quad  (8)
 * \f]
 *
 *   and for the \f$ t \in \{1, ..., T\}  \f$ following inequalities ensure
 *   the relation between active output and primary and secondary spinning
 *   reserves:
 * \f[
 *   p_t^{ac} + p_t^{pr} + p_t^{sc}\leq \bar{p}_t  u_t  \quad (9)
 *
 * \f]
 * \f[
 *
 *   \underline{p}_t  u_t \leq  p_t^{ac} - p_t^{pr} - p_t^{sc} \quad (10)
 * \f]
 *   considering an additional variable \f$ p_t \in R^{|T|}\f$ representing
 *   the power injected in to the grid by the power plant which differ from
 *   \f$ p_t^{ac} \f$, because the power plant in consuming a given power
 *   when it is off. we have the following relation for each
 *   \f$ t \in \{1, ..., T\}  \f$:
 *
 * \f[
 * //TODO this constraint goes to UCBlock
 *   {p}_t = (1 - u_t )P_t{au} + p_t^{ac} \quad (11)
 * \f]
 *   where \f$ P_t{au} \f$ denotes the fixed consumption of the power plant
 *   when it is off.
 *
 * - Active power relation with primary and secondary spinning reserves:
 *   for each
 *   \f$ t \in \{1, ..., T\}  \f$:
 * \f[
 *
 *   p_t^{pr} \leq \rho_t^{pr} p_t^{ac} \quad (12)
 * \f]
 *   and
 *
 * \f[
 *
 *   p_t^{sc} \leq \rho_t^{sc} p_t^{ac} \quad (13)
 * \f]
 *
 *
 * - Time independent Start-Up Costs:
 *   Here we just consider time independent star up cost which is simply equal
 *   to
 *
 * \f[
 *            s \sum_{ t \in T } v_t
 * \f]
 *   where \f$ s \f$ denotes the start up cost value of the unit when the
 *   start up thermal variable \f$ v_t = 1 \f$ at time \f$ t \f$.
 *
 * - objective function:
 * Given the constants and variables defined above, the objective function of
 * the unit commitment representing the total power production cost to be
 * minimized has the form:
 *
 * \f[
 *     min \quad ( s \sum_{ t \in T } v_t +
 *                 \sum_{ t \in T } (a_t p_t^2 + b_t p_t + c_t u_t) )
 * \f]
 *
 * where \f$ s \sum_{ t \in T } v_t \f$ is the total start-up cost of
 * the unit and \f$ a_t \f$, \f$ b_t \f$, and \f$ c_t \f$ are,
 * respectively, the quadratic, linear, and constant terms of the
 * power cost function of the unit at time period \f$t\f$.
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

/** Constructor of ThermalUnitBlock, taking possibly a pointer of its
 * father Block. */

 ThermalUnitBlock( Block * flbock = nullptr ): UnitBlock( flbock ) { }

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

 virtual void load( std::istream &input ) override { };

/*--------------------------------------------------------------------------*/
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the ThermalUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 * - the dimension "TimeHorizon" containing the time horizon;
 *
 * - the variable "MinPower", of type double and indexed over the dimension
 *   "NumberIntervals"; each entry of the variable is assumed to contain the
 *   minimum power output value of the unit for the corresponding time steps;
 *   it must be that MinPower[ i ] >= 0 for all i;
 *
 * - the variable "MaxPower", of type double and indexed over the dimension
 *   "NumberIntervals"; each entry of the variable is assumed to contain the
 *   maximum power output value of the unit for the corresponding time steps;
 *   it must be that MinPower[ i ] <= MaxPower[ i ] for all i;
 *
 *
 * - the variable "DeltaRampUp", of type double and indexed over the dimension
 *   "NumberIntervals"; each entry of the variable is assumed to contain the
 *   increases of power production value of the unit for the corresponding time
 *   steps in same interval whereas it may change(or not) in the other intervals
 *   (if exist any); this variable is optional, if it is not provided then it
 *   is assumed that DeltaRampUp == MaxPower, i.e., the unit can ramp up by
 *   an arbitrary amount, i.e., there are no ramp-up constraints;
 *
 * - the variable "DeltaRampDown", of type double and indexed over the
 *   dimension "NumberIntervals"; each entry of the variable is assumed to
 *   contain the decreases of power production value of the unit for the
 *   corresponding time steps in same interval whereas it may change(or not)
 *   in the other intervals  (if exist any); this variable is optional, if it
 *   is not provided then it is assumed that DeltaRampDown == MaxPower, i.e.,
 *   the unit can ramp down by an arbitrary amount, i.e., there are no
 *   ramp-down constraints;
 *
 * - the variable "PrimaryRho", of type double and indexed over the dimension
 *   "NumberIntervals"; each entry of the variable is assumed to contain the
 *   maximum fraction factor between primary power and the active power of
 *   the unit for the corresponding time steps in same interval whereas it may
 *   change(or not) in the other intervals (if exist any);
 *
 * - the variable "SecondaryRho", of type double and indexed over the
 *   dimension "NumberIntervals"; each entry of the variable is assumed to
 *   contain the maximum fraction factor between secondary power and the
 *   active power of the unit for the corresponding time steps in same
 *   interval whereas it may change(or not) in the other intervals (if exist
 *   any);
 *
 * - the variable "QuadTerm", of type double and indexed over the dimension
 *   "NumberIntervals"; each entry of the variable is assumed to contain the
 *   quadratic term of power cost function of the unit for the corresponding
 *   time steps in same interval whereas it may change(or not) in the other
 *   intervals (if exist any);
 *
 * - the scalar variable "StartUpCost", of type double and not indexed over
 *   any dimension and indicates the start up cost in this unit;
 *
 * - the variable "LinearTerm", of type double and indexed over the
 *   dimension "NumberIntervals"; each entry of the variable is assumed to
 *   contain the linear term of power cost function of the unit for the
 *   corresponding time steps in same interval whereas it may change(or not)
 *   in the other intervals (if exist any);
 *
 * - the variable "ConstTerm", of type double and indexed over the
 *   dimension "NumberIntervals"; each entry of the variable is assumed to
 *   contain the constant term of power cost function of the unit for the
 *   corresponding time steps if the unit is on;
 *
 * - the scalar variable "InitUpDownTime", of type UInt64 and not indexed over
 *   any dimension and indicates the initial time to generating the unit;
 *   if InitUpDownTime > 0, this means that the unit has been on for
 *   InitUpDownTime time stamps prior to time stamp 0 (the beginning of the
 *   horizon); if, instead, InitUpDownTime <= 0, this means that the unit
 *   has been off for - InitUpDownTime time stamps prior to time stamp 0;
 *   note that InitUpDownTime == 0 means that the unit has been just shut
 *   down at the end of time instant -1, i.e., the beginning of time
 *   instant 0;
 *
 * - the scalar variable "InitialPower", of type double and not
 *   indexed over any dimension; if InitUpDownTime > 0, it means that
 *   the unit was on at time instant -1 (prior to the beginning of the
 *   horizon), and then InitialPower indicates the amount of the power
 *   that the unit was producing at time instant -1; if InitUpDownTime
 *   <= 0 then this variable need not be defined since it is not
 *   loaded, if the variable is provided then it must be that MaxPower
 *   >= its value >= MinPower;
 *
 * - the scalar variable "MinUpTime", of type UInt64 and not indexed over
 *   any dimension and indicates the minimum allowed down time in this unit;
 *   this variable is optional, if it is not provided it is taken to be
 *   MinUpTime == 0, which mean that the unit can shut down in the very
 *   same time stamp in which it starts up;
 *
 * - the scalar variable "MinDownTime", of type UInt64 and not indexed over
 *   any dimension and indicates the minimum allowed up time in this unit;
 *   this variable is optional, if it is not provided it is taken to be
 *   MinDownTime == 0, which mean that the unit can start up in the very
 *   same time stamp in which it starts up;
 *
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
   override;

/// generate the abstract variables of the ThermalUnit
/** Method that generates the abstract variables of the ThermalUnitBlock.
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
 *  Note2: for the entry a = 0, ..., init_t - 1 when commitment variables fix
 *  to 0 then active power variables fix to 0, but when commitment variables
 *  fix to 1 the bound constraints, a std::vector<LB0Constraint> with exactly
 *  init_t entries, the entry a = 0, ..., init_t - 1 being the bound
 *  constraints of the ColVariable corresponding to the active power
 *  production of the unit.
 *
 * */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the ThermalUnit
/** Method that generates the static constraint of the ThermalUnitBlock.
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
 *
 * - Power output Constraints, a std::vector<FRowConstraint> with
 *   exactly ((f_time_horizon) - (init_t)) entries, the entry a = init_t, ...,
 *   (f_time_horizon) - 1 being the power output constraints at time t;
 *
 * - Time dependent Start-Up Cost Constraints
 * //TODO
 *
 * - Other Constraints,
 * //TODO
 */

 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
   override ;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the ThermalUnitBlock
/** Method that generates the objective of the ThermalUnitBlock.
 *
 */
 virtual void generate_objective( Configuration *objc = nullptr )
   override;

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
 * details of the format of the created netCDF group.
 *
 * */

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
/// print the ThermalUnitBlock on an ostream with the given verbosity
/** Protected method to print information about the ThermalUnitBlock;
 * */

// virtual void print( std::ostream &output ) const override ;

/*--------------------------------------------------------------------------*/
/// loads the ThermalUnit instance from standard .dat file format
/** Protected method for loading a ThermalUnitBlock out of a std::istream
 *
 *      //TODO
 */

//virtual void load( std::istream &input ) override final;


// void load(std::istream& inStream);

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

 /// the StartUpCost value
 double f_StartUpCost;

 /// the InitialPower value
 double f_initial_power;

 /* TODO Add the following to the comments of deserialize:
  *
  * - f_initial_min_power
  * - f_initial_delta_ramp_up
  * - f_initial_delta_ramp_down
  */

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
