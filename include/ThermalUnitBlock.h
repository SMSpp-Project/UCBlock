/*--------------------------------------------------------------------------*/
/*------------------------- File ThermalUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class ThermalUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define the thermal unit of EDF
 * unit commitment Problem.
 *
 * A ThermalUnitBlock class is designed in order to give mathematical
 * formulation to describe the operation of large conventional power plants
 * (such as Nuclear, Hard coal, Gas turbine, Gas, Combined cycle, Oil, ...)
 * which directly connected to the transmission grid. The technical and
 * physical constraints are mainly divided in  for different categories as
 * bellow:
 *
 * - Maximum and Minimum power output constraints
 * - Ramp Up/Down rate constraints
 * - Min Up/Down time constraints
 * - Start-up cost constraints
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
 *
 *
 * \version 0.11
 *
 * \date 18 - 04 - 2019
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
 *
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
 * \f[
 *     min \quad ( s \sum_{ t \in T } v_t + \sum_{ t \in T } (a_t p_t^2 + b_t p_t
 *                                          + c_t u_t))
 * \f]
 * where \f$ s \sum_{ t \in T } v_t \f$ is the start-up cost of the unit and
 * \f$ a_t \f$, \f$ b_t \f$, and  \f$ c_t \f$  are the quadratic, linear and
 * constant term of power cost function of the unit at time period \f$t\f$
 * respectively.
 * */



 class ThermalUnitBlock : public UnitBlock {


/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
    public:
/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * ThermalUnitBlock defines several main public types:
 *
 * - Index, the type of parameters indices;
 *
 *
 @{ */


/*--------------------------------------------------------------------------*/
typedef unsigned int Index;                 ///< index of parameters
typedef const Index c_Index;                ///< a read-only Index

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/** Constructor of ThermalUnitBlock, taking possibly a pointer of its father
  *  Block, alongside with setting the unit with ramp_constraints */

ThermalUnitBlock( Block * flbock = nullptr ): UnitBlock( flbock ) {


}

/*--------------------------------------------------------------------------*/

/// destructor of ThermalUnitBlock

virtual ~ThermalUnitBlock() { } ;



/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// loads the ThermalUnitBlock instance from memory
/** Loads the ThermalUnitBlock instance from memory.
 * Like load( std::istream & ), if there is any Solver attached to this
 * ThermalUnitBlock then a NBModification (the "nuclear option") is issued.
 * */

virtual void load( );


/*--------------------------------------------------------------------------*/
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the ThermalUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 *
 * - the dimension "TimeHorizon" containing the number of time steps in
 *   this unit;
 *
 * - let's suppose some variables (such as MinPower, MaxPower , ...) may
 *   change independently, and may have different values for some intervals
 *   along the "TimeHorizon".  Without loss of generality, let's collect the
 *   intersection of all the intervals between each two parameters separately.
 *   It means, at the end we may have some intervals such as:
 *   \f$ [0 , a] , [a+1 , b], ... , [j+1 , k] , [k+1, T] \f$ (where
 *   \f$ a, b, ... , k < T \f$ and are positive integer numbers in which
 *   \f$T\f$ is the length of the TimeHorizon). In each interval one parameter
 *   may change or not(if not, copy the corresponding value of its' previous
 *   interval).
 *
 * - the dimension "NumberValues" which is a subset of "TimeHorizon" and
 *   indicates the number of above intervals \f$([0 , a] , [a+1 , b], ... ,
 *   [j+1 , k] , [k+1, T])\f$ where the variables change. The dimension is
 *   optional, if it is not provided than it is taken to be 1.
 *
 *   Three scenarios may happen:
 *
 *    i). In the simplest case scenario, "NumberValues = 1" which means
 *        in all the period of "TimeHorizon"; all the variables are fixed.
 *
 *   ii). In the average case scenario, "1 < NumberValues < TimeHorizon" which
 *        means in some time steps of defined "TimeHorizon"; some of the
 *        variables are changing.
 *
 *  iii). In the worst case scenario, "NumberValues = TimeHorizon" which means
 *         in every time step of defined "TimeHorizon"; all the variables are
 *         changing.
 *
 * - the variable "ChangeIntervals", of type integer and indexed over the
 *   dimension "NumberValues"; the \f$t_{th}\f$ entry of the variable
 *   indicates the positive number of \f$ a, b, ..., k, T\f$ on the above
 *   example.
 *
 * - the variable "MinPower", of type double and indexed over the dimension
 *   "TimeHorizon"; each entry of the variable is assumed to contain the
 *   minimum power output of the unit at time t, and it must be equal to other
 *   entries in the same interval whereas it may change(or not) in the other
 *   intervals (if exist any);
 *
 * - the variable "MaxPower", of type double and indexed over the dimension
 *   "TimeHorizon"; each entry of the variable is assumed to contain the
 *   maximum power output of the unit at time t, and it must be equal to other
 *   entries in the same interval whereas it may change(or not) in the other
 *   intervals (if exist any);
 *
 * - the variable "FixedConsPower", of type UInt64 and not indexed over
 *   any dimension and indicates the fixed consumption of the power plant when
 *   it is off in this unit;
 *
 * - the variable "DeltaRampUp", of type double and indexed over the dimension
 *   "TimeHorizon"; each entry of the variable is assumed to contain the
 *   increases of power production of the unit at time t, and it must be equal
 *   to other entries in the same interval whereas it may change(or not) in the
 *   other intervals (if exist any);
 *
 * - the variable "DeltaRampDown", of type double and indexed over the
 *   dimension "TimeHorizon"; each entry of the variable is assumed to contain
 *   the decreases of power production of the unit at time t, and it must be
 *   equal to other entries in the same interval whereas it may change(or not)
 *   in the other intervals (if exist any);
 *
 * - the variable "PrimaryRho", of type double and indexed over the dimension
 *   "TimeHorizon"; each entry of the variable is assumed to contain the
 *   maximum fraction factor between primary power and the active power for
 *   the unit at time t, and it must be equal to other entries in the same
 *   interval whereas it may change(or not) in the other intervals
 *   (if exist any);
 *
 * - the variable "SecondaryRho", of type double and indexed over the
 *   dimension "TimeHorizon"; each entry of the variable is assumed to contain
 *   the maximum fraction factor between secondary power and the active power
 *   for the unit at time t, and it must be equal to other entries in the same
 *   interval whereas it may change(or not) in the other intervals
 *   (if exist any);
 *
 * - the variable "QuadTerm", of type double and indexed over the
 *   dimension "TimeHorizon"; each entry of the variable is assumed to contain
 *   the quadratic term of power cost function for the unit at time t, and it
 *   must be equal to other entries in the same interval whereas it may change
 *   (or not) in the other intervals (if exist any);
 *
 * - the scalar variable "StartUpCost", of type UInt64 and not indexed over
 *   any dimension and indicates the start up cost in this unit;
 *
 * - the variable "LinearTerm", of type double and indexed over the
 *   dimension "TimeHorizon"; each entry of the variable is assumed to contain
 *   the linear term of power cost function for the unit at time t, and it
 *   must be equal to other entries in the same interval whereas it may change
 *   (or not) in the other intervals (if exist any);
 *
 * - the variable "ConstTerm", of type double and indexed over the
 *   dimension "TimeHorizon"; each entry of the variable is assumed to contain
 *   the constant term of power cost function for the unit at time t, and it
 *   must be equal to other entries in the same interval whereas it may change
 *   (or not) in the other intervals (if exist any);
 *
 * - the scalar variable "PZero", of type UInt64 and not indexed over any
 *   dimension and indicates the initiate amount of the power in this unit;
 *
 * - the scalar variable "MinUpTime", of type UInt64 and not indexed over
 *   any dimension and indicates the minimum allowed down time in this unit;
 *
 * - the scalar variable "MinDownTime", of type UInt64 and not indexed over
 *   any dimension and indicates the minimum allowed up time in this unit;
 *
 * - the scalar variable "InitUpDownTime", of type UInt64 and not indexed over
 *   any dimension and indicates the initial time to generating the unit;
 *
 */

virtual void deserialize( netCDF::NcGroup & group ,
        Block *father = nullptr ) override;

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
 *      (f_time_horizon) - 1 corresponding start_up_thermal and
 *      shut_down_thermal variables.
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
 * - Start-Up Cost Constraints,
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
 * //TODO
 *
*/
virtual void generate_objective( Configuration *objc = nullptr )
        override ;
/*@} -----------------------------------------------------------------------*/
/*----------- Methods for reading the data of the ThermalUnitBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the ThermalUnitBlock
 *  @{  */

/// the number of values
Index f_number_values;

/// the vector of change interval
std::vector < double >  f_change_interval;

/// the vector of MinPower
std::vector < double >  f_MinPower;

/// the vector of MaxPower
std::vector < double >  f_MaxPower;

/// the vector of PrimaryRho
std::vector< double >  f_PrimaryRho;

/// the vector of SecondaryRho
std::vector< double >  f_SecondaryRho;

/// the vector of InitPower
std::vector < double >  f_StartUpLim;

/// the vector of RampUp
std::vector < double >  f_DeltaRampUp;

/// the vector of RampDown
std::vector < double >  f_DeltaRampDown;

/// the vector of QuadTerm
std::vector < double >  f_QuadTerm;

/// the vector of LinearTerm
std::vector < double >  f_LinearTerm;

/// the vector of ConstTerm
std::vector < double >  f_ConstTerm;

/// the FixedConsPower value
double  f_FixedConsPower;

/// the StartUpCost value
double StratUpCost_val;

/// the PZero value
double PZero_val;

/// the MinUpTime value
double MinUpTime_val;

/// the MinDownTime value
double  MinDownTime_val;

/// the InitUpDownTime value
double InitUpDownTime_val;

///< variable denoting the time-steps unit is subjected to initial conditions
int init_t;



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
/*-------------------------- PROTECTED FRIENDS -----------------------------*/
/*--------------------------------------------------------------------------*/

int number_values;

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
* @{ */
/// print the ThermalUnitBlock on an ostream with the given verbosity
/** Protected method to print information about the ThermalUnitBlock;
 * */

virtual void print( std::ostream &output ) const override ;

/*--------------------------------------------------------------------------*/
/// loads the ThermalUnit instance from standard .dat file format
/** Protected method for loading a ThermalUnitBlock out of a std::istream
 *
 *      //TODO
 */

//virtual void load( std::istream &input ) override final;


// void load(std::istream& inStream);

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*-----------------------------variables------------------------------------*/

    ///< the start up binary variables
    std::vector < ColVariable > v_start_up_thermal;

    ///< the shut down binary variables
    std::vector < ColVariable > v_shut_down_thermal;


/*----------------------------constraints-----------------------------------*/

    ///< the connection min up and down time constraints
    std::vector < FRowConstraint > UVW_Const;

     ///< the TURN ON min up and down time constraints
    std::vector < FRowConstraint > UV_Const;

    ///< the SHUT DOWN min up and down time constraints
    std::vector < FRowConstraint > UW_Const;

    ///< the RampUp time constraints
    std::vector < FRowConstraint > RampUp_Const;

    ///< the RampDown time constraints
    std::vector < FRowConstraint > RampDown_Const;

    ///< the PrimaryRho fraction constraints
    std::vector < FRowConstraint > PrimaryRho_Const;

    ///< the SecondaryRho fraction constraints
     std::vector < FRowConstraint > SecondaryRho_Const;

    ///< the PowerOutput constraints

    std::vector < FRowConstraint > PMin_Const;

    std::vector < FRowConstraint > PMax_Const;

    std::vector < FRowConstraint > PowerInjected_Const;

    std::vector <LB0Constraint> PowerFix_Const;

    ///< the (linear) objective function
    FRealObjective object;

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
