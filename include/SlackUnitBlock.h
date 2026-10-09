/*--------------------------------------------------------------------------*/
/*-------------------------- File SlackUnitBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class SlackUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h] and implements a "slack" unit; a (typically, fictitious)
 * unit capable of producing (typically, a large amount of) active power
 * and/or primary/secondary reserve and/or inertia at any time period
 * completely independently from each other and from all other time periods,
 * albeit at a (typically, huge) cost. Such a unit is typically added to a
 * Unit Commitment problem to ensure that it has a (fictitious) feasible
 * solution, which may help solution methods. At the very least such a
 * modified UC would produce a "least unfeasible" solution which can be used
 * to identify the parts of the system that lack capacity/resources.
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

#ifndef __SlackUnitBlock
 #define __SlackUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OneVarConstraint.h"

#include "UnitBlock.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SlackUnitBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the UnitBlock concept for a "slack" unit
/** SlackUnitBlock derives from UnitBlock and implements the concept of
 * "slack" unit, i.e., a (typically fictitious) unit that can give active
 * power, primary and secondary reserve and inertia at each instant,
 * independently of each other and of the other instants, at a (typically very
 * high) cost. Such a unit is added to a unit commitment problem to ensure
 * that it always has a feasible solution, which may help the solution
 * methods. At the optimum, the slack unit is used only where the real units
 * cannot meet the demands, and its use indicates the parts of the system that
 * lack capacity. A slack unit with a negative bound on the active power is a
 * "dump" unit, which absorbs the power that cannot be used otherwise. It has
 * exactly one generator, and its place in the complete model is described in
 * \ref ucblock_model.
 *
 * \par Rows and Objective
 * With \f$ t \in \mathcal{T} = \{ 0 , \ldots , T - 1 \} \f$ and the symbols
 * of deserialize(), the rows are, for the variables that exist,
 * \f{align*}{
 *   & \min\{ P^{mx}_t , 0 \} \leq p^{ac}_t \leq \max\{ P^{mx}_t , 0 \}
 *     \tag{1} \\
 *   & 0 \leq p^{pr}_t \leq P^{pr}_t , \qquad 0 \leq p^{sc}_t \leq P^{sc}_t
 *     \tag{2} \\
 *   & Q^{mn}_t \leq q_t \leq Q^{mx}_t \tag{3} \\
 *   & q_t \leq a^q_t , \qquad - q_t \leq a^q_t \tag{4}
 * \f}
 * i.e., (1) is \f$ 0 \leq p^{ac}_t \leq P^{mx}_t \f$ if
 * \f$ P^{mx}_t \geq 0 \f$ and \f$ P^{mx}_t \leq p^{ac}_t \leq 0 \f$ if
 * \f$ P^{mx}_t < 0 \f$. In (3), \f$ Q^{mn}_t = - Q^{mx}_t \f$ if
 * "MinReactivePower" is absent, and (3)-(4) exist whenever the reactive power
 * does, while an absent bound is 0 (as get_min_reactive_power() and
 * get_max_reactive_power() return it); hence, with no bound the reactive
 * power is 0. Moreover, the commitment \f$ u_t \in [ 0 , 1 ] \f$ is a
 * continuous variable: the unit gives the inertia \f$ h^u_t u_t \f$, any
 * fraction of its largest one \f$ h^u_t \f$. The Objective is
 * \f[
 *   \sum_{ t \in \mathcal{T} } \bigl( b_t p^{ac}_t + c^{pr}_t p^{pr}_t +
 *   c^{sc}_t p^{sc}_t + c^u_t h^u_t u_t + C^q b_t a^q_t \bigr) , \tag{5}
 * \f]
 * where the reactive power costs the constant fraction \f$ C^q = 0.7 \f$
 * (REACTIVE_COST_FACTOR) of the price of the active power, and it is not
 * multiplied by a scale factor (the unit has none, i.e., its scale is 1).
 * Note that the term \f$ b_t p^{ac}_t \f$ is a cost of the power injected if
 * \f$ b_t > 0 \f$ and \f$ P^{mx}_t \geq 0 \f$. For a dump unit
 * (\f$ P^{mx}_t < 0 \f$, hence \f$ p^{ac}_t \leq 0 \f$), instead, the power
 * absorbed costs \f$ - b_t \f$ per unit; hence, a penalty on dumping is given
 * by \f$ b_t < 0 \f$, while \f$ b_t > 0 \f$ would pay the unit for absorbing.
 *
 * \par Features not modeled
 * The slack unit is bounded: "MaxPower" has to be given and large enough to
 * cover any imbalance, since its default 0 means no slack at all, and there
 * is no unbounded slack. Also, its cost is linear, and therefore a cost
 * increasing more than linearly with the imbalance is not represented. */

class SlackUnitBlock : public UnitBlock
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC CONSTANTS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// the cost of the reactive power relative to that of the active power
 /** The factor \f$ C^{q} \f$ that multiplies "ActivePowerCost" in the cost
  * of the absolute value of the reactive power [see generate_objective()]. */

 static constexpr double REACTIVE_COST_FACTOR = 0.7;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father block
 /** Constructor of SlackUnitBlock, taking possibly a pointer of its father
  * Block. */

 explicit SlackUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of SlackUnitBlock

 virtual ~SlackUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the SlackUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The symbols are
 * those of the class description. All the data below are optional and
 * time-indexed: such a variable, of type netCDF::NcDouble, is either a
 * scalar, or indexed over the dimension "NumberIntervals" (one value per
 * interval, see "ChangeIntervals"), or indexed over "TimeHorizon", and it
 * is expanded into one value per instant as described in
 * UnitBlock::deserialize(). The netCDF::NcGroup may contain:
 *
 * - the variable "MaxPower", the datum \f$ P^{mx}_t \f$ of (1), whose sign
 *   gives the direction of the unit: if \f$ P^{mx}_t \geq 0 \f$ the unit
 *   injects at most \f$ P^{mx}_t \f$, if \f$ P^{mx}_t < 0 \f$ it absorbs at
 *   most \f$ - P^{mx}_t \f$ (a "dump" unit); if it is absent,
 *   \f$ P^{mx}_t = 0 \f$ and the unit has no active power at all;
 *
 * - the variables "MaxPrimaryPower" and "MaxSecondaryPower", the
 *   nonnegative bounds \f$ P^{pr}_t \f$ and \f$ P^{sc}_t \f$ of (2); the
 *   reserve variables exist only if the datum is given;
 *
 * - the variable "MaxInertia", the coefficient \f$ h^u_t \geq 0 \f$ of the
 *   commitment \f$ u_t \in [ 0 , 1 ] \f$ in the inertia rows of UCBlock,
 *   i.e., the largest inertia the unit gives; the commitment exists only if
 *   the datum is given;
 *
 * - the variable "ActivePowerCost", the price \f$ b_t \f$ of (5), 0 if
 *   absent (a setting that makes the slack unit as cheap as any other,
 *   while it is normally meant to be used only when nothing else is
 *   possible);
 *
 * - the variables "PrimaryCost" and "SecondaryCost", the costs
 *   \f$ c^{pr}_t \f$ and \f$ c^{sc}_t \f$ of a unit of primary and secondary
 *   reserve in (5), 0 if absent;
 *
 * - the variable "InertiaCost", the cost \f$ c^u_t \f$ in (5) of a unit of
 *   inertia, i.e., the cost of \f$ u_t \f$ is \f$ c^u_t h^u_t \f$; 0 if
 *   absent;
 *
 * - the variables "MinReactivePower" and "MaxReactivePower", the bounds
 *   \f$ Q^{mn}_t \f$ and \f$ Q^{mx}_t \f$ of (3), used only if the
 *   enclosing UCBlock asks for the reactive power; a vector of zeros is
 *   treated as absent. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /* extends UnitBlock::expected_dims()
  * not necessary, no new dimensions

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends UnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the SlackUnitBlock
 /** Generates the static Variable of the SlackUnitBlock, with the symbols
  * of the class description; each group is a std::vector< ColVariable > of
  * size get_time_horizon(), or empty if it is not generated:
  *
  * - the active power \f$ p^{ac}_t \f$ ("p_slack"), continuous, its sign
  *   being given by (1);
  *
  * - the commitment \f$ u_t \in [ 0 , 1 ] \f$ ("u_slack"), continuous,
  *   only if the enclosing UCBlock asks for the inertia [see
  *   set_reserve_vars()] and "MaxInertia" is given [see has_commitment()];
  *
  * - the primary and secondary reserves \f$ p^{pr}_t \f$ ("pr_slack") and
  *   \f$ p^{sc}_t \f$ ("sr_slack"), nonnegative, each only if the enclosing
  *   UCBlock asks for that reserve and "MaxPrimaryPower", respectively
  *   "MaxSecondaryPower", is given;
  *
  * - the reactive power \f$ q_t \f$ ("q_slack"), continuous, and its
  *   absolute value \f$ a^q_t \f$ ("q_a_slack"), nonnegative, only if the
  *   enclosing UCBlock asks for the reactive power [see
  *   set_reactive_power()].
  *
  * The parameter \p stvv is only passed to
  * UnitBlock::generate_abstract_variables(). */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraints of the SlackUnitBlock
 /** Generates the static Constraint of the SlackUnitBlock, i.e., the rows
  * (1)-(4) of the class description:
  *
  * - (1) is "ActivePowerBound_Slack", a vector of BoxConstraint;
  *
  * - (2) is "PrimarySpinningReserveBound_Slack" and
  *   "SecondarySpinningReserveBound_Slack", two vectors of BoxConstraint,
  *   each only if the corresponding reserve exists;
  *
  * - (3) is "ReactivePowerBound_thermal", a vector of BoxConstraint, if the
  *   reactive power exists;
  *
  * - (4) is "Lin_of_Abs_Reactive", a vector of 2 T FRowConstraint (the rows
  *   \f$ q_t - a^q_t \leq 0 \f$ first), if the reactive power exists;
  *
  * - "Inertia_bound_Slack", a vector of ZOConstraint on \f$ u_t \f$, only
  *   if the commitment exists and \p stcc (or, if it is nullptr, the
  *   f_static_constraints_Configuration of the BlockConfig) is a
  *   SimpleConfiguration< int > with a nonzero value: these rows repeat the
  *   bounds of the variables, which some approaches need in order to have a
  *   dual value for them. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the SlackUnitBlock
 /** Generates the Objective (5) of the class description, to be minimized:
  * a LinearFunction with coefficient \f$ b_t \f$ on \f$ p^{ac}_t \f$,
  * \f$ c^{pr}_t \f$ and \f$ c^{sc}_t \f$ on the reserves,
  * \f$ c^{u}_t h^u_t \f$ on \f$ u_t \f$ and \f$ C^q b_t \f$ on
  * \f$ a^q_t \f$, for the variables that exist, an absent cost being 0. The
  * parameter \p objc is not used. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------- Methods for checking the SlackUnitBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the SlackUnitBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this SlackUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this SlackUnitBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this SlackUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or f_BlockConfig->f_is_feasible_Configuration and they are determined as
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
  * - Otherwise, if both f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of the
  *   parameters are obtained analogously as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation
  *   is considered.
  *
  * This function considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed (Constraint::compute()).
  *
  * @param useabstract This parameter is ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


/** @} ---------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE SlackUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the SlackUnitBlock
 * @{ */

 /// returns the upper bound of the active power at instant t
 /** Returns \f$ \max\{ P^{mx}_t , 0 \} \f$, i.e., "MaxPower" if it is
  * nonnegative and 0 otherwise (also if it is absent). The methods below
  * that return the vector of a time-indexed datum return it either empty
  * (the datum is absent) or with one value per instant [see
  * deserialize()]. */

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( ( v_MaxPower.size() > t ) ?
	  ( ( v_MaxPower[ t ] >= 0 ) ? v_MaxPower[ t ] : 0 ) : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the lower bound of the active power at instant t
 /** Returns \f$ \min\{ P^{mx}_t , 0 \} \f$, i.e., "MaxPower" if it is
  * negative and 0 otherwise (also if it is absent). */

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( ( v_MaxPower.size() > t ) ?
	  ( ( v_MaxPower[ t ] >= 0 ) ? 0 : v_MaxPower[ t ] ) : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum reactive power of \p generator at time \p t

 double get_min_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MinReactivePower.size() > t ) ? v_MinReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of \p generator at time \p t

 double get_max_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MaxReactivePower.size() > t ) ? v_MaxReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the bounds of the primary reserve ("MaxPrimaryPower")

 const std::vector< double > & get_max_primary_power( void ) const {
  return( v_MaxPrimaryPower );
  }

/*--------------------------------------------------------------------------*/
 /// returns the prices of the active power ("ActivePowerCost")

 const std::vector< double > & get_active_power_cost( void ) const {
  return( v_ActivePowerCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the bounds of the secondary reserve ("MaxSecondaryPower")

 const std::vector< double > & get_max_secondary_power( void ) const {
  return( v_MaxSecondaryPower );
  }

/*--------------------------------------------------------------------------*/
 /// returns the costs of the primary reserve ("PrimaryCost")

 const std::vector< double > & get_primary_cost( void ) const {
  return( v_PrimaryCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the costs of the secondary reserve ("SecondaryCost")

 const std::vector< double > & get_secondary_cost( void ) const {
  return( v_SecondaryCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the coefficients of the commitment in the inertia rows
 /** Returns a pointer to the get_time_horizon() coefficients \f$ h^u_t \f$
  * ("MaxInertia") of the commitment \f$ u_t \f$ in the inertia rows of
  * UCBlock, nullptr if "MaxInertia" is absent. */

 const double * get_inertia_commitment( Index generator ) const override {
  if( v_MaxInertia.empty() )
   return( nullptr );
  return( &( v_MaxInertia.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the costs of the inertia ("InertiaCost")

 const std::vector< double > & get_inertia_cost( void ) const {
  return( v_InertiaCost );
  }

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Variable OF THE SlackUnitBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the SlackUnitBlock
 *
 * These methods allow to read the each group of Variable that any
 * SlackUnitBlock in principle has (although some may not):
 *
 * - commitment variables;
 *
 * - active_power variables;
 *
 * - primary_spinning_reserve variables;
 *
 * - secondary_spinning_reserve variables;
 * @{ */

 /// returns the vector of commitment variables

 /// the slack unit is committed only if it has an inertia reserve to give
 /** The commitment of a slack unit is there only to say whether the unit is
  * producing inertia reserve, hence it exists only if the enclosing UCBlock
  * asks for the inertia reserve and the unit can produce some. */

 bool has_commitment( void ) const override {
  return( ( reserve_vars & 4u ) && ( ! v_MaxInertia.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

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

 /// the unit provides the reserve only if it has any to give
 /** The enclosing UCBlock asking for the reserve is not enough, the unit has
  * to have some to give: this is the very condition with which the Variable
  * are generated, said in terms of the data alone. */

 bool has_primary_reserve( void ) const override {
  return( ( reserve_vars & 1u ) && ( ! v_MaxPrimaryPower.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_secondary_reserve( void ) const override {
  return( ( reserve_vars & 2u ) && ( ! v_MaxSecondaryPower.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

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
 /// returns the vector of reactive power variables

 ColVariable * get_reactive_power( Index generator ) override {
  if( v_reactive_power.empty() )
   return( nullptr );
  return( &( v_reactive_power.front() ) );
  }

/** @} ---------------------------------------------------------------------*/
/*----------------- METHODS FOR SAVING THE SlackUnitBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the SlackUnitBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * SlackUnitBlock. See SlackUnitBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE SlackUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the SlackUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "SlackUnitBlock::load() not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
/// set the active power cost values (Subset overload)
/** Update v_ActivePowerCost at the time indices in @p subset to the values
 *  pointed by @p values. If the objective has already been generated, the
 *  affected LinearFunction coefficients (the active-power coefficient at the
 *  matching @p t, and the reactive-power-abs coefficient if reactive power
 *  is enabled) are updated in sync. Issues a SlackUnitBlockSbstMod with type
 *  SlackUnitBlockMod::eSetActPCost depending on @p issuePMod. */

 void set_active_power_cost( MF_dbl_it values ,
                             Subset && subset ,
                             bool ordered = false ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/// set the active power cost values (Range overload)
/** Update v_ActivePowerCost at the time indices in @p rng to the values
 *  pointed by @p values. Same Modification dispatch as the Subset overload,
 *  but issues a SlackUnitBlockRngdMod instead. */

 void set_active_power_cost( MF_dbl_it values ,
                             Range rng = Range( 0 , Inf< Index >() ) ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/// set the active power cost to \p value at every instant

 void set_active_power_cost( double value ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector( get_time_horizon() , value );
  set_active_power_cost( vector.cbegin() ,
                         Range( 0 , Inf< Index >() ) ,
                         issuePMod , issueAMod );
  }

/*--------------------------------------------------------------------------*/

 static void static_initialization( void )
 {
  register_method< SlackUnitBlock , MF_dbl_it , Subset && , bool >(
   "SlackUnitBlock::set_active_power_cost" ,
   & SlackUnitBlock::set_active_power_cost );

  register_method< SlackUnitBlock , MF_dbl_it , Range >(
   "SlackUnitBlock::set_active_power_cost" ,
   & SlackUnitBlock::set_active_power_cost );
  }

/*--------------------------------------------------------------------------*/
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

 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the vector of ActivePowerCost
 std::vector< double > v_ActivePowerCost;

  /// the vector of MinReactivePower
 std::vector< double > v_MinReactivePower;

 /// the vector of MaxReactivePower
 std::vector< double > v_MaxReactivePower;

 /// the vector of MaxPrimaryPower
 std::vector< double > v_MaxPrimaryPower;

 /// the vector of PrimaryCost
 std::vector< double > v_PrimaryCost;

 /// the vector of MaxSecondaryPower
 std::vector< double > v_MaxSecondaryPower;

 /// the vector of SecondaryCost
 std::vector< double > v_SecondaryCost;

 /// the vector of MaxInertia
 std::vector< double > v_MaxInertia;

 /// the vector of InertiaCost
 std::vector< double > v_InertiaCost;

/*-------------------------------- variables -------------------------------*/

 /// the commitment variables
 std::vector< ColVariable > v_commitment;

 /// the active power variables
 std::vector< ColVariable > v_active_power;

 /// the reactive power variables
 std::vector< ColVariable > v_reactive_power;

 /// the absolute value of v_reactive_power for the cost function
 std::vector< ColVariable > v_abs_reactive_power;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

/*------------------------------- constraints ------------------------------*/

 /// the active power bound constraints
 std::vector< BoxConstraint > ActivePower_Bound_Const;

 /// the primary spinning reserve bound constraints
 std::vector< LB0Constraint > Primary_Spinning_Reserve_Bound_Const;

 /// the secondary spinning reserve bound constraints
 std::vector< LB0Constraint > Secondary_Spinning_Reserve_Bound_Const;

 /// the inertia variables bound constraints
 std::vector< ZOConstraint > Inertia_Bound_Const;

 /// the reactive power bound constraints
 std::vector< BoxConstraint > ReactivePower_Bound_Const;

 /// Linearization of the v_reactive_power
 std::vector< FRowConstraint > Abs_of_Reactive;

 /*!! Q <= P
 std::vector< FRowConstraint > Reactive_2_Active_Const;
 !!*/

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



};  // end( class( SlackUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SlackUnitBlockMod --------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from UnitBlockMod for changes to a SlackUnitBlock

class SlackUnitBlockMod : public UnitBlockMod
{
 public:

 /// public enum for the types of SlackUnitBlockMod
 enum SUB_mod_type
 {
  eSetActPCost = eUBModLastParam , ///< set active power cost values
  eSUBModLastParam  ///< first allowed parameter for derived classes
  };

 /// constructor, takes the SlackUnitBlock and the type
 SlackUnitBlockMod( SlackUnitBlock * const fblock , const int type )
  : UnitBlockMod( fblock , type ) , f_Block( fblock ) {}

 /// destructor, does nothing
 virtual ~SlackUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the SlackUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "SlackUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetActPCost ):
    output << "Set active power cost values ";
    break;
   default:;
  }
  }

 SlackUnitBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

};  // end( class( SlackUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS SlackUnitBlockRngdMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// derived from SlackUnitBlockMod for "ranged" modifications
class SlackUnitBlockRngdMod : public SlackUnitBlockMod
{
 public:

 /// constructor: takes the SlackUnitBlock, the type, and the range
 SlackUnitBlockRngdMod( SlackUnitBlock * const fblock , const int type ,
                        const Block::Range & rng )
  : SlackUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~SlackUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the SlackUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  SlackUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

};  // end( class( SlackUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS SlackUnitBlockSbstMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// derived from SlackUnitBlockMod for "subset" modifications
class SlackUnitBlockSbstMod : public SlackUnitBlockMod
{
 public:

 /// constructor: takes the SlackUnitBlock, the type, and the subset
 SlackUnitBlockSbstMod( SlackUnitBlock * const fblock , const int type ,
                        Block::Subset && nms )
  : SlackUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~SlackUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the SlackUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  SlackUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

 Block::Subset f_nms;  ///< the subset

};  // end( class( SlackUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __SlackUnitBlock */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SlackUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
