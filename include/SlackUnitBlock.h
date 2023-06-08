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
/** The SlackUnitBlock class derives from UnitBlock and implements the concept
 * of "slack" unit; a (typically, fictitious) unit capable of producing
 * (typically, a large amount of) active power and/or primary/secondary
 * reserve and/or inertia at any time period completely independently from
 * each other and from all other time periods, albeit at a (typically, huge)
 * cost. Such a unit is typically added to a Unit Commitment problem to ensure
 * that it has a (fictitious) feasible solution, which may help solution
 * methods. At the very least such a modified UC would produce a "least
 * unfeasible" solution which can be used to identify the parts of the system
 * that lack capacity/resources. */

class SlackUnitBlock : public UnitBlock
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
 /** Constructor of SlackUnitBlock, taking possibly a pointer of its father
  * Block. */

 explicit SlackUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of SlackUnitBlock

 virtual ~SlackUnitBlock() override;

/**@} ----------------------------------------------------------------------*/
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
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The variable "MaxPower", of type netCDF::NcDouble and either of size 1
 *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
 *   not provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector MxP[ t ] that, for
 *   each time instant t, contains the maximum active power output value of
 *   the unit for the corresponding time step. If "MaxPower" has length 1
 *   then MxP[ t ] contains the same value for all t. Otherwise, MaxPower[ i ]
 *   is the fixed value of MxP[ t ] for all t in the interval [
 *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. This variable is optional, if is not
 *   provided then MxP[ t ] == 0 for all t. Note that it must be MxP[ t ] >= 0
 *   for all t. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon,
 *   then the mapping clearly does not require "ChangeIntervals", which in
 *   fact is not loaded.
 *
 * - The variable "MaxPrimaryPower", of type netCDF::NcDouble and either of
 *   size 1 or indexed over the dimension "NumberIntervals" (if
 *   "NumberIntervals" is not provided, then this variable can also be
 *   indexed over "TimeHorizon"). This is meant to represent the vector
 *   MaxPP[ t ] that, for each time instant t, contains the maximum amount of
 *   primary reserve that the unit can produce in the corresponding time step.
 *   If "MaxPrimaryPower" has length 1 then MaxPP[ t ] contains the same value
 *   for all t. Otherwise, MaxPrimaryPower[ i ] is the fixed value of MaxPP[ t
 *   ] for all t in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[
 *   i ] ], with the assumption that ChangeIntervals[ - 1 ] = 0. This variable
 *   is optional, if is not provided then MaxPP[ t ] == 0 for all t. If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxSecondaryPower", of type netCDF::NcDouble and either of
 *   size 1 or indexed over the dimension "NumberIntervals" (if
 *   "NumberIntervals" is not provided, then this variable can also be
 *   indexed over "TimeHorizon"). This is meant to represent the vector
 *   MaxSP[ t ] that, for each time instant t, contains the maximum amount of
 *   secondary reserve that the unit can produce in the corresponding time
 *   step. If "MaxSecondaryPower" has length 1 then MaxSP[ t ] contains the
 *   same value for all t. Otherwise, MaxSecondaryPower[ i ] is the fixed
 *   value of MaxSP[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
 *   0. This variable is optional, if is not provided then MaxSP[ t ] == 0 for
 *   all t. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
 *   the mapping clearly does not require "ChangeIntervals", which in fact is
 *   not loaded.
 *
 * - The variable "MaxInertia", of type netCDF::NcDouble and either of size 1
 *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals" is
 *   not provided, then this variable can also be indexed over "TimeHorizon").
 *   This is meant to represent the vector MaxI[ t ] which, for each time
 *   instant t, contains the maximum "amount of inertia" (contribution that
 *   the SlackUnit can give to the inertia constraint) at time t. The
 *   variable is optional; if it is not defined, MaxI[ t ] == 0
 *   for all time instants. If it has size 1, then MaxI[ t ] == MaxInertia[ 0
 *   ] for all t, regardless to what "NumberIntervals" says. Otherwise,
 *   MaxInertia[ i ] is the fixed value of MaxI[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "ActivePowerCost", of type netCDF::NcDouble and either of
 *   size 1 or indexed over the dimension "NumberIntervals" (if
 *   "NumberIntervals" is not provided, then this variable can also be
 *   indexed over "TimeHorizon"). This is meant to represent the vector
 *   APC[ t ] that, for each time instant t, contains the cost of producing
 *   one unit of active power at the corresponding time step. This variable
 *   is optional, if it is not provided then it's taken to be zero (although
 *   this is a very strange setting, as it would typically imply that all the
 *   demand, or at least as much as possible of it, is satisfied by the
 *   fictitious SlackUnit rather than from "real" ones). If "ActivePowerCost"
 *   has length 1 then APC[ t ] contains the same value for t. Otherwise,
 *   ActivePowerCost[ i ] is the fixed value of APC[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ] with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "PrimaryCost", of type netCDF::NcDouble and either of size
 *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
 *   is not provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector PC[ t ] that, for
 *   each time instant t, contains the cost of producing one unit of primary
 *   reserve at the corresponding time step. This variable is optional; if it
 *   is not provided then it's taken to be zero (but this is a very strange
 *   setting, cf. the discussion in ActivePowerCost). If "PrimaryCost" has
 *   length 1 then PC[ t ] contains the same value for t. Otherwise,
 *   PrimaryCost[ i ] is the fixed value of PC[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "SecondaryCost", of type netCDF::NcDouble and either of
 *   size 1 or indexed over the dimension "NumberIntervals" (if
 *   "NumberIntervals" is not provided, then this variable can also be
 *   indexed over "TimeHorizon"). This is meant to represent the vector
 *   SC[ t ] that, for each time instant t, contains the cost of producing
 *   one unit of secondary reserve at the corresponding time step. This
 *   variable is optional; if it is not provided then it's taken to be zero
 *   (but this is a very strange setting, cf. the discussion in
 *   ActivePowerCost). If "SecondaryCost" has length 1 then SC[ t ] contains
 *   the same value for t. Otherwise, SecondaryCost[ i ] is the fixed value
 *   of SC[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "InertiaCost", of type netCDF::NcDouble and either of size
 *   1 or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
 *   is not provided, then this variable can also be indexed over
 *   "TimeHorizon"). This is meant to represent the vector IC[ t ] that, for
 *   each time instant t, contains the "the cost of producing "one unit" of
 *   inertia. Since the inertia-producing variable is u[ t ] which is in the
 *   interval [ 0 , 1 ], the cost of u[ t ] is MaxI[ t ] * IC[ t ]; in other
 *   words, u[ t ] represents the fraction the maximum possible amount of
 *   inertia (MaxI[ t ]) that can be produced at time step t This variable is
 *   optional; if it is not provided then it's taken to be zero (but this is a
 *   very strange setting, cf. the discussion in ActivePowerCost). If
 *   "InertiaCost" has length 1 then IC[ t ] contains the same value for
 *   t. Otherwise, InertiaCost[ i ] is the fixed value of IC[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the SlackUnitBlock
 /** The SlackUnitBlock class has several different variables which are:
  *
  * - the commitment variables which takes the continues values between
  *   1 and zero.
  *
  * - the primary spinning reserve variables;
  *
  * - the secondary spinning reserve variables;
  *
  * - the active power variables.
  *
  * Note that of these variables are optional, and it is also possible to
  * restrict which of the subsets are generated without using the parameter
  * stvv. In other word, each group of above variables as binary commitment, or
  * primary or secondary spinning reserve, or active power variables must be
  * defined if and only if the MaxInertia or MaxPrimaryPower or
  * MaxSecondaryPower or MaxPower is defined in the
  * deserialize( netCDF::NcGroup ) respectively. Otherwise, the corresponding
  * variable is not to be needed to generate. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
/// Generate the static constraint of the SlackUnitBlock
/** This method generates the abstract constraints of the SlackUnitBlock.
 *
 * The operations of the slack generating unit are described on a discrete
 * time horizon as dictated by the UnitBlock interface. In this description
 * we indicate it with \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$.
 * This unit just contains the bounds constraint on the ActivePower, Primary
 * and Secondary spinning reserve variables as below:
 *
 * \f[
 *      0 \leq p^{ac}_t \leq P^{mx}_t \quad t \in \mathcal{T}         \quad (1)
 * \f]
 *
 * \f[
 *      0 \leq p^{pr}_t \leq P^{mxP}_t \quad t \in \mathcal{T}        \quad (2)
 * \f]
 *
 * \f[
 *      0 \leq p^{sc}_t \leq P^{mxS}_t \quad t \in \mathcal{T}        \quad (3)
 * \f]
 *
 * Note that the inertia is "produced" by the commitment variable u_t, which
 * is consider as a kPosUnitary and therefore has "implicit" lower and upper
 * bounds 0 and 1, and thus it does not need BoxConstraint. However, other
 * variables are restricted by a BoxConstraint as defined above. The first
 * group is about active power bounds the second one is about the primary
 * spinning reserve, and the last one for the secondary spinning reserve
 * variables. */

 void generate_abstract_constraints( Configuration *stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the SlackUnitBlock
 /** Method that generates the objective of the SlackUnitBlock.
  *
  * - Objective function: the objective function of the SlackUnitBlock
  *   representing the total power production cost to be minimized has the
  *   form:
  *
  *   \f[
  *     \min ( \sum_{ t \in  [t_0 , \mathcal{T}] } C^{ac}_t p^{ac}_t +
  *     C^{pr}_t p^{pr}_t +C^{sc}_t p^{sc}_t +
  *     (P^{MaxI}_t * C^{i}_t) u_t )
  *   \f]
  *
  *   where \f$ C^{ac}_t \f$, \f$ C^{pr}_t \f$, \f$ C^{sc}_t \f$,
  *   \f$ C^{i}_t \f$, and \f$ P^{MaxI}_t \f$ for each time step t are
  *   defined as the ActivePowerCost, PrimaryCost, and SecondaryCost,
  *   InertiaCost, and the MaxInertia respectively.
  *
  *  The objective of the SlackUnitBlock would seem to be an exceedingly simple
  *  object, there is still a nontrivial decision to be made about it, and it
  *  is also possible to restrict. If objc is not nullptr and it is a
  *  SimpleConfiguration< double > or if
  *  f_BlockConfig->f_objective_Configuration is not nullptr and it is a
  *  SimpleConfiguration< double >, then the f_value of the
  *  SimpleConfiguration< int > is taken the objective function. */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
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
/*----------- METHODS FOR READING THE DATA OF THE SlackUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the SlackUnitBlock
 * @{ */

/// returns the vector of maximum power
/** The returned vector contains to maximum power at time t. There are three
 * possible cases:
 *
 * - if the vector is empty, then the maximum power of the unit is 0;
 *
 * - if the vector has only one element, then the maximum power of the unit
 *   for all time horizon;
 *
 * - otherwise, the vector must have size get_time_horizon() and each element
 *   of vector represents the maximum power value at time t. */

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum primary power
 /** The returned vector contains to maximum primary power at time t. There are
  * three possible cases:
  *
  * - if the vector is empty, then the maximum primary power of the unit is 0;
  *
  * - if the vector has only one element, then the maximum primary power of
  *   the unit for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the maximum primary power value at time t. */

 const std::vector< double > & get_max_primary_power( void ) const {
  return( v_MaxPrimaryPower );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active power cost
 /** The returned vector contains to active power cost at time t. There are
  * three possible cases:
  *
  * - if the vector is empty, then the active power cost of the unit is 0;
  *
  * - if the vector has only one element, then the active power cost of the
  *   unit for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the active power cost value at time t. */

 const std::vector< double > & get_active_power_cost( void ) const {
  return( v_ActivePowerCost );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum secondary power
 /** The returned vector contains to maximum secondary power at time t. There
  * are three possible cases:
  *
  * - if the vector is empty, then the maximum secondary power of the unit is
  *   0;
  *
  * - if the vector has only one element, then the maximum secondary power of
  *   the unit for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the maximum secondary power value at time t. */

 const std::vector< double > & get_max_secondary_power( void ) const {
  return( v_MaxSecondaryPower );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary cost
 /** The returned vector contains to primary cost at time t. There are three
  * possible cases:
  *
  * - if the vector is empty, then the primary cost of the unit is 0;
  *
  * - if the vector has only one element, then the primary cost of the unit for
  *   all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the primary cost value at time t. */

 const std::vector< double > & get_primary_cost( void ) const {
  return( v_PrimaryCost );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary cost
 /** The returned vector contains to secondary cost at time t. There are three
  * possible cases:
  *
  * - if the vector is empty, then the secondary cost of the unit is 0;
  *
  * - if the vector has only one element, then the secondary cost of the unit
  *   for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the secondary cost value at time t. */

 const std::vector< double > & get_secondary_cost( void ) const {
  return( v_SecondaryCost );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inertia commitment
 /** The returned value U = get_inertia_commitment() contains the contribution
  * to inertia (basically, the constants to be multiplied by the commitment
  * variables returned by get_commitment()) of all the generators at all time
  * instants. There are three possible cases:
  *
  * - if the vector is empty, then the inertia commitment is always 0 and this
  *   function returns nullptr;
  *
  * - if the vector only has one element, then the inertia commitment for the
  *   fixed consumption of the unit for all t;
  *
  * - otherwise, the vector must have size get_time_horizon(), and each element
  *   of vector represents the inertia commitment at time t. */

 const double * get_inertia_commitment( Index generator ) const override {
  if( v_MaxInertia.empty() )
   return( nullptr );
  return( &( v_MaxInertia.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of inertia cost
 /** The returned vector contains to inertia cost at time t. There are three
  * possible cases:
  *
  * - if the vector is empty, then the inertia cost of the unit is 0;
  *
  * - if the vector has only one element, then the inertia cost of the unit
  *   for all time horizon;
  *
  * - otherwise, the vector must have size get_time_horizon() and each element
  *   of vector represents the inertia cost value at time t. */

 const std::vector< double > & get_inertia_cost( void ) const {
  return( v_InertiaCost );
 }

/**@} ----------------------------------------------------------------------*/
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

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

/*------------------------------- constraints ------------------------------*/

 /// the active power bound constraints
 std::vector< LB0Constraint > ActivePower_Bound_Const;

 /// the primary spinning reserve bound constraints
 std::vector< LB0Constraint > Primary_Spinning_Reserve_Bound_Const;

 /// the secondary spinning reserve bound constraints
 std::vector< LB0Constraint > Secondary_Spinning_Reserve_Bound_Const;


 /// the inertia variables bound constraints
 std::vector< ZOConstraint > Inertia_Bound_Const;


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
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SlackUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SlackUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
