/*--------------------------------------------------------------------------*/
/*----------------------- File IntermittentUnitBlock.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class IntermittentUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define a Unit representing
 * Intermittent Generation in the Unit Commitment Problem.
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

#ifndef __IntermittentUnitBlock
 #define __IntermittentUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColVariable.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "UnitBlock.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS IntermittentUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Intermittent Generation unit
/** The IntermittentUnitBlock class implements the Block concept [see Block.h]
 * for a units representing generation (be ir centralized or distributed) by
 * intermittent (= unreliable) sources in the unit commitment problem, such
 * as wind farms, solar parks and run-of-the-river hydroelectricity. Each
 * unit is supposed to be connected to a specific node of the clustered
 * network (which means that the "distributed" case refers to "distributed in
 * a small region", where of course "small" depends on the granularity of the
 * network description. The model relies mainly on historical data of local
 * generation of wind and solar at each node of the grid; these data are
 * used to develop normalized generation profiles associated with wind and
 * solar generators. Intermittent generators are supposed to be able to
 * contribute to primary and secondary reserves. Contribution to the system
 * inertia concerns more specifically run of river generators. The potential
 * contribution of solar or wind generation to inertia is still the subject
 * of active research. Reserve requirements are specified in order to be
 * symmetrically available to increase or decrease power injected into the
 * grid. Then the technical and physical constraints are mainly divided in
 * three different categories:
 *
 * - the active power bounds;
 *
 * - the maximum power output constraints according to primary and secondary
 *   spinning reserves;
 *
 * - the minimum power output constraints according to primary and secondary
 *   spinning reserves. */

class IntermittentUnitBlock : public UnitBlock
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
 /** Constructor of IntermittentUnitBlock, taking possibly a pointer of its
  * father Block. */

 explicit IntermittentUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of IntermittentUnitBlock

 virtual ~IntermittentUnitBlock() override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the IntermittentUnitBlock. Besides the mandatory "type" attribute of any
  * :Block, the group must contain all the data required by the base
  * UnitBlock, as described in the comments to UnitBlock::deserialize(
  * netCDF::NcGroup ). In particular, we refer to that description for the
  * crucial dimensions "TimeHorizon", "NumberIntervals" and
  * "ChangeIntervals". The netCDF::NcGroup must then also contain:
  *
  * - The variable "MinPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MinP[ t ] that,
  *   for each time instant t, contains the minimum potential production value
  *   of the unit for the corresponding time step. If "MinPower" has length 1
  *   then MinP[ t ] contains the same value for all t. Otherwise, MinPower[ i
  *   ] is the fixed value of MinP[ t ] for all t in the interval [
  *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded. Note that it
  *   must be MnP[ t ] >= 0 for all t.
  *
  * - The variable "MaxPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector MaxP[ t ] that,
  *   for each time instant t, contains the maximum potential production value
  *   of the unit for the corresponding time step. If "MaxPower" has length 1
  *   then MaxP[ t ] contains the same value for all t. Otherwise, MaxPower[ i
  *   ] is the fixed value of MaxP[ t ] for all t in the interval [
  *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
  *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
  *   require "ChangeIntervals", which in fact is not loaded. Note that it
  *   must be MxP[ t ] >= MnP[ t ] [>= 0] for all t. Yet, MxP[ t ] == MnP[ t ]
  *   is possible: it means that (at time instant t) the unit cannot be
  *   curtailed and cannot provide any reserve.
  *
  * - The variable "InertiaPower", of type netCDF::NcDouble and either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   IP[ t ] which, for each time instant t, contains the contribution that
  *   the unit can give to the inertia constraint which depends on the active
  *   power that it is currently generating (basically, the constant to be
  *   multiplied to the active power variable) at time t for this unit. The
  *   variable is optional; if it is not defined, IP[ t ] == 0 for each time
  *   instants t. If it has size 1 then the entry IP[ 0 ] is assumed to
  *   contain the inertia power value for this unit and all time instants t.
  *   Otherwise, InertiaPower[ i ] is the fixed value of IP[ t ] for all t in
  *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
  *   assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
  *   NumberIntervals >= TimeHorizon then the mapping clearly does not require
  *   "ChangeIntervals", which in fact is not loaded.
  *
  * - The scalar variable "Gamma", of type netCDF::NcDouble and not indexed
  *   over any dimension. This variable is used to take into account an
  *   uncertainty on the maximal potential production. Note that it must be
  *   0 <= Gamma <= 1; when Gamma == 0, the unit does not provide any reserve.
  *
  * - The scalar variable "Kappa", of type netCDF::NcDouble and not indexed
  *   over any dimension. This variable is used to multiply to the minimum
  *   and maximum power at each time instant t. This variable is optional, if
  *   it is not provided it is taken to be Kappa == 1. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the IntermittentUnitBlock
 /** The IntermittentUnitBlock class has three different variables which are:
  *
  * - the primary spinning reserve variables;
  *
  * - the secondary spinning reserve variables;
  *
  * - the active power variables.
  *
   of those variables are optional except the active power variables in
  * the sense that the model may just not have them and whenever a group of
  * above variables is created, its size will be the time horizon.
  *
  * In the design scenario of the UC problem, i.e., if an investment cost
  * is given for this IntermittentUnitBlock, an additional binary variable
  * is needed in order to let the model infer how much capacity to install. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraints of the IntermittentUnitBlock
 /** Method that generates the static constraints of the IntermittentUnitBlock.
  * These are the:
  *
  * - maximum and minimum power output constraints according to primary and
  *   secondary spinning reserves are presented in (1)-(2). Each of them is a
  *   std::vector< FRowConstraint >, with the dimension of get_time_horizon(),
  *   where the entry t, for t in \f$ \mathcal{T} = \f$ {0, ...,
  *   get_time_horizon() - 1}, being the maximum and minimum power output value
  *   according to the primary and the secondary spinning reserves at time
  *   t. These constraints ensure the maximum (or minimum) amount of energy
  *   that unit can produce (or use) when it is on (or off).
  *
  *   \f[
  *       p^{pr}_t + p^{sc}_t \leq \gamma ( \kappa P^{mx}_t - p^{ac}_t )
  *                                           \quad t \in \mathcal{T} \quad (1)
  *   \f]
  *
  *   \f[
  *       p^{pr}_t + p^{sc}_t \leq  p^{ac}_t - ( \kappa P^{mn}_t)
  *                                           \quad t \in \mathcal{T} \quad (2)
  *   \f]
  *
  *   where \f$ P^{mx}_t \f$ and \f$ P^{mn}_t \f$ are the maximum and
  *   minimum power output parameters for each time t in \f$ \mathcal{T} \f$,
  *   respectively.
  *
  * - the active power bounds:
  *
  *   \f[
  *    p^{ac}_t \in [ \kappa P^{mn}_t , \kappa P^{mx}_t]
  *                                          \quad t \in \mathcal{T} \quad (3a)
  *   \f]
  *
  *   which, in the design scenario of the UC problem, become:
  *
  *   \f[
  *    x ( \kappa P^{mn}_t ) \leq p^{ac}_t \leq x ( \kappa P^{mx}_t )
  *                                          \quad t \in \mathcal{T} \quad (3b)
  *   \f]
  *
  *   where \f$ x \f$ is the design variable. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the IntermittentUnitBlock
 /** Method that generates the objective of the IntermittentUnitBlock. The
  * objective function of the IntermittentUnitBlock in the design scenario of
  * the UC problem is given as follow:
  *
  * \f[
  *   \min ( I x )
  * \f]
  *
  * where \f$ I \f$ is the investment cost and \f$ x \f$ is the design variable.
  * Otherwise, the objective function of the IntermittentUnitBlock is "empty"
  * (a FRealObjective with a LinearFunction inside with no active variables).
  */

 void generate_objective( Configuration * objc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// setting the BlockConfig
 /** This method sets the BlockConfig of this IntermittentUnitBlock. Besides
  * the Configuration for the is_feasible() function, the
  * IntermittentUnitBlock also considers the extra Configuration of the
  * BlockConfig. If the extra Configuration is a non-null pointer to a
  * SimpleConfiguration< double >, then the value, let us call it epsilon,
  * stored in that Configuration will replace any zero value that may appear
  * as maximum power at any time instant.
  *
  * For instance, if the maximum power provided during deserialization (see
  * IntermittentUnitBlock::deserialize( netCDF::NcGroup )) is zero for some time
  * instant t, then it will become epsilon for that time instant. Moreover, if
  * any zero value is provided to set_maximum_power() for some time instant t,
  * then the maximum power for time instant t will become epsilon.
  *
  * When epsilon > 0, this can be used to prevent the maximum power from being
  * zero. Notice, however, that the actual maximum power may become zero even
  * if epsilon > 0 if the kappa constant is zero (see set_kappa()).
  *
  * The reason behind this is that some Solver may not be able to handle
  * modifications in the maximum power if it is initially zero and become
  * nonzero after a modification. By setting epsilon > 0, this issue is
  * avoided.
  *
  * Please see the comments to Block::set_BlockConfig() for more details about
  * the BlockConfig. */

 void set_BlockConfig( BlockConfig * newBC = nullptr ,
                       bool deleteold = true ) override;

/**@} ----------------------------------------------------------------------*/
/*------------- Methods for checking the IntermittentUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the IntermittentUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this IntermittentUnitBlock is
  * approximately feasible within the given tolerance. That is, a solution is
  * considered feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this IntermittentUnitBlock is not
  *      greater than the tolerance.
  *
  * Every Constraint of this IntermittentUnitBlock is a RowConstraint and its
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
/*------- METHODS FOR READING THE DATA OF THE IntermittentUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the IntermittentUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of Intermittent Generation units
 * @{ */

 /// returns the gamma value
 double get_gamma( void ) const { return( f_gamma ); }

 /// returns the kappa value
 double get_kappa( void ) const { return( f_kappa ); }

 /// returns the investment cost
 double get_investment_cost( void ) const { return( f_InvestmentCost ); }

 /// returns the maximum installable capacity by the user
 double get_max_capacity( void ) const { return( f_MaxCapacity ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of minimum power
 /** The method returned a std::vector< double > V and each element of V
  * contains the minimum power at time t. There are three possible cases:
  *
  * - if the vector is empty, then the minimum power of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the minimum power of
  *   the unit for all time horizon;
  *
  * - otherwise, the std::vector< double > V must have size get_time_horizon()
  *   and each V[ t ] represents the minimum power value at time t. */

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( v_MinPower[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum power
 /** The method returned a std::vector< double > V and each element of V
  * contains the maximum power at time t. There are three possible cases:
  *
  * - if the vector is empty, then the maximum power of the unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] is the maximum power of
  *   the unit for all time horizon;
  *
  * - otherwise, the std::vector< double > V must have size get_time_horizon()
  *   and each V[ t ] represents the maximum power value at time t. */

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower[ t ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inertia power
 /** The returned value U = get_inertia_power() contains the contribution to
  * inertia (basically, the constants to be multiplied by the active power
  * variables returned by get_active_power()) of all the generators at all
  * time instants. There are four possible cases:
  *
  * - if the matrix is empty, then the inertia power is always 0 and this
  *   function returns nullptr;
  *
  * - if the matrix only has one row (i.e., the first dimension has size 1),
  *   then the inertia power for each generator g is U[ 0 , g ] for all t
  *   which means that the second dimension has size get_number_generators();
  *
  * - if the matrix only has one column with size get_time_horizon() (i.e.,
  *   the second dimension has size 1), then the InertiaPower[ t , 0 ] gives
  *   the inertia power for the problem at time t. Since in this unit there is
  *   only one electrical generator, this case should happen by assumption;
  *
  * - otherwise, the matrix has size get_time_horizon() per
  *   get_number_generators(), then the InertiaPower[ t , g ] represents the
  *   inertia power for the problem at time t for each electrical generator
  *   g. */

 const double * get_inertia_power( Index generator ) const override {
  if( v_InertiaPower.empty() )
   return( nullptr );
  return( &( v_InertiaPower.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the scale factor

 double get_scale( void ) const override { return( f_scale ); }

/**@} ----------------------------------------------------------------------*/
/*------ METHODS FOR READING THE Variable OF THE IntermittentUnitBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the IntermittentUnitBlock
 *
 * These methods allow to read the each group of Variable that any
 * IntermittentUnitBlock in principle has (although some may not):
 *
 * - active_power variables;
 *
 * - primary_spinning_reserve variables;
 *
 * - secondary_spinning_reserve variables.
 * @{ */

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
 /// returns the design variable

 ColVariable & get_design( void ) {
  return( design );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum total power constraints

 const std::vector< FRowConstraint > & get_min_power_constraints( void ) const {
  return( min_power_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum total power constraints

 const std::vector< FRowConstraint > & get_max_power_constraints( void ) const {
  return( max_power_Const );
 }

/*--------------------------------------------------------------------------*/
 /// returns the bound constraints on the active power

 const std::vector< BoxConstraint > &
 get_active_power_bound_constraints( void ) const {
  return( active_power_bounds_Const );
 }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR SAVING THE IntermittentUnitBlock---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the IntermittentUnitBlock
 * @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * IntermittentUnitBlock. See IntermittentUnitBlock::deserialize(
 * netCDF::NcGroup ) for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*----------- METHODS FOR INITIALIZING THE IntermittentUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the IntermittentUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error(
   "IntermittentUnitBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_maximum_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_maximum_power( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
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
 /// set the kappa constant
 /** This function sets the kappa constant, which multiplies the minimum and
  * maximum power in the constraints of this IntermittentUnitBlock.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param subset If non-empty, the kappa constant is set to the value pointed
  *               by \p values. If empty, no operation is performed.
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
  * maximum power in the constraints of this IntermittentUnitBlock.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param rng If non-empty, the kappa constant is set to the value pointed by
  *            \p values. If empty, no operation is performed.
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
  * maximum power in the constraints of this IntermittentUnitBlock.
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

 /// the vector of MinPower
 std::vector< double > v_MinPower;

 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the matrix of inertia power of generators
 std::vector< double > v_InertiaPower;


 /// the investment cost
 double f_InvestmentCost{};

 /// the maximum installable capacity by the user
 double f_MaxCapacity{};

 /// the gamma value
 double f_gamma{};

 /// the kappa value
 double f_kappa = 1;

 /// the scale factor
 double f_scale = 1;

 /// this is the value that will replace any zero value in maximum power
 double f_max_power_epsilon{};

/*-------------------------------- variables -------------------------------*/

 /// the active power variables
 std::vector< ColVariable > v_active_power;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

 /// the design variable
 ColVariable design;

/*------------------------------- constraints ------------------------------*/

 /// the active power lower bound constraints
 std::vector< FRowConstraint > min_power_Const;

 /// the active power upper bound constraints
 std::vector< FRowConstraint > max_power_Const;


 /// the active power bounds design constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_design_Const;


 /// the active power bounds constraints
 std::vector< BoxConstraint > active_power_bounds_Const;


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

 /// updates the constraints for the current maximum power
 /** This function updates the right-hand side of the "maximum power" and the
  * "active power bounds" constraints associated with the time instants given
  * in \p time. */

 void update_max_power_in_cnstrs( const Block::Subset & time ,
                                  c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current maximum power
 /** This function updates the right-hand side of the "maximum power" and the
  * "active power bounds" constraints associated with the time instants given
  * in \p time. */

 void update_max_power_in_cnstrs( const Block::Range & time ,
                                  c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this IntermittentUnitBlock is consistent
 /** This function checks whether the data in this IntermittentUnitBlock is
  * consistent. The data is consistent if all of the following conditions are
  * met.
  *
  * - The maximum power is greater than or equal to the minimum power.
  *
  * - The minimum power is nonnegative.
  *
  * - Gamma is between 0 and 1.
  *
  * - Kappa is nonnegative.
  *
  * - The inertia power is nonnegative.
  *
  * If any of the above conditions are not met, an exception is thrown. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/

 static void static_initialization( void ) {
  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args(),
   *
   * register_method< IntermittentUnitBlock >
   *                ( "IntermittentUnitBlock::set_maximum_power",
   *                  &IntermittentUnitBlock::set_maximum_power,
   *                  MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::set_maximum_power" ,
   &IntermittentUnitBlock::set_maximum_power );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::set_maximum_power" ,
   &IntermittentUnitBlock::set_maximum_power );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::scale" ,
   &IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::scale" ,
   &IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::set_kappa" ,
   &IntermittentUnitBlock::set_kappa );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::set_kappa" ,
   &IntermittentUnitBlock::set_kappa );
 }

};  // end( class( IntermittentUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS IntermittentUnitBlockMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a IntermittentUnitBlock
class IntermittentUnitBlockMod : public UnitBlockMod
{

 public:

 /// public enum for the types of IntermittentUnitBlockMod
 enum IUB_mod_type
 {
  eSetMaxP = eUBModLastParam , ///< set max power values
  eSetKappa ,                  ///< set the kappa constant
  eIUBModLastParam             ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to
   * extend the set of types of IntermittentUnitBlockMod. */
 };

 /// constructor, takes the IntermittentUnitBlock and the type
 IntermittentUnitBlockMod( IntermittentUnitBlock * const fblock ,
                           const int type )
  : UnitBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~IntermittentUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the IntermittentUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "IntermittentUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   default:
    output << "Set max power values ";
  }
 }

 IntermittentUnitBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

};  // end( class( IntermittentUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------- CLASS IntermittentUnitBlockRngdMod -------------------*/
/*--------------------------------------------------------------------------*/

/// derived from IntermittentUnitBlockMod for "ranged" modifications
class IntermittentUnitBlockRngdMod : public IntermittentUnitBlockMod
{

 public:

 /// constructor: takes the IntermittentUnitBlock, the type, and the range
 IntermittentUnitBlockRngdMod( IntermittentUnitBlock * const fblock ,
                               const int type , Block::Range rng )
  : IntermittentUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~IntermittentUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the IntermittentUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  IntermittentUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng;  ///< the range

};  // end( class( IntermittentUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------- CLASS IntermittentUnitBlockSbstMod -------------------*/
/*--------------------------------------------------------------------------*/

/// derived from IntermittentUnitBlockMod for "subset" modifications
class IntermittentUnitBlockSbstMod : public IntermittentUnitBlockMod
{

 public:

 /// constructor: takes the IntermittentUnitBlock, the type, and the subset
 IntermittentUnitBlockSbstMod( IntermittentUnitBlock * const fblock ,
                               const int type , Block::Subset && nms )
  : IntermittentUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~IntermittentUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the IntermittentUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  IntermittentUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

};  // end( class( IntermittentUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* IntermittentUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File IntermittentUnitBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
