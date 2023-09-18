/*--------------------------------------------------------------------------*/
/*------------------------- File NuclearUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class NuclearUnitBlock, which derives from
 * ThermalUnitBlock [see ThermalUnitBlock.h] and therefore from UnitBlock
 * [see UnitBlock.h], in order to define a "nuclear unit", i.e., "a thermal
 * unit subject to modulation constraints".
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NuclearUnitBlock
 #define __NuclearUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*----------------------- CLASS NuclearUnitBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a nuclear (thermal) unit
/** The NuclearUnitBlock class derives from ThermalUnitBlock (and therefore
 * from UnitBlock) and implements the concept of a "nuclear unit", i.e.,
 * a thermal unit (as implemented in ThermalUnitBlock) subject to extra
 * "modulation constraints" that limit the number of time instants in which
 * the energy production can (significantly) change.
 *
 * The class makes some assumptions about the behavour and implementation of
 * the base class:
 *
 * - whatever formulation is implemented in ThermalUnitBlock, it comprises
 *   the start-up and shut-down variables (v_t and w_t in the standard
 *   notation), where
 *
 *    v_t = v_start_up[ t - init_t ]  ,  w_t = v_shut_down[ t - init_t ]
 *
 *   with init_t a field of the class that indicates the first time instant
 *   in which the commitment variables are free to be chosen (while in all
 *   instants 0 <= t < init_t, if any, they are fixed to either 0 or 1
 *   depending on the state of the unit prior to the beginning of time
 *   (instant 0) due to the minumum up- or down-time constraints
 *
 * - variables_generated() and constraints_generated() can be used to assess
 *   whether or not (static) Variable and Constraint have already been
 *   generated
 *
 * - standard ramp-up and ramp-down deltas are defined, i.e.,
 *   v_DeltaRampUp and v_DeltaRampDown are nonempty
 *
 * Current significant limitations of the class are:
 *
 * - the duration of the modulation period and the initial state of
 *   modulation cannot be changed at all
 *
 * - the modulation ramp-up and ramp-down deltas can be changed, but only
 *   from the physical representation: doing that from the abstract
 *   representation is not allowed (it requires changing coefficients in
 *   the corresponding constraints, which is not supported already in
 *   ThermalUnitBlock, so one would get an exception)
 */

class NuclearUnitBlock : public ThermalUnitBlock {

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
 /** Constructor of NuclearUnitBlock, possibly taking a pointer of its
  * father Block (presumably, but not necessarily, a UCBlock). */

 explicit NuclearUnitBlock( Block * f_block = nullptr ) :
  ThermalUnitBlock( f_block ) , f_initial_modulation( 0 ) ,
  f_modulation_interval( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of NuclearUnitBlock

 virtual ~NuclearUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the NuclearUnitBlock, i.e., that of ThermalUnitBlock plus the extra
 * information needed to handle modulation constraints:
 *
 * - The scalar variable "ModulationTime" of type netCDF::NcUint for the
 *   modulationinterval, i.e., the number of consecutive time instants
 *   within which at most one modulation is allowed. The value must be >= 2,
 *   as otherwise the modulation constraint is useless and one can use an
 *   original ThermalUnitBlock. Note that, like with "MinUpTime" and
 *   "MinDownTime", there is no way to change the value of "ModulationTime"
 *   once the NuclearUnitBlock is loaded, so useless once, useless always.
 *   This variable is optional. If it is not provided, then the value is
 *   assumed to be 2.
 *
 * - The scalar variable "InitModulation" of type netCDF::NcUint for the
 *   initial modulation state, i.e., the number of instants before the first
 *   (0) when the unit last modulated. The value must be >= 1, as a value of
 *   0 would mean that the unit modulated at the first interval, which is
 *   nonsensical since it may be forced to be off, and even if it is on is
 *   should be free to decide whether or not modulating then: "the history
 *   is written only for instants in the past, i.e, before 0". This variable
 *   is optional. If it is not provided, then  the value is assumed to be
 *   equal to that of "ModulationTime", which means that lhe unit last
 *   modulated "a long time ago" (possibly never) and it is allowed to start
 *   modulating immediately (provided it is on, see "InitUpDownTime" in
 *   the original ThermalUnitBlock).
 *
 * - The variable "ModulationDeltaRampUp" of type double and either of size 1
 *   or indexed over the dimension "NumberIntervals"; if "NumberIntervals" is
 *   not provided, then this variable can also be indexed over "TimeHorizon".
 *   This is meant to represent the vector MDP[ t ] that, for each time
 *   instant t, contains the modulation ramp-up value of the unit for the
 *   corresponding time step, i.e., the maximum possible increase of active
 *   power production w.r.t. the power that had been produced in time instant
 *   t - 1 *unless a modulation is performed* (note that start-ups are not
 *   modulations, hence for a modulation to be performed the unit necessarily
 *   had to be producing power at t - 1). This variable is optional; if it is
 *   not provided then it is assumed that MDP[ t ] == 0, i.e., the unit cannot
 *   increase its power output unless a modulation is perofrmed. Note that the
 *   value of the variable must always be such that 0 <= MDP[ t ] <= DP[ t ],
 *   where DP[ t ] is the vector of original ramp-up values fot the unit (see
 *   "DeltaRampUp" in the original ThermalUnitBlock). If
 *   "ModulationDeltaRampUp" has length 1, then MDP[ t ] contains the same
 *   value for all t. Otherwise, ModulationDeltaRampUp[ i ] is the fixed
 *   value of MDP[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *    = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
 *   the mapping clearly does not require "ChangeIntervals", which in fact
 *   is not loaded.
 *
 * - The variable "ModulationDeltaRampDown" of type double and either of size
 *   1 or indexed over the dimension "NumberIntervals"; if "NumberIntervals"
 *   is not provided, then this variable can also be indexed over
 *   "TimeHorizon". This is meant to represent the vector MDM[ t ] that, for
 *   each time instant t, contains the modulation ramp-down value of the unit
 *   for the corresponding time step, i.e., the maximum possible decrease of
 *   active power production w.r.t. the power that had been produced in time
 *   instant t - 1 *unless a modulation is performed* (note that start-ups
 *   are not modulations, hence for a modulation to be performed the unit
 *   necessarily had to be producing power at t - 1). This variable is
 *   optional; if it is not provided then it is assumed that MDM[ t ] == 0,
 *   i.e., the unit cannot decrease its power output unless a modulation is
 *   perofrmed. Note that the value of the variable must always be such that
 *   0 <= MDM[ t ] <= DM[ t ], where DM[ t ] is the vector of original
 *   ramp-down values fot the unit (see "DeltaRampDown" in the original
 *   ThermalUnitBlock). If "ModulationDeltaRampDown" has length 1, then
 *   MDM[ t ] contains the same value for all t. Otherwise,
 *   ModulationDeltaRampDown[ i ] is the fixed value of MDM[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ]  = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * Important: nuclear units *must* have ramp-up and ramp-down constraints,
 * hence DeltaRampUp and DeltaRampDown must be provided in the \p group
 * (although they are in principle optional for ThermalUnitBlock). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the NuclearUnitBlock
/** Besides those of ThermalUnitBlock (depending on the exact chosen
 * formulation, see the Configuration parameter in the ThermalUnitBlock
 * class), NuclearUnitBlock has just one more vector of variables: the
 * std::vector of size TimeHorizon containing the binary variables m[ t ]
 * indicating whether or not the unit is performing a modulation at time
 * instant t. These variables are mandatory.
 *
 * The Configuration parameter has no use here and it is just passed up to
 * the method of the base ThermalUnitBlock class that is called first thing
 * here inside. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
/// Generate the static constraint of the NuclearUnitBlock
/** Besides those of ThermalUnitBlock (depending on the exact chosen
 * formulation, see the Configuration parameter in the ThermalUnitBlock
 * class), this method generates the abstract constraints of the
 * NuclearUnitBlock.
 *
 * The method constructs five groups of constraints, oll of which are
 * mandatory and instantiated for all time instants t:
 *
 * - the modulation ramp-up constraint
 *   \f[
 *     p_t - p_{t-1} \leq \Delta^M_{t+} u_{t-1} +
 *     ( \Delta_{t+} - \Delta^M_{t+} ) m_t + \bar{l}_t v_t
 *   \f]
 *
 *   with \f$ \Delta^M_{t+} \f$ the modulation ramp-up value at time t,
 *   \f$ \Delta_{t+} \f$ the original ramp-up value at time t, and
 *   \f$ \bar{l}_t \f$ the start-up limit at time t
 *
 * - the modulation ramp-down constraint
 *   \f[
 *     p_{t-1} - p_t \leq \Delta^M_{t-} u_t +
 *     ( \Delta_{t-} - \Delta^M_{t-} ) m_t + \bar{u}_t w_t
 *   \f]
 *
 *   with \f$ \Delta^M_{t-} \f$ the modulation ramp-down value at time t,
 *   \f$ \Delta_{t-} \f$ the original ramp-down value at time t, and
 *   \f$ \bar{u}_t \f$ the shut-down limit at time t
 *
 * - the logical constraints \f$ m_t \leq u_t \f$ (modulations cannot
 *   happen when the unit is down)
 *
 * - the logical constraints \f$ m_t \leq ( 1 - v_t ) \f$ (start-ups are
 *   not modulations, so the unit cannot modulate while starting up=
 *
 * - the modulation constraint proper
 *   \f[
 *     \sum_{h = \max\{ 0 , t - \tau^M + 1 \}}^t m_h \leq 1
 *   \f]
 *
 *   where \f$ \tau^M \f$ is the modulation interval.
 *
 * Note that another group of constraints is needed to fix to 0 certain
 * modulation variables depending on the initial modulation value; this is
 * done by changing the bounds (and fixing them) and therefore it does not
 * result in a separate group of constraints.
 *
 * The Configuration parameter has no use here and it is just passed up to
 * the method of the base ThermalUnitBlock class that is called first thing
 * here inside. */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// generate the Objective of the NuclearUnitBlock
/* Not necessary, NuclearUnitBlock does not change the Objective of 
 * ThermalUnitBlock. */

// void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- Methods for checking the NuclearUnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the NuclearUnitBlock
 *  @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this NuclearUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  *   -# is is (approximately) feasible for the original ThermalUnitBlock;
  *
  *   -# each new ColVariable of NuclearUnitBlock is feasible;
  *
  *   -# the violation of each new Constraint of this NuclearUnitBlock is not
  *      greater than the tolerance.
  *
  * Every Constraint of this NuclearUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or #f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  *   - If \p fsbc is not a nullptr and it is a pointer to a
  *     SimpleConfiguration< double >, then the tolerance is the value present
  *     in that SimpleConfiguration and the relative violation is considered.
  *
  *   - If \p fsbc is not nullptr and it is a
  *     SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *     fsbc->f_value.first and the type of violation is determined by
  *     fsbc->f_value.second (any nonzero number for relative violation and
  *     zero for absolute violation);
  *
  *   - Otherwise, if both #f_BlockConfig and
  *     f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *     latter is a pointer to either a SimpleConfiguration< double > or to a
  *     SimpleConfiguration< std::pair< double , int > >, then the values of the
  *     parameters are obtained analogously as above;
  *
  *   - Otherwise, by default, the tolerance is 0 and the relative violation
  *     is considered.
  *
  * This function currently considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * currently ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed
  * (Constraint::compute()).
  *
  * @param useabstract This parameter is currently ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *        and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE NuclearUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NuclearUnitBlock
 * @{ */

 /// returns the initial modulation value
 double get_initial_modulation( void ) const {
  return( f_initial_modulation );
  }

 /// returns the minimum allowed modulation interval
 Index get_modulation_interval( void ) const {
  return( f_modulation_interval );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of modulation delta ramp-up
 /** The returned vector contains the modulation delta ramp-up at each time.
  * The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_modulation_ramp_up( void ) const {
  return( v_modulation_ramp_up );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of modulation delta ramp-down
 /** The returned vector contains the modulation delta ramp-up at each time.
  * The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_modulation_ramp_down( void ) const {
  return( v_modulation_ramp_down );
  }

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Variable OF THE NuclearUnitBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the NuclearUnitBlock
 * @{ */

 /// returns the vector of modulation variables

 ColVariable * get_modulation( void ) {
  if( v_modulation.empty() )
   return( nullptr );
  return &( v_modulation.front() );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE NuclearUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NuclearUnitBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * NuclearUnitBlock (extending that of ThermalUnitBlock). See
 * NuclearUnitBlock::deserialize( netCDF::NcGroup ) (and therefore
 * ThermalUnitBlock::deserialize( netCDF::NcGroup )) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR INITIALIZING THE NuclearUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the NuclearUnitBlock
 *  @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "NuclearUnitBlock::load() not implemented yet" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
 /* Method for handling Modification.
  *
  * This method has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the NuclearUnitBlock, and
  * "translate" them into both changes of the actual data structures and
  * corresponding "physical Modification". These Modification are those
  * for which Modification::concerns_Block() is true.
  *
  * Not needed: NuclearUnitBlock does not allow changes via (its part of)
  * the "abstract representation", any Modification produced in that way
  * will pass through ThermalUnitBlock::add_Modification(), not be
  * recognised as valid, and throw exception. */

 // void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// update a subset of the modulation ramp up values of the unit
 /** This method updates the modulation ramp up values of the unit. The
  * \p subset parameter contains a list of time instants and \p values
  * contains the modulation ramp up values of the unit at those time
  * instants. The modulation ramp up values of the unit at time subset[ i ]
  * is given by std::next( values , i ) for each i in
  * { 0 , ..., subset.size() - 1 }.
  *
  * Recall that the modulation ramp up values must be non-negative and not
  * larger than the original ramp up values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_up( MF_dbl_it values ,
			      Subset && subset , bool ordered = false ,
			      ModParam issuePMod = eNoBlck ,
			      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update a range of the modulation ramp up values of the unit
 /** This method updates the modulation ramp up values of the unit. The
  * \p rng parameter contains a range of time instants and \p values contains
  * the modulation ramp up values of the unit at those time instants. The
  * modulation ramp up values of the unit at time rng.first + i is given by
  * std::next( values , i ) for each i in {0, ...,
  * ( std::min( rng.second, get_time_horizon() ) - rng.first - 1 )}.
  *
  * Recall that the modulation ramp up values must be non-negative and not
  * larger than the original ramp up values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_up( MF_dbl_it values , Range rng = INFRange ,
			      ModParam issuePMod = eNoBlck ,
			      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update a subset of the modulation ramp down values of the unit
 /** This method updates the modulation ramp down values of the unit. The
  * \p subset parameter contains a list of time instants and \p values
  * contains the modulation ramp down values of the unit at those time
  * instants. The modulation ramp down values of the unit at time subset[ i ]
  * is given by std::next( values , i ) for each i in
  * { 0 , ..., subset.size() - 1 }.
  *
  * Recall that the modulation ramp down values must be non-negative and not
  * larger than the original ramp down values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_down( MF_dbl_it values ,
				Subset && subset , bool ordered = false ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update a range of the modulation ramp down values of the unit
 /** This method updates the modulation ramp down values of the unit. The
  * \p rng parameter contains a range of time instants and \p values contains
  * the modulation ramp down values of the unit at those time instants. The
  * modulation ramp down values of the unit at time rng.first + i is given by
  * std::next( values , i ) for each i in {0, ...,
  * ( std::min( rng.second, get_time_horizon() ) - rng.first - 1 )}.
  *
  * Recall that the modulation ramp down values must be non-negative and not
  * larger than the original ramp down values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_down( MF_dbl_it values , Range rng = INFRange ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*------------------- PROTECTED METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 // void guts_of_add_Modification( p_Mod mod , ChnlName chnl );
 // not required, changes via the "abstract representation" are not supported

 void check_modulation_consistency( void ) const;

 void update_initial_power_in_cnstrs( c_ModParam issueAMod = eNoBlck )
  override;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

//-------------------------------- data -------------------------------------

 /// the vector of modulation delta ramp up
 std::vector< double > v_modulation_ramp_up;

 /// the vector of modulation delta ramp down
 std::vector< double > v_modulation_ramp_down;

 /// the initial modulation value
 int f_initial_modulation;

 /// the minimum allowed modulation interval
 int f_modulation_interval;

//----------------------------- Variable ------------------------------------

 /// the modulation (binary) variables
 std::vector< ColVariable > v_modulation;

//---------------------------- Constraint -----------------------------------

 /// the Modulation RampUp time constraints
 std::vector< FRowConstraint > Modulation_RampUp_Constraints;

 /// the Modulation RampDown time constraints
 std::vector< FRowConstraint > Modulation_RampDown_Constraints;

 /// the NoDownModulation constraints
 std::vector< FRowConstraint > NoDownModulation;

 /// the NoStartUpModulation constraints
 std::vector< FRowConstraint > NoStartUpModulation;

 /// the Modulation constraints constraints proper
 std::vector< FRowConstraint > ModulationConst;

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
 /// register the methods in the methods factory

 static void static_initialization( void ) {
  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args(),
   *
   * register_method< NuclearUnitBlock >(
   *           "NuclearUnitBlock::set_modulation_ramp_up" ,
   *           & NuclearUnitBlock::set_availability , MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< NuclearUnitBlock , MF_dbl_it , Subset && , bool >(
                               "NuclearUnitBlock::set_modulation_ramp_up" ,
                               & NuclearUnitBlock::set_modulation_ramp_up );

  register_method< NuclearUnitBlock , MF_dbl_it , Range >(
                               "NuclearUnitBlock::set_modulation_ramp_up" ,
                               & NuclearUnitBlock::set_modulation_ramp_up );

  register_method< NuclearUnitBlock , MF_dbl_it , Subset && , bool >(
                             "NuclearUnitBlock::set_modulation_ramp_down" ,
                             & NuclearUnitBlock::set_modulation_ramp_down );

  register_method< NuclearUnitBlock , MF_dbl_it , Range >(
                             "NuclearUnitBlock::set_modulation_ramp_down" ,
                             & NuclearUnitBlock::set_modulation_ramp_down );
  }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( class( NuclearUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS NuclearUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// Derived class from Modification for modifications to a NuclearUnitBlock
class NuclearUnitBlockMod : public ThermalUnitBlockMod {

 public:

 /// Public enum for the types of ThermalUnitBlockMod
 enum NUB_mod_type {
  eSetModDP = eTUBBModLastParam , ///< set modulation delta ramp up
  eSetModDM ,                     ///< set modulation delta ramp down
  eNUBModLastParam  ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of NuclearUnitBlockMod. */
  };

 /// Constructor, takes the NuclearUnitBlock and the type
 NuclearUnitBlockMod( NuclearUnitBlock * fblock , int type )
  : ThermalUnitBlockMod( fblock , type ) {}

 ///< Destructor, does nothing
 virtual ~NuclearUnitBlockMod() override = default;

 /// prints the NuclearUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "NuclearUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case eSetMaxP:
    output << "set max power values";
    break;
   case eSetInitP:
    output << "set initial power values";
    break;
   case eSetInitUD:
    output << "Set initial up/down times";
    break;
   case eSetAv:
    output << "Set availability";
    break;
   case eSetSUC:
    output << "Set startup costs";
    break;
   case eSetLinT:
    output << "Set linear term";
    break;
   case eSetQuadT:
    output << "Set quad term";
    break;
   case eSetConstT:
    output << "Set constant term";
    break;
   case eSetModDP:
    output << "set modulation delta ramp up";
    break;
   case eSetModDM:
    output << "set modulation delta ramp down";
    break;
   default:;
   }
  }
 }; // end( class( NuclearUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS NuclearUnitBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from NuclearUnitBlockMod for "ranged" changes in a nuclear

class NuclearUnitBlockRngdMod : public NuclearUnitBlockMod {

 public:

 /// constructor: takes the NuclearUnitBlock, the type, and the range
 NuclearUnitBlockRngdMod( NuclearUnitBlock * fblock , int type ,
                          Block::Range rng )
  : NuclearUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~NuclearUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the NuclearUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  NuclearUnitBlockRngdMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

 };  // end( class( ThermalUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS NuclearUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from NuclearUnitBlockMod for "subset" changes in a nuclear

class NuclearUnitBlockSbstMod : public NuclearUnitBlockMod {

 public:

 /// constructor: takes the NuclearUnitBlock, the type, and the subset
 NuclearUnitBlockSbstMod( NuclearUnitBlock * fblock , int type ,
                          Block::Subset && nms )
  : NuclearUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~NuclearUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the NuclearUnitBlockSbstMod
 void print( std::ostream &output ) const override {
  NuclearUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

 Block::Subset f_nms;  ///< the subset

 };  // end( class( NuclearUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* NuclearUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File NuclearUnitBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
