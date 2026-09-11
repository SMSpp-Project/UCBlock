/*--------------------------------------------------------------------------*/
/*-------------------------- File UnitBlock.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class UnitBlock, which derives from the Block, in order
 * to define a base class for any possible "unit" in a UCBlock. A unit is in
 * general a set of electrical generators tied together by some technical
 * constraints, although many units actually correspond to only one generator.
 * The base UnitBlock class only has very basic information that can
 * characterize almost any different kind of unit, which includes the length
 * of the time horizon, the number of electrical generators in the unit (one
 * by default), the default implementation of four sets of Variables: active
 * power variables, commitment variables, primary and secondary spinning
 * reserve variables, defined for each time instant in the time horizon and
 * each generator in the unit. It also outputs some general information
 * regarding how the active power and/or commitment status of the generators
 * in the unit at a given time instant impact the unit's capability of
 * satisfying inertia constraints, and the fixed consumption (if any) of the
 * generators in the unit when they are off.
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __UnitBlock
 #define __UnitBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "ColVariable.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------------- FORWARD DECLARATIONS ---------------------------*/
/*--------------------------------------------------------------------------*/

 class UnitBlockSolution;  // forward definition of UnitBlockSolution

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS UnitBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for "a generic unit" in UC
/** The class UnitBlock, which derives from the Block, defines a base class
 * for any possible "unit" that can be attached to a UCBlock. A unit is in
 * general a set of electrical generators tied together by some technical
 * constraints, although many units actually correspond to only one generator.
 * The base UnitBlock class only has very basic information that can
 * characterize almost any different kind of unit:
 *
 * - The time horizon of the problem;
 *
 * - The number of generators in the unit (1 by default, see
 *   get_number_generators());
 *
 * - Whether or not the unit produces also reactive power in addition to
 *   active power.
 *
 * - The default implementation of four variables which are assumed that the
 *   Variable (of each type) for each generator are organised in arrays of
 *   size get_time_horizon() which are:
 *
 *     (i)   the commitment of the generators in the unit;
 *
 *     (ii)  the primary spinning reserve of the generators in the unit;
 *
 *     (iii) the secondary spinning reserve of the generators in the unit;
 *
 *     (iv)  the active power produced by the generators in the unit;
 *
 *     (v)   the reactive power produced by the generators in the unit.
 *
 * Note that not all units produce all types of variables, except for the
 * active power ones.
 *
 * The class also outputs some general information regarding how the active
 * power and / or commitment status of each generator of the unit at a given
 * time instant impact the unit's capability of satisfying inertia
 * constraints, and the fixed consumption (if any) of each generator in the
 * unit when it is off. */

class UnitBlock : public Block
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
 /** Constructor of UnitBlock, taking possibly a pointer of its father
  * Block and doing nothing. */

 explicit UnitBlock( Block * father = nullptr ) : Block( father ) ,
  f_time_horizon( 0 ) , f_number_intervals( 0 ) , reserve_vars( 0 ) ,
  f_reactive_power( false ) , AR( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of UnitBlock

 virtual ~UnitBlock() override {
  for( auto & block : v_Block )
   delete( block );
  v_Block.clear();
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the UnitBlock. Besides the mandatory "type" attribute of any :Block, the
  * group should contain the following:
  *
  * - The dimension "TimeHorizon" containing the time horizon. The dimension
  *   is optional because the same information may be passed via the method
  *   set_time_horizon(), or directly retrieved from the father if it is a
  *   UCBlock; see the comments to set_time_horizon() for details.
  *
  * - The dimension "NumberIntervals", that is provided to allow that all
  *   time-dependent data in the UnitBlock can only change at a subset of
  *   the time instants of the time interval, being therefore
  *   piecewise-constant (possibly, constant). "NumberIntervals" should
  *   therefore be <= "TimeHorizon", with four distinct cases:
  *
  *   i)   1 < "NumberIntervals" < "TimeHorizon", which means that at some
  *        time instants, *but not all of them*, the values of some of the
  *        relevant data are changing; the intervals are then described in
  *        variable "ChangeIntervals".
  *
  *   ii) "NumberIntervals" == 1, which means that the value of each relevant
  *        data in the UnitBlock is the same for each time instant 0, ...,
  *        "TimeHorizon" - 1 in the time horizon. In this case, the variable
  *        "ChangeIntervals" (see below) is ignored.
  *
  *   iii) "NumberIntervals" == "TimeHorizon", which means that values of the
  *        relevant data changes at every time interval (in principle; of
  *        course there is nothing preventing the same value to be repeated in
  *        the netCDF input). Also in this case the variable "ChangeIntervals"
  *        is ignored, since it is useless.
  *
  *   iv)  The dimension "NumberIntervals" is not provided, which means that
  *        the values of the relevant data may be the same for each time
  *        instant (as in case ii above) or indexed over "TimeHorizon" (as in
  *        case iii above). Also in this case, of course, "ChangeIntervals"
  *        (see below) is ignored, and therefore it can (and should) not be
  *        present.
  *
  *   Note that this (together with "ChangeIntervals", if defined) obviously
  *   sets the "maximum frequency" at which data can change; if some data
  *   changes less frequently (say, it is constant), then the same value
  *   will have to be repeated. Individual data can also have specific
  *   provisions for the case where the data is all equal despite
  *   "NumberIntervals" saying differently.
  *
  * - The variable "ChangeIntervals", of type integer and indexed over the
  *   dimension "NumberIntervals". The time horizon is subdivided into
  *   NumberIntervals = k of the form [ 0 , i_0 ], [ i_0 + 1 , i_1 ], ... [
  *   i_{k-2} + 1 , "TimeHorizon" - 1 ]; "ChangeIntervals" then has to contain
  *   [ i_0 , i_1 , ... , i_{k-2} ] as the first k-1 elements. Note that,
  *   since the upper endpoint of the last interval must necessarily be
  *   "TimeHorizon" - 1, the last element of "ChangeIntervals", namely
  *   ChangeIntervals[ NumberIntervals - 1 ], is ignored and does not need
  *   to be set (although the variable has actually "NumberIntervals"
  *   elements). Anyway, the whole variable is ignored if either
  *   "NumberIntervals" <= 1 (such as if it is not defined), or
  *   "NumberIntervals" >= "TimeHorizon". */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends Block::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends Block::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE UnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of electrical generation units, i.e.:
 *
 * - the time horizon
 *
 * - fixed consumption when the unit is off;
 *
 * - the contribution to the inertia depending on the commitment status;
 *
 * - the contribution to the inertia depending on the active power produced.
 *
 * Note that most of this data is optional, and only "few" of the UnitBlock
 * actually have it. Hence, no data structures are defined by the base
 * UnitBlock class: the methods have a default implementation returning an
 * empty vector, and derived classes will have to handle their own data (if
 * any).
 * @{ */

 /// returns the time horizon of the problem
 Index get_time_horizon( void ) const { return( f_time_horizon ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of electrical generators of each unit in the problem
 /** Returns the number of electrical generators for this UnitBlock. Since in
  * most of the cases, each unit has only one electrical generator, this
  * method in the base UnitBlock class returns to one by default. Therefore,
  * for all the units that have only one generator, the implementation of
  * this method is already done right in the base UnitBlock class. The units
  * that have more than one electrical generator (tied together by technical
  * constraints) will have to handle this number by their-self. */

 virtual Index get_number_generators( void ) const { return( 1 ); }

/*--------------------------------------------------------------------------*/
 /// returns the fixed consumption of the given generator
 /** This method returns a pointer to the array containing the fixed
  * consumption (basically, the constants to be multiplied by the commitment
  * variables returned by get_commitment()) of the given \p generator at all
  * time instants. Being C the value returned by this method, C[t] is the
  * inertia commitment of the given \p generator at the time instant t for
  * each t in {0, ..., time_horizon - 1}.
  *
  * The default implementation of the method returns nullptr; and derived
  * classes will have to handle their own data (if any).
  *
  * @param generator The index of the generator whose fixed consumption is
  *                  desired. */

 virtual const double * get_fixed_consumption( Index generator ) const {
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia commitment of the given generator
 /** This method returns a pointer to the array containing the contribution to
  * inertia (basically, the constants to be multiplied by the commitment
  * variables returned by get_commitment()) of the given \p generator at all
  * time instants. Being C the value returned by this method, C[t] is the
  * inertia commitment of the given \p generator at the time instant t for
  * each t in {0, ..., time_horizon - 1}.
  *
  * The default implementation of the method returns nullptr; and derived
  * classes will have to handle their own data (if any).
  *
  * @param generator The index of the generator whose inertia commitment is
  *                  desired. */

 virtual const double * get_inertia_commitment( Index generator ) const {
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia power of the given generator
 /** This method returns a pointer to the array of inertia power (basically,
  * the constants to be multiplied by the active power variables returned by
  * get_active_power()) of the given \p generator at all time instants. Being
  * P the value returned by this method, P[t] is the inertia power of the
  * given \p generator at the time instant t for each t in {0, ...,
  * time_horizon - 1}.
  *
  * The default implementation of the method returns nullptr; and derived
  * classes will have to handle their own data (if any).
  *
  * @param generator The index of the generator whose inertia power is
  *                  desired. */

 virtual const double * get_inertia_power( Index generator ) const {
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power of \p generator at time \p t

 virtual double get_min_power( Index t , Index generator = 0 ) const {
  return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power of \p generator at time \p t

 virtual double get_max_power( Index t , Index generator = 0 ) const {
  return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum reactive power of \p generator at time \p t

 virtual double get_min_reactive_power( Index t , Index generator = 0 )
  const {
  return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of \p generator at time \p t

 virtual double get_max_reactive_power( Index t , Index generator = 0 )
  const {
  return( 0 );
  }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE Variable OF THE UnitBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the UnitBlock
 *
 * These methods allow to read the four groups of Variable that any UnitBlock
 * in principle has (although some may not):
 *
 * - commitment variables;
 *
 * - primary spinning reserve variables;
 *
 * - secondary spinning reserve variables;
 *
 * - active power variables.
 *
 * - reactive power variables.
 * @{ */

 /// returns the array of commitment variables
 /** This method returns a pointer to the array containing the commitment
  * variables of the given \p generator at all time instants. Being C the
  * value returned by this method, C[ t ] is the commitment variable at time t
  * for each t in {0, ..., time_horizon - 1}.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the commitment variable (if any).
  *
  * @param generator The index of the generator whose commitment variables are
  *                  desired. */

 virtual ColVariable * get_commitment( Index generator ) {
  return( nullptr );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like get_commitment(), but returns a const * so that it can be const

 const ColVariable * get_const_commitment( Index generator ) const
 {
  // this is a dirty trick: casting away const-ness to be able to call the
  // standard version of the method; however, this allows to avoid to
  // redefine the *const_* version in derived classes, assuming of course
  // that their methods will do nothing except returning the pointer

  return( const_cast< UnitBlock * >( this )->get_commitment( generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the array of primary spinning reserve variables
 /** This method returns a pointer to the array containing the primary
  * spinning reserve variables of the given \p generator at all time
  * instants. Being R the value returned by this method, R[t] is the primary
  * spinning reserve variable at time t for each t in { 0 , ... ,
  * time_horizon - 1 }.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the primary spinning reserve variable (if
  * any).
  *
  * @param generator The index of the generator whose primary spinning reserve
  *                  variables are desired. */

 virtual ColVariable * get_primary_spinning_reserve( Index generator ) {
  return( nullptr );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like get_primary_spinning_reserve(), but returns a const *

 const ColVariable * get_const_primary_spinning_reserve( Index generator )
  const {
  // this is a dirty trick: casting away const-ness to be able to call the
  // standard version of the method; however, this allows to avoid to
  // redefine the *const_* version in derived classes, assuming of course
  // that their methods will do nothing except returning the pointer

  return( const_cast< UnitBlock * >( this )->get_primary_spinning_reserve(
							       generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the array of secondary spinning reserve variables
 /** This method returns a pointer to the array containing the secondary
  * spinning reserve variables of the given \p generator at all time
  * instants. Being R the value returned by this method, R[t] is the secondary
  * spinning reserve variable at time t for each t in { 0 , ... ,
  * time_horizon - 1 }.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the secondary spinning reserve variable (if
  * any).
  *
  * @param generator The index of the generator whose secondary spinning
  *                  reserve variables are desired. */

 virtual ColVariable * get_secondary_spinning_reserve( Index generator ) {
  return( nullptr );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like get_secondary_spinning_reserve(), but returns a const *

 const ColVariable * get_const_secondary_spinning_reserve( Index generator )
  const {
  // this is a dirty trick: casting away const-ness to be able to call the
  // standard version of the method; however, this allows to avoid to
  // redefine the *const_* version in derived classes, assuming of course
  // that their methods will do nothing except returning the pointer

  return( const_cast< UnitBlock * >( this )->get_secondary_spinning_reserve(
							       generator ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the array of active power variables
 /** This method returns a pointer to the array containing the active power
  * variables of the given \p generator at all time instants. Being P the
  * value returned by this method, P[t] is the active power variable at time t
  * for each t in { 0 , ... , time_horizon - 1 }.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the active power variable (if any).
  *
  * @param generator The index of the generator whose active power variables
  *                  are desired. */

 virtual ColVariable * get_active_power( Index generator ) {
  return( nullptr );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like get_active_power(), but returns a const * so that it can be const

 const ColVariable * get_const_active_power( Index generator ) const {
  // this is a dirty trick: casting away const-ness to be able to call the
  // standard version of the method; however, this allows to avoid to
  // redefine the *const_* version in derived classes, assuming of course
  // that their methods will do nothing except returning the pointer

  return( const_cast< UnitBlock * >( this )->get_active_power( generator ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the units of this UnitBlock are committed
 /** Returns true if this UnitBlock has commitment Variable, i.e., if the
  * units it describes are switched on and off rather than being always on.
  *
  * This is a property of the *data*, not of the abstract representation:
  * unlike get_commitment(), which is a Variable and therefore only exists
  * once the latter has been generated, this can be asked at any time, which
  * is what whoever has to decide what a Solution of this UnitBlock is made
  * of needs [see get_Solution()]. The base class has no commitment, a
  * :UnitBlock whose units are committed says so. */

 virtual bool has_commitment( void ) const { return( false ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if this UnitBlock has primary spinning reserve Variable
 /** Returns true if this UnitBlock has primary spinning reserve Variable,
  * i.e., if the enclosing UCBlock asked for them [see set_reserve_vars()]
  * and the units can provide them: as has_commitment(), this is a property
  * of the data and can be asked before the abstract representation is
  * generated. A :UnitBlock that cannot provide the reserve, whatever it is
  * asked, says so by redefining this. */

 virtual bool has_primary_reserve( void ) const {
  return( reserve_vars & 1u );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if this UnitBlock has secondary spinning reserve Variable

 virtual bool has_secondary_reserve( void ) const {
  return( reserve_vars & 2u );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if this UnitBlock has reactive power Variable
 /** Returns true if this UnitBlock has reactive power Variable, i.e., if the
  * enclosing UCBlock asked for them [see set_reactive_power()] and the units
  * can produce them; as has_commitment(), it is a property of the data. A
  * :UnitBlock that produces no reactive power, whatever it is asked, says so
  * by redefining this. */

 virtual bool has_reactive_power( void ) const {
  return( f_reactive_power );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of reactive power variables
 /** The base UnitBlock class does not handle reactive power variables, so
  *  this method always returns nullptr. */

 virtual ColVariable * get_reactive_power( Index generator ) {
  return( nullptr );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like get_reactive_power(), but returns a const * so that it can be const

 const ColVariable * get_const_reactive_power( Index generator ) const { 
  // this is a dirty trick: casting away const-ness to be able to call the
  // standard version of the method; however, this allows to avoid to
  // redefine the *const_* version in derived classes, assuming of course
  // that their methods will do nothing except returning the pointer

  return( const_cast< UnitBlock * >( this )->get_reactive_power( generator )
	  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the scale factor of this UnitBlock
 /** This method returns the scale factor of this UnitBlock. Since not every
  * UnitBlock may support the notion of scaling, this method has a default
  * implementation that returns 1. Derived classes that support scaling must
  * override this method. See UnitBlock::scale() for more details about the
  * scaling of a UnitBlock.
  *
  * @return The scaling factor number of this UnitBlock. */

 virtual double get_scale( void ) const { return( 1 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the upper bound of the design variable
 /** This method returns the upper bound of the design variable of this
  * UnitBlock when the unit is in design mode. The intended use is for
  * pre-computing fleet-level bounds on node injection in the embedding
  * UCBlock: the per-time-step contribution of one unit to the node
  * injection is bounded by `get_scale() * get_design_ub() *
  * get_max_power( t , g )`.
  *
  * UnitBlock provides a default implementation returning 1 (which is
  * correct for blocks whose design variable is binary, like
  * ThermalUnitBlock). UnitBlocks that support a design upper bound
  * greater than 1 (IntermittentUnitBlock / BatteryUnitBlock via
  * MaxCapacityDesign) override it.
  *
  * @return The upper bound of the design variable (1 if no design is
  *         in use or the design is binary). */

 virtual double get_design_ub( void ) const { return( 1 ); }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution representing the current solution of this UnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this UnitBlock. This may
  * either be a UnitBlockSolution or a further derived class containing more
  * specific solution information for derived :UnitBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - bit 0 (& 1) means "store the power", which surely is the active one and
  *         also the reactive one if any of the generators of the unit
  *         produces it (the others will produce 0 reactive power)
  *
  * - bit 1 (& 2) means "store the commitment"
  *
  * - bit 2 (& 4) means "store the primary reserve"
  *
  * - bit 3 (& 8) means "store the secondary reserve"
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 15 (save everything).
  *
  * Note, however, that the UnitBlockSolution may not contain some of the
  * data that ws implies since it may just not be present in the :UnitBlock.
  * This is true for commitment, primary and secondary reserve that are
  * optional, but not for active power which is mandatory. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" UnitBlockSolution
 /** Small virtual method that just returns an "empty" UnitBlockSolution
  * object. It is used by get_Solution(), with the idea that derived classes
  * can override it to make it return a :UnitBlockSolution better suited for
  * the specific :UnitBlock at hand. */

 virtual UnitBlockSolution * new_Solution( void ) const;

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the UnitBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * UnitBlock. See deserialize( const netCDF::NcGroup & ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*---------------- METHODS FOR MODIFYING THE UnitBlock ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the UnitBlock
 * @{ */

 /// sets the time horizon method
 /** This method can be called *before* that deserialize() is called to
  * provide the UnitBlock with the time horizon. This allows the information
  * not to be duplicated in the netCDF group that describes the unit, since
  * usually (but not necessarily) a UnitBlock is deserialized inside a
  * UCBlock, and all units have the same time horizon, that can therefore be
  * read once and for all by the father UCBlock.
  *
  * If this method is *not* called, which means that f_time_horizon is at its
  * initial value of 0 (not initialized), then when deserialize() is called
  * the information has to be available by other means, i.e.:
  *
  * (i)  If there is no dimension TimeHorizon in netCDF input, then the
  *      UnitBlock must have a father, which must be a UCBlock: the
  *      time horizon is then taken to be that of the father. If the
  *      UnitBlock does not have a father (or it is not a UCBlock), then
  *      exception is thrown.
  *
  * (ii) If the dimension TimeHorizon is present in the netCDF input of
  *      UnitBlock, the value provided there is used with no check that
  *      the UnitBlock has a father at all, that the father is a UCBlock,
  *      or that the two time horizon agree.
  *
  * If this method *is* called, which has to happen before that deserialize()
  * is called, and f_time_horizon is set at a value != 0, then if the
  * dimension TimeHorizon is present in netCDF input, then the two values must
  * agree. If the dimension TimeHorizon is not present, then the value set by
  * this method is used. Note that, of course, the data in the netCDF file (if
  * the unit has any data induced over the time horizon) has to agree with the
  * value set by this method.
  *
  * If this method is called *after* that deserialize() is called, this is
  * taken to mean that the UnitBlock is being "reset", and that
  * immediately after deserialize() will be called again. The same rules as
  * above are to be followed for that subsequent call to deserialize(). */

 void set_time_horizon( Index t ) { f_time_horizon = t; }

/*--------------------------------------------------------------------------*/
 /// sets which reserve variables should be there
 /** This method can be called *after* that deserialize() and before
  * generate_abstract_variables() and generate_abstract_constraints(). This
  * is called to provide the UnitBlock with information about which of the
  * possible reserve variables the enclosing UCBlock (if any) requires, and
  * therfore should actually be created. This is because the reason for
  * these variables existing is that there be constraints about the total
  * amount of available reserve (of some type), which is something that the
  * enclosing UCBlock (if any) knows about but the UnitBlock cannot possibly
  * know unless it is informed by calling this method.
  *
  * The \p what parameter is a bitwise value that allows to specify which
  * reserve variables (and, of course, the corresponding constraints) the
  * UnitBlock should define:
  *
  * - bit 0 (& 1): define primary spinning reserve variables
  *
  * - bit 1 (& 2): define secondary spinning reserve variables
  *
  * - bit 2 (& 4): define inertia reserve variables
  *
  * If the method is not called, 0 (no reserve variables) is assumed. Note
  * that calling this method will not guarantee that the variables will be
  * there, since some units do not provide reserve of some type; hence, the
  * method says "do define the variables if you can". */

 virtual void set_reserve_vars( unsigned char what = 0 ) {
  reserve_vars = what;
  }

/*--------------------------------------------------------------------------*/
 /// sets whether or not reactive power variables should be there
 /** This method can be called *after* that deserialize() and before
  * generate_abstract_variables() and generate_abstract_constraints(). This
  * is called to provide the UnitBlock with information about whether or not
  * the enclosing UCBlock (if any) requires reactive power variables to be
  * defined. This is because the reason for these variables existing is that
  * the interconnection network has a treatment for both active and reactive
  * power, which is something that the enclosing UCBlock (if any) knows about
  * but the UnitBlock cannot possibly know unless it is informed by calling
  * this method.
  *
  * If the method is not called, false (no reactive power variables) is
  * assumed. Note that calling this method will not guarantee that the
  * variables will be there, since some units do not provide reactive power;
  * hence, the method says "do define the variables if you can". */

 virtual void set_reactive_power( bool reactive = false ) {
  f_reactive_power = reactive;
  }

/*--------------------------------------------------------------------------*/
 /// sets the scale factor of this UnitBlock
 /** Some situations may require the presence of multiple identical units. By
  * identical units we mean units that represent the exactly same mathematical
  * model: they have the same type, data, variables, constraints, etc. In
  * short, these are units that are equivalent in every aspect. An immediate
  * way of considering multiple identical units could be simply to have
  * multiple instances of the same unit. In particular cases, however, it may
  * be possible to have a compact and computational efficient representation
  * of this set of identical units. One of these cases occurs when all units
  * behave exactly as each other, as if they were synchronized and performing
  * the same tasks simultaneously. This case is implemented by stating that a
  * UnitBlock can be scaled.
  *
  * A scaled UnitBlock must be interpreted as follows. The four sets of
  * Variable considered by the UnitBlock (namely, active power, commitment,
  * and primary and secondary spinning reserves) represent what happens to the
  * UnitBlock independently of the scale factor. For instance, consider the
  * active power variables. These variables represent the active power
  * produced by the generators of the UnitBlock. If we denote by \f$ P(g,t)
  * \f$ the value of the Variable representing the active power produced by
  * generator g at time t (see get_active_power()) and by \f$ S \f$ the scale
  * factor of this UnitBlock, then the g-th generator of this UnitBlock
  * <b>must be interpreted</b> as if it would produce \f$ S P(g,t) \f$ at time
  * t. That is, the scale factor does not affect the values of the Variable of
  * the UnitBlock. Any other object that uses the values of these active power
  * variables must explicitly multiply them by the scale factor in order to
  * obtain the actual amount of active power produced by the unit. The
  * treatment of the primary and secondary spinning reserves variables is
  * similar. Notice, however, that the values of the commitment variables
  * should not be multiplied by the scale factor, as they simply indicate
  * whether each generator of the unit is committed or not.
  *
  * The Objective of the UnitBlock, on the other hand, has a different
  * relation with the scale factor than that of Variable. If the UnitBlock has
  * an Objective, then this Objective must take into account the scale factor.
  * That is, differently from what happens to the Variable of the UnitBlock,
  * no action is required on the part of an object that uses the Objective of
  * this UnitBlock: the value of the Objective already considers the scale
  * factor. For instance, suppose that the Objective of a UnitBlock represents
  * the cost of that unit, say, \f$ \sum_{g,t} C(g,t) P(g,t) \f$, where \f$
  * C(g,t) \f$ is the cost of generating one unit of power by the g-th
  * generator at time t. Then, being \f$ S \f$ the scale factor of the
  * UnitBlock, its Objective must actually be \f$ S \sum_{g,t} C(g,t) P(g,t)
  * \f$.
  *
  * Since not every UnitBlock may support the notion of scaling, this method
  * has a default empty implementation. Derived classes that support scaling
  * must override this method.
  *
  * If the UnitBlock is really modified, a UnitBlockMod with type
  * UnitBlockMod::eScale is issued depending on the value of \p issuePMod.
  *
  * @param values An iterator to a vector containing the scale factor.
  *
  * @param subset If non-empty, the scale factor must be set to the value
  *               pointed by \p values. If empty, no operation must be
  *               performed.
  *
  * @param ordered This parameter is ignored.
  *
  * @param issuePMod Controls how physical Modification are issued.
  *
  * @param issueAMod Controls how abstract Modification are issued. */

 virtual void scale( MF_dbl_it values ,
                     Subset && subset ,
                     const bool ordered = false ,
                     c_ModParam issuePMod = eNoBlck ,
                     c_ModParam issueAMod = eNoBlck ) { }

/*--------------------------------------------------------------------------*/
 /// sets the scale factor of this UnitBlock
 /** This method sets the scale factor of this UnitBlock. A default
  * implementation is provided which simply call the Subset version of this
  * method. See UnitBlock::scale() for the semantics of scaling a
  * UnitBlock. If the UnitBlock is really modified, a UnitBlockMod with type
  * UnitBlockMod::eScale is issued depending on the value of \p issuePMod.
  *
  * @param values An iterator to a vector containing the scale factor.
  *
  * @param rng If non-empty, the scale factor is set to the value pointed by
  *            \p values. If empty, no operation is performed.
  *
  * @param issuePMod Controls how physical Modification are issued.
  *
  * @param issueAMod Controls how abstract Modification are issued. */

 virtual void scale( MF_dbl_it values ,
                     Range rng = Range( 0, Inf< Index >() ) ,
                     c_ModParam issuePMod = eNoBlck ,
                     c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the scale factor of this UnitBlock
 /** This method sets the scale factor of this UnitBlock. A default
  * implementation is provided which simply call the Subset version of this
  * method. See UnitBlock::scale() for the semantics of scaling a
  * UnitBlock. If the UnitBlock is really modified, a UnitBlockMod with type
  * UnitBlockMod::eScale is issued depending on the value of \p issuePMod.
  *
  * @param scale_factor The factor by which this UnitBlock should be scaled.
  *
  * @param issuePMod Controls how physical Modification are issued.
  *
  * @param issueAMod Controls how abstract Modification are issued. */

 virtual void scale( double scale_factor , c_ModParam issuePMod = eNoBlck ,
                     c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the linearization coefficient of the kappa-parametrized objective
 /** When the investment in this UnitBlock is represented by its kappa constant
  * (rather than by its scale factor), the kappa multiplies the right-hand side
  * of a set of the UnitBlock's own Constraints. By the envelope theorem, the
  * derivative of the optimal objective value of this UnitBlock with respect to
  * kappa is the sum, over those Constraints, of the dual value times the base
  * coefficient that kappa multiplies. This method computes that derivative
  * from the current (dual) solution of this UnitBlock, so that it can be used
  * as the linearization coefficient of the value function with respect to the
  * Variable associated with the kappa of this UnitBlock.
  *
  * Since not every UnitBlock supports the notion of kappa, this method has no
  * meaningful default: UnitBlocks whose investment is represented by kappa
  * (IntermittentUnitBlock, BatteryUnitBlock) override it; the default throws.
  *
  * @return The derivative of the optimal objective value of this UnitBlock
  *         with respect to its kappa constant. */

 virtual double get_kappa_linearization( void ) const {
  throw( std::logic_error( "UnitBlock::get_kappa_linearization: kappa is not "
                           "supported by this UnitBlock" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "UnitBlock::load() not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED TYPES OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 using MAdouble = boost::multi_array< double , 2 >;
 using MAdouble_ext = MAdouble::extent_gen;
 
 using MACV = boost::multi_array< ColVariable , 2 >;
 using MACV_ext = MAdouble::extent_gen;

 using MABC = boost::multi_array< BoxConstraint , 2 >;
 using MABC_ext = MAdouble::extent_gen;

 using MAFRC = boost::multi_array< FRowConstraint , 2 >;
 using MAFRC_ext = MAdouble::extent_gen;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// deserializes the time horizon from a netCDF group
 void deserialize_time_horizon( const netCDF::NcGroup & group );

 /// deserializes the change intervals vector from a netCDF group
 void deserialize_change_intervals( const netCDF::NcGroup & group );

 /// states that the Variable of the UnitBlock have been generated
 void set_variables_generated( void ) { AR |= HasVar; }

 /// states that the Constraint of the UnitBlock have been generated
 void set_constraints_generated( void ) { AR |= HasCst; }

 /// states that the Objective of the UnitBlock has been generated
 void set_objective_generated( void ) { AR |= HasObj; }

 /// indicates whether the Variable of the UnitBlock have been generated
 bool variables_generated( void ) const { return( AR & HasVar ); }

 /// indicates whether the Constraint of the UnitBlock have been generated
 bool constraints_generated( void ) const { return( AR & HasCst ); }

 /// indicates whether the Objective of the UnitBlock has been generated
 bool objective_generated( void ) const { return( AR & HasObj ); }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// the time horizon of the problem
 Index f_time_horizon;

 /// the number of intervals
 Index f_number_intervals;

 /// the vector of change intervals
 std::vector< std::size_t > v_change_intervals;

 /// bit-wise coded: which reserve variables generate
 unsigned char reserve_vars;

 /// bool: whether or not define reactive power variables
 bool f_reactive_power;

 ///< bit-wise coded: what abstract is there
 unsigned char AR;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 static constexpr unsigned char HasVar = 16;
 ///< 5th bit of AR == 1 if the Variables have been constructed

 static constexpr unsigned char HasCst = 32;
 ///< 6th bit of AR == 1 if the Constraints have been constructed

 static constexpr unsigned char HasObj = 64;
 ///< 7th bit of AR == 1 if the Objective has been constructed

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/



};  // end( class( UnitBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS UnitBlockMod ----------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a UnitBlock
class UnitBlockMod : public Modification {

public:

 /// public enum for the types of UnitBlockMod
 enum UB_mod_type {
  eScale = 0 ,    ///< set the scale factor
                  /**< This indicates that the scale factor of the UnitBlock
                   * has been modified. See UnitBlock::scale(). */
  eUBModLastParam ///< first allowed parameter value for derived classes
                  /**< Convenience value to easily allow derived classes to
                   * extend the set of types of UnitBlockMod. */
 };

 /// constructor, takes the UnitBlock and the type
 UnitBlockMod( UnitBlock * const fblock, const int type )
  : f_Block( fblock ), f_type( type ) {}

 /// destructor, default version
 virtual ~UnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 /// accessor to the type of modification
 int type( void ) const { return( f_type ); }

 protected:

 /// prints the UnitBlockMod
 void print( std::ostream & output ) const override {
  output << "UnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eScale ):
    output << "Set the scale factor";
    break;
   default:;
  }
 }

 /// pointer to the Block to which the Modification refers
 UnitBlock * f_Block{};

 int f_type;  ///< type of modification

};  // end( class( UnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS UnitBlockSolution --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a Solution of a UnitBlock
/** The UnitBlockSolution class, derived from Solution, represents a solution
 * of a "generic" UnitBlock, i.e., the values of
 *
 * - active [and possibly reactive] power variables;
 *
 * - [possibly] commitment variables;
 *
 * - [possibly] primary spinning reserve variables;
 *
 * - [possibly] secondary spinning reserve variables
 *
 * for every generator of the unit. UnitBlockSolution is not thought to be
 * "final", since :UnitBlock may want to define and handle their derived
 * :UnitBlockSolution to store unit-specific solution information. */

class UnitBlockSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 using Index = Block::Index;  // "import" Index

/*------------------------------- FRIENDS ----------------------------------*/

 friend UnitBlock;  ///< make UnitBlock friend

/*------------- CONSTRUCTING AND DESTRUCTING UnitBlockSolution -------------*/

 explicit UnitBlockSolution( void ) : f_time_horizon( 0 ) ,
  f_number_generators( 0 ) {}  /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~UnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*---------- METHODS DESCRIBING THE BEHAVIOR OF A UnitBlockSolution --------*/

 void read( const Block * block ) override;

 void write( Block * block ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a UnitBlockSolution into a netCDF::NcGroup
 /** Serialize a UnitBlockSolution into a netCDF::NcGroup, with the following
  * format:
  *
  * - The dimension "TimeHorizon" containing the number of time steps in the
  *   problem. It is mandatory.
  *
  * - The dimension "NumberGenerators" containing the number of generators
  *   in the unit; the dimension is optional, if it is missing then 1 (one)
  *   generator is assumed.
  *
  * - The variable "ActivePower", of type netCDF::NcDouble. If
  *   "NumberGenerators" is defined then it is indexed both over the
  *   dimensions "NumberGenerators" and "TimeHorizon", otherwise only
  *   over the dimension "TimeHorizon". ActivePower[ i , t ] is assumed to
  *   contain the optimal active power for generator i at the time t. The
  *   variable is optional.
  *
  * - The variable "ReactivePower", of type netCDF::NcDouble. If
  *   "NumberGenerators" is defined then it is indexed both over the
  *   dimensions "NumberGenerators" and "TimeHorizon", otherwise only
  *   over the dimension "TimeHorizon". ReactivePower[ i , t ] is assumed to
  *   contain the optimal active power for generator i at the time t. The
  *   variable is optional: it is in principle there if "ActivePower" is
  *   defined, but only if at least one of the generators of the unit
  *   actually produces it (if only a subset of the generators do, the
  *   variable is still defined for all of ones and will be filled with
  *   zeros for those who do not).
  *
  * - The variable "Commitment", of type netCDF::NcDouble (note that
  *   commitment variables are generally integer valued, in fact binary,
  *   but one may want to save the values of continuous relaxations).
  *   If "NumberGenerators" is defined then it is indexed both over the
  *   dimensions "NumberGenerators" and "TimeHorizon", otherwise only
  *   over the dimension "TimeHorizon". Commitment[ i , t ] is assumed to
  *   contain the optimal active power for generator i at the time t. The
  *   variable is optional.
  *
  * - The variable "PrimaryReserve", of type netCDF::NcDouble. If
  *   "NumberGenerators" is defined then it is indexed both over the
  *   dimensions "NumberGenerators" and "TimeHorizon", otherwise only
  *   over the dimension "TimeHorizon". PrimaryReserve[ i , t ] is assumed
  *   to contain the optimal active power for generator i at the time t. The
  *   variable is optional.
  *
  * - The variable "SecondaryReserve", of type netCDF::NcDouble. If
  *   "NumberGenerators" is defined then it is indexed both over the
  *   dimensions "NumberGenerators" and "TimeHorizon", otherwise only
  *   over the dimension "TimeHorizon". SecondaryReserve[ i , t ] is assumed
  *   to contain the optimal active power for generator i at the time t. The
  *   variable is optional. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 UnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

/*----------- METHODS FOR READING AND WRITING THE SOLUTION -----------------*/
/** @name Reading and writing the solution
 *
 * The parts of the solution this UnitBlockSolution saves, each of which is
 * empty if it does not save it [see UnitBlock::get_Solution()]. The setters
 * are what a Solver fills the Solution with directly out of its own data
 * structures, rather than writing the solution in the Variable of the
 * UnitBlock and having it read back from there, which requires the Variable
 * to exist at all.
 *  @{ */

 /// returns the active power saved in this UnitBlockSolution

 [[nodiscard]] const boost::multi_array< double , 2 > & get_active_power(
  void ) const { return( v_active_power ); }

 /// returns the reactive power saved in this UnitBlockSolution

 [[nodiscard]] const boost::multi_array< double , 2 > & get_reactive_power(
  void ) const { return( v_reactive_power ); }

 /// returns the commitment saved in this UnitBlockSolution

 [[nodiscard]] const boost::multi_array< double , 2 > & get_commitment(
  void ) const { return( v_commitment ); }

 /// returns the primary spinning reserve saved in this UnitBlockSolution

 [[nodiscard]] const boost::multi_array< double , 2 > &
  get_primary_spinning_reserve( void ) const {
  return( v_primary_reserve );
  }

 /// returns the secondary spinning reserve saved in this UnitBlockSolution

 [[nodiscard]] const boost::multi_array< double , 2 > &
  get_secondary_spinning_reserve( void ) const {
  return( v_secondary_reserve );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the number of generators and the time horizon of the solution
 /** Sets the size of the solution, which the setters below must agree with;
  * it is what read() takes from the UnitBlock. */

 void set_dimensions( Index number_generators , Index time_horizon ) {
  f_number_generators = number_generators;
  f_time_horizon = time_horizon;
  }

 /// sets the active power saved in this UnitBlockSolution

 void set_active_power( boost::multi_array< double , 2 > && ap ) {
  v_active_power.resize( boost::extents[ ap.shape()[ 0 ] ]
                                  [ ap.shape()[ 1 ] ] );
  v_active_power = ap;
  }

 /// sets the reactive power saved in this UnitBlockSolution

 void set_reactive_power( boost::multi_array< double , 2 > && rp ) {
  v_reactive_power.resize( boost::extents[ rp.shape()[ 0 ] ]
                                  [ rp.shape()[ 1 ] ] );
  v_reactive_power = rp;
  }

 /// sets the commitment saved in this UnitBlockSolution

 void set_commitment( boost::multi_array< double , 2 > && cm ) {
  v_commitment.resize( boost::extents[ cm.shape()[ 0 ] ]
                                  [ cm.shape()[ 1 ] ] );
  v_commitment = cm;
  }

 /// sets the primary spinning reserve saved in this UnitBlockSolution

 void set_primary_spinning_reserve( boost::multi_array< double , 2 > && pr ) {
  v_primary_reserve.resize( boost::extents[ pr.shape()[ 0 ] ]
                                  [ pr.shape()[ 1 ] ] );
  v_primary_reserve = pr;
  }

 /// sets the secondary spinning reserve saved in this UnitBlockSolution

 void set_secondary_spinning_reserve( boost::multi_array< double , 2 > && sr )
 {
  v_secondary_reserve.resize( boost::extents[ sr.shape()[ 0 ] ]
                                  [ sr.shape()[ 1 ] ] );
  v_secondary_reserve = sr;
  }

/** @} ---------------------------------------------------------------------*/

 UnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "UnitBlockSolution [" << this << "]: " << std::endl;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// do the heavy lifting of cloning a non-empty UnitBlockSolution
 /** This method does the actually copying of the fields for an already
  * existing :UnitBlockSolution; this is provided to make life easier to
  * the clone() of derived classes. */

 void guts_of_clone( UnitBlockSolution * sol ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// do the heavy lifting of scaling a non-empty UnitBlockSolution
 /** This method does the actually scaling of the fields for an already
  * existing :UnitBlockSolution; this is provided to make life easier to
  * the scale() of derived classes. */

 void guts_of_scale( UnitBlockSolution * sol , double factor ) const;

/*-------------------------- PROTECTED FIELDS ------------------------------*/

 Index f_time_horizon;            ///< the time horizon
 Index f_number_generators;       ///< the number of generators

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 boost::multi_array< double , 2 > v_active_power;
 ///< v_active_power[ i ][ t ] = active power of generator i at time t

 boost::multi_array< double , 2 > v_reactive_power;
 ///< v_reactive_power[ i ][ t ] = reactive power of generator i at time t

 boost::multi_array< double , 2 > v_commitment;
 ///< v_commitment[ i ][ t ] = commitment of generator i at time t

 boost::multi_array< double , 2 > v_primary_reserve;
 ///< v_primary_reserve[ i ][ t ] = primary reserve of generator i at time t

 boost::multi_array< double , 2 > v_secondary_reserve;
 ///< v_secondary_reserve[ i ][ t ] = secondary reserve of generator i at time t

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( UnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __UnitBlock */

/*--------------------------------------------------------------------------*/
/*------------------------ End File UnitBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
