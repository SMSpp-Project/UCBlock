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
 *                   Rafael Durbano Lobato
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

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS UnitBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of the Block concept for "a generic unit" in UC
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
 *     (iv)  the active power produced by the generators in the unit.
 *
 * The class also outputs some general information regarding how the active
 * power and/or commitment status of each generator of the unit at a given
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
 *  @{ */

 /// Constructor, takes the father and the time horizon
 /** Constructor of UnitBlock, taking possibly a pointer of its father
  * Block and the time horizon. By default the time horizon is initialized to
  * 0, which means "not set yet". */

 explicit UnitBlock( Block * father_block = nullptr , Index t = 0 );

/*--------------------------------------------------------------------------*/
 /// Destructor of UnitBlock

 virtual ~UnitBlock() override {
  for( auto & block : v_Block )
   delete block;
  v_Block.clear();
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
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
 *   NumberIntervals = k of the form [ 0 , i_0 ], [ i_0 + 1 , i_1 ], ...  [
 *   i_{k-2} + 1 , "TimeHorizon" - 1 ]; "ChangeIntervals" then has to contain
 *   [ i_0 , i_1 , ... , i_{k-2} ] as the first k-1 elements. Note that, since
 *   the upper endpoint of the last interval must necessarily be "TimeHorizon"
 *   - 1, the last element of "ChangeIntervals", namely ChangeIntervals[
 *   NumberIntervals - 1 ], is ignored and does not need to be set (although
 *   the variable has actually "NumberIntervals" elements). Anyway, the whole
 *   variable is ignored if either "NumberIntervals" <= 1 (such as if it is
 *   not defined), or "NumberIntervals" >= "TimeHorizon". */
 void deserialize( const netCDF::NcGroup & group ) override;

/**@} ----------------------------------------------------------------------*/
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
 Index get_time_horizon() const { return f_time_horizon; }

/*--------------------------------------------------------------------------*/
 /// returns the number of electrical generators of each unit in the problem
 /** Returns the number of electrical generators for this UnitBlock. Since in
  *  most of the cases each unit has only one electrical generator, this
  *  method in the base UnitBlock class returns to one by default. Therefore,
  *  for all the units that have only one generator, the implementation of
  *  this method is already done right in the base UnitBlock class. The units
  *  that have more than one electrical generator (tied together by technical
  *  constraints) will have to handle this number by their-self. */

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
  *        desired. */

 virtual double * get_fixed_consumption( Index generator ) {
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
  *        desired. */

 virtual double * get_inertia_commitment( Index generator ) {
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
  *        desired. */

 virtual double * get_inertia_power( Index generator ) {
  return( nullptr );
 }

/**@} ----------------------------------------------------------------------*/
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
 * @{ */

 /// returns the array of commitment variables
 /** This method returns a pointer to the array containing the commitment
  * variables of the given \p generator at all time instants. Being C the
  * value returned by this method, C[t] is the commitment variable at time t
  * for each t in {0, ..., time_horizon - 1}.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the commitment variable (if any).
  *
  * @param generator The index of the generator whose commitment variables are
  *        desired. */

 virtual ColVariable * get_commitment( Index generator ) {
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the array of primary spinning reserve variables
 /** This method returns a pointer to the array containing the primary
  * spinning reserve variables of the given \p generator at all time
  * instants. Being R the value returned by this method, R[t] is the primary
  * spinning reserve variable at time t for each t in {0, ..., time_horizon -
  * 1}.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the primary spinning reserve variable (if
  * any).
  *
  * @param generator The index of the generator whose primary spinning reserve
  *        variables are desired. */

 virtual ColVariable * get_primary_spinning_reserve( Index generator ) {
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the array of secondary spinning reserve variables
 /** This method returns a pointer to the array containing the secondary
  * spinning reserve variables of the given \p generator at all time
  * instants. Being R the value returned by this method, R[t] is the secondary
  * spinning reserve variable at time t for each t in {0, ..., time_horizon -
  * 1}.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the secondary spinning reserve variable (if
  * any).
  *
  * @param generator The index of the generator whose secondary spinning
  *        reserve variables are desired. */

 virtual ColVariable * get_secondary_spinning_reserve( Index generator ) {
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the array of active power variables
 /** This method returns a pointer to the array containing the active power
  * variables of the given \p generator at all time instants. Being P the
  * value returned by this method, P[t] is the active power variable at time t
  * for each t in {0, ..., time_horizon - 1}.
  *
  * The default implementation of this method returns nullptr; and derived
  * classes will have to handle the active power variable (if any).
  *
  * @param generator The index of the generator whose active power variables
  *        are desired. */

 virtual ColVariable * get_active_power( Index generator ) {
  return( nullptr );
 }
/**@} ----------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */
 /// returns a Solution representing the current solution of this UnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this UnitBlock. The base
  * UnitBlock class defaults to ColVariableSolution, RowConstraintSolution,
  * and ColRowSolution, but :UnitBlock may make different choices.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value. If this value is
  *
  * - 1, then a RowConstraintSolution is returned;
  *
  * - 2, then a ColRowSolution is returned;
  *
  * - any other value, then a ColVariable Solution is returned.
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
  * - otherwise, it is 0. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the UnitBlock
 *  @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  *  UnitBlock. See UnitBlock::deserialize( netCDF::NcGroup ) for
  *  details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*---------------- METHODS FOR MODIFYING THE UnitBlock ---------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for modifying the UnitBlock
  *  @{ */

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
 /// sets reserve vars method
 /** This method can be called *after* that deserialize() and before
  * generate_abstract_variables() and generate_abstract_constraints(). This is
  * called to provide the UCBlock with the reserve variables if it's needed.
  * The input parameter is a bitwise value that allows to specify which unit
  * could have the reserve variables:
  *
  * - 1 the unit could have primary spinning reserve variables
  * - 2 the unit could have secondary spinning reserve variables
  * - 4 the unit could have inertia reserve variables.
  *
  * Note: this method is only to "destroy" the (primary, secondary and inertia)
  * reserve variables; it cannot create them if they are not there.*/

 virtual void set_reserve_vars( unsigned char what ) {
  reserve_vars = what;
 }

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
 /** @name Handling the data of the UnitBlock
    @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "UnitBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// deserializes the time horizon from a netCDF group
 void deserialize_time_horizon( const netCDF::NcGroup & group );

 /// deserializes the change intervals vector from a netCDF group
 void deserialize_change_intervals( const netCDF::NcGroup & group );

 /// states that the Variable of the UnitBlock have been generated
 void set_variables_generated() { AR |= HasVar; }

 /// states that the Constraint of the UnitBlock have been generated
 void set_constraints_generated() { AR |= HasCst; }

 /// states that the Objective of the UnitBlock has been generated
 void set_objective_generated() { AR |= HasObj; }

 /// indicates whether the Variable of the UnitBlock have been generated
 bool variables_generated() const { return( AR & HasVar ); }

 /// indicates whether the Constraint of the UnitBlock have been generated
 bool constraints_generated() const { return( AR & HasCst ); }

 /// indicates whether the Objective of the UnitBlock has been generated
 bool objective_generated() const { return( AR & HasObj ); }

 /// Resizes a vector to time_horizon by using change_intervals
 template< typename T >
 void decompress_vector( std::vector< T > & v );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// the time horizon of the problem
 Index f_time_horizon;

 /// the number of intervals
 Index f_number_intervals;

 /// the vector of change intervals
 std::vector< Index > v_change_intervals;

 unsigned char reserve_vars{};
 ///< bit-wise coded: which reserve variables generate

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 unsigned char AR{}; ///< bit-wise coded: what abstract is there

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variables have been constructed
 static constexpr unsigned char HasCst = 2;
 ///< second bit of AR == 1 if the Constraints have been constructed
 static constexpr unsigned char HasObj = 4;
 ///< third bit of AR == 1 if the Objective has been constructed

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( class( UnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* UnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File UnitBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
