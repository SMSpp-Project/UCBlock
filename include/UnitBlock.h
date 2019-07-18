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
 * by default), four sets of Variables: active power variables, commitment
 * variables, primary and secondary spinning reserve variables, defined for
 * each time instant in the time horizon and each generator in the unit. It
 * also outputs some general information regarding how the active power and/or
 * commitment status of the generators in the unit at a given time instant
 * impact the unit's capability of satisfying inertia constraints, and the
 * fixed consumption (if any) of the generators in the unit when they are off.
 *
 * \version 0.11
 *
 * \date 08 - 07 - 2019
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __UnitBlock
#define __UnitBlock   /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

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
 * - Four boost::multi_array< ColVariable, 2 > objects, that are used to
 *   store the information regarding:
 *
 *     (i)   the commitment of the generators in the unit;
 *
 *     (ii)  the primary spinning reserve of the generators in the unit;
 *
 *     (iii) the secondary spinning reserve of the generators in the unit;
 *
 *     (iv)  the active power produced by the generators in the unit.
 *
 *   Each of the boost::multi_array< ColVariable, 2 > has as first dimension
 *   the time horizon and as second dimension the number of generators (as
 *   returned get_number_generators()). There are two possible cases:
 *
 *   - if some boost::multi_array<ColVariable, 2> is empty(), then the
 *     corresponding variable does not exist (for instance, the generator
 *     in the unit may not have reserve);
 *
 *   - otherwise, the boost::multi_array< ColVariable, 2 >, say M, must have
 *     f_time_horizon rows and get_number_generators() columns: each element
 *     M[ t , g ] gives the variable for time step t of generator g.
 *
 * The class also outputs some general information regarding how the active
 * power and/or commitment status of each generator of the unit at a given
 * time instant impact the unit's capability of satisfying inertia
 * constraints, and the fixed consumption (if any) of each generator in the
 * unit when it is off. */

class UnitBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * UnitBlock defines the following main public type:
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
 /** Constructor of UnitBlock, taking possibly a pointer of its father
  * Block and the time horizon. By default the time horizon is initialized to
  * 0, which means "not set yet". */

 explicit UnitBlock( Block * father_block = nullptr, Index t = 0 );

/*--------------------------------------------------------------------------*/
 /// Destructor of UnitBlock: it is virtual, and empty

 ~UnitBlock() override = default;

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
 *   therefore be <= "TimeHorizon", with three distinct cases:
 *
 *    i)  "NumberIntervals" <= 1, which is taken to mean "NumberIntervals"
 *        == 1; this is what is assumed if the dimension, that is optional,
 *        is not there. This means that the value of each relevant data in
 *        the UnitBlock (see e.g. "FixedConsumption", "InertiaCommitment"
 *        and "InertiaPower" below) is the same for each time instant
 *        0, ..., "TimeHorizon" - 1 in the time horizon. In this case, the
 *        variable "ChangeIntervals" (see below) is ignored.
 *
 *   ii)  1 < "NumberIntervals" < "TimeHorizon", which means that in some
 *        time instants, *but not all of them*, the values of some of the
 *        relevant data are changing; the intervals are then described in
 *        variable "ChangeIntervals".
 *
 *   iii) "NumberIntervals" == "TimeHorizon",  which means that values of
 *        the relevant data changes at every time interval (in principle;
 *	       of course there is nothing preventing the same value to be
 *        repeated in the netCDF input). Also in this case the variable
 *        "ChangeIntervals" is ignored, since it is useless.
 *
 *   Note that this (together with "ChangeIntervals", if defined) obviously
 *   sets the "maximum frequency" at which data can change; if some data
 *   changes less frequently (say, it is constant), then the same value
 *   will have to be repeated. Individual data can also have specific
 *   provisions for the case where the data is all equal despite
 *   "NumberIntervals" saying differently, see e.g. "FixedConsumption".
 *
 * - The variable "ChangeIntervals", of type integer and indexed over the
 *   dimension "NumberIntervals". The time horizon is subdivided into
 *   NumberIntervals = k of the form [ 0 , i_1 ], [ i_1 + 1 , i_2 ], ...
 *   [ i_{k-1} + 1 , "TimeHorizon" - 1 ]; "ChangeIntervals" then has to
 *   contain [ i_1 , i_2 , ... i_{k-1} ]. Note that, therefore,
 *   "ChangeIntervals" has one significant value less than
 *   "NumberIntervals", which means that
 *   ChangeIntervals[ NumberIntervals - 1 ] is ignored. Anyway, the whole
 *   variable is ignored if either "NumberIntervals" <= 1 (such as if it
 *   is not defined), or "NumberIntervals" >= "TimeHorizon".
 */
 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// Generates the static variables of UnitBlock
/** The base UnitBlock class has four different "groups" of variables:
 *
 * - the commitment variables;
 *
 * - the primary spinning reserve variables;
 *
 * - the secondary spinning reserve variables;
 *
 * - the active power variables.
 *
 * All of these variables are optional, except the active power variables, in
 * the sense that the model may just not have them (say, because the unit
 * does not have 0-1 commitment decisions, or it cannot generate spinning
 * reserve). However, it is also possible to restrict which of the subsets
 * are generated with the parameter stvv.
 *
 * If stvv is not nullptr and it is a SimpleConfiguration<int>, or if
 * f_BlockConfig->f_static_variables_Configuration is not nullptr and it is
 * a SimpleConfiguration<int>, then the f_value (an int) indicates whether
 * each of the optional variables should be created. If the Configuration is
 * not available, the default value is taken to be 0. The value of the int
 * is interpreted bit-wise, with commitment variables being bit 0, primary
 * spinning reserve variables being bit 1, secondary spinning reserve
 * variables being bit 2, and active power variables being bit 3. If the bit
 * associated with a variable is 0 then the variable (assuming the model
 * actually has it) *is* created, otherwise it is *not*; hence, the default
 * value of 0 means that all the variables (that the model has) are created.
 *
 * Whenever a group of variables is created, its size will be the time
 * horizon.
 *
 * Note that derived classes are free to use the other bits of the int to
 * similarly encode for creation of their own specific groups of variables.
 */

 void generate_abstract_variables( Configuration * stvv ) override;

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE UnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
 /** @name Reading the data of the UnitBlock
  *
  * These methods allow to read data that must be common to (in principle) all
  * the kind of electrical generation units, i.e.:
  *
  * - fixed consumption when the unit is off;
  *
  * - the contribution to the inertia depending on the commitment status;
  *
  * - the contribution to the inertia depending on the active power produced.
  *
  * Note that this data is optional, and only "few" of the UnitBlock actually
  * have it. Hence, no data structures are defined by the base UnitBlock
  * class: the methods have a default implementation returning an empty
  * vector, and derived classes will have to handle their own data (if any).
  * @{ */

 /// Returns the time horizon of the problem
 Index get_time_horizon() const {
  return f_time_horizon;
 }
/*--------------------------------------------------------------------------*/
 /// Returns the number of electrical generators of each unit in the problem
 /** Returns the number of electrical generators for this UnitBlock. Since in
  *  most of the cases each unit has only one electrical generator, this
  *  method in the base UnitBlock class returns to one by default. Therefore,
  *  for all the units that have only one generator, the implementation of
  *  this method is already done right in the base UnitBlock class. The units
  *  that have more than one electrical generator (tied together by technical
  *  constraints) will have to handle this number by their-self. */

 virtual Index get_number_generators( void ) const {
  return( 1 );
 }
/*--------------------------------------------------------------------------*/
 /// Returns the empty matrix of fixed consumption
 /** The returned value U = get_fixed_consumption() contains the contribution
  *  to fixed consumption (basically, the constants to be multiplied by the
  *  commitment variables returned by get_commitment()) of all the generators
  *  at all time instants.
  *
  * The default implementation of the methods returns an empty matrix; and
  * derived classes will have to handle their own data (if any). */

 virtual const boost::multi_array< double , 2 > & get_fixed_consumption()
 const { const static boost::multi_array< double , 2 > _fc {};
  return ( _fc );
  }

/*--------------------------------------------------------------------------*/
 /// It returns the empty matrix of inertia commitment
 /** The returned value U = get_inertia_commitment() contains the contribution
  *  to inertia (basically, the constants to be multiplied by the commitment
  *  variables returned by get_commitment()) of all the generators at all time
  *  instants.
  *
  * The default implementation of the methods returns an empty matrix; and
  * derived classes will have to handle their own data (if any). */

 virtual const boost::multi_array< double , 2 > & get_inertia_commitment()
 const { const static boost::multi_array< double , 2 > _ic {};
  return ( _ic );
 }

/*--------------------------------------------------------------------------*/
 /// Returns the empty matrix of inertia power
 /** The returned value U = get_inertia_power() contains the contribution to
  *  inertia (basically, the constants to be multiplied by the active power
  *  variables returned by get_active_power()) of all the generators at all
  *  time instants.
  *
  * The default implementation of the methods returns an empty matrix; and
  * derived classes will have to handle their own data (if any).  */

 virtual const boost::multi_array< double , 2 > & get_inertia_power() const {
  const static boost::multi_array< double , 2 > _ip {};
  return ( _ip );
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
 *
 * All these groups of variables are (if not empty)
 * boost::multi_array< ColVariable , 2 > with first dimension time horizon
 * and second dimension number of generators.
 * @{ */

 /// Returns the matrix of commitment variables
 /** The returned boost::multi_array< ColVariable , 2 >, say U, contains the
  * commitment variables and is indexed over the dimensions time horizon and
  * get_number_generators(). There are two possible cases:
  *
  *  - if U is empty(), then these variables are not defined;
  *
  *  - otherwise, U must have f_time_horizon rows and get_number_generators()
  *    columns, and M[ t , g ] is the commitment variable for time step t of
  *    generator g. */

 const boost::multi_array< ColVariable , 2 > & get_commitment() const {
  return v_commitment;
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of primary spinning reserve variables
 /** The returned boost::multi_array< ColVariable , 2 >, say PR, contains the
  * primary spinning reserve variables and is indexed over the dimensions
  * time horizon and get_number_generators(). There are two possible cases:
  *
  *  - if PR is empty(), then these variables are not defined;
  *
  *  - otherwise, PR must have f_time_horizon rows and get_number_generators()
  *    columns, and PR[ t , g ] is the primary spinning reserve variable for
  *    time step t of generator g. */

 const boost::multi_array< ColVariable , 2 > & get_primary_spinning_reserve()
  const { return v_primary_spinning_reserve; }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of secondary spinning reserve variables
 /** The returned boost::multi_array< ColVariable , 2 >, say SR, contains the
  * secondary spinning reserve variables and is indexed over the dimensions
  * time horizon and get_number_generators(). There are two possible cases:
  *
  *  - if SR is empty(), then these variables are not defined;
  *
  *  - otherwise, SR must have f_time_horizon rows and get_number_generators()
  *    columns, and SR[ t , g ] is the secondary spinning reserve variable for
  *    time step t of generator g. */

 const boost::multi_array< ColVariable , 2 > &
  get_secondary_spinning_reserve() const {
  return v_secondary_spinning_reserve;
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of active power variables
 /** The returned boost::multi_array< ColVariable , 2 >, say P, contains the
  * active power variables and is indexed over the dimensions time horizon
  * and get_number_generators(). There are two possible cases:
  *
  *  - if P is empty(), then these variables are not defined;
  *
  *  - otherwise, P must have f_time_horizon rows and get_number_generators()
  *    columns, and P[ t , g ] is the secondary spinning reserve variable for
  *    time step t of generator g. */

 const boost::multi_array< ColVariable , 2 > & get_active_power() const {
  return v_active_power;
  }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Methods for loading, printing & saving the UnitBlock
  *  @{ */

 /// Extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  *  UnitBlock. See UnitBlock::deserialize( netCDF::NcGroup ) for
  *  details of the format of the created netCDF group.
  */
 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*---------------- METHODS FOR MODIFYING THE UnitBlock ---------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Methods for modifying the UnitBlock
  *  @{ */

 /// Sets the time horizon method
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

/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
 /** @name Handling the data of the UnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "UnitBlock::load() not implemented yet" ) );
 };

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// Utility method for resetting all the Variables
 void guts_of_destructor();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// the time horizon of the problem
 Index f_time_horizon;

 /// the number of intervals
 Index f_number_intervals;

 /// the vector of change intervals
 std::vector< Index > v_change_intervals;

 /// the matrix of commitment variables of generators
 boost::multi_array< ColVariable , 2 > v_commitment;

 /// the matrix of active power variables of generators
 boost::multi_array< ColVariable , 2 >v_active_power;

 /// the matrix of primary spinning reserve variables of generators
 boost::multi_array< ColVariable , 2 >v_primary_spinning_reserve;

 /// the matrix of secondary spinning reserve variables of generators
 boost::multi_array< ColVariable , 2 > v_secondary_spinning_reserve;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Returns which variables must be generated
 /** This method returns an int that indicates which variables of
  * UnitBlock must be generated by the generate_abstract_variables()
  * method. This value may be given in stvv as explained in
  * generate_abstract_variables(). If this value is not given in stvv, then
  * this method returns the appropriate value according to what is specified
  * in the generate_abstract_variables() method.
  */
 unsigned int get_variables_to_be_generated( Configuration * stvv );

 /// Deserializes the time horizon from a netCDF group
 void deserialize_time_horizon( netCDF::NcGroup & group );

 /// Deserializes the change intervals vector from a netCDF group
 void deserialize_change_intervals( netCDF::NcGroup & group );

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
