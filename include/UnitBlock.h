/*--------------------------------------------------------------------------*/
/*--------------------------- File UnitBlock.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class UnitBlock, which derives from the Block, in
 * order to define a base class for any possible unit that can be attached to
 * a UCBlock. It has very basic information that can characterize almost any
 * different kind of unit, which includes the lenght of the time horizon and
 * four sets of Variables: active power variables, commitment variables,
 * primary and secondary spinning reserve variables. It also outputs some
 * general information regarding how the active power and/or commitment
 * status of the unit at a given time instant impact the unit's capability
 * of satisfying inertia constraints, and the fixed consumption of the unit
 * (if any) when it is off.
 *
 * \version 0.11
 *
 * \date 19 - 06 - 2019
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
 * Copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
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

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS UnitBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// implementation of the Block concept for "a generic unit" in UC
/** The class UnitBlock, which derives from the Block, defines a base class
 * for any possible unit that can be attached to a UCBlock. It has very basic
 * information that can characterize almost any different kind of unit, which
 * includes four sets of Variables: power variables, commitment variables,
 * primary and secondary spinning reserve variables. This class has thus been
 * constructed having the following elements:
 *
 * - The time horizon of the problem.
 *
 * - Four vectors of ColVariable objects, that are used to store the
 *   information regarding:
 *
 *     (i)   the commitment of the unit;
 *
 *     (ii)  the primary spinning reserve of the unit;
 *
 *     (iii) the secondary spinning reserve of the unit;
 *
 *     (iv)  the active power produced by the unit;
 *
 *   Each of these vectors either have size equal to the time horizon
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the unit may not have reserve).
 *
 * The class also outputs some general information regarding how the active
 * power and/or commitment status of the unit at a given time instant impact
 * the unit's capability of satisfying inertia constraints, and the fixed
 * consumption of the unit (if any) when it is off. */

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
 * UnitBlock defines the following main public types:
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

/// constructor, takes the father and the time horizon
/** Constructor of UnitBlock, taking possibly a pointer of its father
 * Block and the time horizon. By default the time horizon is initialized
 * to 0, which means "not set yet". */
 explicit UnitBlock( Block * father_block = nullptr , Index t = 0 )
   : Block( father_block ), f_time_horizon( t ) {}

/*--------------------------------------------------------------------------*/

/// destructor of UnitBlock: it is virtual, and empty
 ~UnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the UnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
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
 *	  of course there is nothing preventing the same value to be
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
 *
 * - The variable "FixedConsumption", of type double and either indexed
 *   over the dimension "NumberIntervals", or having size 1. This is meant
 *   to represent the vector FC[ t ] which, for each time instant t,
 *   contains the fixed consumption of the power plant if it is OFF at time
 *   t. The variable is optional; if it is not defined, FC[ t ] == 0 for all
 *   time instants. If it is defined, it can either have size 1 or size
 *   "NumberIntervals". If it has size 1, then FC[ t ] ==
 *   FixedConsumption[ 0 ] for all t, regardless to what "NumberIntervals"
 *   says. Otherwise, FixedConsumption[ i ] is the fixed value of FC[ t ] for
 *   all t in the interval [ ChangeIntervals[ i - 1 ], ChangeIntervals[ i ] ],
 *   with the assumption that ChangeIntervals[ - 1 ] = 0.
 *
 * - The variable "InertiaCommitment", of type double and either indexed over
 *   the dimension "NumberIntervals" or has size 1. This is meant to
 *   represent the vector IC[ t ] which, for each time instant t, contains
 *   the contribution that the unit can give to the inertia constraint for
 *   the sole fact that is is on (basically, the constant to be multiplied to
 *   the commitment variable) at time t. The variable is optional; if it is
 *   not defined, IC[ t ] == 0 for all time instants. If it is defined, it
 *   can either have size 1 or size "NumberIntervals". If it has size 1, then
 *   IC[ t ] == InertiaCommitment[ 0 ] for all t, regardless to what
 *   "NumberIntervals" says. Otherwise, InertiaCommitment[ i ] is the
 *   fixed value of IC[ t ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0.
 *
 * - The variable "InertiaPower", of type double and either is indexed over
 *   the dimension "NumberIntervals" or has size 1. This is meant to
 *   represent the vector IP[ t ] which, for each time instant t, contains
 *   the contribution that the unit can give to the inertia constraint which
 *   depends on the active power that it is currently generating (basically,
 *   the constant to be multiplied to the active power variable) at time t.
 *   The variable is optional; if it is not defined, IP[ t ] == 0 for all
 *   time instants. If it is defined, it can either have size 1 or size
 *   "NumberIntervals". If it has size 1, then IP[ t ] == InertiaPower[ 0 ]
 *   for all t, regardless to what "NumberIntervals" says. Otherwise,
 *   InertiaPower[ i ] is the fixed value of IP[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static variables of UnitBlock
 /** Method that generates the static variables of this UnitBlock. The base
  * UnitBlock class has four different "groups" of variables:
  *
  * - the commitment variables;
  *
  * - the primary spinning reserve variables;
  *
  * - the secondary spinning reserve variables;
  *
  * - the active power variables.
  *
  * All of these variables are optional, except the active power variables,
  * in the sense that the model may just not have them (say, because the
  * unit does not have 0-1 commitment decisions, or it cannot generate
  * spinning reserve). However, it is also possible to restrict which of
  * the subsets are generated with the parameter stvv.
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

 void generate_abstract_variables( Configuration *stvv ) override;

/**@} ----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE UnitBlock -------------*/
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
 * @{ */

 /// method for returning the time horizon
 Index get_time_horizon() const { return( f_time_horizon ); }

/*--------------------------------------------------------------------------*/
 /// method for returning the vector of fixed consumption
 /** Method for returning the vector of fixed consumption. There are three
  * possible cases:
  *
  * - if the vector is empty, then the fixed consumption is 0;
  *
  * - if the vector only has one element, then the fixed consumption is
  *   always equal to the value of that element;
  *
  * - otherwise the vector must have the size of the time horizon, and the
  *   t-th entry gives the fixed consumption at time instant t.
  */
 const std::vector< double > & get_fixed_consumption() const {
  return( v_fixed_consumption );
  }

/*--------------------------------------------------------------------------*/
 /// returns the fixed consumption at time t
 /** Method for returning the fixed consumption at time t, for t
  * between 0 and time horizon minus 1.
  */
 double get_fixed_consumption( Index t ) const {
  return( v_fixed_consumption.empty() ? 0 :
	  v_fixed_consumption[ std::min( v_fixed_consumption.size() - 1 ,
					 t ) ] );
  }

/*--------------------------------------------------------------------------*/
 /// method for returning the vector of inertia commitment
 /** Method for returning the vector of "inertia commitment", i.e., the
  * contribution to inertia that the unit at each time t if the unit is on
  * (basically, the constants to be multiplied by the commitment variables).
  * There are three possible cases:
  *
  * - if the vector is empty, then the inertia commitment is 0;
  *
  * - if the vector only has one element, then the inertia commitment is
  *   always equal to the value of that element;
  *
  * - otherwise the vector must have the size of the time horizon, and the
  *   t-th entry gives the inertia commitment at time instant t.
  */

 const std::vector< double > & get_inertia_commitment() const {
  return( v_inertia_commitment );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia commitment at time t
 /** Method for returning the contribution to the inertia that the unit
  * can give at time t, for t between 0 and time horizon minus 1, if the
  * unit is on (basically, the constant to be multiplied by the commitment
  * variable at time t). */

 double get_inertia_commitment( Index t ) const {
  return( v_inertia_commitment.empty() ? 0 :
	  v_inertia_commitment[ std::min( v_inertia_commitment.size() - 1 ,
					  t ) ] );
  }

/*--------------------------------------------------------------------------*/
 /// method for returning the vector of inertia power
 /** Method for returning the vector of "inertia power", i.e., the
  * contribution to inertia that the unit at each time t that is proportional
  * to the active power generated (basically, the constants to be multiplied
  * by the active power variables). There are three possible cases:
  *
  * - if the vector is empty, then the inertia power is 0;
  *
  * - if the vector only has one element, then the inertia power is
  *   always equal to the value of that element;
  *
  * - otherwise the vector must have the size of the time horizon, and the
  *   t-th entry gives the inertia power at time instant t.
  */

 const std::vector< double > & get_inertia_power() const {
  return( v_inertia_power );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia power at time t
 /** Method for returning the contribution to the inertia that the unit
  * can give at time t, for t between 0 and time horizon minus 1, that is
  * proportional to the active power generated at that time (basically, the
  * constant to be multiplied by the active power variable at time t). */

 double get_inertia_power( Index t ) const {
  return( v_inertia_power.empty() ? 0 :
	  v_inertia_power[ std::min( v_inertia_power.size() - 1 , t ) ] );
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
 * @{ */

/// method for returning the vector of commitment variables

 const std::vector<ColVariable> & get_commitment() const {
  return( v_commitment );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// method for returning a pointer to the commitment variable at time t

 ColVariable & get_commitment( Index t )  {
  return( v_commitment[ t ] );
  }

/*--------------------------------------------------------------------------*/
/// method for returning the vector of primary spinning reserve variables
 const std::vector<ColVariable> & get_primary_spinning_reserve() const {
  return( v_primary_spinning_reserve );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// method for returning a pointer to the primary spinning reserve at time t

 ColVariable & get_primary_spinning_reserve( Index t )  {
  return( v_primary_spinning_reserve[ t ] );
  }

/*--------------------------------------------------------------------------*/
/// method for returning the vector of secondary reserve variables
 const std::vector<ColVariable> & get_secondary_spinning_reserve()
   const { return( v_secondary_spinning_reserve ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// method for returning a pointer to the secondary spinning reserve at time t

 ColVariable & get_secondary_spinning_reserve( Index t )  {
  return( v_secondary_spinning_reserve[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// method for returning the vector of power variables
 const std::vector<ColVariable> & get_active_power() const {
  return( v_active_power );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// Method for returning the pointer to the power variable at time t

 ColVariable & get_active_power( Index t ) { return( v_active_power[ t ] ); }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UnitBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the UnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * UnitBlock. See UnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*----------------- METHODS FOR MODIFYING THE UnitBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the UnitBlock
 *  @{ */

 /// set the time horizon method
 /** Set the time horizon method.
  *
  * This method can be called *before* that deserialize() is called to
  * provide the UnitBlock with the time horizon. This allows the information
  * not to be duplicated in the netCDF group that describes the unit, since
  * usually (bit not necessarily) a UnitBlock is deserialized inside a
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
  * the unit has any data indiced over the time horizon) has to agree with the
  * value set by this method.
  *
  * If this method is called *after* that deserialize() is called, this is
  * taken to mean that the UnitBlock is being "reset", and that immediately
  * after deserialize() will be called again. The same rules as above are to
  * be followed for that subsequent call to deserialize(). */

 void set_time_horizon( Index t ) { f_time_horizon = t; }

/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UnitBlock
    @{ */

 void load( std::istream &input ) override {
  throw( std::logic_error( "UnitBlock::load() not implemented yet" ) );
  };

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// The time horizon of the problem
 Index f_time_horizon;

 /// the number of intervals
 Index f_number_intervals{};

 /// the vector of change intervals
 std::vector< Index > v_change_intervals;

 /// Vector of fixed consumption
 std::vector< double > v_fixed_consumption;

 /// Vector of inertia commitment
 std::vector< double > v_inertia_commitment;

 /// Vector of inertia power
 std::vector< double > v_inertia_power;

 /* Each of the following vectors of Variables should either have size
  * f_time_horizon, meaning that there is one Variable for each time
  * step, or be empty, in which case the variables simply do not
  * exist. */

 /// Vector of commitment variables
 std::vector<ColVariable> v_commitment;

 /// Vector of power variables
 std::vector<ColVariable> v_active_power;

 /// Vector of primary spinning reserve variables
 std::vector<ColVariable> v_primary_spinning_reserve;

 /// Vector of secondary spinning reserve variables
 std::vector<ColVariable> v_secondary_spinning_reserve;

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

 // returns which variables must be generated
 /* This method returns an int that indicates which variables of UnitBlock
  * must be generated by the generate_abstract_variables() method. This value
  * may be given in stvv as explained in generate_abstract_variables(). If
  * this value is not given in stvv, then this method returns the appropriate
  * value according to what is specified in the generate_abstract_variables()
  * method. */

 unsigned int get_variables_to_be_generated( Configuration *stvv );

 void deserialize_time_horizon( netCDF::NcGroup & group );

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
/*------------------------- End File UnitBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
