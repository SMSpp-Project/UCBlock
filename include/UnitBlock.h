/*--------------------------------------------------------------------------*/
/*--------------------------- File UnitBlock.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the *derived* class UnitBlock, which derives from
 * the Block, in order to define a base for any possible unit that can
 * be attached to a UCBlock. It has very basic information that can
 * characterize almost any different kind of unit, which includes four
 * sets of Variables: power variables, commitment variables, and the
 * primary and secondary spinning reserve variables. This class has
 * thus been constructed having the following elements:
 *
 * - A virtual public method that is used to initialize and read the
 *   data of any possible derived UnitBlock class.
 *
 * - The time horizon of the problem.
 *
 * - Four vectors of ColVariable objects, that are used to store the
 *   information regarding:
 *
 *     (i)   the commitment of the unit;
 *
 *     (ii) the primary spinning reserve of the unit;
 *
 *     (iii)  the secondary spinning reserve of the unit;
 *
 *     (iv)  the active power produced by the unit.
 *
 *   Each of these vectors either have size equal to the time horizon
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the unit may not have reserve).
 *
 * \version 0.11
 *
 * \date 06 - 06 - 2019
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __UnitBlock
#define __UnitBlock
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "ColVariable.h"
#include <map>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

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

typedef unsigned int Index;                 ///< index of parameters

/*@}------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/** Constructor of UnitBlock, taking possibly a pointer of its father
 * Block and the time horizon. */
 UnitBlock( Block * father_block = nullptr , Index t = 0 )
   : Block( father_block ), f_time_horizon( t ) {}

/*--------------------------------------------------------------------------*/

 /// destructor of UnitBlock: it is virtual, and empty
 virtual ~UnitBlock() { }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the UnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 * - the dimension "TimeHorizon" containing the time horizon.
 *
 * Note 1: consider set time horizon \f$\{0, \dots,
 * "TimeHorizon-1"\}\f$ with dimension "TimeHorizon". It can be
 * presented as the union of some intervals.
 *
 * Note 2: let's suppose the values of each variable may change independently,
 * and may have different values for some intervals along the set time horizon
 * \f$\{0, \dots, "TimeHorizon-1"\}\f$. It means each variable has its own
 * changes along some intervals independently. Without loss of generality,
 * let's take the union of the intersection of all the intervals between each
 * two variables separately. It means, at the end we may have a set of
 * intervals such as: \f$ [0 , a], [a+1 , b], \dots, [j+1 , k], [k+1 ,
 * TimeHorizon-1] \f$which where they cover all the changes for all the
 * variables. Then in each interval, one variable may change or not(if not,
 * copy the corresponding value of its' previous interval).
 *
 * - the dimension "NumberIntervals" which is a subset of \f$ \{1, ...,
 *   "TimeHorizon"\}\f$ and indicates the number of above intervals \f$([0
 *   , a] , [a+1 , b], ... , [j+1 , k] , [k+1, TimeHorizon-1])\f$ where the
 *   variables change. This dimension is optional. If it is not
 *   provided then it is taken to be 1.
 *
 *   Three scenarios may happen:
 *
 *    i). In the simplest case scenario, "NumberIntervals = 1" which
 *        means that the value of each variable does not change, i.e.,
 *        it is the same for each period in \f$ \{1, \dots ,
 *        "TimeHorizon"\}\f$.
 *
 *   ii). In the average case scenario, "1 < NumberIntervals <
 *        TimeHorizon", which means that in some time steps the values
 *        of some of the variables are changing.
 *
 *  iii). In the worst case scenario, "NumberIntervals = TimeHorizon",
 *        which means that in every time step, the value of each
 *        variable may change.
 *
 * - the variable "ChangeIntervals", of type integer and indexed over
 *   the dimension "NumberIntervals"; the \f$t_{th}\f$ entry of the
 *   variable indicates the positive number of \f$ a, b, ..., k,
 *   TimeHorizon-1\f$ on the above example.
 *
 * - the variable "FixedConsumption", of type double and indexed over the
 *   dimension "TimeHorizon"; the entry  FixedConsumption[ t ] shows the
 *   fixed consumption of the power plant when it is OFF; this means that when
 *   the unit is off, it has a fixed consumption at time t which is a negative
 *   value; the variable is optional; if not defined,
 *   FixedConsumption[ t ] == 0 for all t; if it is defined it can either have
 *   size 1 or size TimeHorizon if it has size 1 then the value is the same
 *   for all t;
 *
 * - the variable "InertiaCommitment", of type double and indexed over the
 *   dimension "TimeHorizon"; the entry  InertiaCommitment[ t ] shows the
 *   value of inertia commitment parameter for each thermal unit; this means
 *   that the unit gives a contribution to the inertia at time t which is
 *   get_commitment()[ t ] * InertiaCommitment[ t ]; the variable is optional.
 *   if not defined, InertiaCommitment[ t ] == 0 for all t. if it is defined
 *   it can either have size 1 or size TimeHorizon if it has size 1 then the
 *   value is the same for all t;
 *
 * - the variable "InertiaPower", of type double and indexed over the
 *   dimension "TimeHorizon"; the entry InertiaPower[ t ] shows the amount of
 *   inertia power value for the unit at time t; this means that the unit
 *   gives a contribution to the inertia at time t which is
 *   get_active_power()[ t ] * InertiaCommitment[ t ]; the variable is
 *   optional. if not defined, InertiaCommitment[ t ] == 0 for all t. if it is
 *   defined it can either have size 1 or size TimeHorizon if it has size 1
 *   then the value is the same for all t
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 /// generate the static variables of UnitBlock
 /** Method that generates the static variables of this
  * UnitBlock. These may be:
  *
  * - the commitment variables [bit 0]
  *
  * - the primary spinning reserve variables [bit 1]
  *
  * - the secondary spinning reserve variables [bit 2]
  *
  * - the active power variables [bit 3]
  *
  * All of these variables are optional, except the active power
  * variables. The parameter stvv is used to decide which of the
  * optional variables should be used. If stvv is not nullptr and it
  * is a SimpleConfiguration<int> or if
  * f_BlockConfig->f_static_variables_Configuration is not nullptr
  * and it is a SimpleConfiguration<int>, then the f_value (an int of
  * this (the first possible) Configuration indicates whether each of
  * the optional variables should be used. This is done according to
  * the corresponding bit of f_value. If the bit associated with a
  * variable is 1, then the variable should be used; otherwise, it
  * should not. The first bit is associated with the commitment
  * variables, the second one with the primary spinning reserve
  * variables and so on according to the order the variables are
  * listed above. Whenever a group of variables should be used, its
  * size will be the time horizon. In any other case (i.e., if a
  * SimpleConfiguration<int> is not present), none of the variables
  * above is considered. */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
   override ;

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE UnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UnitBlock
    @{ */

 /// Method for returning the time horizon
 Index get_time_horizon( void ) const { return f_time_horizon; }

 /// Method for returning the vector of fixed consumption
 const std::vector< double > & get_fixed_consumption( void ) const {
   return v_fixed_consumption;
 }

 /// Method for returning the vector of inertia commitment
 const std::vector< double > & get_inertia_commitment( void ) const {
   return v_inertia_commitment;
 }

 /// Method for returning the vector of inertia power
 const std::vector< double > & get_inertia_power( void ) const {
   return v_inertia_power;
 }

 /// Method for returning the vector of commitment variables
 const std::vector<ColVariable> & get_commitment( void ) const {
   return v_commitment;
 }

 /// Method for returning the pointer to the commitment variable at time t
 ColVariable * get_commitment( int t ) { return & ( v_commitment[ t ] ); }

 /// Method for returning the vector of primary spinning reserve variables
  const std::vector<ColVariable> & get_primary_spinning_reserve( void ) const {
   return v_primary_spinning_reserve;
 }

 /// Method for returning the vector of secondary reserve variables
 const std::vector<ColVariable> & get_secondary_spinning_reserve( void ) const {
   return v_secondary_spinning_reserve;
 }

 /// Method for returning the vector of power variables
 const std::vector<ColVariable> & get_active_power( void ) const {
   return v_active_power;
 }

 /// Method for returning the pointer to the power variable at time t
 ColVariable * get_active_power( int i ) { return & ( v_active_power[i] ); }

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UnitBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the UnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * UnitBlock. See UnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 */

 virtual void serialize( netCDF::NcGroup & group ) const override;

/*@} -----------------------------------------------------------------------*/
/*----------------- METHODS FOR MODIFYING THE UnitBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the UnitBlock
 *  @{ */

 /// Set the time horizon \\TODO
 /** Better comment again when this is supposed to be called and why: should
  * be called by the father just before calling deserialize(), it says tha
  * the object is just going to be deserialized and it should use the value
  * passed here if it does not find one in the netCDF file. If there is a
  * value in the netCDF and the two disagre ...
  */

 void set_time_horizon( Index t );

/*@} -----------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UnitBlock
    @{ */

 virtual void load( std::istream &input ) override {};


/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// The time horizon of the problem
 Index f_time_horizon;

 /// the number of intervals
 Index f_number_intervals;

 /// the vector of change intervals
 std::vector< int > v_change_intervals;

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

 /// returns which variables must be generated
 /** This method returns an int that indicates which variables of
  * UnitBlock must be generated by the generate_abstract_variables()
  * method. This value may be given in stvv as explained in
  * generate_abstract_variables(). If this value is not given in stvv,
  * then this method returns the appropriate value according to what
  * is specified in the generate_abstract_variables() method. */
 int get_variables_to_be_generated( Configuration *stvv );

 void deserialize_time_horizon( netCDF::NcGroup & group );

 void deserialize_change_intervals( netCDF::NcGroup & group );

};  // end( class( UnitBlock ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* UnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File UnitBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
