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
 *     (ii)  the power injected into the grid;
 *
 *     (iii) the primary spinning reserve of the unit;
 *
 *     (iv)  the secondary spinning reserve of the unit;
 *
 *     (v)  the active power produced by the unit.
 *
 *   Each of these vectors either have size equal to the time horizon
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the unit may not have reserve).
 *
 * \version 0.11
 *
 * \date 17 - 05 - 2019
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
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/** Constructor of UnitBlock, taking possibly a pointer of its father
 * Block and the time horizon. */
 UnitBlock( Block * father_block = nullptr , int t = 0 )
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
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 /// generate the static variables of UnitBlock
 /** Method that generates the static variables of this
  * UnitBlock. These may be:
  *
  * - the commitment variables
  *
  * - the power injected variables
  *
  * - the primary spinning reserve variables
  *
  * - the secondary spinning reserve variables
  *
  * - the active power variables
  *
  * All of these variables are optional, except the active power
  * variables. The parameter stvv is used to decide which of the
  * optional variables should be used. If stvv is not nullptr and it
  * is a SimpleConfiguration<int> or if
  * f_BlockConfig->f_static_constraints_Configuration is not nullptr
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
 int get_time_horizon( void ) const { return f_time_horizon; }

 /// Method for returning the vector of commitment variables
 const std::vector<ColVariable> & get_commitment( void ) const {
   return v_commitment;
 }

 /// Method for returning the vector of power injected variables
 const std::vector<ColVariable> & get_power_injected( void ) const {
   return v_power_injected;
 }

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

 /// Set the time horizon
 void set_time_horizon( int t );

/*@} -----------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UnitBlock
    @{ */

 virtual void load( std::istream &input ) override {};

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

 /// The time horizon of the problem
 int f_time_horizon;

 /* Each of the following vectors of Variables should either have size
  * f_time_horizon, meaning that there is one Variable for each time
  * step, or be empty, in which case the variables simply do not
  * exist. */

 /// Vector of commitment variables
 std::vector<ColVariable> v_commitment;

 /// Vector of power injected into the grid
 std::vector<ColVariable> v_power_injected;

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

};  // end( class( UnitBlock ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* UnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File UnitBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
