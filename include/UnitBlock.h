/*--------------------------------------------------------------------------*/
/*--------------------------- File UnitBlock.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the *derived* class UnitBlock, which derives from
 * the Block, in order to define a base for any possible unit that can
 * be attached to a UCBlock. It has very basic information that can
 * characterize almost any different kind of unit, which includes five
 * sets of Variables: power variables, commitment variables, primary
 * and secondary spinning reserve variables, and the heat
 * variables. This class has thus been constructed having the
 * following elements:
 *
 * - A virtual public method is used in order to initialize and read
 *   the data of any possible derived UnitBlock class.
 *
 * - A factory that is used in order for any possible derived class to
 *   be able to "automatically" register itself and be
 *   initialized. The factory is defined by a static method that
 *   initializes the static map that is used in order to store all the
 *   different possible derived classes that are linked with a unique
 *   string.
 *
 * - The time horizon of the problem.
 *
 * - Five vectors of ColVariable objects, that are used to store the
 *   information regarding:
 *
 *     (i)   the commitment of the unit;
 *
 *     (ii)  the power produced by the unit;
 *
 *     (iii) the primary spinning reverve of the unit;
 *
 *     (iv)  the secondary spinning reserve of the unit;
 *
 *     (v)   the heat produced by the unit.
 *
 *   Each of these vectors either have size equal to the time horizon
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the unit may not have reserve).
 *
 * \version 0.11
 *
 * \date 05 - 03 - 2019
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
#define __UnitBlock  /* self-identification: #endif at the end
                      * of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "boost/function.hpp"
#include "Block.h"
#include "ColVariable.h"
#include <map>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class UCBlock; ///< forward declaration of UCBlock

class UnitBlock : public Block {

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED TYPES OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/** Definition of the UnitFactory, used to properly initialize all the
 * different possible derived classes of UnitBlock */
 typedef boost::function <UnitBlock * (UCBlock *)> UnitFactory;

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

 /// generate the static variables of UnitBlock
 /** Method that generates the static variables of this
  * UnitBlock. These may be:
  *
  * - the commitment variables
  *
  * - the power variables
  *
  * - the primary spinning reserve variables
  *
  * - the secondary spinning reserve variables
  *
  * - the heat variables
  *
  * All of these variables are optional. The parameter stvv is used to
  * decide which of these variables should be used. If stvv is not
  * nullptr and it is a SimpleConfiguration<int> or if
  * f_BlockConfig->f_static_constraints_Configuration is not nullptr
  * and it is a SimpleConfiguration<int>, then the f_value (an int of
  * this (the first possible) Configuration indicates whether each of
  * the variables should be used. This is done according to the
  * corresponding bit of f_value. If the bit associated with a
  * variable is 1, then the variable should be used; otherwise, it
  * should not. The first bit is associated with the commitment
  * variables, the second one with the power variables and so on
  * according to the order the variables are listed above. Whenever a
  * group of variables should be used, its size will be the time
  * horizon. In any other case (i.e., if a SimpleConfiguration<int> is
  * not present), none of the variables above are considered. */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
   override;

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE UnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UnitBlock
    @{ */

 /// Method for returning the time horizon
 int get_time_horizon( void ) const { return f_time_horizon; }

 /// Method for returning the vector of commitment variables
 std::vector<ColVariable>* get_commitment( void ) { return &v_commitment; }

 /// Method for returning the i-th commitment variables
 ColVariable * get_commitment( int i ) { return &v_commitment[i]; }

 /// Method for returning the total number of commitment variable
 int get_commitment_size( void ) { return v_commitment.size(); }

 /// Method for returning the vector of power variables
 std::vector<ColVariable> * get_power( void ) { return &v_power; }

 /// Method for returning the i-th power variable
 ColVariable * get_power( int i ) { return &v_power[i]; }

 /// Method for returning the total number of power variables
 int get_power_size( void ) { return v_power.size(); }

 /// Method for returning the vector of primary spinning reserve variables
 std::vector<ColVariable> * get_primary_spinning_reserve( void ) {
   return &v_primary_spinning_reserve;
 }

 /// Method for returning the i-th primary reserve variable
 ColVariable * get_primary_spinning_reserve( int i ) {
   return &v_primary_spinning_reserve[i];
 }

 /// Method for returning the total number of primary reserve variables
 int get_primary_spinning_reserve_size( void ) {
   return v_primary_spinning_reserve.size();
 }

 /// Method for returning the vector of secondary reserve variables
 std::vector<ColVariable> * get_secondary_spinning_reserve( void ) {
   return &v_secondary_spinning_reserve;
 }

 /// Method for returning the i-th secondary reserve variable
 ColVariable * get_secondary_spinning_reserve( int i ) {
   return &v_secondary_spinning_reserve[i];
 }

 /// Method for returning the total number of secondary reserve variables
 int get_secondary_spinning_reserve_size( void ) {
   return v_secondary_spinning_reserve.size();
 }

 /// Method for returning the vector of heat variables
 std::vector<ColVariable> * get_heat( void ) { return &v_heat; }

 /// Method for returning the i-th heat variable
 ColVariable * get_heat( int i ) { return &v_heat[i]; }

 /// Method for returning the total number of heat variables
 int get_heat_size( void ) { return v_heat.size(); }

/*@} -----------------------------------------------------------------------*/
/*----------------- METHODS FOR MODIFYING THE UnitBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the UnitBlock
 *  @{ */

 /// Set the time horizon
 void set_time_horizon( int t ) { f_time_horizon = t; }

/*@} -----------------------------------------------------------------------*/
/*------------------ METHODS FOR INITIALIZING THE UnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UnitBlock
    @{ */

/** Static member method for initializing the map of the factory
 * of UnitBlock */
 static std::map<std::string, UnitBlock::UnitFactory>& U_factory();

 virtual void load( std::istream &input ) {};

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

 /// Vector of power variables
 std::vector<ColVariable> v_power;

 /// Vector of primary spinning reserve variables
 std::vector<ColVariable> v_primary_spinning_reserve;

 /// Vector of secondary spinning reserve variables
 std::vector<ColVariable> v_secondary_spinning_reserve;

 /// Vector of heat variables
 std::vector<ColVariable> v_heat;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

 SMSpp_insert_in_factory_h;

};  // end( class( UnitBlock ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* UnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File UnitBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
