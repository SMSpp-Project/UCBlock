/*--------------------------------------------------------------------------*/
/*--------------------------- File UCBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class UCBlock, which implements the
 * base class Block, in order to define a very basic Unit Commitment
 * Block, that will be able to fit and be used as a base for almost
 * any different variation of the Unit Commitement Problem. As a
 * result of this, the UCBlock class is characterized by the following
 * ingredients:
 *
 * - A time horizon of the optimization problem.
 *
 * - A Network, that defines the topology of the network.
 *
 * - A set of units (UnitBlock; that may be referring to any different
 *   kind of unit, such as thermal, hydro, etc).
 *
 * - A set of networks (NetworkBlock; that can potentially refer to
 *   demand satisfaction, transmission line capacities, etc.). This
 *   set is either empty, which means that there is no network in the
 *   model, or has size equals the time horizon.
 *
 * - A vector of node injection ColVariables.
 *
 * - A multi_array of node injection constraints.
 *
 * - Some basic public methods to read the data and initialize the
 *   optimisation problem.
 *
 * \version 0.11
 *
 * \date 07 - 03 - 2019
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
#ifndef __UCBlock
#define __UCBlock
/* self-identification: #endif at the end of the file */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <boost/multi_array.hpp>
#include <vector>
#include "Block.h"
#include "ColVariable.h"
#include "Network.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class FRowConstraint; ///< forward declaration of FRowConstraint
class NetworkBlock;   ///< forward declaration of NetworkBlock
class UnitBlock;      ///< forward declaration of UnitBlock

class UCBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of UCBlock, taking possibly a pointer of its father Block
 UCBlock( Block *father = nullptr ) : Block( father ) {}

/*--------------------------------------------------------------------------*/

 virtual ~UCBlock() {};   ///< destructor of UCBlock: it is virtual, and empty

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
   override;

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
 *  @{ */

 /// Method that initiazes the instance and passes all the needed data
 void instance( std::istream& inStream );

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
    @{ */

 /// returns the time horizon of the problem
 int get_time_horizon( void ) const { return f_time_horizon; }

 /// sets the number of units of the problem
 void set_units_size( int units ) { f_units_size = units; }

 /// returns the vector of (pointers to) NetworkBlocks
 const std::vector<NetworkBlock *> & get_network_blocks( void ) {
   return v_network_blocks;
 }

 /// returns the Network
 const Network & get_network( void ) { return f_network; }

/*@} -----------------------------------------------------------------------*/
/*------------------ METHODS FOR MODIFYING THE UCBlock ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the UCBlock
 *  @{ */

 /// sets the time horizon of the problem
 void set_time_horizon( int t ) { f_time_horizon = t; }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

 /// The time horizon of the problem
 int f_time_horizon;

 /// The number of units of the problem
 int f_units_size;

 /** Vector of pointers to the NetworkBlocks. This vector either is
  * empty or has size f_time_horizon. If it is empty, it means there
  * is no network. If it has positive size, then the NetworkBlock at
  * position i in this vector refers to the network at the i-th time
  * step. */
 std::vector<NetworkBlock *> v_network_blocks;

 /// The network
 Network f_network;

 /// Node injection constraints at each time and for
 boost::multi_array<FRowConstraint *, 2> v_node_injection_constraints;

};   // end( class( UCBlock ) )

} /* namespace SMSpp_di_unipi_it */

#endif /* UCBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File UCBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
