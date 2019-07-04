/*--------------------------------------------------------------------------*/
/*-------------------------- File BusNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BusNetworkBlock, which derives from NetworkBlock
 * [see NetworkBlock.h] in order to define a "bus" transmission network in
 * the Unit Commitment problem.
 *
 * \version 0.11
 *
 * \date 01 - 07 - 2019
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

#ifndef __BusNetworkBlock
#define __BusNetworkBlock
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BusNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// A NetworkBlock with only one node, i.e., a "bus" transmission network
/** The BusNetworkBlock class, which derives from NetworkBlock [see
 *  NetworkBlock.h] implements the Block concept [see Block.h] in order to
 *  define a "bus" transmission network in the Unit Commitment problem for a
 *  given instant in the time horizon. */

class BusNetworkBlock : public NetworkBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// Constructor of BusNetworkBlock, taking possibly a pointer to its father

 explicit BusNetworkBlock( Block * f_block = nullptr ) :
  NetworkBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// Destructor of BusNetworkBlock, (understandably) does nothing

 ~BusNetworkBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Loads the BusNetworkBlock instance from memory
/** Like load( std::istream & ), if there is any Solver attached to this
 *  BusNetworkBlock then a NBModification (the "nuclear option") is issued.
 */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "BusNetworkBlock::load() not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
/// Generates the static variables of BusNetworkBlock
/** The base BusNetworkBlock class has just the node injection variables.
 * Since a "bus" network has just one node, and therefore a single value D for
 * the demand and a single injection variable s, which can hardly be called a
 * variable since the only possible way to satisfy the constraints is by
 * having s = D which in fact makes the variable a constant.*/

 void generate_abstract_variables( Configuration * stvv ) override;

/*--------------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE BusNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the BusNetworkBlock
 *  @{ */

 /// Method to set the NetworkData object
 /** This method does nothing in BusNetworkBlock because by definition
  * the network is made by only one node. */

 void set_NetworkData( NetworkBlock::NetworkData * nd ) override {}

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};   // end( class( BusNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  /* end namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* BusNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BusNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
