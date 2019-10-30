/*--------------------------------------------------------------------------*/
/*-------------------- File BusNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BusNetworkBlock class.
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
 * \copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>
#include "NetworkBlock.h"
#include "BusNetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register NetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( BusNetworkBlock );

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_abstract_variables( Configuration * stvv ) {

 // In BusNetworkBlock, number_nodes = 1
 auto active_demand = get_active_demand()[ 0 ];
 v_node_injection[ 0 ].set_value( active_demand );
 v_node_injection[ 0 ].is_fixed( true );
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE BusNetworkBlock ----*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*--------------------- End File BusNetworkBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
