/*--------------------------------------------------------------------------*/
/*-------------------- File BusNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BusNetworkBlock class.
 *
 * \version 0.11
 *
 * \date 22 - 06 - 2019
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>
#include "LinearFunction.h"
#include "NetworkBlock.h"
#include "BusNetworkBlock.h"
#include "FRowConstraint.h"
#include "UCBlock.h"

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

void BusNetworkBlock::deserialize( netCDF::NcGroup & group ) {


  unsigned int  number_nodes =  f_NetworkData ? f_NetworkData->
      f_number_nodes : 1;

  if (number_nodes == 1) {

    ::deserialize(group, "ActiveDemand", &f_active_demand);
  }
}  // end( BusNetworkBlock::deserialize )
/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_abstract_variables( Configuration *stvv ) {

  unsigned int number_nodes =  f_NetworkData ? f_NetworkData->
      f_number_nodes : 1;

  if (number_nodes == 1) {

    v_node_injection[number_nodes].set_value(f_active_demand);
    v_node_injection[number_nodes].is_fixed(true);
  } else {
    throw (std::logic_error("BusNetworkBlock has not define"));
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE BusNetworkBlock ----*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::serialize( netCDF::NcGroup & group ) const {

	group.putAtt( "type" , "BusNetworkBlock" );

  ::serialize( group, "ActiveDemand",
               netCDF::NcDouble(), f_active_demand );
}    // end( BusNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------------- End File BusNetworkBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
