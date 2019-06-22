/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 21 - 06 - 2019
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>
#include "LinearFunction.h"
#include "NetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register NetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( NetworkBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( netCDF::NcGroup & group )
{
  ::deserialize_dim( group, "NumberNodes",   f_NetworkData->f_number_nodes );

  ::deserialize(group, "ActiveDemand", f_NetworkData->f_number_nodes,
                v_active_demand);

 }  // end( NetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

  void NetworkBlock::serialize(netCDF::NcGroup &group) const {

    group.putAtt("type", "NetworkBlock");
  auto dim_number_nodes = group.addDim("NumberNodes",
                          f_NetworkData ? f_NetworkData->f_number_nodes : 1);

    ::serialize(group, "ActiveDemand", netCDF::NcDouble(),
                {dim_number_nodes}, v_active_demand);

  }    // end( NetworkBlock::serialize )
/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
