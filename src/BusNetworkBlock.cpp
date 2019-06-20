/*--------------------------------------------------------------------------*/
/*-------------------- File BusNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BusNetworkBlock class.
 *
 * \version 0.11
 *
 * \date 05 - 06 - 2019
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

	// Default values for optional dimensions
	f_number_nodes = 1;

	::deserialize_dim(group, "NumberNodes", f_number_nodes);

	::deserialize(group, "ActiveDemand", f_number_nodes, v_active_demand);

}  // end( BusNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_abstract_constraints( Configuration *stcc ) {

	//TODO
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE BusNetworkBlock ----*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::serialize( netCDF::NcGroup & group ) const {

	group.putAtt( "type" , "BusNetworkBlock" );

	auto dim_number_nodes = group.addDim( "NumberNodes", f_number_nodes );

	::serialize( group, "ActiveDemand", netCDF::NcDouble(),
							 { dim_number_nodes}, v_active_demand);
}    // end( BusNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------------- End File BusNetworkBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
