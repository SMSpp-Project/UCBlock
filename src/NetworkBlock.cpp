/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 21 - 03 - 2020
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
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

NetworkBlock::NetworkData::NetworkData() {
 f_number_lines = 0;
 f_number_nodes = 0;
}

void NetworkBlock::NetworkData::deserialize( netCDF::NcGroup & group ) {

 if( ! ::deserialize_dim( group, "NumberNodes", f_number_nodes, true ) )
  f_number_nodes = 1;

 if ( f_number_nodes > 1  ) {  // DCNetworkBlock

  ::deserialize_dim( group, "NumberLines", f_number_lines, false );

  ::deserialize( group, "StartLine", f_number_lines, v_start_line, false, true );

  ::deserialize( group, "EndLine", f_number_lines, v_end_line, false, true );

  ::deserialize( group, "MinPowerFlow", f_number_lines, v_min_power_flow, true, true );

  ::deserialize( group, "MaxPowerFlow", f_number_lines, v_max_power_flow, true, true );

  ::deserialize( group, "Susceptance", f_number_lines, v_susceptance, true, true );
 }

}

/*--------------------------------------------------------------------------*/
void NetworkBlock::deserialize( netCDF::NcGroup & group ) {

 auto network_data = new NetworkBlock::NetworkData();
 network_data->deserialize( group );

 auto dim_number_nodes = group.getDim( "NumberNodes" );

 if( !dim_number_nodes.isNull() )
  ::deserialize( group, "ActiveDemand", dim_number_nodes.getSize(),
                 v_active_demand );

 Block::deserialize( group );
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::serialize( netCDF::NcGroup & group ) const {

 auto dim_number_nodes = group.addDim( "NumberNodes", f_number_nodes );

 if( f_number_nodes > 1 ) { // DCNetworkBlock
  auto dim_number_lines = group.addDim( "NumberLines", f_number_lines );

  ::serialize( group, "StartLine", netCDF::NcUint64(),
               { dim_number_lines }, v_start_line );

  ::serialize( group, "EndLine", netCDF::NcUint64(),
               { dim_number_lines }, v_end_line );

  ::serialize( group, "MinPowerFlow", netCDF::NcDouble(),
               { dim_number_lines }, v_min_power_flow );

  ::serialize( group, "MaxPowerFlow", netCDF::NcDouble(),
               { dim_number_lines }, v_max_power_flow );

  ::serialize( group, "Susceptance", netCDF::NcDouble(),
               { dim_number_lines }, v_susceptance );
 }

}

/*--------------------------------------------------------------------------*/
void NetworkBlock::serialize( netCDF::NcGroup & group ) const {

 Block::serialize( group );

 auto dim_number_nodes = group.getDim( "NumberNodes" );

 if( !dim_number_nodes.isNull() )
  ::serialize( group, "ActiveDemand", netCDF::NcDouble(),
               { dim_number_nodes }, v_active_demand );

} // end( NetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
