/*--------------------------------------------------------------------------*/
/*------------------------ File NetworkData.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkData class.
 *
 * \version 0.11
 *
 * \date 17 - 06 - 2019
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
 * Copyright &copy by Antonio Frangioni and Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <random>
#include <map>

#include "NetworkData.h"
#include "UCBlock.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register NetworkData the Block factory

SMSpp_insert_in_factory_cpp_1( NetworkData );

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF NetworkData ------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void NetworkData::deserialize( netCDF::NcGroup & group ) {

  // Default values for optional dimensions
  f_number_nodes           = 1;
  f_number_lines           = 0;

  ::deserialize_dim( group, "NumberNodes",          f_number_nodes );
  ::deserialize_dim( group, "NumberLines",          f_number_lines );

  if( f_number_nodes > 1 ) {
    ::deserialize( group, "StartLine", f_number_nodes, v_start_line );
  }

  if( f_number_nodes > 1 ) {
    ::deserialize( group, "EndLine", f_number_nodes, v_end_line);
  }


  if( f_number_nodes > 1 ) {
    ::deserialize( group, "MinPowerFlow", f_number_lines , v_min_power_flow );
  }

  if( f_number_nodes > 1 ) {
    ::deserialize( group, "MaxPowerFlow", f_number_lines, v_max_power_flow );
  }

  if( f_number_nodes > 1 ) {
    ::deserialize( group, "Susceptance", f_number_lines, v_susceptance);
  }

}  // end( NetworkData::deserialize )

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkData ---------*/
/*--------------------------------------------------------------------------*/

void NetworkData::serialize( netCDF::NcGroup & group ) const {


  group.putAtt( "type" , "NetworkData" );

  auto dim_number_nodes = group.addDim( "NumberNodes", f_number_nodes );
  auto dim_number_lines = group.addDim( "NumberLines", f_number_lines );

  if( f_number_nodes > 1 ) {
    ::serialize( group, "StartLine", netCDF::NcUint64(),
                 { dim_number_nodes }, v_start_line );

    ::serialize( group, "EndLine", netCDF::NcUint64(),
                 { dim_number_nodes }, v_end_line );

    ::serialize( group, "MinPowerFlow", netCDF::NcDouble(),
                 { dim_number_lines }, v_min_power_flow );

    ::serialize( group, "MaxPowerFlow", netCDF::NcDouble(),
                 { dim_number_lines }, v_max_power_flow );

    ::serialize( group, "Susceptance", netCDF::NcDouble(),
                 { dim_number_lines }, v_susceptance );
  }

}  // end( NetworkData::serialize )

/*--------------------------------------------------------------------------*/
/*---------------------- End File NetworkData.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
