/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 14 - 06 - 2019
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
 ::deserialize(group, "ActiveDemand", f_number_nodes, v_active_demand);

 }  // end( NetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::generate_abstract_variables(Configuration *stvv)
{
 if( f_number_nodes < 0 ) {
      throw (std::logic_error("NetworkBlock::generate_abstract_variables: "
                              "number of nodes of NetworkBlock is not set"));
    }

    if (v_node_injection.size() != f_number_nodes) {
      assert(v_node_injection.size() == 0); // this should only happen once
      v_node_injection.resize(f_number_nodes);
      add_static_variable(v_node_injection);
    }
  }
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

  void NetworkBlock::serialize(netCDF::NcGroup &group) const {

    group.putAtt("type", "NetworkBlock");

    auto dim_number_nodes = group.addDim("NumberNodes", f_number_nodes);
    auto dim_number_lines = group.addDim("NumberLines", f_number_lines);

    ::serialize(group, "ActiveDemand", netCDF::NcDouble(),
                {dim_number_nodes}, v_active_demand);

    if (f_number_nodes > 1) {
      ::serialize(group, "StartLine", netCDF::NcUint64(),
                  {dim_number_nodes}, v_start_line);

      ::serialize(group, "EndLine", netCDF::NcUint64(),
                  {dim_number_nodes}, v_end_line);

      ::serialize(group, "MinPowerFlow", netCDF::NcDouble(),
                  {dim_number_lines}, v_min_power_flow);

      ::serialize(group, "MaxPowerFlow", netCDF::NcDouble(),
                  {dim_number_lines}, v_max_power_flow);

      ::serialize(group, "Susceptance", netCDF::NcDouble(),
                  {dim_number_lines}, v_susceptance);
    }

  }    // end( NetworkBlock::serialize )
/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
