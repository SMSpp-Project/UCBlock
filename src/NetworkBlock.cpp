/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 23 - 05 - 2019
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

NetworkBlock::NetworkBlock( Block * block ) : Block( block ) { }

/*--------------------------------------------------------------------------*/

NetworkBlock::~NetworkBlock() { }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( netCDF::NcGroup & group ) {

  std::vector < size_t > start = { 0 };
  std::vector < size_t > count_nodes = { (size_t) f_number_nodes };
  std::vector < size_t > count_lines = { (size_t) f_number_lines };

  // Read the number of nodes and lines

  netCDF::NcDim number_nodes_NcDim = group.getDim( "NumberNodes" );
  if( number_nodes_NcDim.isNull() )
    throw( std::logic_error( "NumberNodes dimension is required" ) );
  f_number_nodes = number_nodes_NcDim.getSize();

  netCDF::NcDim number_lines_NcDim = group.getDim( "NumberLines" );
  if( number_lines_NcDim.isNull() )
    throw( std::logic_error( "NumberLines dimension is required" ) );
  f_number_lines = number_lines_NcDim.getSize();

  // Read starting lines

  netCDF::NcVar start_line_NcVar = group.getVar( "StartLine" );
  if( start_line_NcVar.isNull() )
    throw( std::logic_error( "StartLine not found" ) );

  v_startline.resize( f_number_nodes );
  start_line_NcVar.getVar( start , count_nodes , v_startline.data() );

  // Read ending lines

  netCDF::NcVar end_line_NcVar = group.getVar( "EndLine" );
  if( end_line_NcVar.isNull() )
    throw( std::logic_error( "EndLine not found" ) );

  v_endline.resize( f_number_nodes );
  start_line_NcVar.getVar( start , count_nodes , v_endline.data() );

  // Read active demand

  netCDF::NcVar active_demand_NcVar = group.getVar( "ActiveDemand" );
  if( active_demand_NcVar.isNull() )
    throw( std::logic_error( "ActiveDemand not found" ) );

  v_active_demand.resize( f_number_nodes );
  active_demand_NcVar.getVar( start , count_nodes , v_active_demand.data() );

  // Read susceptance

  netCDF::NcVar susceptance_NcVar = group.getVar( "Susceptance" );
  if( susceptance_NcVar.isNull() )
    throw( std::logic_error( "Susceptance not found" ) );

  v_susceptance.resize( f_number_lines );
  susceptance_NcVar.getVar( start , count_lines , v_susceptance.data() );

  // Read minimum power flow

  netCDF::NcVar min_power_flow_NcVar = group.getVar( "MinPowerFlow" );
  if( min_power_flow_NcVar.isNull() )
    throw( std::logic_error( "MinPowerFlow not found" ) );

  v_minimum_power_flow.resize( f_number_lines );
  min_power_flow_NcVar.getVar( start , count_lines , v_minimum_power_flow.data() );

  // Read maximum power flow

  netCDF::NcVar max_power_flow_NcVar = group.getVar( "MaxPowerFlow" );
  if( max_power_flow_NcVar.isNull() )
    throw( std::logic_error( "MaxPowerFlow not found" ) );

  v_maximum_power_flow.resize( f_number_lines );
  max_power_flow_NcVar.getVar( start , count_lines , v_maximum_power_flow.data() );

  // Issue Modification. Note: this is a NBModification, the "nuclear
  // option"

  if( anyone_there() )
    add_Modification( std::make_shared<NBModification>( this ) );

}  // end( NetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::generate_abstract_variables( Configuration *stvv ) {

  if( f_number_nodes < 0 ) {
    throw( std::logic_error( "NetworkBlock::generate_abstract_variables: "
                             "number of nodes of NetworkBlock is not set" ) );
  }

  if( v_node_injection.size() != f_number_nodes ) {
    assert( v_node_injection.size() == 0 ); // this should only happen once
    v_node_injection.resize( f_number_nodes );
    add_static_variable( v_node_injection );
  }
}

/*--------------------------------------------------------------------------*/

void NetworkBlock::generate_abstract_constraints( Configuration *stcc ) {

  if( f_number_lines < 0 ) {
    throw( std::logic_error( "NetworkBlock::generate_abstract_constraints: "
                             "number of lines of NetworkBlock is not set" ) );
  }

  if( v_flow_limit_constraints.size() != f_number_lines ) {
    // this should only happen once
    assert( v_flow_limit_constraints.size() == 0 );
    v_flow_limit_constraints.resize( f_number_lines );
  }

  // Flow limit constraints

  // TODO Put these constraints in the DCNetworkBlock when (and if) it
  // is created.

  for( int line_id = 0; line_id < f_number_lines; ++line_id ) {

    auto linear_function = new LinearFunction();
    double constant_term = 0;

    for( int node_id = 0; node_id < v_node_injection.size(); ++node_id ) {

      double coefficient = 0.0; // TODO Compute the Power Transfer
                                // Distribution Factor Matrix

      if( coefficient == 0.0 )
        continue;

      linear_function->add_variable
        ( & v_node_injection[node_id] , coefficient );

      constant_term -= coefficient * v_active_demand[node_id];

    } // for each node

    // Set the function of the constraint

    v_flow_limit_constraints[line_id].set_function( linear_function );

    // Set the left- and right-hand sides

    v_flow_limit_constraints[line_id].set_lhs
      ( v_minimum_power_flow[line_id] - constant_term );

    v_flow_limit_constraints[line_id].set_rhs
      ( v_maximum_power_flow[line_id] - constant_term );

  } // for each line

  add_static_constraint( v_flow_limit_constraints );
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::serialize( netCDF::NcGroup & group ) const {

 group.putAtt( "type" , "NetworkBlock" );

 netCDF::NcDim number_nodes_NcDim =
   group.addDim( "NumberNodes" , f_number_nodes );

 netCDF::NcDim number_lines_NcDim =
   group.addDim( "NumberLines" , f_number_lines );

 std::vector < size_t > start = { 0 };
 std::vector < size_t > count_nodes = { (size_t) f_number_nodes };
 std::vector < size_t > count_lines = { (size_t) f_number_lines };


  if( v_startline.size() )
    ( group.addVar( "StartLine" , netCDF::NcUint64() , number_nodes_NcDim )
    ).putVar( start , count_nodes , v_startline.data() );

  if( v_endline.size() )
    ( group.addVar( "EndLine" , netCDF::NcUint64() , number_nodes_NcDim )
    ).putVar( start , count_nodes , v_endline.data() );

 if( v_active_demand.size() )
   ( group.addVar( "ActiveDemand" , netCDF::NcDouble() , number_nodes_NcDim )
     ).putVar( start , count_nodes , v_active_demand.data() );

 if( v_susceptance.size() )
   ( group.addVar( "Susceptance" , netCDF::NcDouble() , number_lines_NcDim )
     ).putVar( start , count_lines , v_susceptance.data() );

 if( v_minimum_power_flow.size() )
   ( group.addVar( "MinPowerFlow" , netCDF::NcDouble() , number_lines_NcDim )
     ).putVar( start , count_lines , v_minimum_power_flow.data() );

 if( v_maximum_power_flow.size() )
   ( group.addVar( "MaxPowerFlow" , netCDF::NcDouble() , number_lines_NcDim )
     ).putVar( start , count_lines , v_maximum_power_flow.data() );

}    // end( NetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
