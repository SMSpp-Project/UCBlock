/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 04 - 03 - 2019
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
#include "UCBlock.h"

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

void NetworkBlock::generate_abstract_constraints( Configuration *stcc ) {

  if( ! f_network ) {
    throw( std::logic_error( "NetworkBlock::generate_abstract_constraints: "
			     "f_network of NetworkBlock is not set" ) );
  }

  auto num_lines = f_network->get_num_lines();

  if( v_flow_limit_constraints.size() != num_lines ) {
    // this should only happen once
    assert( v_flow_limit_constraints.size() == 0 );
    v_flow_limit_constraints.resize( num_lines );
  }

  // Flow limit constraints

  auto num_nodes = f_network->get_num_nodes();

  for( int line_id = 0; line_id < num_lines; ++line_id ) {

    auto linear_function = new LinearFunction();
    double constant_term = 0;

    for( int node_id = 0; node_id < num_nodes; ++node_id ) {

      double coefficient = 0.0; // TODO Compute the Power Transfer
				// Distribution Factor Matrix"

      if( coefficient == 0.0 )
	continue;

      auto unit_blocks = f_network->get_node( node_id )->get_unit_blocks();

      for( int unit_block_id = 0; unit_block_id < unit_blocks.size();
	   ++unit_block_id ) {
	linear_function.add_variable
	  ( unit_blocks[unit_block_id]->get_power( f_time, coefficient ) );
      }

      constant_term -= coefficient * v_demand[node_id];

    } // for each node

    // Set the function of the constraint

    v_flow_limit_constraints[line_id].set_function(linear_function);

    // Set the left- and right-hand sides

    v_flow_limit_constraints[line_id].set_lhs
      ( v_minimum_power_flow[line_id] - constant_term );

    v_flow_limit_constraints[line_id].set_rhs
      ( v_maximum_power_flow[line_id] - constant_term );

  } // for each line

  add_static_constraint( v_flow_limit_constraints );
}

/*--------------------------------------------------------------------------*/
