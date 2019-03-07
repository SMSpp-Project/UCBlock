/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
 * \version 0.11
 *
 * \date 07 - 03 - 2019
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

#include <iostream>
#include <vector>
#include "FRowConstraint.h"
#include "NetWorkBlock.h"
#include "UCBlock.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_variables( Configuration *stvv ) {

  if( ! f_network ) {
    throw( std::logic_error( "UCBlock::generate_abstract_variables: "
			     "f_network of UCBlock is not set" ) );
  }

  auto num_nodes = f_network->get_num_nodes();

  if( v_node_injection.size() != num_nodes ) {
    assert( v_node_injection.size() == 0 ); // this should only happen once
    v_node_injection.resize( num_nodes );
    add_static_variable( v_node_injection );
  }
}

/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_constraints( Configuration *stcc ) {

  if( ! f_network ) {
    throw( std::logic_error( "UCBlock::generate_abstract_constraints: "
			     "f_network of UCBlock is not set" ) );
  }

  auto num_nodes = f_network->get_num_nodes();

  if( v_node_injection_constraints.size() != f_time_horizon ) {
    // this should only happen once
    assert( v_node_injection_constraints.size() == 0 );

    v_node_injection_constraints.resize
      ( boost::multi_array<FRowConstraint *, 2>::
	extent_gen()[f_time_horizon][num_nodes] );
  }

  // Node injection constraints.

  int node_id = 0;
  for( int t = 0; t < f_time_horizon; ++t ) {
    for( auto node : f_network->get_nodes() ) {

      auto linear_function = new LinearFunction();

      for( auto unit_block : node->get_unit_blocks() )
	linear_function.add_variable( unit_block->get_power(t), 1.0);

      linear_function.add_variable
	( v_network_blocks[t].get_node_injection( node_id ), - 1.0);

      v_node_injection_constraints[t][node_id].set_both(0.0);
      v_node_injection_constraints[t][node_id].set_function(linear_function);

      ++node_id;
    }

    add_static_constraint( v_node_injection_constraints[t] );
  }
}

/*--------------------------------------------------------------------------*/
