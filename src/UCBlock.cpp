/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
 * \version 0.11
 *
 * \date 19 - 03 - 2019
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
#include "LinearFunction.h"
#include "NetworkBlock.h"
#include "NetworkNode.h"
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

void UCBlock::generate_abstract_constraints( Configuration *stcc ) {

  auto num_nodes = f_network.get_num_nodes();

  if( v_node_injection_constraints.size() != f_time_horizon ) {
    // this should only happen once
    assert( v_node_injection_constraints.size() == 0 );

    v_node_injection_constraints.resize
      ( boost::multi_array<FRowConstraint *, 2>::
        extent_gen()[f_time_horizon][num_nodes] );
  }

  // Node injection constraints.

  for( int t = 0; t < f_time_horizon; ++t ) {

    auto node_injection = v_network_blocks[t]->get_node_injection();

    for( int node_id = 0; node_id < num_nodes; ++node_id ) {

      auto linear_function = new LinearFunction();

      for( auto unit_block : f_network.get_node( node_id )->get_unit_blocks() )
        linear_function->add_variable( unit_block->get_power( t ), 1.0 );

      linear_function->add_variable( & node_injection[ node_id ], - 1.0 );

      v_node_injection_constraints[t][node_id]->set_both( 0.0 );
      v_node_injection_constraints[t][node_id]->set_function( linear_function );

    }
  }

  add_static_constraint( v_node_injection_constraints );
}

/*--------------------------------------------------------------------------*/
