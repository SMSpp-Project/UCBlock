/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCNetworkBlock class.
 *
 * \version 0.11
 *
 * \date 24 - 06 - 2019
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
 * Copyright &copy by Antonio Frangioni, and Ali Ghezelsoflu
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
#include "DCNetworkBlock.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DCNetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( DCNetworkBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration *stcc ) {

  if( f_NetworkData->get_number_lines() < 0 ) {
    throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                             "number of lines of DCNetworkBlock is not set"));
  }

  if( v_flow_limit_constraints.size() != f_NetworkData->get_number_lines() ) {
    // this should only happen once
    assert( v_flow_limit_constraints.size() == 0 );
    v_flow_limit_constraints.resize( f_NetworkData->get_number_lines() );
  }

  // Flow limit constraints

  // TODO

  for( Index line_id = 0; line_id < f_NetworkData->get_number_lines() ;
       ++line_id ) {

    auto min_power_flow = f_NetworkData-> get_min_power_flow()[line_id];
    auto max_power_flow = f_NetworkData-> get_max_power_flow()[line_id];
    auto linear_function = new LinearFunction();
    double constant_term = 0;

    for( Index node_id = 0; node_id < v_node_injection.size(); ++node_id ) {

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
        ( min_power_flow - constant_term );

    v_flow_limit_constraints[line_id].set_rhs
        ( max_power_flow - constant_term );

  } // for each line

  add_static_constraint( v_flow_limit_constraints );

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock -----*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
