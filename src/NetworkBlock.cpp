/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 24 - 04 - 2019
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



/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */



//void NetworkBlock::load( )
//{


//}


/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( netCDF::NcGroup & group , Block * father ) {



// read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -


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
				// Distribution Factor Matrix"

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

void NetworkBlock::serialize( netCDF::NcGroup & group ) const
{

}    // end( NetworkBlock::serialize )
/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/