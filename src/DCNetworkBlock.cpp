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

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv ) {

 int number_nodes = f_NetworkData->get_number_nodes();

 int number_lines = f_NetworkData->get_number_lines();

if (number_nodes > 1) {
 // the node injection variables
 if (v_node_injection.size() != number_nodes ){
  assert( v_node_injection.empty()); // this should only happen once
  v_node_injection.resize( number_nodes );
  int n = 0;
  for( auto & i : v_node_injection ) {
   i.set_type( ColVariable::kContinuous );
  }
  add_static_variable( v_node_injection, "S" );

 }
}

 if ( number_lines > 0 ) {

  // the power flow Variable

  if( v_power_flow.size() != number_lines ) {
   assert( v_power_flow.empty()); // this should only happen once
   v_power_flow.resize( number_lines );
   int n = 0;
   for( auto & i : v_power_flow ) {
    i.set_type( ColVariable::kContinuous );
   }
   add_static_variable( v_power_flow, "pf" );
  }
 }
 AR |= HasVar;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc ) {

 if( f_NetworkData->get_number_nodes() > 1 ) {

  if( f_NetworkData->get_number_lines() <= 0 ) {
   throw ( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                             "number of lines of DCNetworkBlock is not set" ));
  }

  // initial condition of number lines
  int number_lines = f_NetworkData->get_number_lines();

  // initial condition of vector StartLine
  std::vector< Index > StartLine = f_NetworkData->get_start_line();

  // initial condition of vector EndLine
  std::vector< Index > EndLine = f_NetworkData->get_end_line();

  // initial condition of minimum power flow
  std::vector< double > MinPowerFlow = f_NetworkData->get_min_power_flow();
  if( MinPowerFlow.size() == 1 ) {
   MinPowerFlow.resize( number_lines, MinPowerFlow[0] );

  }
  // initial condition of maximum power flow
  std::vector< double > MaxPowerFlow = f_NetworkData->get_max_power_flow();
  if( MaxPowerFlow.size() == 1 ) {
   MaxPowerFlow.resize( number_lines, MaxPowerFlow[0] );

  }
  // initial condition of Susceptance
  std::vector< double > Susceptance = f_NetworkData->get_susceptance();
  if( Susceptance.size() == 1 ) {
   Susceptance.resize( number_lines, Susceptance[0] );

  }

  //  Net Transfer Capacity (NTC) model
/*--------------------------------------------------------------------------*/
  // HVDC power flow limit
  if( v_HVDC_power_flow_limit_constraints.size() != f_NetworkData->get_number_lines()) {

   assert( v_HVDC_power_flow_limit_constraints.empty());
   v_HVDC_power_flow_limit_constraints.resize( f_NetworkData->get_number_lines());
  }

  for( Index line_id = 0; line_id < f_NetworkData->get_number_lines();
       ++line_id ) {

   if( Susceptance.empty() || Susceptance[line_id] == 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_power_flow[line_id], 1.0 );

    v_HVDC_power_flow_limit_constraints[line_id].set_lhs( MinPowerFlow[line_id] );
    v_HVDC_power_flow_limit_constraints[line_id].set_rhs( MaxPowerFlow[line_id] );
    v_HVDC_power_flow_limit_constraints[line_id].set_function( linear_function );

   }
  }
  add_static_constraint( v_HVDC_power_flow_limit_constraints, "HVDC_power_flow_limit" );

/*--------------------------------------------------------------------------*/

  // HVDC power flow and node injection constraints
  if( v_power_flow_injection_constraints.size() != f_NetworkData->get_number_lines()) {

   assert( v_power_flow_injection_constraints.empty());
   v_power_flow_injection_constraints.resize( f_NetworkData->get_number_lines());
  }

  Index end = 0;

  for( Index n = 0; n < f_NetworkData->get_number_nodes(); ++n ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_node_injection[n], -1.0 );

   end += n;

   for( Index line_id = 0; line_id < end; ++line_id ) {

    if( Susceptance.empty() || Susceptance[line_id] == 0 ) {

     if( StartLine[line_id] == n ) {

      linear_function->add_variable( &v_power_flow[line_id], -1.0 );
     } else {
      linear_function->add_variable( &v_power_flow[line_id], 1.0 );

     }
     v_power_flow_injection_constraints[ line_id ]
      .set_both( v_active_demand[ n ] );
     v_power_flow_injection_constraints[line_id].set_function( linear_function );

    }
   }
  }
  add_static_constraint( v_power_flow_injection_constraints, "HVDC_power_flow_injection" );

/*--------------------------------------------------------------------------*/

// TODO implementation of AC and AC-HVDC lines is not ready

  // AC power flow limit
  if( v_AC_power_flow_limit_constraints.size() != f_NetworkData->get_number_lines()) {
   // this should only happen once
   assert( v_AC_power_flow_limit_constraints.empty());
   v_AC_power_flow_limit_constraints.resize( f_NetworkData->get_number_lines());
  }

  // Flow limit constraints

  // TODO
  for( Index line_id = 0; line_id < f_NetworkData->get_number_lines();
       ++line_id ) {
   if( Susceptance[line_id] != 0 ) { //  just AC lines model


    auto min_power_flow = f_NetworkData->get_min_power_flow()[line_id];
    auto max_power_flow = f_NetworkData->get_max_power_flow()[line_id];
    auto linear_function = new LinearFunction();
    double constant_term = 0;

    for( Index node_id = 0; node_id < v_node_injection.size(); ++node_id ) {

     double coefficient = 0.0; // TODO Compute the Power Transfer
     // Distribution Factor Matrix

     if( coefficient == 0.0 )
      continue;

     linear_function->add_variable
             ( &v_node_injection[node_id], coefficient );

     constant_term -= coefficient * v_active_demand[node_id];

    } // for each node

    // Set the function of the constraint

    v_AC_power_flow_limit_constraints[line_id].set_function( linear_function );

    // Set the left- and right-hand sides

    v_AC_power_flow_limit_constraints[line_id].set_lhs
            ( min_power_flow - constant_term );

    v_AC_power_flow_limit_constraints[line_id].set_rhs
            ( max_power_flow - constant_term );

   } // for each line

  }
  add_static_constraint( v_AC_power_flow_limit_constraints, "AC_power_low_limits" );


  // AC-HVDC power flow limit

  if( v_AC_HVDC_power_flow_constraints.size() != f_NetworkData->get_number_lines()) {
   // this should only happen once
   assert( v_AC_HVDC_power_flow_constraints.empty());
   v_AC_HVDC_power_flow_constraints.resize( f_NetworkData->get_number_lines());
  }
  // TODO
  // TODO
  // TODO

  add_static_constraint( v_AC_HVDC_power_flow_constraints, "AC/HVDC_power_flow_limits" );

 }
 AR |= HasCst;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock -----*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
