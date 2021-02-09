/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCNetworkBlock class.
 *
 * \version 0.11
 *
 * \date 11 - 10 - 2020
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
#include "OneVarConstraint.h"
#include "FRealObjective.h"


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
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

DCNetworkBlock::~DCNetworkBlock() {
 for( auto & constraint : v_AC_power_flow_limit_constraints )
  constraint.clear();

 for( auto & constraint : v_HVDC_power_flow_limit_constraints )
  constraint.clear();

 for( auto & constraint : v_AC_HVDC_power_flow_limit_constraints )
  constraint.clear();

 for( auto & constraint : v_power_flow_injection_constraints )
  constraint.clear();

 for( auto & constraint : v_AC_HVDC_power_flow_constraints )
  constraint.clear();

 for( auto & constraint : v_power_flow_auxiliary_variable_one_constraints)
  constraint.clear();

 for( auto & constraint : v_power_flow_auxiliary_variable_two_constraints)
  constraint.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 int number_nodes = f_NetworkData->get_number_nodes();

 int number_lines = f_NetworkData->get_number_lines();

 if( number_nodes > 1 ) {
  // the node injection variables
  v_node_injection.resize( number_nodes );
  for( auto & var : v_node_injection )
   var.set_type( ColVariable::kContinuous );
  add_static_variable( v_node_injection, "S" );
 }

 if( number_lines > 0 ) {
  // the power flow Variable
  v_power_flow.resize( number_lines );
  for( auto & var : v_power_flow )
   var.set_type( ColVariable::kContinuous );
  add_static_variable( v_power_flow, "F_power_flow" );

  if( !f_NetworkData->get_network_cost().empty() ) {
   // the auxiliary Variable
   v_auxiliary_variable.resize( number_lines );
   for( auto & var : v_auxiliary_variable )
    var.set_type( ColVariable::kContinuous );
   add_static_variable( v_auxiliary_variable, "V_auxiliary" );
  }
 }
 set_variables_generated();

}

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc ) {

 if( constraints_generated() )
  return; // constraints have already been generated

 if( f_NetworkData->get_number_nodes() > 1 ) {

  if( f_NetworkData->get_number_lines() <= 0 ) {
   throw ( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                             "number of lines of DCNetworkBlock is not set" ) );
  }

  // initial condition of number lines
  auto number_lines = f_NetworkData->get_number_lines();

  // initial condition of vector StartLine
  auto & StartLine = f_NetworkData->get_start_line();

  // initial condition of vector EndLine
  auto & EndLine = f_NetworkData->get_end_line();

  // initial condition of minimum power flow
  std::vector< double > MinPowerFlow = f_NetworkData->get_min_power_flow();
  if( MinPowerFlow.size() == 1 ) {
   MinPowerFlow.resize( number_lines , MinPowerFlow[ 0 ] );
  }

  // initial condition of maximum power flow
  std::vector< double > MaxPowerFlow = f_NetworkData->get_max_power_flow();
  if( MaxPowerFlow.size() == 1 ) {
   MaxPowerFlow.resize( number_lines , MaxPowerFlow[ 0 ] );
  }

  // initial condition of Susceptance
  std::vector< double > Susceptance = f_NetworkData->get_susceptance();
  if( Susceptance.size() == 1 ) {
   Susceptance.resize( number_lines , Susceptance[ 0 ] );
  }

  //  Net Transfer Capacity (NTC) model
/*--------------------------------------------------------------------------*/
// Initial check on network

  auto lines_type = f_NetworkData->get_lines_type();

/*--------------------------------------------------------------------------*/

  if( lines_type == kHVDC ) {   // HVDC power flow limit

   if( v_HVDC_power_flow_limit_constraints.size() !=
       f_NetworkData->get_number_lines() ) {

    assert( v_HVDC_power_flow_limit_constraints.empty() );
    v_HVDC_power_flow_limit_constraints.resize
     ( f_NetworkData->get_number_lines() );
   }

   for( Index line_id = 0; line_id < f_NetworkData->get_number_lines();
        ++line_id ) {

    v_HVDC_power_flow_limit_constraints[line_id].set_lhs(MinPowerFlow[line_id]);
    v_HVDC_power_flow_limit_constraints[line_id].set_rhs(MaxPowerFlow[line_id]);
    v_HVDC_power_flow_limit_constraints[line_id].set_variable(&v_power_flow[line_id]);
   }

   add_static_constraint( v_HVDC_power_flow_limit_constraints ,
                          "HVDC_power_flow_limit" );

/*--------------------------------------------------------------------------*/

   // HVDC power flow and node injection constraints

   if( v_power_flow_injection_constraints.size() !=
       f_NetworkData->get_number_nodes() ) {

    assert( v_power_flow_injection_constraints.empty() );
    v_power_flow_injection_constraints.resize
     ( f_NetworkData->get_number_nodes() );
   }

   for( Index n = 0 ; n < f_NetworkData->get_number_nodes() ; ++n ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_node_injection[n], -1.0 );

    for( Index line_id = 0; line_id < number_lines; ++line_id ) {
     if( StartLine[line_id] == n )
      linear_function->add_variable( &v_power_flow[line_id], 1.0 );
     if( EndLine[line_id] == n )
      linear_function->add_variable( &v_power_flow[line_id], -1.0 );
    }
    v_power_flow_injection_constraints[n].set_both( -v_active_demand[n] );
    v_power_flow_injection_constraints[n].set_function( linear_function );
   }

   add_static_constraint( v_power_flow_injection_constraints ,
                          "HVDC_power_flow_injection" );

/*--------------------------------------------------------------------------*/
   if( !f_NetworkData->get_network_cost().empty() ) {

    if( v_power_flow_auxiliary_variable_one_constraints.size() !=
        f_NetworkData->get_number_lines()) {

     assert( v_power_flow_auxiliary_variable_one_constraints.empty());
     v_power_flow_auxiliary_variable_one_constraints.resize
             ( f_NetworkData->get_number_lines());
    }
    auto linear_function = new LinearFunction();

    for( Index line_id = 0; line_id < f_NetworkData->get_number_lines();
         ++line_id ) {
     linear_function->add_variable( &v_power_flow[line_id], -1.0 );
     linear_function->add_variable( &v_auxiliary_variable[line_id], 1.0 );
     v_power_flow_auxiliary_variable_one_constraints[line_id].set_lhs( 0.0 );
     v_power_flow_auxiliary_variable_one_constraints[line_id].set_function( linear_function );
    }

    add_static_constraint( v_power_flow_auxiliary_variable_one_constraints,
                           "power_flow_auxiliary_variable_one" );

/*--------------------------------------------------------------------------*/
    if( v_power_flow_auxiliary_variable_two_constraints.size() !=
        f_NetworkData->get_number_lines()) {

     assert( v_power_flow_auxiliary_variable_two_constraints.empty());
     v_power_flow_auxiliary_variable_two_constraints.resize
             ( f_NetworkData->get_number_lines());
    }
    auto linear_f = new LinearFunction();

    for( Index line_id = 0; line_id < f_NetworkData->get_number_lines();
         ++line_id ) {
     linear_function->add_variable( &v_power_flow[line_id], 1.0 );
     linear_function->add_variable( &v_auxiliary_variable[line_id], 1.0 );
     v_power_flow_auxiliary_variable_two_constraints[line_id].set_lhs( 0.0 );
     v_power_flow_auxiliary_variable_two_constraints[line_id].set_function( linear_f );
    }

    add_static_constraint( v_power_flow_auxiliary_variable_two_constraints,
                           "power_flow_auxiliary_variable_two" );

   }
  } // end HVDC_Lines constraints
/*--------------------------------------------------------------------------*/
// TODO implementation of AC and AC-HVDC lines is not ready

  if( lines_type == kAC ) {    // AC power flow limit
/*
  if( v_AC_power_flow_limit_constraints.size() != f_NetworkData->get_number_lines()) {
   // this should only happen once
   assert( v_AC_power_flow_limit_constraints.empty());
   v_AC_power_flow_limit_constraints.resize( f_NetworkData->get_number_lines());
  }

  // Flow limit constraints

  for( Index line_id = 0; line_id < f_NetworkData->get_number_lines();
       ++line_id ) {


    auto linear_function = new LinearFunction();
    double constant_term = 0;

    for( Index node_id = 0; node_id < f_NetworkData->get_number_nodes(); ++node_id ) {

     double coefficient = 0.0;
     // Distribution Factor Matrix

     linear_function->add_variable
             ( &v_node_injection[node_id], coefficient );

     constant_term -= coefficient * v_active_demand[node_id];

    } // for each node


    // Set the function of the constraint

    v_AC_power_flow_limit_constraints[line_id].set_function( linear_function );

    // Set the left- and right-hand sides

    v_AC_power_flow_limit_constraints[line_id].set_lhs
            ( MinPowerFlow[line_id] - constant_term );

    v_AC_power_flow_limit_constraints[line_id].set_rhs
            ( MaxPowerFlow[line_id] - constant_term );

  }
  add_static_constraint( v_AC_power_flow_limit_constraints, "AC_power_low_limits" );
*/
  } // end AC_Lines constraints
/*--------------------------------------------------------------------------*/

  if( lines_type == kAC_HVDC ) { // AC-HVDC power flow limit


/*
  if( v_AC_HVDC_power_flow_constraints.size() != f_NetworkData->get_number_lines()) {
   // this should only happen once
   assert( v_AC_HVDC_power_flow_constraints.empty());
   v_AC_HVDC_power_flow_constraints.resize( f_NetworkData->get_number_lines());
  }
  // TODO
  // TODO
  // TODO

  add_static_constraint( v_AC_HVDC_power_flow_constraints, "AC/HVDC_power_flow_limits" );
*/
  } // end AC-HVDC constraints
 }
 set_constraints_generated();
}
/*--------------------------------------------------------------------------*/
void DCNetworkBlock::generate_objective( Configuration * objc ) {

// Initial check on network

 auto lines_type = f_NetworkData->get_lines_type();

/*--------------------------------------------------------------------------*/

 if( lines_type == kHVDC ) {   // HVDC power flow limit

  if( objective_generated())
   return; // Objective has already been generated

  if( get_objective() != nullptr )  // an objective is there already
   return;                         // cowardly (and silently) return

  if( !f_NetworkData->get_network_cost().empty() ) { // empty objective function

   auto linear_function = new LinearFunction();
   objective.set_function( linear_function );
  } else {

   auto linear_function = new LinearFunction();
   for( Index l = 0; l < f_NetworkData->get_number_lines(); ++l ) {
    linear_function->add_variable( &v_auxiliary_variable[ l ] ,
                                  f_NetworkData->get_network_cost()[l] ,
                                  0.0 );
    objective.set_function( linear_function );
    objective.set_sense( Objective::eMin );

   }
  }
  // Set Block objective
  this->set_objective( &objective );
 }
// TODO The implementation of objective function for AC and AC-HVDC lines is
//  not ready

 if( lines_type == kAC ) {    // AC power flow limit

  //TODO
 }

 if( lines_type == kAC_HVDC ) { // AC-HVDC power flow limit

  //TODO

 }


 set_objective_generated();

}  // end( DCNetworkBlock::generate_objective )
/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void
NetworkBlock::set_active_demand( std::vector< double >::const_iterator values,
                                 Block::Subset && subset,
                                 const bool ordered,
                                 c_ModParam issuePMod,
                                 c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_active_demand.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = * std::max_element( std::begin( subset ),
                                        std::end( subset ) );
  v_active_demand.assign( max_index, 0 );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_active_demand.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( v_active_demand[ i ] != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   v_active_demand[ i ] = *( values++ );
  }

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   // FIXME: This is correct only for BusNetworkBlock

  }
 }

}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void
DCNetworkBlock::set_active_demand( std::vector< double >::const_iterator values,
                                   Block::Subset && subset,
                                   const bool ordered,
                                   c_ModParam issuePMod,
                                   c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_active_demand.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) { return cst == 0; } ) ) {
   return;
  }

  Index max_index = * std::max_element( std::begin( subset ),
                                        std::end( subset ) );
  v_active_demand.assign( max_index, 0 );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_active_demand.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( v_active_demand[ i ] != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   v_active_demand[ i ] = *( values++ );
  }

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   auto s = f_NetworkData->get_susceptance();

   if( s.empty() || std::all_of( s.begin(), s.end(),
                                 []( double i ) { return i == 0; } ) ) {
    // NTC Model
    // for( auto i : subset ) {
    //  v_power_flow_injection_constraints[ i ]
    //   .set_both( v_active_demand[ i ], issueAMod );
    // }

   } else if( std::all_of( s.begin(), s.end(),
                           []( double i ) { return i != 0; } ) ) {
    // Just DC lines
    // TODO
   } else {
    // HVDC/DC
    // TODO
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< NetworkBlockSbstMod >( this,
                                            NetworkBlockMod::eSetActD,
                                            std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

void
DCNetworkBlock::set_active_demand( std::vector< double >::const_iterator values,
                                   Block::Range rng,
                                   c_ModParam issuePMod,
                                   c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, get_number_nodes() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_active_demand.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = rng.second;
  v_active_demand.assign( max_index, 0 );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_active_demand.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_active_demand.begin() + rng.first );

  if( constraints_generated() ) {
   // Change the abstract representation

   auto s = f_NetworkData->get_susceptance();

   if( s.empty() || std::all_of( s.begin(), s.end(),
                                 []( double i ) { return i == 0; } ) ) {
    // NTC Model
    // for( Index i = rng.first; i < rng.second; ++i ) {
    //  v_power_flow_injection_constraints[ i ]
    //   .set_both( v_active_demand[ i ], issueAMod );
    // }

   } else if( std::all_of( s.begin(), s.end(),
                           []( double i ) { return i != 0; } ) ) {
    // Just DC lines
    // TODO
   } else {
    // HVDC/DC
    // TODO
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< NetworkBlockRngdMod >( this,
                                            NetworkBlockMod::eSetActD,
                                            rng ),
   Observer::par2chnl( issuePMod ) );
 }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock -----*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
