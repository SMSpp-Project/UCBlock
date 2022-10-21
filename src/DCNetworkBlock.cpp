/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCNetworkBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                   Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
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
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

DCNetworkBlock::~DCNetworkBlock()
{
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

 for( auto & constraint : v_power_flow_relax_abs_1 )
  constraint.clear();

 for( auto & constraint : v_power_flow_relax_abs_2 )
  constraint.clear();

 objective.clear();

 // Delete the NetworkData if it is local.
 if( f_local_NetworkData )
  delete f_NetworkData;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
 auto NumberNodes = group.getDim( "NumberNodes" );

 if( ! NumberNodes.isNull() ) {
  // Since the dimension "NumberNodes" has been provided, it means that a
  // NetworkData has been provided. Thus, the NetworkData is deserialized and
  // it is marked as being local.
  delete f_NetworkData;
  f_NetworkData = new NetworkData();
  f_NetworkData->deserialize( group );
  f_local_NetworkData = true;
 }

 NetworkBlock::deserialize( group );
 }

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )
  return; // variables have already been generated

 const auto number_nodes = get_number_nodes();
 const auto number_lines = get_number_lines();

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

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )
  return; // constraints have already been generated

 const auto number_nodes = get_number_nodes();

 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 ) {
  throw ( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                            "number of lines of DCNetworkBlock is not set" ) );
 }

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 const auto lines_type = f_NetworkData->get_lines_type();

 if( lines_type == kHVDC ) {   // HVDC power flow limit
  //  Net Transfer Capacity (NTC) model

  // Flow limit constraints

  if( v_HVDC_power_flow_limit_constraints.size() != number_lines ) {
   assert( v_HVDC_power_flow_limit_constraints.empty() );
   v_HVDC_power_flow_limit_constraints.resize( number_lines );
  }

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {

   const auto kappa = get_kappa( line_id );
   v_HVDC_power_flow_limit_constraints[ line_id ].set_lhs
    ( kappa * get_min_power_flow( line_id ) );
   v_HVDC_power_flow_limit_constraints[ line_id ].set_rhs
    ( kappa * get_max_power_flow( line_id ) );
   v_HVDC_power_flow_limit_constraints[ line_id ].set_variable
    ( & v_power_flow[ line_id ] );
  }

  add_static_constraint( v_HVDC_power_flow_limit_constraints ,
                         "HVDC_power_flow_limit" );

  // HVDC power flow and node injection constraints

  if( v_power_flow_injection_constraints.size() != number_nodes ) {
   assert( v_power_flow_injection_constraints.empty() );
   v_power_flow_injection_constraints.resize( number_nodes );
  }

  for( Index n = 0 ; n < number_nodes ; ++n ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( & v_node_injection[ n ] , -1.0 );

   for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
    if( start_line[ line_id ] == n )
     linear_function->add_variable( & v_power_flow[ line_id ] , 1.0 );
    if( end_line[ line_id ] == n )
     linear_function->add_variable( & v_power_flow[ line_id ] , -1.0 );
   }
   v_power_flow_injection_constraints[ n ].set_both( -v_active_demand[ n ] );
   v_power_flow_injection_constraints[ n ].set_function( linear_function );
  }

  add_static_constraint( v_power_flow_injection_constraints ,
                         "HVDC_power_flow_injection" );

/*--------------------------------------------------------------------------*/
  if( ! f_NetworkData->get_network_cost().empty() ) {

   // Auxiliary constraints for relaxing the absolute value

   // F_l <= V_l

   if( v_power_flow_relax_abs_1.size() != number_lines ) {
    assert( v_power_flow_relax_abs_1.empty() );
    v_power_flow_relax_abs_1.resize( number_lines );
   }

   for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_power_flow[ line_id ] , -1.0 );
    linear_function->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
    v_power_flow_relax_abs_1[ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs_1[ line_id ].set_rhs( Inf< double >() );

    v_power_flow_relax_abs_1[ line_id ].set_function( linear_function );
   }

   add_static_constraint( v_power_flow_relax_abs_1 , "power_flow_relax_abs_1" );

   // - F_l <= V_l

   if( v_power_flow_relax_abs_2.size() != number_lines ) {
    assert( v_power_flow_relax_abs_2.empty() );
    v_power_flow_relax_abs_2.resize( number_lines );
   }

   for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {

    auto linear_f = new LinearFunction();

    linear_f->add_variable( &v_power_flow[ line_id ] , 1.0 );
    linear_f->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
    v_power_flow_relax_abs_2[ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs_2[ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs_2[ line_id ].set_function( linear_f );
   }

   add_static_constraint( v_power_flow_relax_abs_2 , "power_flow_relax_abs_2" );
  } // end( cost not empty )
 } // end( HVDC_Lines constraints )

/*--------------------------------------------------------------------------*/
// TODO implementation of AC and AC-HVDC lines is not ready

 else if( lines_type == kAC ) {    // AC power flow limit
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
            ( get_min_power_flow(line_id) - constant_term );

    v_AC_power_flow_limit_constraints[line_id].set_rhs
            ( get_max_power_flow(line_id) - constant_term );

  }
  add_static_constraint( v_AC_power_flow_limit_constraints, "AC_power_low_limits" );
*/
 } // end AC_Lines constraints
/*--------------------------------------------------------------------------*/

 else if( lines_type == kAC_HVDC ) { // AC-HVDC power flow limit


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

  if( !f_NetworkData->get_network_cost().empty() ) {

   auto linear_function = new LinearFunction();
   for( Index l = 0; l < f_NetworkData->get_number_lines(); ++l ) {
    linear_function->add_variable( &v_auxiliary_variable[ l ] ,
                                   f_NetworkData->get_network_cost()[l] ,
                                   0.0 );
    objective.set_function( linear_function );
    objective.set_sense( Objective::eMin );

   }

  } else { // empty objective function
   auto linear_function = new LinearFunction();
   objective.set_function( linear_function );
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
/*----------------- METHODS FOR CHECKING THE DCNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/

bool DCNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc ) {

 // Retrieve the tolerance and the type of violation.

 double tolerance = 0;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ & tolerance , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tolerance = tc->f_value;
   return true;
  }
  if( auto tc = dynamic_cast< SimpleConfiguration<
      std::pair< double , int > > * >( c ) ) {
   tolerance = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return true;
  }
  return false;
 };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 auto is_feasible = [ tolerance , rel_viol ]( auto & constraints ) {
  return RowConstraint::is_feasible( constraints , tolerance , rel_viol );
 };

 return NetworkBlock::is_feasible( useabstract )
  && ColVariable::is_feasible( v_node_injection , tolerance )
  && ColVariable::is_feasible( v_power_flow , tolerance )
  && ColVariable::is_feasible( v_auxiliary_variable , tolerance )
  && is_feasible( v_AC_power_flow_limit_constraints )
  && is_feasible( v_HVDC_power_flow_limit_constraints )
  && is_feasible( v_AC_HVDC_power_flow_limit_constraints )
  && is_feasible( v_power_flow_injection_constraints )
  && is_feasible( v_AC_HVDC_power_flow_constraints )
  && is_feasible( v_power_flow_relax_abs_1 )
  && is_feasible( v_power_flow_relax_abs_2 );

} // end( DCNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_active_demand
( std::vector< double >::const_iterator values , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 if( v_active_demand.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return cst == 0; } ) ) {
   return;
  }

  v_active_demand.assign( get_number_nodes() , 0 );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_active_demand.size() )
   throw ( std::invalid_argument( "DCNetworkBlock::set_active_demand: "
                                  "invalid value in subset" ) );
  auto demand = *( values++ );
  if( v_active_demand[ i ] != demand ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_active_demand[ i ] = demand;
  }
 }
 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) && not_dry_run( issueAMod ) &&
     constraints_generated() ) {

  // Change the abstract representation

  switch( f_NetworkData->get_lines_type() ) {
   case( kHVDC ): {
    for( auto i : subset )
     v_power_flow_injection_constraints[ i ].set_both( -v_active_demand[ i ] ,
                                                       issueAMod );
    break;
   }
   case( kAC ):
    // TODO
    break;
   case( kAC_HVDC ):
    // TODO
    break;
   default: break;
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered )
   std::sort( subset.begin(), subset.end() );

  Block::add_Modification( std::make_shared< NetworkBlockSbstMod >
                           ( this , NetworkBlockMod::eSetActD ,
                             std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_active_demand
( std::vector< double >::const_iterator values , Block::Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second , get_number_nodes() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_active_demand.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return ( cst == 0 ); } ) ) {
   return;
  }

  v_active_demand.assign( get_number_nodes() , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_active_demand.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_active_demand.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   switch( f_NetworkData->get_lines_type() ) {
    case( kHVDC ): {
     for( Index i = rng.first ; i < rng.second ; ++i )
      v_power_flow_injection_constraints[ i ].set_both( -v_active_demand[ i ] ,
                                                        issueAMod );
     break;
    }
    case( kAC ):
     // TODO
     break;
    case( kAC_HVDC ):
     // TODO
     break;
    default: break;
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< NetworkBlockRngdMod >
                           ( this , NetworkBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa
( std::vector< double >::const_iterator values , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 if( v_kappa.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return cst == 1; } ) )
   return;

  v_kappa.assign( get_number_lines() , 1 );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_kappa.size() )
   throw ( std::invalid_argument( "DCNetworkBlock::set_kappa: "
                                  "invalid value in subset: " +
                                  std::to_string( i ) + "." ) );
  const auto kappa = *( values++ );
  if( v_kappa[ i ] != kappa ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_kappa[ i ] = kappa;
  }
 }
 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) && not_dry_run( issueAMod ) &&
     constraints_generated() ) {

  // Change the abstract representation

  switch( f_NetworkData->get_lines_type() ) {
   case( kHVDC ): {
    for( auto i : subset ) {
     v_HVDC_power_flow_limit_constraints[ i ].set_lhs
      ( v_kappa[ i ] * get_min_power_flow( i ) , issueAMod );

     v_HVDC_power_flow_limit_constraints[ i ].set_rhs
      ( v_kappa[ i ] * get_max_power_flow( i ) , issueAMod );
    }
    break;
   }
   case( kAC ):
    // TODO
    break;
   case( kAC_HVDC ):
    // TODO
    break;
   default: break;
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< DCNetworkBlockSbstMod >
                           ( this , DCNetworkBlockMod::eSetKappa ,
                             std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa
( std::vector< double >::const_iterator values , Block::Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second , get_number_lines() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_kappa.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return ( cst == 1 ); } ) ) {
   return;
  }

  v_kappa.assign( get_number_lines() , 1 );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_kappa.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_kappa.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   switch( f_NetworkData->get_lines_type() ) {
    case( kHVDC ): {
     for( Index i = rng.first ; i < rng.second ; ++i ) {
      v_HVDC_power_flow_limit_constraints[ i ].set_lhs
       ( v_kappa[ i ] * get_min_power_flow( i ) , issueAMod );

      v_HVDC_power_flow_limit_constraints[ i ].set_rhs
       ( v_kappa[ i ] * get_max_power_flow( i ) , issueAMod );
     }
     break;
    }
    case( kAC ):
     // TODO
     break;
    case( kAC_HVDC ):
     // TODO
     break;
    default: break;
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< DCNetworkBlockRngdMod >
                           ( this , DCNetworkBlockMod::eSetKappa , rng ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
