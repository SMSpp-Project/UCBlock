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

DCNetworkBlock::~DCNetworkBlock() {

 clear_constraints( v_AC_power_flow_limit_constraints );
 clear_constraints( v_AC_HVDC_power_flow_limit_constraints );
 clear_constraints( v_power_flow_injection_constraints );
 clear_constraints( v_AC_HVDC_power_flow_constraints );
 clear_constraints( v_power_flow_aux_var_one_constraints );
 clear_constraints( v_power_flow_aux_var_two_constraints );

 clear_constraints( v_HVDC_power_flow_limit_constraints );

 objective.clear();

 // Delete the DCNetworkData if it is local.
 if( f_local_NetworkData )
  delete f_NetworkData;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::DCNetworkData::deserialize(
 const netCDF::NcGroup & group ) {

 NetworkBlock::NetworkData::deserialize( group );

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberLines" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "StartLine" , "EndLine" ,
                                                     "MinPowerFlow" ,
                                                     "MaxPowerFlow" ,
                                                     "Susceptance" ,
                                                     "NetworkCost" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Optional variables

 if( !::deserialize_dim( group , "NumberNodes" , f_number_nodes , true ) )
  f_number_nodes = 1;

 if( f_number_nodes > 1 ) {

  ::deserialize( group , "MinPowerFlow" , f_number_lines , v_min_power_flow ,
                 true , true );

  ::deserialize( group , "MaxPowerFlow" , f_number_lines , v_max_power_flow ,
                 true , true );

  ::deserialize( group , "Susceptance" , f_number_lines , v_susceptance ,
                 true , true );

  ::deserialize( group , "NetworkCost" , f_number_lines , v_network_cost ,
                 true , true );
 }
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::deserialize( const netCDF::NcGroup & group ) {

 NetworkBlock::deserialize( group );

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "ConstTerm" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Optional variables

 Index NumberNodes;
 if( ::deserialize_dim( group , "NumberNodes" , NumberNodes , true ) ) {
  // Since the dimension "NumberNodes" has been provided, it means that a
  // DCNetworkData has been provided. Thus, the DCNetworkData is deserialized,
  // and it is marked as being local.
  delete f_NetworkData;
  f_NetworkData = new DCNetworkData();
  f_NetworkData->deserialize( group );
  f_local_NetworkData = true;
  // A DCNetworkData has been provided. So, the size of the given vector of
  // active demand must be equal to the number of nodes.
  ::deserialize( group , "ActiveDemand" , NumberNodes ,
                 v_active_demand , true );
 } else {
  // A DCNetworkData has not been provided. However, the active demand may still
  // have been provided.
  auto ActiveDemand = group.getVar( "ActiveDemand" );

  if( !ActiveDemand.isNull() ) {
   // The active demand has indeed been provided.

   if( ActiveDemand.getDimCount() != 1 )
    // The active demand must be a one-dimensional array.
    throw ( std::invalid_argument(
     "DCNetworkBlock::deserialize(): ActiveDemand should have one dimension, "
     "but it has " + std::to_string( ActiveDemand.getDimCount() ) ) );

   // Retrieve the number of nodes from the size of the given netCDF variable.
   const auto number_nodes = ActiveDemand.getDim( 0 ).getSize();

   // Resize the vector of active demand.
   v_active_demand.resize( number_nodes );

   // Retrieve the active demand from the netCDF variable.
   ActiveDemand.getVar( v_active_demand.data() );
  }
 }
}

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::DCNetworkData::serialize( netCDF::NcGroup & group ) const {

 NetworkBlock::NetworkData::serialize( group );

 if( f_number_nodes > 1 ) {
  auto NumberLines = group.getDim( "NumberLines" );

  ::serialize( group , "MinPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_min_power_flow );

  ::serialize( group , "MaxPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_max_power_flow );

  ::serialize( group , "Susceptance" , netCDF::NcDouble() , NumberLines ,
               v_susceptance );

  ::serialize( group , "NetworkCost" , netCDF::NcDouble() , NumberLines ,
               v_network_cost );
 }
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::serialize( netCDF::NcGroup & group ) const {

 NetworkBlock::serialize( group );

 if( auto network_data = get_NetworkData() )
  // If a DCNetworkData is present, serialize it.
  network_data->serialize( group );

 if( !v_active_demand.empty() ) {
  // This DCNetworkBlock has active demand, so it is serialized.

  auto NumberNodes = group.getDim( "NumberNodes" );

  if( NumberNodes.isNull() ) {
   /* The dimension "NumberNodes" is not present in the group (which means
    * that a DCNetworkData is not present). However, the number of nodes can
    * still be obtained from the size of the active demand vector. Notice that
    * the name "NumberNodes" is not used for this new dimension, because it
    * would indicate that a DCNetworkData is present (which is not the
    * case). Therefore, we create an alternative dimension in order to be able
    * to serialize the active demand. */
   NumberNodes = group.addDim( "__NumberNodes__" , v_active_demand.size() );
  }

  // Finally, serialize the active demand.
  ::serialize( group , "ActiveDemand" , netCDF::NcDouble() ,
               NumberNodes , v_active_demand );
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv ) {
 if( variables_generated() )
  return; // variables have already been generated

 Index number_nodes = get_number_nodes();

 Index number_lines = get_number_lines();

 if( number_nodes > 1 ) {
  // the node injection variables
  v_node_injection.resize( number_nodes );
  for( auto & var : v_node_injection )
   var.set_type( ColVariable::kContinuous );
  add_static_variable( v_node_injection , "S" );
 }

 if( number_lines > 0 ) {
  // the power flow Variable
  v_power_flow.resize( number_lines );
  for( auto & var : v_power_flow )
   var.set_type( ColVariable::kContinuous );
  add_static_variable( v_power_flow , "F_power_flow" );

  if( !f_NetworkData->get_network_cost().empty() ) {
   // the auxiliary Variable
   v_auxiliary_variable.resize( number_lines );
   for( auto & var : v_auxiliary_variable )
    var.set_type( ColVariable::kContinuous );
   add_static_variable( v_auxiliary_variable , "V_auxiliary" );
  }
 }

 set_variables_generated();
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc ) {
 if( constraints_generated() )
  return; // constraints have already been generated

 if( get_number_nodes() > 1 ) {

  if( get_number_lines() <= 0 ) {
   throw ( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                             "number of lines of DCNetworkBlock is not set" ) );
  }

  // initial condition of number lines
  auto number_lines = get_number_lines();

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

   if( v_HVDC_power_flow_limit_constraints.size() != get_number_lines() ) {

    assert( v_HVDC_power_flow_limit_constraints.empty() );
    v_HVDC_power_flow_limit_constraints.resize( get_number_lines() );
   }

   for( Index line_id = 0 ; line_id < get_number_lines() ; ++line_id ) {

    v_HVDC_power_flow_limit_constraints[ line_id ].set_lhs(
     MinPowerFlow[ line_id ] );
    v_HVDC_power_flow_limit_constraints[ line_id ].set_rhs(
     MaxPowerFlow[ line_id ] );
    v_HVDC_power_flow_limit_constraints[ line_id ].set_variable(
     &v_power_flow[ line_id ] );
   }

   add_static_constraint( v_HVDC_power_flow_limit_constraints ,
                          "HVDC_power_flow_limit" );

/*--------------------------------------------------------------------------*/

   // HVDC power flow and node injection constraints

   if( v_power_flow_injection_constraints.size() != get_number_nodes() ) {

    assert( v_power_flow_injection_constraints.empty() );
    v_power_flow_injection_constraints.resize( get_number_nodes() );
   }

   for( Index n = 0 ; n < get_number_nodes() ; ++n ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_node_injection[ n ] , -1.0 );

    for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
     if( StartLine[ line_id ] == n )
      linear_function->add_variable( &v_power_flow[ line_id ] , 1.0 );
     if( EndLine[ line_id ] == n )
      linear_function->add_variable( &v_power_flow[ line_id ] , -1.0 );
    }
    v_power_flow_injection_constraints[ n ].set_both( -v_active_demand[ n ] );
    v_power_flow_injection_constraints[ n ].set_function( linear_function );
   }

   add_static_constraint( v_power_flow_injection_constraints ,
                          "HVDC_power_flow_injection" );

/*--------------------------------------------------------------------------*/

   if( !f_NetworkData->get_network_cost().empty() ) {

    if( v_power_flow_aux_var_one_constraints.size() != get_number_lines() ) {

     assert( v_power_flow_aux_var_one_constraints.empty() );
     v_power_flow_aux_var_one_constraints.resize( get_number_lines() );
    }

    for( Index line_id = 0 ; line_id < get_number_lines() ; ++line_id ) {
     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_power_flow[ line_id ] , -1.0 );
     linear_function->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
     v_power_flow_aux_var_one_constraints[ line_id ].set_lhs( 0.0 );
     v_power_flow_aux_var_one_constraints[ line_id ].set_rhs(
      Inf< double >() );

     v_power_flow_aux_var_one_constraints[ line_id ].set_function(
      linear_function );
    }

    add_static_constraint( v_power_flow_aux_var_one_constraints ,
                           "power_flow_auxiliary_variable_one" );

/*--------------------------------------------------------------------------*/

    if( v_power_flow_aux_var_two_constraints.size() != get_number_lines() ) {

     assert( v_power_flow_aux_var_two_constraints.empty() );
     v_power_flow_aux_var_two_constraints.resize( get_number_lines() );
    }

    for( Index line_id = 0 ; line_id < get_number_lines() ; ++line_id ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_power_flow[ line_id ] , 1.0 );
     linear_function->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
     v_power_flow_aux_var_two_constraints[ line_id ].set_lhs( 0.0 );
     v_power_flow_aux_var_two_constraints[ line_id ].set_rhs(
      Inf< double >() );
     v_power_flow_aux_var_two_constraints[ line_id ].set_function(
      linear_function );
    }

    add_static_constraint( v_power_flow_aux_var_two_constraints ,
                           "power_flow_auxiliary_variable_two" );

   }
  } // end HVDC_Lines constraints

/*--------------------------------------------------------------------------*/

// TODO implementation of AC and AC-HVDC lines is not ready

  if( lines_type == kAC ) {    // AC power flow limit
/*
  if( v_AC_power_flow_limit_constraints.size() != get_number_lines() ) {
   // this should only happen once
   assert( v_AC_power_flow_limit_constraints.empty() );
   v_AC_power_flow_limit_constraints.resize( get_number_lines() );
  }

  // Flow limit constraints

  for( Index line_id = 0; line_id < get_number_lines(); ++line_id ) {

    auto linear_function = new LinearFunction();
    double constant_term = 0;

    for( Index node_id = 0; node_id < get_number_nodes(); ++node_id ) {

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
  if( v_AC_HVDC_power_flow_constraints.size() != get_number_lines() ) {
   // this should only happen once
   assert( v_AC_HVDC_power_flow_constraints.empty() );
   v_AC_HVDC_power_flow_constraints.resize( get_number_lines() );
  }
  // TODO

  add_static_constraint( v_AC_HVDC_power_flow_constraints, "AC/HVDC_power_flow_limits" );
*/
  } // end AC-HVDC constraints
 }
 set_constraints_generated();
}

/*--------------------------------------------------------------------------*/

/// verifies whether the current solution is feasible for the given constraints
/** This function checks whether the relative violation of each RowConstraint
 * in the given group of RowConstraint is not greater than the provided
 * tolerance.
 *
 * @return This function returns true if and only if the relative violation of
 *         each RowConstraint in the given group is not greater than the given
 *         tolerance. */

template< class C >
static std::enable_if_t< std::is_base_of_v< RowConstraint , C > , bool >
is_feasible( std::vector< C > & constraints , double tolerance ) {
 for( auto & constraint : constraints ) {
  if( constraint.is_relaxed() )
   continue;
  constraint.compute();
  if( constraint.rel_viol() > tolerance )
   return false;
 }
 return true;
}

/*--------------------------------------------------------------------------*/

bool DCNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc ) {

 // Retrieve the tolerance.

 auto config = dynamic_cast< SimpleConfiguration< double > * >( fsbc );

 if( ( !config ) && f_BlockConfig )
  config = dynamic_cast< SimpleConfiguration< double > * >
  ( f_BlockConfig->f_is_feasible_Configuration );

 // If a tolerance has not been provided, use the default tolerance.
 const auto tolerance = config ? config->f_value : 1.0e-8;

 return NetworkBlock::is_feasible( useabstract )
        && ::is_feasible( v_AC_power_flow_limit_constraints , tolerance )
        && ::is_feasible( v_HVDC_power_flow_limit_constraints , tolerance )
        && ::is_feasible( v_AC_HVDC_power_flow_limit_constraints , tolerance )
        && ::is_feasible( v_power_flow_injection_constraints , tolerance )
        && ::is_feasible( v_AC_HVDC_power_flow_constraints , tolerance )
        && ::is_feasible( v_power_flow_aux_var_one_constraints ,
                          tolerance )
        && ::is_feasible( v_power_flow_aux_var_two_constraints ,
                          tolerance );

}  // end( DCNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_objective( Configuration * objc ) {

 // Initial check on network

 auto lines_type = f_NetworkData->get_lines_type();

/*--------------------------------------------------------------------------*/

 // HVDC power flow limit
 if( lines_type == kHVDC ) {

  if( objective_generated() )
   return; // Objective has already been generated

  if( get_objective() != nullptr )  // an objective is there already
   return;                         // cowardly (and silently) return

  auto linear_function = new LinearFunction();

  if( !f_NetworkData->get_network_cost().empty() ) {
   for( Index line_id = 0 ; line_id < get_number_lines() ; ++line_id ) {
    linear_function->add_variable( &v_auxiliary_variable[ line_id ] ,
                                   f_NetworkData->get_network_cost()[ line_id ] ,
                                   0.0 );
   }
   linear_function->set_constant_term( f_const_term );
   objective.set_function( linear_function );
   objective.set_sense( Objective::eMin );

  } else { // empty objective function
   linear_function->set_constant_term( f_const_term );
   objective.set_function( linear_function );
  }
  // Set Block objective
  this->set_objective( &objective );
 }

 // TODO The implementation of objective function for
 //  AC and AC-HVDC lines is not ready

 // AC power flow limit
 if( lines_type == kAC ) {

  // TODO
 }

 // AC-HVDC power flow limit
 if( lines_type == kAC_HVDC ) {

  // TODO
 }

 set_objective_generated();

}

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

  Index max_index = *std::max_element( std::begin( subset ) ,
                                       std::end( subset ) );
  assert( max_index < get_number_nodes() );
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
   case ( kHVDC ): {
    for( auto i : subset )
     v_power_flow_injection_constraints[ i ].set_both( -v_active_demand[ i ] );
    break;
   }
   case ( kAC ):
    // TODO
    break;
   case ( kAC_HVDC ):
    // TODO
    break;
   default:
    break;
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered )
   std::sort( subset.begin() , subset.end() );

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
    case ( kHVDC ): {
     for( Index i = rng.first ; i < rng.second ; ++i )
      v_power_flow_injection_constraints[ i ].set_both( -v_active_demand[ i ] );
     break;
    }
    case ( kAC ):
     // TODO
     break;
    case ( kAC_HVDC ):
     // TODO
     break;
    default:
     break;
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
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
