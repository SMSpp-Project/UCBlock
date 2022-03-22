/*--------------------------------------------------------------------------*/
/*--------------------- File ECNetworkBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ECNetworkBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NetworkBlock.h"
#include "ECNetworkBlock.h"
#include "LinearFunction.h"
#include "BatteryUnitBlock.h"
#include "IntermittentUnitBlock.h"
#include "ThermalUnitBlock.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( ECNetworkBlock );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

ECNetworkBlock::~ECNetworkBlock() {

 for( auto & constraint : micro_power_balance_constraints )
  constraint.clear();

 clear_constraints( power_balance_constraints );
 clear_constraints( power_flow_limit_constraints );

 objective.clear();

 delete f_NetworkData;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 static std::vector< std::string > expected_dims =
  { "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars =
  { "BuyPrice" , "ConsumptionPrice" , "SellPrice" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Optional variables

 if( !::deserialize_dim( group , "NumberIntervals" ,
                         f_number_intervals , false ) )
  f_number_intervals = 1;
 else {
  if( ( f_number_intervals < 1 ) )
   throw ( std::invalid_argument(
    "ECNetworkBlock::::deserialize: NumberIntervals must be > 0." ) );
 }

 // Mandatory variables

 Index number_nodes;
 ::deserialize_dim( group , "NumberNodes" , number_nodes , false )
 // Since the dimension "NumberNodes" must be provided since there not could
 // be an Energy Community with only one node, i.e., only one user, it means
 // that a NetworkData has always been provided. Thus, the NetworkData is
 // deserialized, and it is marked as being local.
 delete f_NetworkData;
 f_NetworkData = new NetworkData();
 f_NetworkData->deserialize( group );

 ::deserialize( group , "BuyPrice" , f_number_intervals ,
                v_buy_price , false , true );
 ::deserialize( group , "ConsumptionPrice" , f_number_intervals ,
                v_consumption_price , false , true );
 ::deserialize( group , "SellPrice" , f_number_intervals ,
                v_sell_price , false , true );

 // Deserialize data from the base class
 NetworkBlock::deserialize( group );
}

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_variables(
 Configuration * stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 auto number_nodes = f_NetworkData->get_number_nodes();

 // the power injected variables
 v_micro_power_injection.resize( number_nodes );
 for( auto & var : v_micro_power_injection )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_micro_power_injection , "micro_power_injection" );

 // the power absorbed variables
 v_micro_power_absorption.resize( number_nodes );
 for( auto & var : v_micro_power_absorption )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_micro_power_absorption , "micro_power_absorption" );

 // the power injected variables
 v_public_power_injection.resize( number_nodes );
 for( auto & var : v_public_power_injection )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_public_power_injection , "public_power_injection" );

 // the power absorbed variables
 v_public_power_absorption.resize( number_nodes );
 for( auto & var : v_public_power_absorption )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_public_power_absorption , "public_power_absorption" );

 // the max power variables
 v_max_power.resize( number_nodes );
 for( auto & var : v_max_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_max_power , "max_power" );

 set_variables_generated();
}

void ECNetworkBlock::generate_abstract_constraints(
 Configuration * stcc ) {

 if( constraints_generated() )
  return; // constraints have already been generated

 auto number_nodes = f_NetworkData->get_number_nodes();
 auto number_intervals = f_NetworkData->get_number_intervals();

 // inequality constraints

 // set that the hourly dispatch cannot go beyond the maximum dispatch
 // of the corresponding peak power period, i.e.:
 //
 //    P^{max} >= P^{POD,+} - P^{POD,-}       for all u, t       (1)
 // => P^{POD,+} - P^{POD,-} - P^{max} <= 0   for all u, t
 //
 //    P^{max} >= - [ P^{POD,+} - P^{POD,-} ]   for all u, t     (2)
 // => - P^{POD,+} + P^{POD,-} - P^{max} <= 0   for all u, t
 //
 // where:
 //
 //    P^{POD,+} = P^{P,+} + P^{M,+}
 //    P^{POD,-} = P^{P,-} + P^{M,-}

 power_flow_limit_constraints.resize(
  boost::multi_array< FRowConstraint , 3 >::extent_gen()
  [ number_nodes ][ number_intervals ][ 2 ] ); // 2 dims, i.e., the sign (+/-)

 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

  for( Index t = 0 ; t < number_intervals ; ++t ) {

   LinearFunction::v_coeff_pair vars_p;
   LinearFunction::v_coeff_pair vars_n;

   // P^{max} vars also depends from P^{M+} and P^{M-}
   // vars as specified in the paper
   if( default_config ) {
    // case (1)
    vars_p.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                      1.0 ) );
    vars_p.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                      -1.0 ) );
    // case (2)
    vars_n.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                      -1.0 ) );
    vars_n.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                      1.0 ) );
   }

   // case (1)
   vars_p.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                     1.0 ) );
   vars_p.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                     -1.0 ) );
   vars_p.push_back( std::make_pair( &v_max_power[ node_id ] , -1.0 ) );
   // case (2)
   vars_n.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                     -1.0 ) );
   vars_n.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                     1.0 ) );
   vars_n.push_back( std::make_pair( &v_max_power[ node_id ] , -1.0 ) );

   // case (1)
   power_flow_limit_constraints[ node_id ][ t ][ 0 ].set_rhs( 0.0 );
   power_flow_limit_constraints[ node_id ][ t ][ 0 ].set_lhs(
    -Inf< double >() );
   power_flow_limit_constraints[ node_id ][ t ][ 0 ].set_function(
    new LinearFunction( std::move( vars_p ) ) );

   // case (2)
   power_flow_limit_constraints[ node_id ][ t ][ 1 ].set_rhs( 0.0 );
   power_flow_limit_constraints[ node_id ][ t ][ 1 ].set_lhs(
    -Inf< double >() );
   power_flow_limit_constraints[ node_id ][ t ][ 1 ].set_function(
    new LinearFunction( std::move( vars_n ) ) );
  }
 }

 add_static_constraint( power_flow_limit_constraints ,
                        "power_flow_limit_constraints" );

/*--------------------------------------------------------------------------*/

 // equality constraints

 // set the power balance within the microgrid market/network, i.e., the
 // flows within the microgrid to have sum equal to zero:
 //
 //    P^{M,+} = P^{M,-}       for all t
 // => P^{M,+} - P^{M,-} = 0   for all t

 if( micro_power_balance_constraints.size() != number_intervals ) {
  assert( micro_power_balance_constraints.empty() );
  micro_power_balance_constraints.resize( number_intervals );
 }

 for( Index t = 0 ; t < number_intervals ; ++t ) {

  LinearFunction::v_coeff_pair vars;

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   vars.push_back(
    std::make_pair( &v_micro_power_injection[ node_id ] , 1.0 ) );
   vars.push_back(
    std::make_pair( &v_micro_power_absorption[ node_id ] , -1.0 ) );
  }

  micro_power_balance_constraints[ t ].set_both( 0.0 );
  micro_power_balance_constraints[ t ].set_function(
   new LinearFunction( std::move( vars ) ) );
 }

 add_static_constraint( micro_power_balance_constraints ,
                        "micro_power_balance_constraints" );

/*--------------------------------------------------------------------------*/

 // set the power balance, i.e.:
 //
 //    P^{POD,+} - P^{POD,-} - node_injection = - active_demand   for all t, u
 //
 // where:
 //
 //    P^{POD,+} = P^{P,+} + P^{M,+}
 //    P^{POD,-} = P^{P,-} + P^{M,-}

 power_balance_constraints.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ number_nodes ][ number_intervals ] );

 for( Index t = 0 ; t < number_intervals ; ++t ) {

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   LinearFunction::v_coeff_pair vars;

   vars.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                   -1.0 ) );
   vars.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                   -1.0 ) );
   vars.push_back( std::make_pair( &v_node_injection[ node_id ] , -1.0 ) );
   power_balance_constraints[ node_id ][ t ].set_both(
    -v_active_demand[ node_id ][ t ] );
   power_balance_constraints[ node_id ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }
 }

 add_static_constraint( power_balance_constraints ,
                        "power_balance_constraints" );
}

void ECNetworkBlock::generate_objective( Configuration * objc ) {

 if( objective_generated() )
  return; // objective has already been generated

 if( get_objective() != nullptr )  // an objective is there already
  return;                          // cowardly (and silently) return

 LinearFunction::v_coeff_pair vars;

 for( Index node_id = 0 ;
      node_id < f_NetworkData->get_number_nodes() ;
      ++node_id ) {

  for( Index t = 0 ;
       t < f_NetworkData->get_number_intervals() ;
       ++t ) {

   // net economic balance wrt the public market
   vars.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                   v_sell_price[ t ] ) );
   // v_sell_price, i.e.:
   // ( v_energy_weight[ t ] * v_time_resolution[ t ] * v_sell_price[ t ] ) /
   // pow( ( 1 + f_discount_rate ) , f_project_lifetime ) )
   vars.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                   v_sell_price[ t ] ) );
   // v_sell_price, i.e.:
   //( v_energy_weight[ t ] * v_time_resolution[ t ] * v_sell_price[ t ] ) /
   //pow( ( 1 + f_discount_rate ) , f_project_lifetime ) )
   vars.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                   -v_buy_price[ t ] ) );
   // v_buy_price, i.e.:
   // -( v_energy_weight[ t ] * v_time_resolution[ t ] * v_buy_price[ t ] ) /
   // pow( ( 1 + f_discount_rate ) , f_project_lifetime ) ) );
   vars.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                   -v_buy_price[ t ] ) );
   // v_buy_price, i.e.:
   // -( v_energy_weight[ t ] * v_time_resolution[ t ] * v_buy_price[ t ] ) /
   // pow( ( 1 + f_discount_rate ) , f_project_lifetime ) )
  }

  // the costs due to the peak power
  vars.push_back( std::make_pair( &v_max_power[ node_id ] , f_tariff );
  // f_tariff, i.e.: f_weight * f_tariff
 }

 auto lf = new LinearFunction( std::move( vars ) );
 // f_constant_term, i.e.:
 // -( v_energy_weight[ t ] * v_time_resolution[ t ] *
 //    v_consumption_price[ t ] * v_active_demand[ node_id ][ t ] ) /
 //  pow( ( 1 + f_discount_rate ) , f_project_lifetime ) );
 lf->set_constant_term( f_constant_term );
 objective.set_function( lf );
 objective.set_sense( Objective::eMax );

 // set block objective
 this->set_objective( &objective );
 set_objective_generated();
}

/*--------------------------------------------------------------------------*/

template< unsigned long T >
void ECNetworkBlock::clear_constraints(
 boost::multi_array< FRowConstraint , T > & constraints ) {
 auto constraint = constraints.data();
 auto n = constraints.num_elements();
 for( decltype( n ) i = 0 ; i < n ; ++i , ++constraint )
  constraint->clear();
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File ECNetworkBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/